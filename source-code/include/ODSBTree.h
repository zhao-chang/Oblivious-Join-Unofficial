#ifndef __ODS_BTREE_H__
#define __ODS_BTREE_H__

#include "Util.h"
#include "ORAM.h"
#include "Schema.h"
#include "App.h"
#include <cstdio>
#include <cstdint>
#include <algorithm>
#include <string>
#include <queue>
#include <unordered_map>
#include <chrono>
#include <time.h>

using namespace std;

const Schema* btree_sch = NULL;
const CMP* btree_cmp = NULL;

int btree_compare(const void* entry1, const void* entry2) {
    int8_t* p1 = (int8_t *)entry1;
    int8_t* p2 = (int8_t *)entry2;
    
    assert(btree_cmp->nCMPs == 1);
    uint32_t attrID = btree_cmp->attrID[0];
    ATTR_TYPE attrType = btree_sch->attrType[attrID];
    uint8_t order = btree_cmp->order[0];
    if (attrType == CHAR) {
        int8_t* char1 = p1 + 2 * sizeof(int32_t);
        int8_t* char2 = p2 + 2 * sizeof(int32_t);
        int8_t res = *char1 - *char2;
        if (res != 0) {
            if (order == 0) return res;
            else if (order == 1) return -res;
            else return 0;
        }
    }
    else if (attrType == INTEGER) {
        int8_t* pos1 = p1 + 2 * sizeof(int32_t);
        int8_t* pos2 = p2 + 2 * sizeof(int32_t);
        int32_t val1;
        int32_t val2;
        memcpy(&val1, pos1, sizeof(int32_t));
        memcpy(&val2, pos2, sizeof(int32_t));
        if (val1 != val2) {
            if (order == 0) return val1 - val2;
            else if (order == 1) return val2 - val1;
            else return 0;
        }
    }
    else if (attrType == DOUBLE) {
        int8_t* pos1 = p1 + 2 * sizeof(int32_t);
        int8_t* pos2 = p2 + 2 * sizeof(int32_t);
        double val1;
        double val2;
        memcpy(&val1, pos1, sizeof(double));
        memcpy(&val2, pos2, sizeof(double));
        if (val1 != val2) {
            if (order == 0) {
                if (val1 - val2 < 0.0) return -1;
                else return 1;
            }
            else if (order == 1) {
                if (val2 - val1 < 0.0) return -1;
                else return 1;
            }
            else return 0;
        }
    }
    else if (attrType == STRING || attrType == TINYTEXT) {
        char* str1 = (char *)p1 + 2 * sizeof(int32_t);
        char* str2 = (char *)p2 + 2 * sizeof(int32_t);
        int32_t res = strncmp(str1, str2, btree_sch->attrSize[attrID]);
        if (res != 0) {
            if (order == 0) return res;
            else if (order == 1) return -res;
            else return 0;
        }
    }
    
    int32_t blockID1;
    int32_t blockID2;
    memcpy(&blockID1, p1, sizeof(int32_t));
    memcpy(&blockID2, p2, sizeof(int32_t));
    if (order == 0) return blockID1 - blockID2;
    else if (order == 1) return blockID2 - blockID1;
    return 0;
}

template<class T>
class ODSBTree {
public:
    ODSBTree(const uint32_t mode, const std::string& in_file, const Schema* in_sch, const std::string& in_dst, const std::string& in_prefix) {
        // Load Meta-Data
        sch = in_sch;
        btree_sch = in_sch;
        generateIndexInfo(mode, in_file, index_num, index_ids);
        root_ids = new uint32_t[index_num];
        root_pos = new uint32_t[index_num];
        height = new uint32_t[index_num];
        outsourced_height = new uint32_t[index_num];
        for (uint32_t i = 0; i < index_num; ++i)
            outsourced_height[i] = 1;
        cached_height = new uint32_t[index_num];
        
        // Calculate n_blocks needed
        printf("\n-------------------------------------\n");
        uint32_t real_item_num = 0;
        getInputItemNum(in_file, real_item_num);
        uint32_t item_size = sch->item_size;
        uint32_t item_per_blk = sch->item_per_blk;
        uint32_t block_num = (uint32_t)ceil((double)real_item_num / item_per_blk);
        uint32_t item_num = block_num * item_per_blk;
        
        const uint32_t blk_cont_len = B - aes_block_size - sizeof(uint32_t);
        const uint32_t blk_entry_len = blk_cont_len - sizeof(char) - 2 * sizeof(uint32_t);
        
        uint32_t res = block_num;
        for (uint32_t h = 0; h < index_num; ++h) {
            uint32_t attr_id = index_ids[h];
            ATTR_TYPE attr_type = sch->attrType[attr_id];
            uint32_t attr_offset = sch->attrOffset[attr_id];
            uint32_t attr_size = sch->attrSize[attr_id];
            
            const uint32_t ele_len = 2 * sizeof(uint32_t) + attr_size; // uint32_t pid; uint32_t pos; attr_type key;
            uint32_t max_entry_per_node = (uint32_t)(blk_entry_len * btree_node_rate / ele_len);
            
            uint32_t n_data = real_item_num;
            for (uint32_t outsourced_level = 0; outsourced_level < outsourced_height[h]; ++outsourced_level) {
                if (n_data <= max_entry_per_node) {
                    ++res; break;
                }
                uint32_t cur_level = (uint32_t)ceil((double)n_data / max_entry_per_node);
                res += cur_level;
                n_data = cur_level;
            }
        }
        n_blocks = (uint32_t)1 << (uint32_t)ceil(log2((double)res));
        printf("Total blocks required: %d...\n", res);
        printf("n_blocks: %d...\n", n_blocks);
        printf("\n-------------------------------------\n");
        
        // Data Blocks
        std::unordered_map <uint32_t, std::string> blocks;
        importData(in_file, sch, blocks, true, n_blocks);
        assert(block_num == blocks.size());
        
        // Build B-Tree
        printf("\n-------------------------------------\n");
        printf("# of indexed columns: %u\n", index_num);
        
        uint32_t block_cnt = block_num;
        for (uint32_t h = 0; h < index_num; ++h) {
            uint32_t attr_id = index_ids[h];
            ATTR_TYPE attr_type = sch->attrType[attr_id];
            uint32_t attr_offset = sch->attrOffset[attr_id];
            uint32_t attr_size = sch->attrSize[attr_id];
            
            const uint32_t ele_len = 2 * sizeof(uint32_t) + attr_size; // uint32_t pid; uint32_t pos; attr_type key;
            uint32_t max_entry_per_node = (uint32_t)(blk_entry_len * btree_node_rate / ele_len);
            
            printf("Index #%u, Column ID #%u\n", h, attr_id);
            printf("Fanout of internal and leaf index blocks: %u\n", max_entry_per_node);
            
            uint32_t level_entry_num = 0;
            char* level_entry_array = new char[item_num * ele_len];
            char* pEntry = level_entry_array;
            for (uint32_t i = 0; i < block_num; ++i) {
                std::string block_value = blocks[i];
                const char* pPos = block_value.c_str();
                const char* pItem = pPos + 2 * sizeof(uint32_t);
                for (uint32_t j = 0; j < item_per_blk; ++j) {
                    const char* pAttr = pItem + attr_offset;
                    if (*pItem == 'r') {
                        ++level_entry_num;
                        memcpy(pEntry, &i, sizeof(uint32_t));
                        memcpy(pEntry + sizeof(uint32_t), pPos, sizeof(uint32_t));
                        memcpy(pEntry + 2 * sizeof(uint32_t), pAttr, attr_size);
                        pEntry += ele_len;
                    }
                    if (j < item_per_blk - 1)
                        pItem += item_size;
                }
            }
            assert(level_entry_num == real_item_num);
            
            // Sort B-tree Entries
            CMP cmp {1, {attr_id}, {0}};
            btree_cmp = &cmp;
            qsort(level_entry_array, level_entry_num, ele_len, btree_compare);
            
            // Index Blocks
            char isLeaf = 1;
            uint32_t curHeight = 0;
            while (level_entry_num > max_entry_per_node) {
                uint32_t next_level_entry_num = (uint32_t)ceil((double)level_entry_num / max_entry_per_node);
                char* next_level_entry_array = new char[next_level_entry_num * ele_len];
                pEntry = next_level_entry_array;
                char* pArray = level_entry_array;
                uint32_t now = 0;
                while (now < level_entry_num) {
                    uint32_t key = block_cnt++;
                    uint32_t tmp_cnt = std::min(max_entry_per_node, level_entry_num - now);
                    uint32_t current_pos = rand_int(n_blocks);
                    std::string value((const char *)&current_pos, sizeof(uint32_t));
                    value += std::string(&isLeaf, sizeof(char));
                    value += std::string((const char *)&tmp_cnt, sizeof(uint32_t));
                    value += std::string(pArray, tmp_cnt * ele_len);
                    now += tmp_cnt;
                    pArray += tmp_cnt * ele_len;
                    
                    assert(blk_cont_len >= value.length());
                    uint32_t tmp_len = blk_cont_len - value.length();
                    if (tmp_len > 0) value += generate_random_block(tmp_len);
                    if (curHeight < outsourced_height[h]) blocks[key] = value;
                    else cached_index[key] = value;
                    
                    memcpy(pEntry, &key, sizeof(uint32_t));
                    memcpy(pEntry + sizeof(uint32_t), &current_pos, sizeof(uint32_t));
                    memcpy(pEntry + 2 * sizeof(uint32_t), pArray - attr_size, attr_size);
                    pEntry += ele_len;
                }
                delete[] level_entry_array;
                level_entry_array = next_level_entry_array;
                level_entry_num = next_level_entry_num;
                ++curHeight;
                if (isLeaf == 1) isLeaf = 0;
                
                printf("# of index blocks in this level: %u\n", level_entry_num);
            }
            
            // Root Block
            uint32_t key = block_cnt++;
            uint32_t tmp_cnt = level_entry_num;
            uint32_t current_pos = rand_int(n_blocks);
            std::string value((const char *)&current_pos, sizeof(uint32_t));
            value += std::string(&isLeaf, sizeof(char));
            value += std::string((const char *)&tmp_cnt, sizeof(uint32_t));
            value += std::string(level_entry_array, tmp_cnt * ele_len);
            
            assert(blk_cont_len >= value.length());
            uint32_t tmp_len = blk_cont_len - value.length();
            if (tmp_len > 0) value += generate_random_block(tmp_len);
            if (curHeight < outsourced_height[h]) blocks[key] = value;
            else cached_index[key] = value;
            root_ids[h] = key;
            root_pos[h] = current_pos;
            delete[] level_entry_array;
            ++curHeight;
            height[h] = curHeight;
            outsourced_height[h] = std::min(height[h], outsourced_height[h]);
            cached_height[h] = height[h] - outsourced_height[h];
            
            printf("# of index block in this level: 1\n");
            printf("Index #%u, Root block ID #%u\n", h, root_ids[h]);
            printf("B-tree height: %u\n", height[h]);
        }
        
        // Build ORAM
        printf("Total # of index and data blocks required: %u\n", block_cnt);
        uint32_t outsourced_block_cnt = blocks.size();
        assert(outsourced_block_cnt == res);
        printf("Total # of cached index blocks: %u\n", block_cnt - outsourced_block_cnt);
        printf("Total # of index and data blocks on the server: %u\n", outsourced_block_cnt);
        
        oram = new T(blocks, in_dst.c_str(), (in_prefix + "_oram.txt").c_str(), false);
        oram->resetCommSize();
        max_cache_size = 0;
        
        // Write Information
        std::string fnode = in_prefix + "_node.txt";
        FILE* fp = fopen(fnode.c_str(), "wb");
        fwrite(&index_num, sizeof(uint32_t), 1, fp);
        fwrite(index_ids, sizeof(uint32_t), index_num, fp);
        fwrite(root_ids, sizeof(uint32_t), index_num, fp);
        fwrite(root_pos, sizeof(uint32_t), index_num, fp);
        fwrite(&n_blocks, sizeof(uint32_t), 1, fp);
        fwrite(height, sizeof(uint32_t), index_num, fp);
        fwrite(outsourced_height, sizeof(uint32_t), index_num, fp);
        fwrite(cached_height, sizeof(uint32_t), index_num, fp);
        uint32_t cached_index_size = cached_index.size();
        fwrite(&cached_index_size, sizeof(uint32_t), 1, fp);
        for (auto it = cached_index.begin(); it != cached_index.end(); ++it) {
            uint32_t bid = it->first;
            fwrite(&bid, sizeof(uint32_t), 1, fp);
            fwrite(it->second.c_str(), sizeof(char), blk_cont_len, fp);
        }
        fclose(fp);
        printf("Finish building B-Tree...\n");
        printf("-------------------------------------\n");
    }
    
    ODSBTree(const Schema* in_sch, const std::string& in_prefix) {
        sch = in_sch;
        btree_sch = in_sch;
        
        // Load Meta-Data
        std::string fnode = in_prefix + "_node.txt";
        FILE* fp = fopen(fnode.c_str(), "rb");
        fread(&index_num, sizeof(uint32_t), 1, fp);
        index_ids = new uint32_t[index_num];
        fread(index_ids, sizeof(uint32_t), index_num, fp);
        root_ids = new uint32_t[index_num];
        fread(root_ids, sizeof(uint32_t), index_num, fp);
        root_pos = new uint32_t[index_num];
        fread(root_pos, sizeof(uint32_t), index_num, fp);
        fread(&n_blocks, sizeof(uint32_t), 1, fp);
        height = new uint32_t[index_num];
        fread(height, sizeof(uint32_t), index_num, fp);
        outsourced_height = new uint32_t[index_num];
        fread(outsourced_height, sizeof(uint32_t), index_num, fp);
        cached_height = new uint32_t[index_num];
        fread(cached_height, sizeof(uint32_t), index_num, fp);
        uint32_t cached_index_size;
        fread(&cached_index_size, sizeof(uint32_t), 1, fp);
        const uint32_t blk_cont_len = B - aes_block_size - sizeof(uint32_t);
        uint32_t bid;
        char bvalue[blk_cont_len];
        for (uint32_t i = 0; i < cached_index_size; ++i) {
            fread(&bid, sizeof(uint32_t), 1, fp);
            fread(bvalue, sizeof(char), blk_cont_len, fp);
            cached_index[bid] = std::string(bvalue, blk_cont_len);
        }
        fclose(fp);
        
        printf("Start Loading Mongo....\n");
        oram = new T((in_prefix + "_oram.txt").c_str());
        oram->resetCommSize();
        max_cache_size = 0;
        printf("Finish Loading Mongo...\n");
        for (uint32_t h = 0; h < index_num; ++h)
            printf("Index #%u, Root block ID #%u\n", h, root_ids[h]);
    }
    
    ~ODSBTree() {
        if (oram != NULL) delete oram;
        if (index_ids != NULL) delete[] index_ids;
        if (root_ids != NULL) delete[] root_ids;
        if (root_pos != NULL) delete[] root_pos;
        if (height != NULL) delete[] height;
        if (outsourced_height != NULL) delete[] outsourced_height;
        if (cached_height != NULL) delete[] cached_height;
    }
    
    // get the height of a B-tree
    uint32_t getHeight(uint32_t attrID) {
        for (uint32_t h = 0; h < index_num; ++h) {
            if (attrID == index_ids[h])
                return height[h];
        }
        return 0;
    }
    
    // get the outsourced height of a B-tree
    uint32_t getOutsourcedHeight(uint32_t attrID) {
        for (uint32_t h = 0; h < index_num; ++h) {
            if (attrID == index_ids[h])
                return outsourced_height[h];
        }
        return 0;
    }
    
    // for index based sort-merge join and block-nested loop join
    void getDummyTuple() {
        int32_t dummyID = -1;
        int32_t dummyPos = -1;
        std::string dummyValue = oram->getAndRemove(dummyID, dummyPos);
        oram->put(dummyID, dummyValue, dummyPos);
    }
    
    // for the first table in index based sort-merge join and block-nested loop join
    void getFirstTuple(bool& next, uint32_t attr_id, char* leaf, char* data, int32_t& data_bid, int32_t* all_index, int32_t& end_level) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t entryLen = 2 * sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - 2 * sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        data_bid = -1;
        end_level = -1;
        memset(all_index, -1, sizeof(int32_t) * (height[h] + 1));
        
        // index blocks
        uint32_t access_num = 0;
        int32_t now_id = root_ids[h];
        int32_t now_pos = root_pos[h];
        int32_t now_new_pos = rand_int(n_blocks);
        root_pos[h] = now_new_pos;
        char* pKey = NULL;
        
        std::string value;
        for (int32_t i = 0; i < height[h]; ++i) {
            if (i < cached_height[h]) value = cached_index[now_id];
            else {
                value = oram->getAndRemove(now_id, now_pos);
                ++access_num;
            }
            value.replace(0, sizeof(uint32_t), (const char *)&now_new_pos, sizeof(uint32_t));
            
            const char* block = value.c_str() + sizeof(uint32_t);
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            if (isLeaf) memcpy(leaf, block, blk_cont_len);
            int32_t tmp_cnt;
            memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            if (tmp_cnt <= 0) {
                if (i < cached_height[h]) cached_index[now_id] = value;
                else oram->put(now_id, value, now_new_pos);
                while (access_num < outsourced_height[h] + 1) {
                    getDummyTuple();
                    ++access_num;
                }
                return;
            }
            else {
                all_index[i] = 0;
                if (all_index[i] < tmp_cnt - 1) end_level = i;
                
                const char* first_id = block + offset;
                const char* first_pos = first_id + sizeof(uint32_t);
                int32_t child_id, child_pos, child_new_pos;
                memcpy(&child_id, first_id, sizeof(int32_t));
                memcpy(&child_pos, first_pos, sizeof(int32_t));
                child_new_pos = rand_int(n_blocks);
                value.replace(first_pos - value.c_str(), sizeof(uint32_t), (const char *)&child_new_pos, sizeof(uint32_t));
                if (i < cached_height[h]) cached_index[now_id] = value;
                else oram->put(now_id, value, now_new_pos);
                
                now_id = child_id;
                now_pos = child_pos;
                now_new_pos = child_new_pos;
                assert(now_id >= 0);
                
                if (isLeaf) {
                    const char* first_key = first_pos + sizeof(uint32_t);
                    pKey = leaf + (first_key - block);
                }
            }
        }
        
        // data block
        value = oram->getAndRemove(now_id, now_pos);
        value.replace(0, sizeof(uint32_t), (const char *)&now_new_pos, sizeof(uint32_t));
        oram->put(now_id, value, now_new_pos);
        const char* block = value.c_str() + sizeof(uint32_t);
        memcpy(data, block, blk_cont_len);
        data_bid = now_id;
        
        char* pAttr = (char *)(block + sizeof(uint32_t) + attrOffset);
        for (uint32_t j = 0; j < itemPerBlk; ++j) {
            if (attrType == TINYTEXT && strncmp(pKey, pAttr, attrSize) == 0 ||
                attrType != TINYTEXT && memcmp(pKey, pAttr, attrSize) == 0) {
                next = true;
                all_index[height[h]] = j;
                break;
            }
            if (j < itemPerBlk - 1)
                pAttr += itemSize;
        }
    }
    
    // for other tables in index based sort-merge join and block-nested loop join
    void getFirstTuple(bool& next, char* p_attr, ATTR_TYPE attr_type, uint32_t attr_size, uint32_t attr_id, char* leaf, char* data, int32_t& data_bid, int32_t* all_index, int32_t& end_level) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t entryLen = 2 * sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - 2 * sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        data_bid = -1;
        end_level = -1;
        memset(all_index, -1, sizeof(int32_t) * (height[h] + 1));
        
        // index blocks
        uint32_t access_num = 0;
        int32_t now_id = root_ids[h];
        int32_t now_pos = root_pos[h];
        int32_t now_new_pos = rand_int(n_blocks);
        root_pos[h] = now_new_pos;
        char* pKey = NULL;
        
        std::string value;
        for (int32_t i = 0; i < height[h]; ++i) {
            if (i < cached_height[h]) value = cached_index[now_id];
            else {
                value = oram->getAndRemove(now_id, now_pos);
                ++access_num;
            }
            value.replace(0, sizeof(uint32_t), (const char *)&now_new_pos, sizeof(uint32_t));
            
            const char* block = value.c_str() + sizeof(uint32_t);
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            if (isLeaf) memcpy(leaf, block, blk_cont_len);
            int32_t tmp_cnt;
            memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            if (tmp_cnt <= 0) {
                if (i < cached_height[h]) cached_index[now_id] = value;
                else oram->put(now_id, value, now_new_pos);
                while (access_num < outsourced_height[h] + 1) {
                    getDummyTuple();
                    ++access_num;
                }
                return;
            }
            else {
                bool find = false;
                const char* pEntry = block + offset;
                for (int32_t j = 0; j < tmp_cnt; ++j) {
                    const char* pAttr = pEntry + 2 * sizeof(int32_t);
                    int cmpres = attrCompare(p_attr, attr_type, attr_size, pAttr, attrType, attrSize);
                    if (cmpres <= 0) {
                        find = true;
                        const char* matched_id = pEntry;
                        const char* matched_pos = matched_id + sizeof(uint32_t);
                        int32_t child_id, child_pos, child_new_pos;
                        memcpy(&child_id, matched_id, sizeof(int32_t));
                        memcpy(&child_pos, matched_pos, sizeof(int32_t));
                        child_new_pos = rand_int(n_blocks);
                        value.replace(matched_pos - value.c_str(), sizeof(uint32_t), (const char *)&child_new_pos, sizeof(uint32_t));
                        if (i < cached_height[h]) cached_index[now_id] = value;
                        else oram->put(now_id, value, now_new_pos);
                        
                        now_id = child_id;
                        now_pos = child_pos;
                        now_new_pos = child_new_pos;
                        assert(now_id >= 0);
                        
                        if (isLeaf) pKey = leaf + (pAttr - block);
                        all_index[i] = j;
                        if (j < tmp_cnt - 1) end_level = i;
                        break;
                    }
                    if (j < tmp_cnt - 1)
                        pEntry += entryLen;
                }
                if (!find) {
                    if (i < cached_height[h]) cached_index[now_id] = value;
                    else oram->put(now_id, value, now_new_pos);
                    while (access_num < outsourced_height[h] + 1) {
                        getDummyTuple();
                        ++access_num;
                    }
                    return;
                }
            }
        }
        
        // data block
        value = oram->getAndRemove(now_id, now_pos);
        value.replace(0, sizeof(uint32_t), (const char *)&now_new_pos, sizeof(uint32_t));
        oram->put(now_id, value, now_new_pos);
        const char* block = value.c_str() + sizeof(uint32_t);
        memcpy(data, block, blk_cont_len);
        
        data_bid = now_id;
        char* pAttr = (char *)(block + sizeof(uint32_t) + attrOffset);
        for (int32_t j = 0; j < itemPerBlk; ++j) {
            if (attrType == TINYTEXT && strncmp(pKey, pAttr, attrSize) == 0 ||
                attrType != TINYTEXT && memcmp(pKey, pAttr, attrSize) == 0) {
                next = true;
                all_index[height[h]] = j;
                break;
            }
            if (j < itemPerBlk - 1)
                pAttr += itemSize;
        }
    }
    
    // for index based sort-merge join and block-nested loop join
    void getNextTuple(bool& next, uint32_t attr_id, char* leaf, char* data, int32_t& data_bid, int32_t* all_index, int32_t& end_level) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t entryLen = 2 * sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - 2 * sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        int32_t new_end_level = -1;
        if (end_level == -1) {
            for (int32_t i = 0; i <= outsourced_height[h]; ++i)
                getDummyTuple();
            data_bid = -1;
            memset(all_index, -1, sizeof(int32_t) * (height[h] + 1));
            return;
        }
        
        // index blocks
        char* test_key[2] = {NULL, NULL};
        test_key[0] = new char[attrSize];
        memcpy(test_key[0], leaf + offset + 2 * sizeof(uint32_t) + all_index[height[h] - 1] * entryLen, attrSize);
        
        uint32_t access_num = 0;
        int32_t now_id = root_ids[h];
        int32_t now_pos = root_pos[h];
        int32_t now_new_pos = rand_int(n_blocks);
        root_pos[h] = now_new_pos;
        
        std::string value;
        for (int32_t i = 0; i < height[h]; ++i) {
            if (i < cached_height[h]) value = cached_index[now_id];
            else {
                value = oram->getAndRemove(now_id, now_pos);
                ++access_num;
            }
            value.replace(0, sizeof(uint32_t), (const char *)&now_new_pos, sizeof(uint32_t));
            
            const char* block = value.c_str() + sizeof(uint32_t);
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            if (isLeaf) memcpy(leaf, block, blk_cont_len);
            int32_t tmp_cnt;
            memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            assert(tmp_cnt > 0);
            
            assert(all_index[i] >= 0);
            if (i == end_level) ++all_index[i];
            else if (i > end_level) all_index[i] = 0;
            if (all_index[i] < tmp_cnt - 1) new_end_level = i;
            assert(all_index[i] < tmp_cnt);
            
            const char* matched_id = block + offset + all_index[i] * entryLen;
            const char* matched_pos = matched_id + sizeof(uint32_t);
            int32_t child_id, child_pos, child_new_pos;
            memcpy(&child_id, matched_id, sizeof(int32_t));
            memcpy(&child_pos, matched_pos, sizeof(int32_t));
            child_new_pos = rand_int(n_blocks);
            value.replace(matched_pos - value.c_str(), sizeof(uint32_t), (const char *)&child_new_pos, sizeof(uint32_t));
            if (i < cached_height[h]) cached_index[now_id] = value;
            else oram->put(now_id, value, now_new_pos);
            
            now_id = child_id;
            now_pos = child_pos;
            now_new_pos = child_new_pos;
            assert(now_id >= 0);
            
            if (isLeaf) {
                const char* matched_key = matched_pos + sizeof(uint32_t);
                test_key[1] = leaf + (matched_key - block);
            }
        }
        end_level = new_end_level;
        
        // data block
        value = oram->getAndRemove(now_id, now_pos);
        value.replace(0, sizeof(uint32_t), (const char *)&now_new_pos, sizeof(uint32_t));
        oram->put(now_id, value, now_new_pos);
        const char* block = value.c_str() + sizeof(uint32_t);
        memcpy(data, block, blk_cont_len);
        
        // check whether the recent two keys are equal
        bool equal = false;
        if (attrType == TINYTEXT) {
            uint32_t attrLen[2];
            for (uint32_t i = 0; i < 2; ++i)
                attrLen[i] = std::min((uint32_t)strlen(test_key[i]), attrSize);
            if (attrLen[0] == attrLen[1] && memcmp(test_key[0], test_key[1], attrLen[0]) == 0)
                equal = true;
        }
        else if (memcmp(test_key[0], test_key[1], attrSize) == 0)
            equal = true;
        
        int32_t j = 0;
        char* pAttr = (char *)(block + sizeof(uint32_t) + attrOffset);
        if (equal && data_bid == now_id) {
            j = all_index[height[h]] + 1;
            assert(j < itemPerBlk);
            pAttr += j * itemSize;
        }
        else data_bid = now_id;
        
        char* pKey = test_key[1];
        for (; j < itemPerBlk; ++j) {
            if (attrType == TINYTEXT && strncmp(pKey, pAttr, attrSize) == 0 ||
                attrType != TINYTEXT && memcmp(pKey, pAttr, attrSize) == 0) {
                next = true;
                all_index[height[h]] = j;
                break;
            }
            if (j < itemPerBlk - 1)
                pAttr += itemSize;
        }
        delete[] test_key[0];
    }
    
    /***********************/
    size_t getServerSize() const {
        return oram->getServerSize();
    }
    
    size_t getClientSize() const {
        return oram->getClientSize() + cached_index.size() * B + max_cache_size * B;
    }
    
    size_t getCommSize() const {
        return oram->getCommSize();
    }
    
    void resetCommSize() {
        oram->resetCommSize();
    }
    
    size_t getAccessCount() const {
        return oram->getAccessCount();
    }
    
    double getReadTime() const {
        return oram->getReadTime();
    }
    
    double getWriteTime() const {
        return oram->getWriteTime();
    }
    
    double getEncDecTime() const {
        return oram->getEncDecTime();
    }
    
    double getORAMTime() const {
        return oram->getORAMTime();
    }
    /***********************/
    
private:
    // for index based block-nested loop join
    int attrCompare(const char* pAttr1, const ATTR_TYPE attrType1, const uint32_t attrSize1, const char* pAttr2, const ATTR_TYPE attrType2, const uint32_t attrSize2) {
        if (attrType1 == CHAR) {
            if (*pAttr1 < *pAttr2) return -1;
            else if (*pAttr1 > *pAttr2) return 1;
            return 0;
        }
        else if (attrType1 == INTEGER || attrType1 == DOUBLE) {
            int32_t tmp_int_1;
            double tmp_double_1;
            if (attrType1 == INTEGER) memcpy(&tmp_int_1, pAttr1, sizeof(int32_t));
            else memcpy(&tmp_double_1, pAttr1, sizeof(double));
            if (attrType2 == INTEGER) {
                int32_t tmp_int_2;
                memcpy(&tmp_int_2, pAttr2, sizeof(int32_t));
                if (attrType1 == INTEGER) {
                    if (tmp_int_1 < tmp_int_2) return -1;
                    else if (tmp_int_1 > tmp_int_2) return 1;
                }
                else {
                    if (tmp_double_1 + 1e-6 < tmp_int_2) return -1;
                    else if (tmp_double_1 - 1e-6 > tmp_int_2) return 1;
                }
            }
            else {
                double tmp_double_2;
                memcpy(&tmp_double_2, pAttr2, sizeof(double));
                if (attrType1 == INTEGER) {
                    if (tmp_int_1 < tmp_double_2 - 1e-6) return -1;
                    else if (tmp_int_1 > tmp_double_2 + 1e-6) return 1;
                }
                else {
                    if (tmp_double_1 < tmp_double_2 - 1e-6) return -1;
                    else if (tmp_double_1 > tmp_double_2 + 1e-6) return 1;
                }
            }
            return 0;
        }
        else if (attrType1 == STRING || attrType1 == TINYTEXT) {
            uint32_t attrLen1;
            if (attrType1 == TINYTEXT)
                attrLen1 = std::min((uint32_t)strlen(pAttr1), attrSize1);
            else attrLen1 = attrSize1;
            uint32_t attrLen2;
            if (attrType2 == TINYTEXT)
                attrLen2 = std::min((uint32_t)strlen(pAttr2), attrSize2);
            else attrLen2 = attrSize2;
            
            int res = strncmp(pAttr1, pAttr2, std::min(attrLen1, attrLen2));
            if (res == 0) {
                if (attrLen1 < attrLen2) res = -1;
                else if (attrLen1 > attrLen2) res = 1;
            }
            return res;
        }
        return 0;
    }
    
    const Schema* sch;
    uint32_t index_num;
    uint32_t* index_ids = NULL;
    uint32_t* root_ids = NULL;
    uint32_t* root_pos = NULL;
    uint32_t n_blocks;
    size_t max_cache_size;
    uint32_t* height = NULL;
    uint32_t* outsourced_height = NULL;
    uint32_t* cached_height = NULL;
    std::unordered_map <uint32_t, std::string> cached_index;
    ORAM* oram = NULL;
};

#endif // __ODS_BTREE_H__
