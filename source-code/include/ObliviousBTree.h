#ifndef __OBLIVIOUS_BTREE_H__
#define __OBLIVIOUS_BTREE_H__

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
        int8_t* char1 = p1 + sizeof(int32_t);
        int8_t* char2 = p2 + sizeof(int32_t);
        int8_t res = *char1 - *char2;
        if (res != 0) {
            if (order == 0) return res;
            else if (order == 1) return -res;
            else return 0;
        }
    }
    else if (attrType == INTEGER) {
        int8_t* pos1 = p1 + sizeof(int32_t);
        int8_t* pos2 = p2 + sizeof(int32_t);
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
        int8_t* pos1 = p1 + sizeof(int32_t);
        int8_t* pos2 = p2 + sizeof(int32_t);
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
        char* str1 = (char *)p1 + sizeof(int32_t);
        char* str2 = (char *)p2 + sizeof(int32_t);
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
class ObliviousBTree {
public:
    ObliviousBTree(const uint32_t mode, const std::string scale, const uint32_t table_id, const Schema* in_sch, const uint32_t o_height, const std::string& in_dst, const std::string& in_prefix) {
        // Load Meta-Data
        sch = in_sch;
        btree_sch = in_sch;
        generateIndexInfo(mode, table_id, index_num, index_ids);
        root_ids = new uint32_t[index_num];
        height = new uint32_t[index_num];
        outsourced_height = new uint32_t[index_num];
        for (uint32_t i = 0; i < index_num; ++i)
            outsourced_height[i] = o_height;
        cached_height = new uint32_t[index_num];
        cached_index = new std::unordered_map <uint32_t, std::string> [index_num];
        index_oram = new ORAM*[index_num];
        
        // Import Data Blocks
        printf("\n-------------------------------------\n");
        std::unordered_map <uint32_t, std::string> blocks;
        importData(mode, scale, table_id, sch, blocks);
        block_num = blocks.size();
        printf("Total # of data blocks required: %u\n", block_num);
        
        // Build ORAM for Data Blocks
        std::unordered_map <uint32_t, std::string> data_blocks = blocks;
        std::string data_in_dst = in_dst + "_data";
        std::string data_in_prefix = in_prefix + "_data";
        data_oram = new T(data_blocks, data_in_dst.c_str(), (data_in_prefix + "_oram.txt").c_str());
        
        // Build B-Tree
        // The format of each entry is as follows:
        // entry: uint32_t pid; attr_type key;
        printf("\n-------------------------------------\n");
        printf("# of indexed columns: %u\n\n", index_num);
        
        const uint32_t blk_cont_len = B - aes_block_size - sizeof(uint32_t);
        const uint32_t blk_internal_entry_len = blk_cont_len - sizeof(char) - sizeof(uint32_t);
        const uint32_t blk_leaf_entry_len = blk_internal_entry_len - sizeof(int32_t);
        uint32_t item_size = sch->item_size;
        uint32_t item_per_blk = sch->item_per_blk;
        uint32_t item_num = block_num * item_per_blk;
        
        for (uint32_t h = 0; h < index_num; ++h) {
            std::unordered_map <uint32_t, std::string> index_blocks;
            uint32_t block_cnt = 0;
            
            uint32_t attr_id = index_ids[h];
            ATTR_TYPE attr_type = sch->attrType[attr_id];
            uint32_t attr_offset = sch->attrOffset[attr_id];
            uint32_t attr_size = sch->attrSize[attr_id];
            
            const uint32_t ele_len = sizeof(uint32_t) + attr_size;
            uint32_t max_internal_per_node = (uint32_t)(blk_internal_entry_len * btree_node_rate / ele_len);
            uint32_t max_leaf_per_node = (uint32_t)(blk_leaf_entry_len * btree_node_rate / ele_len);
            
            printf("Index #%u, Column ID #%u\n", h, attr_id);
            printf("Fanout of internal index blocks: %u\n", max_internal_per_node);
            printf("Fanout of leaf index blocks: %u\n", max_leaf_per_node);
            
            // Generate B-tree Leaf Entries
            uint32_t level_entry_num = 0;
            char* level_entry_array = new char[item_num * ele_len];
            char* pEntry = level_entry_array;
            for (uint32_t i = 0; i < block_num; ++i) {
                std::string block_value = blocks[i];
                const char* pItem = block_value.c_str() + sizeof(uint32_t);
                for (uint32_t j = 0; j < item_per_blk; ++j) {
                    const char* pAttr = pItem + attr_offset;
                    if (*pItem == 'r') {
                        ++level_entry_num;
                        memcpy(pEntry, &i, sizeof(uint32_t));
                        memcpy(pEntry + sizeof(uint32_t), pAttr, attr_size);
                        pEntry += ele_len;
                    }
                    if (j < item_per_blk - 1)
                        pItem += item_size;
                }
            }
            
            // Sort B-tree Leaf Entries
            CMP cmp {1, {attr_id}, {0}};
            btree_cmp = &cmp;
            qsort(level_entry_array, level_entry_num, ele_len, btree_compare);
            
            // Index Blocks
            char isLeaf = 1;
            uint32_t curHeight = 0;
            uint32_t max_entry_per_node = max_leaf_per_node;
            while (level_entry_num > max_entry_per_node) {
                uint32_t next_level_entry_num = (uint32_t)ceil((double)level_entry_num / max_entry_per_node);
                char* next_level_entry_array = new char[next_level_entry_num * ele_len];
                pEntry = next_level_entry_array;
                char* pArray = level_entry_array;
                uint32_t now = 0;
                while (now < level_entry_num) {
                    uint32_t key = block_cnt++;
                    uint32_t tmp_cnt = std::min(max_entry_per_node, level_entry_num - now);
                    std::string value(&isLeaf, sizeof(char));
                    if (isLeaf == 1) {
                        int32_t next = tmp_cnt + now < level_entry_num ? block_cnt : -1;
                        value += std::string((const char *)&next, sizeof(int32_t));
                    }
                    value += std::string((const char *)&tmp_cnt, sizeof(uint32_t));
                    value += std::string(pArray, tmp_cnt * ele_len);
                    now += tmp_cnt;
                    pArray += tmp_cnt * ele_len;
                    assert(blk_cont_len >= value.length());
                    uint32_t tmp_len = blk_cont_len - value.length();
                    if (tmp_len > 0) value += generate_random_block(tmp_len);
                    
                    if (curHeight < outsourced_height[h]) index_blocks[key] = value;
                    else cached_index[h][key] = value;
                    
                    memcpy(pEntry, &key, sizeof(uint32_t));
                    memcpy(pEntry + sizeof(uint32_t), pArray - attr_size, attr_size);
                    pEntry += ele_len;
                }
                delete[] level_entry_array;
                level_entry_array = next_level_entry_array;
                level_entry_num = next_level_entry_num;
                ++curHeight;
                if (isLeaf == 1) {
                    isLeaf = 0;
                    max_entry_per_node = max_internal_per_node;
                }
                printf("# of index blocks in this level: %u\n", level_entry_num);
            }
            
            // Root Block
            uint32_t key = block_cnt++;
            uint32_t tmp_cnt = level_entry_num;
            std::string value(&isLeaf, sizeof(char));
            if (isLeaf == 1) {
                int32_t next = -1;
                value += std::string((const char *)&next, sizeof(int32_t));
            }
            value += std::string((const char *)&tmp_cnt, sizeof(uint32_t));
            value += std::string(level_entry_array, tmp_cnt * ele_len);
            assert(blk_cont_len >= value.length());
            uint32_t tmp_len = blk_cont_len - value.length();
            if (tmp_len > 0) value += generate_random_block(tmp_len);
            delete[] level_entry_array;
            
            if (curHeight < outsourced_height[h]) index_blocks[key] = value;
            else cached_index[h][key] = value;
            root_ids[h] = key;
            ++curHeight;
            height[h] = curHeight;
            outsourced_height[h] = std::min(height[h], outsourced_height[h]);
            cached_height[h] = height[h] - outsourced_height[h];
            
            printf("# of index block in this level: 1\n");
            printf("Index #%u, Root block ID #%u\n", h, root_ids[h]);
            printf("B-tree height: %u\n", height[h]);
            
            // Build ORAM
            printf("Total # of index blocks required: %u\n", block_cnt);
            uint32_t outsourced_block_cnt = index_blocks.size();
            printf("Total # of cached index blocks: %u\n", block_cnt - outsourced_block_cnt);
            printf("Total # of index blocks on the server: %u\n", outsourced_block_cnt);
            
            std::string index_in_dst = in_dst + "_index_" + std::to_string(h);
            std::string index_in_prefix = in_prefix + "_index_" + std::to_string(h);
            index_oram[h] = new T(index_blocks, index_in_dst.c_str(), (index_in_prefix + "_oram.txt").c_str());
        }
        resetCommSize();
        
        // Write Information
        std::string fnode = in_prefix + "_node.txt";
        FILE* fp = fopen(fnode.c_str(), "wb");
        fwrite(&block_num, sizeof(uint32_t), 1, fp);
        fwrite(&index_num, sizeof(uint32_t), 1, fp);
        fwrite(index_ids, sizeof(uint32_t), index_num, fp);
        fwrite(root_ids, sizeof(uint32_t), index_num, fp);
        fwrite(height, sizeof(uint32_t), index_num, fp);
        fwrite(outsourced_height, sizeof(uint32_t), index_num, fp);
        fwrite(cached_height, sizeof(uint32_t), index_num, fp);
        for (uint32_t h = 0; h < index_num; ++h) {
            uint32_t cached_index_size = cached_index[h].size();
            fwrite(&cached_index_size, sizeof(uint32_t), 1, fp);
            for (auto it = cached_index[h].begin(); it != cached_index[h].end(); ++it) {
                uint32_t bid = it->first;
                fwrite(&bid, sizeof(uint32_t), 1, fp);
                fwrite(it->second.c_str(), sizeof(char), blk_cont_len, fp);
            }
        }
        fclose(fp);
        printf("Finish building B-Tree...\n");
        printf("-------------------------------------\n");
    }
    
    ObliviousBTree(const Schema* in_sch, const std::string& in_prefix) {
        sch = in_sch;
        btree_sch = in_sch;
        
        // Load Meta-Data
        std::string fnode = in_prefix + "_node.txt";
        FILE* fp = fopen(fnode.c_str(), "rb");
        fread(&block_num, sizeof(uint32_t), 1, fp);
        fread(&index_num, sizeof(uint32_t), 1, fp);
        index_ids = new uint32_t[index_num];
        fread(index_ids, sizeof(uint32_t), index_num, fp);
        root_ids = new uint32_t[index_num];
        fread(root_ids, sizeof(uint32_t), index_num, fp);
        height = new uint32_t[index_num];
        fread(height, sizeof(uint32_t), index_num, fp);
        outsourced_height = new uint32_t[index_num];
        fread(outsourced_height, sizeof(uint32_t), index_num, fp);
        cached_height = new uint32_t[index_num];
        fread(cached_height, sizeof(uint32_t), index_num, fp);
        cached_index = new std::unordered_map <uint32_t, std::string> [index_num];
        const uint32_t blk_cont_len = B - aes_block_size - sizeof(uint32_t);
        for (uint32_t h = 0; h < index_num; ++h) {
            uint32_t cached_index_size;
            fread(&cached_index_size, sizeof(uint32_t), 1, fp);
            uint32_t bid;
            char bvalue[blk_cont_len];
            for (uint32_t i = 0; i < cached_index_size; ++i) {
                fread(&bid, sizeof(uint32_t), 1, fp);
                fread(bvalue, sizeof(char), blk_cont_len, fp);
                cached_index[h][bid] = std::string(bvalue, blk_cont_len);
            }
        }
        fclose(fp);
        
        printf("Start Loading Mongo....\n");
        std::string data_in_prefix = in_prefix + "_data";
        data_oram = new T((data_in_prefix + "_oram.txt").c_str());
        index_oram = new ORAM*[index_num];
        for (uint32_t h = 0; h < index_num; ++h) {
            std::string index_in_prefix = in_prefix + "_index_" + std::to_string(h);
            index_oram[h] = new T((index_in_prefix + "_oram.txt").c_str());
        }
        resetCommSize();
        printf("Finish Loading Mongo...\n");
    }
    
    ~ObliviousBTree() {
        if (data_oram != NULL) delete data_oram;
        if (index_oram != NULL) {
            for (uint32_t i = 0; i < index_num; ++i)
                if (index_oram[i] != NULL)
                    delete index_oram[i];
            delete[] index_oram;
        }
        if (index_ids != NULL) delete[] index_ids;
        if (root_ids != NULL) delete[] root_ids;
        if (height != NULL) delete[] height;
        if (outsourced_height != NULL) delete[] outsourced_height;
        if (cached_height != NULL) delete[] cached_height;
        if (cached_index != NULL) delete[] cached_index;
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
    
    void getDummyIndexBlock(uint32_t attrID) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attrID == index_ids[h]) break;
        
        int32_t dummyID = -1;
        index_oram[h]->get(dummyID);
    }
    
    void getDummyDataBlock() {
        int32_t dummyID = -1;
        data_oram->get(dummyID);
    }
    
    // for oblivious sort-merge join
    void getFirstTuple(bool& next, uint32_t attrID, char* leaf, int32_t& leafIndex, char* data, int32_t& dataBID, int32_t& dataIndex) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attrID == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attrID];
        uint32_t attrSize = sch->attrSize[attrID];
        uint32_t attrOffset = sch->attrOffset[attrID];
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        leafIndex = -1; dataBID = -1; dataIndex = -1;
        
        // index blocks
        int32_t now = root_ids[h];
        for (int32_t i = 0; i < height[h]; ++i) {
            std::string value;
            if (i < cached_height[h]) value = cached_index[h][now];
            else value = index_oram[h]->get(now);
            const char* block = value.c_str();
            
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            int32_t tmp_cnt;
            if (isLeaf) {
                memcpy(leaf, block, blk_cont_len);
                memcpy(&tmp_cnt, block + offset, sizeof(int32_t));
            }
            else memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            assert(tmp_cnt > 0);
            if (!isLeaf) memcpy(&now, block + offset, sizeof(int32_t));
        }
        
        // data block
        offset += sizeof(uint32_t);
        char* entry = leaf + offset;
        memcpy(&now, entry, sizeof(int32_t));
        if (now < 0) getDummyDataBlock();
        else {
            leafIndex = 0;
            std::string value = data_oram->get(now);
            const char* block = value.c_str();
            memcpy(data, block, blk_cont_len);
            
            uint32_t j = 0;
            char* key = entry + sizeof(uint32_t);
            char* attr = (char *)(block + sizeof(uint32_t) + attrOffset);
            if (dataBID == now) {
                j = dataIndex + 1;
                assert(j < itemPerBlk);
                attr += j * itemSize;
            }
            else dataBID = now;
            for (; j < itemPerBlk; ++j) {
                if (attrType == TINYTEXT && strncmp(key, attr, attrSize) == 0 ||
                    attrType != TINYTEXT && memcmp(key, attr, attrSize) == 0) {
                    next = true;
                    dataIndex = j;
                    break;
                }
                if (j < itemPerBlk - 1)
                    attr += itemSize;
            }
            if (!next) dataIndex = -1;
        }
    }
    
    // for oblivious sort-merge join
    void getNextTuple(bool& next, uint32_t attrID, char* leaf, int32_t& leafIndex, char* data, int32_t& dataBID, int32_t& dataIndex) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attrID == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attrID];
        uint32_t attrSize = sch->attrSize[attrID];
        uint32_t attrOffset = sch->attrOffset[attrID];
        uint32_t entryLen = sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        next = false;
        char* pEntry = NULL;
        bool equal = false;
        
        char* test_key[2];
        char* p_leaf_offset = leaf + offset;
        test_key[0] = p_leaf_offset + 2 * sizeof(uint32_t) + leafIndex * entryLen;
        
        int32_t entryNum;
        memcpy(&entryNum, p_leaf_offset, sizeof(int32_t));
        assert(leafIndex >= 0);
        ++leafIndex;
        if (leafIndex >= entryNum) {
            int32_t next_leaf_id;
            memcpy(&next_leaf_id, leaf + sizeof(char), sizeof(int32_t));
            if (next_leaf_id >= 0) {
                std::string value;
                if (outsourced_height[h] > 0) value = index_oram[h]->get(next_leaf_id);
                else value = cached_index[h][next_leaf_id];
                const char* next_leaf = value.c_str();
                const char* p_next_leaf_offset = next_leaf + offset;
                
                // check whether the recent two keys are equal
                test_key[1] = (char *)p_next_leaf_offset + 2 * sizeof(uint32_t);
                if (attrType == TINYTEXT && strncmp(test_key[0], test_key[1], attrSize) == 0 ||
                    attrType != TINYTEXT && memcmp(test_key[0], test_key[1], attrSize) == 0)
                    equal = true;
                
                memcpy(leaf, next_leaf, blk_cont_len);
                int32_t tmp_cnt;
                memcpy(&tmp_cnt, p_leaf_offset, sizeof(int32_t));
                assert(tmp_cnt > 0);
                
                leafIndex = 0;
                pEntry = p_leaf_offset + sizeof(uint32_t);
            }
            else {
                if (outsourced_height[h] > 0) getDummyIndexBlock(attrID);
                getDummyDataBlock();
                leafIndex = -1; dataBID = -1; dataIndex = -1;
                return;
            }
        }
        else {
            if (outsourced_height[h] > 0) getDummyIndexBlock(attrID);
            pEntry = test_key[0] + attrSize;
            
            // check whether the recent two keys are equal
            test_key[1] = test_key[0] + entryLen;
            if (attrType == TINYTEXT && strncmp(test_key[0], test_key[1], attrSize) == 0 ||
                attrType != TINYTEXT && memcmp(test_key[0], test_key[1], attrSize) == 0)
                equal = true;
        }
        
        int32_t now;
        memcpy(&now, pEntry, sizeof(int32_t));
        if (now < 0) {
            getDummyDataBlock();
            leafIndex = -1; dataBID = -1; dataIndex = -1;
        }
        else {
            std::string value = data_oram->get(now);
            const char* block = value.c_str();
            memcpy(data, block, blk_cont_len);
            
            uint32_t j = 0;
            char* key = pEntry + sizeof(uint32_t);
            char* attr = (char *)(block + sizeof(uint32_t) + attrOffset);
            if (equal && dataBID == now) {
                j = dataIndex + 1;
                assert(j < itemPerBlk);
                attr += j * itemSize;
            }
            else dataBID = now;
            for (; j < itemPerBlk; ++j) {
                if (attrType == TINYTEXT && strncmp(key, attr, attrSize) == 0 ||
                    attrType != TINYTEXT && memcmp(key, attr, attrSize) == 0) {
                    next = true;
                    dataIndex = j;
                    break;
                }
                if (j < itemPerBlk - 1)
                    attr += itemSize;
            }
            if (!next) dataIndex = -1;
        }
    }
    
    // for the first table in oblivious index nested-loop join
    void getFirstTuple(bool& next, char* data, int32_t& dataBID, int32_t& dataIndex) {
        // initialize return values
        next = false;
        dataBID = -1; dataIndex = -1;
        
        assert(block_num > 0);
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        std::string value = data_oram->get(0);
        const char* block = value.c_str();
        memcpy(data, block, blk_cont_len);
        
        int32_t itemNum;
        memcpy(&itemNum, data, sizeof(int32_t));
        if (itemNum > 0) {
            char* pItem = data + sizeof(uint32_t);
            assert(*pItem == 'r');
            next = true;
            dataBID = 0;
            dataIndex = 0;
        }
    }
    
    // for the first table in oblivious index nested-loop join
    void getNextTuple(bool& next, char* data, int32_t& dataBID, int32_t& dataIndex) {
        // initialize return values
        next = false;
        
        int32_t itemNum;
        memcpy(&itemNum, data, sizeof(int32_t));
        ++dataIndex;
        if (dataIndex < itemNum) {
            getDummyDataBlock();
            
            uint32_t itemSize = sch->item_size;
            char* pItem = data + sizeof(uint32_t) + dataIndex * itemSize;
            assert(*pItem == 'r');
            next = true;
        }
        else {
            ++dataBID;
            if (dataBID >= block_num) getDummyDataBlock();
            else {
                const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
                std::string value = data_oram->get(dataBID);
                const char* block = value.c_str();
                memcpy(data, block, blk_cont_len);
                
                memcpy(&itemNum, data, sizeof(int32_t));
                if (itemNum > 0) {
                    char* pItem = data + sizeof(uint32_t);
                    assert(*pItem == 'r');
                    next = true;
                    dataIndex = 0;
                }
            }
        }
    }
    
    // TODO: "comparator" is equal to '=', 'b', or '<'
    // for the second table in oblivious index nested-loop join
    void getFirstTuple(bool& next, const char* p_attr, const ATTR_TYPE attr_type, const uint32_t attr_size, char comparator, uint32_t attr_id, char* leaf, char* data, int32_t& data_bid, int32_t* all_index, int32_t& end_level, double* band_range = NULL) {
        assert(comparator == '=' || comparator == 'b' || comparator == '<');
        
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t entryLen = sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        data_bid = -1; end_level = -1;
        memset(all_index, -1, sizeof(int32_t) * (height[h] + 1));
        
        // index blocks
        uint32_t access_num = 0;
        int32_t now = root_ids[h];
        char* pKey = NULL;
        for (int32_t i = 0; i < height[h]; ++i) {
            std::string value;
            if (i < cached_height[h]) value = cached_index[h][now];
            else {
                value = index_oram[h]->get(now);
                ++access_num;
            }
            const char* block = value.c_str();
            
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            int32_t tmp_cnt;
            if (isLeaf) {
                memcpy(leaf, block, blk_cont_len);
                memcpy(&tmp_cnt, block + offset, sizeof(int32_t));
            }
            else memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            assert(tmp_cnt > 0);
            
            bool find = false;
            char* pEntry = NULL;
            if (isLeaf) pEntry = leaf + offset + sizeof(uint32_t);
            else pEntry = (char *)(block + offset);
            for (int32_t j = 0; j < tmp_cnt; ++j) {
                char* pAttr = pEntry + sizeof(int32_t);
                int cmpres = attrCompare(p_attr, attr_type, attr_size, pAttr, attrType, attrSize, band_range);
                if (comparator == '=' || comparator == 'b') find = (cmpres <= 0);
                else if (comparator == '<') find = (cmpres < 0);
                if (find) {
                    memcpy(&now, pEntry, sizeof(int32_t));
                    if (isLeaf) pKey = pEntry + sizeof(uint32_t);
                    all_index[i] = j;
                    if (j < tmp_cnt - 1) end_level = i;
                    break;
                }
                if (j < tmp_cnt - 1)
                    pEntry += entryLen;
            }
            if (!find || now < 0) {
                while (access_num < outsourced_height[h]) {
                    getDummyIndexBlock(attr_id);
                    ++access_num;
                }
                getDummyDataBlock();
                return;
            }
        }
        
        // data block
        assert(now >= 0);
        std::string value = data_oram->get(now);
        const char* block = value.c_str();
        memcpy(data, block, blk_cont_len);
        
        data_bid = now;
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
    
    // for the second table in oblivious index nested-loop join
    void getNextTuple(bool& next, uint32_t attr_id, char* leaf, char* data, int32_t& data_bid, int32_t* all_index, int32_t& end_level) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t entryLen = sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        int32_t new_end_level = -1;
        if (end_level == -1) {
            for (int32_t i = 0; i < outsourced_height[h]; ++i)
                getDummyIndexBlock(attr_id);
            getDummyDataBlock();
            data_bid = -1;
            memset(all_index, -1, sizeof(int32_t) * (height[h] + 1));
            return;
        }
        
        // index blocks
        bool equal = false;
        char* test_key[2] = {NULL, NULL};
        test_key[0] = new char[attrSize];
        memcpy(test_key[0], leaf + offset + 2 * sizeof(uint32_t) + all_index[height[h] - 1] * entryLen, attrSize);
        
        int32_t now = root_ids[h];
        for (int32_t i = 0; i < height[h]; ++i) {
            std::string value;
            if (i < cached_height[h]) value = cached_index[h][now];
            else value = index_oram[h]->get(now);
            const char* block = value.c_str();
            
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            int32_t tmp_cnt;
            if (isLeaf) {
                memcpy(leaf, block, blk_cont_len);
                memcpy(&tmp_cnt, block + offset, sizeof(int32_t));
            }
            else memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            assert(tmp_cnt > 0);
            
            assert(all_index[i] >= 0);
            if (i == end_level) ++all_index[i];
            else if (i > end_level) all_index[i] = 0;
            if (all_index[i] < tmp_cnt - 1) new_end_level = i;
            assert(all_index[i] < tmp_cnt);
            
            char* pEntry = NULL;
            if (isLeaf) pEntry = leaf + offset + sizeof(uint32_t) + all_index[i] * entryLen;
            else pEntry = (char *)(block + offset + all_index[i] * entryLen);
            memcpy(&now, pEntry, sizeof(int32_t));
            assert(now >= 0);
            
            if (isLeaf) {
                test_key[1] = pEntry + sizeof(uint32_t);
                // check whether the recent two keys are equal
                equal = (attrType == TINYTEXT && strncmp(test_key[0], test_key[1], attrSize) == 0 ||
                    attrType != TINYTEXT && memcmp(test_key[0], test_key[1], attrSize) == 0);
                delete[] test_key[0];
            }
        }
        end_level = new_end_level;
        
        // data block
        assert(now >= 0);
        std::string value = data_oram->get(now);
        const char* block = value.c_str();
        memcpy(data, block, blk_cont_len);
        
        int32_t j = 0;
        char* pAttr = (char *)(block + sizeof(uint32_t) + attrOffset);
        if (equal && data_bid == now) {
            j = all_index[height[h]] + 1;
            assert(j < itemPerBlk);
            pAttr += j * itemSize;
        }
        else data_bid = now;
        
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
    }
    
    /***********************/
    uint32_t getDataBlockNum() const {
        return block_num;
    }
    
    size_t getServerSize() const {
        size_t server_size = data_oram->getServerSize();
        for (uint32_t h = 0; h < index_num; ++h)
            server_size += index_oram[h]->getServerSize();
        return server_size;
    }
    
    size_t getClientSize() const {
        size_t client_size = data_oram->getClientSize();
        for (uint32_t h = 0; h < index_num; ++h)
            client_size += index_oram[h]->getClientSize() + cached_index[h].size() * B;
        return client_size;
    }
    
    size_t getCommSize() const {
        size_t comm_size = data_oram->getCommSize();
        for (uint32_t h = 0; h < index_num; ++h)
            comm_size += index_oram[h]->getCommSize();
        return comm_size;
    }
    
    void resetCommSize() {
        data_oram->resetCommSize();
        for (uint32_t h = 0; h < index_num; ++h)
            index_oram[h]->resetCommSize();
    }
    
    size_t getAccessCount() const {
        size_t access_count = data_oram->getAccessCount();
        for (uint32_t h = 0; h < index_num; ++h)
            access_count += index_oram[h]->getAccessCount();
        return access_count;
    }
    
    double getReadTime() const {
        double read_time = data_oram->getReadTime();
        for (uint32_t h = 0; h < index_num; ++h)
            read_time += index_oram[h]->getReadTime();
        return read_time;
    }
    
    double getWriteTime() const {
        double write_time = data_oram->getWriteTime();
        for (uint32_t h = 0; h < index_num; ++h)
            write_time += index_oram[h]->getWriteTime();
        return write_time;
    }
    
    double getEncDecTime() const {
        double enc_dec_time = data_oram->getEncDecTime();
        for (uint32_t h = 0; h < index_num; ++h)
            enc_dec_time += index_oram[h]->getEncDecTime();
        return enc_dec_time;
    }
    
    double getORAMTime() const {
        double oram_time = data_oram->getORAMTime();
        for (uint32_t h = 0; h < index_num; ++h)
            oram_time += index_oram[h]->getORAMTime();
        return oram_time;
    }
    /***********************/
    
private:
    // for oblivious index nested-loop join
    int attrCompare(const char* pAttr1, const ATTR_TYPE attrType1, const uint32_t attrSize1, const char* pAttr2, const ATTR_TYPE attrType2, const uint32_t attrSize2, double* bRange = NULL) {
        if (attrType1 == CHAR) {
            if (*pAttr1 < *pAttr2) return -1;
            else if (*pAttr1 > *pAttr2) return 1;
            return 0;
        }
        else if (attrType1 == INTEGER || attrType1 == DOUBLE) {
            double tmp_double_1;
            if (attrType1 == INTEGER) {
                int32_t tmp_int_1;
                memcpy(&tmp_int_1, pAttr1, sizeof(int32_t));
                tmp_double_1 = tmp_int_1;
            }
            else memcpy(&tmp_double_1, pAttr1, sizeof(double));
            if (bRange != NULL)
                tmp_double_1 -= bRange[0];
            
            double tmp_double_2;
            if (attrType2 == INTEGER) {
                int32_t tmp_int_2;
                memcpy(&tmp_int_2, pAttr2, sizeof(int32_t));
                tmp_double_2 = tmp_int_2;
            }
            else memcpy(&tmp_double_2, pAttr2, sizeof(double));
            
            if (tmp_double_1 < tmp_double_2 - 1e-6) return -1;
            else if (tmp_double_1 > tmp_double_2 + 1e-6) return 1;
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
    uint32_t block_num;
    uint32_t index_num;
    uint32_t* index_ids = NULL;
    uint32_t* root_ids = NULL;
    uint32_t* height = NULL;
    uint32_t* outsourced_height = NULL;
    uint32_t* cached_height = NULL;
    std::unordered_map <uint32_t, std::string>* cached_index;
    ORAM* data_oram = NULL;
    ORAM** index_oram = NULL;
};

#endif //__OBLIVIOUS_BTREE_H__
