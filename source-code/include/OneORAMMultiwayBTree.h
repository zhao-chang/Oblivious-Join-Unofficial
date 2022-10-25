#ifndef __ONE_ORAM_MULTIWAY_BTREE_H__
#define __ONE_ORAM_MULTIWAY_BTREE_H__

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

const uint32_t internal_pre_len = 1;
const uint32_t leaf_flag_len = 2;
const uint32_t leaf_pre_len = 4;

const Schema* btree_sch = NULL;
const CMP* btree_cmp = NULL;

int btree_compare(const void* entry1, const void* entry2) {
    int8_t* tID1 = ((int8_t *)entry1) + leaf_flag_len;
    int8_t* tID2 = ((int8_t *)entry2) + leaf_flag_len;
    int8_t* bID1 = ((int8_t *)entry1) + leaf_pre_len;
    int8_t* bID2 = ((int8_t *)entry2) + leaf_pre_len;
    int8_t* attr1 = bID1 + sizeof(int32_t);
    int8_t* attr2 = bID2 + sizeof(int32_t);
    
    assert(btree_cmp->nCMPs == 1);
    uint32_t attrID = btree_cmp->attrID[0];
    ATTR_TYPE attrType = btree_sch->attrType[attrID];
    uint8_t order = btree_cmp->order[0];
    if (attrType == CHAR) {
        int8_t res = *attr1 - *attr2;
        if (res != 0) {
            if (order == 0) return res;
            else if (order == 1) return -res;
            else return 0;
        }
    }
    else if (attrType == INTEGER) {
        int32_t val1;
        int32_t val2;
        memcpy(&val1, attr1, sizeof(int32_t));
        memcpy(&val2, attr2, sizeof(int32_t));
        if (val1 != val2) {
            if (order == 0) return val1 - val2;
            else if (order == 1) return val2 - val1;
            else return 0;
        }
    }
    else if (attrType == DOUBLE) {
        double val1;
        double val2;
        memcpy(&val1, attr1, sizeof(double));
        memcpy(&val2, attr2, sizeof(double));
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
        char* str1 = (char *)attr1;
        char* str2 = (char *)attr2;
        int32_t res = strncmp(str1, str2, btree_sch->attrSize[attrID]);
        if (res != 0) {
            if (order == 0) return res;
            else if (order == 1) return -res;
            else return 0;
        }
    }
    
    int32_t blockID1;
    int32_t blockID2;
    memcpy(&blockID1, bID1, sizeof(int32_t));
    memcpy(&blockID2, bID2, sizeof(int32_t));
    if (blockID1 != blockID2) {
        if (order == 0) return blockID1 - blockID2;
        else if (order == 1) return blockID2 - blockID1;
        else return 0;
    }
    else {
        int16_t tupleID1;
        int16_t tupleID2;
        memcpy(&tupleID1, tID1, sizeof(int16_t));
        memcpy(&tupleID2, tID2, sizeof(int16_t));
        if (order == 0) return tupleID1 - tupleID2;
        else if (order == 1) return tupleID2 - tupleID1;
        else return 0;
    }
    return 0;
}

template<class T>
class OneORAMMultiwayBTree {
public:
    OneORAMMultiwayBTree(const uint32_t mode, const std::string scale, const uint32_t table_id, const Schema* in_sch, const uint32_t o_height, std::unordered_map <uint32_t, std::string>& blocks, uint32_t& block_id, const std::string& in_prefix) {
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
        
        // Import Data Blocks
        printf("\n-------------------------------------\n");
        start_data_block_id = block_id;
        importData(mode, scale, table_id, sch, blocks, block_id);
        end_data_block_id = block_id;
        block_num = end_data_block_id - start_data_block_id;
        printf("Total # of data blocks required: %u\n", block_num);
        
        // Build B-Tree
        // The format of each entry is as follows:
        // leaf entry: bool enabled; bool next; uint16_t eid; uint32_t pid; attr_type key;
        // internal entry: bool enabled; uint32_t pid; attr_type key;
        printf("\n-------------------------------------\n");
        printf("# of indexed columns: %u\n", index_num);
        
        const uint32_t blk_cont_len = B - aes_block_size - sizeof(uint32_t);
        const uint32_t blk_internal_entry_len = blk_cont_len - sizeof(char) - sizeof(uint32_t);
        const uint32_t blk_leaf_entry_len = blk_internal_entry_len - sizeof(int32_t);
        uint32_t item_size = sch->item_size;
        uint32_t item_per_blk = sch->item_per_blk;
        uint32_t item_num = block_num * item_per_blk;
        
        for (uint32_t h = 0; h < index_num; ++h) {
            uint32_t block_cnt = 0;
            
            uint32_t attr_id = index_ids[h];
            ATTR_TYPE attr_type = sch->attrType[attr_id];
            uint32_t attr_offset = sch->attrOffset[attr_id];
            uint32_t attr_size = sch->attrSize[attr_id];
            
            const uint32_t internal_ele_len = internal_pre_len + sizeof(uint32_t) + attr_size;
            const uint32_t leaf_ele_len = leaf_pre_len + sizeof(uint32_t) + attr_size;
            uint32_t max_internal_per_node = (uint32_t)(blk_internal_entry_len * btree_node_rate / internal_ele_len);
            uint32_t max_leaf_per_node = (uint32_t)(blk_leaf_entry_len * btree_node_rate / leaf_ele_len);
            
            printf("Index #%u, Column ID #%u\n", h, attr_id);
            printf("Fanout of internal index blocks: %u\n", max_internal_per_node);
            printf("Fanout of leaf index blocks: %u\n", max_leaf_per_node);
            
            // Generate B-tree Leaf Entries
            uint32_t level_entry_num = 0;
            char* level_entry_array = new char[item_num * leaf_ele_len];
            char* pEntry = level_entry_array;
            for (uint32_t i = start_data_block_id; i < end_data_block_id; ++i) {
                std::string block_value = blocks[i];
                const char* pItem = block_value.c_str() + sizeof(uint32_t);
                for (uint32_t j = 0; j < item_per_blk; ++j) {
                    const char* pAttr = pItem + attr_offset;
                    if (*pItem == 'r') {
                        uint16_t tupleID = (uint16_t)j;
                        char* pPos = pEntry + leaf_flag_len;
                        memcpy(pPos, &tupleID, sizeof(uint16_t));
                        pPos = pEntry + leaf_pre_len;
                        memcpy(pPos, &i, sizeof(uint32_t));
                        pPos += sizeof(uint32_t);
                        memcpy(pPos, pAttr, attr_size);
                        pEntry += leaf_ele_len;
                        ++level_entry_num;
                    }
                    if (j < item_per_blk - 1)
                        pItem += item_size;
                }
            }
            
            // Sort B-tree Entries
            CMP cmp {1, {attr_id}, {0}};
            btree_cmp = &cmp;
            qsort(level_entry_array, level_entry_num, leaf_ele_len, btree_compare);
            
            // Initialize Flags
            pEntry = level_entry_array;
            bool enabled = true;
            for (uint32_t i = 0; i < level_entry_num; ++i) {
                memcpy(pEntry, &enabled, sizeof(bool));
                
                char* pNextEntry = NULL;
                bool next = false;
                if (i < level_entry_num - 1) {
                    pNextEntry = pEntry + leaf_ele_len;
                    const char* pAttr = pEntry + leaf_pre_len + sizeof(uint32_t);
                    const char* pNextAttr = pNextEntry + leaf_pre_len + sizeof(uint32_t);
                    if (attr_type == TINYTEXT && strncmp(pAttr, pNextAttr, attr_size) == 0 ||
                        attr_type != TINYTEXT && memcmp(pAttr, pNextAttr, attr_size) == 0)
                        next = true;
                }
                char* pPos = pEntry + sizeof(bool);
                memcpy(pPos, &next, sizeof(bool));
                
                if (i < level_entry_num - 1)
                    pEntry = pNextEntry;
            }
            
            // Index Blocks
            char isLeaf = 1;
            uint32_t curHeight = 0;
            uint32_t max_entry_per_node = max_leaf_per_node;
            uint32_t ele_len = leaf_ele_len;
            while (level_entry_num > max_entry_per_node) {
                uint32_t next_level_entry_num = (uint32_t)ceil((double)level_entry_num / max_entry_per_node);
                char* next_level_entry_array = new char[next_level_entry_num * internal_ele_len];
                pEntry = next_level_entry_array;
                char* pArray = level_entry_array;
                uint32_t now = 0;
                while (now < level_entry_num) {
                    uint32_t key = block_id++;
                    block_cnt++;
                    uint32_t tmp_cnt = std::min(max_entry_per_node, level_entry_num - now);
                    std::string value(&isLeaf, sizeof(char));
                    if (isLeaf == 1) {
                        int32_t next = tmp_cnt + now < level_entry_num ? block_id : -1;
                        value += std::string((const char *)&next, sizeof(int32_t));
                    }
                    value += std::string((const char *)&tmp_cnt, sizeof(uint32_t));
                    value += std::string(pArray, tmp_cnt * ele_len);
                    now += tmp_cnt;
                    pArray += tmp_cnt * ele_len;
                    assert(blk_cont_len >= value.length());
                    uint32_t tmp_len = blk_cont_len - value.length();
                    if (tmp_len > 0) value += generate_random_block(tmp_len);
                    
                    if (curHeight < outsourced_height[h]) blocks[key] = value;
                    else cached_index[h][key] = value;
                    
                    char* pPos = pEntry;
                    memcpy(pPos, &enabled, sizeof(bool));
                    pPos += sizeof(bool);
                    memcpy(pPos, &key, sizeof(uint32_t));
                    pPos += sizeof(uint32_t);
                    memcpy(pPos, pArray - attr_size, attr_size);
                    pEntry += internal_ele_len;
                }
                delete[] level_entry_array;
                level_entry_array = next_level_entry_array;
                level_entry_num = next_level_entry_num;
                ++curHeight;
                if (isLeaf == 1) {
                    isLeaf = 0;
                    max_entry_per_node = max_internal_per_node;
                    ele_len = internal_ele_len;
                }
                printf("# of index blocks in this level: %u\n", level_entry_num);
            }
            
            // Root Block
            uint32_t key = block_id++;
            block_cnt++;
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
            
            if (curHeight < outsourced_height[h]) blocks[key] = value;
            else cached_index[h][key] = value;
            root_ids[h] = key;
            ++curHeight;
            height[h] = curHeight;
            outsourced_height[h] = std::min(height[h], outsourced_height[h]);
            cached_height[h] = height[h] - outsourced_height[h];
            
            printf("# of index block in this level: 1\n");
            printf("Index #%u, Root block ID #%u\n", h, root_ids[h]);
            printf("B-tree height: %u\n", height[h]);
            
            printf("Total # of index blocks required: %u\n", block_cnt);
            uint32_t cached_block_cnt = cached_index[h].size();
            printf("Total # of cached index blocks: %u\n", cached_block_cnt);
            printf("Total # of index blocks on the server: %u\n", block_cnt - cached_block_cnt);
        }
        
        // Write Information
        std::string fnode = in_prefix + "_node.txt";
        FILE* fp = fopen(fnode.c_str(), "wb");
        fwrite(&start_data_block_id, sizeof(uint32_t), 1, fp);
        fwrite(&end_data_block_id, sizeof(uint32_t), 1, fp);
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
    
    OneORAMMultiwayBTree(const Schema* in_sch, const std::string& in_prefix) {
        sch = in_sch;
        btree_sch = in_sch;
        
        // Load Meta-Data
        std::string fnode = in_prefix + "_node.txt";
        FILE* fp = fopen(fnode.c_str(), "rb");
        fread(&start_data_block_id, sizeof(uint32_t), 1, fp);
        fread(&end_data_block_id, sizeof(uint32_t), 1, fp);
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
    }
    
    ~OneORAMMultiwayBTree() {
        if (index_ids != NULL) delete[] index_ids;
        if (root_ids != NULL) delete[] root_ids;
        if (height != NULL) delete[] height;
        if (outsourced_height != NULL) delete[] outsourced_height;
        if (cached_height != NULL) delete[] cached_height;
        if (cached_index != NULL) delete[] cached_index;
    }
    
    void initORAM(ORAM* oram) {
        one_oram = oram;
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
    
    void getDummyBlock() {
        int32_t dummyID = -1;
        one_oram->get(dummyID);
    }
    
    // for the first table in multiway join
    void getFirstTuple(bool& next, char* data, int32_t& dataBID, int32_t& dataIndex) {
        // initialize return values
        next = false;
        dataBID = -1; dataIndex = -1;
        
        assert(block_num > 0);
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        std::string value = one_oram->get(start_data_block_id);
        const char* block = value.c_str();
        memcpy(data, block, blk_cont_len);
        
        int32_t itemNum;
        memcpy(&itemNum, data, sizeof(int32_t));
        if (itemNum > 0) {
            char* pItem = data + sizeof(uint32_t);
            assert(*pItem == 'r');
            next = true;
            dataBID = start_data_block_id;
            dataIndex = 0;
        }
    }
    
    // for the first table in multiway join
    void getNextTuple(bool& next, char* data, int32_t& dataBID, int32_t& dataIndex) {
        // initialize return values
        next = false;
        
        int32_t itemNum;
        memcpy(&itemNum, data, sizeof(int32_t));
        ++dataIndex;
        if (dataIndex < itemNum) {
            getDummyBlock();
            
            uint32_t itemSize = sch->item_size;
            char* pItem = data + sizeof(uint32_t) + dataIndex * itemSize;
            assert(*pItem == 'r');
            next = true;
        }
        else {
            ++dataBID;
            if (dataBID >= end_data_block_id) getDummyBlock();
            else {
                const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
                std::string value = one_oram->get(dataBID);
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
    
    // for the other tables in multiway join
    void getFirstTuple(bool& next, const char* p_attr, const ATTR_TYPE attr_type, const uint32_t attr_size, uint32_t attr_id, char** all_block, int32_t* all_bid, int32_t* all_index, bool* all_last_enabled) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t internalEntryLen = internal_pre_len + sizeof(uint32_t) + attrSize;
        uint32_t leafEntryLen = leaf_pre_len + sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        
        // index blocks
        uint32_t access_num = 0;
        int32_t now = root_ids[h];
        int16_t tupleID = -1;
        char* pKey = NULL;
        for (int32_t i = 0; i < height[h]; ++i) {
            std::string value;
            if (i < cached_height[h]) value = cached_index[h][now];
            else {
                value = one_oram->get(now);
                ++access_num;
            }
            const char* block = value.c_str();
            all_bid[i] = now;
            memcpy(all_block[i], block, blk_cont_len);
            
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            int32_t tmp_cnt;
            if (isLeaf) memcpy(&tmp_cnt, block + offset, sizeof(int32_t));
            else memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            assert(tmp_cnt > 0);
            
            bool find = false;
            char* pEntry = NULL;
            if (isLeaf) pEntry = all_block[i] + offset + sizeof(uint32_t);
            else pEntry = all_block[i] + offset;
            for (int32_t j = 0; j < tmp_cnt; ++j) {
                bool enabled;
                memcpy(&enabled, pEntry, sizeof(bool));
                if (enabled) {
                    char* pBID = NULL;
                    if (isLeaf) pBID = pEntry + leaf_pre_len;
                    else pBID = pEntry + internal_pre_len;
                    char* pAttr = pBID + sizeof(int32_t);
                    int cmpres = attrCompare(p_attr, attr_type, attr_size, pAttr, attrType, attrSize);
                    if (cmpres <= 0) {
                        find = true;
                        memcpy(&now, pBID, sizeof(int32_t));
                        all_index[i] = j;
                        char* pNextEntry = pEntry;
                        bool isLastEnabled = true;
                        for (int32_t l = j + 1; l < tmp_cnt; ++l) {
                            if (isLeaf) pNextEntry += leafEntryLen;
                            else pNextEntry += internalEntryLen;
                            bool next_enabled;
                            memcpy(&next_enabled, pNextEntry, sizeof(bool));
                            if (next_enabled) {
                                isLastEnabled = false;
                                break;
                            }
                        }
                        all_last_enabled[i] = isLastEnabled;
                        if (isLeaf) {
                            pKey = pAttr;
                            memcpy(&tupleID, pEntry + leaf_flag_len, sizeof(int16_t));
                        }
                        break;
                    }
                }
                if (j < tmp_cnt - 1) {
                    if (isLeaf) pEntry += leafEntryLen;
                    else pEntry += internalEntryLen;
                }
            }
            if (find) assert(now >= 0);
            else {
                while (access_num <= outsourced_height[h]) {
                    getDummyBlock();
                    ++access_num;
                }
                return;
            }
        }
        
        // data block
        assert(now >= 0 && tupleID >= 0);
        std::string value = one_oram->get(now);
        const char* block = value.c_str();
        all_bid[height[h]] = now;
        memcpy(all_block[height[h]], block, blk_cont_len);
        
        char* pAttr = all_block[height[h]] + sizeof(uint32_t) + attrOffset + tupleID * itemSize;
        assert(attrType == TINYTEXT && strncmp(pKey, pAttr, attrSize) == 0 ||
            attrType != TINYTEXT && memcmp(pKey, pAttr, attrSize) == 0);
        
        next = true;
        all_index[height[h]] = tupleID;
    }
    
    // the boolean return value indicates whether block "pre_leaf" must be updated
    bool disableCurrentTuple(uint32_t attr_id, char** all_block, int32_t* all_bid, int32_t* all_index, char* pre_leaf, int32_t& pre_leaf_bid) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t internalEntryLen = internal_pre_len + sizeof(uint32_t) + attrSize;
        uint32_t leafEntryLen = leaf_pre_len + sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        bool updated = false;
        
        // search the preceding enabled leaf entry and update its "next" flag, if the "next" flag of current leaf entry is "false"
        int32_t entryPos = all_index[height[h] - 1];
        char* pEntry = all_block[height[h] - 1] + offset + sizeof(uint32_t) + entryPos * leafEntryLen;
        memset(pEntry, false, sizeof(bool));
        char* pKey = pEntry + leaf_pre_len + sizeof(uint32_t);
        char* pNext = pEntry + sizeof(bool);
        bool bNext;
        memcpy(&bNext, pNext, sizeof(bool));
        if (!bNext) {
            bool exist = true;
            bool find = false;
            char* pPreEntry = pEntry;
            for (int32_t j = entryPos - 1; j >= 0; --j) {
                pPreEntry -= leafEntryLen;
                char* pPreKey = pPreEntry + leaf_pre_len + sizeof(uint32_t);
                char* pPreNext = pPreEntry + sizeof(bool);
                bool enabled;
                memcpy(&enabled, pPreEntry, sizeof(bool));
                bool equal = (attrType == TINYTEXT && strncmp(pPreKey, pKey, attrSize) == 0) ||
                    (attrType != TINYTEXT && memcmp(pPreKey, pKey, attrSize) == 0);
                if (equal && enabled) {
                    bool bPreNext;
                    memcpy(&bPreNext, pPreNext, sizeof(bool));
                    assert(bPreNext);
                    memset(pPreNext, false, sizeof(bool));
                    find = true;
                    break;
                }
                else if (!equal) {
                    exist = false;
                    break;
                }
            }
            if (pre_leaf_bid >= 0 && exist && !find) {
                int32_t entryNum;
                memcpy(&entryNum, pre_leaf + offset, sizeof(int32_t));
                pPreEntry = pre_leaf + offset + sizeof(uint32_t) + (entryNum - 1) * leafEntryLen;
                for (int32_t j = entryNum - 1; j >= 0; --j) {
                    char* pPreKey = pPreEntry + leaf_pre_len + sizeof(uint32_t);
                    char* pPreNext = pPreEntry + sizeof(bool);
                    bool enabled;
                    memcpy(&enabled, pPreEntry, sizeof(bool));
                    bool equal = (attrType == TINYTEXT && strncmp(pPreKey, pKey, attrSize) == 0) ||
                        (attrType != TINYTEXT && memcmp(pPreKey, pKey, attrSize) == 0);
                    if (equal && enabled) {
                        bool bPreNext;
                        memcpy(&bPreNext, pPreNext, sizeof(bool));
                        assert(bPreNext);
                        memset(pPreNext, false, sizeof(bool));
                        updated = true;
                        break;
                    }
                    else if (!equal) break;
                    if (j > 0) pPreEntry -= leafEntryLen;
                }
            }
        }
        
        // check if any entry in the block is still enabled
        for (int32_t i = height[h] - 1; i >= 1; --i) {
            bool find = false;
            char* pEntry = NULL;
            int32_t entryNum;
            if (i == height[h] - 1) {
                pEntry = all_block[i] + offset + sizeof(uint32_t);
                memcpy(&entryNum, all_block[i] + offset, sizeof(int32_t));
            }
            else {
                pEntry = all_block[i] + offset;
                memcpy(&entryNum, all_block[i] + sizeof(char), sizeof(int32_t));
            }
            for (int32_t j = 0; j < entryNum; ++j) {
                bool enabled = false;
                memcpy(&enabled, pEntry, sizeof(bool));
                if (enabled) {
                    find = true;
                    break;
                }
                if (j < entryNum - 1) {
                    if (i == height[h] - 1) pEntry += leafEntryLen;
                    else pEntry += internalEntryLen;
                }
            }
            if (find) break;
            else {
                int32_t entryPos = all_index[i - 1];
                char* pEntry = all_block[i - 1] + offset + entryPos * internalEntryLen;
                memset(pEntry, false, sizeof(bool));
            }
        }
        
        // write out a B-tree path
        for (int32_t i = 0; i < height[h]; ++i) {
            std::string value(all_block[i], blk_cont_len);
            if (i < cached_height[h]) cached_index[h][all_bid[i]] = value;
            else one_oram->put(all_bid[i], value);
        }
        getDummyBlock();
        return updated;
    }
    
    // write the updated "pre_leaf" block to the cloud
    void updatePreLeaf(uint32_t attr_id, char* pre_leaf, int32_t& pre_leaf_bid) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        // write out a B-tree path
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        for (int32_t i = 0; i < height[h]; ++i) {
            if (i == height[h] - 1 && pre_leaf_bid >= 0) {
                std::string value(pre_leaf, blk_cont_len);
                if (i < cached_height[h]) cached_index[h][pre_leaf_bid] = value;
                else one_oram->put(pre_leaf_bid, value);
            }
            else if (i >= cached_height[h]) getDummyBlock();
        }
        getDummyBlock();
    }
    
    bool hasNextTuple(uint32_t attr_id, char** all_block, int32_t* all_index) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t leafEntryLen = leaf_pre_len + sizeof(uint32_t) + attrSize;
        
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        char* pEntry = all_block[height[h] - 1] + offset + sizeof(uint32_t) + all_index[height[h] - 1] * leafEntryLen;
        char* pNext = pEntry + sizeof(bool);
        bool hasNext;
        memcpy(&hasNext, pNext, sizeof(bool));
        return hasNext;
    }
    
    // for the other tables in multiway join
    void getNextTuple(bool& next, uint32_t attr_id, char** all_block, int32_t* all_bid, int32_t* all_index, bool* all_last_enabled, char* pre_leaf, int32_t& pre_leaf_bid) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t attrOffset = sch->attrOffset[attr_id];
        uint32_t internalEntryLen = internal_pre_len + sizeof(uint32_t) + attrSize;
        uint32_t leafEntryLen = leaf_pre_len + sizeof(uint32_t) + attrSize;
        uint32_t itemSize = sch->item_size;
        uint32_t itemPerBlk = sch->item_per_blk;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        // initialize return values
        next = false;
        
        int32_t nextI = -1;
        for (int32_t i = height[h] - 1; i >= 0; --i) {
            if (!all_last_enabled[i]) {
                nextI = i;
                break;
            }
        }
        assert(nextI >= 0);
        
        int32_t now = root_ids[h];
        int16_t tupleID = -1;
        char* pKey = NULL;
        for (int32_t i = 0; i < height[h]; ++i) {
            if (i == height[h] - 1 && now != all_bid[i]) {
                bool find = false;
                char* pPreAttr = all_block[height[h]] + sizeof(uint32_t) + attrOffset + all_index[height[h]] * itemSize;
                uint32_t entryNum;
                memcpy(&entryNum, all_block[i] + offset, sizeof(uint32_t));
                char* pPreEntry = all_block[i] + offset + sizeof(uint32_t) + (entryNum - 1) * leafEntryLen;
                for (int32_t j = entryNum - 1; j >= 0; --j) {
                    bool enabled;
                    memcpy(&enabled, pPreEntry, sizeof(bool));
                    char* pPreKey = pPreEntry + leaf_pre_len + sizeof(uint32_t);
                    bool equal = (attrType == TINYTEXT && strncmp(pPreKey, pPreAttr, attrSize) == 0) ||
                        (attrType != TINYTEXT && memcmp(pPreKey, pPreAttr, attrSize) == 0);
                    if (equal && enabled) {
                        find = true;
                        break;
                    }
                    else if (!equal) break;
                    if (j > 0) pPreEntry -= leafEntryLen;
                }
                if (find) {
                    pre_leaf_bid = all_bid[i];
                    memcpy(pre_leaf, all_block[i], blk_cont_len);
                }
            }
            
            std::string value;
            if (i < cached_height[h]) value = cached_index[h][now];
            else value = one_oram->get(now);
            const char* block = value.c_str();
            all_bid[i] = now;
            memcpy(all_block[i], block, blk_cont_len);
            
            char isLeaf = 0;
            memcpy(&isLeaf, block, sizeof(char));
            int32_t tmp_cnt;
            if (isLeaf) memcpy(&tmp_cnt, block + offset, sizeof(int32_t));
            else memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
            assert(tmp_cnt > 0);
            
            uint32_t begin;
            if (i < nextI) begin = all_index[i];
            else if (i == nextI) begin = all_index[i] + 1;
            else begin = 0;
            char* pEntry = NULL;
            if (isLeaf) pEntry = all_block[i] + offset + sizeof(uint32_t) + begin * leafEntryLen;
            else pEntry = all_block[i] + offset + begin * internalEntryLen;
            bool find = false;
            for (uint32_t j = begin; j < tmp_cnt; ++j) {
                bool enabled;
                memcpy(&enabled, pEntry, sizeof(bool));
                if (enabled) {
                    find = true;
                    if (i < nextI) assert(j == all_index[i]);
                    else {
                        all_index[i] = j;
                        char* pNextEntry = pEntry;
                        bool isLastEnabled = true;
                        for (int32_t l = j + 1; l < tmp_cnt; ++l) {
                            if (isLeaf) pNextEntry += leafEntryLen;
                            else pNextEntry += internalEntryLen;
                            bool next_enabled;
                            memcpy(&next_enabled, pNextEntry, sizeof(bool));
                            if (next_enabled) {
                                isLastEnabled = false;
                                break;
                            }
                        }
                        all_last_enabled[i] = isLastEnabled;
                    }
                    if (isLeaf) {
                        memcpy(&now, pEntry + leaf_pre_len, sizeof(int32_t));
                        memcpy(&tupleID, pEntry + leaf_flag_len, sizeof(int16_t));
                        pKey = pEntry + leaf_pre_len + sizeof(int32_t);
                    }
                    else memcpy(&now, pEntry + internal_pre_len, sizeof(int32_t));
                    break;
                }
                if (j < tmp_cnt - 1) {
                    if (isLeaf) pEntry += leafEntryLen;
                    else pEntry += internalEntryLen;
                }
            }
            assert(find && now >= 0);
        }
        
        // data block
        assert(now >= 0 && tupleID >= 0);
        std::string value = one_oram->get(now);
        const char* block = value.c_str();
        all_bid[height[h]] = now;
        memcpy(all_block[height[h]], block, blk_cont_len);
        
        char* pAttr = all_block[height[h]] + sizeof(uint32_t) + attrOffset + tupleID * itemSize;
        assert(attrType == TINYTEXT && strncmp(pKey, pAttr, attrSize) == 0 ||
            attrType != TINYTEXT && memcmp(pKey, pAttr, attrSize) == 0);
        
        next = true;
        all_index[height[h]] = tupleID;
    }
    
    void resetFlags(uint32_t attr_id) {
        uint32_t h;
        for (h = 0; h < index_num; ++h)
            if (attr_id == index_ids[h]) break;
        
        ATTR_TYPE attrType = sch->attrType[attr_id];
        uint32_t attrSize = sch->attrSize[attr_id];
        uint32_t internalEntryLen = internal_pre_len + sizeof(uint32_t) + attrSize;
        uint32_t leafEntryLen = leaf_pre_len + sizeof(uint32_t) + attrSize;
        
        const uint32_t blk_cont_len = plain_len - sizeof(int32_t);
        uint32_t offset = sizeof(char) + sizeof(uint32_t);
        
        char block[B];
        char nextBlock[B];
        bool enabled = true;
        int32_t lBlock = root_ids[h];
        int32_t rBlock = lBlock;
        //non-leaf nodes
        for (int32_t i = 0; i < height[h] - 1; ++i) {
            int32_t nextL, nextR;
            for (int32_t now = lBlock; now <= rBlock; ++now) {
                std::string value;
                if (i < cached_height[h]) value = cached_index[h][now];
                else value = one_oram->get(now);
                memcpy(block, value.c_str(), blk_cont_len);
                char isLeaf;
                memcpy(&isLeaf, block, sizeof(char));
                assert(!isLeaf);
                int32_t tmp_cnt;
                memcpy(&tmp_cnt, block + sizeof(char), sizeof(int32_t));
                assert(tmp_cnt > 0);
                
                char* pEntry = block + offset;
                for (int32_t j = 0; j < tmp_cnt; ++j) {
                    memcpy(pEntry, &enabled, sizeof(bool));
                    if (now == lBlock && j == 0) {
                        char* pBID = pEntry + internal_pre_len;
                        memcpy(&nextL, pBID, sizeof(int32_t));
                    }
                    if (now == rBlock && j == tmp_cnt - 1) {
                        char* pBID = pEntry + internal_pre_len;
                        memcpy(&nextR, pBID, sizeof(int32_t));
                    }
                    
                    if (j < tmp_cnt - 1)
                        pEntry += internalEntryLen;
                }
                value = std::string(block, blk_cont_len);
                if (i < cached_height[h]) cached_index[h][now] = value;
                else one_oram->put(now, value);
            }
            lBlock = nextL;
            rBlock = nextR;
        }
        
        //leaf nodes
        for (int32_t now = lBlock; now <= rBlock; ++now) {
            if (now == lBlock) {
                std::string value;
                if (outsourced_height[h] <= 0) value = cached_index[h][now];
                else value = one_oram->get(now);
                memcpy(block, value.c_str(), blk_cont_len);
            }
            else memcpy(block, nextBlock, B);
            char isLeaf;
            memcpy(&isLeaf, block, sizeof(char));
            assert(isLeaf);
            int32_t tmp_cnt;
            memcpy(&tmp_cnt, block + offset, sizeof(int32_t));
            assert(tmp_cnt > 0);
            
            char* pEntry = block + offset + sizeof(uint32_t);
            for (int32_t j = 0; j < tmp_cnt; ++j) {
                memcpy(pEntry, &enabled, sizeof(bool));
                bool next = false;
                char* pNextEntry = NULL;
                if (j < tmp_cnt - 1)
                    pNextEntry = pEntry + leafEntryLen;
                else if (now < rBlock) {
                    std::string nextValue;
                    if (outsourced_height[h] <= 0) nextValue = cached_index[h][now + 1];
                    else nextValue = one_oram->get(now + 1);
                    memcpy(nextBlock, nextValue.c_str(), blk_cont_len);
                    
                    char isNextLeaf;
                    memcpy(&isNextLeaf, nextBlock, sizeof(char));
                    assert(isNextLeaf);
                    int32_t tmp_next_cnt;
                    memcpy(&tmp_next_cnt, nextBlock + offset, sizeof(int32_t));
                    assert(tmp_next_cnt > 0);
                    pNextEntry = nextBlock + offset + sizeof(uint32_t);
                }
                if (pNextEntry != NULL) {
                    char* pAttr = pEntry + leaf_pre_len + sizeof(uint32_t);
                    char* pNextAttr = pNextEntry + leaf_pre_len + sizeof(uint32_t);
                    if (attrType == TINYTEXT && strncmp(pAttr, pNextAttr, attrSize) == 0 ||
                        attrType != TINYTEXT && memcmp(pAttr, pNextAttr, attrSize) == 0)
                        next = true;
                }
                char* pPos = pEntry + sizeof(bool);
                memcpy(pPos, &next, sizeof(bool));
                
                if (j < tmp_cnt - 1)
                    pEntry += leafEntryLen;
            }
            std::string value = std::string(block, blk_cont_len);
            if (outsourced_height[h] <= 0) cached_index[h][now] = value;
            else one_oram->put(now, value);
        }
    }
    
    /***********************/
    uint32_t getDataBlockNum() const {
        return block_num;
    }
    
    size_t getClientSize() const {
        size_t client_size = 0;
        for (uint32_t h = 0; h < index_num; ++h)
            client_size += cached_index[h].size() * B;
        return client_size;
    }
    /***********************/
    
private:
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
    uint32_t block_num;
    uint32_t index_num;
    uint32_t* index_ids = NULL;
    uint32_t* root_ids = NULL;
    uint32_t* height = NULL;
    uint32_t* outsourced_height = NULL;
    uint32_t* cached_height = NULL;
    std::unordered_map <uint32_t, std::string>* cached_index;
    
    ORAM* one_oram = NULL;
    uint32_t start_data_block_id;
    uint32_t end_data_block_id;
};

#endif //__ONE_ORAM_MULTIWAY_BTREE_H__
