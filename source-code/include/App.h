#ifndef __APP_H__
#define __APP_H__

#include "Schema.h"
#include "Util.h"
#include "FileSimulator.h"
#include "Basic.h"
#include "ObliviousSort.h"
#include "ObliviousCompact.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <assert.h>
#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_set>
#include <unordered_map>

// mode: 0: rankings&uservisits; 1: TPC-H; 2: twitter
uint32_t getTableNum(const uint32_t mode) {
    const uint32_t n_mode = 3;
    const uint32_t table_num[n_mode] = {3, 9, 4};
    return table_num[mode];
}

std::string getCollectionPrefix(uint32_t mode, const std::string method, const std::string scale) {
    const uint32_t n_mode = 3;
    const std::string file_prefix[n_mode] = {"rankings_uservisits", "tpch", "twitter"};
    std::string predix = file_prefix[mode] +  '_' + method + '_' + scale;
    return predix;
}

std::string getCollectionPrefix(const uint32_t mode, const std::string method, const std::string scale, const uint32_t outsourced_height) {
    const uint32_t n_mode = 3;
    const std::string file_prefix[n_mode] = {"rankings_uservisits", "tpch", "twitter"};
    std::string predix = file_prefix[mode] + '_' + method + '_' + scale + '_' + std::to_string(outsourced_height);
    return predix;
}

void getFileName(const uint32_t mode, const std::string scale, const uint32_t table_id, std::string &input_file) {
    const uint32_t n_mode = 3;
    const std::string data_prefix = "data/";
    const std::string file_prefix[n_mode] = {"rankings_uservisits/", "tpch/", "twitter/"};
    const std::string file_suffix[n_mode] = {"csv", "tbl", "txt"};
    const std::string file_name[n_mode][MAX_TBLS] = {{"rankings", "uservisits", "test"},
        {"supplier", "part", "partsupp", "customer", "orders", "lineitem", "nation", "region", "nation"},
        {"inactive_user", "inactive_user", "popular_user", "normal_user"}};
    input_file = data_prefix + file_prefix[mode] + file_name[mode][table_id] + '_' + scale + '.' + file_suffix[mode];
}

void generateInputSchema(const uint32_t mode, const uint32_t table_id, Schema* sch, const bool ods = false) {
    uint32_t nAttrs[MAX_TBLS];
    ATTR_TYPE inAttr[MAX_TBLS][MAX_COLS];
    uint32_t inLen[MAX_TBLS][MAX_COLS];
    const uint32_t tot_file_num = getTableNum(mode);
    if (mode == 0) {
        const uint32_t n_attrs[tot_file_num] {4, 10, 10};
        ATTR_TYPE in_attr[tot_file_num][MAX_COLS] {{CHAR, TINYTEXT, INTEGER, INTEGER},
            {CHAR, TINYTEXT, TINYTEXT, INTEGER, DOUBLE, TINYTEXT, STRING, STRING, TINYTEXT, INTEGER},
            {CHAR, TINYTEXT, TINYTEXT, INTEGER, DOUBLE, TINYTEXT, STRING, STRING, TINYTEXT, INTEGER}};
        uint32_t in_len[tot_file_num][MAX_COLS] {{1, 127, 4, 4},
            {1, 15, 127, 4, 8, 255, 3, 5, 32, 4}, {1, 15, 127, 4, 8, 255, 3, 5, 32, 4}};
        memcpy(nAttrs, n_attrs, tot_file_num * sizeof(uint32_t));
        memcpy(inAttr[0], in_attr[0], tot_file_num * MAX_COLS * sizeof(ATTR_TYPE));
        memcpy(inLen[0], in_len[0], tot_file_num * MAX_COLS * sizeof(uint32_t));
    }
    else if (mode == 1) {
        const uint32_t n_attrs[tot_file_num] {8, 10, 6, 9, 10, 17, 5, 4, 5};
        ATTR_TYPE in_attr[tot_file_num][MAX_COLS] {{CHAR, INTEGER, STRING, TINYTEXT, INTEGER, STRING, DOUBLE, TINYTEXT},
            {CHAR, INTEGER, TINYTEXT, STRING, STRING, TINYTEXT, INTEGER, TINYTEXT, DOUBLE, TINYTEXT},
            {CHAR, INTEGER, INTEGER, INTEGER, DOUBLE, TINYTEXT},
            {CHAR, INTEGER, STRING, TINYTEXT, INTEGER, STRING, DOUBLE, TINYTEXT, TINYTEXT},
            {CHAR, INTEGER, INTEGER, CHAR, DOUBLE, STRING, TINYTEXT, STRING, CHAR, TINYTEXT},
            {CHAR, INTEGER, INTEGER, INTEGER, INTEGER, INTEGER, DOUBLE, DOUBLE, DOUBLE, CHAR, CHAR, STRING, STRING, STRING, TINYTEXT, TINYTEXT, TINYTEXT},
            {CHAR, INTEGER, TINYTEXT, INTEGER, TINYTEXT}, {CHAR, INTEGER, TINYTEXT, TINYTEXT},
            {CHAR, INTEGER, TINYTEXT, INTEGER, TINYTEXT}};
        uint32_t in_len[tot_file_num][MAX_COLS] {{1, 4, 18, 40, 4, 15, 8, 100},
            {1, 4, 54, 14, 8, 25, 4, 10, 8, 22}, {1, 4, 4, 4, 8, 198},
            {1, 4, 18, 40, 4, 15, 8, 10, 116}, {1, 4, 4, 1, 8, 10, 15, 15, 1, 78},
            {1, 4, 4, 4, 4, 4, 8, 8, 8, 1, 1, 10, 10, 10, 17, 7, 43},
            {1, 4, 14, 4, 114}, {1, 4, 11, 115}, {1, 4, 14, 4, 114}};
        memcpy(nAttrs, n_attrs, tot_file_num * sizeof(uint32_t));
        memcpy(inAttr[0], in_attr[0], tot_file_num * MAX_COLS * sizeof(ATTR_TYPE));
        memcpy(inLen[0], in_len[0], tot_file_num * MAX_COLS * sizeof(uint32_t));
    }
    else if (mode == 2) {
        const uint32_t n_attrs[tot_file_num] {3, 3, 3, 3};
        ATTR_TYPE in_attr[tot_file_num][MAX_COLS] {{CHAR, INTEGER, INTEGER},
            {CHAR, INTEGER, INTEGER}, {CHAR, INTEGER, INTEGER}, {CHAR, INTEGER, INTEGER}};
        uint32_t in_len[tot_file_num][MAX_COLS] {{1, 4, 4}, {1, 4, 4}, {1, 4, 4}, {1, 4, 4}};
        memcpy(nAttrs, n_attrs, tot_file_num * sizeof(uint32_t));
        memcpy(inAttr[0], in_attr[0], tot_file_num * MAX_COLS * sizeof(ATTR_TYPE));
        memcpy(inLen[0], in_len[0], tot_file_num * MAX_COLS * sizeof(uint32_t));
    }
    
    uint32_t index = table_id;
    sch->nAttrs = nAttrs[index];
    memcpy(sch->attrType, inAttr[index], sizeof(ATTR_TYPE) * nAttrs[index]);
    uint32_t totSize = 0;
    for (uint32_t j = 0; j < nAttrs[index]; ++j) {
        uint32_t curSize = 0;
        if (inAttr[index][j] == CHAR)
            curSize = 1;
        else if (inAttr[index][j] == INTEGER)
            curSize = 4;
        else if (inAttr[index][j] == DOUBLE)
            curSize = 8;
        else if (inAttr[index][j] == STRING || inAttr[index][j] == TINYTEXT)
            curSize = inLen[index][j];
        
        sch->attrName[j] = (index << 16) + j;
        sch->attrOffset[j] = totSize;
        sch->attrSize[j] = curSize;
        totSize += curSize;
    }
    sch->item_size = totSize;
    if (ods) sch->item_per_blk = (plain_len - META_BLOCK_SIZE - sizeof(uint32_t)) / totSize;
    else sch->item_per_blk = (plain_len - META_BLOCK_SIZE) / totSize;
}

void getInputItemNum(const uint32_t mode, const std::string scale, const uint32_t table_id, uint32_t& real_item_num) {
    std::string first = std::to_string(table_id) + '_' + scale;
    std::unordered_map <std::string, uint32_t> input_item_num;
    if (mode == 0) {
        input_item_num["0_1x"] = 360000;
        input_item_num["1_1x"] = 350000;
        input_item_num["1_0.0001x"] = 35;
        input_item_num["2_1x"] = 35;
    }
    else if (mode == 1) {
        input_item_num["0_0.01x"] = 100;
        input_item_num["1_0.01x"] = 2000;
        input_item_num["2_0.01x"] = 8000;
        input_item_num["3_0.01x"] = 1500;
        input_item_num["4_0.01x"] = 15000;
        input_item_num["5_0.01x"] = 60175;
        input_item_num["6_0.01x"] = 25;
        input_item_num["7_0.01x"] = 5;
        input_item_num["8_0.01x"] = 25;
        
        input_item_num["0_0.02x"] = 200;
        input_item_num["1_0.02x"] = 4000;
        input_item_num["2_0.02x"] = 16000;
        input_item_num["3_0.02x"] = 3000;
        input_item_num["4_0.02x"] = 30000;
        input_item_num["5_0.02x"] = 120515;
        input_item_num["6_0.02x"] = 25;
        input_item_num["7_0.02x"] = 5;
        input_item_num["8_0.02x"] = 25;
        
        input_item_num["0_0.05x"] = 500;
        input_item_num["1_0.05x"] = 10000;
        input_item_num["2_0.05x"] = 40000;
        input_item_num["3_0.05x"] = 7500;
        input_item_num["4_0.05x"] = 75000;
        input_item_num["5_0.05x"] = 299814;
        input_item_num["6_0.05x"] = 25;
        input_item_num["7_0.05x"] = 5;
        input_item_num["8_0.05x"] = 25;
        
        input_item_num["0_0.1x"] = 1000;
        input_item_num["1_0.1x"] = 20000;
        input_item_num["2_0.1x"] = 80000;
        input_item_num["3_0.1x"] = 15000;
        input_item_num["4_0.1x"] = 150000;
        input_item_num["5_0.1x"] = 600572;
        input_item_num["6_0.1x"] = 25;
        input_item_num["7_0.1x"] = 5;
        input_item_num["8_0.1x"] = 25;
        
        input_item_num["0_0.2x"] = 2000;
        input_item_num["1_0.2x"] = 40000;
        input_item_num["2_0.2x"] = 160000;
        input_item_num["3_0.2x"] = 30000;
        input_item_num["4_0.2x"] = 300000;
        input_item_num["5_0.2x"] = 1199969;
        input_item_num["6_0.2x"] = 25;
        input_item_num["7_0.2x"] = 5;
        input_item_num["8_0.2x"] = 25;
        
        input_item_num["0_0.5x"] = 5000;
        input_item_num["1_0.5x"] = 100000;
        input_item_num["2_0.5x"] = 400000;
        input_item_num["3_0.5x"] = 75000;
        input_item_num["4_0.5x"] = 750000;
        input_item_num["5_0.5x"] = 2999671;
        input_item_num["6_0.5x"] = 25;
        input_item_num["7_0.5x"] = 5;
        input_item_num["8_0.5x"] = 25;
        
        input_item_num["0_1x"] = 10000;
        input_item_num["1_1x"] = 200000;
        input_item_num["2_1x"] = 800000;
        input_item_num["3_1x"] = 150000;
        input_item_num["4_1x"] = 1500000;
        input_item_num["5_1x"] = 6001215;
        input_item_num["6_1x"] = 25;
        input_item_num["7_1x"] = 5;
        input_item_num["8_1x"] = 25;
    }
    else if (mode == 2) {
        input_item_num["0_5k"] = 117;
        input_item_num["1_5k"] = 117;
        input_item_num["2_5k"] = 29383;
        input_item_num["3_5k"] = 62087;
        
        input_item_num["0_10k"] = 239;
        input_item_num["1_10k"] = 239;
        input_item_num["2_10k"] = 43931;
        input_item_num["3_10k"] = 112635;
        
        input_item_num["0_20k"] = 429;
        input_item_num["1_20k"] = 429;
        input_item_num["2_20k"] = 82779;
        input_item_num["3_20k"] = 226780;
        
        input_item_num["0_50k"] = 1007;
        input_item_num["1_50k"] = 1007;
        input_item_num["2_50k"] = 185972;
        input_item_num["3_50k"] = 607131;
        
        input_item_num["0_100k"] = 2005;
        input_item_num["1_100k"] = 2005;
        input_item_num["2_100k"] = 483905;
        input_item_num["3_100k"] = 1188739;
        
        input_item_num["0_200k"] = 4001;
        input_item_num["1_200k"] = 4001;
        input_item_num["2_200k"] = 1626569;
        input_item_num["3_200k"] = 2327154;
        
        input_item_num["0_500k"] = 10187;
        input_item_num["1_500k"] = 10187;
        input_item_num["2_500k"] = 12021336;
        input_item_num["3_500k"] = 5662487;
    }
    real_item_num = input_item_num[first];
}

void generateIndexInfo(const uint32_t mode, const uint32_t table_id, uint32_t& index_num, uint32_t* &index_ids) {
    uint32_t index = table_id;
    const uint32_t tot_file_num = getTableNum(mode);
    if (mode == 0) {
        const uint32_t indexNum[tot_file_num] = {2, 2, 2};
        const uint32_t indexIDs[tot_file_num][2] = {{1, 2}, {2, 3}, {2, 3}};
        
        index_num = indexNum[index];
        index_ids = new uint32_t[index_num];
        memcpy(index_ids, indexIDs[index], sizeof(uint32_t) * index_num);
    }
    else if (mode == 1) {
        const uint32_t indexNum[tot_file_num] = {3, 2, 2, 2, 2, 3, 2, 1, 2};
        const uint32_t indexIDs[tot_file_num][3] = {{1, 4, 6}, {1, 8}, {1, 2}, {1, 4}, {1, 2}, {1, 2, 3}, {1, 3}, {1}, {1, 3}};
        
        index_num = indexNum[index];
        index_ids = new uint32_t[index_num];
        memcpy(index_ids, indexIDs[index], sizeof(uint32_t) * index_num);
    }
    if (mode == 2) {
        const uint32_t indexNum[tot_file_num] = {2, 2, 2, 2};
        const uint32_t indexIDs[tot_file_num][2] = {{1, 2}, {1, 2}, {1, 2}, {1, 2}};
        
        index_num = indexNum[index];
        index_ids = new uint32_t[index_num];
        memcpy(index_ids, indexIDs[index], sizeof(uint32_t) * index_num);
    }
}

void importData(const std::string& in_file, const Schema* sch, std::unordered_map <uint32_t, std::string>& blocks, const bool ods = false, const uint32_t n_blocks = 0) {
    char line[plain_len];
    char item[plain_len];
    item[0] = 'r';
    char* pItem = &(item[1]);
    
    uint32_t block_id = 0;
    uint32_t item_num = 0;
    std::string value = "";
    int32_t value_size = plain_len - sizeof(uint32_t);
    
    char split_char = ',';
    char end_char = '\n';
    const std::string suffix = in_file.substr(in_file.length() - 3);
    if (suffix == std::string("csv")) {
        split_char = ',';
        end_char = '\n';
    }
    else if (suffix == std::string("tbl")) {
        split_char = '|';
        end_char = '|';
    }
    else if (suffix == std::string("txt")) {
        split_char = '\t';
        end_char = '\n';
    }
    
    FILE* fp = fopen(in_file.c_str(), "r");
    while (fgets(line, plain_len, fp) != NULL) {
        char* start = line;
        char* pAttr = pItem;
        for (uint32_t l = 1; l < sch->nAttrs; ++l) {
            char* next = NULL;
            if (l < sch->nAttrs - 1)
                next = strchr(start, split_char);
            else next = strchr(start, end_char);
            
            if (sch->attrType[l] == CHAR) {
                *pAttr = *start;
                ++pAttr;
            }
            else if (sch->attrType[l] == INTEGER) {
                int32_t val;
                sscanf(start, "%d", &val);
                memcpy(pAttr, &val, sizeof(int32_t));
                pAttr += sizeof(int32_t);
            }
            else if (sch->attrType[l] == DOUBLE) {
                double val;
                sscanf(start, "%lf", &val);
                memcpy(pAttr, &val, sizeof(double));
                pAttr += sizeof(double);
            }
            else if (sch->attrType[l] == STRING || sch->attrType[l] == TINYTEXT) {
                uint32_t textLen = next - start;
                memcpy(pAttr, start, textLen);
                if (textLen < sch->attrSize[l])
                    pAttr[textLen] = '\0';
                pAttr += sch->attrSize[l];
            }
            
            if (l < sch->nAttrs - 1) {
                ++next;
                start = next;
            }
        }
        value += std::string(item, sch->item_size);
        ++item_num;
        if (item_num >= sch->item_per_blk) {
            value = std::string((const char*)&item_num, sizeof(uint32_t)) + value;
            if (ods) {
                uint32_t current_pos = rand_int(n_blocks);
                value = std::string((const char*)&current_pos, sizeof(uint32_t)) + value;
            }
            int32_t tmp_len = value_size - value.length();
            assert(tmp_len >= 0);
            if (tmp_len > 0) value += generate_random_block(tmp_len);
            blocks[block_id] = value;
            ++block_id;
            item_num = 0;
            value = "";
        }
    }
    fclose(fp);
    
    if (item_num > 0) {
        item[0] = 'd';
        for (uint32_t j = item_num; j < sch->item_per_blk; ++j) {
            std::string rand_block = generate_random_block(sch->item_size - 1);
            memcpy(pItem, rand_block.c_str(), sch->item_size - 1);
            value += std::string(item, sch->item_size);
        }
        value = std::string((const char*)&item_num, sizeof(uint32_t)) + value;
        if (ods) {
            uint32_t current_pos = rand_int(n_blocks);
            value = std::string((const char*)&current_pos, sizeof(uint32_t)) + value;
        }
        int32_t tmp_len = value_size - value.length();
        assert(tmp_len >= 0);
        if (tmp_len > 0) value += generate_random_block(tmp_len);
        blocks[block_id] = value;
        ++block_id;
        item_num = 0;
        value = "";
    }
}

void importData(const std::string& in_file, const Schema* sch, std::unordered_map <uint32_t, std::string>& blocks, uint32_t& block_id, const bool ods = false, const uint32_t n_blocks = 0) {
    char line[plain_len];
    char item[plain_len];
    item[0] = 'r';
    char* pItem = &(item[1]);
    
    uint32_t item_num = 0;
    std::string value = "";
    int32_t value_size = plain_len - sizeof(uint32_t);
    
    char split_char = ',';
    char end_char = '\n';
    const std::string suffix = in_file.substr(in_file.length() - 3);
    if (suffix == std::string("csv")) {
        split_char = ',';
        end_char = '\n';
    }
    else if (suffix == std::string("tbl")) {
        split_char = '|';
        end_char = '|';
    }
    else if (suffix == std::string("txt")) {
        split_char = '\t';
        end_char = '\n';
    }
    
    FILE* fp = fopen(in_file.c_str(), "r");
    while (fgets(line, plain_len, fp) != NULL) {
        char* start = line;
        char* pAttr = pItem;
        for(uint32_t l = 1; l < sch->nAttrs; ++l) {
            char* next = NULL;
            if (l < sch->nAttrs - 1)
                next = strchr(start, split_char);
            else next = strchr(start, end_char);
            
            if (sch->attrType[l] == CHAR) {
                *pAttr = *start;
                ++pAttr;
            }
            else if (sch->attrType[l] == INTEGER) {
                int32_t val;
                sscanf(start, "%d", &val);
                memcpy(pAttr, &val, sizeof(int32_t));
                pAttr += sizeof(int32_t);
            }
            else if (sch->attrType[l] == DOUBLE) {
                double val;
                sscanf(start, "%lf", &val);
                memcpy(pAttr, &val, sizeof(double));
                pAttr += sizeof(double);
            }
            else if (sch->attrType[l] == STRING || sch->attrType[l] == TINYTEXT) {
                uint32_t textLen = next - start;
                memcpy(pAttr, start, textLen);
                if (textLen < sch->attrSize[l])
                    pAttr[textLen] = '\0';
                pAttr += sch->attrSize[l];
            }
            
            if (l < sch->nAttrs - 1) {
                ++next;
                start = next;
            }
        }
        value += std::string(item, sch->item_size);
        ++item_num;
        if (item_num >= sch->item_per_blk) {
            value = std::string((const char*)&item_num, sizeof(uint32_t)) + value;
            if (ods) {
                uint32_t current_pos = rand_int(n_blocks);
                value = std::string((const char*)&current_pos, sizeof(uint32_t)) + value;
            }
            int32_t tmp_len = value_size - value.length();
            assert(tmp_len >= 0);
            if (tmp_len > 0) value += generate_random_block(tmp_len);
            blocks[block_id] = value;
            ++block_id;
            item_num = 0;
            value = "";
        }
    }
    fclose(fp);
    
    if (item_num > 0) {
        item[0] = 'd';
        for (uint32_t j = item_num; j < sch->item_per_blk; ++j) {
            std::string rand_block = generate_random_block(sch->item_size - 1);
            memcpy(pItem, rand_block.c_str(), sch->item_size - 1);
            value += std::string(item, sch->item_size);
        }
        value = std::string((const char*)&item_num, sizeof(uint32_t)) + value;
        if (ods) {
            uint32_t current_pos = rand_int(n_blocks);
            value = std::string((const char*)&current_pos, sizeof(uint32_t)) + value;
        }
        int32_t tmp_len = value_size - value.length();
        assert(tmp_len >= 0);
        if (tmp_len > 0) value += generate_random_block(tmp_len);
        blocks[block_id] = value;
        ++block_id;
        item_num = 0;
        value = "";
    }
}

void importData(const uint32_t mode, const std::string scale, const uint32_t table_id, const Schema* sch, std::unordered_map <uint32_t, std::string>& blocks, const bool ods = false, const uint32_t n_blocks = 0) {
    std::string in_file;
    getFileName(mode, scale, table_id, in_file);
    importData(in_file, sch, blocks, ods, n_blocks);
}

void importData(const uint32_t mode, const std::string scale, const uint32_t table_id, const Schema* sch, std::unordered_map <uint32_t, std::string>& blocks, uint32_t& block_id, const bool ods = false, const uint32_t n_blocks = 0) {
    std::string in_file;
    getFileName(mode, scale, table_id, in_file);
    importData(in_file, sch, blocks, block_id, ods, n_blocks);
}

// for OQP paper
void generateOutputSchema(const Schema* inSchema, const uint32_t nInputs, const int32_t* attrID, Schema* resSchema, bool equi = true) {
    uint32_t resNAttrs = inSchema[0].nAttrs;
    for (uint32_t i = 1; i < nInputs; ++i) {
        if (equi) resNAttrs += inSchema[i].nAttrs - 2; // flag column, join column
        else resNAttrs += inSchema[i].nAttrs - 1; // flag column
    }
    resSchema->nAttrs = resNAttrs;
    memcpy(resSchema->attrType, inSchema[0].attrType, sizeof(ATTR_TYPE) * inSchema[0].nAttrs);
    memcpy(resSchema->attrSize, inSchema[0].attrSize, sizeof(uint32_t) * inSchema[0].nAttrs);
    memcpy(resSchema->attrOffset, inSchema[0].attrOffset, sizeof(uint32_t) * inSchema[0].nAttrs);
    uint32_t index = inSchema[0].nAttrs;
    uint32_t resTotSize = inSchema[0].item_size;
    for (uint32_t i = 1; i < nInputs; ++i) {
        for (uint32_t j = 1; j < inSchema[i].nAttrs; ++j) {
            if (equi && j == attrID[i]) continue;
            resSchema->attrType[index] = inSchema[i].attrType[j];
            resSchema->attrSize[index] = inSchema[i].attrSize[j];
            resSchema->attrOffset[index] = resTotSize;
            resTotSize += resSchema->attrSize[index];
            ++index;
        }
    }
    resSchema->item_size = resTotSize;
    resSchema->item_per_blk = (plain_len - META_BLOCK_SIZE) / resTotSize;
}

// for queries without projection (generating intermediate result)
void generateOutputSchema(const Schema** inSchema, const uint32_t nJoinTables, const int32_t* joinAttrID, Schema* resSchema, bool equi = true) {
    uint32_t resNAttrs = inSchema[0]->nAttrs;
    for (uint32_t i = 1; i < nJoinTables; ++i) {
        if (equi) resNAttrs += inSchema[i]->nAttrs - 2; // flag column, join column
        else resNAttrs += inSchema[i]->nAttrs - 1; // flag column
    }
    resSchema->nAttrs = resNAttrs;
    memcpy(resSchema->attrType, inSchema[0]->attrType, sizeof(ATTR_TYPE) * inSchema[0]->nAttrs);
    memcpy(resSchema->attrName, inSchema[0]->attrName, sizeof(uint32_t) * inSchema[0]->nAttrs);
    memcpy(resSchema->attrSize, inSchema[0]->attrSize, sizeof(uint32_t) * inSchema[0]->nAttrs);
    memcpy(resSchema->attrOffset, inSchema[0]->attrOffset, sizeof(uint32_t) * inSchema[0]->nAttrs);
    uint32_t index = inSchema[0]->nAttrs;
    uint32_t resTotSize = inSchema[0]->item_size;
    for (uint32_t i = 1; i < nJoinTables; ++i) {
        for (uint32_t j = 1; j < inSchema[i]->nAttrs; ++j) {
            if (equi && j == joinAttrID[i]) continue;
            resSchema->attrType[index] = inSchema[i]->attrType[j];
            resSchema->attrName[index] = inSchema[i]->attrName[j];
            resSchema->attrSize[index] = inSchema[i]->attrSize[j];
            resSchema->attrOffset[index] = resTotSize;
            resTotSize += resSchema->attrSize[index];
            ++index;
        }
    }
    resSchema->item_size = resTotSize;
    resSchema->item_per_blk = (plain_len - META_BLOCK_SIZE) / resTotSize;
}

// for queries with projection (generating intermediate result)
void generateOutputSchema(const Schema** inSchema, const uint32_t nJoinTables, const uint32_t nTables, const uint32_t* iTable, const int32_t* attrID, const int32_t* iParent, const int32_t* pAttrID, const uint32_t nProjCols, uint32_t* projID[], Schema* resSchema) {
    std::unordered_set <uint32_t> candAttrName;
    for (uint32_t index = 0; index < nProjCols; ++index) {
        uint32_t iTableIndex = projID[0][index];
        uint32_t iAttrIndex = projID[1][index];
        uint32_t tableID = iTable[iTableIndex];
        uint32_t attrName = (tableID << 16) + iAttrIndex;
        candAttrName.insert(attrName);
    }
    for (uint32_t index = 1; index < nTables; ++index) {
        uint32_t tableID = iTable[index];
        uint32_t iAttrIndex = attrID[index];
        uint32_t attrName = (tableID << 16) + iAttrIndex;
        candAttrName.insert(attrName);
    }
    for (uint32_t index = 1; index < nTables; ++index) {
        uint32_t tableID = iParent[index];
        uint32_t iAttrIndex = pAttrID[index];
        uint32_t attrName = (tableID << 16) + iAttrIndex;
        candAttrName.insert(attrName);
    }
    
    resSchema->attrType[0] = CHAR;
    resSchema->attrName[0] = 0;
    resSchema->attrSize[0] = 1;
    resSchema->attrOffset[0] = 0;
    uint32_t resNAttrs = 1;
    uint32_t resTotSize = resSchema->attrSize[0];
    for (uint32_t i = 0; i < nJoinTables; ++i) {
        for (uint32_t j = 1; j < inSchema[i]->nAttrs; ++j) {
            uint32_t attrName = inSchema[i]->attrName[j];
            auto got = candAttrName.find(attrName);
            if (got != candAttrName.end()) {
                resSchema->attrType[resNAttrs] = inSchema[i]->attrType[j];
                resSchema->attrName[resNAttrs] = inSchema[i]->attrName[j];
                resSchema->attrSize[resNAttrs] = inSchema[i]->attrSize[j];
                resSchema->attrOffset[resNAttrs] = resTotSize;
                resTotSize += resSchema->attrSize[resNAttrs];
                ++resNAttrs;
            }
        }
    }
    resSchema->nAttrs = resNAttrs;
    resSchema->item_size = resTotSize;
    resSchema->item_per_blk = (plain_len - META_BLOCK_SIZE) / resTotSize;
}

// for queries without projection (generating final result)
void generateOutputSchema(const Schema* inSchema, const uint32_t nTables, const uint32_t* iTable, const int32_t* attrID, Schema* resSchema, bool equi = true) {
    uint32_t tableID = iTable[0];
    uint32_t resNAttrs = inSchema[tableID].nAttrs;
    for (uint32_t i = 1; i < nTables; ++i) {
        tableID = iTable[i];
        if (equi) resNAttrs += inSchema[tableID].nAttrs - 2; // flag column, join column
        else resNAttrs += inSchema[tableID].nAttrs - 1; // flag column
    }
    resSchema->nAttrs = resNAttrs;
    
    tableID = iTable[0];
    memcpy(resSchema->attrType, inSchema[tableID].attrType, sizeof(ATTR_TYPE) * inSchema[tableID].nAttrs);
    memcpy(resSchema->attrName, inSchema[tableID].attrName, sizeof(uint32_t) * inSchema[tableID].nAttrs);
    memcpy(resSchema->attrSize, inSchema[tableID].attrSize, sizeof(uint32_t) * inSchema[tableID].nAttrs);
    memcpy(resSchema->attrOffset, inSchema[tableID].attrOffset, sizeof(uint32_t) * inSchema[tableID].nAttrs);
    uint32_t index = inSchema[tableID].nAttrs;
    uint32_t resTotSize = inSchema[tableID].item_size;
    for (uint32_t i = 1; i < nTables; ++i) {
        tableID = iTable[i];
        for (uint32_t j = 1; j < inSchema[tableID].nAttrs; ++j) {
            if (equi && j == attrID[i]) continue;
            resSchema->attrType[index] = inSchema[tableID].attrType[j];
            resSchema->attrName[index] = inSchema[tableID].attrName[j];
            resSchema->attrSize[index] = inSchema[tableID].attrSize[j];
            resSchema->attrOffset[index] = resTotSize;
            resTotSize += resSchema->attrSize[index];
            ++index;
        }
    }
    resSchema->item_size = resTotSize;
    resSchema->item_per_blk = (plain_len - META_BLOCK_SIZE) / resTotSize;
}

// for queries with projection (generating final result)
void generateOutputSchema(const Schema* inSchema, const uint32_t nTables, const uint32_t* iTable, const uint32_t nProjCols, uint32_t* projID[], Schema* resSchema) {
    resSchema->nAttrs = nProjCols + 1;
    resSchema->attrType[0] = CHAR;
    resSchema->attrName[0] = 0;
    resSchema->attrSize[0] = 1;
    resSchema->attrOffset[0] = 0;
    uint32_t resTotSize = resSchema->attrSize[0];
    for (uint32_t index = 1; index <= nProjCols; ++index) {
        uint32_t iTableIndex = projID[0][index - 1];
        uint32_t iAttrIndex = projID[1][index - 1];
        int32_t tableID = iTable[iTableIndex];
        resSchema->attrType[index] = inSchema[tableID].attrType[iAttrIndex];
        resSchema->attrName[index] = inSchema[tableID].attrName[iAttrIndex];
        resSchema->attrSize[index] = inSchema[tableID].attrSize[iAttrIndex];
        resSchema->attrOffset[index] = resTotSize;
        resTotSize += resSchema->attrSize[index];
    }
    resSchema->item_size = resTotSize;
    resSchema->item_per_blk = (plain_len - META_BLOCK_SIZE) / resTotSize;
}

void filterOutput(ServerConnector* out_conn, uint32_t & out_tot_block, const Schema* out_sch, const CMP* out_cmp, uint32_t & access_num) {
    // initialize the buffer
    initBuffer();

    // SIGMOD 2022: obliviously sort the output table
    // ObliviousSort* oblisort = new ObliviousSort(out_conn, out_tot_block, out_sch, out_cmp);
    // access_num = oblisort->BatcherSort();
    // delete oblisort;
    
    // Journal Version: obliviously compact the output table
    ObliviousCompact* oblicompact = new ObliviousCompact(out_conn, out_tot_block, out_sch);
    access_num = oblicompact->Compact();
    delete oblicompact;
    
    // obliviously filter the output table
    /****************/
    uint32_t test_real_item_num = 0;
    /****************/
    uint32_t lBlock = 0;
    uint32_t rBlock;
    uint32_t remain_block = out_tot_block;
    uint32_t block_cnt = 0;
    while (remain_block > 0) {
        uint32_t pickup_block;
        if (remain_block >= two_m_block)
            pickup_block = two_m_block;
        else pickup_block = remain_block;
        remain_block -= pickup_block;
        rBlock = lBlock + pickup_block;
        readBlock(out_conn, lBlock, rBlock, buffer);
        
        char* pBlock = buffer + META_BLOCK_SIZE;
        bool find = false;
        uint32_t indexI, indexJ;
        for (indexI = 0; indexI < pickup_block; ++indexI) {
            char* pItem = pBlock;
            for (indexJ = 0; indexJ < out_sch->item_per_blk; ++indexJ) {
                if (*pItem == 'd') {
                    find = true;
                    break;
                }
                assert(*pItem == 'r');
                /****************/
                ++test_real_item_num;
                /****************/
                if (indexJ < out_sch->item_per_blk - 1)
                    pItem += out_sch->item_size;
            }
            if (find) break;
            if (indexI < pickup_block - 1)
                pBlock += B;
        }
        if (find) {
            if (indexJ == 0) block_cnt += indexI;
            else block_cnt += indexI + 1;
            if (block_cnt < out_tot_block)
                resizeFile(out_conn, block_cnt);
            break;
        }
        block_cnt += pickup_block;
        lBlock = rBlock;
    }
    out_tot_block = block_cnt;
    access_num += block_cnt * 2;
    /****************/
    printf("\nThe number of real items in result table: %u\n\n", test_real_item_num);
    /****************/
    
    // destroy the buffer
    destroyBuffer();
}

void printResult(ServerConnector* conn, const uint32_t nBlocks, const Schema* resSchema, const uint32_t attrLen, const int32_t* attrID) {
    // initialize the buffer
    initBuffer();
    
    uint32_t test_real_item_num = 0;
    uint32_t lBlock = 0;
    uint32_t rBlock;
    uint32_t remain_block = nBlocks;
    uint32_t block_cnt = 0;
    while (remain_block > 0) {
        uint32_t pickup_block;
        if (remain_block >= two_m_block)
            pickup_block = two_m_block;
        else pickup_block = remain_block;
        remain_block -= pickup_block;
        rBlock = lBlock + pickup_block;
        readBlock(conn, lBlock, rBlock, buffer);
        
        char* pBlock = buffer + META_BLOCK_SIZE;
        uint32_t indexI, indexJ, indexL;
        for (indexI = 0; indexI < pickup_block; ++indexI) {
            char* pItem = pBlock;
            for (indexJ = 0; indexJ < resSchema->item_per_blk; ++indexJ) {
                if (*pItem == 'r') {
                    /****************
                    for (indexL = 0; indexL < attrLen; ++indexL) {
                        uint32_t curAttrID = attrID[indexL];
                        uint32_t offset = resSchema->attrOffset[curAttrID];
                        ATTR_TYPE attrType = resSchema->attrType[curAttrID];
                        uint32_t attrSize = resSchema->attrSize[curAttrID];
                        
                        char* pAttr = pItem + offset;
                        if (attrType == CHAR) {
                            char attrValue = *pAttr;
                            printf("%c ", attrValue);
                        }
                        else if (attrType == INTEGER) {
                            int32_t attrValue;
                            memcpy(&attrValue, pAttr, sizeof(int32_t));
                            printf("%d ", attrValue);
                        }
                        else if (attrType == DOUBLE) {
                            double attrValue;
                            memcpy(&attrValue, pAttr, sizeof(double));
                            printf("%lf ", attrValue);
                        }
                        else if (attrType == STRING || attrType == TINYTEXT) {
                            char attrValue[attrSize + 1];
                            memcpy(attrValue, pAttr, attrSize);
                            attrValue[attrSize] = '\0';
                            printf("%s ", attrValue);
                        }
                    }
                    printf("\n");
                    /****************/
                    ++test_real_item_num;
                }
                if (indexJ < resSchema->item_per_blk - 1)
                    pItem += resSchema->item_size;
            }
            if (indexI < pickup_block - 1)
                pBlock += B;
        }
        block_cnt += pickup_block;
        lBlock = rBlock;
    }
    printf("\nThe number of real items in result table: %u\n\n", test_real_item_num);
    
    // destroy the buffer
    destroyBuffer();
}

#endif //__APP_H__
