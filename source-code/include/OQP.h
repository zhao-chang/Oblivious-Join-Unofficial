#ifndef __OQP_H__
#define __OQP_H__

#include "Schema.h"
#include "Util.h"
#include "App.h"
#include "Basic.h"
#include "ObliviousSort.h"
#include <algorithm>
#include <memory>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <queue>

using namespace std;

// ICDT 2014: Oblivious Query Processing
class OQP {
public:
    OQP(const uint32_t mode, const uint32_t n_inputs, const std::string* input_file, const std::string& dst, const std::string& prefix) : nInputs(n_inputs) {
        in_sch = new Schema[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i)
            generateInputSchema(mode, i, &(in_sch[i]));
        
        dst_prefix = dst;
        std::string key_name = prefix + std::string("_key.txt");
        initKey(key_name.c_str(), true);
        
        in_conn = new ServerConnector*[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i) {
            std::string cur_dst = dst + "_input_" + std::to_string(i);
            in_conn[i] = new FileSimulator(server_host, cur_dst, true);
        }
        
        nBlocks = new uint32_t[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i) {
            std::unordered_map <uint32_t, std::string> blocks;
            importData(input_file[i], &(in_sch[i]), blocks);
            nBlocks[i] = blocks.size();
            std::vector <std::pair<uint32_t, std::string> > insert_buffer;
            for (auto j = blocks.begin(); j != blocks.end(); ++j) {
                std::string b_id = std::string((const char *)(&(j->first)), sizeof(uint32_t));
                std::string cipher;
                enc_engine.encrypt(b_id + j->second, enc_key, cipher);
                insert_buffer.emplace_back(std::make_pair(j->first, cipher));
                if (insert_buffer.size() >= 10000) {
                    in_conn[i]->insert(insert_buffer);
                    insert_buffer.clear();
                }
            }
            in_conn[i]->insert(insert_buffer);
            insert_buffer.clear();
        }
    }
    
    ~OQP() {
        if (in_sch != NULL) {
            delete[] in_sch;
            in_sch = NULL;
        }
        if (in_conn != NULL) {
            for (uint32_t i = 0; i < nInputs; ++i)
                delete in_conn[i];
            delete[] in_conn;
            in_conn = NULL;
        }
        if (nBlocks != NULL) {
            delete[] nBlocks;
            nBlocks = NULL;
        }
    }
    
    void ObliJoin(const uint32_t* attr_id) {
        attrID = new int32_t[nInputs];
        memcpy(attrID, attr_id, nInputs * sizeof(uint32_t));
        std::string out_dst = dst_prefix + "_output";
        for (uint32_t i = 0; i < nInputs; ++i) {
            out_dst += "_" + std::to_string(attrID[i]);
        }
        out_conn = new FileSimulator(server_host, out_dst.c_str(), true);
        out_sch = new Schema;
        generateOutputSchema(in_sch, nInputs, attrID, out_sch);
        out_cmp = new CMP({1, {attrID[0]}, {0}});
        
        // initialize the buffer
        initBuffer();
        
        // create input tables with join degrees
        ServerConnector** weight_conn = new ServerConnector*[nInputs];
        Schema* weight_sch = new Schema[nInputs];
        uint32_t* weight_n_blocks = new uint32_t[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i) {
            std::string weight_dst = dst_prefix + "_weight_" + std::to_string(i);
            weight_conn[i] = new FileSimulator(server_host, weight_dst.c_str(), true);
        }
        
        // compute the join degrees of input tables
        computeDegree(weight_conn, weight_sch, weight_n_blocks);
        
        // create expanded input tables
        ServerConnector** expand_conn = new ServerConnector*[nInputs];
        Schema* expand_sch = new Schema[nInputs];
        uint32_t* expand_n_blocks = new uint32_t[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i) {
            std::string expand_dst = dst_prefix + "_expand_" + std::to_string(i);
            expand_conn[i] = new FileSimulator(server_host, expand_dst.c_str(), true);
        }
        
        // expand the input tables based on the degrees
        generateJId(weight_conn[1], &(weight_sch[1]), weight_n_blocks[1], attrID[1]);
        expandTable(weight_conn, weight_sch, weight_n_blocks, expand_conn, expand_sch, expand_n_blocks, 0);
        generateJId(expand_conn[0], &(expand_sch[0]), expand_n_blocks[0], expand_sch[0].nAttrs - 1);
        expandTable(weight_conn, weight_sch, weight_n_blocks, expand_conn, expand_sch, expand_n_blocks, 1);
        
        for (uint32_t i = 0; i < nInputs; ++i) {
            delete weight_conn[i];
        }
        delete[] weight_conn;
        delete[] weight_sch;
        delete[] weight_n_blocks;
        
        // stitch expanded input tables
        stitchTable(expand_conn, expand_sch, expand_n_blocks);
        
        for (uint32_t i = 0; i < nInputs; ++i) {
            delete expand_conn[i];
        }
        delete[] expand_conn;
        delete[] expand_sch;
        delete[] expand_n_blocks;
        
        // destroy the buffer
        destroyBuffer();
        
        if (attrID != NULL) {
            delete[] attrID;
            attrID = NULL;
        }
        if (out_conn != NULL) {
            delete out_conn;
            out_conn = NULL;
        }
        if (out_sch != NULL) {
            delete out_sch;
            out_sch = NULL;
        }
        if (out_cmp != NULL) {
            delete out_cmp;
            out_cmp = NULL;
        }
    }
    
private:
    void computeDegree(ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks) {
        // prepare tables for join degree computing
        ServerConnector** prep_conn = new ServerConnector*[nInputs]; // prepared input tables
        Schema* prep_sch = new Schema[nInputs]; // schema of prepared input tables
        uint32_t* prep_n_blocks = new uint32_t[nInputs]; // # of blocks in prepared input tables
        for (uint32_t i = 0; i < nInputs; ++i) {
            std::string prep_dst = dst_prefix + "_prep_" + std::to_string(i);
            prep_conn[i] = new FileSimulator(server_host, prep_dst.c_str(), true);
        }
        
        for (uint32_t i = 0; i < nInputs; ++i)
            prepTable(prep_conn, prep_sch, prep_n_blocks, i);
        
        // unify prepared tables
        std::string union_dst = dst_prefix + "_union";
        ServerConnector* union_conn = new FileSimulator(server_host, union_dst.c_str(), true);
        Schema* union_sch = new Schema; // schema of unified input table
        uint32_t union_n_blocks; // # of blocks in unified input table
        
        unifyTable(prep_conn, prep_sch, prep_n_blocks, union_conn, union_sch, union_n_blocks);
        
        for (uint32_t i = 0; i < nInputs; ++i) {
            delete prep_conn[i];
        }
        delete[] prep_conn;
        delete[] prep_sch;
        delete[] prep_n_blocks;
        
        // compute group running sums
        for (uint32_t i = 0; i < nInputs; ++i)
            rsumTable(union_conn, union_sch, union_n_blocks, i);
        
        // generate the input tables with join degrees
        for (uint32_t i = 0; i < nInputs; ++i)
            generateWTable(union_conn, union_sch, union_n_blocks, weight_conn, weight_sch, weight_n_blocks, i);
        
        delete union_conn;
        delete union_sch;
        
        // filter out dummy items
        for (uint32_t i = 0; i < nInputs; ++i)
            filterWTable(weight_conn, weight_sch, weight_n_blocks, i);
    }
    
    void prepTable(ServerConnector** prep_conn, Schema* prep_sch, uint32_t* prep_n_blocks, const uint32_t tableID) {
        // [0]: in_sch[tableID]; [1]: prep_sch[tableID];
        const uint32_t append_num = 4;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = in_sch[tableID].nAttrs;
        nAttrs[1] = nAttrs[0] + append_num;
        itemSize[0] = in_sch[tableID].item_size;
        itemSize[1] = itemSize[0] + append_num * sizeof(uint32_t); // N <- 1, N_R/N_S, Src <- tableID, Id <- ID
        itemPerBlk[0] = in_sch[tableID].item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        // initialize the schema of the prepared table
        Schema* prep_sch_tid = &(prep_sch[tableID]);
        prep_sch_tid->nAttrs = nAttrs[1];
        prep_sch_tid->item_size = itemSize[1];
        prep_sch_tid->item_per_blk = itemPerBlk[1];
        memcpy(prep_sch_tid->attrType, in_sch[tableID].attrType, sizeof(ATTR_TYPE) * nAttrs[0]);
        memcpy(prep_sch_tid->attrSize, in_sch[tableID].attrSize, sizeof(uint32_t) * nAttrs[0]);
        memcpy(prep_sch_tid->attrOffset, in_sch[tableID].attrOffset, sizeof(uint32_t) * nAttrs[0]);
        uint32_t curOffset = itemSize[0];
        for (uint32_t i = nAttrs[0]; i < nAttrs[1]; ++i) {
            prep_sch_tid->attrType[i] = INTEGER;
            prep_sch_tid->attrSize[i] = sizeof(uint32_t);
            prep_sch_tid->attrOffset[i] = curOffset;
            curOffset += sizeof(uint32_t);
        }
        
        // construct the prepared table
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t append_cols[append_num]; // N <- 1, N_R/N_S, Src <- tableID, Id <- ID
        append_cols[0] = 1;
        append_cols[1] = 0;
        append_cols[2] = tableID;
        append_cols[3] = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = nBlocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(in_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        ++append_cols[3];
                        memcpy(oPItem, iPItem, itemSize[0]);
                        memcpy(oPItem + itemSize[0], append_cols, append_num * sizeof(uint32_t));
                        
                        ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(prep_conn[tableID], oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(prep_conn[tableID], oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        prep_n_blocks[tableID] = oBlockID;
    }
    
    void unifyTable(ServerConnector** prep_conn, const Schema* prep_sch, const uint32_t* prep_n_blocks, ServerConnector* union_conn, Schema* union_sch, uint32_t & union_n_blocks) {
        // [0]: prep_sch[0]; [1]: prep_sch[1]; [2]: union_sch[2];
        uint32_t nAttrs[3];
        uint32_t itemSize[3];
        uint32_t itemPerBlk[3];
        for (uint32_t tableID = 0; tableID < 2; ++tableID) {
            nAttrs[tableID] = prep_sch[tableID].nAttrs;
            itemSize[tableID] = prep_sch[tableID].item_size;
            itemPerBlk[tableID] = prep_sch[tableID].item_per_blk;
        }
        
        // initialize the schema of the unified table
        const uint32_t append_num = 4;
        union_sch->nAttrs = nAttrs[0] + nAttrs[1] - append_num - 2; // flag column, joined column
        uint32_t union_index = nAttrs[0] - append_num;
        memcpy(union_sch->attrType, prep_sch[0].attrType, sizeof(ATTR_TYPE) * union_index);
        memcpy(union_sch->attrSize, prep_sch[0].attrSize, sizeof(uint32_t) * union_index);
        memcpy(union_sch->attrOffset, prep_sch[0].attrOffset, sizeof(uint32_t) * union_index);
        uint32_t union_item_size = prep_sch[0].attrOffset[union_index];
        for (uint32_t j = 1; j < nAttrs[1]; ++j) {
            if (j == attrID[1]) continue;
            union_sch->attrType[union_index] = prep_sch[1].attrType[j];
            union_sch->attrSize[union_index] = prep_sch[1].attrSize[j];
            union_sch->attrOffset[union_index] = union_item_size;
            union_item_size += union_sch->attrSize[union_index];
            ++union_index;
        }
        union_sch->item_size = union_item_size;
        union_sch->item_per_blk = (plain_len - META_BLOCK_SIZE) / union_item_size;
        
        nAttrs[2] = union_sch->nAttrs;
        itemSize[2] = union_sch->item_size;
        itemPerBlk[2] = union_sch->item_per_blk;
        
        // construct the unified table
        uint32_t attr_size[nInputs];
        uint32_t offset[nInputs];
        uint32_t append_attr_size[nInputs];
        uint32_t append_offset[nInputs];
        for (uint32_t tableID = 0; tableID < 2; ++tableID) {
            uint32_t curAttrID = attrID[tableID];
            attr_size[tableID] = prep_sch[tableID].attrSize[curAttrID];
            offset[tableID] = prep_sch[tableID].attrOffset[curAttrID];
            curAttrID = nAttrs[tableID] - append_num;
            append_attr_size[tableID] = append_num * sizeof(uint32_t);
            append_offset[tableID] = prep_sch[tableID].attrOffset[curAttrID];
        }
        
        uint32_t iItemSegPos[2];
        iItemSegPos[0] = itemSize[0] - append_attr_size[0];
        iItemSegPos[1] = offset[1] + attr_size[1];
        uint32_t iZeroSegSize[2];
        iZeroSegSize[0] = itemSize[2] - itemSize[0];
        iZeroSegSize[1] = iItemSegPos[0] - sizeof(char);
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = offset[1] - sizeof(char);
        iItemSegSize[1] = itemSize[1] - iItemSegPos[1];
        uint32_t oItemSegPos[2];
        oItemSegPos[0] = itemSize[2] - append_attr_size[0];
        oItemSegPos[1] = iItemSegPos[0] + iItemSegSize[0];
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        for (uint32_t tableID = 0; tableID < 2; ++tableID) {
            uint32_t iLeft = 0;
            uint32_t iRight;
            uint32_t remain_block = prep_n_blocks[tableID];
            while (remain_block > 0) {
                uint32_t pickup_block;
                if (remain_block >= two_m_block)
                    pickup_block = two_m_block;
                else pickup_block = remain_block;
                remain_block -= pickup_block;
                iRight = iLeft + pickup_block;
                readBlock(prep_conn[tableID], iLeft, iRight, buffer);
                
                char* iPBlock = buffer + META_BLOCK_SIZE;
                for (uint32_t i = 0; i < pickup_block; ++i) {
                    char* iPItem = iPBlock;
                    for (uint32_t j = 0; j < itemPerBlk[tableID]; ++j) {
                        if (*iPItem == 'r') {
                            if (tableID == 0) {
                                memcpy(oPItem, iPItem, iItemSegPos[0]);
                                memset(oPItem + iItemSegPos[0], 0, iZeroSegSize[0]);
                                memcpy(oPItem + oItemSegPos[0], iPItem + iItemSegPos[0], append_attr_size[0]);
                            }
                            else if (tableID == 1) {
                                memcpy(oPItem, iPItem, sizeof(char));
                                memset(oPItem + sizeof(char), 0, iZeroSegSize[1]);
                                memcpy(oPItem + offset[0], iPItem + offset[1], attr_size[1]);
                                memcpy(oPItem + iItemSegPos[0], iPItem + sizeof(char), iItemSegSize[0]);
                                memcpy(oPItem + oItemSegPos[1], iPItem + iItemSegPos[1], iItemSegSize[1]);
                            }
                            ++oRItemCnt;
                            ++oItemCnt;
                            if (oItemCnt < itemPerBlk[2])
                                oPItem += itemSize[2];
                            else {
                                memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                                memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                                insertBlock(union_conn, oBlockID, oBlock);
                                ++oBlockID;
                                oItemCnt = 0;
                                oRItemCnt = 0;
                                oPItem = oBlock + META_BLOCK_SIZE;
                            }
                        }
                        if (j < itemPerBlk[tableID] - 1)
                            iPItem += itemSize[tableID];
                    }
                    if (i < pickup_block - 1)
                        iPBlock += B;
                }
                iLeft = iRight;
            }
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[2] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[2]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[2]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[2])
                    oPItem += itemSize[2];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(union_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        union_n_blocks = oBlockID;
    }
    
    void rsumTable(ServerConnector* union_conn, const Schema* union_sch, const uint32_t union_n_blocks, const uint32_t src) {
        // obliviously sort the unified table
        assert(src == 0 || src == 1);
        uint8_t order = 1 - src;
        uint32_t nAttrs = union_sch->nAttrs;
        CMP union_cmp {3, {attrID[0], nAttrs - 2, nAttrs - 1}, {0, order, 0}};
        
        ObliviousSort* oblisort = new ObliviousSort(union_conn, union_n_blocks, union_sch, &union_cmp);
        oblisort->BatcherSort();
        delete oblisort;
        
        uint32_t itemSize = union_sch->item_size;
        uint32_t itemPerBlk = union_sch->item_per_blk;
        uint32_t attrType = union_sch->attrType[attrID[0]];
        uint32_t offset = union_sch->attrOffset[attrID[0]];
        uint32_t attrSize = union_sch->attrSize[attrID[0]];
        
        // [0]: N; [1]: N_R/N_S; [2]: Src
        const uint32_t append_num = 4;
        uint32_t append_offset[3];
        for (uint32_t i = 0, j = nAttrs - append_num; i < 3; ++i, ++j)
            append_offset[i] = union_sch->attrOffset[j];
        
        // compute the group running sums
        int32_t pre_nrs_val = -1;
        uint32_t pre_attr_len = 0;
        char pre_attr[attrSize];
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = union_n_blocks;
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(union_conn, iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk; ++j) {
                    if (*iPItem == 'r') {
                        uint32_t cur_tid;
                        memcpy(&cur_tid, iPItem + append_offset[2], sizeof(uint32_t));
                        uint32_t cur_n_val;
                        memcpy(&cur_n_val, iPItem + append_offset[0], sizeof(uint32_t));
                        assert(cur_tid == src || cur_tid == 1 - src);
                        assert(cur_n_val == 1);
                        
                        char* cur_attr = iPItem + offset;
                        uint32_t cur_attr_len = 0;
                        if (attrType == TINYTEXT)
                            cur_attr_len = std::min((uint32_t)strlen(cur_attr), attrSize);
                        else cur_attr_len = attrSize;
                        
                        bool equal;
                        if (pre_nrs_val == -1) // first item in the sorted unified table
                            equal = false;
                        else if (pre_attr_len == cur_attr_len && memcmp(pre_attr, cur_attr, cur_attr_len) == 0)
                            equal = true;
                        else equal = false;
                        
                        uint32_t cur_nrs_val = 0;
                        if (equal) {
                            if (cur_tid == src)
                                cur_nrs_val = pre_nrs_val;
                            else if (cur_tid == 1 - src)
                                cur_nrs_val = pre_nrs_val + cur_n_val;
                        }
                        else {
                            if (cur_tid == src)
                                cur_nrs_val = 0;
                            else if (cur_tid == 1 - src)
                                cur_nrs_val = cur_n_val;
                        }
                        if (cur_tid == src)
                            memcpy(iPItem + append_offset[1], &cur_nrs_val, sizeof(uint32_t));
                        
                        pre_nrs_val = cur_nrs_val;
                        pre_attr_len = cur_attr_len;
                        memcpy(pre_attr, cur_attr, cur_attr_len);
                    }
                    if (j < itemPerBlk - 1)
                        iPItem += itemSize;
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            updateBlock(union_conn, iLeft, iRight, buffer);
            iLeft = iRight;
        }
    }
    
    void generateWTable(ServerConnector* union_conn, const Schema* union_sch, const uint32_t union_n_blocks, ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks, const uint32_t tableID) {
        // initialize the schema of the table with join degrees
        uint32_t weight_index = in_sch[tableID].nAttrs;
        weight_sch[tableID].nAttrs = weight_index + 2;
        memcpy(weight_sch[tableID].attrType, in_sch[tableID].attrType, sizeof(ATTR_TYPE) * weight_index);
        memcpy(weight_sch[tableID].attrSize, in_sch[tableID].attrSize, sizeof(uint32_t) * weight_index);
        memcpy(weight_sch[tableID].attrOffset, in_sch[tableID].attrOffset, sizeof(uint32_t) * weight_index);
        uint32_t weight_item_size = in_sch[tableID].item_size;
        for (uint32_t i = 0; i < 2; ++i) {
            weight_sch[tableID].attrType[weight_index] = INTEGER;
            weight_sch[tableID].attrSize[weight_index] = sizeof(uint32_t);
            weight_sch[tableID].attrOffset[weight_index] = weight_item_size;
            weight_item_size += sizeof(uint32_t);
            ++weight_index;
        }
        weight_sch[tableID].item_size = weight_item_size;
        weight_sch[tableID].item_per_blk = (plain_len - META_BLOCK_SIZE) / weight_item_size;
        
        // [0]: union_sch; [1]: weight_sch[tableID];
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = union_sch->nAttrs;
        nAttrs[1] = weight_sch[tableID].nAttrs;
        itemSize[0] = union_sch->item_size;
        itemSize[1] = weight_sch[tableID].item_size;
        itemPerBlk[0] = union_sch->item_per_blk;
        itemPerBlk[1] = weight_sch[tableID].item_per_blk;
        
        // [0]: in_sch[0]; [1]: in_sch[1];
        uint32_t inItemSize[2];
        uint32_t offset[2];
        uint32_t attrSize[2];
        for (uint32_t i = 0; i < 2; ++i) {
            inItemSize[i] = in_sch[i].item_size;
            offset[i] = in_sch[i].attrOffset[attrID[i]];
            attrSize[i] = in_sch[i].attrSize[attrID[i]];
        }
        
        // [0]: N_R/N_S; [1]: Src; [2]: Id;
        uint32_t append_offset[3];
        for (uint32_t i = 0, j = nAttrs[0] - 3; i < 3; ++i, ++j)
            append_offset[i] = union_sch->attrOffset[j];
        
        uint32_t iItemSegSize[2];
        uint32_t oItemSegPos[2];
        iItemSegSize[0] = offset[1] - sizeof(char);
        uint32_t iItemSegPos = inItemSize[0] + iItemSegSize[0];
        oItemSegPos[0] = sizeof(char) + iItemSegSize[0];
        oItemSegPos[1] = oItemSegPos[0] + attrSize[0];
        iItemSegSize[1] = inItemSize[1] - oItemSegPos[1];
        uint32_t oItemLastPos[2];
        oItemLastPos[0] = inItemSize[0] + sizeof(uint32_t);
        oItemLastPos[1] = inItemSize[1] + sizeof(uint32_t);
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        char oItem[plain_len];
        char* oSItem = oItem + sizeof(char);
        uint32_t oItemContentSize = itemSize[1] - sizeof(char);
        
        // generate the input tables with join degrees
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = union_n_blocks;
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(union_conn, iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        uint32_t tid;
                        memcpy(&tid, iPItem + append_offset[1], sizeof(uint32_t));
                        assert(tid == 0 || tid == 1);
                        if (tid == tableID) {
                            if (tableID == 0) {
                                memcpy(oItem, iPItem, inItemSize[0]);
                                memcpy(oItem + inItemSize[0], iPItem + append_offset[0], sizeof(uint32_t));
                                memcpy(oItem + oItemLastPos[0], iPItem + append_offset[2], sizeof(uint32_t));
                            }
                            else if (tableID == 1) {
                                memcpy(oItem, iPItem, sizeof(char));
                                memcpy(oItem + sizeof(char), iPItem + inItemSize[0], iItemSegSize[0]);
                                memcpy(oItem + oItemSegPos[0], iPItem + offset[0], attrSize[0]);
                                memcpy(oItem + oItemSegPos[1], iPItem + iItemSegPos, iItemSegSize[1]);
                                memcpy(oItem + inItemSize[1], iPItem + append_offset[0], sizeof(uint32_t));
                                memcpy(oItem + oItemLastPos[1], iPItem + append_offset[2], sizeof(uint32_t));
                            }
                            ++oRItemCnt;
                        }
                        else {
                            oItem[0] = 'd';
                            std::string rnd_str = generate_random_block(oItemContentSize);
                            memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                        }
                        memcpy(oPItem, oItem, itemSize[1]);
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(weight_conn[tableID], oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(weight_conn[tableID], oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        weight_n_blocks[tableID] = oBlockID;
    }
    
    void filterWTable(ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks, const uint32_t tableID) {
        // obliviously sort the table with join degrees
        uint32_t nAttrs = weight_sch[tableID].nAttrs;
        CMP weight_cmp {2, {attrID[tableID], nAttrs - 1}, {0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(weight_conn[tableID], weight_n_blocks[tableID], &(weight_sch[tableID]), &weight_cmp);
        oblisort->BatcherSort();
        delete oblisort;
        
        uint32_t itemSize = weight_sch[tableID].item_size;
        uint32_t itemPerBlk = weight_sch[tableID].item_per_blk;
        
        // obliviously filter the table with join degrees
        uint32_t lBlock = 0;
        uint32_t rBlock;
        uint32_t remain_block = weight_n_blocks[tableID];
        uint32_t block_cnt = 0;
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            rBlock = lBlock + pickup_block;
            readBlock(weight_conn[tableID], lBlock, rBlock, buffer);
            
            char* pBlock = buffer + META_BLOCK_SIZE;
            bool find = false;
            uint32_t indexI, indexJ;
            for (indexI = 0; indexI < pickup_block; ++indexI) {
                char* pItem = pBlock;
                for (indexJ = 0; indexJ < itemPerBlk; ++indexJ) {
                    if (*pItem == 'd') {
                        find = true;
                        break;
                    }
                    if (indexJ < itemPerBlk - 1)
                        pItem += itemSize;
                }
                if (find) break;
                if (indexI < pickup_block - 1)
                    pBlock += B;
            }
            if (find) {
                if (indexJ == 0) block_cnt += indexI;
                else block_cnt += indexI + 1;
                resizeFile(weight_conn[tableID], block_cnt);
                break;
            }
            block_cnt += pickup_block;
            lBlock = rBlock;
        }
        weight_n_blocks[tableID] = block_cnt;
    }
    
    // generate JId using S.(JId <- ID_J) or R_{exp}.(JId <- ID_{Id})
    void generateJId(ServerConnector* & jid_conn, Schema* jid_sch, uint32_t & jid_n_blocks, const uint32_t sort_col) {
        // [0]: old schema; [1]: new schema;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = jid_sch->nAttrs;
        nAttrs[1] = jid_sch->nAttrs + 1;
        itemSize[0] = jid_sch->item_size;
        itemSize[1] = jid_sch->item_size + sizeof(uint32_t);
        itemPerBlk[0] = jid_sch->item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        uint32_t attrType = jid_sch->attrType[sort_col];
        uint32_t offset = jid_sch->attrOffset[sort_col];
        uint32_t attrSize = jid_sch->attrSize[sort_col];
        
        // obliviously sort the old table
        CMP jid_cmp {2, {sort_col, nAttrs[0] - 1}, {0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(jid_conn, jid_n_blocks, jid_sch, &jid_cmp);
        oblisort->BatcherSort();
        delete oblisort;
        
        // unify prepared tables
        std::string tmp_dst = jid_conn->getCollectionName() + "_tmp_jid";
        ServerConnector* tmp_conn = new FileSimulator(server_host, tmp_dst.c_str(), true);
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        // assign the JId values
        int32_t pre_jid = -1;
        uint32_t pre_attr_len = 0;
        char pre_attr[attrSize];
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = jid_n_blocks;
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(jid_conn, iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        char* cur_attr = iPItem + offset;
                        uint32_t cur_attr_len;
                        if (attrType == TINYTEXT)
                            cur_attr_len = std::min((uint32_t)strlen(cur_attr), attrSize);
                        else cur_attr_len = attrSize;
                        
                        bool equal;
                        if (pre_jid == -1) // first item in the sorted unified table
                            equal = false;
                        else if (pre_attr_len == cur_attr_len && memcmp(pre_attr, cur_attr, cur_attr_len) == 0)
                            equal = true;
                        else equal = false;
                        
                        uint32_t cur_jid;
                        if (equal) cur_jid = pre_jid + 1;
                        else cur_jid = 1;
                        
                        pre_jid = cur_jid;
                        pre_attr_len = cur_attr_len;
                        memcpy(pre_attr, cur_attr, cur_attr_len);
                        
                        memcpy(oPItem, iPItem, itemSize[0]);
                        memcpy(oPItem + itemSize[0], &cur_jid, sizeof(uint32_t));
                        
                        ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(tmp_conn, oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(tmp_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        jid_n_blocks = oBlockID;
        
        // update the schema of the table with JId column
        jid_sch->nAttrs = nAttrs[1];
        jid_sch->attrType[nAttrs[0]] = INTEGER;
        jid_sch->attrSize[nAttrs[0]] = sizeof(uint32_t);
        jid_sch->attrOffset[nAttrs[0]] = itemSize[0];
        jid_sch->item_size = itemSize[1];
        jid_sch->item_per_blk = itemPerBlk[1];
        
        // update the connector of the table with JId column
        delete jid_conn;
        jid_conn = tmp_conn;
    }
    
    void expandTable(ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks, ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID) {
        double wDistAvg;
        uint32_t wDistLen = 30;
        uint32_t wDist[wDistLen][2]; // 0, 1, 2, 4, 8, ...
        generateWDist(weight_conn, weight_sch, weight_n_blocks, expand_conn, expand_sch, expand_n_blocks, tableID, wDistLen, wDist);
        generateHId(expand_conn, expand_sch, expand_n_blocks, tableID, wDistLen, wDist, wDistAvg);
        findPos(expand_conn, expand_sch, expand_n_blocks, tableID, wDistLen, wDist, wDistAvg);
        reorderPos(expand_conn, expand_sch, expand_n_blocks, tableID);
        makeCopy(expand_conn, expand_sch, expand_n_blocks, tableID, wDistAvg);
        attachFlag(expand_conn, expand_sch, expand_n_blocks, tableID);
        filterDummy(expand_conn, expand_sch, expand_n_blocks, tableID);
    }
    
    void generateWDist(ServerConnector** weight_conn, const Schema* weight_sch, const uint32_t* weight_n_blocks, ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID, uint32_t & wDistLen, uint32_t wDist[][2]) {
        // [0]: weight_sch[tableID]; [1]: expand_sch[tableID];
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = weight_sch[tableID].nAttrs;
        nAttrs[1] = nAttrs[0] + 1;
        itemSize[0] = weight_sch[tableID].item_size;
        itemSize[1] = itemSize[0] + sizeof(uint32_t);
        itemPerBlk[0] = weight_sch[tableID].item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        // initialize the schema of the table with rounded weights
        expand_sch[tableID].nAttrs = nAttrs[1];
        uint32_t expand_index = 0;
        if (tableID == 0) expand_index = nAttrs[0] - 1;
        else if (tableID == 1) expand_index = nAttrs[0] - 2;
        memcpy(expand_sch[tableID].attrType, weight_sch[tableID].attrType, sizeof(ATTR_TYPE) * expand_index);
        memcpy(expand_sch[tableID].attrSize, weight_sch[tableID].attrSize, sizeof(uint32_t) * expand_index);
        memcpy(expand_sch[tableID].attrOffset, weight_sch[tableID].attrOffset, sizeof(uint32_t) * expand_index);
        uint32_t curOffset = weight_sch[tableID].attrOffset[expand_index];
        expand_sch[tableID].attrType[expand_index] = INTEGER;
        expand_sch[tableID].attrSize[expand_index] = sizeof(uint32_t);
        expand_sch[tableID].attrOffset[expand_index] = curOffset;
        curOffset += sizeof(uint32_t);
        for (uint32_t j = expand_index; j < nAttrs[0]; ++j) {
            uint32_t i = j + 1;
            expand_sch[tableID].attrType[i] = weight_sch[tableID].attrType[j];
            expand_sch[tableID].attrSize[i] = weight_sch[tableID].attrSize[j];
            expand_sch[tableID].attrOffset[i] = curOffset;
            curOffset += expand_sch[tableID].attrSize[i];
        }
        expand_sch[tableID].item_size = itemSize[1];
        expand_sch[tableID].item_per_blk = itemPerBlk[1];
        
        // generate the table with rounded weights
        uint32_t offset = weight_sch[tableID].attrOffset[expand_index - 1];
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = weight_sch[tableID].attrOffset[expand_index];
        iItemSegSize[1] = itemSize[0] - iItemSegSize[0];
        uint32_t oItemSegPos = iItemSegSize[0] + sizeof(uint32_t);
        
        uint32_t sumW = 0;
        uint32_t sumTW = 0;
        memset(wDist[0], 0, sizeof(uint32_t) * wDistLen * 2);
        wDistLen = 0;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        uint32_t oItemID = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = weight_n_blocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(weight_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        char* curAttr = iPItem + offset;
                        uint32_t curW;
                        memcpy(&curW, curAttr, sizeof(uint32_t));
                        uint32_t curTW;
                        if (curW <= 2) {
                            curTW = curW;
                            ++wDist[curTW][1];
                            if (wDistLen < curTW + 1)
                                wDistLen = curTW + 1;
                        }
                        else {
                            uint32_t logTW = (uint32_t)ceil(log2((double)curW));
                            curTW = 1 << logTW;
                            ++wDist[logTW + 1][1];
                            if (wDistLen < logTW + 2)
                                wDistLen = logTW + 2;
                        }
                        sumW += curW;
                        sumTW += curTW;
                        
                        memcpy(oPItem, iPItem, iItemSegSize[0]);
                        memcpy(oPItem + iItemSegSize[0], &curTW, sizeof(uint32_t));
                        memcpy(oPItem + oItemSegPos, iPItem + iItemSegSize[0], iItemSegSize[1]);
                        
                        ++oRItemCnt;
                        ++oItemCnt;
                        ++oItemID;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(expand_conn[tableID], oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        
        // padded with one more dummy item with TW = 2 * sumW - sumTW
        uint32_t wLast = 2 * sumW - sumTW;
        char oItem[plain_len];
        char* oSItem = oItem + sizeof(char);
        oItem[0] = 'r';
        std::string rnd_str = generate_random_block(offset - sizeof(char));
        memcpy(oSItem, rnd_str.c_str(), offset - sizeof(char));
        memset(oItem + offset, 0, sizeof(uint32_t));
        offset += sizeof(uint32_t);
        memcpy(oItem + offset, &wLast, sizeof(uint32_t));
        offset += sizeof(uint32_t);
        ++oItemID;
        memcpy(oItem + offset, &oItemID, sizeof(uint32_t));
        memcpy(oPItem, oItem, itemSize[1]);
        
        ++oRItemCnt;
        ++oItemCnt;
        if (oItemCnt < itemPerBlk[1])
            oPItem += itemSize[1];
        else {
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(expand_conn[tableID], oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
            oPItem = oBlock + META_BLOCK_SIZE;
        }
        
        // update the weight distribution
        for (uint32_t i = 1; i < wDistLen; ++i)
            wDist[i][0] = 1 << (i - 1);
        bool isPowof2 = (wLast == 0) || (wLast & (wLast - 1) == 0);
        bool find = isPowof2 && (wLast < (1 << (wDistLen - 1)));
        if (find) {
            if (wLast <= 2) ++wDist[wLast][1];
            else {
                uint32_t logTW = (uint32_t)log2((double)wLast);
                ++wDist[logTW + 1][1];
            }
        }
        else {
            wDist[wDistLen][0] = wLast;
            wDist[wDistLen][1] = 1;
            ++wDistLen;
        }
        
        // padded with dummy items
        if (oItemCnt > 0) {
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(expand_conn[tableID], oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        expand_n_blocks[tableID] = oBlockID;
    }
    
    void generateHId(ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID, const uint32_t wDistLen, const uint32_t wDist[][2], double & wDistAvg) {
        // compute the value of wDistAvg
        uint32_t wDistSum = 0;
        uint32_t wDistNum = 0;
        for (uint32_t i = 0; i < wDistLen; ++i) {
            wDistNum += wDist[i][1];
            wDistSum += wDist[i][0] * wDist[i][1];
        }
        wDistAvg = (double)wDistSum / wDistNum;
        
        /*****************************************************/
        // generate H column
        
        // [0]: old schema; [1]: new schema;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = expand_sch[tableID].nAttrs;
        nAttrs[1] = nAttrs[0] + 1;
        itemSize[0] = expand_sch[tableID].item_size;
        itemSize[1] = itemSize[0] + sizeof(uint32_t);
        itemPerBlk[0] = expand_sch[tableID].item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        // create a temporary output table
        std::string tmp_dst = dst_prefix + "_tmp_h_" + std::to_string(tableID);
        ServerConnector* tmp_conn = new FileSimulator(server_host, tmp_dst.c_str(), true);
        
        // assign H values
        uint32_t expand_index = 0;
        if (tableID == 0) expand_index = nAttrs[0] - 1;
        else if (tableID == 1) expand_index = nAttrs[0] - 2;
        uint32_t offset = expand_sch[tableID].attrOffset[expand_index - 1];
        
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = expand_sch[tableID].attrOffset[expand_index];
        iItemSegSize[1] = itemSize[0] - iItemSegSize[0];
        uint32_t oItemSegPos = iItemSegSize[0] + sizeof(uint32_t);
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = expand_n_blocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(expand_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        char* curAttr = iPItem + offset;
                        uint32_t curTW;
                        memcpy(&curTW, curAttr, sizeof(uint32_t));
                        uint32_t curH;
                        if ((double)curTW > wDistAvg) curH = 1;
                        else curH = 0;
                        
                        memcpy(oPItem, iPItem, iItemSegSize[0]);
                        memcpy(oPItem + iItemSegSize[0], &curH, sizeof(uint32_t));
                        memcpy(oPItem + oItemSegPos, iPItem + iItemSegSize[0], iItemSegSize[1]);

                        ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(tmp_conn, oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(tmp_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        expand_n_blocks[tableID] = oBlockID;
        
        // update the schema
        expand_sch[tableID].nAttrs = nAttrs[1];
        for (uint32_t i = nAttrs[0]; i > expand_index; --i) {
            uint32_t j = i - 1;
            expand_sch[tableID].attrType[i] = expand_sch[tableID].attrType[j];
            expand_sch[tableID].attrSize[i] = expand_sch[tableID].attrSize[j];
            expand_sch[tableID].attrOffset[i] = expand_sch[tableID].attrOffset[j] + sizeof(uint32_t);
        }
        expand_sch[tableID].attrType[expand_index] = INTEGER;
        expand_sch[tableID].attrSize[expand_index] = sizeof(uint32_t);
        expand_sch[tableID].item_size = itemSize[1];
        expand_sch[tableID].item_per_blk = itemPerBlk[1];
        
        // update the connector
        delete expand_conn[tableID];
        expand_conn[tableID] = tmp_conn;
        
        /*****************************************************/
        // generate HId column
        
        // obliviously sort the old table
        CMP expand_cmp {3, {expand_index, expand_index - 1, expand_index + 1}, {1, 1, 0}};
        ObliviousSort* oblisort = new ObliviousSort(expand_conn[tableID], expand_n_blocks[tableID], &(expand_sch[tableID]), &expand_cmp);
        oblisort->BatcherSort();
        delete oblisort;
        
        // [0]: old schema; [1]: new schema;
        nAttrs[0] = expand_sch[tableID].nAttrs;
        nAttrs[1] = nAttrs[0] + 1;
        itemSize[0] = expand_sch[tableID].item_size;
        itemSize[1] = itemSize[0] + sizeof(uint32_t);
        itemPerBlk[0] = expand_sch[tableID].item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        // create a temporary output table
        tmp_dst = dst_prefix + "_tmp_hid_" + std::to_string(tableID);
        tmp_conn = new FileSimulator(server_host, tmp_dst.c_str(), true);
        
        // assign HId values
        if (tableID == 0) expand_index = nAttrs[0] - 1;
        else if (tableID == 1) expand_index = nAttrs[0] - 2;
        offset = expand_sch[tableID].attrOffset[expand_index - 1];
        
        int32_t preH = -1;
        uint32_t preHId = 0;
        
        iItemSegSize[0] = expand_sch[tableID].attrOffset[expand_index];
        iItemSegSize[1] = itemSize[0] - iItemSegSize[0];
        oItemSegPos = iItemSegSize[0] + sizeof(uint32_t);
        
        oPItem = oBlock + META_BLOCK_SIZE;
        oBlockID = 0;
        oItemCnt = 0;
        oRItemCnt = 0;
        
        iLeft = 0;
        remain_block = expand_n_blocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(expand_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        char* curAttr = iPItem + offset;
                        int32_t curH;
                        memcpy(&curH, curAttr, sizeof(uint32_t));
                        assert(curH == 0 || curH == 1);
                        uint32_t curHId;
                        if (preH == -1) { // first item in the table
                            curHId = 1;
                            preH = curH;
                        }
                        else if (curH == preH)
                            curHId = preHId + 1;
                        else {
                            assert(preH == 1 && curH == 0);
                            curHId = 1;
                            preH = curH;
                        }
                        preHId = curHId;
                        
                        memcpy(oPItem, iPItem, iItemSegSize[0]);
                        memcpy(oPItem + iItemSegSize[0], &curHId, sizeof(uint32_t));
                        memcpy(oPItem + oItemSegPos, iPItem + iItemSegSize[0], iItemSegSize[1]);
                        
                        ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(tmp_conn, oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(tmp_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        expand_n_blocks[tableID] = oBlockID;
        
        // update the schema
        expand_sch[tableID].nAttrs = nAttrs[1];
        for (uint32_t i = nAttrs[0]; i > expand_index; --i) {
            uint32_t j = i - 1;
            expand_sch[tableID].attrType[i] = expand_sch[tableID].attrType[j];
            expand_sch[tableID].attrSize[i] = expand_sch[tableID].attrSize[j];
            expand_sch[tableID].attrOffset[i] = expand_sch[tableID].attrOffset[j] + sizeof(uint32_t);
        }
        expand_sch[tableID].attrType[expand_index] = INTEGER;
        expand_sch[tableID].attrSize[expand_index] = sizeof(uint32_t);
        expand_sch[tableID].item_size = itemSize[1];
        expand_sch[tableID].item_per_blk = itemPerBlk[1];
        
        // update the connector
        delete expand_conn[tableID];
        expand_conn[tableID] = tmp_conn;
    }
    
    void findPos(ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID, const uint32_t wDistLen, const uint32_t wDist[][2], const double wDistAvg) {
        // generate W_Attr information
        struct W_Attr {
            int32_t ID;
            int32_t TW;
            int32_t N;
            int32_t SumN;
            double Delta;
            double SumDelta;
        };
        
        uint32_t wHDistLen = 0, wLDistLen = 0;
        W_Attr wHDist[wDistLen];
        W_Attr wLDist[wDistLen];
        for (uint32_t i = 0; i < wDistLen; ++i) {
            if ((double)(wDist[i][0]) > wDistAvg) {
                wHDist[wHDistLen].TW = wDist[i][0];
                wHDist[wHDistLen].N = wDist[i][1];
                ++wHDistLen;
            }
            else {
                wLDist[wLDistLen].TW = wDist[i][0];
                wLDist[wLDistLen].N = wDist[i][1];
                ++wLDistLen;
            }
        }
        
        auto cmp = [](const void *arg1, const void *arg2)->int {
            return ((W_Attr *)arg2)->TW - ((W_Attr *)arg1)->TW;
        };
        if (wHDistLen > 0) std::qsort(wHDist, wHDistLen, sizeof(W_Attr), cmp);
        if (wLDistLen > 0) std::qsort(wLDist, wLDistLen, sizeof(W_Attr), cmp);
        
        uint32_t sumN = 0;
        double sumDelta = 0.0;
        for (uint32_t i = 0; i < wHDistLen; ++i) {
            wHDist[i].ID = i + 1;
            sumN += wHDist[i].N;
            wHDist[i].SumN = sumN;
            wHDist[i].Delta = wHDist[i].N * (wHDist[i].TW - wDistAvg);
            sumDelta += wHDist[i].Delta;
            wHDist[i].SumDelta = sumDelta;
        }
        sumN = 0;
        sumDelta = 0.0;
        for (uint32_t i = 0; i < wLDistLen; ++i) {
            wLDist[i].ID = i + 1;
            sumN += wLDist[i].N;
            wLDist[i].SumN = sumN;
            wLDist[i].Delta = wLDist[i].N * (wDistAvg - wLDist[i].TW);
            sumDelta += wLDist[i].Delta;
            wLDist[i].SumDelta = sumDelta;
        }
        
        // [0]: old schema; [1]: new schema;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = expand_sch[tableID].nAttrs;
        nAttrs[1] = nAttrs[0] + 1;
        itemSize[0] = expand_sch[tableID].item_size;
        itemSize[1] = itemSize[0] + sizeof(uint32_t);
        itemPerBlk[0] = expand_sch[tableID].item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        // create a temporary output table
        std::string tmp_dst = dst_prefix + "_tmp_find_pos_" + std::to_string(tableID);
        ServerConnector* tmp_conn = new FileSimulator(server_host, tmp_dst.c_str(), true);
        
        // generate Position column
        uint32_t expand_index = 0;
        if (tableID == 0) expand_index = nAttrs[0] - 1;
        else if (tableID == 1) expand_index = nAttrs[0] - 2;
        uint32_t offsetH = expand_sch[tableID].attrOffset[expand_index - 2];
        uint32_t offsetHId = expand_sch[tableID].attrOffset[expand_index - 1];
        
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = expand_sch[tableID].attrOffset[expand_index];
        iItemSegSize[1] = itemSize[0] - iItemSegSize[0];
        uint32_t oItemSegPos = iItemSegSize[0] + sizeof(uint32_t);
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = expand_n_blocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(expand_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        char* curAttr = iPItem + offsetH;
                        uint32_t curH;
                        memcpy(&curH, curAttr, sizeof(uint32_t));
                        assert(curH == 0 || curH == 1);
                        curAttr = iPItem + offsetHId;
                        uint32_t curHId;
                        memcpy(&curHId, curAttr, sizeof(uint32_t));
                        uint32_t curPos;
                        
                        if (curH == 1) {
                            uint32_t hIndex = 0;
                            for (hIndex = 0; hIndex < wHDistLen; ++hIndex)
                                if (wHDist[hIndex].SumN >= curHId)
                                    break;
                            uint32_t curP = curHId - 1 + wHDist[hIndex].N - wHDist[hIndex].SumN;
                            double curDelta = wHDist[hIndex].SumDelta - wHDist[hIndex].Delta + curP * (wHDist[hIndex].TW - wDistAvg);
                            uint32_t lIndex = 0;
                            for (lIndex = 0; lIndex < wLDistLen; ++lIndex)
                                if (wLDist[lIndex].SumDelta >= curDelta)
                                    break;
                            curPos = curHId + (uint32_t)((curDelta + wLDist[lIndex].Delta - wLDist[lIndex].SumDelta) / (wDistAvg - wLDist[lIndex].TW)) + wLDist[lIndex].SumN - wLDist[lIndex].N;
                        }
                        else if (curH == 0) {
                            uint32_t lIndex = 0;
                            for (lIndex = 0; lIndex < wLDistLen; ++lIndex)
                                if (wLDist[lIndex].SumN >= curHId)
                                    break;
                            uint32_t curP = curHId + wLDist[lIndex].N - wLDist[lIndex].SumN;
                            double curDelta = wLDist[lIndex].SumDelta - wLDist[lIndex].Delta + curP * (wDistAvg - wLDist[lIndex].TW);
                            uint32_t hIndex = 0;
                            for (hIndex = 0; hIndex < wHDistLen; ++hIndex)
                                if (wHDist[hIndex].SumDelta >= curDelta)
                                    break;
                            curPos = curHId + (uint32_t)ceil((curDelta + wHDist[hIndex].Delta - wHDist[hIndex].SumDelta) / (wHDist[hIndex].TW - wDistAvg)) + wHDist[hIndex].SumN - wHDist[hIndex].N;
                        }
                        
                        memcpy(oPItem, iPItem, iItemSegSize[0]);
                        memcpy(oPItem + iItemSegSize[0], &curPos, sizeof(uint32_t));
                        memcpy(oPItem + oItemSegPos, iPItem + iItemSegSize[0], iItemSegSize[1]);
                        
                        ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(tmp_conn, oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(tmp_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        expand_n_blocks[tableID] = oBlockID;
        
        // update the schema
        expand_sch[tableID].nAttrs = nAttrs[1];
        for (uint32_t i = nAttrs[0]; i > expand_index; --i) {
            uint32_t j = i - 1;
            expand_sch[tableID].attrType[i] = expand_sch[tableID].attrType[j];
            expand_sch[tableID].attrSize[i] = expand_sch[tableID].attrSize[j];
            expand_sch[tableID].attrOffset[i] = expand_sch[tableID].attrOffset[j] + sizeof(uint32_t);
        }
        expand_sch[tableID].attrType[expand_index] = INTEGER;
        expand_sch[tableID].attrSize[expand_index] = sizeof(uint32_t);
        expand_sch[tableID].item_size = itemSize[1];
        expand_sch[tableID].item_per_blk = itemPerBlk[1];
        
        // update the connector
        delete expand_conn[tableID];
        expand_conn[tableID] = tmp_conn;
    }
    
    void reorderPos(ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID) {
        // [0]: old schema; [1]: new schema;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = expand_sch[tableID].nAttrs;
        nAttrs[1] = nAttrs[0] - 3;
        itemSize[0] = expand_sch[tableID].item_size;
        itemSize[1] = itemSize[0] - 3 * sizeof(uint32_t);
        itemPerBlk[0] = expand_sch[tableID].item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        // obliviously sort the old table based on Position column
        uint32_t pos_index = 0;
        if (tableID == 0) pos_index = nAttrs[0] - 2;
        else if (tableID == 1) pos_index = nAttrs[0] - 3;
        CMP expand_cmp {1, {pos_index}, {0}};
        ObliviousSort* oblisort = new ObliviousSort(expand_conn[tableID], expand_n_blocks[tableID], &(expand_sch[tableID]), &expand_cmp);
        oblisort->BatcherSort();
        delete oblisort;
        
        // create a temporary output table
        std::string tmp_dst = dst_prefix + "_tmp_reorder_pos_" + std::to_string(tableID);
        ServerConnector* tmp_conn = new FileSimulator(server_host, tmp_dst.c_str(), true);
        
        // projection
        uint32_t iItemSegPos = expand_sch[tableID].attrOffset[pos_index + 1];
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = expand_sch[tableID].attrOffset[pos_index - 2];
        iItemSegSize[1] = itemSize[0] - iItemSegPos;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = expand_n_blocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(expand_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        memcpy(oPItem, iPItem, iItemSegSize[0]);
                        memcpy(oPItem + iItemSegSize[0], iPItem + iItemSegPos, iItemSegSize[1]);
                        
                        ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(tmp_conn, oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(tmp_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        expand_n_blocks[tableID] = oBlockID;
        
        // update the schema
        expand_sch[tableID].nAttrs = nAttrs[1];
        uint32_t curOffset = expand_sch[tableID].attrOffset[pos_index - 2];
        for (uint32_t i = pos_index - 2; i < nAttrs[1]; ++i) {
            uint32_t j = i + 3;
            expand_sch[tableID].attrType[i] = expand_sch[tableID].attrType[j];
            expand_sch[tableID].attrSize[i] = expand_sch[tableID].attrSize[j];
            expand_sch[tableID].attrOffset[i] = curOffset;
            curOffset += expand_sch[tableID].attrSize[i];
        }
        expand_sch[tableID].item_size = itemSize[1];
        expand_sch[tableID].item_per_blk = itemPerBlk[1];
        
        // update the connector
        delete expand_conn[tableID];
        expand_conn[tableID] = tmp_conn;
    }
    
    void makeCopy(ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID, const double wDistAvg) {
        // create a temporary output table
        std::string tmp_dst = dst_prefix + "_tmp_make_copy_" + std::to_string(tableID);
        ServerConnector* tmp_conn = new FileSimulator(server_host, tmp_dst.c_str(), true);
        
        uint32_t nAttrs = expand_sch[tableID].nAttrs;
        uint32_t itemSize = expand_sch[tableID].item_size;
        uint32_t itemPerBlk = expand_sch[tableID].item_per_blk;
        uint32_t offset = 0;
        if (tableID == 0) offset = expand_sch[tableID].attrOffset[nAttrs - 2];
        else if (tableID == 1) offset = expand_sch[tableID].attrOffset[nAttrs - 3];
        
        uint32_t totIn = 0;
        int32_t preTotOut = 0;
        struct Item {
            uint32_t itemID;
            int32_t TW;
            char* content;
        };
        std::queue <Item> leftItem;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        // make copy
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = expand_n_blocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(expand_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk; ++j) {
                    if (*iPItem == 'r') {
                        ++totIn;
                        int32_t curTotOut = (int32_t)(totIn * wDistAvg);
                        int32_t curOut = curTotOut - preTotOut;
                        preTotOut = curTotOut;
                        
                        char* pAttr = iPItem + offset;
                        int32_t curTW;
                        memcpy(&curTW, pAttr, sizeof(int32_t));
                        if (curTW <= curOut) {
                            for (uint32_t l = 0; l < curTW; ++l) {
                                memcpy(oPItem, iPItem, itemSize);
                                ++oRItemCnt;
                                ++oItemCnt;
                                if (oItemCnt < itemPerBlk)
                                    oPItem += itemSize;
                                else {
                                    memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                                    memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                                    insertBlock(tmp_conn, oBlockID, oBlock);
                                    ++oBlockID;
                                    oItemCnt = 0;
                                    oRItemCnt = 0;
                                    oPItem = oBlock + META_BLOCK_SIZE;
                                }
                            }
                            curOut -= curTW;
                        }
                        else {
                            Item newItem {totIn, curTW, NULL};
                            newItem.content = new char[itemSize];
                            memcpy(newItem.content, iPItem, itemSize);
                            leftItem.push(newItem);
                        }
                        while (curOut > 0) {
                            Item curItem = leftItem.front();
                            uint32_t copyNum = std::min(curItem.TW, curOut);
                            for (uint32_t l = 0; l < copyNum; ++l) {
                                memcpy(oPItem, curItem.content, itemSize);
                                ++oRItemCnt;
                                ++oItemCnt;
                                if (oItemCnt < itemPerBlk)
                                    oPItem += itemSize;
                                else {
                                    memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                                    memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                                    insertBlock(tmp_conn, oBlockID, oBlock);
                                    ++oBlockID;
                                    oItemCnt = 0;
                                    oRItemCnt = 0;
                                    oPItem = oBlock + META_BLOCK_SIZE;
                                }
                            }
                            if (curItem.TW > curOut) {
                                leftItem.front().TW -= curOut;
                                curOut = 0;
                            }
                            else {
                                curOut -= curItem.TW;
                                delete[] curItem.content;
                                leftItem.pop();
                            }
                        }
                    }
                    if (j < itemPerBlk - 1)
                        iPItem += itemSize;
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk)
                    oPItem += itemSize;
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(tmp_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        expand_n_blocks[tableID] = oBlockID;
        
        // update the connector
        delete expand_conn[tableID];
        expand_conn[tableID] = tmp_conn;
    }
    
    void attachFlag(ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID) {
        uint32_t nAttrs = expand_sch[tableID].nAttrs;
        uint32_t itemSize = expand_sch[tableID].item_size;
        uint32_t itemPerBlk = expand_sch[tableID].item_per_blk;
        
        // obliviously sort the old table based on Id column
        uint32_t id_index = 0;
        if (tableID == 0) id_index = nAttrs - 1;
        else if (tableID == 1) id_index = nAttrs - 2;
        CMP expand_cmp {1, {id_index}, {0}};
        ObliviousSort* oblisort = new ObliviousSort(expand_conn[tableID], expand_n_blocks[tableID], &(expand_sch[tableID]), &expand_cmp);
        oblisort->BatcherSort();
        delete oblisort;
        
        uint32_t w_offset = expand_sch[tableID].attrOffset[id_index - 2];
        uint32_t flag_offset = expand_sch[tableID].attrOffset[id_index - 1];
        uint32_t id_offset = expand_sch[tableID].attrOffset[id_index];
        
        // attach Flag column
        int32_t preID = -1;
        uint32_t preNum = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = expand_n_blocks[tableID];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(expand_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk; ++j) {
                    if (*iPItem == 'r') {
                        char* pAttr = iPItem + id_offset;
                        int32_t curID;
                        memcpy(&curID, pAttr, sizeof(int32_t));
                        uint32_t curNum;
                        if (curID == preID) curNum = preNum + 1;
                        else {
                            curNum = 1;
                            preID = curID;
                        }
                        preNum = curNum;
                        pAttr = iPItem + w_offset;
                        uint32_t curW;
                        memcpy(&curW, pAttr, sizeof(uint32_t));
                        
                        uint32_t curFlag;
                        if (curNum <= curW) curFlag = 1;
                        else curFlag = 0;
                        pAttr = iPItem + flag_offset;
                        memcpy(pAttr, &curFlag, sizeof(uint32_t));
                    }
                    if (j < itemPerBlk - 1)
                        iPItem += itemSize;
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            updateBlock(expand_conn[tableID], iLeft, iRight, buffer);
            iLeft = iRight;
        }
    }
    
    void filterDummy(ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t tableID) {
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = expand_sch[tableID].nAttrs;
        nAttrs[1] = nAttrs[0] - 2;
        itemSize[0] = expand_sch[tableID].item_size;
        itemSize[1] = itemSize[0] - 2 * sizeof(uint32_t);
        itemPerBlk[0] = expand_sch[tableID].item_per_blk;
        itemPerBlk[1] = (plain_len - META_BLOCK_SIZE) / itemSize[1];
        
        // obliviously sort the old table based on Id column
        uint32_t flag_index = 0;
        if (tableID == 0) flag_index = nAttrs[0] - 2;
        else if (tableID == 1) flag_index = nAttrs[0] - 3;
        CMP expand_cmp {2, {flag_index, flag_index + 1}, {1, 0}};
        ObliviousSort* oblisort = new ObliviousSort(expand_conn[tableID], expand_n_blocks[tableID], &(expand_sch[tableID]), &expand_cmp);
        oblisort->BatcherSort();
        delete oblisort;
        
        // create a temporary output table
        std::string tmp_dst = dst_prefix + "_tmp_filter_dummy_" + std::to_string(tableID);
        ServerConnector* tmp_conn = new FileSimulator(server_host, tmp_dst.c_str(), true);
        
        // projection
        uint32_t flag_offset = expand_sch[tableID].attrOffset[flag_index];
        
        uint32_t iItemSegPos = expand_sch[tableID].attrOffset[flag_index + 1];
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = expand_sch[tableID].attrOffset[flag_index - 1];
        iItemSegSize[1] = itemSize[0] - iItemSegPos;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = expand_n_blocks[tableID];
        bool iEnd = false;
        while (remain_block > 0) {
            if (iEnd) break;
            
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(expand_conn[tableID], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        char* iPAttr = iPItem + flag_offset;
                        uint32_t flag;
                        memcpy(&flag, iPAttr, sizeof(uint32_t));
                        assert(flag == 0 || flag == 1);
                        if (flag == 0) {
                            iEnd = true;
                            break;
                        }
                        memcpy(oPItem, iPItem, iItemSegSize[0]);
                        memcpy(oPItem + iItemSegSize[0], iPItem + iItemSegPos, iItemSegSize[1]);
                        ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(tmp_conn, oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[0] - 1)
                        iPItem += itemSize[0];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
                if (iEnd) break;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[1] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[1]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[1]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[1])
                    oPItem += itemSize[1];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(tmp_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        expand_n_blocks[tableID] = oBlockID;
        
        // update the schema
        expand_sch[tableID].nAttrs = nAttrs[1];
        uint32_t curOffset = expand_sch[tableID].attrOffset[flag_index - 1];
        for (uint32_t i = flag_index - 1; i < nAttrs[1]; ++i) {
            uint32_t j = i + 2;
            expand_sch[tableID].attrType[i] = expand_sch[tableID].attrType[j];
            expand_sch[tableID].attrSize[i] = expand_sch[tableID].attrSize[j];
            expand_sch[tableID].attrOffset[i] = curOffset;
            curOffset += expand_sch[tableID].attrSize[i];
        }
        expand_sch[tableID].item_size = itemSize[1];
        expand_sch[tableID].item_per_blk = itemPerBlk[1];
        
        // update the connector
        delete expand_conn[tableID];
        expand_conn[tableID] = tmp_conn;
    }
    
    void stitchTable(ServerConnector** expand_conn, const Schema* expand_sch, const uint32_t* expand_n_blocks) {
        // [0, ... , nInputs - 1]: expand_sch[0, ... , nInputs - 1]; [nInputs]: out_sch
        uint32_t nAttrs[nInputs + 1];
        uint32_t itemSize[nInputs + 1];
        uint32_t itemPerBlk[nInputs + 1];
        uint32_t attrSize[nInputs];
        uint32_t offset[nInputs];
        for (uint32_t j = 0; j < nInputs; ++j) {
            nAttrs[j] = expand_sch[j].nAttrs;
            itemSize[j] = expand_sch[j].item_size;
            itemPerBlk[j] = expand_sch[j].item_per_blk;
            uint32_t curAttrID = attrID[j];
            attrSize[j] = in_sch[j].attrSize[curAttrID];
            offset[j] = in_sch[j].attrOffset[curAttrID];
        }
        nAttrs[nInputs] = out_sch->nAttrs;
        itemSize[nInputs] = out_sch->item_size;
        itemPerBlk[nInputs] = out_sch->item_per_blk;
        
        // obliviously sort the expanded input tables based on join column and JId
        for (uint32_t j = 0; j < nInputs; ++j) {
            CMP expand_cmp {3, {attrID[j], nAttrs[j] - 1, nAttrs[j] - 2}, {0, 0, 0}};
            ObliviousSort* oblisort = new ObliviousSort(expand_conn[j], expand_n_blocks[j], &(expand_sch[j]), &expand_cmp);
            oblisort->BatcherSort();
            delete oblisort;
        }
        
        char iBlock[nInputs][B];
        char* iPItem[nInputs];
        uint32_t iBlockCnt[nInputs];
        uint32_t iItemCnt[nInputs];
        memset(iBlockCnt, 0, sizeof(uint32_t) * nInputs);
        memset(iItemCnt, 0, sizeof(uint32_t) * nInputs);
        for (uint32_t j = 0; j < nInputs; ++j) {
            readBlock(expand_conn[j], iBlockCnt[j], iBlock[j]);
            ++iBlockCnt[j];
            iPItem[j] = iBlock[j] + META_BLOCK_SIZE;
        }
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iItemSegPos[nInputs];
        uint32_t iItemSegSize[nInputs][2];
        uint32_t oItemSegPos[nInputs][2];
        iItemSegSize[0][0] = itemSize[0] - 2 * sizeof(uint32_t);
        iItemSegSize[0][1] = 0;
        oItemSegPos[0][0] = 0;
        oItemSegPos[0][1] = iItemSegSize[0][0];
        for (uint32_t j = 1; j < nInputs; ++j) {
            iItemSegPos[j] = offset[j] + attrSize[j];
            iItemSegSize[j][0] = offset[j] - sizeof(char);
            iItemSegSize[j][1] = itemSize[j] - iItemSegPos[j] - 2 * sizeof(uint32_t);
            oItemSegPos[j][0] = oItemSegPos[j - 1][1] + iItemSegSize[j - 1][1];
            oItemSegPos[j][1] = oItemSegPos[j][0] + iItemSegSize[j][0];
        }
        
        // stitch the expanded tables
        /****************/
        uint32_t test_real_item_num = 0;
        uint32_t out_offset = out_sch->attrOffset[attrID[0]];
        /****************/
        bool iEnd = false;
        while (!iEnd) {
            if (*iPItem[0] == 'r') {
                for (uint32_t j = 1; j < nInputs; ++j)
                    assert(*iPItem[j] == 'r');
                
                memcpy(oPItem, iPItem[0], iItemSegSize[0][0]);
                for (uint32_t j = 1; j < nInputs; ++j) {
                    memcpy(oPItem + oItemSegPos[j][0], iPItem[j] + sizeof(char), iItemSegSize[j][0]);
                    memcpy(oPItem + oItemSegPos[j][1], iPItem[j] + iItemSegPos[j], iItemSegSize[j][1]);
                }
                /****************/
                char testString[256];
                memcpy(testString, oPItem + out_offset, 256);
                testString[255] = '\0';
                printf("%s\n", testString);
                ++test_real_item_num;
                /****************/
                ++oRItemCnt;
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[nInputs])
                    oPItem += itemSize[nInputs];
                else {
                    memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                    memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                    insertBlock(out_conn, oBlockID, oBlock);
                    ++oBlockID;
                    oItemCnt = 0;
                    oRItemCnt = 0;
                    oPItem = oBlock + META_BLOCK_SIZE;
                }
                
                for (uint32_t j = 0; j < nInputs; ++j) {
                    if (iItemCnt[j] < itemPerBlk[j] - 1) {
                        iPItem[j] += itemSize[j];
                        ++iItemCnt[j];
                    }
                    else {
                        if (iBlockCnt[j] >= expand_n_blocks[j]) {
                            iEnd = true;
                            break;
                        }
                        readBlock(expand_conn[j], iBlockCnt[j], iBlock[j]);
                        ++iBlockCnt[j];
                        iItemCnt[j] = 0;
                        iPItem[j] = iBlock[j] + META_BLOCK_SIZE;
                    }
                }
            }
            else iEnd = true;
        }
        if (oItemCnt > 0) {
            char oItem[plain_len];
            char* oSItem = oItem + sizeof(char);
            uint32_t oItemContentSize = itemSize[nInputs] - sizeof(char);
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[nInputs]) {
                std::string rnd_str = generate_random_block(oItemContentSize);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
                memcpy(oPItem, oItem, itemSize[nInputs]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[nInputs])
                    oPItem += itemSize[nInputs];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(out_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        out_tot_block = oBlockID;
        /****************/
        printf("\nThe number of real items in result table: %u\n\n", test_real_item_num);
        /****************/
    }
    
    // # of input tables
    const uint32_t nInputs = 2;
    // server connectors of input tables
    ServerConnector** in_conn = NULL;
    // # of blocks in input tables
    uint32_t* nBlocks = NULL;
    // schema of input tables
    Schema* in_sch = NULL;
    // attribute IDs to be joined
    int32_t* attrID = NULL;
    
    // output table (file)
    ServerConnector* out_conn = NULL;
    // prefix of output collection name
    std::string dst_prefix;
    // schema of output table
    Schema* out_sch = NULL;
    // cmp of output table (only for ORDER BY)
    CMP* out_cmp = NULL;
    // # of blocks in output table
    uint32_t out_tot_block;
};

#endif //__OQP_H__

