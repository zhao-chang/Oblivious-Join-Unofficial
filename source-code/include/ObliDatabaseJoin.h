#ifndef __OBLI_DATABASE_JOIN_H__
#define __OBLI_DATABASE_JOIN_H__

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

class ObliDatabaseJoin {
public:
    // load
    ObliDatabaseJoin(const uint32_t mode, const std::string scale, const std::string& dst, const std::string& prefix) {
        nInputs = getTableNum(mode);
        iItemNum = new uint32_t[nInputs];
        in_sch = new Schema[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i) {
            generateInputSchema(mode, i, &(in_sch[i]));
            getInputItemNum(mode, scale, i, iItemNum[i]);
        }
        
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
            importData(mode, scale, i, &(in_sch[i]), blocks);
            nBlocks[i] = blocks.size();
            std::vector <std::pair<uint32_t, std::string>> insert_buffer;
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
        initBufferSize(2);
        
        /***********************/
        comm_size = 0;
        /***********************/
        
        // write info
        std::string odj_name = prefix + std::string("_odj_info.txt");
        FILE* fp = fopen(odj_name.c_str(), "w");
        fprintf(fp, "%d", nBlocks[0]);
        for (uint32_t i = 1; i < nInputs; ++i)
            fprintf(fp, " %d", nBlocks[i]);
        fprintf(fp, "\n");
        fprintf(fp, "%d", iItemNum[0]);
        for (uint32_t i = 1; i < nInputs; ++i)
            fprintf(fp, " %d", iItemNum[i]);
        fprintf(fp, "\n");
        fprintf(fp, "%s\n", dst_prefix.c_str());
        fclose(fp);
    }
    
    // run
    ObliDatabaseJoin(const uint32_t mode, const std::string& prefix) {
        /***********************/
        comm_size = 0;
        /***********************/
        
        nInputs = getTableNum(mode);
        in_sch = new Schema[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i)
            generateInputSchema(mode, i, &(in_sch[i]));
        
        // read info
        nBlocks = new uint32_t[nInputs];
        iItemNum = new uint32_t[nInputs];
        std::string odj_name = prefix + std::string("_odj_info.txt");
        FILE* fp = fopen(odj_name.c_str(), "r");
        for (uint32_t i = 0; i < nInputs; ++i)
            fscanf(fp, "%d", nBlocks + i);
        for (uint32_t i = 0; i < nInputs; ++i)
            fscanf(fp, "%d", iItemNum + i);
        char out_dst_prefix[100];
        fscanf(fp, "%s", out_dst_prefix);
        dst_prefix = std::string(out_dst_prefix);
        fclose(fp);
        
        std::string key_name = prefix + std::string("_key.txt");
        initKey(key_name.c_str(), false);
        
        in_conn = new ServerConnector*[nInputs];
        for (uint32_t i = 0; i < nInputs; ++i) {
            std::string cur_dst = dst_prefix + "_input_" + std::to_string(i);
            in_conn[i] = new FileSimulator(server_host, cur_dst, false);
        }
        initBufferSize(2);
    }
    
    ~ObliDatabaseJoin() {
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
        if (iItemNum != NULL) {
            delete[] iItemNum;
            iItemNum = NULL;
        }
    }
    
    void ObliEquiJoin(const uint32_t n_tables, const uint32_t* table_id, const int32_t* parent_id, const int32_t* attr_id, const int32_t* parent_attr_id, const uint32_t n_proj_cols, const uint32_t proj_id[][MAX_COLS]) {
        assert(n_tables == 2);
        nTables = n_tables;
        iTable = new uint32_t[nTables];
        attrID = new int32_t[nTables];
        iParent = new int32_t[nTables];
        pAttrID = new int32_t[nTables];
        memcpy(iTable, table_id, nTables * sizeof(uint32_t));
        memcpy(attrID, attr_id, nTables * sizeof(int32_t));
        memcpy(iParent, parent_id, nTables * sizeof(int32_t));
        memcpy(pAttrID, parent_attr_id, nTables * sizeof(int32_t));
        attrID[0] = pAttrID[1];
        nProjCols = n_proj_cols;
        if (nProjCols > 0) {
            for (uint32_t i = 0; i < 2; ++i) {
                projID[i] = new uint32_t[nProjCols];
                memcpy(projID[i], proj_id[i], nProjCols * sizeof(uint32_t));
            }
        }
        
        std::string out_dst = dst_prefix + "_output";
        for (uint32_t i = 0; i < nTables; ++i)
            out_dst += "_" + std::to_string(iTable[i]);
        for (uint32_t i = 1; i < nTables; ++i) {
            out_dst += "_" + std::to_string(attrID[i]);
            out_dst += "_" + std::to_string(pAttrID[i]);
        }
        out_conn = new FileSimulator(server_host, out_dst.c_str(), true);
        out_sch = new Schema;
        if (nProjCols == 0) generateOutputSchema(in_sch, nTables, iTable, attrID, out_sch);
        else generateOutputSchema(in_sch, nTables, iTable, nProjCols, projID, out_sch);
        // TODO: support ORDER-BY operator
        if (nProjCols == 0) out_cmp = new CMP({1, {pAttrID[1]}, {0}});
        else out_cmp = new CMP({1, {1}, {0}});
        
        // initialize the buffer
        initBuffer();
        
        // create input tables with join degrees
        ServerConnector** weight_conn = new ServerConnector*[nTables];
        Schema* weight_sch = new Schema[nTables];
        uint32_t* weight_n_blocks = new uint32_t[nTables];
        for (uint32_t i = 0; i < nTables; ++i) {
            std::string weight_dst = dst_prefix + "_weight_" + std::to_string(i);
            weight_conn[i] = new FileSimulator(server_host, weight_dst.c_str(), true);
        }
        
        // compute the join degrees of input tables
        augmentTables(weight_conn, weight_sch, weight_n_blocks);
        
        if (est_real_item_num > 0) {
            // create expanded input tables
            ServerConnector** expand_conn = new ServerConnector*[nTables];
            Schema* expand_sch = new Schema[nTables];
            uint32_t* expand_n_blocks = new uint32_t[nTables];
            for (uint32_t i = 0; i < nTables; ++i) {
                std::string expand_dst = dst_prefix + "_expand_" + std::to_string(i);
                expand_conn[i] = new FileSimulator(server_host, expand_dst.c_str(), true);
            }
            
            // expand the input tables based on the degrees
            for (uint32_t i = 0; i < nTables; ++i)
                expandTable(weight_conn, weight_sch, weight_n_blocks, expand_conn, expand_sch, expand_n_blocks, i);
            
            // stitch expanded input tables
            stitchTable(expand_conn, expand_sch, expand_n_blocks);
            
            for (uint32_t i = 0; i < nTables; ++i) {
                delete expand_conn[i];
            }
            delete[] expand_conn;
            delete[] expand_sch;
            delete[] expand_n_blocks;
        }
        else out_tot_block = 0;
        
        for (uint32_t i = 0; i < nTables; ++i) {
            delete weight_conn[i];
        }
        delete[] weight_conn;
        delete[] weight_sch;
        delete[] weight_n_blocks;
        
        // destroy the buffer
        destroyBuffer();
        
        /***********************/
        // print join result
        if (nProjCols == 0) {
            printResult(out_conn, out_tot_block, out_sch, 1, &(pAttrID[1]));
        }
        else {
            int32_t printAttrID[nProjCols];
            for (int32_t i = 0; i < nProjCols; ++i)
                printAttrID[i] = i + 1;
            printResult(out_conn, out_tot_block, out_sch, nProjCols, printAttrID);
        }
        /***********************/
        
        if (iTable != NULL) {
            delete[] iTable;
            iTable = NULL;
        }
        if (attrID != NULL) {
            delete[] attrID;
            attrID = NULL;
        }
        if (iParent != NULL) {
            delete[] iParent;
            iParent = NULL;
        }
        if (pAttrID != NULL) {
            delete[] pAttrID;
            pAttrID = NULL;
        }
        if (nProjCols > 0) {
            for (uint32_t i = 0; i < 2; ++i) {
                if (projID[i] != NULL) {
                    delete[](projID[i]);
                    projID[i] = NULL;
                }
            }
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
    
    void ObliBandJoin(const uint32_t n_tables, const uint32_t* table_id, const int32_t* parent_id, const int32_t* attr_id, const int32_t* parent_attr_id, const double* band_range, const uint32_t n_proj_cols, const uint32_t proj_id[][MAX_COLS]) {
        assert(n_tables == 2);
        nTables = n_tables;
        iTable = new uint32_t[nTables];
        attrID = new int32_t[nTables];
        iParent = new int32_t[nTables];
        pAttrID = new int32_t[nTables];
        bRange = new double[2];
        memcpy(iTable, table_id, nTables * sizeof(uint32_t));
        memcpy(attrID, attr_id, nTables * sizeof(int32_t));
        memcpy(iParent, parent_id, nTables * sizeof(int32_t));
        memcpy(pAttrID, parent_attr_id, nTables * sizeof(int32_t));
        attrID[0] = pAttrID[1];
        memcpy(bRange, band_range, 2 * sizeof(double));
        nProjCols = n_proj_cols;
        if (nProjCols > 0) {
            for (uint32_t i = 0; i < 2; ++i) {
                projID[i] = new uint32_t[nProjCols];
                memcpy(projID[i], proj_id[i], nProjCols * sizeof(uint32_t));
            }
        }
        
        std::string out_dst = dst_prefix + "_output";
        for (uint32_t i = 0; i < nTables; ++i)
            out_dst += "_" + std::to_string(iTable[i]);
        for (uint32_t i = 1; i < nTables; ++i) {
            out_dst += "_" + std::to_string(attrID[i]);
            out_dst += "_" + std::to_string(pAttrID[i]);
        }
        for (uint32_t i = 0; i < 2; ++i)
            out_dst += "_" + std::to_string(bRange[i]);
        out_conn = new FileSimulator(server_host, out_dst.c_str(), true);
        out_sch = new Schema;
        if (nProjCols == 0) generateOutputSchema(in_sch, nTables, iTable, attrID, out_sch, false);
        else generateOutputSchema(in_sch, nTables, iTable, nProjCols, projID, out_sch);
        // TODO: support ORDER-BY operator
        if (nProjCols == 0) out_cmp = new CMP({1, {pAttrID[1]}, {0}});
        else out_cmp = new CMP({1, {1}, {0}});
        
        // initialize the buffer
        initBuffer();
        
        // create input tables with join degrees
        ServerConnector** weight_conn = new ServerConnector*[nTables];
        Schema* weight_sch = new Schema[nTables];
        uint32_t* weight_n_blocks = new uint32_t[nTables];
        for (uint32_t i = 0; i < nTables; ++i) {
            std::string weight_dst = dst_prefix + "_weight_" + std::to_string(i);
            weight_conn[i] = new FileSimulator(server_host, weight_dst.c_str(), true);
        }
        
        // compute the join degrees of input tables
        augmentBandTables(weight_conn, weight_sch, weight_n_blocks);
        
        if (est_real_item_num > 0) {
            // create expanded input tables
            ServerConnector** expand_conn = new ServerConnector*[nTables];
            Schema* expand_sch = new Schema[nTables];
            uint32_t* expand_n_blocks = new uint32_t[nTables];
            for (uint32_t i = 0; i < nTables; ++i) {
                std::string expand_dst = dst_prefix + "_expand_" + std::to_string(i);
                expand_conn[i] = new FileSimulator(server_host, expand_dst.c_str(), true);
            }
            
            // expand the input tables based on the degrees
            for (uint32_t i = 0; i < nTables; ++i)
                expandBandTable(weight_conn, weight_sch, weight_n_blocks, expand_conn, expand_sch, expand_n_blocks, i);
            
            // stitch expanded input tables
            stitchBandTable(expand_conn, expand_sch, expand_n_blocks);
            
            for (uint32_t i = 0; i < nTables; ++i) {
                delete expand_conn[i];
            }
            delete[] expand_conn;
            delete[] expand_sch;
            delete[] expand_n_blocks;
        }
        else out_tot_block = 0;
        
        for (uint32_t i = 0; i < nTables; ++i) {
            delete weight_conn[i];
        }
        delete[] weight_conn;
        delete[] weight_sch;
        delete[] weight_n_blocks;
        
        // destroy the buffer
        destroyBuffer();
        
        /***********************/
        // print join result
        if (nProjCols == 0) {
            printResult(out_conn, out_tot_block, out_sch, 1, &(pAttrID[1]));
        }
        else {
            int32_t printAttrID[nProjCols];
            for (int32_t i = 0; i < nProjCols; ++i)
                printAttrID[i] = i + 1;
            printResult(out_conn, out_tot_block, out_sch, nProjCols, printAttrID);
        }
        /***********************/
        
        if (iTable != NULL) {
            delete[] iTable;
            iTable = NULL;
        }
        if (attrID != NULL) {
            delete[] attrID;
            attrID = NULL;
        }
        if (iParent != NULL) {
            delete[] iParent;
            iParent = NULL;
        }
        if (pAttrID != NULL) {
            delete[] pAttrID;
            pAttrID = NULL;
        }
        if (bRange != NULL) {
            delete[] bRange;
            bRange = NULL;
        }
        if (nProjCols > 0) {
            for (uint32_t i = 0; i < 2; ++i) {
                if (projID[i] != NULL) {
                    delete[] (projID[i]);
                    projID[i] = NULL;
                }
            }
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
    
    /***********************/
    uint32_t getDataBlockNum() const {
        uint32_t block_num = 0;
        for (uint32_t i = 0; i < nInputs; ++i)
            block_num += nBlocks[i];
        return block_num;
    }
    
    size_t getServerSize() const {
        size_t server_size = (size_t)getDataBlockNum() * (size_t)B;
        return server_size;
    }
    
    size_t getClientSize() const {
        size_t client_size = oblisort_mem;
        return client_size;
    }
    
    size_t getCommSize() const {
        size_t total_comm_size = comm_size * (size_t)B;
        return total_comm_size;
    }
    
    void resetCommSize() {
        comm_size = 0;
    }
    /***********************/
    
private:
    /***********************/
    //The followings are for band joins
    void augmentBandTables(ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks) {
        for (uint32_t index = 0; index < nTables; ++index) {
            // obliviously sort input tables
            uint32_t tableID = iTable[index];
            CMP in_cmp {1, {attrID[index]}, {0}};
            ObliviousSort* oblisort = new ObliviousSort(in_conn[tableID], nBlocks[tableID], &(in_sch[tableID]), &in_cmp);
            comm_size += oblisort->BatcherSort();
            delete oblisort;
            
            // unify input tables
            std::string union_dst = dst_prefix + "_union";
            ServerConnector* union_conn = new FileSimulator(server_host, union_dst.c_str(), true);
            Schema* union_sch = new Schema; // schema of unified input table
            uint32_t union_n_blocks; // # of blocks in unified input table
            
            unifyTables(union_conn, union_sch, union_n_blocks, index);
            fillPositions(union_conn, union_sch, union_n_blocks, index);
            
            // separate the unified table
            ServerConnector** separate_conn = new ServerConnector*[nTables];
            Schema* separate_sch = new Schema; // schema of separate table
            uint32_t separate_n_blocks; // # of blocks in separate table
            for (uint32_t j = 0; j < nTables; ++j) {
                std::string separate_dst = dst_prefix + "_separate_" + std::to_string(j);
                separate_conn[j] = new FileSimulator(server_host, separate_dst.c_str(), true);
            }
            separateTables(union_conn, union_sch, union_n_blocks, separate_conn, separate_sch, separate_n_blocks, index);
            
            delete union_conn;
            delete union_sch;
            
            // generate the input table with join degrees
            genWTables(separate_conn, separate_sch, separate_n_blocks, weight_conn, weight_sch, weight_n_blocks, index);
            
            for (uint32_t j = 0; j < nTables; ++j) {
                delete separate_conn[j];
            }
            delete[] separate_conn;
            delete separate_sch;
        }
    }
    
    void unifyTables(ServerConnector* union_conn, Schema* union_sch, uint32_t & union_n_blocks, const uint32_t union_index) {
        // [0]: in_sch[iTable[0]]; [1]: in_sch[iTable[1]]; [2]: union_sch;
        uint32_t nAttrs[3];
        uint32_t itemSize[3];
        uint32_t itemPerBlk[3];
        for (uint32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
            nAttrs[index] = in_sch[tableID].nAttrs;
            itemSize[index] = in_sch[tableID].item_size;
            itemPerBlk[index] = in_sch[tableID].item_per_blk;
        }
        
        // initialize the schema of the unified table
        const uint32_t append_num = 3; // id, tid, pos
        union_sch->nAttrs = nAttrs[0] + nAttrs[1] + append_num - 2; // flag column, joined column
        uint32_t attr_index = nAttrs[0];
        memcpy(union_sch->attrType, in_sch[iTable[0]].attrType, sizeof(ATTR_TYPE) * attr_index);
        memcpy(union_sch->attrSize, in_sch[iTable[0]].attrSize, sizeof(uint32_t) * attr_index);
        memcpy(union_sch->attrOffset, in_sch[iTable[0]].attrOffset, sizeof(uint32_t) * attr_index);
        uint32_t union_item_size = itemSize[0];
        for (uint32_t j = 1; j < nAttrs[1]; ++j) {
            if (j == attrID[1]) continue;
            union_sch->attrType[attr_index] = in_sch[iTable[1]].attrType[j];
            union_sch->attrSize[attr_index] = in_sch[iTable[1]].attrSize[j];
            union_sch->attrOffset[attr_index] = union_item_size;
            union_item_size += union_sch->attrSize[attr_index];
            ++attr_index;
        }
        for (uint32_t j = 0; j < append_num; ++j) {
            union_sch->attrType[attr_index] = INTEGER;
            union_sch->attrSize[attr_index] = sizeof(uint32_t);
            union_sch->attrOffset[attr_index] = union_item_size;
            union_item_size += sizeof(uint32_t);
            ++attr_index;
        }
        union_sch->item_size = union_item_size;
        union_sch->item_per_blk = (plain_len - META_BLOCK_SIZE) / union_item_size;
        
        nAttrs[2] = union_sch->nAttrs;
        itemSize[2] = union_sch->item_size;
        itemPerBlk[2] = union_sch->item_per_blk;
        
        // construct the unified table
        uint32_t attr_size[nTables];
        ATTR_TYPE attr_type[nTables];
        uint32_t offset[nTables];
        uint32_t id_offset = itemSize[2] - append_num * sizeof(uint32_t);
        uint32_t tid_offset = id_offset + sizeof(uint32_t);
        for (uint32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
            uint32_t curAttrID = attrID[index];
            attr_size[index] = in_sch[tableID].attrSize[curAttrID];
            attr_type[index] = in_sch[tableID].attrType[curAttrID];
            offset[index] = in_sch[tableID].attrOffset[curAttrID];
            assert(attr_type[index] == INTEGER || attr_type[index] == DOUBLE);
        }
        
        uint32_t iItemSegPos = offset[1] + attr_size[1];
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = offset[1] - sizeof(char);
        iItemSegSize[1] = itemSize[1] - iItemSegPos;
        uint32_t oItemSegPos = itemSize[0] + iItemSegSize[0];
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        for (int32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
            uint32_t iLeft = 0;
            uint32_t iRight;
            uint32_t remain_block = nBlocks[tableID];
            uint32_t item_id = 0;
            uint32_t oitems;
            if (index == union_index) oitems = 2;
            else oitems = 1;
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
                    for (uint32_t j = 0; j < itemPerBlk[index]; ++j) {
                        if (*iPItem == 'r') {
                            for (int32_t l = 0; l < oitems; ++l) {
                                if (index == 0) {
                                    memcpy(oPItem, iPItem, itemSize[0]);
                                }
                                else if (index == 1) {
                                    oPItem[0] = 'r';
                                    memcpy(oPItem + offset[0], iPItem + offset[1], attr_size[1]);
                                    memcpy(oPItem + itemSize[0], iPItem + sizeof(char), iItemSegSize[0]);
                                    memcpy(oPItem + oItemSegPos, iPItem + iItemSegPos, iItemSegSize[1]);
                                }
                                if (index == union_index) {
                                    if(attr_type[index] == INTEGER) {
                                        int32_t join_key;
                                        memcpy(&join_key, oPItem + offset[0], sizeof(int32_t));
                                        if (l == 0) join_key -= (int32_t)(bRange[index]);
                                        else if (l == 1) join_key += (int32_t)(bRange[1 - index]);
                                        memcpy(oPItem + offset[0], &join_key, sizeof(int32_t));
                                    }
                                    else if(attr_type[index] == DOUBLE) {
                                        double join_key;
                                        memcpy(&join_key, oPItem + offset[0], sizeof(double));
                                        if (l == 0) join_key -= bRange[index];
                                        else if (l == 1) join_key += bRange[1 - index];
                                        memcpy(oPItem + offset[0], &join_key, sizeof(double));
                                    }
                                }
                                memcpy(oPItem + id_offset, &item_id, sizeof(uint32_t));
                                if (index == union_index) {
                                    int32_t tid;
                                    if (l == 0) tid = -nTables;
                                    else if (l == 1) tid = nTables;
                                    memcpy(oPItem + tid_offset, &tid, sizeof(int32_t));
                                }
                                else memcpy(oPItem + tid_offset, &index, sizeof(int32_t));
                                
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
                            ++item_id;
                        }
                        if (j < itemPerBlk[index] - 1)
                            iPItem += itemSize[index];
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
        
        /***********************/
        for (uint32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
            comm_size += nBlocks[tableID];
        }
        comm_size += union_n_blocks;
        /***********************/
        
        // obliviously sort the unified table
        CMP union_cmp {3, {attrID[0], nAttrs[2] - 2, nAttrs[2] - 3}, {0, 0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(union_conn, union_n_blocks, union_sch, &union_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
    }
    
    void fillPositions(ServerConnector* union_conn, Schema* union_sch, const uint32_t union_n_blocks, const uint32_t union_index) {
        uint32_t nAttrs = union_sch->nAttrs;
        uint32_t itemSize = union_sch->item_size;
        uint32_t itemPerBlk = union_sch->item_per_blk;
        ATTR_TYPE attrType = union_sch->attrType[attrID[0]];
        uint32_t offset = union_sch->attrOffset[attrID[0]];
        uint32_t attrSize = union_sch->attrSize[attrID[0]];
        
        // [0]: id, [1]: tid, [2]: pos
        const uint32_t append_num = 3;
        uint32_t append_offset[append_num];
        for (uint32_t i = 0, j = nAttrs - append_num; i < append_num; ++i, ++j)
            append_offset[i] = union_sch->attrOffset[j];
        
        int32_t position = 0;
        
        uint32_t iLeft, iRight;
        uint32_t remain_block = union_n_blocks;
        iLeft = 0;
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
                        int32_t cur_tid;
                        memcpy(&cur_tid, iPItem + append_offset[1], sizeof(int32_t));
                        if (cur_tid == 1 - union_index) ++position;
                        else {
                            if (cur_tid > 0) cur_tid = -1;
                            memcpy(iPItem + append_offset[1], &cur_tid, sizeof(int32_t));
                            memcpy(iPItem + append_offset[2], &position, sizeof(int32_t));
                        }
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
        /***********************/
        comm_size += 2 * union_n_blocks;
        /***********************/
        
        // obliviously sort the unified table
        CMP union_cmp {2, {nAttrs - 2, nAttrs - 3}, {0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(union_conn, union_n_blocks, union_sch, &union_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
    }
    
    void separateTables(ServerConnector* union_conn, Schema* union_sch, const uint32_t union_n_blocks, ServerConnector** separate_conn, Schema* separate_sch, uint32_t & separate_n_blocks, const uint32_t union_index) {
        // initialize the schema of separate tables
        uint32_t append_num = 2; // [0]: id, [1]: pos
        uint32_t tableID = iTable[union_index];
        uint32_t attr_index = in_sch[tableID].nAttrs;
        separate_sch->nAttrs = attr_index + append_num;
        memcpy(separate_sch->attrType, in_sch[tableID].attrType, sizeof(ATTR_TYPE) * attr_index);
        memcpy(separate_sch->attrSize, in_sch[tableID].attrSize, sizeof(uint32_t) * attr_index);
        memcpy(separate_sch->attrOffset, in_sch[tableID].attrOffset, sizeof(uint32_t) * attr_index);
        uint32_t separate_item_size = in_sch[tableID].item_size;
        for (uint32_t i = 0; i < append_num; ++i) {
            separate_sch->attrType[attr_index] = INTEGER;
            separate_sch->attrSize[attr_index] = sizeof(uint32_t);
            separate_sch->attrOffset[attr_index] = separate_item_size;
            separate_item_size += sizeof(uint32_t);
            ++attr_index;
        }
        separate_sch->item_size = separate_item_size;
        separate_sch->item_per_blk = (plain_len - META_BLOCK_SIZE) / separate_item_size;
        
        // [0]: union_sch; [1]: separate_sch;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = union_sch->nAttrs;
        itemSize[0] = union_sch->item_size;
        itemPerBlk[0] = union_sch->item_per_blk;
        nAttrs[1] = separate_sch->nAttrs;
        itemSize[1] = separate_sch->item_size;
        itemPerBlk[1] = separate_sch->item_per_blk;
        
        // [0]: in_sch[iTable[0]]; [1]: in_sch[iTable[1]];
        uint32_t inItemSize[2];
        uint32_t offset[2];
        uint32_t attrSize[2];
        for (uint32_t i = 0; i < nTables; ++i) {
            int32_t iTableID = iTable[i];
            inItemSize[i] = in_sch[iTableID].item_size;
            offset[i] = in_sch[iTableID].attrOffset[attrID[i]];
            attrSize[i] = in_sch[iTableID].attrSize[attrID[i]];
        }
        
        append_num = 3; // [0]: id, [1]: tid, [2]: pos
        uint32_t append_offset[append_num];
        for (uint32_t i = 0, j = nAttrs[0] - append_num; i < append_num; ++i, ++j)
            append_offset[i] = union_sch->attrOffset[j];
        
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = offset[1] - sizeof(char);
        uint32_t iItemSegPos = inItemSize[0] + iItemSegSize[0];
        uint32_t oItemSegPos = offset[1] + attrSize[0];
        iItemSegSize[1] = inItemSize[1] - oItemSegPos;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        char oItem[plain_len];
        char* oSItem = oItem + sizeof(char);
        uint32_t oItemContentSize = itemSize[1] - sizeof(char);
        
        // generate the separate tables
        bool iEnd = false;
        int32_t pre_tid = -nTables;
        uint32_t separate_index = 0;
        char* iPBlock = buffer + META_BLOCK_SIZE;
        uint32_t union_block_id = 0;
        for (union_block_id = 0; union_block_id < union_n_blocks; ++union_block_id) {
            readBlock(union_conn, union_block_id, buffer);
            char* iPItem = iPBlock;
            for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                if (*iPItem == 'r') {
                    int32_t cur_tid;
                    memcpy(&cur_tid, iPItem + append_offset[1], sizeof(int32_t));
                    if (cur_tid != pre_tid) {
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
                            insertBlock(separate_conn[separate_index], oBlockID, oBlock);
                            ++oBlockID;
                        }
                        separate_n_blocks = oBlockID;
                        oPItem = oBlock + META_BLOCK_SIZE;
                        oBlockID = 0;
                        oItemCnt = 0;
                        oRItemCnt = 0;
                        
                        if (cur_tid >= 0) {
                            iEnd = true;
                            break;
                        }
                        else {
                            assert(cur_tid == -1);
                            separate_index = 1;
                            pre_tid = cur_tid;
                        }
                    }
                    if (union_index == 0)
                        memcpy(oItem, iPItem, inItemSize[0]);
                    else if (union_index == 1) {
                        oItem[0] = 'r';
                        memcpy(oItem + sizeof(char), iPItem + inItemSize[0], iItemSegSize[0]);
                        memcpy(oItem + offset[1], iPItem + offset[0], attrSize[0]);
                        memcpy(oItem + oItemSegPos, iPItem + iItemSegPos, iItemSegSize[1]);
                    }
                    memcpy(oItem + inItemSize[union_index], iPItem + append_offset[0], sizeof(uint32_t));
                    memcpy(oItem + inItemSize[union_index] + sizeof(uint32_t), iPItem + append_offset[2], sizeof(uint32_t));
                    ++oRItemCnt;
                    memcpy(oPItem, oItem, itemSize[1]);
                    ++oItemCnt;
                    if (oItemCnt < itemPerBlk[1])
                        oPItem += itemSize[1];
                    else {
                        memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                        memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                        insertBlock(separate_conn[separate_index], oBlockID, oBlock);
                        ++oBlockID;
                        oItemCnt = 0;
                        oRItemCnt = 0;
                        oPItem = oBlock + META_BLOCK_SIZE;
                    }
                }
                if (j < itemPerBlk[0] - 1)
                    iPItem += itemSize[0];
            }
            if (iEnd) break;
        }
        
        /***********************/
        comm_size += union_block_id + 1 + 2 * separate_n_blocks;
        /***********************/
    }
    
    void genWTables(ServerConnector** separate_conn, Schema* separate_sch, const uint32_t separate_n_blocks, ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks, const uint32_t weight_index) {
        // initialize the schema of input table with join degrees
        uint32_t append_num = 1; // [0]: alpha
        uint32_t attr_index = separate_sch->nAttrs;
        weight_sch[weight_index].nAttrs = attr_index + append_num;
        memcpy(weight_sch[weight_index].attrType, separate_sch->attrType, sizeof(ATTR_TYPE) * attr_index);
        memcpy(weight_sch[weight_index].attrSize, separate_sch->attrSize, sizeof(uint32_t) * attr_index);
        memcpy(weight_sch[weight_index].attrOffset, separate_sch->attrOffset, sizeof(uint32_t) * attr_index);
        uint32_t weight_item_size = separate_sch->item_size;
        for (uint32_t i = 0; i < append_num; ++i) {
            weight_sch[weight_index].attrType[attr_index] = INTEGER;
            weight_sch[weight_index].attrSize[attr_index] = sizeof(uint32_t);
            weight_sch[weight_index].attrOffset[attr_index] = weight_item_size;
            weight_item_size += sizeof(uint32_t);
            ++attr_index;
        }
        weight_sch[weight_index].item_size = weight_item_size;
        weight_sch[weight_index].item_per_blk = (plain_len - META_BLOCK_SIZE) / weight_item_size;
        
        // [0]: separate_sch; [1]: weight_sch[weight_index];
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = separate_sch->nAttrs;
        itemSize[0] = separate_sch->item_size;
        itemPerBlk[0] = separate_sch->item_per_blk;
        nAttrs[1] = weight_sch[weight_index].nAttrs;
        itemSize[1] = weight_sch[weight_index].item_size;
        itemPerBlk[1] = weight_sch[weight_index].item_per_blk;
        
        int32_t iTableID = iTable[weight_index];
        uint32_t curAttrID = attrID[weight_index];
        ATTR_TYPE attr_type = in_sch[iTableID].attrType[curAttrID];
        uint32_t offset = in_sch[iTableID].attrOffset[curAttrID];
        uint32_t attrSize = in_sch[iTableID].attrSize[curAttrID];
        uint32_t pos_offset = itemSize[0] - sizeof(uint32_t);
        
        char* iPBlock[2];
        char* iPItem[2];
        for (uint32_t j = 0; j < 2; ++j)
            iPBlock[j] = buffer + j * B;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        char oItem[plain_len];
        char* oSItem = oItem + sizeof(char);
        uint32_t oItemContentSize = itemSize[1] - sizeof(char);
        
        // generate input table with join degrees
        est_real_item_num = 0;
        for (uint32_t i = 0; i < separate_n_blocks; ++i) {
            for (uint32_t j = 0; j < 2; ++j) {
                readBlock(separate_conn[j], i, iPBlock[j]);
                iPItem[j] = iPBlock[j] + META_BLOCK_SIZE;
            }
            for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                if (*(iPItem[0]) == 'r') {
                    assert(*(iPItem[1]) == 'r');
                    memcpy(oItem, iPItem[0], itemSize[0]);
                    
                    if(attr_type == INTEGER) {
                        int32_t join_key;
                        memcpy(&join_key, oItem + offset, sizeof(int32_t));
                        join_key += (int32_t)(bRange[weight_index]);
                        memcpy(oItem + offset, &join_key, sizeof(int32_t));
                    }
                    else if(attr_type == DOUBLE) {
                        double join_key;
                        memcpy(&join_key, oItem + offset, sizeof(double));
                        join_key += bRange[weight_index];
                        memcpy(oItem + offset, &join_key, sizeof(double));
                    }
                    
                    int32_t rpos, spos;
                    memcpy(&rpos, iPItem[0] + pos_offset, sizeof(int32_t));
                    memcpy(&spos, iPItem[1] + pos_offset, sizeof(int32_t));
                    int32_t alpha = spos - rpos;
                    assert(alpha >= 0);
                    memcpy(oItem + itemSize[0], &alpha, sizeof(int32_t));
                    est_real_item_num += alpha;
                    
                    ++oRItemCnt;
                    memcpy(oPItem, oItem, itemSize[1]);
                    ++oItemCnt;
                    if (oItemCnt < itemPerBlk[1])
                        oPItem += itemSize[1];
                    else {
                        memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                        memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                        insertBlock(weight_conn[weight_index], oBlockID, oBlock);
                        ++oBlockID;
                        oItemCnt = 0;
                        oRItemCnt = 0;
                        oPItem = oBlock + META_BLOCK_SIZE;
                    }
                }
                if (j < itemPerBlk[0] - 1) {
                    iPItem[0] += itemSize[0];
                    iPItem[1] += itemSize[0];
                }
            }
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
            insertBlock(weight_conn[weight_index], oBlockID, oBlock);
            ++oBlockID;
        }
        weight_n_blocks[weight_index] = oBlockID;
        
        /***********************/
        comm_size += 2 * separate_n_blocks + weight_n_blocks[weight_index];
        /***********************/
    }
    
    void expandBandTable(ServerConnector** weight_conn, Schema* weight_sch, const uint32_t* weight_n_blocks, ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t expand_index) {
        // input table with f values
        std::string fval_dst = dst_prefix + "_fval";
        ServerConnector* fval_conn = new FileSimulator(server_host, fval_dst.c_str(), true);
        Schema* fval_sch = new Schema; // schema of input table with f values
        uint32_t fval_n_blocks; // # of blocks in input table with f values
        
        genBandFValues(weight_conn, weight_sch, weight_n_blocks, fval_conn, fval_sch, fval_n_blocks, expand_index);
        
        // obliviously sort the input table with f values
        uint32_t nAttrs = fval_sch->nAttrs;
        CMP fval_cmp {1, {nAttrs - 1}, {0}};
        ObliviousSort* oblisort = new ObliviousSort(fval_conn, fval_n_blocks, fval_sch, &fval_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
        
        // if (result size < input size), then resize the input file with f values
        uint32_t est_real_block_num = (uint32_t)ceil((double)est_real_item_num / fval_sch->item_per_blk);
        if (est_real_block_num < fval_n_blocks) {
            resizeFile(fval_conn, est_real_block_num);
            comm_size += 2 * est_real_block_num;
        }
        
        //obliviously distribute the items
        distItems(fval_conn, fval_sch, fval_n_blocks);
        
        //fill in missing items
        fillBandItems(fval_conn, fval_sch, fval_n_blocks, expand_conn, expand_sch, expand_n_blocks, expand_index);
        
        delete fval_conn;
        delete fval_sch;
    }
    
    void genBandFValues(ServerConnector** weight_conn, Schema* weight_sch, const uint32_t* weight_n_blocks, ServerConnector* fval_conn, Schema* fval_sch, uint32_t & fval_n_blocks, const uint32_t weight_index) {
        // initialize the schema of input table with f values
        const uint32_t append_num = 1; // f value
        fval_sch->nAttrs = weight_sch[weight_index].nAttrs + append_num;
        uint32_t fval_index = weight_sch[weight_index].nAttrs;
        memcpy(fval_sch->attrType, weight_sch[weight_index].attrType, sizeof(ATTR_TYPE) * fval_index);
        memcpy(fval_sch->attrSize, weight_sch[weight_index].attrSize, sizeof(uint32_t) * fval_index);
        memcpy(fval_sch->attrOffset, weight_sch[weight_index].attrOffset, sizeof(uint32_t) * fval_index);
        uint32_t fval_item_size = weight_sch[weight_index].item_size;
        for (uint32_t j = 0; j < append_num; ++j) {
            fval_sch->attrType[fval_index] = INTEGER;
            fval_sch->attrSize[fval_index] = sizeof(uint32_t);
            fval_sch->attrOffset[fval_index] = fval_item_size;
            fval_item_size += sizeof(uint32_t);
            ++fval_index;
        }
        fval_sch->item_size = fval_item_size;
        fval_sch->item_per_blk = (plain_len - META_BLOCK_SIZE) / fval_item_size;
        
        // [0]: weight_sch[weight_index]; [1]: fval_sch;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = weight_sch[weight_index].nAttrs;
        nAttrs[1] = fval_sch->nAttrs;
        itemSize[0] = weight_sch[weight_index].item_size;
        itemSize[1] = fval_sch->item_size;
        itemPerBlk[0] = weight_sch[weight_index].item_per_blk;
        itemPerBlk[1] = fval_sch->item_per_blk;
        
        // current f value
        int32_t fval = 0;
        uint32_t alpha_offset = weight_sch[weight_index].attrOffset[nAttrs[0] - 1];
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = weight_n_blocks[weight_index];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(weight_conn[weight_index], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        int32_t alpha;
                        memcpy(&alpha, iPItem + alpha_offset, sizeof(int32_t));
                        memcpy(oPItem, iPItem, itemSize[0]);
                        memcpy(oPItem + itemSize[0], &fval, sizeof(int32_t));
                        fval += alpha;
                        
                        if(alpha == 0) *oPItem = 'd';
                        else ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(fval_conn, oBlockID, oBlock);
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
        
        uint32_t est_real_block_num = (uint32_t)ceil((double)est_real_item_num / itemPerBlk[1]);
        while (oBlockID < est_real_block_num || oItemCnt > 0) {
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
            insertBlock(fval_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
            oPItem = oBlock + META_BLOCK_SIZE;
        }
        fval_n_blocks = oBlockID;
        
        /***********************/
        comm_size += weight_n_blocks[weight_index] + fval_n_blocks;
        /***********************/
    }
    
    void fillBandItems(ServerConnector* fval_conn, Schema* fval_sch, const uint32_t fval_n_blocks, ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t expand_index) {
        // initialize the schema of the expanded table
        expand_sch[expand_index].nAttrs = fval_sch->nAttrs - 2;
        memcpy(expand_sch[expand_index].attrType, fval_sch->attrType, sizeof(ATTR_TYPE) * expand_sch[expand_index].nAttrs);
        memcpy(expand_sch[expand_index].attrSize, fval_sch->attrSize, sizeof(uint32_t) * expand_sch[expand_index].nAttrs);
        memcpy(expand_sch[expand_index].attrOffset, fval_sch->attrOffset, sizeof(uint32_t) * expand_sch[expand_index].nAttrs);
        expand_sch[expand_index].item_size = fval_sch->item_size - 2 * sizeof(uint32_t);
        expand_sch[expand_index].item_per_blk = (plain_len - META_BLOCK_SIZE) / expand_sch[expand_index].item_size;
        
        // [0]: fval_sch; [1]: expand_sch[expand_index];
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = fval_sch->nAttrs;
        nAttrs[1] = expand_sch[expand_index].nAttrs;
        itemSize[0] = fval_sch->item_size;
        itemSize[1] = expand_sch[expand_index].item_size;
        itemPerBlk[0] = fval_sch->item_per_blk;
        itemPerBlk[1] = expand_sch[expand_index].item_per_blk;
        
        uint32_t pos_offset = fval_sch->attrOffset[nAttrs[0] - 3];
        
        char pre_item[itemSize[0]];
        int32_t item_index = -1;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = fval_n_blocks;
        uint32_t real_item_num = 0;
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(fval_conn, iLeft, iRight, buffer);
            
            char* iPBlock = buffer;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock + META_BLOCK_SIZE;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (real_item_num < est_real_item_num) {
                        if (*iPItem == 'r') {
                            memcpy(pre_item, iPItem, itemSize[0]);
                            item_index = 1;
                        }
                        else ++item_index;
                        memcpy(oPItem, pre_item, itemSize[1]);
                        
                        int32_t position;
                        memcpy(&position, pre_item + pos_offset, sizeof(int32_t));
                        position += item_index;
                        memcpy(oPItem + pos_offset, &position, sizeof(int32_t));

                        ++oRItemCnt;
                        ++real_item_num;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(expand_conn[expand_index], oBlockID, oBlock);
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
            insertBlock(expand_conn[expand_index], oBlockID, oBlock);
            ++oBlockID;
        }
        expand_n_blocks[expand_index] = oBlockID;
        /***********************/
        comm_size += fval_n_blocks + expand_n_blocks[expand_index];
        /***********************/
    }
    
    void stitchBandTable(ServerConnector** expand_conn, const Schema* expand_sch, const uint32_t* expand_n_blocks) {
        uint32_t nAttrs[3];
        uint32_t itemSize[3];
        uint32_t itemPerBlk[3];
        uint32_t attrSize[2];
        uint32_t offset[2];
        for (uint32_t j = 0; j < 2; ++j) {
            uint32_t tableID = iTable[j];
            nAttrs[j] = expand_sch[j].nAttrs;
            itemSize[j] = expand_sch[j].item_size;
            itemPerBlk[j] = expand_sch[j].item_per_blk;
            uint32_t curAttrID = attrID[j];
            attrSize[j] = in_sch[tableID].attrSize[curAttrID];
            offset[j] = in_sch[tableID].attrOffset[curAttrID];
        }
        nAttrs[2] = out_sch->nAttrs;
        itemSize[2] = out_sch->item_size;
        itemPerBlk[2] = out_sch->item_per_blk;
        
        // obliviously sort the expanded input table based on position and item_id
        CMP expand_cmp {2, {nAttrs[1] - 1, nAttrs[1] - 2}, {0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(expand_conn[1], expand_n_blocks[1], &(expand_sch[1]), &expand_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
        
        char* iBlock[2];
        char* iPItem[2];
        uint32_t iBlockCnt[2];
        uint32_t iItemCnt[2];
        memset(iBlockCnt, 0, sizeof(uint32_t) * 2);
        memset(iItemCnt, 0, sizeof(uint32_t) * 2);
        for (uint32_t j = 0; j < 2; ++j) {
            iBlock[j] = buffer + j * B;
            readBlock(expand_conn[j], iBlockCnt[j], iBlock[j]);
            ++iBlockCnt[j];
            iPItem[j] = iBlock[j] + META_BLOCK_SIZE;
        }
        
        char oBlock[B];
        char oItem[itemSize[2]];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        // stitch the expanded tables
        bool iEnd = false;
        while (!iEnd) {
            if (*iPItem[0] == 'r') {
                assert(*iPItem[1] == 'r');
                writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
                
                for (uint32_t j = 0; j < 2; ++j) {
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
        while (oItemCnt != 0) {
            writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
        }
        out_tot_block = oBlockID;
        
        /***********************/
        for (uint32_t j = 0; j < 2; ++j)
            comm_size += expand_n_blocks[j];
        comm_size += out_tot_block;
        /***********************/
    }
    
    /***********************/
    //The followings are for equi-joins
    void augmentTables(ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks) {
        // unify input tables
        std::string union_dst = dst_prefix + "_union";
        ServerConnector* union_conn = new FileSimulator(server_host, union_dst.c_str(), true);
        Schema* union_sch = new Schema; // schema of unified input table
        uint32_t union_n_blocks; // # of blocks in unified input table
        
        unifyTables(union_conn, union_sch, union_n_blocks);
        
        // compute the join degrees
        for (uint32_t i = 0; i < nTables; ++i)
            fillDimensions(union_conn, union_sch, union_n_blocks, i);
        
        // generate the input tables with join degrees
        genWTables(union_conn, union_sch, union_n_blocks, weight_conn, weight_sch, weight_n_blocks);
        
        delete union_conn;
        delete union_sch;
    }
    
    void unifyTables(ServerConnector* union_conn, Schema* union_sch, uint32_t & union_n_blocks) {
        // [0]: in_sch[iTable[0]]; [1]: in_sch[iTable[1]]; [2]: union_sch;
        uint32_t nAttrs[3];
        uint32_t itemSize[3];
        uint32_t itemPerBlk[3];
        for (uint32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
            nAttrs[index] = in_sch[tableID].nAttrs;
            itemSize[index] = in_sch[tableID].item_size;
            itemPerBlk[index] = in_sch[tableID].item_per_blk;
        }
        
        // initialize the schema of the unified table
        const uint32_t append_num = 3; // tid, alpha1, alpha2
        union_sch->nAttrs = nAttrs[0] + nAttrs[1] + append_num - 2; // flag column, joined column
        uint32_t union_index = nAttrs[0];
        memcpy(union_sch->attrType, in_sch[iTable[0]].attrType, sizeof(ATTR_TYPE) * union_index);
        memcpy(union_sch->attrSize, in_sch[iTable[0]].attrSize, sizeof(uint32_t) * union_index);
        memcpy(union_sch->attrOffset, in_sch[iTable[0]].attrOffset, sizeof(uint32_t) * union_index);
        uint32_t union_item_size = itemSize[0];
        for (uint32_t j = 1; j < nAttrs[1]; ++j) {
            if (j == attrID[1]) continue;
            union_sch->attrType[union_index] = in_sch[iTable[1]].attrType[j];
            union_sch->attrSize[union_index] = in_sch[iTable[1]].attrSize[j];
            union_sch->attrOffset[union_index] = union_item_size;
            union_item_size += union_sch->attrSize[union_index];
            ++union_index;
        }
        for (uint32_t j = 0; j < append_num; ++j) {
            union_sch->attrType[union_index] = INTEGER;
            union_sch->attrSize[union_index] = sizeof(uint32_t);
            union_sch->attrOffset[union_index] = union_item_size;
            union_item_size += sizeof(uint32_t);
            ++union_index;
        }
        union_sch->item_size = union_item_size;
        union_sch->item_per_blk = (plain_len - META_BLOCK_SIZE) / union_item_size;
        
        nAttrs[2] = union_sch->nAttrs;
        itemSize[2] = union_sch->item_size;
        itemPerBlk[2] = union_sch->item_per_blk;
        
        // construct the unified table
        uint32_t attr_size[nTables];
        uint32_t offset[nTables];
        uint32_t append_offset = itemSize[2] - append_num * sizeof(uint32_t);
        for (uint32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
            uint32_t curAttrID = attrID[index];
            attr_size[index] = in_sch[tableID].attrSize[curAttrID];
            offset[index] = in_sch[tableID].attrOffset[curAttrID];
        }
        
        uint32_t iItemSegPos = offset[1] + attr_size[1];
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = offset[1] - sizeof(char);
        iItemSegSize[1] = itemSize[1] - iItemSegPos;
        uint32_t oItemSegPos = itemSize[0] + iItemSegSize[0];
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        for (uint32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
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
                    for (uint32_t j = 0; j < itemPerBlk[index]; ++j) {
                        if (*iPItem == 'r') {
                            if (index == 0) {
                                memcpy(oPItem, iPItem, itemSize[0]);
                            }
                            else if (index == 1) {
                                oPItem[0] = 'r';
                                memcpy(oPItem + offset[0], iPItem + offset[1], attr_size[1]);
                                memcpy(oPItem + itemSize[0], iPItem + sizeof(char), iItemSegSize[0]);
                                memcpy(oPItem + oItemSegPos, iPItem + iItemSegPos, iItemSegSize[1]);
                            }
                            memcpy(oPItem + append_offset, &index, sizeof(uint32_t));
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
                        if (j < itemPerBlk[index] - 1)
                            iPItem += itemSize[index];
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
        
        /***********************/
        for (uint32_t index = 0; index < nTables; ++index) {
            uint32_t tableID = iTable[index];
            comm_size += nBlocks[tableID];
        }
        comm_size += union_n_blocks;
        /***********************/
        
        // obliviously sort the unified table
        CMP union_cmp {2, {attrID[0], nAttrs[2] - append_num}, {0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(union_conn, union_n_blocks, union_sch, &union_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
    }
    
    int joinCompare(char* attr0, ATTR_TYPE attrType0, uint32_t attrSize0, char* attr1, ATTR_TYPE attrType1, uint32_t attrSize1, double* bRange = NULL) {
        if (attrType0 == CHAR) {
            char tmp_char_0 = *attr0;
            char tmp_char_1 = *attr1;
            if (tmp_char_0 < tmp_char_1) return -1;
            else if (tmp_char_0 > tmp_char_1) return 1;
            return 0;
        }
        else if (attrType0 == INTEGER || attrType0 == DOUBLE) {
            double tmp_double_0;
            if (attrType0 == INTEGER) {
                int32_t tmp_int_0;
                memcpy(&tmp_int_0, attr0, sizeof(int32_t));
                tmp_double_0 = tmp_int_0;
            }
            else memcpy(&tmp_double_0, attr0, sizeof(double));
            
            double tmp_double_1;
            if (attrType1 == INTEGER) {
                int32_t tmp_int_1;
                memcpy(&tmp_int_1, attr1, sizeof(int32_t));
                tmp_double_1 = tmp_int_1;
            }
            else memcpy(&tmp_double_1, attr1, sizeof(double));
            
            if (bRange != NULL) {
                if (tmp_double_0 + bRange[1] < tmp_double_1 - 1e-6) return -1;
                else if (tmp_double_0 - bRange[0] > tmp_double_1 + 1e-6) return 1;
            }
            else {
                if (tmp_double_0 < tmp_double_1 - 1e-6) return -1;
                else if (tmp_double_0 > tmp_double_1 + 1e-6) return 1;
            }
            return 0;
        }
        else if (attrType0 == STRING || attrType0 == TINYTEXT) {
            uint32_t attrLen0 = attrSize0;
            if (attrType0 == TINYTEXT)
                attrLen0 = std::min((uint32_t)strlen(attr0), attrSize0);
            uint32_t attrLen1 = attrSize1;
            if (attrType1 == TINYTEXT)
                attrLen1 = std::min((uint32_t)strlen(attr1), attrSize1);
            
            int res = strncmp(attr0, attr1, std::min(attrLen0, attrLen1));
            if (res == 0) {
                if (attrLen0 < attrLen1) res = -1;
                else if (attrLen0 > attrLen1) res = 1;
            }
            return res;
        }
        return 0;
    }
    
    void fillDimensions(ServerConnector* union_conn, const Schema* union_sch, const uint32_t union_n_blocks, const uint32_t scan_order) {
        uint32_t nAttrs = union_sch->nAttrs;
        uint32_t itemSize = union_sch->item_size;
        uint32_t itemPerBlk = union_sch->item_per_blk;
        ATTR_TYPE attrType = union_sch->attrType[attrID[0]];
        uint32_t offset = union_sch->attrOffset[attrID[0]];
        uint32_t attrSize = union_sch->attrSize[attrID[0]];
        
        // [0]: tid, [1]: alpha1, [2]: alpha2
        const uint32_t append_num = 3;
        uint32_t append_offset[append_num];
        for (uint32_t i = 0, j = nAttrs - append_num; i < append_num; ++i, ++j)
            append_offset[i] = union_sch->attrOffset[j];
        
        // compute the group running sums
        bool begin = true;
        int32_t alpha[2];
        uint32_t pre_tid;
        uint32_t pre_attr_len = 0;
        char pre_attr[attrSize];
        
        uint32_t iLeft, iRight;
        uint32_t remain_block = union_n_blocks;
        est_real_item_num = 0;
        if (scan_order == 0) iLeft = 0;
        else iRight = union_n_blocks;
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            if (scan_order == 0) iRight = iLeft + pickup_block;
            else iLeft = iRight - pickup_block;
            readBlock(union_conn, iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            if (scan_order != 0) iPBlock += (pickup_block - 1) * B + (itemPerBlk - 1) * itemSize;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk; ++j) {
                    if (*iPItem == 'r') {
                        uint32_t cur_tid;
                        memcpy(&cur_tid, iPItem + append_offset[0], sizeof(uint32_t));
                        assert(cur_tid == 0 || cur_tid == 1);
                        char* cur_attr = iPItem + offset;
                        uint32_t cur_attr_len = 0;
                        if (attrType == TINYTEXT)
                            cur_attr_len = std::min((uint32_t)strlen(cur_attr), attrSize);
                        else cur_attr_len = attrSize;
                        
                        bool equal;
                        if (!begin) {
                            int cmpres = joinCompare(pre_attr, attrType, pre_attr_len, cur_attr, attrType, cur_attr_len);
                            equal = (cmpres == 0);
                        }
                        else {
                            equal = false;
                            begin = false;
                        }
                        if (scan_order == 0) {
                            if (equal) {
                                if (cur_tid == 0) {
                                    ++alpha[0];
                                    alpha[1] = 0;
                                }
                                else if (cur_tid == 1) {
                                    if (pre_tid == 0) alpha[1] = 1;
                                    else ++alpha[1];
                                }
                            }
                            else {
                                if (cur_tid == 0) {
                                    alpha[0] = 1;
                                    alpha[1] = 0;
                                }
                                else if (cur_tid == 1) {
                                    alpha[0] = 0;
                                    alpha[1] = 1;
                                }
                            }
                            memcpy(iPItem + append_offset[1], alpha, 2 * sizeof(uint32_t));
                        }
                        else {
                            if (equal) memcpy(iPItem + append_offset[1], alpha, 2 * sizeof(uint32_t));
                            else {
                                memcpy(alpha, iPItem + append_offset[1], 2 * sizeof(uint32_t));
                                est_real_item_num += alpha[0] * alpha[1];
                            }
                        }
                        pre_tid = cur_tid;
                        pre_attr_len = cur_attr_len;
                        memcpy(pre_attr, cur_attr, cur_attr_len);
                    }
                    if (j < itemPerBlk - 1) {
                        if (scan_order == 0) iPItem += itemSize;
                        else iPItem -= itemSize;
                    }
                }
                if (i < pickup_block - 1) {
                    if (scan_order == 0) iPBlock += B;
                    else iPBlock -= B;
                }
            }
            updateBlock(union_conn, iLeft, iRight, buffer);
            if (scan_order == 0) iLeft = iRight;
            else iRight = iLeft;
        }
        /***********************/
        if (scan_order == 1) {
            printf("\nThe number of estimated real items in result table: %u\n\n", est_real_item_num);
        }
        comm_size += 2 * union_n_blocks;
        /***********************/
    }
    
    void genWTables(ServerConnector* union_conn, const Schema* union_sch, const uint32_t union_n_blocks, ServerConnector** weight_conn, Schema* weight_sch, uint32_t* weight_n_blocks) {
        // obliviously sort the unified table with join degrees
        CMP union_cmp {2, {union_sch->nAttrs - 3, attrID[0]}, {0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(union_conn, union_n_blocks, union_sch, &union_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
        
        // initialize the schema of input tables with join degrees
        const uint32_t append_num = 2; // [0]: alpha1, [1]: alpha2
        for (uint32_t weight_index = 0; weight_index < nTables; ++weight_index) {
            uint32_t tableID = iTable[weight_index];
            uint32_t attr_index = in_sch[tableID].nAttrs;
            weight_sch[weight_index].nAttrs = attr_index + append_num;
            memcpy(weight_sch[weight_index].attrType, in_sch[tableID].attrType, sizeof(ATTR_TYPE) * attr_index);
            memcpy(weight_sch[weight_index].attrSize, in_sch[tableID].attrSize, sizeof(uint32_t) * attr_index);
            memcpy(weight_sch[weight_index].attrOffset, in_sch[tableID].attrOffset, sizeof(uint32_t) * attr_index);
            uint32_t weight_item_size = in_sch[tableID].item_size;
            for (uint32_t i = 0; i < append_num; ++i) {
                weight_sch[weight_index].attrType[attr_index] = INTEGER;
                weight_sch[weight_index].attrSize[attr_index] = sizeof(uint32_t);
                weight_sch[weight_index].attrOffset[attr_index] = weight_item_size;
                weight_item_size += sizeof(uint32_t);
                ++attr_index;
            }
            weight_sch[weight_index].item_size = weight_item_size;
            weight_sch[weight_index].item_per_blk = (plain_len - META_BLOCK_SIZE) / weight_item_size;
        }
        
        // [0]: weight_sch[0]; [1]: weight_sch[1]; [2]: union_sch;
        uint32_t nAttrs[3];
        uint32_t itemSize[3];
        uint32_t itemPerBlk[3];
        for (uint32_t weight_index = 0; weight_index < nTables; ++weight_index) {
            nAttrs[weight_index] = weight_sch[weight_index].nAttrs;
            itemSize[weight_index] = weight_sch[weight_index].item_size;
            itemPerBlk[weight_index] = weight_sch[weight_index].item_per_blk;
        }
        nAttrs[2] = union_sch->nAttrs;
        itemSize[2] = union_sch->item_size;
        itemPerBlk[2] = union_sch->item_per_blk;
        
        // [0]: in_sch[iTable[0]]; [1]: in_sch[iTable[1]];
        uint32_t inItemSize[2];
        uint32_t offset[2];
        uint32_t attrSize[2];
        for (uint32_t i = 0; i < nTables; ++i) {
            int32_t iTableID = iTable[i];
            inItemSize[i] = in_sch[iTableID].item_size;
            offset[i] = in_sch[iTableID].attrOffset[attrID[i]];
            attrSize[i] = in_sch[iTableID].attrSize[attrID[i]];
        }
        
        // [0]: tid, [1]: alpha1, [2]: alpha2
        uint32_t append_offset[3];
        for (uint32_t i = 0, j = nAttrs[2] - 3; i < 3; ++i, ++j)
            append_offset[i] = union_sch->attrOffset[j];
        
        uint32_t iItemSegSize[2];
        iItemSegSize[0] = offset[1] - sizeof(char);
        uint32_t iItemSegPos = inItemSize[0] + iItemSegSize[0];
        uint32_t oItemSegPos = offset[1] + attrSize[0];
        iItemSegSize[1] = inItemSize[1] - oItemSegPos;
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        char oItem[plain_len];
        char* oSItem = oItem + sizeof(char);
        uint32_t oItemContentSize[2];
        for (uint32_t i = 0; i < nTables; ++i)
            oItemContentSize[i] = itemSize[i] - sizeof(char);
        
        // generate the input tables with join degrees
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = union_n_blocks;
        uint32_t weight_index = 0;
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
                for (uint32_t j = 0; j < itemPerBlk[2]; ++j) {
                    if (*iPItem == 'r') {
                        uint32_t tid;
                        memcpy(&tid, iPItem + append_offset[0], sizeof(uint32_t));
                        assert(tid == 0 || tid == 1);
                        if (tid != weight_index) {
                            assert(weight_index == 0 && tid == 1);
                            if (oItemCnt > 0) {
                                oItem[0] = 'd';
                                while (oItemCnt < itemPerBlk[weight_index]) {
                                    std::string rnd_str = generate_random_block(oItemContentSize[weight_index]);
                                    memcpy(oSItem, rnd_str.c_str(), oItemContentSize[weight_index]);
                                    memcpy(oPItem, oItem, itemSize[weight_index]);
                                    ++oItemCnt;
                                    if (oItemCnt < itemPerBlk[weight_index])
                                        oPItem += itemSize[weight_index];
                                }
                                memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                                memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                                insertBlock(weight_conn[weight_index], oBlockID, oBlock);
                                ++oBlockID;
                            }
                            weight_n_blocks[weight_index] = oBlockID;
                            weight_index = 1;
                            oPItem = oBlock + META_BLOCK_SIZE;
                            oBlockID = 0;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                        }
                        assert(tid == weight_index);
                        if (weight_index == 0)
                            memcpy(oItem, iPItem, inItemSize[0]);
                        else if (weight_index == 1) {
                            oItem[0] = 'r';
                            memcpy(oItem + sizeof(char), iPItem + inItemSize[0], iItemSegSize[0]);
                            memcpy(oItem + offset[1], iPItem + offset[0], attrSize[0]);
                            memcpy(oItem + oItemSegPos, iPItem + iItemSegPos, iItemSegSize[1]);
                        }
                        memcpy(oItem + inItemSize[weight_index], iPItem + append_offset[1], 2 * sizeof(uint32_t));
                        ++oRItemCnt;
                        memcpy(oPItem, oItem, itemSize[weight_index]);
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[weight_index])
                            oPItem += itemSize[weight_index];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(weight_conn[weight_index], oBlockID, oBlock);
                            ++oBlockID;
                            oItemCnt = 0;
                            oRItemCnt = 0;
                            oPItem = oBlock + META_BLOCK_SIZE;
                        }
                    }
                    if (j < itemPerBlk[2] - 1)
                        iPItem += itemSize[2];
                }
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            iLeft = iRight;
        }
        if (oItemCnt > 0) {
            oItem[0] = 'd';
            while (oItemCnt < itemPerBlk[weight_index]) {
                std::string rnd_str = generate_random_block(oItemContentSize[weight_index]);
                memcpy(oSItem, rnd_str.c_str(), oItemContentSize[weight_index]);
                memcpy(oPItem, oItem, itemSize[weight_index]);
                ++oItemCnt;
                if (oItemCnt < itemPerBlk[weight_index])
                    oPItem += itemSize[weight_index];
            }
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(weight_conn[weight_index], oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
        weight_n_blocks[weight_index] = oBlockID;
        
        /***********************/
        comm_size += union_n_blocks;
        for (uint32_t index = 0; index < nTables; ++index)
            comm_size += weight_n_blocks[index];
        /***********************/
    }
    
    void expandTable(ServerConnector** weight_conn, Schema* weight_sch, const uint32_t* weight_n_blocks, ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t expand_index) {
        // input table with f values
        std::string fval_dst = dst_prefix + "_fval";
        ServerConnector* fval_conn = new FileSimulator(server_host, fval_dst.c_str(), true);
        Schema* fval_sch = new Schema; // schema of input table with f values
        uint32_t fval_n_blocks; // # of blocks in input table with f values
        
        genFValues(weight_conn, weight_sch, weight_n_blocks, fval_conn, fval_sch, fval_n_blocks, expand_index);
        
        // obliviously sort the input table with f values
        uint32_t nAttrs = fval_sch->nAttrs;
        CMP fval_cmp {1, {nAttrs - 1}, {0}};
        ObliviousSort* oblisort = new ObliviousSort(fval_conn, fval_n_blocks, fval_sch, &fval_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
        
        // if (result size < input size), then resize the input file with f values
        uint32_t est_real_block_num = (uint32_t)ceil((double)est_real_item_num / fval_sch->item_per_blk);
        if (est_real_block_num < fval_n_blocks) {
            resizeFile(fval_conn, est_real_block_num);
            comm_size += 2 * est_real_block_num;
        }
        
        //obliviously distribute the items
        distItems(fval_conn, fval_sch, fval_n_blocks);
        
        //fill in missing items
        fillItems(fval_conn, fval_sch, fval_n_blocks, expand_conn, expand_sch, expand_n_blocks, expand_index);
        
        delete fval_conn;
        delete fval_sch;
    }
    
    void genFValues(ServerConnector** weight_conn, Schema* weight_sch, const uint32_t* weight_n_blocks, ServerConnector* fval_conn, Schema* fval_sch, uint32_t & fval_n_blocks, const uint32_t weight_index) {
        // initialize the schema of input table with f values
        const uint32_t append_num = 1; // f value
        fval_sch->nAttrs = weight_sch[weight_index].nAttrs + append_num;
        uint32_t fval_index = weight_sch[weight_index].nAttrs;
        memcpy(fval_sch->attrType, weight_sch[weight_index].attrType, sizeof(ATTR_TYPE) * fval_index);
        memcpy(fval_sch->attrSize, weight_sch[weight_index].attrSize, sizeof(uint32_t) * fval_index);
        memcpy(fval_sch->attrOffset, weight_sch[weight_index].attrOffset, sizeof(uint32_t) * fval_index);
        uint32_t fval_item_size = weight_sch[weight_index].item_size;
        for (uint32_t j = 0; j < append_num; ++j) {
            fval_sch->attrType[fval_index] = INTEGER;
            fval_sch->attrSize[fval_index] = sizeof(uint32_t);
            fval_sch->attrOffset[fval_index] = fval_item_size;
            fval_item_size += sizeof(uint32_t);
            ++fval_index;
        }
        fval_sch->item_size = fval_item_size;
        fval_sch->item_per_blk = (plain_len - META_BLOCK_SIZE) / fval_item_size;
        
        // [0]: weight_sch[weight_index]; [1]: fval_sch;
        uint32_t nAttrs[2];
        uint32_t itemSize[2];
        uint32_t itemPerBlk[2];
        nAttrs[0] = weight_sch[weight_index].nAttrs;
        nAttrs[1] = fval_sch->nAttrs;
        itemSize[0] = weight_sch[weight_index].item_size;
        itemSize[1] = fval_sch->item_size;
        itemPerBlk[0] = weight_sch[weight_index].item_per_blk;
        itemPerBlk[1] = fval_sch->item_per_blk;
        
        // current f value
        int32_t fval = 0;
        uint32_t alpha_offset = weight_sch[weight_index].attrOffset[nAttrs[0] - 1 - weight_index];
        
        char oBlock[B];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = weight_n_blocks[weight_index];
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(weight_conn[weight_index], iLeft, iRight, buffer);
            
            char* iPBlock = buffer + META_BLOCK_SIZE;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                char* iPItem = iPBlock;
                for (uint32_t j = 0; j < itemPerBlk[0]; ++j) {
                    if (*iPItem == 'r') {
                        int32_t alpha;
                        memcpy(&alpha, iPItem + alpha_offset, sizeof(int32_t));
                        memcpy(oPItem, iPItem, itemSize[0]);
                        memcpy(oPItem + itemSize[0], &fval, sizeof(int32_t));
                        fval += alpha;
                        
                        if(alpha == 0) *oPItem = 'd';
                        else ++oRItemCnt;
                        ++oItemCnt;
                        if (oItemCnt < itemPerBlk[1])
                            oPItem += itemSize[1];
                        else {
                            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
                            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
                            insertBlock(fval_conn, oBlockID, oBlock);
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
        
        uint32_t est_real_block_num = (uint32_t)ceil((double)est_real_item_num / itemPerBlk[1]);
        while (oBlockID < est_real_block_num || oItemCnt > 0) {
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
            insertBlock(fval_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
            oPItem = oBlock + META_BLOCK_SIZE;
        }
        fval_n_blocks = oBlockID;
        
        /***********************/
        comm_size += weight_n_blocks[weight_index] + fval_n_blocks;
        /***********************/
    }
    
    void distItems(ServerConnector* fval_conn, Schema* fval_sch, const uint32_t fval_n_blocks) {
        char oItem[fval_sch->item_size];
        oItem[0] = 'd';
        char* oSItem = oItem + sizeof(char);
        uint32_t oItemContentSize = fval_sch->item_size - sizeof(char);
        std::string rnd_str = generate_random_block(oItemContentSize);
        memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
        
        uint32_t fval_offset = fval_sch->item_size - sizeof(int32_t);
        int32_t max_level = (int32_t)ceil(log2(est_real_item_num));
        for (int32_t level = max_level - 1; level >= 0; --level) {
            int32_t skip = 1 << level;
            int32_t dst_tot_item_id = est_real_item_num - 1;
            int32_t src_tot_item_id = dst_tot_item_id - skip;
            int32_t dst_block_id = dst_tot_item_id / fval_sch->item_per_blk;
            int32_t dst_item_id = dst_tot_item_id % fval_sch->item_per_blk;
            int32_t src_block_id = src_tot_item_id / fval_sch->item_per_blk;
            int32_t src_item_id = src_tot_item_id % fval_sch->item_per_blk;
            char* src_block = buffer;
            char* dst_block = buffer + B;
            readBlock(fval_conn, src_block_id, src_block);
            /***********************/
            ++comm_size;
            /***********************/
            if (src_block_id == dst_block_id)
                memcpy(dst_block, src_block, B);
            else {
                readBlock(fval_conn, dst_block_id, dst_block);
                /***********************/
                ++comm_size;
                /***********************/
            }
            while (src_tot_item_id >= 0) {
                char* src_item = src_block + META_BLOCK_SIZE + src_item_id * fval_sch->item_size;
                char* dst_item = dst_block + META_BLOCK_SIZE + dst_item_id * fval_sch->item_size;
                if (*src_item == 'r') {
                    int32_t fval;
                    memcpy(&fval, src_item + fval_offset, sizeof(int32_t));
                    if (dst_tot_item_id <= fval) {
                        memcpy(dst_item, src_item, fval_sch->item_size);
                        memcpy(src_item, oItem, fval_sch->item_size);
                        if (src_block_id == dst_block_id) {
                            src_item = dst_block + META_BLOCK_SIZE + src_item_id * fval_sch->item_size;
                            memcpy(src_item, oItem, fval_sch->item_size);
                        }
                    }
                }
                
                --src_tot_item_id;
                --dst_tot_item_id;
                if (src_item_id > 0) --src_item_id;
                else {
                    if (src_block_id != dst_block_id) {
                        updateBlock(fval_conn, src_block_id, src_block);
                        /***********************/
                        ++comm_size;
                        /***********************/
                    }
                    --src_block_id;
                    src_item_id = fval_sch->item_per_blk - 1;
                    if (src_block_id >= 0) {
                        readBlock(fval_conn, src_block_id, src_block);
                        /***********************/
                        ++comm_size;
                        /***********************/
                    }
                    else {
                        updateBlock(fval_conn, dst_block_id, dst_block);
                        /***********************/
                        ++comm_size;
                        /***********************/
                        break;
                    }
                }
                if (dst_item_id > 0) --dst_item_id;
                else {
                    updateBlock(fval_conn, dst_block_id, dst_block);
                    /***********************/
                    ++comm_size;
                    /***********************/
                    --dst_block_id;
                    dst_item_id = fval_sch->item_per_blk - 1;
                    assert(dst_block_id >= 0);
                    if (src_block_id == dst_block_id)
                        memcpy(dst_block, src_block, B);
                    else {
                        readBlock(fval_conn, dst_block_id, dst_block);
                        /***********************/
                        ++comm_size;
                        /***********************/
                    }
                }
            }
        }
    }
    
    void fillItems(ServerConnector* fval_conn, Schema* fval_sch, const uint32_t fval_n_blocks, ServerConnector** expand_conn, Schema* expand_sch, uint32_t* expand_n_blocks, const uint32_t expand_index) {
        uint32_t nAttrs = fval_sch->nAttrs;
        uint32_t itemSize = fval_sch->item_size;
        uint32_t itemPerBlk = fval_sch->item_per_blk;
        ATTR_TYPE attrType = fval_sch->attrType[attrID[expand_index]];
        uint32_t offset = fval_sch->attrOffset[attrID[expand_index]];
        uint32_t attrSize = fval_sch->attrSize[attrID[expand_index]];
        
        // initialize the schema of the expanded table
        expand_sch[expand_index].nAttrs = nAttrs;
        memcpy(expand_sch[expand_index].attrType, fval_sch->attrType, sizeof(ATTR_TYPE) * nAttrs);
        memcpy(expand_sch[expand_index].attrSize, fval_sch->attrSize, sizeof(uint32_t) * nAttrs);
        memcpy(expand_sch[expand_index].attrOffset, fval_sch->attrOffset, sizeof(uint32_t) * nAttrs);
        expand_sch[expand_index].item_size = itemSize;
        expand_sch[expand_index].item_per_blk = itemPerBlk;
        
        // [0]: alpha1, [1]: alpha2, [2]: ii
        const uint32_t append_num = 3;
        uint32_t append_offset[append_num];
        for (uint32_t i = 0, j = nAttrs - append_num; i < append_num; ++i, ++j)
            append_offset[i] = fval_sch->attrOffset[j];
        
        int32_t alpha[2];
        char pre_item[itemSize];
        uint32_t pre_attr_len = 0;
        char pre_attr[attrSize];
        int32_t attr_index = -1;
        
        uint32_t iLeft = 0;
        uint32_t iRight;
        uint32_t remain_block = fval_n_blocks;
        uint32_t real_item_num = 0;
        while (remain_block > 0) {
            uint32_t pickup_block;
            if (remain_block >= two_m_block)
                pickup_block = two_m_block;
            else pickup_block = remain_block;
            remain_block -= pickup_block;
            iRight = iLeft + pickup_block;
            readBlock(fval_conn, iLeft, iRight, buffer);
            
            char* iPBlock = buffer;
            for (uint32_t i = 0; i < pickup_block; ++i) {
                uint32_t oRItemCnt = 0;
                char* iPItem = iPBlock + META_BLOCK_SIZE;
                for (uint32_t j = 0; j < itemPerBlk; ++j) {
                    if (real_item_num < est_real_item_num) {
                        if (*iPItem == 'r')
                            memcpy(pre_item, iPItem, itemSize);
                        else memcpy(iPItem, pre_item, itemSize);
                        
                        if (expand_index == 1) {
                            char* cur_attr = iPItem + offset;
                            uint32_t cur_attr_len = 0;
                            if (attrType == TINYTEXT)
                                cur_attr_len = std::min((uint32_t)strlen(cur_attr), attrSize);
                            else cur_attr_len = attrSize;
                            
                            bool equal;
                            if (real_item_num == 0)
                                equal = false;
                            else {
                                int cmpres = joinCompare(pre_attr, attrType, pre_attr_len, cur_attr, attrType, cur_attr_len);
                                equal = (cmpres == 0);
                            }
                            if (equal) ++attr_index;
                            else attr_index = 0;
                            
                            memcpy(alpha, iPItem + append_offset[0], 2 * sizeof(int32_t));
                            int32_t ii = attr_index / alpha[0] + (attr_index % alpha[0]) * alpha[1];
                            memcpy(iPItem + append_offset[2], &ii, sizeof(int32_t));
                            
                            pre_attr_len = cur_attr_len;
                            memcpy(pre_attr, cur_attr, cur_attr_len);
                        }
                        ++oRItemCnt;
                        ++real_item_num;
                    }
                    if (j < itemPerBlk - 1)
                        iPItem += itemSize;
                }
                memcpy(iPBlock + sizeof(uint32_t), &oRItemCnt, sizeof(uint32_t));
                if (i < pickup_block - 1)
                    iPBlock += B;
            }
            insertBlock(expand_conn[expand_index], iLeft, iRight, buffer);
            iLeft = iRight;
        }
        expand_n_blocks[expand_index] = fval_n_blocks;
        /***********************/
        comm_size += fval_n_blocks + expand_n_blocks[expand_index];
        /***********************/
    }
    
    void writeItem(char type, char** iPItem, char* oBlock, char* oItem, uint32_t& oItemCnt, uint32_t& oRItemCnt, uint32_t& oBlockID, bool equi = true) {
        char* oSItem = oItem + sizeof(char);
        /*********************
        if (type == 'r') printf("real result\n");
        else printf("dummy result\n");
        *********************/
        if (type == 'r') {
            if (nProjCols == 0) {
                uint32_t attrSize[nTables];
                uint32_t offset[nTables];
                uint32_t itemSize[nTables];
                for (uint32_t j = 0; j < nTables; ++j) {
                    int32_t tableID = iTable[j];
                    if (j > 0) {
                        uint32_t curAttrID = attrID[j];
                        attrSize[j] = in_sch[tableID].attrSize[curAttrID];
                        offset[j] = in_sch[tableID].attrOffset[curAttrID];
                    }
                    itemSize[j] = in_sch[tableID].item_size;
                }
                
                uint32_t iItemSegPos[nTables];
                uint32_t iItemSegSize[nTables][2];
                uint32_t oItemSegPos[nTables][2];
                iItemSegSize[0][0] = itemSize[0] - sizeof(char);
                iItemSegSize[0][1] = 0;
                oItemSegPos[0][0] = 0;
                oItemSegPos[0][1] = iItemSegSize[0][0];
                if (equi) {
                    for (uint32_t j = 1; j < nTables; ++j) {
                        iItemSegPos[j] = offset[j] + attrSize[j];
                        iItemSegSize[j][0] = offset[j] - sizeof(char);
                        iItemSegSize[j][1] = itemSize[j] - iItemSegPos[j];
                        oItemSegPos[j][0] = oItemSegPos[j - 1][1] + iItemSegSize[j - 1][1];
                        oItemSegPos[j][1] = oItemSegPos[j][0] + iItemSegSize[j][0];
                    }
                    oItem[0] = 'r';
                    memcpy(oSItem, iPItem[0] + sizeof(char), iItemSegSize[0][0]);
                    for (uint32_t j = 1; j < nTables; ++j) {
                        memcpy(oSItem + oItemSegPos[j][0], iPItem[j] + sizeof(char), iItemSegSize[j][0]);
                        memcpy(oSItem + oItemSegPos[j][1], iPItem[j] + iItemSegPos[j], iItemSegSize[j][1]);
                    }
                    ++oRItemCnt;
                }
                else {
                    for (uint32_t j = 1; j < nTables; ++j) {
                        iItemSegSize[j][0] = itemSize[j] - sizeof(char);
                        oItemSegPos[j][0] = oItemSegPos[j - 1][0] + iItemSegSize[j - 1][0];
                    }
                    oItem[0] = 'r';
                    memcpy(oSItem, iPItem[0] + sizeof(char), iItemSegSize[0][0]);
                    for (uint32_t j = 1; j < nTables; ++j)
                        memcpy(oSItem + oItemSegPos[j][0], iPItem[j] + sizeof(char), iItemSegSize[j][0]);
                    ++oRItemCnt;
                }
            }
            else {
                oItem[0] = 'r';
                uint32_t oItemPos = 0;
                for (uint32_t j = 0; j < nProjCols; ++j) {
                    uint32_t iTableIndex = projID[0][j];
                    int32_t tableID = iTable[iTableIndex];
                    uint32_t curAttrID = projID[1][j];
                    uint32_t attrSize = in_sch[tableID].attrSize[curAttrID];
                    uint32_t offset = in_sch[tableID].attrOffset[curAttrID];
                    memcpy(oSItem + oItemPos, iPItem[iTableIndex] + offset, attrSize);
                    oItemPos += attrSize;
                }
                ++oRItemCnt;
            }
        }
        else {
            oItem[0] = 'd';
            uint32_t oItemContentSize = out_sch->item_size - sizeof(char);
            std::string rnd_str = generate_random_block(oItemContentSize);
            memcpy(oSItem, rnd_str.c_str(), oItemContentSize);
        }
        
        char* oPItem = oBlock + META_BLOCK_SIZE + oItemCnt * out_sch->item_size;
        memcpy(oPItem, oItem, out_sch->item_size);
        ++oItemCnt;
        if (oItemCnt >= out_sch->item_per_blk) {
            memcpy(oBlock, &oBlockID, sizeof(uint32_t));
            memcpy(&(oBlock[sizeof(uint32_t)]), &oRItemCnt, sizeof(uint32_t));
            insertBlock(out_conn, oBlockID, oBlock);
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
    }
    
    void stitchTable(ServerConnector** expand_conn, const Schema* expand_sch, const uint32_t* expand_n_blocks) {
        uint32_t nAttrs[3];
        uint32_t itemSize[3];
        uint32_t itemPerBlk[3];
        uint32_t attrSize[2];
        uint32_t offset[2];
        for (uint32_t j = 0; j < 2; ++j) {
            uint32_t tableID = iTable[j];
            nAttrs[j] = expand_sch[j].nAttrs;
            itemSize[j] = expand_sch[j].item_size;
            itemPerBlk[j] = expand_sch[j].item_per_blk;
            uint32_t curAttrID = attrID[j];
            attrSize[j] = in_sch[tableID].attrSize[curAttrID];
            offset[j] = in_sch[tableID].attrOffset[curAttrID];
        }
        nAttrs[2] = out_sch->nAttrs;
        itemSize[2] = out_sch->item_size;
        itemPerBlk[2] = out_sch->item_per_blk;
        
        // obliviously sort the expanded input table based on join column and ii
        CMP expand_cmp {2, {attrID[1], nAttrs[1] - 1}, {0, 0}};
        ObliviousSort* oblisort = new ObliviousSort(expand_conn[1], expand_n_blocks[1], &(expand_sch[1]), &expand_cmp);
        comm_size += oblisort->BatcherSort();
        delete oblisort;
        
        char* iBlock[2];
        char* iPItem[2];
        uint32_t iBlockCnt[2];
        uint32_t iItemCnt[2];
        memset(iBlockCnt, 0, sizeof(uint32_t) * 2);
        memset(iItemCnt, 0, sizeof(uint32_t) * 2);
        for (uint32_t j = 0; j < 2; ++j) {
            iBlock[j] = buffer + j * B;
            readBlock(expand_conn[j], iBlockCnt[j], iBlock[j]);
            ++iBlockCnt[j];
            iPItem[j] = iBlock[j] + META_BLOCK_SIZE;
        }
        
        char oBlock[B];
        char oItem[itemSize[2]];
        char* oPItem = oBlock + META_BLOCK_SIZE;
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        // stitch the expanded tables
        bool iEnd = false;
        while (!iEnd) {
            if (*iPItem[0] == 'r') {
                assert(*iPItem[1] == 'r');
                writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                
                for (uint32_t j = 0; j < 2; ++j) {
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
        while (oItemCnt != 0) {
            writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        }
        out_tot_block = oBlockID;
        
        /***********************/
        for (uint32_t j = 0; j < 2; ++j)
            comm_size += expand_n_blocks[j];
        comm_size += out_tot_block;
        /***********************/
    }
    
    // # of input tables
    uint32_t nInputs;
    // schema of input tables
    Schema* in_sch = NULL;
    
    // server connectors of input tables
    ServerConnector** in_conn = NULL;
    // # of blocks in input tables
    uint32_t* nBlocks = NULL;
    // # of items in input tables
    uint32_t* iItemNum = NULL;
    
    // join query info
    uint32_t nTables;
    uint32_t* iTable = NULL;
    int32_t* attrID = NULL;
    int32_t* iParent = NULL;
    int32_t* pAttrID = NULL;
    // range for band join
    double* bRange = NULL;
    // # of columns in the result after projection
    // 0 means all columns remain in the result
    uint32_t nProjCols = 0;
    uint32_t* projID[2] = {NULL, NULL};
    
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
    
    // # of estimated join records
    uint32_t est_real_item_num;
    
    /***********************/
    size_t comm_size = 0;
    /***********************/
};

#endif //__OBLI_DATABASE_JOIN_H__
