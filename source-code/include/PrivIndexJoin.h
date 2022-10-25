#ifndef __PRIV_INDEX_JOIN_H__
#define __PRIV_INDEX_JOIN_H__

#include "Schema.h"
#include "PathORAM.h"
#include "ObliviousSort.h"
#include "App.h"

using namespace std;

class PrivIndexJoin {
public:
    // load
    PrivIndexJoin(const uint32_t mode, const std::string& dst, const std::string& prefix) {
        nInputs = getTableNum(mode);
        in_sch = new Schema[nInputs];
        dst_prefix = dst;
        
        std::string key_name = prefix + std::string("_key.txt");
        initKey(key_name.c_str(), true);
        
        /***********************/
        res_time = 0.0;
        /***********************/
        
        // write info
        std::string oij_name = prefix + std::string("_pij_info.txt");
        FILE* fp = fopen(oij_name.c_str(), "w");
        fprintf(fp, "%s\n", dst_prefix.c_str());
        fclose(fp);
    }
    
    // run
    PrivIndexJoin(const uint32_t mode, const std::string& prefix) {
        nInputs = getTableNum(mode);
        in_sch = new Schema[nInputs];
        
        // read info
        std::string oij_name = prefix + std::string("_pij_info.txt");
        FILE* fp = fopen(oij_name.c_str(), "r");
        char out_dst_prefix[100];
        fscanf(fp, "%s", out_dst_prefix);
        dst_prefix = std::string(out_dst_prefix);
        fclose(fp);
        
        std::string key_name = prefix + std::string("_key.txt");
        initKey(key_name.c_str(), false);
        
        /***********************/
        res_time = 0.0;
        /***********************/
    }
    
    ~PrivIndexJoin() {
        if (in_sch != NULL) {
            delete[] in_sch;
            in_sch = NULL;
        }
    }
    
    /***********************/
    double getResTime() const {
        return res_time;
    }
    /***********************/
    
protected:
    int joinCompare(char* iPAttr[], ATTR_TYPE attrType[], uint32_t attrSize[], double* bRange = NULL) {
        if (attrType[0] == CHAR) {
            char tmp_char_0 = *(iPAttr[0]);
            char tmp_char_1 = *(iPAttr[1]);
            if (tmp_char_0 < tmp_char_1) return -1;
            else if (tmp_char_0 > tmp_char_1) return 1;
            return 0;
        }
        else if (attrType[0] == INTEGER || attrType[0] == DOUBLE) {
            double tmp_double_0;
            if (attrType[0] == INTEGER) {
                int32_t tmp_int_0;
                memcpy(&tmp_int_0, iPAttr[0], sizeof(int32_t));
                tmp_double_0 = tmp_int_0;
            }
            else memcpy(&tmp_double_0, iPAttr[0], sizeof(double));
            
            double tmp_double_1;
            if (attrType[1] == INTEGER) {
                int32_t tmp_int_1;
                memcpy(&tmp_int_1, iPAttr[1], sizeof(int32_t));
                tmp_double_1 = tmp_int_1;
            }
            else memcpy(&tmp_double_1, iPAttr[1], sizeof(double));
            
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
        else if (attrType[0] == STRING || attrType[0] == TINYTEXT) {
            uint32_t attrLen[2];
            for (uint32_t i = 0; i < 2; ++i) {
                if (attrType[i] == TINYTEXT)
                    attrLen[i] = std::min((uint32_t)strlen(iPAttr[i]), attrSize[i]);
                else attrLen[i] = attrSize[i];
            }
            int res = strncmp(iPAttr[0], iPAttr[1], std::min(attrLen[0], attrLen[1]));
            if (res == 0) {
                if (attrLen[0] < attrLen[1]) res = -1;
                else if (attrLen[0] > attrLen[1]) res = 1;
            }
            return res;
        }
        return 0;
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
            /***********************/
            auto start = std::chrono::high_resolution_clock::now();
            /***********************/
            insertBlock(out_conn, oBlockID, oBlock);
            /***********************/
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
            res_time += (double)(diff.count());
            /***********************/
            ++oBlockID;
            oItemCnt = 0;
            oRItemCnt = 0;
        }
    }
    
    // # of input tables
    uint32_t nInputs;
    // schema of input tables
    Schema* in_sch = NULL;
    
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
    
    // output connector
    ServerConnector* out_conn = NULL;
    // prefix of output collection name
    std::string dst_prefix;
    // schema of output table
    Schema* out_sch = NULL;
    // cmp of output table (only for ORDER BY)
    CMP* out_cmp = NULL;
    // # of blocks in output table
    uint32_t out_tot_block;
    
    /***********************/
    size_t comm_size = 0;
    double res_time;
    /***********************/
};

#endif //__PRIV_INDEX_JOIN_H__
