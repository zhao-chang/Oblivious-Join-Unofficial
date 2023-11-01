#ifndef __ODS_INDEX_JOIN_H__
#define __ODS_INDEX_JOIN_H__

#include "PrivIndexJoin.h"
#include "ODSBTree.h"

using namespace std;

template<class T>
class ODSIndexJoin : public PrivIndexJoin {
public:
    // load
    ODSIndexJoin(const uint32_t mode, const uint32_t n_inputs, const std::string* input_file, const std::string& dst, const std::string& prefix) : PrivIndexJoin (n_inputs, dst, prefix) {
        /***********************/
        btree_time = 0.0;
        /***********************/
        
        std::string oij_name = prefix + std::string("_oij_type2_info.txt");
        FILE* fp = fopen(oij_name.c_str(), "w");
        btree = new ODSBTree <T>* [this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            fprintf(fp, "%s\n", input_file[i].c_str());
            generateInputSchema(mode, input_file[i], &(this->in_sch[i]), true);
            
            std::string cur_dst = dst + "_input_" + std::to_string(i);
            std::string cur_prefix = prefix + "_input_" + std::to_string(i);
            btree[i] = new ODSBTree <T> (mode, input_file[i], &(this->in_sch[i]), cur_dst, cur_prefix);
        }
        fclose(fp);
    }
    
    // run
    ODSIndexJoin(const uint32_t mode, const uint32_t n_inputs, const std::string& prefix) : PrivIndexJoin (n_inputs, prefix) {
        /***********************/
        btree_time = 0.0;
        /***********************/
        
        std::string oij_name = prefix + std::string("_oij_type2_info.txt");
        FILE* fp = fopen(oij_name.c_str(), "r");
        btree = new ODSBTree <T>* [this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            char input_file_name[100];
            fscanf(fp, "%s", input_file_name);
            std::string input_file = std::string(input_file_name);
            generateInputSchema(mode, input_file, &(this->in_sch[i]), true);
            
            std::string cur_prefix = prefix + "_input_" + std::to_string(i);
            btree[i] = new ODSBTree <T> (&(this->in_sch[i]), cur_prefix);
        }
        fclose(fp);
    }
    
    ~ODSIndexJoin() {
        if (btree != NULL) {
            for (uint32_t i = 0; i < this->nInputs; ++i)
                delete btree[i];
            delete[] btree;
            btree = NULL;
        }
    }
    
    // oblivious sort-merge join
    void ObliSMJ(const uint32_t* attr_id) {
        this->attrID = new int32_t[this->nInputs];
        memcpy(this->attrID, attr_id, this->nInputs * sizeof(uint32_t));
        std::string out_dst = this->dst_prefix + "_output";
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            out_dst += "_" + std::to_string(this->attrID[i]);
        }
        this->out_conn = new MongoConnector(server_host, out_dst.c_str(), true);
        this->out_sch = new Schema;
        generateOutputSchema(this->in_sch, this->nInputs, this->attrID, this->out_sch);
        this->out_cmp = new CMP({1, {this->attrID[0]}, {0}});
        
        // perform oblivious join
        performSMJ();
        
        // perform oblivious filter
        uint32_t access_num = 0;
        filterOutput(this->out_conn, this->out_tot_block, this->out_sch, this->out_cmp, access_num);
        this->comm_size += access_num;
        
        if (this->attrID != NULL) {
            delete[] this->attrID;
            this->attrID = NULL;
        }
        if (this->out_conn != NULL) {
            delete this->out_conn;
            this->out_conn = NULL;
        }
        if (this->out_sch != NULL) {
            delete this->out_sch;
            this->out_sch = NULL;
        }
        if (this->out_cmp != NULL) {
            delete this->out_cmp;
            this->out_cmp = NULL;
        }
    }
    
    // oblivious block-nested loop join
    void ObliBNLJ(const uint32_t* attr_id) {
        this->attrID = new int32_t[this->nInputs];
        memcpy(this->attrID, attr_id, this->nInputs * sizeof(uint32_t));
        std::string out_dst = this->dst_prefix + "_output";
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            out_dst += "_" + std::to_string(this->attrID[i]);
        }
        this->out_conn = new MongoConnector(server_host, out_dst.c_str(), true);
        this->out_sch = new Schema;
        generateOutputSchema(this->in_sch, this->nInputs, this->attrID, this->out_sch);
        this->out_cmp = new CMP({1, {this->attrID[0]}, {0}});
        
        // perform oblivious join
        performBNLJ();
        
        // perform oblivious filter
        uint32_t access_num = 0;
        filterOutput(this->out_conn, this->out_tot_block, this->out_sch, this->out_cmp, access_num);
        this->comm_size += access_num;
        
        if (this->attrID != NULL) {
            delete[] this->attrID;
            this->attrID = NULL;
        }
        if (this->out_conn != NULL) {
            delete this->out_conn;
            this->out_conn = NULL;
        }
        if (this->out_sch != NULL) {
            delete this->out_sch;
            this->out_sch = NULL;
        }
        if (this->out_cmp != NULL) {
            delete this->out_cmp;
            this->out_cmp = NULL;
        }
    }
    
    /***********************/
    size_t getServerSize() const {
        size_t server_size = 0;
        for (uint32_t i = 0; i < nInputs; ++i)
            server_size += btree[i]->getServerSize();
        return server_size;
    }
    
    size_t getClientSize() const {
        size_t client_size = 0;
        for (uint32_t i = 0; i < nInputs; ++i)
            client_size += btree[i]->getClientSize();
        return client_size;
    }
    
    size_t getCommSize() const {
        size_t total_comm_size = 0;
        for (uint32_t i = 0; i < nInputs; ++i) {
            printf("B-Tree %u: comm size is %lu\n", i, btree[i]->getCommSize());
            total_comm_size += btree[i]->getCommSize();
        }
        total_comm_size += this->comm_size * (uint64_t)B;
        return total_comm_size;
    }
    
    void resetCommSize() {
        this->comm_size = 0;
        for (uint32_t i = 0; i < nInputs; ++i)
            btree[i]->resetCommSize();
    }
    
    size_t getAccessCount() const {
        size_t access_count = 0;
        for (uint32_t i = 0; i < nInputs; ++i) {
            access_count += btree[i]->getAccessCount();
            printf("Access Count %u: %lu\n", i, btree[i]->getAccessCount());
        }
        return access_count;
    }
    
    double getReadTime() const {
        double read_time = 0.0;
        for (uint32_t i = 0; i < nInputs; ++i) {
            read_time += btree[i]->getReadTime();
            printf("Read Time %u: %lf\n", i, btree[i]->getReadTime());
        }
        return read_time;
    }
    
    double getWriteTime() const {
        double write_time = 0.0;
        for (uint32_t i = 0; i < nInputs; ++i) {
            write_time += btree[i]->getWriteTime();
            printf("Write Time %u: %lf\n", i, btree[i]->getWriteTime());
        }
        return write_time;
    }
    
    double getEncDecTime() const {
        double enc_dec_time = 0.0;
        for (uint32_t i = 0; i < nInputs; ++i) {
            enc_dec_time += btree[i]->getEncDecTime();
            printf("EncDec Time %u: %lf\n", i, btree[i]->getEncDecTime());
        }
        return enc_dec_time;
    }
    
    double getORAMTime() const {
        double oram_time = 0.0;
        for (uint32_t i = 0; i < nInputs; ++i) {
            oram_time += btree[i]->getORAMTime();
            printf("ORAM Time %u: %lf\n", i, btree[i]->getORAMTime());
        }
        return oram_time;
    }
    
    double getBTreeTime() const {
        return btree_time;
    }
    /***********************/
    
private:
    // for index based sort-merge join and block-nested loop join
    void getNextTuples(int32_t iIndex, bool* iNext, char* iLeaf, char* iData, int32_t* iDataBID, int32_t** iAllIndex, int32_t* iEndLevel, int32_t* iOutsourcedTreeHeight) {
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        assert(iIndex == -1 || iIndex > 0);
        if (iIndex == -1) {
            for (uint32_t i = 0; i < this->nInputs; ++i)
                for (uint32_t j = 0; j <= iOutsourcedTreeHeight[i]; ++j)
                    btree[i]->getDummyTuple();
        }
        else if (iIndex > 0) {
            for (uint32_t i = 0; i < iIndex; ++i)
                for (uint32_t j = 0; j <= iOutsourcedTreeHeight[i]; ++j)
                    btree[i]->getDummyTuple();
            btree[iIndex]->getNextTuple(iNext[iIndex], this->attrID[iIndex], iLeaf + iIndex * B, iData + iIndex * B, iDataBID[iIndex], iAllIndex[iIndex], iEndLevel[iIndex]);
            for (uint32_t i = iIndex + 1; i < this->nInputs; ++i)
                for (uint32_t j = 0; j <= iOutsourcedTreeHeight[i]; ++j)
                    btree[i]->getDummyTuple();
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
    }
    
    // index based sort-merge join
    void performSMJ() {
        ATTR_TYPE attrType[this->nInputs];
        uint32_t attrSize[this->nInputs];
        uint32_t offset[this->nInputs];
        uint32_t itemSize[this->nInputs];
        for (uint32_t j = 0; j < this->nInputs; ++j) {
            uint32_t curAttrID = this->attrID[j];
            attrType[j] = this->in_sch[j].attrType[curAttrID];
            attrSize[j] = this->in_sch[j].attrSize[curAttrID];
            offset[j] = this->in_sch[j].attrOffset[curAttrID];
            itemSize[j] = this->in_sch[j].item_size;
        }
        
        if (attrType[0] == CHAR) {
            for (uint32_t j = 1; j < this->nInputs; ++j)
                assert(attrType[j] == CHAR);
        }
        else if (attrType[0] == INTEGER || attrType[0] == DOUBLE) {
            for (uint32_t j = 1; j < this->nInputs; ++j)
                assert(attrType[j] == INTEGER || attrType[j] == DOUBLE);
        }
        else if (attrType[0] == STRING || attrType[0] == TINYTEXT) {
            for (uint32_t j = 1; j < this->nInputs; ++j)
                assert(attrType[j] == STRING || attrType[j] == TINYTEXT);
        }
        
        char iLeaf[this->nInputs][B];
        char iData[this->nInputs][B];
        int32_t iDataBID[this->nInputs];
        char iBackLeaf[this->nInputs][B];
        char iBackData[this->nInputs][B];
        int32_t iBackDataBID[this->nInputs];
        char iForwardLeaf[this->nInputs][B];
        char iForwardData[this->nInputs][B];
        int32_t iForwardDataBID[this->nInputs];
        
        int32_t iTreeHeight[this->nInputs];
        int32_t iOutsourcedTreeHeight[this->nInputs];
        int32_t* iAllIndex[this->nInputs];
        int32_t iEndLevel[this->nInputs];
        int32_t* iBackAllIndex[this->nInputs];
        int32_t iBackEndLevel[this->nInputs];
        int32_t* iForwardAllIndex[this->nInputs];
        int32_t iForwardEndLevel[this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            iTreeHeight[i] = btree[i]->getHeight(this->attrID[i]);
            iOutsourcedTreeHeight[i] = btree[i]->getOutsourcedHeight(this->attrID[i]);
            iAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
            iBackAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
            iForwardAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
        }
        
        bool iNext[this->nInputs];
        bool iForwardNext[this->nInputs];
        uint32_t iJoinNum[this->nInputs];
        uint32_t iCurJoinNum[this->nInputs];
        bool iFindJoinNum[this->nInputs];
        memset(iNext, true, this->nInputs * sizeof(bool));
        memset(iForwardNext, true, this->nInputs * sizeof(bool));
        memset(iJoinNum, 0, this->nInputs * sizeof(uint32_t));
        memset(iCurJoinNum, 0, this->nInputs * sizeof(uint32_t));
        memset(iFindJoinNum, false, this->nInputs * sizeof(bool));
        
        char* iItem[this->nInputs];
        char* iPItem[this->nInputs];
        char* iPAttr[this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i)
            iItem[i] = iData[i] + sizeof(uint32_t);
        
        char oBlock[B];
        char oItem[plain_len];
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        bool match = false;
        uint32_t curCmpNum = 0;
        uint32_t iStartIndex = 0;
        uint32_t iCurIndex = 0;
        
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            btree[i]->getFirstTuple(iNext[i], this->attrID[i], iLeaf[i], iData[i], iDataBID[i], iAllIndex[i], iEndLevel[i]);
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        while (true) {
            bool iEnd = false;
            bool iFirstLastMatch = false;
            bool iOtherLastMatch = false;
            uint32_t iIndex;
            for (iIndex = 0; iIndex < this->nInputs; ++iIndex) {
                if (!iNext[iIndex]) {
                    if (match) {
                        if (iIndex == 0) iFirstLastMatch = true;
                        else iOtherLastMatch = true;
                    }
                    else iEnd = true;
                    break;
                }
            }
            if (iEnd) break;
            
            int cmpres;
            if (iFirstLastMatch) {
                assert(iIndex == 0 && iIndex == iCurIndex);
                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                cmpres = 1;
            }
            else if (iOtherLastMatch) {
                assert(iIndex > 0 && iIndex == iCurIndex);
                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                cmpres = -1;
            }
            else {
                for (uint32_t i = 0; i < this->nInputs; ++i) {
                    iPItem[i] = iItem[i] + iAllIndex[i][iTreeHeight[i]] * itemSize[i];
                    iPAttr[i] = iPItem[i] + offset[i];
                }
                cmpres = this->joinCompare(iPAttr, attrType, attrSize, iStartIndex);
                if (cmpres == 0) this->writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                else this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
            }
            
            //Post-processing
            //TODO: get dummy block once or twice?
            if (cmpres < 0) {
                /*******************
                printf("cmpres < 0\n");
                if (match) printf("match -> match\n");
                else printf("not match -> not match\n");
                *******************/
                
                if (match) {
                    assert(iCurIndex > 0 && !iFindJoinNum[iCurIndex]);
                    ++curCmpNum;
                    iFindJoinNum[iCurIndex] = true;
                    if (iOtherLastMatch) {
                        iNext[iCurIndex] = true;
                        iForwardDataBID[iCurIndex] = -1;
                        memset(iForwardAllIndex[iCurIndex], -1, sizeof(int32_t) * (iTreeHeight[iCurIndex] + 1));
                        iForwardEndLevel[iCurIndex] = -1;
                        iForwardNext[iCurIndex] = false;
                    }
                    else {
                        memcpy(iForwardLeaf[iCurIndex], iLeaf[iCurIndex], B);
                        memcpy(iForwardData[iCurIndex], iData[iCurIndex], B);
                        iForwardDataBID[iCurIndex] = iDataBID[iCurIndex];
                        memcpy(iForwardAllIndex[iCurIndex], iAllIndex[iCurIndex], sizeof(int32_t) * (iTreeHeight[iCurIndex] + 1));
                        iForwardEndLevel[iCurIndex] = iEndLevel[iCurIndex];
                        iForwardNext[iCurIndex] = true;
                    }
                    memcpy(iLeaf[iCurIndex], iBackLeaf[iCurIndex], B);
                    memcpy(iData[iCurIndex], iBackData[iCurIndex], B);
                    iDataBID[iCurIndex] = iBackDataBID[iCurIndex];
                    memcpy(iAllIndex[iCurIndex], iBackAllIndex[iCurIndex], sizeof(int32_t) * (iTreeHeight[iCurIndex] + 1));
                    iEndLevel[iCurIndex] = iBackEndLevel[iCurIndex];
                    iCurJoinNum[iCurIndex] = 1;
                    for (uint32_t i = iCurIndex + 1; i < this->nInputs; ++i)
                        assert(iCurJoinNum[i] == 1);
                    
                    --iCurIndex;
                    getNextTuples(iCurIndex, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                }
                else {
                    assert(iStartIndex > 0);
                    --iStartIndex;
                    getNextTuples(iStartIndex, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                }
            }
            else if (cmpres > 0) {
                /*******************
                printf("cmpres > 0\n");
                if (match) printf("match -> not match\n");
                else printf("not match -> not match\n");
                *******************/
                
                if (match) {
                    assert(iCurIndex == 0 && !iFindJoinNum[iCurIndex]);
                    ++curCmpNum;
                    iFindJoinNum[iCurIndex] = true;
                    if (iFirstLastMatch) iNext[iCurIndex] = false;
                    
                    match = false;
                    for (uint32_t i = 1; i < this->nInputs; ++i) {
                        if (iForwardNext[i]) {
                            memcpy(iLeaf[i], iForwardLeaf[i], B);
                            memcpy(iData[i], iForwardData[i], B);
                            iDataBID[i] = iForwardDataBID[i];
                            memcpy(iAllIndex[i], iForwardAllIndex[i], sizeof(int32_t) * (iTreeHeight[i] + 1));
                            iEndLevel[i] = iForwardEndLevel[i];
                            iNext[i] = true;
                        }
                        else {
                            iDataBID[i] = -1;
                            memset(iAllIndex[i], -1, sizeof(int32_t) * (iTreeHeight[i] + 1));
                            iEndLevel[i] = -1;
                            iNext[i] = false;
                        }
                    }
                    
                    // padding to make "curCmpNum = iJoinNum[0] * ... * iJoinNum[nInputs - 1]
                    // + iJoinNum[0] + ... + iJoinNum[nInputs - 1]"
                    uint32_t totCmpNum = 1;
                    for (uint32_t i = 0; i < this->nInputs; ++i)
                        totCmpNum *= iJoinNum[i];
                    for (uint32_t i = 0; i < this->nInputs; ++i)
                        totCmpNum += iJoinNum[i];
                    
                    assert(curCmpNum <= totCmpNum);
                    while (curCmpNum < totCmpNum) {
                        getNextTuples(-1, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                        this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                        ++curCmpNum;
                    }
                }
                else {
                    getNextTuples(iStartIndex, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                }
            }
            else {
                /*******************
                printf("cmpres == 0\n");
                if (match) printf("match -> match\n");
                else printf("not match -> match\n");
                *******************/
                
                if (match) {
                    int32_t i;
                    for (i = this->nInputs - 1; i >= 0; --i) {
                        if (!iFindJoinNum[i] || iCurJoinNum[i] < iJoinNum[i]) {
                            ++curCmpNum;
                            if (!iFindJoinNum[i]) ++iJoinNum[i];
                            ++iCurJoinNum[i];
                            getNextTuples(i, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                            for (uint32_t j = i + 1; j < this->nInputs; ++j) {
                                assert(iFindJoinNum[j]);
                                memcpy(iLeaf[j], iBackLeaf[j], B);
                                memcpy(iData[j], iBackData[j], B);
                                iDataBID[j] = iBackDataBID[j];
                                memcpy(iAllIndex[j], iBackAllIndex[j], sizeof(int32_t) * (iTreeHeight[j] + 1));
                                iEndLevel[j] = iBackEndLevel[j];
                                iCurJoinNum[j] = 1;
                            }
                            iCurIndex = i;
                            break;
                        }
                    }
                    assert(i >= 0);
                }
                else {
                    match = true;
                    curCmpNum = 1;
                    for (uint32_t i = 0; i < this->nInputs; ++i) {
                        iJoinNum[i] = 1;
                        iCurJoinNum[i] = 1;
                    }
                    memset(iFindJoinNum, false, this->nInputs * sizeof(bool));
                    
                    for (uint32_t i = 0; i < this->nInputs; ++i) {
                        memcpy(iBackLeaf[i], iLeaf[i], B);
                        memcpy(iBackData[i], iData[i], B);
                        iBackDataBID[i] = iDataBID[i];
                        memcpy(iBackAllIndex[i], iAllIndex[i], sizeof(int32_t) * (iTreeHeight[i] + 1));
                        iBackEndLevel[i] = iEndLevel[i];
                    }
                    iCurIndex = this->nInputs - 1;
                    getNextTuples(this->nInputs - 1, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                }
            }
        }
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            while (iNext[i]) {
                getNextTuples(i, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
            }
        }
        while (oItemCnt != 0) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        }
        this->out_tot_block = oBlockID;
        
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            delete[] iAllIndex[i];
            delete[] iBackAllIndex[i];
            delete[] iForwardAllIndex[i];
        }
    }
    
    // index based block-nested loop join
    void performBNLJ() {
        ATTR_TYPE attrType[this->nInputs];
        uint32_t attrSize[this->nInputs];
        uint32_t offset[this->nInputs];
        uint32_t itemSize[this->nInputs];
        for (uint32_t j = 0; j < this->nInputs; ++j) {
            uint32_t curAttrID = this->attrID[j];
            attrType[j] = this->in_sch[j].attrType[curAttrID];
            attrSize[j] = this->in_sch[j].attrSize[curAttrID];
            offset[j] = this->in_sch[j].attrOffset[curAttrID];
            itemSize[j] = this->in_sch[j].item_size;
        }
        
        if (attrType[0] == CHAR) {
            for (uint32_t j = 1; j < this->nInputs; ++j)
                assert(attrType[j] == CHAR);
        }
        else if (attrType[0] == INTEGER || attrType[0] == DOUBLE) {
            for (uint32_t j = 1; j < this->nInputs; ++j)
                assert(attrType[j] == INTEGER || attrType[j] == DOUBLE);
        }
        else if (attrType[0] == STRING || attrType[0] == TINYTEXT) {
            for (uint32_t j = 1; j < this->nInputs; ++j)
                assert(attrType[j] == STRING || attrType[j] == TINYTEXT);
        }
        
        char iLeaf[this->nInputs][B];
        char iData[this->nInputs][B];
        int32_t iDataBID[this->nInputs];
        char iBackLeaf[this->nInputs][B];
        char iBackData[this->nInputs][B];
        int32_t iBackDataBID[this->nInputs];
        
        int32_t iTreeHeight[this->nInputs];
        int32_t iOutsourcedTreeHeight[this->nInputs];
        int32_t* iAllIndex[this->nInputs];
        int32_t iEndLevel[this->nInputs];
        int32_t* iBackAllIndex[this->nInputs];
        int32_t iBackEndLevel[this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            iTreeHeight[i] = btree[i]->getHeight(this->attrID[i]);
            iOutsourcedTreeHeight[i] = btree[i]->getOutsourcedHeight(this->attrID[i]);
            iAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
            iBackAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
        }
        
        bool iNext[this->nInputs];
        uint32_t iJoinNum[this->nInputs];
        uint32_t iCurJoinNum[this->nInputs];
        bool iFindJoinNum[this->nInputs];
        
        char* iItem[this->nInputs];
        char* iPItem[this->nInputs];
        char* iPAttr[this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i)
            iItem[i] = iData[i] + sizeof(uint32_t);
        
        char oBlock[B];
        char oItem[plain_len];
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        btree[0]->getFirstTuple(iNext[0], this->attrID[0], iLeaf[0], iData[0], iDataBID[0], iAllIndex[0], iEndLevel[0]);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        while (iNext[0]) {
            iPItem[0] = iItem[0] + iAllIndex[0][iTreeHeight[0]] * itemSize[0];
            iPAttr[0] = iPItem[0] + offset[0];
            /***********************/
            start = std::chrono::high_resolution_clock::now();
            /***********************/
            for (uint32_t i = 1; i < this->nInputs; ++i) {
                btree[i]->getFirstTuple(iNext[i], iPAttr[0], attrType[0], attrSize[0], this->attrID[i], iLeaf[i], iData[i], iDataBID[i], iAllIndex[i], iEndLevel[i]);
            }
            /***********************/
            end = std::chrono::high_resolution_clock::now();
            diff = end - start;
            btree_time += (double)(diff.count());
            /***********************/
            
            bool match = false;
            uint32_t curCmpNum = 0;
            uint32_t iStartIndex = 0;
            uint32_t iCurIndex = 0;
            memset(iJoinNum, 0, this->nInputs * sizeof(uint32_t));
            memset(iCurJoinNum, 0, this->nInputs * sizeof(uint32_t));
            memset(iFindJoinNum, false, this->nInputs * sizeof(bool));
            
            while (true) {
                bool iEnd = false;
                bool iFirstLastMatch = false;
                bool iOtherLastMatch = false;
                uint32_t iIndex;
                for (iIndex = 1; iIndex < this->nInputs; ++iIndex) {
                    if (!iNext[iIndex]) {
                        if (match) {
                            if (iIndex == 1) iFirstLastMatch = true;
                            else iOtherLastMatch = true;
                        }
                        else iEnd = true;
                        break;
                    }
                }
                if (iEnd) break;
                
                int cmpres;
                if (iFirstLastMatch || iOtherLastMatch) {
                    assert(iIndex >= 1 && iIndex == iCurIndex);
                    this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                    cmpres = -1;
                }
                else {
                    for (uint32_t i = 1; i < this->nInputs; ++i) {
                        iPItem[i] = iItem[i] + iAllIndex[i][iTreeHeight[i]] * itemSize[i];
                        iPAttr[i] = iPItem[i] + offset[i];
                    }
                    cmpres = this->joinCompare(iPAttr, attrType, attrSize, iStartIndex);
                    if (cmpres == 0) this->writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                    else this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                }
                
                //Post-processing
                if (cmpres < 0) {
                    if (match) {
                        ++curCmpNum;
                        assert(!iFindJoinNum[iCurIndex]);
                        iFindJoinNum[iCurIndex] = true;
                        if (iCurIndex > 1) {
                            /***********
                            printf("cmpres < 0 and iCurIndex > 1\n");
                            printf("match -> match\n");
                            ***********/
                            
                            if (iOtherLastMatch) iNext[iCurIndex] = true;
                            memcpy(iLeaf[iCurIndex], iBackLeaf[iCurIndex], B);
                            memcpy(iData[iCurIndex], iBackData[iCurIndex], B);
                            iDataBID[iCurIndex] = iBackDataBID[iCurIndex];
                            memcpy(iAllIndex[iCurIndex], iBackAllIndex[iCurIndex], sizeof(int32_t) * (iTreeHeight[iCurIndex] + 1));
                            iEndLevel[iCurIndex] = iBackEndLevel[iCurIndex];
                            iCurJoinNum[iCurIndex] = 1;
                            for (uint32_t i = iCurIndex + 1; i < this->nInputs; ++i)
                                assert(iCurJoinNum[i] == 1);
                            
                            --iCurIndex;
                            getNextTuples(iCurIndex, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                        }
                        else {
                            /***********
                            printf("cmpres < 0 and iCurIndex == 1\n");
                            printf("match -> not match\n");
                            ***********/
                            
                            match = false;
                            // padding to make "curCmpNum = iJoinNum[1] * ... * iJoinNum[nInputs - 1]
                            // + (nInputs - 1)"
                            uint32_t totCmpNum = 1;
                            for (uint32_t i = 1; i < this->nInputs; ++i)
                                totCmpNum *= iJoinNum[i];
                            totCmpNum += this->nInputs - 1;
                            
                            assert(curCmpNum <= totCmpNum);
                            while (curCmpNum < totCmpNum) {
                                getNextTuples(-1, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                                ++curCmpNum;
                            }
                            break;
                        }
                    }
                    else {
                        /***********
                        printf("cmpres < 0\n");
                        printf("not match -> not match\n");
                        ***********/
                        
                        //padding to (nInputs - 1) = 1 + (nInputs - 2)
                        for (uint32_t i = 0; i < this->nInputs - 2; ++i) {
                            getNextTuples(-1, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                        }
                        break;
                    }
                }
                else {
                    /***********
                    printf("cmpres == 0\n");
                    if (match) printf("match -> match\n");
                    else printf("not match -> match\n");
                    ***********/
                    
                    if (match) {
                        int32_t i;
                        for (i = this->nInputs - 1; i >= 1; --i) {
                            if (!iFindJoinNum[i] || iCurJoinNum[i] < iJoinNum[i]) {
                                ++curCmpNum;
                                if (!iFindJoinNum[i]) ++iJoinNum[i];
                                ++iCurJoinNum[i];
                                getNextTuples(i, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                                
                                for (uint32_t j = i + 1; j < this->nInputs; ++j) {
                                    assert(iFindJoinNum[j]);
                                    memcpy(iLeaf[j], iBackLeaf[j], B);
                                    memcpy(iData[j], iBackData[j], B);
                                    iDataBID[j] = iBackDataBID[j];
                                    memcpy(iAllIndex[j], iBackAllIndex[j], sizeof(int32_t) * (iTreeHeight[j] + 1));
                                    iEndLevel[j] = iBackEndLevel[j];
                                    iCurJoinNum[j] = 1;
                                }
                                iCurIndex = i;
                                break;
                            }
                        }
                        assert(i >= 1);
                    }
                    else {
                        match = true;
                        curCmpNum = 1;
                        for (uint32_t i = 1; i < this->nInputs; ++i) {
                            iJoinNum[i] = 1;
                            iCurJoinNum[i] = 1;
                        }
                        memset(iFindJoinNum, false, this->nInputs * sizeof(bool));
                        
                        for (uint32_t i = 1; i < this->nInputs; ++i) {
                            memcpy(iBackLeaf[i], iLeaf[i], B);
                            memcpy(iBackData[i], iData[i], B);
                            iBackDataBID[i] = iDataBID[i];
                            memcpy(iBackAllIndex[i], iAllIndex[i], sizeof(int32_t) * (iTreeHeight[i] + 1));
                            iBackEndLevel[i] = iEndLevel[i];
                        }
                        iCurIndex = this->nInputs - 1;
                        getNextTuples(this->nInputs - 1, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight);
                    }
                }
            }
            /***********************/
            start = std::chrono::high_resolution_clock::now();
            /***********************/
            btree[0]->getNextTuple(iNext[0], this->attrID[0], iLeaf[0], iData[0], iDataBID[0], iAllIndex[0], iEndLevel[0]);
            /***********************/
            end = std::chrono::high_resolution_clock::now();
            diff = end - start;
            btree_time += (double)(diff.count());
            /***********************/
        }
        while (oItemCnt != 0) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        }
        this->out_tot_block = oBlockID;
        
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            delete[] iAllIndex[i];
            delete[] iBackAllIndex[i];
        }
    }
    
    // private B-tree
    ODSBTree <T>** btree = NULL;
    
    /***********************/
    double btree_time;
    /***********************/
};

#endif //__ODS_INDEX_JOIN_H__
