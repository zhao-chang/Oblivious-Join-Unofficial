#ifndef __OBLI_MULTIWAY_JOIN_H__
#define __OBLI_MULTIWAY_JOIN_H__

#include "PrivIndexJoin.h"
#include "ObliMultiwayBTree.h"

using namespace std;

template<class T>
class ObliMultiwayJoin : public PrivIndexJoin {
public:
    // load
    ObliMultiwayJoin(const uint32_t mode, const std::string scale, const uint32_t outsourced_height, const std::string& dst, const std::string& prefix) : PrivIndexJoin (mode, dst, prefix) {
        /***********************/
        btree_time = 0.0;
        /***********************/
        
        btree = new ObliMultiwayBTree <T>* [this->nInputs];
        iItemNum = new uint32_t[this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            generateInputSchema(mode, i, &(this->in_sch[i]));
            getInputItemNum(mode, scale, i, iItemNum[i]);
            
            std::string cur_dst = dst + "_input_" + std::to_string(i);
            std::string cur_prefix = prefix + "_input_" + std::to_string(i);
            btree[i] = new ObliMultiwayBTree <T> (mode, scale, i, &(this->in_sch[i]), outsourced_height, cur_dst, cur_prefix);
        }
        initBufferPara(getDataBlockNum());
        
        std::string omj_name = prefix + std::string("_omj_type1_info.txt");
        FILE* fp = fopen(omj_name.c_str(), "wb");
        fwrite(iItemNum, sizeof(uint32_t), this->nInputs, fp);
        fclose(fp);
    }
    
    // run
    ObliMultiwayJoin(const uint32_t mode, const std::string& prefix) : PrivIndexJoin (mode, prefix) {
        /***********************/
        btree_time = 0.0;
        /***********************/
        
        iItemNum = new uint32_t[this->nInputs];
        std::string omj_name = prefix + std::string("_omj_type1_info.txt");
        FILE* fp = fopen(omj_name.c_str(), "rb");
        fread(iItemNum, sizeof(uint32_t), this->nInputs, fp);
        fclose(fp);
        
        btree = new ObliMultiwayBTree <T>*[this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            generateInputSchema(mode, i, &(this->in_sch[i]));
            
            std::string cur_prefix = prefix + "_input_" + std::to_string(i);
            btree[i] = new ObliMultiwayBTree <T> (&(this->in_sch[i]), cur_prefix);
        }
        initBufferPara(getDataBlockNum());
    }
    
    ~ObliMultiwayJoin() {
        if (btree != NULL) {
            for (uint32_t i = 0; i < this->nInputs; ++i)
                delete btree[i];
            delete[] btree;
            btree = NULL;
        }
        if (iItemNum != NULL) {
            delete[] iItemNum;
            iItemNum = NULL;
        }
    }
    
    // oblivious block-nested loop join
    void ObliBNLJ(const uint32_t n_tables, const uint32_t* table_id, const int32_t* parent_id, const int32_t* attr_id, const int32_t* parent_attr_id, const uint32_t n_proj_cols, const uint32_t proj_id[][MAX_COLS]) {
        this->nTables = n_tables;
        this->iTable = new uint32_t[this->nTables];
        this->attrID = new int32_t[this->nTables];
        this->iParent = new int32_t[this->nTables];
        this->pAttrID = new int32_t[this->nTables];
        memcpy(this->iTable, table_id, this->nTables * sizeof(uint32_t));
        memcpy(this->attrID, attr_id, this->nTables * sizeof(int32_t));
        memcpy(this->iParent, parent_id, this->nTables * sizeof(int32_t));
        memcpy(this->pAttrID, parent_attr_id, this->nTables * sizeof(int32_t));
        this->nProjCols = n_proj_cols;
        if (this->nProjCols > 0) {
            for (uint32_t i = 0; i < 2; ++i) {
                this->projID[i] = new uint32_t[this->nProjCols];
                memcpy(this->projID[i], proj_id[i], this->nProjCols * sizeof(uint32_t));
            }
        }
        
        std::string out_dst = this->dst_prefix + "_output";
        for (uint32_t i = 0; i < this->nTables; ++i)
            out_dst += "_" + std::to_string(this->iTable[i]);
        for (uint32_t i = 1; i < this->nTables; ++i) {
            out_dst += "_" + std::to_string(this->attrID[i]);
            out_dst += "_" + std::to_string(this->pAttrID[i]);
        }
        this->out_conn = new FileSimulator(server_host, out_dst.c_str(), true);
        this->out_sch = new Schema;
        if (this->nProjCols == 0) generateOutputSchema(this->in_sch, this->nTables, this->iTable, this->attrID, this->out_sch);
        else generateOutputSchema(this->in_sch, this->nTables, this->iTable, this->nProjCols, this->projID, this->out_sch);
        // TODO: support ORDER-BY operator
        if (this->nProjCols == 0) this->out_cmp = new CMP({1, {this->pAttrID[1]}, {0}});
        else this->out_cmp = new CMP({1, {1}, {0}});
        
        // perform oblivious join
        performBNLJ();
        
        // perform oblivious filter
        uint32_t access_num = 0;
        filterOutput(this->out_conn, this->out_tot_block, this->out_sch, this->out_cmp, access_num);
        this->comm_size += access_num;
        
        /***********************/
        // print join result
        if (this->nProjCols == 0) {
            printResult(this->out_conn, this->out_tot_block, this->out_sch, 1, &(this->pAttrID[1]));
        }
        else {
            int32_t printAttrID[this->nProjCols];
            for (int32_t i = 0; i < this->nProjCols; ++i)
                printAttrID[i] = i + 1;
            printResult(this->out_conn, this->out_tot_block, this->out_sch, this->nProjCols, printAttrID);
        }
        /***********************/
        
        if (this->iTable != NULL) {
            delete[] this->iTable;
            this->iTable = NULL;
        }
        if (this->attrID != NULL) {
            delete[] this->attrID;
            this->attrID = NULL;
        }
        if (this->iParent != NULL) {
            delete[] this->iParent;
            this->iParent = NULL;
        }
        if (this->pAttrID != NULL) {
            delete[] this->pAttrID;
            this->pAttrID = NULL;
        }
        if (this->nProjCols > 0) {
            for (uint32_t i = 0; i < 2; ++i) {
                if (this->projID[i] != NULL) {
                    delete[] (this->projID[i]);
                    this->projID[i] = NULL;
                }
            }
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
    uint32_t getDataBlockNum() const {
        uint32_t block_num = 0;
        for (uint32_t i = 0; i < nInputs; ++i)
            block_num += btree[i]->getDataBlockNum();
        return block_num;
    }
    
    size_t getServerSize() const {
        size_t server_size = 0;
        for (uint32_t i = 0; i < nInputs; ++i)
            server_size += btree[i]->getServerSize();
        return server_size;
    }
    
    size_t getClientSize() const {
        size_t client_size = oblisort_mem;
        printf("Additional trusted memory size is %lf MB\n", (double)oblisort_mem / (1 << 20));
        for (uint32_t i = 0; i < nInputs; ++i)
            client_size += btree[i]->getClientSize();
        return client_size;
    }
    
    size_t getCommSize() const {
        size_t total_comm_size = 0;
        for (uint32_t i = 0; i < nInputs; ++i) {
            size_t btree_comm_size = btree[i]->getCommSize();
            printf("B-Tree %u: comm size is %lf MB\n", i, (double)btree_comm_size / (1 << 20));
            total_comm_size += btree_comm_size;
        }
        total_comm_size += this->comm_size * (size_t)B;
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
    
	// private B-tree
	ObliMultiwayBTree <T>** btree = NULL;

private:
    void outputDummyRecord(int32_t iIndex, char** iPItem, char* oBlock, char* oItem, uint32_t& oItemCnt, uint32_t& oRItemCnt, uint32_t& oBlockID) {
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        for (uint32_t i = iIndex + 1; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            int32_t curAttrID = this->attrID[i];
            for (uint32_t j = 0; j < iOutsourcedTreeHeight[i]; ++j)
                btree[tableID]->getDummyIndexBlock(curAttrID);
            btree[tableID]->getDummyDataBlock();
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
    }
    
    void updatePreLeaf(int32_t iIndex) {
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        for (uint32_t i = 0; i < iIndex; ++i) {
            int32_t tableID = this->iTable[i];
            if (i == 0) btree[tableID]->getDummyDataBlock();
            else {
                int32_t curAttrID = this->attrID[i];
                for (uint32_t j = 0; j < iOutsourcedTreeHeight[i]; ++j)
                    btree[tableID]->getDummyIndexBlock(curAttrID);
                btree[tableID]->getDummyDataBlock();
            }
        }
        int32_t tableID = this->iTable[iIndex];
        btree[tableID]->updatePreLeaf(this->attrID[iIndex], iPreLeaf[iIndex], iPreLeafBID[iIndex]);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
    }
    
    // the boolean return value indicates whether block "iPreLeaf[iIndex]" must be updated
    bool disableCurrentTuple(int32_t iIndex) {
        bool updated = false;
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        for (uint32_t i = 0; i < iIndex; ++i) {
            int32_t tableID = this->iTable[i];
            if (i == 0) btree[tableID]->getDummyDataBlock();
            else {
                int32_t curAttrID = this->attrID[i];
                for (uint32_t j = 0; j < iOutsourcedTreeHeight[i]; ++j)
                    btree[tableID]->getDummyIndexBlock(curAttrID);
                btree[tableID]->getDummyDataBlock();
            }
        }
        int32_t tableID = this->iTable[iIndex];
        updated = btree[tableID]->disableCurrentTuple(this->attrID[iIndex], iAllBlock[iIndex], iAllBID[iIndex], iAllIndex[iIndex], iPreLeaf[iIndex], iPreLeafBID[iIndex]);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        return updated;
    }
    
    bool getNextTuples(int32_t iIndex) {
        bool iNext = false;
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        assert(iIndex == -1 || iIndex > 0);
        
        int32_t endIndex;
        if (iIndex == -1) endIndex = this->nTables;
        else endIndex = iIndex;
        
        for (uint32_t i = 0; i < endIndex; ++i) {
            int32_t tableID = this->iTable[i];
            if (i == 0) btree[tableID]->getDummyDataBlock();
            else {
                int32_t curAttrID = this->attrID[i];
                for (uint32_t j = 0; j < iOutsourcedTreeHeight[i]; ++j)
                    btree[tableID]->getDummyIndexBlock(curAttrID);
                btree[tableID]->getDummyDataBlock();
            }
        }
        if (iIndex > 0) {
            int32_t tableID = this->iTable[iIndex];
            btree[tableID]->getNextTuple(iNext, this->attrID[iIndex], iAllBlock[iIndex], iAllBID[iIndex], iAllIndex[iIndex], iAllLastEnabled[iIndex], iPreLeaf[iIndex], iPreLeafBID[iIndex]);
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        return iNext;
    }
    
    void join(int32_t index, bool& return_flag, int32_t& return_index, char** iPItem, char* oBlock, char* oItem, uint32_t& oItemCnt, uint32_t& oRItemCnt, uint32_t& oBlockID) {
        return_flag = false;
        
        int32_t tableID = this->iTable[index];
        int32_t parentID = this->iParent[index];
        int32_t pIndex = -1;
        for (int32_t j = 0; j < index; ++j) {
            if (parentID == this->iTable[j]) {
                pIndex = j;
                break;
            }
        }
        assert(pIndex >= 0);
        
        ATTR_TYPE attrType[2];
        uint32_t attrSize[2];
        uint32_t offset[2];
        uint32_t itemSize[2];
        for (uint32_t j = 0; j < 2; ++j) {
            int32_t curTableID;
            int32_t curAttrID;
            if (j == 0) {
                curTableID = parentID;
                curAttrID = this->pAttrID[index];
            }
            else if (j == 1) {
                curTableID = tableID;
                curAttrID = this->attrID[index];
            }
            attrType[j] = this->in_sch[curTableID].attrType[curAttrID];
            attrSize[j] = this->in_sch[curTableID].attrSize[curAttrID];
            offset[j] = this->in_sch[curTableID].attrOffset[curAttrID];
            itemSize[j] = this->in_sch[curTableID].item_size;
        }
        
        char* iItem[2];
        char* iPAttr[2];
        iItem[0] = iAllBlock[pIndex][iTreeHeight[pIndex]] + sizeof(uint32_t);
        iPAttr[0] = iItem[0] + iAllIndex[pIndex][iTreeHeight[pIndex]] * itemSize[0] + offset[0];
        
        bool iNext;
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        btree[tableID]->getFirstTuple(iNext, iPAttr[0], attrType[0], attrSize[0], this->attrID[index], iAllBlock[index], iAllBID[index], iAllIndex[index], iAllLastEnabled[index]);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        
        if (!iNext) {
            outputDummyRecord(index, iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
            ++iRetrievalNum;
        }
        else {
            iPreLeafBID[index] = -1;
            iItem[1] = iAllBlock[index][iTreeHeight[index]] + sizeof(uint32_t);
            iPItem[index] = iItem[1] + iAllIndex[index][iTreeHeight[index]] * itemSize[1];
            iPAttr[1] = iPItem[index] + offset[1];
            int match = this->joinCompare(iPAttr, attrType, attrSize);
            if (match != 0) {
                outputDummyRecord(index, iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                ++iRetrievalNum;
            }
            else {
                while (true) {
                    if (index + 1 == this->nTables) {
                        return_flag = true;
                        this->writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                        ++iRetrievalNum;
                        ++oResNum;
                    }
                    else {
                        bool next_flag = false;
                        int32_t next_index = -1;
                        join(index + 1, next_flag, next_index, iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                        if (index != next_index) {
                            return_flag = next_flag;
                            return_index = next_index;
                            return;
                        }
                        if (next_flag) return_flag = true;
                        else {
                            bool pre_updated = disableCurrentTuple(index);
                            outputDummyRecord(index, iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                            ++iRetrievalNum;
                            if (pre_updated) {
                                updatePreLeaf(index);
                                outputDummyRecord(index, iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                                ++iRetrievalNum;
                            }
                        }
                    }
                    if (btree[tableID]->hasNextTuple(this->attrID[index], iAllBlock[index], iAllIndex[index])) {
                        iNext = getNextTuples(index);
                        assert(iNext);
                        iPItem[index] = iItem[1] + iAllIndex[index][iTreeHeight[index]] * itemSize[1];
                    }
                    else break;
                }
            }
        }
        
        if (return_flag) return_index = index - 1;
        else return_index = pIndex;
        return;
    }
    
    // index based block-nested loop join
    void performBNLJ() {
        ATTR_TYPE attrType[this->nTables];
        ATTR_TYPE pAttrType[this->nTables];
        uint32_t itemSize[this->nTables];
        for (uint32_t j = 0; j < this->nTables; ++j) {
            int32_t tableID = this->iTable[j];
            if (j >= 1) {
                int32_t parentID = this->iParent[j];
                int32_t curAttrID = this->attrID[j];
                int32_t curPAttrID = this->pAttrID[j];
                attrType[j] = this->in_sch[tableID].attrType[curAttrID];
                pAttrType[j] = this->in_sch[parentID].attrType[curPAttrID];
            }
            itemSize[j] = this->in_sch[tableID].item_size;
        }
        
        for (uint32_t j = 1; j < this->nTables; ++j) {
            assert(attrType[j] == CHAR && pAttrType[j] == CHAR ||
                (attrType[j] == INTEGER || attrType[j] == DOUBLE) && (pAttrType[j] == INTEGER || pAttrType[j] == DOUBLE) ||
                (attrType[j] == STRING || attrType[j] == TINYTEXT) && (pAttrType[j] == STRING || pAttrType[j] == TINYTEXT));
        }
        
        iTreeHeight = new int32_t[this->nTables];
        iOutsourcedTreeHeight = new int32_t[this->nTables];
        iAllBlock = new char**[this->nTables];
        iAllBID = new int32_t*[this->nTables];
        iAllIndex = new int32_t*[this->nTables];
        iAllLastEnabled = new bool*[this->nTables];
        iPreLeaf = new char*[this->nTables];
        iPreLeafBID = new int32_t[this->nTables];
        for (int32_t i = 0; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            if (i == 0) {
                iTreeHeight[i] = 0;
                iOutsourcedTreeHeight[i] = 0;
            }
            else {
                int32_t curAttrID = this->attrID[i];
                iTreeHeight[i] = btree[tableID]->getHeight(curAttrID);
                iOutsourcedTreeHeight[i] = btree[tableID]->getOutsourcedHeight(curAttrID);
            }
            iAllBlock[i] = new char*[iTreeHeight[i] + 1];
            for (int32_t j = 0; j <= iTreeHeight[i]; ++j)
                iAllBlock[i][j] = new char[B];
            iAllBID[i] = new int32_t[iTreeHeight[i] + 1];
            iAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
            iAllLastEnabled[i] = new bool[iTreeHeight[i] + 1];
            iPreLeaf[i] = new char[B];
        }
        
        char* iPItem[this->nTables];
        char oBlock[B];
        char oItem[plain_len];
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        iRetrievalNum = 0;
        oResNum = 0;
        
        bool iNext;
        int32_t tableID = this->iTable[0];
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        btree[tableID]->getFirstTuple(iNext, iAllBlock[0][0], iAllBID[0][0], iAllIndex[0][0]);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        
        bool return_flag = false;
        int32_t return_index = -1;
        char* iItem = iAllBlock[0][0] + sizeof(uint32_t);
        iPItem[0] = iItem + iAllIndex[0][0] * itemSize[0];
        join(1, return_flag, return_index, iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        
        for (uint32_t i = 1; i < iItemNum[tableID]; ++i) {
            /***********************/
            start = std::chrono::high_resolution_clock::now();
            /***********************/
            btree[tableID]->getNextTuple(iNext, iAllBlock[0][0], iAllBID[0][0], iAllIndex[0][0]);
            /***********************/
            end = std::chrono::high_resolution_clock::now();
            diff = end - start;
            btree_time += (double)(diff.count());
            /***********************/
            
            return_flag = false;
            return_index = -1;
            iPItem[0] = iItem + iAllIndex[0][0] * itemSize[0];
            join(1, return_flag, return_index, iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        }
        
        uint32_t upperBound = iItemNum[tableID];
        for (uint32_t j = 1; j < this->nTables; ++j) {
            int32_t tableID = this->iTable[j];
            upperBound += 2 * iItemNum[tableID];
        }
        upperBound += oResNum;
        while (iRetrievalNum < upperBound) {
            iNext = getNextTuples(-1);
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
            ++iRetrievalNum;
        }
        while (oItemCnt != 0) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        }
        this->out_tot_block = oBlockID;
        
        for (uint32_t j = 1; j < this->nTables; ++j) {
            int32_t tableID = this->iTable[j];
            int32_t curAttrID = this->attrID[j];
            /***********************/
            start = std::chrono::high_resolution_clock::now();
            /***********************/
            btree[tableID]->resetFlags(curAttrID);
            /***********************/
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
            btree_time += (double)(diff.count());
            /***********************/
        }
        
        for (uint32_t i = 0; i < this->nTables; ++i) {
            for (uint32_t j = 0; j <= iTreeHeight[i]; ++j)
                delete[] iAllBlock[i][j];
            delete[] iAllBlock[i];
        }
        delete[] iAllBlock;
        for (uint32_t i = 0; i < this->nTables; ++i)
            delete[] iAllBID[i];
        delete[] iAllBID;
        for (uint32_t i = 0; i < this->nTables; ++i)
            delete[] iAllIndex[i];
        delete[] iAllIndex;
        for (uint32_t i = 0; i < this->nTables; ++i)
            delete[] iAllLastEnabled[i];
        delete[] iAllLastEnabled;
        for (uint32_t i = 0; i < this->nTables; ++i)
            delete[] iPreLeaf[i];
        delete[] iPreLeaf;
        delete[] iPreLeafBID;
        delete[] iTreeHeight;
        delete[] iOutsourcedTreeHeight;
    }
    
    uint32_t* iItemNum = NULL;
    int32_t* iTreeHeight = NULL;
    int32_t* iOutsourcedTreeHeight = NULL;
    char*** iAllBlock = NULL;
    int32_t** iAllBID = NULL;
    int32_t** iAllIndex = NULL;
    bool** iAllLastEnabled = NULL;
    char** iPreLeaf = NULL;
    int32_t* iPreLeafBID = NULL;
    
    uint32_t iRetrievalNum;
    uint32_t oResNum;
    

    
    /***********************/
    double btree_time;
    /***********************/
};

#endif //__OBLI_MULTIWAY_JOIN_H__
