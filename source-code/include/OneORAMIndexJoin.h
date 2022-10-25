#ifndef __ONE_ORAM_INDEX_JOIN_H__
#define __ONE_ORAM_INDEX_JOIN_H__

#include "PrivIndexJoin.h"
#include "OneORAMBTree.h"

using namespace std;

template<class T>
class OneORAMIndexJoin : public PrivIndexJoin {
public:
    // load
    OneORAMIndexJoin(const uint32_t mode, const std::string scale, const uint32_t outsourced_height, const std::string& dst, const std::string& prefix) : PrivIndexJoin (mode, dst, prefix) {
        /***********************/
        btree_time = 0.0;
        /***********************/
        
        uint32_t block_id = 0;
        std::unordered_map <uint32_t, std::string> blocks;
        btree = new OneORAMBTree <T>* [this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            generateInputSchema(mode, i, &(this->in_sch[i]));
            
            std::string cur_prefix = prefix + "_input_" + std::to_string(i);
            btree[i] = new OneORAMBTree <T> (mode, scale, i, &(this->in_sch[i]), outsourced_height, blocks, block_id, cur_prefix);
        }
        initBufferPara(getDataBlockNum());
        
        // build the ORAM storage
        std::string in_dst = dst + "_input";
        std::string in_prefix = prefix + "_input";
        one_oram = new T(blocks, in_dst.c_str(), (in_prefix + "_oram.txt").c_str());
        resetCommSize();
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            btree[i]->initORAM(one_oram);
        }
    }
    
    // run
    OneORAMIndexJoin(const uint32_t mode, const std::string& prefix) : PrivIndexJoin (mode, prefix) {
        /***********************/
        btree_time = 0.0;
        /***********************/
        
        btree = new OneORAMBTree <T>* [this->nInputs];
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            generateInputSchema(mode, i, &(this->in_sch[i]));
            
            std::string cur_prefix = prefix + "_input_" + std::to_string(i);
            btree[i] = new OneORAMBTree <T> (&(this->in_sch[i]), cur_prefix);
        }
        initBufferPara(getDataBlockNum());
        
        // build the ORAM storage
        printf("Start Loading Mongo....\n");
        std::string in_prefix = prefix + "_input";
        one_oram = new T((in_prefix + "_oram.txt").c_str());
        resetCommSize();
        for (uint32_t i = 0; i < this->nInputs; ++i) {
            btree[i]->initORAM(one_oram);
        }
        printf("Finish Loading Mongo...\n");
    }
    
    ~OneORAMIndexJoin() {
        if (one_oram != NULL) {
            delete one_oram;
            one_oram = NULL;
        }
        if (btree != NULL) {
            for (uint32_t i = 0; i < this->nInputs; ++i)
                delete btree[i];
            delete[] btree;
            btree = NULL;
        }
    }
    
    // oblivious sort-merge equi-join
    void ObliSMEJ(const uint32_t n_tables, const uint32_t* table_id, const int32_t* parent_id, const int32_t* attr_id, const int32_t* parent_attr_id, const uint32_t n_proj_cols, const uint32_t proj_id[][MAX_COLS]) {
        assert(n_tables == 2);
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
        performSMEJ();
        
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
    
    // oblivious sort-merge band join
    void ObliSMBJ(const uint32_t n_tables, const uint32_t* table_id, const int32_t* parent_id, const int32_t* attr_id, const int32_t* parent_attr_id, const double* band_range, const uint32_t n_proj_cols, const uint32_t proj_id[][MAX_COLS]) {
        assert(n_tables == 2);
        this->nTables = n_tables;
        this->iTable = new uint32_t[this->nTables];
        this->attrID = new int32_t[this->nTables];
        this->iParent = new int32_t[this->nTables];
        this->pAttrID = new int32_t[this->nTables];
        this->bRange = new double[2];
        memcpy(this->iTable, table_id, this->nTables * sizeof(uint32_t));
        memcpy(this->attrID, attr_id, this->nTables * sizeof(int32_t));
        memcpy(this->iParent, parent_id, this->nTables * sizeof(int32_t));
        memcpy(this->pAttrID, parent_attr_id, this->nTables * sizeof(int32_t));
        memcpy(this->bRange, band_range, 2 * sizeof(double));
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
        for (uint32_t i = 0; i < 2; ++i)
            out_dst += "_" + std::to_string(this->bRange[i]);
        this->out_conn = new FileSimulator(server_host, out_dst.c_str(), true);
        this->out_sch = new Schema;
        if (this->nProjCols == 0) generateOutputSchema(this->in_sch, this->nTables, this->iTable, this->attrID, this->out_sch, false);
        else generateOutputSchema(this->in_sch, this->nTables, this->iTable, this->nProjCols, this->projID, this->out_sch);
        // TODO: support ORDER-BY operator
        if (this->nProjCols == 0) this->out_cmp = new CMP({1, {this->pAttrID[1]}, {0}});
        else this->out_cmp = new CMP({1, {1}, {0}});

        // perform oblivious join
        performSMBJ();

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
        if (this->bRange != NULL) {
            delete[] this->bRange;
            this->bRange = NULL;
        }
        if (this->nProjCols > 0) {
            for (uint32_t i = 0; i < 2; ++i) {
                if (this->projID[i] != NULL) {
                    delete[](this->projID[i]);
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
    
    // oblivious index-nested loop equi-join
    void ObliINLEJ(const uint32_t n_tables, const uint32_t* table_id, const int32_t* parent_id, const int32_t* attr_id, const int32_t* parent_attr_id, const uint32_t n_proj_cols, const uint32_t proj_id[][MAX_COLS]) {
        assert(n_tables == 2);
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
        performINLEJ();
        
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
    
    // oblivious index-nested loop band join
    void ObliINLBJ(const uint32_t n_tables, const uint32_t* table_id, const int32_t* parent_id, const int32_t* attr_id, const int32_t* parent_attr_id, const double* band_range, const uint32_t n_proj_cols, const uint32_t proj_id[][MAX_COLS]) {
        assert(n_tables == 2);
        this->nTables = n_tables;
        this->iTable = new uint32_t[this->nTables];
        this->attrID = new int32_t[this->nTables];
        this->iParent = new int32_t[this->nTables];
        this->pAttrID = new int32_t[this->nTables];
        this->bRange = new double[2];
        memcpy(this->iTable, table_id, this->nTables * sizeof(uint32_t));
        memcpy(this->attrID, attr_id, this->nTables * sizeof(int32_t));
        memcpy(this->iParent, parent_id, this->nTables * sizeof(int32_t));
        memcpy(this->pAttrID, parent_attr_id, this->nTables * sizeof(int32_t));
        memcpy(this->bRange, band_range, 2 * sizeof(double));
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
        for (uint32_t i = 0; i < 2; ++i)
            out_dst += "_" + std::to_string(this->bRange[i]);
        this->out_conn = new FileSimulator(server_host, out_dst.c_str(), true);
        this->out_sch = new Schema;
        if (this->nProjCols == 0) generateOutputSchema(this->in_sch, this->nTables, this->iTable, this->attrID, this->out_sch, false);
        else generateOutputSchema(this->in_sch, this->nTables, this->iTable, this->nProjCols, this->projID, this->out_sch);
        // TODO: support ORDER-BY operator
        if (this->nProjCols == 0) this->out_cmp = new CMP({1, {this->pAttrID[1]}, {0}});
        else this->out_cmp = new CMP({1, {1}, {0}});
        
        // perform oblivious join
        performINLBJ();
        
        // perform oblivious filter
        uint32_t access_num = 0;
        filterOutput(this->out_conn, this->out_tot_block, this->out_sch, this->out_cmp, access_num);
        this->comm_size += access_num;
        
        /***********************/
        // print join result
        if (this->nProjCols == 0) {
            int32_t printAttrID[this->nTables];
            printAttrID[0] = this->pAttrID[1];
            printAttrID[1] = this->in_sch->nAttrs + this->attrID[1] - 1;
            printResult(this->out_conn, this->out_tot_block, this->out_sch, 2, printAttrID);
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
        if (this->bRange != NULL) {
            delete[] this->bRange;
            this->bRange = NULL;
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
        size_t server_size = one_oram->getServerSize();
        return server_size;
    }
    
    size_t getClientSize() const {
        size_t client_size = oblisort_mem;
        printf("Additional trusted memory size is %lf MB\n", (double)oblisort_mem / (1 << 20));
        client_size += one_oram->getClientSize();
        for (uint32_t i = 0; i < nInputs; ++i)
            client_size += btree[i]->getClientSize();
        return client_size;
    }
    
    size_t getCommSize() const {
        size_t total_comm_size = one_oram->getCommSize();
        printf("One ORAM comm size is %lf MB\n", (double)total_comm_size / (1 << 20));
        total_comm_size += this->comm_size * (size_t)B;
        return total_comm_size;
    }
    
    void resetCommSize() {
        this->comm_size = 0;
        one_oram->resetCommSize();
    }
    
    size_t getAccessCount() const {
        size_t access_count = one_oram->getAccessCount();
        return access_count;
    }
    
    double getReadTime() const {
        double read_time = one_oram->getReadTime();
        return read_time;
    }
    
    double getWriteTime() const {
        double write_time = one_oram->getWriteTime();
        return write_time;
    }
    
    double getEncDecTime() const {
        double enc_dec_time = one_oram->getEncDecTime();
        return enc_dec_time;
    }
    
    double getORAMTime() const {
        double oram_time = one_oram->getORAMTime();
        return oram_time;
    }
    
    double getBTreeTime() const {
        return btree_time;
    }
    /***********************/
    
private:
    void getDummyBlock() {
        int32_t dummyID = -1;
        one_oram->get(dummyID);
    }
    
    // for index based sort-merge equi-join
    void getNextTuple(int32_t iIndex, bool* iNext, char* iLeaf, int32_t* iLeafIndex, char* iData, int32_t* iDataBID, int32_t* iDataIndex, int32_t* iOutsourcedTreeHeight, int32_t iMaxOutsourcedTreeHeight) {
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        assert(iIndex >= -1);
        if (iIndex == -1) {
            if (iMaxOutsourcedTreeHeight > 0)
                getDummyBlock();
            getDummyBlock();
        }
        else {
            int32_t tableID = this->iTable[iIndex];
            int32_t curAttrID;
            if (iIndex == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[iIndex];
            btree[tableID]->getNextTuple(iNext[iIndex], curAttrID, iLeaf + iIndex * B, iLeafIndex[iIndex], iData + iIndex * B, iDataBID[iIndex], iDataIndex[iIndex]);
            if (iOutsourcedTreeHeight[iIndex] <= 0 && iMaxOutsourcedTreeHeight > 0)
                getDummyBlock();
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
    }
    
    // index based sort-merge equi-join
    void performSMEJ() {
        ATTR_TYPE attrType[this->nTables];
        uint32_t attrSize[this->nTables];
        uint32_t offset[this->nTables];
        uint32_t itemSize[this->nTables];
        for (uint32_t j = 0; j < this->nTables; ++j) {
            int32_t tableID = this->iTable[j];
            int32_t curAttrID;
            if (j == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[j];
            attrType[j] = this->in_sch[tableID].attrType[curAttrID];
            attrSize[j] = this->in_sch[tableID].attrSize[curAttrID];
            offset[j] = this->in_sch[tableID].attrOffset[curAttrID];
            itemSize[j] = this->in_sch[tableID].item_size;
        }
        
        if (attrType[0] == CHAR) {
            for (uint32_t j = 1; j < this->nTables; ++j)
                assert(attrType[j] == CHAR);
        }
        else if (attrType[0] == INTEGER || attrType[0] == DOUBLE) {
            for (uint32_t j = 1; j < this->nTables; ++j)
                assert(attrType[j] == INTEGER || attrType[j] == DOUBLE);
        }
        else if (attrType[0] == STRING || attrType[0] == TINYTEXT) {
            for (uint32_t j = 1; j < this->nTables; ++j)
                assert(attrType[j] == STRING || attrType[j] == TINYTEXT);
        }
        
        char iLeaf[this->nTables][B];
        char iData[this->nTables][B];
        int32_t iLeafIndex[this->nTables];
        int32_t iDataBID[this->nTables];
        int32_t iDataIndex[this->nTables];
        
        char iBackLeaf[this->nTables][B];
        char iBackData[this->nTables][B];
        int32_t iBackLeafIndex[this->nTables];
        int32_t iBackDataBID[this->nTables];
        int32_t iBackDataIndex[this->nTables];
        
        int32_t iOutsourcedTreeHeight[this->nTables];
        int32_t iMaxOutsourcedTreeHeight = 0;
        for (uint32_t i = 0; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            int32_t curAttrID;
            if (i == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[i];
            iOutsourcedTreeHeight[i] = btree[tableID]->getOutsourcedHeight(curAttrID);
            if (iOutsourcedTreeHeight[i] > iMaxOutsourcedTreeHeight)
                iMaxOutsourcedTreeHeight = iOutsourcedTreeHeight[i];
        }
        
        bool iNext[this->nTables];
        memset(iNext, true, this->nTables * sizeof(bool));
        
        char* iItem[this->nTables];
        char* iPItem[this->nTables];
        char* iPAttr[this->nTables];
        for (uint32_t i = 0; i < this->nTables; ++i)
            iItem[i] = iData[i] + sizeof(uint32_t);
        
        char oBlock[B];
        char oItem[plain_len];
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        for (uint32_t i = 0; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            int32_t curAttrID;
            if (i == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[i];
            btree[tableID]->getFirstTuple(iNext[i], curAttrID, iLeaf[i], iLeafIndex[i], iData[i], iDataBID[i], iDataIndex[i]);
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        
        while (iNext[0] || iNext[1]) {
            int cmpres;
            if (!iNext[0]) cmpres = 1;
            else if (!iNext[1]) cmpres = -1;
            else {
                for (uint32_t i = 0; i < this->nTables; ++i) {
                    iPItem[i] = iItem[i] + iDataIndex[i] * itemSize[i];
                    iPAttr[i] = iPItem[i] + offset[i];
                }
                cmpres = this->joinCompare(iPAttr, attrType, attrSize);
            }
            
            if (cmpres == 0) {
                memcpy(iBackLeaf[1], iLeaf[1], B);
                memcpy(iBackData[1], iData[1], B);
                iBackLeafIndex[1] = iLeafIndex[1];
                iBackDataBID[1] = iDataBID[1];
                iBackDataIndex[1] = iDataIndex[1];

                while (cmpres == 0) {
                    /*******************
                    printf("cmpres == 0\n");
                    *******************/
                    this->writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                    getNextTuple(1, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
                    if (!iNext[0]) cmpres = 1;
                    else if (!iNext[1]) cmpres = -1;
                    else {
                        for (uint32_t i = 0; i < this->nTables; ++i) {
                            iPItem[i] = iItem[i] + iDataIndex[i] * itemSize[i];
                            iPAttr[i] = iPItem[i] + offset[i];
                        }
                        cmpres = this->joinCompare(iPAttr, attrType, attrSize);
                    }
                }
                /*******************
                printf("cmpres != 0\n");
                *******************/
                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                
                memcpy(iLeaf[1], iBackLeaf[1], B);
                memcpy(iData[1], iBackData[1], B);
                iLeafIndex[1] = iBackLeafIndex[1];
                iDataBID[1] = iBackDataBID[1];
                iDataIndex[1] = iBackDataIndex[1];
                iNext[1] = true;
                getNextTuple(0, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
            }
            else {
                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                if (cmpres < 0) {
                    /*******************
                    printf("cmpres < 0\n");
                    *******************/
                    getNextTuple(0, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
                }
                else if (cmpres > 0) {
                    /*******************
                    printf("cmpres > 0\n");
                    *******************/
                    getNextTuple(1, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
                }
            }
        }
        while (oItemCnt != 0) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        }
        this->out_tot_block = oBlockID;
    }
    
    // index based sort-merge band join
    void performSMBJ() {
        ATTR_TYPE attrType[this->nTables];
        uint32_t attrSize[this->nTables];
        uint32_t offset[this->nTables];
        uint32_t itemSize[this->nTables];
        for (uint32_t j = 0; j < this->nTables; ++j) {
            int32_t tableID = this->iTable[j];
            int32_t curAttrID;
            if (j == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[j];
            attrType[j] = this->in_sch[tableID].attrType[curAttrID];
            attrSize[j] = this->in_sch[tableID].attrSize[curAttrID];
            offset[j] = this->in_sch[tableID].attrOffset[curAttrID];
            itemSize[j] = this->in_sch[tableID].item_size;
        }
        
        if (attrType[0] == CHAR) {
            for (uint32_t j = 1; j < this->nTables; ++j)
                assert(attrType[j] == CHAR);
        }
        else if (attrType[0] == INTEGER || attrType[0] == DOUBLE) {
            for (uint32_t j = 1; j < this->nTables; ++j)
                assert(attrType[j] == INTEGER || attrType[j] == DOUBLE);
        }
        else if (attrType[0] == STRING || attrType[0] == TINYTEXT) {
            for (uint32_t j = 1; j < this->nTables; ++j)
                assert(attrType[j] == STRING || attrType[j] == TINYTEXT);
        }
        
        char iLeaf[this->nTables][B];
        char iData[this->nTables][B];
        int32_t iLeafIndex[this->nTables];
        int32_t iDataBID[this->nTables];
        int32_t iDataIndex[this->nTables];

        char iBackLeaf[this->nTables][B];
        char iBackData[this->nTables][B];
        int32_t iBackLeafIndex[this->nTables];
        int32_t iBackDataBID[this->nTables];
        int32_t iBackDataIndex[this->nTables];

        int32_t iOutsourcedTreeHeight[this->nTables];
        int32_t iMaxOutsourcedTreeHeight = 0;
        for (uint32_t i = 0; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            int32_t curAttrID;
            if (i == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[i];
            iOutsourcedTreeHeight[i] = btree[tableID]->getOutsourcedHeight(curAttrID);
            if (iOutsourcedTreeHeight[i] > iMaxOutsourcedTreeHeight)
                iMaxOutsourcedTreeHeight = iOutsourcedTreeHeight[i];
        }

        bool iNext[this->nTables];
        memset(iNext, true, this->nTables * sizeof(bool));

        char* iItem[this->nTables];
        char* iPItem[this->nTables];
        char* iPAttr[this->nTables];
        for (uint32_t i = 0; i < this->nTables; ++i)
            iItem[i] = iData[i] + sizeof(uint32_t);

        char oBlock[B];
        char oItem[plain_len];
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;

        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        for (uint32_t i = 0; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            int32_t curAttrID;
            if (i == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[i];
            btree[tableID]->getFirstTuple(iNext[i], curAttrID, iLeaf[i], iLeafIndex[i], iData[i], iDataBID[i], iDataIndex[i]);
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/

        while (iNext[0] || iNext[1]) {
            int cmpres;
            if (!iNext[0]) cmpres = 1;
            else if (!iNext[1]) cmpres = -1;
            else {
                for (uint32_t i = 0; i < this->nTables; ++i) {
                    iPItem[i] = iItem[i] + iDataIndex[i] * itemSize[i];
                    iPAttr[i] = iPItem[i] + offset[i];
                }
                cmpres = this->joinCompare(iPAttr, attrType, attrSize, this->bRange);
            }

            if (cmpres == 0) {
                memcpy(iBackLeaf[1], iLeaf[1], B);
                memcpy(iBackData[1], iData[1], B);
                iBackLeafIndex[1] = iLeafIndex[1];
                iBackDataBID[1] = iDataBID[1];
                iBackDataIndex[1] = iDataIndex[1];

                while (cmpres == 0) {
                    /*******************
                    printf("cmpres == 0\n");
                    *******************/
                    this->writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
                    getNextTuple(1, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
                    if (!iNext[0]) cmpres = 1;
                    else if (!iNext[1]) cmpres = -1;
                    else {
                        for (uint32_t i = 0; i < this->nTables; ++i) {
                            iPItem[i] = iItem[i] + iDataIndex[i] * itemSize[i];
                            iPAttr[i] = iPItem[i] + offset[i];
                        }
                        cmpres = this->joinCompare(iPAttr, attrType, attrSize, this->bRange);
                    }
                }
                /*******************
                printf("cmpres != 0\n");
                *******************/
                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);

                memcpy(iLeaf[1], iBackLeaf[1], B);
                memcpy(iData[1], iBackData[1], B);
                iLeafIndex[1] = iBackLeafIndex[1];
                iDataBID[1] = iBackDataBID[1];
                iDataIndex[1] = iBackDataIndex[1];
                iNext[1] = true;
                getNextTuple(0, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
            }
            else {
                this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
                if (cmpres < 0) {
                    /*******************
                    printf("cmpres < 0\n");
                    *******************/
                    getNextTuple(0, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
                }
                else if (cmpres > 0) {
                    /*******************
                    printf("cmpres > 0\n");
                    *******************/
                    getNextTuple(1, iNext, *iLeaf, iLeafIndex, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
                }
            }
        }
        while (oItemCnt != 0) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
        }
        this->out_tot_block = oBlockID;
    }
    
    // for the first table in oblivious index nested-loop equi-join
    void getNextTuple(int32_t iIndex, bool* iNext, char* iData, int32_t* iDataBID, int32_t& iDataIndex, int32_t* iOutsourcedTreeHeight, int32_t iMaxOutsourcedTreeHeight) {
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        assert(iIndex == 0);
        int32_t tableID = this->iTable[0];
        btree[tableID]->getNextTuple(iNext[0], iData, iDataBID[0], iDataIndex);
        for (uint32_t j = 1; j <= iMaxOutsourcedTreeHeight; ++j)
            getDummyBlock();
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
    }
    
    // for the other table in oblivious index nested-loop band join
    void getNextTuple(int32_t iIndex, bool* iNext, char* iLeaf, char* iData, int32_t* iDataBID, int32_t** iAllIndex, int32_t* iEndLevel, int32_t* iOutsourcedTreeHeight, int32_t iMaxOutsourcedTreeHeight) {
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        assert(iIndex == -1 || iIndex > 0);
        if (iIndex == -1) {
            for (uint32_t j = 0; j <= iMaxOutsourcedTreeHeight; ++j)
                getDummyBlock();
        }
        else {
            int32_t tableID = this->iTable[iIndex];
            btree[tableID]->getNextTuple(iNext[iIndex], this->attrID[iIndex], iLeaf + iIndex * B, iData + iIndex * B, iDataBID[iIndex], iAllIndex[iIndex], iEndLevel[iIndex]);
            for (uint32_t j = iOutsourcedTreeHeight[iIndex] + 1; j <= iMaxOutsourcedTreeHeight; ++j)
                getDummyBlock();
        }
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
    }
    
    // oblivious index nested-loop equi-join
    void performINLEJ() {
        ATTR_TYPE attrType[this->nTables];
        uint32_t attrSize[this->nTables];
        uint32_t offset[this->nTables];
        uint32_t itemSize[this->nTables];
        for (uint32_t j = 0; j < this->nTables; ++j) {
            int32_t tableID = this->iTable[j];
            int32_t curAttrID;
            if (j == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[j];
            attrType[j] = this->in_sch[tableID].attrType[curAttrID];
            attrSize[j] = this->in_sch[tableID].attrSize[curAttrID];
            offset[j] = this->in_sch[tableID].attrOffset[curAttrID];
            itemSize[j] = this->in_sch[tableID].item_size;
        }

        if (attrType[0] == CHAR)
            assert(attrType[1] == CHAR);
        else if (attrType[0] == INTEGER || attrType[0] == DOUBLE)
            assert(attrType[1] == INTEGER || attrType[1] == DOUBLE);
        else if (attrType[0] == STRING || attrType[0] == TINYTEXT)
            assert(attrType[1] == STRING || attrType[1] == TINYTEXT);

        char iLeaf[this->nTables][B];
        char iData[this->nTables][B];
        int32_t iDataBID[this->nTables];

        int32_t iTreeHeight[this->nTables];
        int32_t iOutsourcedTreeHeight[this->nTables];
        int32_t iLeafIndex, iDataIndex;
        int32_t* iAllIndex[this->nTables];
        int32_t iEndLevel[this->nTables];
        for (uint32_t i = 0; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            int32_t curAttrID;
            if (i == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[i];
            iTreeHeight[i] = btree[tableID]->getHeight(curAttrID);
            iOutsourcedTreeHeight[i] = btree[tableID]->getOutsourcedHeight(curAttrID);
            iAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
        }
        int32_t iMaxOutsourcedTreeHeight = 0;
        for (uint32_t i = 1; i < this->nTables; ++i) {
            if (iOutsourcedTreeHeight[i] > iMaxOutsourcedTreeHeight)
                iMaxOutsourcedTreeHeight = iOutsourcedTreeHeight[i];
        }
        if (iMaxOutsourcedTreeHeight <= 0 && iOutsourcedTreeHeight[0] > 0)
            iMaxOutsourcedTreeHeight = 1;

        bool iNext[this->nTables];
        char* iItem[this->nTables];
        char* iPItem[this->nTables];
        char* iPAttr[this->nTables];
        for (uint32_t i = 0; i < this->nTables; ++i)
            iItem[i] = iData[i] + sizeof(uint32_t);

        char oBlock[B];
        char oItem[plain_len];
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;

        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        int32_t tableID = this->iTable[0];
        btree[tableID]->getFirstTuple(iNext[0], iData[0], iDataBID[0], iDataIndex);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        while (iNext[0]) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);

            iPItem[0] = iItem[0] + iDataIndex * itemSize[0];
            iPAttr[0] = iPItem[0] + offset[0];
            /***********************/
            start = std::chrono::high_resolution_clock::now();
            /***********************/
            tableID = this->iTable[1];
            btree[tableID]->getFirstTuple(iNext[1], iPAttr[0], attrType[0], attrSize[0], '=', this->attrID[1], iLeaf[1], iData[1], iDataBID[1], iAllIndex[1], iEndLevel[1]);
            /***********************/
            end = std::chrono::high_resolution_clock::now();
            diff = end - start;
            btree_time += (double)(diff.count());
            /***********************/

            while (iNext[1]) {
                iPItem[1] = iItem[1] + iAllIndex[1][iTreeHeight[1]] * itemSize[1];
                iPAttr[1] = iPItem[1] + offset[1];
                int cmpres = this->joinCompare(iPAttr, attrType, attrSize);
                if (cmpres == 0) this->writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
                else break;
                getNextTuple(1, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
            }
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);

            getNextTuple(0, iNext, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
        }
        while (oItemCnt != 0) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID);
        }
        this->out_tot_block = oBlockID;

        for (uint32_t i = 0; i < this->nTables; ++i)
            delete[] iAllIndex[i];
    }
    
    // oblivious index nested-loop band join
    void performINLBJ() {
        ATTR_TYPE attrType[this->nTables];
        uint32_t attrSize[this->nTables];
        uint32_t offset[this->nTables];
        uint32_t itemSize[this->nTables];
        for (uint32_t j = 0; j < this->nTables; ++j) {
            int32_t tableID = this->iTable[j];
            int32_t curAttrID;
            if (j == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[j];
            attrType[j] = this->in_sch[tableID].attrType[curAttrID];
            attrSize[j] = this->in_sch[tableID].attrSize[curAttrID];
            offset[j] = this->in_sch[tableID].attrOffset[curAttrID];
            itemSize[j] = this->in_sch[tableID].item_size;
        }
        
        if (attrType[0] == CHAR)
            assert(attrType[1] == CHAR);
        else if (attrType[0] == INTEGER || attrType[0] == DOUBLE)
            assert(attrType[1] == INTEGER || attrType[1] == DOUBLE);
        else if (attrType[0] == STRING || attrType[0] == TINYTEXT)
            assert(attrType[1] == STRING || attrType[1] == TINYTEXT);
        
        char iLeaf[this->nTables][B];
        char iData[this->nTables][B];
        int32_t iDataBID[this->nTables];
        
        int32_t iTreeHeight[this->nTables];
        int32_t iOutsourcedTreeHeight[this->nTables];
        int32_t iLeafIndex, iDataIndex;
        int32_t* iAllIndex[this->nTables];
        int32_t iEndLevel[this->nTables];
        for (uint32_t i = 0; i < this->nTables; ++i) {
            int32_t tableID = this->iTable[i];
            int32_t curAttrID;
            if (i == 0) curAttrID = this->pAttrID[1];
            else curAttrID = this->attrID[i];
            iTreeHeight[i] = btree[tableID]->getHeight(curAttrID);
            iOutsourcedTreeHeight[i] = btree[tableID]->getOutsourcedHeight(curAttrID);
            iAllIndex[i] = new int32_t[iTreeHeight[i] + 1];
        }
        int32_t iMaxOutsourcedTreeHeight = 0;
        for (uint32_t i = 1; i < this->nTables; ++i) {
            if (iOutsourcedTreeHeight[i] > iMaxOutsourcedTreeHeight)
                iMaxOutsourcedTreeHeight = iOutsourcedTreeHeight[i];
        }
        if (iMaxOutsourcedTreeHeight <= 0 && iOutsourcedTreeHeight[0] > 0)
            iMaxOutsourcedTreeHeight = 1;
        
        bool iNext[this->nTables];
        char* iItem[this->nTables];
        char* iPItem[this->nTables];
        char* iPAttr[this->nTables];
        for (uint32_t i = 0; i < this->nTables; ++i)
            iItem[i] = iData[i] + sizeof(uint32_t);
        
        char oBlock[B];
        char oItem[plain_len];
        uint32_t oBlockID = 0;
        uint32_t oItemCnt = 0;
        uint32_t oRItemCnt = 0;
        
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        int32_t tableID = this->iTable[0];
        btree[tableID]->getFirstTuple(iNext[0], iData[0], iDataBID[0], iDataIndex);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        btree_time += (double)(diff.count());
        /***********************/
        while (iNext[0]) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
            
            iPItem[0] = iItem[0] + iDataIndex * itemSize[0];
            iPAttr[0] = iPItem[0] + offset[0];
            /***********************/
            start = std::chrono::high_resolution_clock::now();
            /***********************/
            tableID = this->iTable[1];
            btree[tableID]->getFirstTuple(iNext[1], iPAttr[0], attrType[0], attrSize[0], 'b', this->attrID[1], iLeaf[1], iData[1], iDataBID[1], iAllIndex[1], iEndLevel[1], this->bRange);
            /***********************/
            end = std::chrono::high_resolution_clock::now();
            diff = end - start;
            btree_time += (double)(diff.count());
            /***********************/
            
            while (iNext[1]) {
                iPItem[1] = iItem[1] + iAllIndex[1][iTreeHeight[1]] * itemSize[1];
                iPAttr[1] = iPItem[1] + offset[1];
                int cmpres = this->joinCompare(iPAttr, attrType, attrSize, this->bRange);
                if (cmpres == 0) this->writeItem('r', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
                else break;
                getNextTuple(1, iNext, *iLeaf, *iData, iDataBID, iAllIndex, iEndLevel, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
            }
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
            
            getNextTuple(0, iNext, *iData, iDataBID, iDataIndex, iOutsourcedTreeHeight, iMaxOutsourcedTreeHeight);
        }
        while (oItemCnt != 0) {
            this->writeItem('d', iPItem, oBlock, oItem, oItemCnt, oRItemCnt, oBlockID, false);
        }
        this->out_tot_block = oBlockID;
        
        for (uint32_t i = 0; i < this->nTables; ++i)
            delete[] iAllIndex[i];
    }
    
    // private ORAM storage
    ORAM* one_oram = NULL;
    // private B-tree
    OneORAMBTree <T>** btree = NULL;
    
    /***********************/
    double btree_time;
    /***********************/
};

#endif //__ONE_ORAM_INDEX_JOIN_H__
