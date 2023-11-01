#ifndef __OBLIVIOUS_HEAP_H__
#define __OBLIVIOUS_HEAP_H__

#include <unordered_map>
#include <assert.h>
#include "FileSimulator.h"
#include "EncryptionEngine.h"
#include "ORAM.h"
#include "Util.h"
#include "Schema.h"

class ObliviousHeap {
public:
    ObliviousHeap(ServerConnector* _conn, const uint32_t _n_block, const Schema* _sort_sch, const CMP* _sort_cmp);
    ~ObliviousHeap();
    
    /***********************/
    std::string findMin();
    std::pair<int32_t, int64_t> insertItem(int32_t item_id, std::string content);
    std::string deleteItem(int32_t pos, int64_t ts);
    std::string extractMin();
    /***********************/
    
    /***********************/
    uint32_t getDataItemNum() const;
    size_t getServerSize() const;
    size_t getClientSize() const;
    size_t getCommSize() const;
    void resetCommSize();
    /***********************/
    
private:
    /***********************/
    std::string genDummyItem();
    int32_t P(int32_t pos, int32_t level);
    void readBucket(const int32_t bucket_pos, std::vector<std::string> & sbuffer);
    void writeBucket(const int32_t bucket_pos, std::vector<std::string> & sbuffer);
    /***********************/
    
    /***********************/
    int itemCompare(const void* item1, const void* item2);
    std::string onlyFindMin();
    std::string onlyDeleteItem(int32_t pos, int64_t timestamp);
    std::pair<int32_t, int64_t> onlyInsertItem(int32_t item_id, std::string content);
    void insertDummyItem();
    void evictAfterInsert();
    void evictAndUpdateMin(int32_t pos);
    std::string readAndRemove(int32_t pos, int64_t timestamp);
    void updateMin(int32_t pos);
    /***********************/
    
    std::unordered_map<uint32_t, std::string> stash;
    std::vector<std::string> sbuffer;
    
    bool type_hiding_security;
    int64_t timestamp;
    
    byte* key;
    uint32_t n_data_items;
    uint32_t n_items;
    uint32_t enc_item_size;
    uint32_t height;
    
    ServerConnector* conn;
    encryption_engine engine;
    
    const uint32_t n_input_blocks;
    const Schema* sort_sch = NULL;
    const CMP* sort_cmp = NULL;
    
    /***********************/
    size_t n_server_items, comm_size, max_stash_size, access_count;
    /***********************/
};

#endif //__OBLIVIOUS_HEAP_H__
