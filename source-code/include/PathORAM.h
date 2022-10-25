#ifndef __PATHORAM_H__
#define __PATHORAM_H__

#include <unordered_map>
#include <assert.h>
#include "FileSimulator.h"
#include "EncryptionEngine.h"
#include "ORAM.h"
#include "Util.h"

class PathORAM: public ORAM {
public:
    PathORAM(const uint32_t& n, const char* dst, const char* output, const bool& pos_host = true);
    PathORAM(std::unordered_map<uint32_t, std::string>& blocks, const char* dst, const char* output, const bool& pos_host = true);
    PathORAM(const char* input);
    virtual ~PathORAM();
    
    virtual std::string get(const std::string & key);
    virtual void put(const std::string & key, const std::string & value);
    
    virtual std::string get(const int32_t & key);
    virtual std::string getAndRemove(const int32_t & key, const uint32_t & pos);
    virtual void put(const int32_t & key, const std::string & value);
    virtual void put(const int32_t & key, const std::string & value, const uint32_t & pos);
    
    virtual void bulkLoad(std::unordered_map<uint32_t, std::string>& blocks);
    virtual void batchAccess(std::unordered_map<uint32_t, std::string>& blocks);
    
    /***********************/
    virtual uint32_t getDataBlockNum() const;
    virtual size_t getServerSize() const;
    virtual size_t getClientSize() const;
    virtual size_t getCommSize() const;
    virtual void resetCommSize();
    virtual size_t getAccessCount() const;
    
    virtual double getReadTime() const;
    virtual double getWriteTime() const;
    virtual double getEncDecTime() const;
    virtual double getORAMTime() const;
    /***********************/
    
private:
    void access(const char& op, const int32_t& block_id, std::string& data);
    void access(const char& op, const int32_t& block_id, std::string& data, const uint32_t& pos);
    bool check(int x, int y, int l);
    
    void fetchAlongPath(const uint32_t& x, std::string* sbuffer, size_t& length);
    void loadAlongPath(const uint32_t& x, const std::string* sbuffer);
    
    std::unordered_map<uint32_t, std::string> stash;
    std::unordered_map<uint32_t, uint32_t>* pos_map = NULL;
    std::vector< std::pair<uint32_t, std::string> > insert_buffer;
    
    byte* key;
    std::string* sbuffer;
    uint32_t n_data_blocks;
    uint32_t n_blocks;
    uint32_t height;
    ServerConnector* conn;
    encryption_engine engine;
    
    /***********************/
    size_t n_server_blocks, comm_size, max_stash_size, access_count;
    double read_time, write_time, enc_dec_time, oram_time;
    /***********************/
};

#endif //__PATHORAM_H__
