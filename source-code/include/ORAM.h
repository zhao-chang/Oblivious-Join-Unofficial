#ifndef __ORAM_H__
#define __ORAM_H__

#include <string>

class ORAM {
public:
    ORAM() {}
    virtual ~ORAM() {}
    
    virtual std::string get(const int32_t & key) = 0;
    virtual std::string getAndRemove(const int32_t & key, const uint32_t& pos) = 0;
    virtual void put(const int32_t& key, const std::string & value) = 0;
    virtual void put(const int32_t& key, const std::string & value, const uint32_t& pos) = 0;
    
    virtual void bulkLoad(std::unordered_map<uint32_t, std::string>& blocks) = 0;
    virtual void batchAccess(std::unordered_map<uint32_t, std::string>& blocks) = 0;
    
    /***********************/
    virtual uint32_t getDataBlockNum() const = 0;
    virtual size_t getServerSize() const = 0;
    virtual size_t getClientSize() const = 0;
    virtual size_t getCommSize() const = 0;
    virtual void resetCommSize() = 0;
    virtual size_t getAccessCount() const = 0;
    
    virtual double getReadTime() const = 0;
    virtual double getWriteTime() const = 0;
    virtual double getEncDecTime() const = 0;
    virtual double getORAMTime() const = 0;
    /***********************/
};

#endif //__ORAM_H__
