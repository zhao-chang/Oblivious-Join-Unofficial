#ifndef __MONGOCONNECTOR_H__
#define __MONGOCONNECTOR_H__

#include <mongo/client/dbclient.h>
#include "ServerConnector.h"

using namespace mongo;

class MongoConnector: public ServerConnector {
public:
    MongoConnector(const std::string& url, const std::string& collection_name, const bool create = false);
    MongoConnector(const std::string& url);
    virtual ~MongoConnector();
    
    struct iterator: public ServerConnector::iterator {
        iterator(std::unique_ptr<DBClientCursor> c): cursor(std::move(c)) {}

        virtual ~iterator() {}

        virtual bool hasNext() {
            return cursor->more();
        }

        virtual std::string next() {
            BSONObj p = cursor->next();
            int len;
            const char* raw_data = p.getField("data").binData(len);
            return std::string(raw_data, (size_t)len);
        }

        std::unique_ptr<DBClientCursor> cursor;
    };
    
    virtual void clear(const std::string& ns = "");
    virtual void resize(const uint32_t& len, const std::string& ns = "");
    virtual void insert(const uint32_t& id, const std::string& encrypted_block, const std::string& ns = "");
    virtual void insert(const std::vector< std::pair<uint32_t, std::string> >& blocks, const std::string& ns = "");
    virtual void insert(const std::string* sbuffer, const uint32_t& low, const size_t& len, const std::string& ns = "");
    virtual void insert(const char* blocks, const uint32_t& low, const uint32_t& high, const std::string& ns = "");
    virtual void insert(const std::vector< std::pair<std::string, std::string> >& blocks, const std::string& ns = "");
    virtual void insertWithTag(const std::vector< std::pair<std::string, std::string> >& blocks, const std::string& ns = "");
    virtual iterator* scan();
    virtual std::string find(const uint32_t& id, const std::string& ns = "");
    virtual void find(const std::vector<uint32_t>& ids, std::string* sbuffer, size_t& length, const std::string& ns = "");
    virtual std::string fetch(const std::string& id, const std::string& ns = "");
    virtual std::string fetch(const uint32_t& id, const std::string& ns = "");
    virtual void fetch(const uint32_t& low, const uint32_t& high, char* blocks, const std::string& ns = "");
    virtual void find(const uint32_t& low, const uint32_t& high, std::vector<std::string>& blocks, const std::string& ns = "");
    virtual void find(const uint32_t& low, const uint32_t& high, char* blocks, const std::string& ns = "");
    virtual void findByTag(const uint32_t& tag, std::string* sbuffer, size_t& length, const std::string& ns = "");
    virtual void update(const uint32_t& id, const std::string& data, const std::string& ns = "");
    virtual void update(const std::string* sbuffer, const uint32_t& low, const size_t& len, const std::string& ns = "");
    virtual void update(const char* blocks, const uint32_t& low, const uint32_t& high, const std::string& ns = "");
    virtual void update(const std::vector< std::pair<uint32_t, std::string> > blocks, const std::string& ns = "");
    virtual void remove(const uint32_t& id, const std::string& ns = "");
    virtual void remove(const std::string& id, const std::string& ns = "");
    virtual void remove(const uint32_t& low, const uint32_t& high, const std::string& ns = "");
    
    virtual void initialize(const std::string& ns = "");
    virtual void finalize(const std::string& ns = "");
    
    virtual std::string getCollectionName();
    
private:
    DBClientConnection mongo;
    std::string collection_name;
    bool by_tag;
};

#endif //__MONGOCONNECTOR_H__
