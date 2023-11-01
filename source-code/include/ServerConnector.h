#ifndef SERVER_CONNECTOR_H
#define SERVER_CONNECTOR_H

#include <string>
#include <vector>

class ServerConnector {
public:
    ServerConnector() {}
    virtual ~ServerConnector() {}

    struct iterator {
        virtual ~iterator() {}
        virtual bool hasNext() = 0;
        virtual std::string next() = 0;
    };

    virtual void finalize() = 0;
    virtual void clear() = 0;
    virtual void resize(const uint32_t& len) = 0;
    virtual void insert(const uint32_t& id, const std::string& encrypted_block) = 0;
    virtual void insert(const std::vector< std::pair<uint32_t, std::string> >& blocks) = 0;
    virtual void insert(const std::string* sbuffer, const uint32_t& low, const size_t& len) = 0;
    virtual void insert(const char* blocks, const uint32_t& low, const uint32_t& high) = 0;
    virtual void insert(const std::vector<std::string>& blocks, const uint32_t& low, const uint32_t& high) = 0;
    virtual void insert(const std::vector< std::pair<std::string, std::string> >& blocks) = 0;
    virtual void insertWithTag(const std::vector< std::pair<std::string, std::string> >& blocks) = 0;
    virtual iterator* scan() = 0;
    virtual std::string find(const uint32_t& id) = 0;
    virtual std::string find(const std::string& id) = 0;
    virtual void find(const std::vector<uint32_t>& ids, std::string* sbuffer, size_t& length) = 0;
    virtual std::string fetch(const std::string& id) = 0;
    virtual std::string fetch(const uint32_t& id) = 0;
    virtual void fetch(const uint32_t& low, const uint32_t& high, char* blocks) = 0;
    virtual void find(const uint32_t& low, const uint32_t& high, std::vector<std::string>& blocks) = 0;
    virtual void find(const uint32_t& low, const uint32_t& high, char* blocks) = 0;
    virtual void findByTag(const uint32_t& tag, std::string* sbuffer, size_t& length) = 0;
    virtual void update(const uint32_t& id, const std::string& data) = 0;
    virtual void update(const std::string& id, const std::string& data) = 0;
    virtual void update(const std::string* sbuffer, const uint32_t& low, const size_t& len) = 0;
    virtual void update(const std::vector<std::string>& blocks, const uint32_t& low, const uint32_t& high) = 0;
    virtual void update(const char* blocks, const uint32_t& low, const uint32_t& high) = 0;
    virtual void update(const std::vector< std::pair<uint32_t, std::string> > blocks) = 0;
    virtual void remove(const uint32_t& id) = 0;
    virtual void remove(const std::string& id) = 0;
    virtual void remove(const uint32_t& low, const uint32_t& high) = 0;

    virtual std::string getCollectionName() = 0;
};

#endif //SERVER_CONNECTOR_H
