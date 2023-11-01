#include "ObliviousHeap.h"
#include <cmath>
#include <cstring>
#include <chrono>
#include <time.h>
#include <assert.h>
#include <iostream>

ObliviousHeap::ObliviousHeap(ServerConnector* _conn, const uint32_t _n_block, const Schema* _sort_sch, const CMP* _sort_cmp)
        :n_input_blocks(_n_block) {
    sort_sch = _sort_sch;
    sort_cmp = _sort_cmp;
    
    n_data_items = n_input_blocks * sort_sch->item_per_blk;
    
    timestamp = 0;
    type_hiding_security = false;
    
    height = (uint32_t)ceil(log2((double)n_data_items)) + 1;
    n_items = (uint32_t)1 << (height - 1);
    stash.clear();
    sbuffer.clear();
    
    enc_item_size = sort_sch->item_size + 2 * sizeof(int64_t) + aes_block_size;
    std::string dst = "oram.oheap";
    conn = new FileSimulator(server_host, dst, true, enc_item_size);
    
    key = new byte[key_size];
    std::string rnd_str = generate_random_block(key_size);
    memcpy(key, rnd_str.c_str(), key_size);
    
    uint32_t n_buckets = (1 << height) - 1;
    for (uint32_t i = 0; i < n_buckets; ++i) {
        for (uint32_t j = 0; j < ObliviousHeap_Z; ++j) {
            std::string dummy_item = genDummyItem();
            std::string dummy_cipher;
            engine.encrypt(dummy_item, key, dummy_cipher);
            sbuffer.emplace_back(dummy_cipher);
        }
        conn->insert(sbuffer, i * ObliviousHeap_Z, (i + 1) * ObliviousHeap_Z - 1);
        sbuffer.clear();
    }
    /***********************/
    printf("Building ORAM with # of items: %u\n\n", n_items);
    
    n_server_items = ObliviousHeap_Z * ((1 << height) - 1);
    comm_size = 0;
    max_stash_size = 0;
    /***********************/
}

ObliviousHeap::~ObliviousHeap() {
    delete[] key;
    conn->finalize();
    delete conn;
}

std::pair<int32_t, int64_t> ObliviousHeap::insertItem(int32_t item_id, std::string content) {
    if (type_hiding_security) {
        onlyFindMin();
        onlyDeleteItem(-1, -1);
    }
    auto ret = std::pair<int32_t, int64_t>(-1, -1);
    ret = onlyInsertItem(item_id, content);
    ++timestamp;
    return ret;
}

std::string ObliviousHeap::deleteItem(int32_t pos, int64_t ts) {
    if (type_hiding_security)
        onlyFindMin();
    std::string ret = onlyDeleteItem(pos, ts);
    if (type_hiding_security)
        insertDummyItem();
    ++timestamp;
    return ret;
}

std::string ObliviousHeap::findMin() {
    std::string min_item = onlyFindMin();
    std::string ret = min_item.substr(2 * sizeof(int64_t));
    if (type_hiding_security) {
        onlyDeleteItem(-1, -1);
        insertDummyItem();
    }
    return ret;
}

std::string ObliviousHeap::extractMin() {
    std::string min_item = onlyFindMin();
    int32_t pos;
    int64_t timestamp;
    memcpy(&pos, min_item.substr(sizeof(int32_t)).c_str(), sizeof(int32_t));
    memcpy(&timestamp, min_item.substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
    std::string ret = onlyDeleteItem(pos, timestamp);
    if (type_hiding_security)
        insertDummyItem();
    ++timestamp;
    return ret;
}

/***********************
void ObliviousHeap::printItem(const char* pItem) {
    if (*pItem == 'r') {
        uint32_t attrLen = sort_sch->nAttrs;
        for (uint32_t indexL = 0; indexL < attrLen; ++indexL) {
            uint32_t offset = sort_sch->attrOffset[indexL];
            ATTR_TYPE attrType = sort_sch->attrType[indexL];
            uint32_t attrSize = sort_sch->attrSize[indexL];
            
            const char* pAttr = pItem + offset;
            if (attrType == CHAR) {
                char attrValue = *pAttr;
                printf("%c ", attrValue);
            }
            else if (attrType == INTEGER) {
                int32_t attrValue;
                memcpy(&attrValue, pAttr, sizeof(int32_t));
                printf("%d ", attrValue);
            }
            else if (attrType == DOUBLE) {
                double attrValue;
                memcpy(&attrValue, pAttr, sizeof(double));
                printf("%lf ", attrValue);
            }
            else if (attrType == STRING || attrType == TINYTEXT) {
                char attrValue[attrSize + 1];
                memcpy(attrValue, pAttr, attrSize);
                attrValue[attrSize] = '\0';
                printf("%s ", attrValue);
            }
        }
        printf("\n");
    }
    else printf("dummy\n");
}
/***********************/

/***********************
void ObliviousHeap::print() {
    std::cout << "==============================" << std::endl;
    std::cout << "n_data_items: " << n_data_items << std::endl;
    uint32_t n_buckets = (1 << height) - 1;
    for (uint32_t i = 0; i < n_buckets; ++i) {
        std::vector <std::string> sbuffer;
        readBucket(i, sbuffer);
        //std::cout << "-----------" << std::endl;
        for (int32_t j = 0; j < ObliviousHeap_Z; ++j) {
            int32_t item_id;
            memcpy(&item_id, sbuffer[j].c_str(), sizeof(int32_t));
            int32_t item_pos;
            memcpy(&item_pos, sbuffer[j].substr(sizeof(int32_t)).c_str(), sizeof(int32_t));
            int64_t item_timestamp;
            memcpy(&item_timestamp, sbuffer[j].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
            if (item_id != -1 && item_timestamp != -1) {
                std::string item_content = sbuffer[j].substr(2 * sizeof(int64_t));
                std::cout << "bucket_id: " << i << ", block: " << j  << std::endl;
                std::cout << "item_id: " << item_id << ", item_pos: " << item_pos << ", item_timestamp: " << item_timestamp << std::endl;
                std::cout << "item_content: ";
                printItem(item_content.c_str());
            }
        }
        sbuffer.clear();
        //std::cout << "-----------" << std::endl;
    }
    std::cout << "-----------" << std::endl;
    for (auto j = stash.begin(); j != stash.end(); ++j) {
        int32_t item_id = j->first;
        int32_t item_pos;
        memcpy(&item_pos, j->second.c_str(), sizeof(int32_t));
        int64_t item_timestamp;
        memcpy(&item_timestamp, j->second.substr(sizeof(int32_t)).c_str(), sizeof(int64_t));
        if (item_id != -1 && item_timestamp != -1) {
            std::string item_content = j->second.substr(sizeof(int32_t) + sizeof(int64_t));
            std::cout << "stash:" << std::endl;
            std::cout << "item_id: " << item_id << ", item_pos: " << item_pos << ", item_timestamp: " << item_timestamp << std::endl;
            std::cout << "item_content: ";
            printItem(item_content.c_str());
        }
    }
    std::cout << "-----------" << std::endl;
    std::cout << "==============================" << std::endl;
}
/***********************/

/***********************
void ObliviousHeap::check_pos() {
    uint32_t n_buckets = (1 << height) - 1;
    for (uint32_t i = 0; i < n_buckets; ++i) {
        std::vector <std::string> sbuffer;
        readBucket(i, sbuffer);
        for (int32_t j = 0; j < ObliviousHeap_Z; ++j) {
            int32_t item_id;
            memcpy(&item_id, sbuffer[j].c_str(), sizeof(int32_t));
            int32_t item_pos;
            memcpy(&item_pos, sbuffer[j].substr(sizeof(int32_t)).c_str(), sizeof(int32_t));
            int64_t item_timestamp;
            memcpy(&item_timestamp, sbuffer[j].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
            if (item_id != -1 && item_timestamp != -1) {
                int32_t level = (int32_t)log2(i + 1);
                if (i != P(item_pos, level)) {
                    std::string item_content = sbuffer[j].substr(2 * sizeof(int64_t));
                    std::cout << "==============================" << std::endl;
                    std::cout << "bucket_id: " << i << ", block: " << j  << std::endl;
                    std::cout << "item_id: " << item_id << ", item_pos: " << item_pos
                        << ", item_timestamp: " << item_timestamp << ", item_content: " << item_content << std::endl;
                    std::cout << "==============================" << std::endl;
                    
                    print();
                    abort();
                }
            }
        }
        sbuffer.clear();
    }
}
/***********************/

int ObliviousHeap::itemCompare(const void* item1, const void* item2) {
    int8_t* p1 = (int8_t *)item1;
    int8_t* p2 = (int8_t *)item2;
    if (*p1 != *p2) return *p2 - *p1;
    else if (*p1 == 'r' && *p2 == 'r') {
        for (uint32_t index = 0; index < sort_cmp->nCMPs; ++index) {
            uint32_t attrID = sort_cmp->attrID[index];
            ATTR_TYPE attrType = sort_sch->attrType[attrID];
            uint32_t offset = sort_sch->attrOffset[attrID];
            uint8_t order = sort_cmp->order[index];
            if (attrType == CHAR) {
                int8_t* char1 = p1 + offset;
                int8_t* char2 = p2 + offset;
                int8_t res = *char1 - *char2;
                if (res != 0) {
                    if (order == 0) return res;
                    else if (order == 1) return -res;
                }
            }
            else if (attrType == INTEGER) {
                int8_t* pos1 = p1 + offset;
                int8_t* pos2 = p2 + offset;
                int32_t val1;
                int32_t val2;
                memcpy(&val1, pos1, sizeof(int32_t));
                memcpy(&val2, pos2, sizeof(int32_t));
                if (val1 != val2) {
                    if (order == 0) return val1 - val2;
                    else if (order == 1) return val2 - val1;
                }
            }
            else if (attrType == DOUBLE) {
                int8_t* pos1 = p1 + offset;
                int8_t* pos2 = p2 + offset;
                double val1;
                double val2;
                memcpy(&val1, pos1, sizeof(double));
                memcpy(&val2, pos2, sizeof(double));
                if (val1 != val2) {
                    if (order == 0) {
                        if (val1 - val2 < 0.0) return -1;
                        else return 1;
                    }
                    else if (order == 1) {
                        if (val2 - val1 < 0.0) return -1;
                        else return 1;
                    }
                }
            }
            else if (attrType == STRING || attrType == TINYTEXT) {
                char* str1 = (char *)p1 + offset;
                char* str2 = (char *)p2 + offset;
                int32_t res = strncmp(str1, str2, sort_sch->attrSize[attrID]);
                if (res != 0) {
                    if (order == 0) return res;
                    else if (order == 1) return -res;
                }
            }
        }
        return 0;
    }
    return 0;
}

std::string ObliviousHeap::genDummyItem() {
    int32_t dummy_item_id = -1;
    int32_t dummy_item_pos = -1;
    int64_t dummy_item_ts = -1;
    std::string dummy_id_str = std::string((const char *)(& dummy_item_id), sizeof(int32_t));
    std::string dummy_pos_str = std::string((const char *)(& dummy_item_pos), sizeof(int32_t));
    std::string dummy_ts_str = std::string((const char *)(& dummy_item_ts), sizeof(int64_t));
    std::string dummy_content = generate_random_block(sort_sch->item_size);
    std::string dummy_item = dummy_id_str + dummy_pos_str + dummy_ts_str + dummy_content;
    return dummy_item;
}

int32_t ObliviousHeap::P(int32_t pos, int32_t level) {
    if (level == 0) return 0;
    return (1 << level) - 1 + (pos / (1 << (height - 1 - level)));
}

void ObliviousHeap::readBucket(const int32_t bucket_pos, std::vector<std::string> & sbuffer) {
    sbuffer.clear();
    int32_t begin = bucket_pos * ObliviousHeap_Z;
    int32_t end = begin + ObliviousHeap_Z - 1;
    conn->find(begin, end, sbuffer);
    comm_size += ObliviousHeap_Z;
    for (size_t i = 0; i < ObliviousHeap_Z; ++i) {
        std::string plain;
        engine.decrypt(sbuffer[i], key, plain);
        sbuffer[i] = plain;
    }
    assert(sbuffer.size() == ObliviousHeap_Z);
}

void ObliviousHeap::writeBucket(const int32_t bucket_pos, std::vector<std::string> & sbuffer) {
    assert(sbuffer.size() == ObliviousHeap_Z);
    for (size_t i = 0; i < ObliviousHeap_Z; ++i) {
        std::string cipher;
        engine.encrypt(sbuffer[i], key, cipher);
        sbuffer[i] = cipher;
    }
    int32_t begin = bucket_pos * ObliviousHeap_Z;
    int32_t end = begin + ObliviousHeap_Z - 1;
    conn->update(sbuffer, begin, end);
    comm_size += ObliviousHeap_Z;
    sbuffer.clear();
}

std::string ObliviousHeap::onlyFindMin() {
    readBucket(P(0, 0), sbuffer);
    std::string min_item = sbuffer[0];
    if (stash.size() > 0) {
        int64_t min_timestamp;
        memcpy(&min_timestamp, min_item.substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
        std::string min_content = min_item.substr(2 * sizeof(int64_t));
        for (auto j = stash.begin(); j != stash.end(); ++j) {
            int64_t stash_timestamp;
            memcpy(&stash_timestamp, j->second.substr(sizeof(int32_t)).c_str(), sizeof(int64_t));
            std::string stash_content = j->second.substr(sizeof(int32_t) + sizeof(int64_t));
            int cmp = itemCompare(min_content.c_str(), stash_content.c_str());
            if (cmp > 0 || cmp == 0 && min_timestamp > stash_timestamp) {
                std::string item_id = std::string((const char *)(&(j->first)), sizeof(uint32_t));
                min_item = item_id + j->second;
            }
        }
    }
    sbuffer.clear();
    return min_item;
}

std::pair<int32_t, int64_t> ObliviousHeap::onlyInsertItem(int32_t item_id, std::string content) {
    auto ret = std::pair<int32_t, int64_t>();
    if (item_id != -1) {
        int32_t pos = rand_int(n_items);
        std::string pos_str = std::string((const char *)(& pos), sizeof(int32_t));
        std::string timestamp_str = std::string((const char *)(& timestamp), sizeof(int64_t));
        stash[item_id] = pos_str + timestamp_str + content;
        ret.first = pos;
        ret.second = timestamp;
    }
    evictAfterInsert();
    return ret;
}

void ObliviousHeap::insertDummyItem() {
    evictAfterInsert();
}

void ObliviousHeap::evictAfterInsert() {
    int32_t random_pos_one = rand_int(n_items >> 1);
    int32_t random_pos_two = rand_int(n_items >> 1) + (n_items >> 1);
    evictAndUpdateMin(random_pos_one);
    evictAndUpdateMin(random_pos_two);
}

std::string ObliviousHeap::onlyDeleteItem(int32_t pos, int64_t timestamp) {
    auto j = stash.begin();
    for (; j != stash.end(); ++j) {
        int32_t stash_pos;
        int64_t stash_timestamp;
        memcpy(&stash_pos, j->second.c_str(), sizeof(int32_t));
        memcpy(&stash_timestamp, j->second.substr(sizeof(int32_t)).c_str(), sizeof(int64_t));
        if (stash_pos == pos && stash_timestamp == timestamp)
            break;
    }
    
    std::string stash_item;
    int64_t stash_timestamp = -1;
    if (j != stash.end()) {
        stash_item = j->second.substr(sizeof(int32_t) + sizeof(uint64_t));
        memcpy(&stash_timestamp, stash_item.substr(sizeof(int32_t)).c_str(), sizeof(int64_t));
        stash.erase(j);
    }
    
    std::string storage_item = readAndRemove(pos, timestamp);
    evictAndUpdateMin(pos);
    
    if (stash_timestamp != -1)
        return stash_item;
    else return storage_item;
}

void ObliviousHeap::evictAndUpdateMin(int32_t pos) {
    uint32_t access_pos;
    if (pos >= 0) access_pos = pos;
    else access_pos = rand_int(n_items);
    
    for (int32_t l = 0; l < height; ++l) {
        readBucket(P(access_pos, l), sbuffer);
        for (int32_t i = 1; i < ObliviousHeap_Z; ++i) {
            int32_t item_id;
            int64_t item_timestamp;
            memcpy(&item_id, sbuffer[i].c_str(), sizeof(int32_t));
            memcpy(&item_timestamp, sbuffer[i].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
            if (item_id != -1 && item_timestamp != -1)
                stash[item_id] = sbuffer[i].substr(sizeof(int32_t));
        }
        sbuffer.clear();
    }
    
    for (int32_t l = height - 1; l >= 0; --l) {
        sbuffer.emplace_back(genDummyItem());
        int32_t bucket_access_pos = P(access_pos, l);
        int32_t counter = 1;
        std::unordered_map<uint32_t, std::string>::iterator j, tmp;
        for (auto j = stash.begin(); j != stash.end() && counter < ObliviousHeap_Z;) {
            int32_t stash_access_pos;
            memcpy(&stash_access_pos, j->second.c_str(), sizeof(int32_t));
            if (bucket_access_pos == P(stash_access_pos, l)) {
                int32_t item_id = j->first;
                std::string item_id_str = std::string((const char *)(& item_id), sizeof(int32_t));
                sbuffer.emplace_back(item_id_str + j->second);
                tmp = j; ++j; stash.erase(tmp);
                ++counter;
            }
            else ++j;
        }
        while (counter < ObliviousHeap_Z) {
            sbuffer.emplace_back(genDummyItem());
            ++counter;
        }
        writeBucket(bucket_access_pos, sbuffer);
        sbuffer.clear();
    }
    
    updateMin(access_pos);
}

std::string ObliviousHeap::readAndRemove(int32_t pos, int64_t timestamp) {
    if (pos != -1 && timestamp != -1) {
        std::string ret;
        for (int32_t l = 0; l < height; ++l) {
            int32_t bucket_pos = P(pos, l);
            readBucket(bucket_pos, sbuffer);
            for (int32_t i = 1; i < ObliviousHeap_Z; ++i) {
                int64_t item_ts;
                memcpy(&item_ts, sbuffer[i].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
                if (item_ts == timestamp) {
                    ret = sbuffer[i].substr(2 * sizeof(int64_t));
                    sbuffer[i] = genDummyItem();
                    break;
                }
            }
            writeBucket(bucket_pos, sbuffer);
            sbuffer.clear();
        }
        return ret;
    }
    else {
        uint32_t rnd_pos = rand_int(n_items);
        for (int32_t l = 0; l < height; ++l) {
            int32_t bucket_pos = P(rnd_pos, l);
            readBucket(bucket_pos, sbuffer);
            writeBucket(bucket_pos, sbuffer);
            sbuffer.clear();
        }
        return std::string();
    }
}

void ObliviousHeap::updateMin(int32_t pos) {
    std::vector<std::string> child_buffer;
    if (pos != -1) {
        for (int32_t l = height - 1; l >= 0; --l) {
            int32_t bucket_pos = P(pos, l);
            readBucket(bucket_pos, sbuffer);
            
            std::string min_item = "";
            for (int32_t i = 1; i < ObliviousHeap_Z; ++i) {
                int32_t item_id;
                memcpy(&item_id, sbuffer[i].c_str(), sizeof(int32_t));
                int64_t item_timestamp;
                memcpy(&item_timestamp, sbuffer[i].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
                if (item_id != -1 && item_timestamp != -1) {
                    if (min_item.empty())
                        min_item = sbuffer[i];
                    else {
                        std::string min_content = min_item.substr(2 * sizeof(int64_t));
                        std::string item_content = sbuffer[i].substr(2 * sizeof(int64_t));
                        int cmp = itemCompare(min_content.c_str(), item_content.c_str());
                        int64_t min_timestamp;
                        memcpy(&min_timestamp, min_item.substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
                        if (cmp > 0 || cmp == 0 && min_timestamp > item_timestamp)
                            min_item = sbuffer[i];
                    }
                }
            }
            
            std::vector<std::string> potential_mins;
            if (!min_item.empty()) potential_mins.emplace_back(min_item);
            
            if (l != height - 1) {
                for (int32_t j = 1; j <= 2; ++j) {
                    readBucket((bucket_pos << 1) + j, child_buffer);
                    int32_t item_id;
                    memcpy(&item_id, child_buffer[0].c_str(), sizeof(int32_t));
                    int64_t item_timestamp;
                    memcpy(&item_timestamp, child_buffer[0].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
                    if (item_id != -1 && item_timestamp != -1)
                        potential_mins.push_back(child_buffer[0]);
                    child_buffer.clear();
                }
            }
            
            int32_t potential_size = potential_mins.size();
            if (potential_size > 0) {
                int32_t min_index = 0;
                int64_t min_timestamp;
                memcpy(&min_timestamp, potential_mins[0].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
                for (int32_t i = 1; i < potential_size; ++i) {
                    std::string min_content = potential_mins[min_index].substr(2 * sizeof(int64_t));
                    std::string potential_content = potential_mins[i].substr(2 * sizeof(int64_t));
                    int cmp = itemCompare(min_content.c_str(), potential_content.c_str());
                    int64_t potential_timestamp;
                    memcpy(&potential_timestamp, potential_mins[i].substr(sizeof(int64_t)).c_str(), sizeof(int64_t));
                    if (cmp > 0 || cmp == 0 && min_timestamp > potential_timestamp) {
                        min_index = i;
                        min_timestamp = potential_timestamp;
                    }
                }
                sbuffer[0] = potential_mins[min_index];
            }
            else sbuffer[0] = genDummyItem();
            writeBucket(bucket_pos, sbuffer);
            sbuffer.clear();
        }
    }
    else {
        uint32_t rnd_pos = rand_int(n_items);
        for (int32_t l = height - 1; l >= 0; --l) {
            int32_t bucket_pos = P(rnd_pos, l);
            readBucket(bucket_pos, sbuffer);
            if (l != height - 1) {
                for (int32_t j = 1; j <= 2; ++j) {
                    readBucket((bucket_pos << 1) + j, child_buffer);
                    child_buffer.clear();
                }
            }
            writeBucket(bucket_pos, sbuffer);
            sbuffer.clear();
        }
    }
}

/***********************/
uint32_t ObliviousHeap::getDataItemNum() const {
    return n_data_items;
}

size_t ObliviousHeap::getServerSize() const {
    return n_server_items * enc_item_size;
}

size_t ObliviousHeap::getClientSize() const {
    return max_stash_size * enc_item_size + 2 * enc_item_size * ObliviousHeap_Z;
}

size_t ObliviousHeap::getCommSize() const {
    size_t total_comm_size = comm_size * (size_t)enc_item_size;
    return total_comm_size;
}

void ObliviousHeap::resetCommSize() {
    comm_size = 0;
}
/***********************/
