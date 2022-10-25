#include "PathORAM.h"
#include <cmath>
#include <cstring>
#include <chrono>
#include <time.h>

PathORAM::PathORAM(const uint32_t& n, const char* dst, const char* output, const bool& pos_host) {
    height = (uint32_t)ceil(log2((double)n)) + 1;
    n_blocks = (uint32_t)1 << (height - 1);
    stash.clear();
    if (pos_host) {
        pos_map = new std::unordered_map<uint32_t, uint32_t>;
        for (uint32_t i = 0; i < n; ++i)
            pos_map->insert(std::make_pair(i, rand_int(n_blocks)));
    }
    conn = new FileSimulator(server_host, dst, true);
    
    key = new byte[key_size];
    std::string rnd_str = generate_random_block(key_size);
    memcpy(key, rnd_str.c_str(), key_size);
    sbuffer = new std::string[height * PathORAM_Z];
    
    uint32_t n_buckets = (1 << height) - 1;
    for (uint32_t i = 0; i < n_buckets; ++i) {
        for (uint32_t j = 0; j < PathORAM_Z; ++j) {
            std::string tmp  = generate_random_block(B - aes_block_size - sizeof(uint32_t));
            int32_t dummyID = -1;
            std::string dID = std::string((const char *)(& dummyID), sizeof(uint32_t));
            std::string cipher;
            engine.encrypt(dID + tmp, key, cipher);
            conn->insert(i * PathORAM_Z + j, cipher);
        }
    }
    
    /***********************/
    printf("Building ORAM with # of blocks: %u\n\n", n_blocks);
    
    n_data_blocks = n;
    n_server_blocks = PathORAM_Z * ((1 << height) - 1);
    comm_size = 0;
    max_stash_size = 0;
    access_count = 0;
    
    read_time = 0.0;
    write_time = 0.0;
    enc_dec_time = 0.0;
    oram_time = 0.0;
    /***********************/
    
    FILE* fp = fopen(output, "w");
    fprintf(fp, "%s\n", dst);
    fprintf(fp, "%d\n", n);
    fprintf(fp, "%d\n", n_blocks);
    fprintf(fp, "%d\n", pos_host ? 1 : 0);
    if (pos_host) {
        for (auto it = pos_map->begin(); it != pos_map->end(); ++it)
            fprintf(fp, "%d %d\n", it->first, it->second);
    }
    fprintf(fp, "%d", key[0]);
    for (uint32_t i = 1; i < key_size; ++i)
        fprintf(fp, " %d", key[i]);
    fprintf(fp, "\n");
    fclose(fp);
}

PathORAM::PathORAM(std::unordered_map<uint32_t, std::string>& blocks, const char* dst, const char* output, const bool& pos_host) {
    uint32_t n = blocks.size();
    height = (uint32_t)ceil(log2((double)n)) + 1;
    n_blocks = (uint32_t)1 << (height - 1);
    stash.clear();
    conn = new FileSimulator(server_host, dst, true);
    
    key = new byte[key_size];
    std::string rnd_str = generate_random_block(key_size);
    memcpy(key, rnd_str.c_str(), key_size);
    sbuffer = new std::string[height * PathORAM_Z];
    
    if (pos_host) {
        pos_map = new std::unordered_map<uint32_t, uint32_t>;
        for (auto it = blocks.begin(); it != blocks.end(); ++it)
            pos_map->insert(std::make_pair(it->first, rand_int(n_blocks)));
    }
    bulkLoad(blocks);
    
    /***********************/
    printf("Building ORAM with # of blocks: %u\n\n", n_blocks);
    
    n_data_blocks = n;
    n_server_blocks = PathORAM_Z * ((1 << height) - 1);
    comm_size = 0;
    max_stash_size = 0;
    access_count = 0;
    
    read_time = 0.0;
    write_time = 0.0;
    enc_dec_time = 0.0;
    oram_time = 0.0;
    /***********************/
    
    FILE* fp = fopen(output, "w");
    fprintf(fp, "%s\n", dst);
    fprintf(fp, "%d\n", n);
    fprintf(fp, "%d\n", n_blocks);
    fprintf(fp, "%d\n", pos_host ? 1 : 0);
    if (pos_host) {
        for (auto it = pos_map->begin(); it != pos_map->end(); ++it)
            fprintf(fp, "%d %d\n", it->first, it->second);
    }
    fprintf(fp, "%d", key[0]);
    for (uint32_t i = 1; i < key_size; ++i)
        fprintf(fp, " %d", key[i]);
    fprintf(fp, "\n");
    fclose(fp);
}

PathORAM::PathORAM(const char* input) {
    FILE* fp = fopen(input, "r");
    char dst[100];
    fscanf(fp, "%s", dst);
    conn = new FileSimulator(server_host, dst, false);
    
    uint32_t n;
    fscanf(fp, "%d", &n);
    height = (uint32_t)ceil(log2((double)n)) + 1;
    n_blocks = (uint32_t)1 << (height - 1);
    fscanf(fp, "%d", &n_blocks);
    
    uint32_t mode;
    fscanf(fp, "%d", &mode);
    if (mode) {
        pos_map = new std::unordered_map<uint32_t, uint32_t>;
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t block_id, block_pos;
            fscanf(fp, "%d %d", &block_id, &block_pos);
            pos_map->insert(std::make_pair(block_id, block_pos));
        }
    }
    
    key = new byte[key_size];
    printf("key:");
    for (uint32_t i = 0; i < key_size; ++i) {
        fscanf(fp, "%d", key + i);
        printf(" %d", key[i]);
    }
    printf("\n");
    fclose(fp);
    
    stash.clear();
    sbuffer = new std::string[height * PathORAM_Z];
    
    /***********************/
    n_data_blocks = n;
    n_server_blocks = PathORAM_Z * ((1 << height) - 1);
    comm_size = 0;
    max_stash_size = 0;
    access_count = 0;
    
    read_time = 0.0;
    write_time = 0.0;
    enc_dec_time = 0.0;
    oram_time = 0.0;
    /***********************/
}

PathORAM::~PathORAM() {
    delete[] key;
    delete pos_map;
    delete[] sbuffer;
    delete conn;
}

std::string PathORAM::get(const std::string & key) {
    /***********************/
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    std::string res;
    uint32_t int_key;
    sscanf(key.c_str(), "%d", &int_key);
    access('r', int_key, res);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    oram_time += (double)(diff.count());
    /***********************/
    return res;
}

void PathORAM::put(const std::string & key, const std::string & value) {
    /***********************/
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    uint32_t int_key;
    sscanf(key.c_str(), "%d", &int_key);
    std::string value2 = value;
    access('w', int_key, value2);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    oram_time += (double)(diff.count());
    /***********************/
}

std::string PathORAM::get(const int32_t & key) {
    /***********************/
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    std::string res;
    access('r', key, res);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    oram_time += (double)(diff.count());
    /***********************/
    return res;
}

std::string PathORAM::getAndRemove(const int32_t & key, const uint32_t& pos) {
    /***********************/
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    std::string res;
    access('r', key, res, pos);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    oram_time += (double)(diff.count());
    /***********************/
    return res;
}

void PathORAM::put(const int32_t & key, const std::string & value) {
    /***********************/
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    std::string value2 = value;
    access('w', key, value2);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    oram_time += (double)(diff.count());
    /***********************/
}

void PathORAM::put(const int32_t & key, const std::string & value, const uint32_t& pos) {
    /***********************/
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    std::string value2 = value;
    access('w', key, value2, pos);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    oram_time += (double)(diff.count());
    /***********************/
}

void PathORAM::access(const char& op, const int32_t& block_id, std::string& data) {
    /***********************/
    ++access_count;
    /***********************/
    
    uint32_t x;
    if (block_id == -1) x = rand_int(n_blocks);
    else {
        x = pos_map->at(block_id);
        pos_map->at(block_id) = rand_int(n_blocks);
    }
    
    size_t length;
    fetchAlongPath(x, sbuffer, length);
    for (size_t i = 0; i < length; ++i) {
        std::string plain;
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        engine.decrypt(sbuffer[i], key, plain);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        enc_dec_time += (double)(diff.count());
        /***********************/
        
        int32_t b_id;
        memcpy(&b_id, plain.c_str(), sizeof(uint32_t));
        if (b_id != -1) {
            stash[b_id] = plain.substr(sizeof(uint32_t));
        }
    }
    /***********************/
    max_stash_size = std::max(stash.size(), max_stash_size);
    /***********************/
    
    if (block_id != -1) {
        if (op == 'r') data = stash[block_id];
        else stash[block_id] = data;
    }
    
    for (uint32_t i = 0; i < height; ++i) {
        uint32_t tot = 0;
        uint32_t base = i * PathORAM_Z;
        std::unordered_map<uint32_t, std::string>::iterator j, tmp;
        j = stash.begin();
        while (j != stash.end() && tot < PathORAM_Z) {
            if (check(pos_map->at(j->first), x, i)) {
                std::string b_id = std::string((const char *)(&(j->first)), sizeof(uint32_t));
                sbuffer[base + tot] = b_id + j->second;
                tmp = j; ++j; stash.erase(tmp);
                ++tot;
            } else ++j;
        }
        for (uint32_t k = tot; k < PathORAM_Z; ++k) {
            std::string tmp_block = generate_random_block(B - aes_block_size - sizeof(uint32_t));
            int32_t dummyID = -1;
            std::string dID = std::string((const char *)(& dummyID), sizeof(uint32_t));
            sbuffer[base + k] = dID + tmp_block;
        }
    }
    
    for (size_t i = 0; i < height * PathORAM_Z; ++i) {
        std::string cipher;
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        engine.encrypt(sbuffer[i], key, cipher);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        enc_dec_time += (double)(diff.count());
        /***********************/
        sbuffer[i] = cipher;
    }
    loadAlongPath(x, sbuffer);
}

void PathORAM::access(const char& op, const int32_t& block_id, std::string& data, const uint32_t& pos) {
    /***********************/
    ++access_count;
    /***********************/
    
    uint32_t x;
    if (block_id == -1) x = rand_int(n_blocks);
    else {
        if (op == 'r') x = pos;
        else x = rand_int(n_blocks);
    }
    
    size_t length;
    fetchAlongPath(x, sbuffer, length);
    for (size_t i = 0; i < length; ++i) {
        std::string plain;
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        engine.decrypt(sbuffer[i], key, plain);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        enc_dec_time += (double)(diff.count());
        /***********************/
        
        int32_t b_id;
        memcpy(&b_id, plain.c_str(), sizeof(uint32_t));
        if (b_id != -1) {
            stash[b_id] = plain.substr(sizeof(uint32_t));
        }
    }
    
    /***********************/
    max_stash_size = std::max(stash.size(), max_stash_size);
    /***********************/
    
    if (block_id != -1) {
        if (op == 'r') {
            data = stash[block_id];
            stash.erase(block_id);
        } else stash[block_id] = data;
    }
    
    for (uint32_t i = 0; i < height; ++i) {
        uint32_t tot = 0;
        uint32_t base = i * PathORAM_Z;
        std::unordered_map<uint32_t, std::string>::iterator j, tmp;
        j = stash.begin();
        while (j != stash.end() && tot < PathORAM_Z) {
            uint32_t tmp_pos;
            memcpy(&tmp_pos, j->second.c_str(), sizeof(uint32_t));
            if (check(tmp_pos, x, i)) {
                std::string b_id = std::string((const char *)(&(j->first)), sizeof(uint32_t));
                sbuffer[base + tot] = b_id + j->second;
                tmp = j; ++j; stash.erase(tmp);
                ++tot;
            } else ++j;
        }
        for (uint32_t k = tot; k < PathORAM_Z; ++k) {
            std::string tmp_block  = generate_random_block(B - aes_block_size - sizeof(uint32_t));
            int32_t dummyID = -1;
            std::string dID = std::string((const char *)(& dummyID), sizeof(uint32_t));
            sbuffer[base + k] = dID + tmp_block;
        }
    }
    
    for (size_t i = 0; i < height * PathORAM_Z; ++i) {
        std::string cipher;
        /***********************/
        auto start = std::chrono::high_resolution_clock::now();
        /***********************/
        engine.encrypt(sbuffer[i], key, cipher);
        /***********************/
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        enc_dec_time += (double)(diff.count());
        /***********************/
        sbuffer[i] = cipher;
    }
    loadAlongPath(x, sbuffer);
}

bool PathORAM::check(int x, int y, int l) {
    return (x >> l) == (y >> l);
}

void PathORAM::fetchAlongPath(const uint32_t& x, std::string* sbuffer, size_t& length) {
    uint32_t cur_pos = x + (1 << (height - 1));
    std::vector<uint32_t> ids;
    while (cur_pos > 0) {
        for (uint32_t i = 0; i < PathORAM_Z; ++i)
            ids.push_back((cur_pos - 1) * PathORAM_Z + i);
        cur_pos >>= 1;
    }
    /***********************/
    comm_size += ids.size();
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    conn->find(ids, sbuffer, length);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    read_time += (double)(diff.count());
    /***********************/
}

void PathORAM::loadAlongPath(const uint32_t& x, const std::string* sbuffer) {
    uint32_t cur_pos = x + (1 << (height - 1));
    uint32_t offset = 0;
    insert_buffer.clear();
    while (cur_pos > 0) {
        for (uint32_t i = 0; i < PathORAM_Z; ++i)
            insert_buffer.emplace_back(std::make_pair((cur_pos - 1) * PathORAM_Z + i, sbuffer[offset + i]));
        offset += PathORAM_Z;
        cur_pos >>= 1;
    }
    /***********************/
    comm_size += insert_buffer.size();
    auto start = std::chrono::high_resolution_clock::now();
    /***********************/
    conn->update(insert_buffer);
    /***********************/
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    write_time += (double)(diff.count());
    /***********************/
}

//Note that after bulk loading, the content in `blocks` will be clear!!!!
void PathORAM::bulkLoad(std::unordered_map<uint32_t, std::string>& blocks) {
    //printf("Start Bulk Loading.....\n");
    bool flag = (pos_map != NULL);
    insert_buffer.clear();
    uint32_t q = 0;
    std::unordered_map<uint32_t, uint32_t> blk_cnt;
    n_server_blocks = 0;
    for (uint32_t i = 0; i < n_blocks; ++i) blk_cnt[i] = 0;
    for (uint32_t i = 0; i < height; ++i) {
        //printf("Processing Level %d, %lu blocks left.....\n", i, blocks.size());
        std::unordered_map<uint32_t, std::string>::iterator j, tmp;
        j = blocks.begin();
        
        uint32_t cnt = 0;
        while (j != blocks.end()) {
            uint32_t level_pos;
            if (flag) level_pos = pos_map->at(j->first) >> i;
            else {
                memcpy(&level_pos, j->second.c_str(), sizeof(uint32_t));
                level_pos >>= i;
            }
            if (blk_cnt[level_pos] < PathORAM_Z) {
                std::string b_id = std::string((const char *)(&(j->first)), sizeof(uint32_t));
                std::string cipher;
                engine.encrypt(b_id + j->second, key, cipher);
                uint32_t tree_pos = level_pos + (1 << (height - i - 1));
                insert_buffer.emplace_back(std::make_pair((tree_pos - 1) * PathORAM_Z + blk_cnt[level_pos], cipher));
                if (insert_buffer.size() >= 10000) {
                    n_server_blocks += insert_buffer.size();
                    conn->insert(insert_buffer);
                    insert_buffer.clear();
                }
                blk_cnt[level_pos]++;
                tmp = j; ++j; blocks.erase(tmp);
            } else ++j;
        }
        
        for (auto x : blk_cnt) {
            for (uint32_t k = x.second; k < PathORAM_Z; ++k) {
                std::string tmp_block  = generate_random_block(B - aes_block_size - sizeof(uint32_t));
                int32_t dummyID = -1;
                std::string dID = std::string((const char *)(& dummyID), sizeof(uint32_t));
                std::string cipher;
                engine.encrypt(dID + tmp_block, key, cipher);
                insert_buffer.emplace_back(std::make_pair((x.first + (1 << (height - i - 1)) - 1) * PathORAM_Z + k, cipher));
                if (insert_buffer.size() >= 10000) {
                    n_server_blocks += insert_buffer.size();
                    conn->insert(insert_buffer);
                    insert_buffer.clear();
                }
            }
        }
        
        blk_cnt.clear();
        for (uint32_t j = 0; j < (n_blocks >> (i + 1)); ++j) blk_cnt[j] = 0;
        n_server_blocks += insert_buffer.size();
        conn->insert(insert_buffer);
        insert_buffer.clear();
    }
    
    assert(blocks.empty());
    printf("Bulk loading finished... # blocks on server: %lu\n", n_server_blocks);
}

/***********************/
uint32_t PathORAM::getDataBlockNum() const {
    return n_data_blocks;
}

size_t PathORAM::getServerSize() const {
    return n_server_blocks * B;
}

size_t PathORAM::getClientSize() const {
    size_t pos_map_size = (pos_map == NULL) ? 0 : n_data_blocks * sizeof(uint32_t);
    return max_stash_size * B + pos_map_size + height * B * PathORAM_Z;
}

size_t PathORAM::getCommSize() const {
    return comm_size * B;
}

void PathORAM::resetCommSize() {
    comm_size = 0;
}

size_t PathORAM::getAccessCount() const {
    return access_count;
}

double PathORAM::getReadTime() const {
    return read_time;
}

double PathORAM::getWriteTime() const {
    return write_time;
}

double PathORAM::getEncDecTime() const {
    return enc_dec_time;
}

double PathORAM::getORAMTime() const {
    return oram_time;
}
/***********************/

// Batch access, if (blocks[i] == "") then read block from ORAM to unordered map,
// Otherwise, write the contents stored in blocks[i] to ORAM...

void PathORAM::batchAccess(std::unordered_map<uint32_t, std::string>& blocks) {
    std::vector<uint32_t> keys;
    for (auto x : blocks) keys.push_back(x.first);
    std::sort(keys.begin(), keys.end(), [this](uint32_t a, uint32_t b) { return this->pos_map->at(a) < this->pos_map->at(b); });
    for (auto key : keys) {
        if (blocks[key] == "") blocks[key] = get(key);
        else put(key, blocks[key]);
    }
}
