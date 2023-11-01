#include "FileSimulator.h"
#include <filesystem>
namespace fs = std::filesystem;

FileSimulator::FileSimulator(const std::string& url, const std::string& collection_name, const bool create, const uint32_t block_size)
        :url(url), collection_name(collection_name), block_size(block_size) {
    if (create) {
        FILE* fp = fopen(collection_name.c_str(), "wb");
        if (fp == NULL) {
            printf("Error: file cannot be created successfully!");
            return;
        }
        fclose(fp);
        block_num = 0;
    }
    else {
        FILE* fp = fopen(collection_name.c_str(), "rb");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        //TODO: fseek_64 for Windows Platform
        fseek(fp, 0, SEEK_END);
        block_num = (uint64_t)ftell(fp) / (uint64_t)block_size;
        fclose(fp);
    }
}

FileSimulator::~FileSimulator() {
    remove(collection_name.c_str());
}

std::string FileSimulator::getCollectionName() {
    return collection_name;
}

// remove file
void FileSimulator::finalize() {
    remove(collection_name.c_str());
}

// remove and create new file
void FileSimulator::clear() {
    remove(collection_name.c_str());

    FILE* fp = fopen(collection_name.c_str(), "wb");
    if (fp == NULL) {
        printf("Error: file cannot be created successfully!");
        return;
    }
    fclose(fp);
}

void FileSimulator::resize(const uint32_t& len) {
    fs::path fp = collection_name.c_str();
    uint64_t file_size = (uint64_t)len * (uint64_t)block_size;
    fs::resize_file(fp, file_size);
}

std::string FileSimulator::find(const uint32_t& id) {
    char block[block_size];
    FILE* fp = fopen(collection_name.c_str(), "rb");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return NULL;
    }
    if (id != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)id * (uint64_t)block_size;
        fseek(fp, filePos, SEEK_SET);
    }
    fread(block, 1, block_size, fp);
    fclose(fp);
    return std::string(block, block_size);
}

std::string FileSimulator::find(const std::string& id) {
    uint32_t block_id = atoi(id.c_str());
    return find(block_id);
}

void FileSimulator::find(const std::vector<uint32_t>& ids, std::string* sbuffer, size_t& length) {
    length = ids.size();
    for (uint32_t i = 0; i < length; ++i) {
        sbuffer[i] = find(ids[i]);
    }
}

void FileSimulator::find(const uint32_t& low, const uint32_t& high, std::vector<std::string>& blocks) {
    char block[block_size];
    FILE* fp = fopen(collection_name.c_str(), "rb");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)block_size;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fread(block, 1, block_size, fp);
        blocks.push_back(std::string(block, block_size));
    }
    fclose(fp);
}

void FileSimulator::find(const uint32_t& low, const uint32_t& high, char* blocks) {
    char block[block_size];
    FILE* fp = fopen(collection_name.c_str(), "rb");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)block_size;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    char* pblock = blocks;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fread(block, 1, block_size, fp);
        memcpy(pblock, block, block_size);
        pblock += block_size;
    }
    fclose(fp);
}

std::string FileSimulator::fetch(const uint32_t& id) {
    //TODO:
    return find(id);
}

std::string FileSimulator::fetch(const std::string& id) {
    //TODO:
    return find(id);
}

void FileSimulator::fetch(const uint32_t& low, const uint32_t& high, char* blocks) {
    //TODO:
    find(low, high, blocks);
}

void FileSimulator::update(const uint32_t& id, const std::string& data) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }
    if (id != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)id * (uint64_t)block_size;
        fseek(fp, filePos, SEEK_SET);
    }
    fwrite(data.c_str(), 1, block_size, fp);
    fclose(fp);
}

void FileSimulator::update(const std::string& id, const std::string& data) {
    uint32_t block_id = atoi(id.c_str());
    update(block_id, data);
}

void FileSimulator::update(const std::string* sbuffer, const uint32_t& low, const size_t& len) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)block_size;
        fseek(fp, filePos, SEEK_SET);
    }

    for (uint32_t i = 0; i < len; ++i) {
        fwrite(sbuffer[i].c_str(), 1, block_size, fp);
    }
    fclose(fp);
}

void FileSimulator::update(const std::vector<std::string>& blocks, const uint32_t& low, const uint32_t& high) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)block_size;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fwrite(blocks[i].c_str(), 1, block_size, fp);
    }
    fclose(fp);
}

void FileSimulator::update(const char* blocks, const uint32_t& low, const uint32_t& high) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)block_size;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    const char* pblock = blocks;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fwrite(pblock, 1, block_size, fp);
        pblock += block_size;
    }
    fclose(fp);
}

void FileSimulator::update(const std::vector< std::pair<uint32_t, std::string> > blocks) {
    for (uint32_t i = 0; i < blocks.size(); ++i) {
        update(blocks[i].first, blocks[i].second);
    }
}

void FileSimulator::insert(const uint32_t& id, const std::string& encrypted_block) {
    if (id < block_num)
        update(id, encrypted_block);
    else {
        char block[block_size];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = id - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, block_size, fp);
        fwrite(encrypted_block.c_str(), 1, block_size, fp);
        fclose(fp);
        block_num = id + 1;
    }
}

void FileSimulator::insert(const std::vector< std::pair<uint32_t, std::string> >& blocks) {
    for (uint32_t i = 0; i < blocks.size(); ++i) {
        insert(blocks[i].first, blocks[i].second);
    }
}

void FileSimulator::insert(const std::string* sbuffer, const uint32_t& low, const size_t& len) {
    if (low + len <= block_num) {
        update(sbuffer, low, len);
    }
    else if (low >= block_num) {
        char block[block_size];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = low - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, block_size, fp);
        for (uint32_t i = 0; i < len; ++i)
            fwrite(sbuffer[i].c_str(), 1, block_size, fp);
        fclose(fp);
        block_num = low + len;
    }
    else {
        uint32_t update_len = block_num - low;
        update(sbuffer, low, update_len);

        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        for (uint32_t i = update_len; i < len; ++i)
            fwrite(sbuffer[i].c_str(), 1, block_size, fp);
        fclose(fp);
        block_num = low + len;
    }
}

void FileSimulator::insert(const std::vector<std::string>& blocks, const uint32_t& low, const uint32_t& high) {
    if (high < block_num) {
        update(blocks, low, high);
    }
    else if (low >= block_num) {
        char block[block_size];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = low - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, block_size, fp); 
        uint32_t block_cnt = high - low + 1;
        for (uint32_t i = 0; i < block_cnt; ++i)
            fwrite(blocks[i].c_str(), 1, block_size, fp);
        fclose(fp);
        block_num = high + 1;
    }
    else {
        update(blocks, low, block_num - 1);
        
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t update_len = block_num - low;
        uint32_t block_cnt = high - low + 1;
        for (uint32_t i = update_len; i < block_cnt; ++i)
            fwrite(blocks[i].c_str(), 1, block_size, fp);
        fclose(fp);
        block_num = high + 1;
    }
}

void FileSimulator::insert(const char* blocks, const uint32_t& low, const uint32_t& high) {
    if (high < block_num) {
        update(blocks, low, high);
    }
    else if (low >= block_num) {
        char block[block_size];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = low - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, block_size, fp);
        uint32_t block_cnt = high - low + 1;
        const char* pblock = blocks;
        for (uint32_t i = 0; i < block_cnt; ++i) {
            fwrite(pblock, 1, block_size, fp);
            pblock += block_size;
        }
        fclose(fp);
        block_num = high + 1;
    }
    else {
        update(blocks, low, block_num - 1);

        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t update_len = block_num - low;
        uint32_t block_cnt = high - low + 1;
        const char* pblock = blocks + update_len * block_size;
        for (uint32_t i = update_len; i < block_cnt; ++i) {
            fwrite(pblock, 1, block_size, fp);
            pblock += block_size;
        }
        fclose(fp);
        block_num = high + 1;
    }
}

void FileSimulator::insert(const std::vector< std::pair<std::string, std::string> >& blocks) {
    for (uint32_t i = 0; i < blocks.size(); ++i) {
        uint32_t block_id = atoi(blocks[i].first.c_str());
        insert(block_id, blocks[i].second);
    }
}

void FileSimulator::remove(const uint32_t& id) {
    //TODO:
    find(id);
}

void FileSimulator::remove(const std::string& id) {
    //TODO:
    find(id);
}

void FileSimulator::remove(const uint32_t& low, const uint32_t& high) {
    //TODO:
    std::vector<std::string> blocks;
    find(low, high, blocks);
}
