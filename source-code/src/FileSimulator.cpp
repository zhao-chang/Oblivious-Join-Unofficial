#include "FileSimulator.h"
#include <filesystem>
namespace fs = std::filesystem;

FileSimulator::FileSimulator(const std::string& url, const std::string& collection_name, const bool create) 
        :url(url), collection_name(collection_name) {
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
        block_num = (uint64_t)ftell(fp) / (uint64_t)B;
        fclose(fp);
    }
}

FileSimulator::~FileSimulator() {
    remove(collection_name.c_str());
}

std::string FileSimulator::getCollectionName() {
    return collection_name;
}

// remove and create new file
void FileSimulator::clear(const std::string& ns) {
    remove(collection_name.c_str());

    FILE* fp = fopen(collection_name.c_str(), "wb");
    if (fp == NULL) {
        printf("Error: file cannot be created successfully!");
        return;
    }
    fclose(fp);
}

void FileSimulator::resize(const uint32_t& len, const std::string& ns) {
    fs::path fp = collection_name.c_str();
    uint64_t file_size = (uint64_t)len * (uint64_t)B;
    fs::resize_file(fp, file_size);
}

std::string FileSimulator::find(const uint32_t& id, const std::string& ns) {
    char block[B];
    FILE* fp = fopen(collection_name.c_str(), "rb");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return NULL;
    }
    if (id != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)id * (uint64_t)B;
        fseek(fp, filePos, SEEK_SET);
    }
    fread(block, 1, B, fp);
    fclose(fp);
    return std::string(block, B);
}

std::string FileSimulator::find(const std::string& id, const std::string& ns) {
    uint32_t block_id = atoi(id.c_str());
    return find(block_id, ns);
}

void FileSimulator::find(const std::vector<uint32_t>& ids, std::string* sbuffer, size_t& length, const std::string& ns) {
    length = ids.size();
    for (uint32_t i = 0; i < length; ++i) {
        sbuffer[i] = find(ids[i], ns);
    }
}

void FileSimulator::find(const uint32_t& low, const uint32_t& high, std::vector<std::string>& blocks, const std::string& ns) {
    char block[B];
    FILE* fp = fopen(collection_name.c_str(), "rb");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)B;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fread(block, 1, B, fp);
        blocks.push_back(std::string(block, B));
    }
    fclose(fp);
}

void FileSimulator::find(const uint32_t& low, const uint32_t& high, char* blocks, const std::string& ns) {
    char block[B];
    FILE* fp = fopen(collection_name.c_str(), "rb");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)B;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    char* pblock = blocks;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fread(block, 1, B, fp);
        memcpy(pblock, block, B);
        pblock += B;
    }
    fclose(fp);
}

std::string FileSimulator::fetch(const uint32_t& id, const std::string& ns) {
    //TODO:
    return find(id, ns);
}

std::string FileSimulator::fetch(const std::string& id, const std::string& ns) {
    //TODO:
    return find(id, ns);
}

void FileSimulator::fetch(const uint32_t& low, const uint32_t& high, char* blocks, const std::string& ns) {
    //TODO:
    find(low, high, blocks, ns);
}

void FileSimulator::update(const uint32_t& id, const std::string& data, const std::string& ns) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }
    if (id != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)id * (uint64_t)B;
        fseek(fp, filePos, SEEK_SET);
    }
    fwrite(data.c_str(), 1, B, fp);
    fclose(fp);
}

void FileSimulator::update(const std::string& id, const std::string& data, const std::string& ns) {
    uint32_t block_id = atoi(id.c_str());
    update(block_id, data, ns);
}

void FileSimulator::update(const std::string* sbuffer, const uint32_t& low, const size_t& len, const std::string& ns) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)B;
        fseek(fp, filePos, SEEK_SET);
    }

    for (uint32_t i = 0; i < len; ++i) {
        fwrite(sbuffer[i].c_str(), 1, B, fp);
    }
    fclose(fp);
}

void FileSimulator::update(const std::vector<std::string>& blocks, const uint32_t& low, const uint32_t& high, const std::string& ns) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)B;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fwrite(blocks[i].c_str(), 1, B, fp);
    }
    fclose(fp);
}

void FileSimulator::update(const char* blocks, const uint32_t& low, const uint32_t& high, const std::string& ns) {
    FILE* fp = fopen(collection_name.c_str(), "rb+");
    if (fp == NULL) {
        printf("Error: file cannot be opened successfully!");
        return;
    }

    assert(high >= low);
    if (low != 0) {
        //TODO: fseek_64 for Windows Platform
        uint64_t filePos = (uint64_t)low * (uint64_t)B;
        fseek(fp, filePos, SEEK_SET);
    }

    uint32_t block_cnt = high - low + 1;
    const char* pblock = blocks;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        fwrite(pblock, 1, B, fp);
        pblock += B;
    }
    fclose(fp);
}

void FileSimulator::update(const std::vector< std::pair<uint32_t, std::string> > blocks, const std::string& ns) {
    for (uint32_t i = 0; i < blocks.size(); ++i) {
        update(blocks[i].first, blocks[i].second, ns);
    }
}

void FileSimulator::insert(const uint32_t& id, const std::string& encrypted_block, const std::string& ns) {
    if (id < block_num)
        update(id, encrypted_block, ns);
    else {
        char block[B];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = id - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, B, fp);
        fwrite(encrypted_block.c_str(), 1, B, fp);
        fclose(fp);
        block_num = id + 1;
    }
}

void FileSimulator::insert(const std::vector< std::pair<uint32_t, std::string> >& blocks, const std::string& ns) {
    for (uint32_t i = 0; i < blocks.size(); ++i) {
        insert(blocks[i].first, blocks[i].second, ns);
    }
}

void FileSimulator::insert(const std::string* sbuffer, const uint32_t& low, const size_t& len, const std::string& ns) {
    if (low + len <= block_num) {
        update(sbuffer, low, len, ns);
    }
    else if (low >= block_num) {
        char block[B];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = low - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, B, fp);
        for (uint32_t i = 0; i < len; ++i)
            fwrite(sbuffer[i].c_str(), 1, B, fp);
        fclose(fp);
        block_num = low + len;
    }
    else {
        uint32_t update_len = block_num - low;
        update(sbuffer, low, update_len, ns);

        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        for (uint32_t i = update_len; i < len; ++i)
            fwrite(sbuffer[i].c_str(), 1, B, fp);
        fclose(fp);
        block_num = low + len;
    }
}

void FileSimulator::insert(const std::vector<std::string>& blocks, const uint32_t& low, const uint32_t& high, const std::string& ns) {
    if (high < block_num) {
        update(blocks, low, high, ns);
    }
    else if (low >= block_num) {
        char block[B];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = low - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, B, fp);
        uint32_t block_cnt = high - low + 1;
        for (uint32_t i = 0; i < block_cnt; ++i)
            fwrite(blocks[i].c_str(), 1, B, fp);
        fclose(fp);
        block_num = high + 1;
    }
    else {
        update(blocks, low, block_num - 1, ns);

        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t update_len = block_num - low;
        uint32_t block_cnt = high - low + 1;
        for (uint32_t i = update_len; i < block_cnt; ++i)
            fwrite(blocks[i].c_str(), 1, B, fp);
        fclose(fp);
        block_num = high + 1;
    }
}

void FileSimulator::insert(const char* blocks, const uint32_t& low, const uint32_t& high, const std::string& ns) {
    if (high < block_num) {
        update(blocks, low, high, ns);
    }
    else if (low >= block_num) {
        char block[B];
        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t append_num = low - block_num;
        for (uint32_t i = 0; i < append_num; ++i)
            fwrite(block, 1, B, fp);
        uint32_t block_cnt = high - low + 1;
        const char* pblock = blocks;
        for (uint32_t i = 0; i < block_cnt; ++i) {
            fwrite(pblock, 1, B, fp);
            pblock += B;
        }
        fclose(fp);
        block_num = high + 1;
    }
    else {
        update(blocks, low, block_num - 1, ns);

        FILE* fp = fopen(collection_name.c_str(), "ab");
        if (fp == NULL) {
            printf("Error: file cannot be opened successfully!");
            return;
        }
        uint32_t update_len = block_num - low;
        uint32_t block_cnt = high - low + 1;
        const char* pblock = blocks + update_len * B;
        for (uint32_t i = update_len; i < block_cnt; ++i) {
            fwrite(pblock, 1, B, fp);
            pblock += B;
        }
        fclose(fp);
        block_num = high + 1;
    }
}

void FileSimulator::insert(const std::vector< std::pair<std::string, std::string> >& blocks, const std::string& ns) {
    for (uint32_t i = 0; i < blocks.size(); ++i) {
        uint32_t block_id = atoi(blocks[i].first.c_str());
        insert(block_id, blocks[i].second, ns);
    }
}

void FileSimulator::remove(const uint32_t& id, const std::string& ns) {
    //TODO:
    find(id, ns);
}

void FileSimulator::remove(const std::string& id, const std::string& ns) {
    //TODO:
    find(id, ns);
}

void FileSimulator::remove(const uint32_t& low, const uint32_t& high, const std::string& ns) {
    //TODO:
    std::vector<std::string> blocks;
    find(low, high, blocks, ns);
}
