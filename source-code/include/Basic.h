#ifndef __BASIC_H__
#define __BASIC_H__

#include "Util.h"
#include "FileSimulator.h"
#include "EncryptionEngine.h"
#include <string>

// The length of plaintext in each block
const uint32_t plain_len = B - aes_block_size;

// Additional trusted memory size M (for oblivious external sort)
uint32_t two_m_block = 1 << 10;
uint32_t m_block = two_m_block >> 1;
uint32_t oblisort_mem = two_m_block * B;

// Encryption engine
encryption_engine enc_engine;

// Encryption key
byte enc_key[key_size];

// The pointer of trusted memory
char* buffer = NULL;

void initKey(const char* fname, const bool create) {
    if (create) {
        std::string enc_str = generate_random_block(key_size);
        memcpy(enc_key, enc_str.c_str(), key_size);
        
        FILE* fp = fopen(fname, "w");
        fprintf(fp, "%d", enc_key[0]);
        for (uint32_t i = 1; i < key_size; ++i)
            fprintf(fp, " %d", enc_key[i]);
        fprintf(fp, "\n");
        fclose(fp);
    }
    else {
        FILE* fp = fopen(fname, "r");
        printf("enc_key:");
        for (uint32_t i = 0; i < key_size; ++i) {
            fscanf(fp, "%d", enc_key + i);
            printf(" %d", enc_key[i]);
        }
        printf("\n");
        fclose(fp);
    }
}

void initBufferPara(const uint32_t block_num) {
    uint32_t cache_size = threshold_amplifier * log2((double)block_num);
    m_block = cache_size >> 1;
    two_m_block = m_block << 1;
    oblisort_mem = two_m_block * B;
}

void initBufferSize(const uint32_t block_num) {
    uint32_t cache_size = block_num;
    if (cache_size < 2) cache_size = 2;
    m_block = cache_size >> 1;
    two_m_block = m_block << 1;
    oblisort_mem = two_m_block * B;
}

void initBuffer() {
    buffer = new char[oblisort_mem];
}

void destroyBuffer() {
    if (buffer != NULL) {
        delete[] buffer;
        buffer = NULL;
    }
}

void readBlock(ServerConnector* conn, const uint32_t blockID, char* block) {
    std::string cipher = conn->find(blockID);
    std::string plain;
    enc_engine.decrypt(cipher, enc_key, plain);
    memcpy(block, plain.c_str(), plain_len);
}

void fetchBlock(ServerConnector* conn, const uint32_t blockID, char* block) {
    std::string cipher = conn->fetch(blockID);
    std::string plain;
    enc_engine.decrypt(cipher, enc_key, plain);
    memcpy(block, plain.c_str(), plain_len);
}

void insertBlock(ServerConnector* conn, const uint32_t blockID, const char* block) {
    std::string plain = std::string(block, plain_len);
    std::string cipher;
    enc_engine.encrypt(plain, enc_key, cipher);
    conn->insert(blockID, cipher);
}

void updateBlock(ServerConnector* conn, const uint32_t blockID, const char* block) {
    std::string plain = std::string(block, plain_len);
    std::string cipher;
    enc_engine.encrypt(plain, enc_key, cipher);
    conn->update(blockID, cipher);
}

void removeBlock(ServerConnector* conn, const uint32_t blockID) {
    conn->remove(blockID);
}

void readBlock(ServerConnector* conn, const uint32_t lowID, const uint32_t highID, char* buffer) {
    assert(lowID < highID);
    conn->find(lowID, highID - 1, buffer);
    char* pBuf = buffer;
    uint32_t block_cnt = highID - lowID;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        std::string cipher = std::string(pBuf, B);
        std::string plain;
        enc_engine.decrypt(cipher, enc_key, plain);
        memcpy(pBuf, plain.c_str(), plain_len);
        pBuf += B;
    }
}

void fetchBlock(ServerConnector* conn, const uint32_t lowID, const uint32_t highID, char* buffer) {
    assert(lowID < highID);
    conn->fetch(lowID, highID - 1, buffer);
    char* pBuf = buffer;
    uint32_t block_cnt = highID - lowID;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        std::string cipher = std::string(pBuf, B);
        std::string plain;
        enc_engine.decrypt(cipher, enc_key, plain);
        memcpy(pBuf, plain.c_str(), plain_len);
        pBuf += B;
    }
}

void insertBlock(ServerConnector* conn, const uint32_t lowID, const uint32_t highID, char* buffer) {
    assert(lowID < highID);
    char* pBuf = buffer;
    uint32_t block_cnt = highID - lowID;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        std::string plain = std::string(pBuf, plain_len);
        std::string cipher;
        enc_engine.encrypt(plain, enc_key, cipher);
        memcpy(pBuf, cipher.c_str(), B);
        pBuf += B;
    }
    conn->insert(buffer, lowID, highID - 1);
}

void updateBlock(ServerConnector* conn, const uint32_t lowID, const uint32_t highID, char* buffer) {
    assert(lowID < highID);
    char* pBuf = buffer;
    uint32_t block_cnt = highID - lowID;
    for (uint32_t i = 0; i < block_cnt; ++i) {
        std::string plain = std::string(pBuf, plain_len);
        std::string cipher;
        enc_engine.encrypt(plain, enc_key, cipher);
        memcpy(pBuf, cipher.c_str(), B);
        pBuf += B;
    }
    conn->update(buffer, lowID, highID - 1);
}

void resizeFile(ServerConnector* conn, const uint32_t highID) {
    conn->resize(highID);
}

#endif //__BASIC_H__
