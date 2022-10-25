#ifndef __ENCRYPTION_ENGINE_H__
#define __ENCRYPTION_ENGINE_H__

#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <string>
#include <random>
#include <algorithm>

#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/sha.h>
#include <openssl/rand.h>
#include <openssl/hmac.h>
#include "Config.h"
#include "Util.h"

#define UNUSED(x) ((void)x)
const char hn[] = "SHA256";
const int aes_block_size = 16;
//const int key_size = 32;
const int key_size = 16;

class encryption_engine {
public:
    encryption_engine();
    encryption_engine(const encryption_engine& enc_engine);
    void encrypt(const std::string &plain_text, const byte* enc_key, std::string &cipher_text);
    void decrypt(const std::string &cipher_text, const byte* enc_key, std::string &plain_text);
    std::string hmac(const std::string &key);

private:
    void handle_errors();
    int sign_it(const byte* msg, size_t mlen, byte** sig, size_t* slen, EVP_PKEY* pkey);
    int encrypt(const unsigned char *plaintext, int plaintext_len, const unsigned char *key,
                unsigned char *ciphertext);
    int decrypt(const unsigned char *ciphertext, int ciphertext_len, const unsigned char *key,
                unsigned char *plaintext);
    int verify_it(const byte* msg, size_t mlen, const byte* sig, size_t slen, EVP_PKEY* pkey);
    void print_it(const char* label, const byte* buff, size_t len);
    int make_keys(EVP_PKEY** skey, EVP_PKEY** vkey);
    int hmac_it(const byte* msg, size_t mlen, byte** val, size_t* vlen, EVP_PKEY* pkey);
    
    EVP_PKEY * skey_, * vkey_;
};

#endif //__ENCRYPTION_ENGINE_H__
