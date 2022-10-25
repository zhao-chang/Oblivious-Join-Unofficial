#include "Util.h"

uint32_t rand_uint32(const uint32_t &min, const uint32_t &max) {
    static thread_local std::random_device random_device_;
    static thread_local std::mt19937 generator_(random_device_());
    std::uniform_int_distribution<uint32_t> distribution(min, max);
    return distribution(generator_);
};

double rand_real(const double &min, const double &max) {
    static thread_local std::random_device random_device_;
    static thread_local std::mt19937 generator_(random_device_());
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(generator_);
};

size_t rand_int(const size_t& n) {
    static thread_local std::random_device random_device_;
    static thread_local std::mt19937 generator_(random_device_());
    return std::uniform_int_distribution<size_t>(0, n - 1)(generator_);
}

std::string generate_random_block(const size_t& length) {
    static const char alphanum[] = "0123456789";
    std::string ret;
    ret.resize(length);
    for (size_t i = 0; i < length; ++i) {
        ret[i] = alphanum[rand_uint32(0, 500) % (sizeof(alphanum) - 1)];
    }
    return ret;
};
