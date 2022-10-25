#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdint.h>
#include <string>
#include <random>

typedef unsigned char byte;

uint32_t rand_uint32(const uint32_t &min, const uint32_t &max);
double rand_real(const double &min, const double &max);
size_t rand_int(const size_t& n);
std::string generate_random_block(const size_t& length);

#endif //__UTIL_H__
