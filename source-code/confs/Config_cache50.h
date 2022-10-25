#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cstdint>
#include <string>

#define MAX_TBLS 10
#define MAX_COLS 35
#define MAX_CMPS 3

const uint32_t B = 4096;
const std::string server_host = "ravenserv1-i.cs.utah.edu";
//const std::string server_host = "localhost";

/* Path ORAM */
const uint32_t PathORAM_Z = 4;
const double btree_node_rate = 1.0;

/* Cache Reserved Factor */
const uint32_t threshold_amplifier = 50;

#endif //__CONFIG_H__
