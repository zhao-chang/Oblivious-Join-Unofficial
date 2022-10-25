#ifndef __SCHEMA_H__
#define __SCHEMA_H__

#include "Config.h"

//Block format: uint32_t blockID; uint32_t real_item_num; Item 1; Item 2; ... ; Dummy 1; Dummy 2; ... (no empty item)
//Item format: char flag ('r' or 'd' means "real" or "dummy"); item content ...

const uint32_t META_BLOCK_SIZE = 2 * sizeof(uint32_t);

typedef enum _ATTR_TYPE {
    CHAR,       //1 byte
    INTEGER,    //4 bytes
    DOUBLE,     //8 bytes
    STRING,     //with the same specific length
    TINYTEXT    //with variable length (no larger than a given maximum length)
} ATTR_TYPE;

struct CMP {
    uint32_t nCMPs;
    uint32_t attrID[MAX_CMPS];
    uint8_t order[MAX_CMPS];  //0: ascending; 1: descending
};

struct Schema {
    uint32_t nAttrs;
    uint32_t item_size;
    uint32_t item_per_blk;
    uint32_t attrName[MAX_COLS]; // Table ID (higher 16 bits) + Attribute ID (lower 16 bits)
    uint32_t attrOffset[MAX_COLS];
    uint32_t attrSize[MAX_COLS];
    ATTR_TYPE attrType[MAX_COLS];
};

#endif //__SCHEMA_H__
