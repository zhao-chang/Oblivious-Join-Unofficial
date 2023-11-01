#ifndef __OBLIVIOUS_SORT_H__
#define __OBLIVIOUS_SORT_H__

#include "Schema.h"
#include "Util.h"
#include "ServerConnector.h"
#include "ObliviousHeap.h"
#include "Basic.h"

#include <algorithm>
#include <memory>
#include <cstring>

using namespace std;

const Schema* sort_sch = NULL;
const CMP* sort_cmp = NULL;

int sort_compare(const void* item1, const void* item2) {
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

class ObliviousSort {
public:
    ObliviousSort(ServerConnector* _conn, const uint32_t _n_block, const Schema* _sort_sch, const CMP* _sort_cmp) : n_block(_n_block) {
        conn = _conn;
        sort_sch = _sort_sch;
        sort_cmp = _sort_cmp;
    }
    
    ~ObliviousSort() {}
    
    uint32_t BatcherSort() {
        /*********/
        comm_size = 0;
        /*********/
        uint32_t len = (uint32_t)ceil(n_block / (double)m_block);
        if (len <= 2) internal_sort(0, 2);
        else {
            uint32_t padLen = 1 << (uint32_t)(ceil(log2(len)));
            oddeven_merge_sort_range(0, padLen - 1, len);
        }
        /*********/
        printf("finish Batcher sort successfully!\n");
        return comm_size;
        /*********/
    }
    
    uint32_t HeapSort() {
        /*********/
        comm_size = 0;
        /*********/
        ObliviousHeap* oheap = new ObliviousHeap(conn, n_block, sort_sch, sort_cmp);     
        char* src = buffer + META_BLOCK_SIZE;
        int32_t item_id = 0;
        for (uint32_t i = 0; i < n_block; ++i) {
            readBlock(conn, i, buffer);
            char* cur = src;
            for (uint32_t j = 0; j < sort_sch->item_per_blk; ++j) {
                std::string content(cur, sort_sch->item_size);
                oheap->insertItem(item_id, content);
                ++item_id;
                cur += sort_sch->item_size;
            }
        }
        
        for (uint32_t i = 0; i < n_block; ++i) {
            char* cur = src;
            uint32_t real_count = 0;
            for (uint32_t j = 0; j < sort_sch->item_per_blk; ++j) {
                std::string min_item = oheap->extractMin();
                memcpy(cur, min_item.c_str(), sort_sch->item_size);
                cur += sort_sch->item_size;
                if (min_item[0] == 'r')
                    ++real_count;
            }
            memcpy(buffer, &i, sizeof(uint32_t));
            memcpy(&(buffer[sizeof(uint32_t)]), &real_count, sizeof(uint32_t));
            updateBlock(conn, i, buffer);
        }
        comm_size += oheap->getCommSize();
        delete oheap;
        
        /*********/
        printf("finish heap sort successfully!\n");
        return comm_size;
        /*********/
    }
    
private:
    void internal_sort(const uint32_t lo, const uint32_t m_block_num) {
        uint32_t left = lo * m_block;
        uint32_t right = std::min((lo + m_block_num) * m_block, n_block);
        assert(left < right);
        uint32_t block_num = right - left;
        uint32_t item_num = block_num * sort_sch->item_per_blk;
        
        //read block
        readBlock(conn, left, right, buffer);
        /*********/
        comm_size += (size_t)(right - left) * (size_t)B;
        /*********/
        
        //temporarily keep meta-data for each block during sorting items
        char* meta = new char[block_num * META_BLOCK_SIZE];
        char* src = buffer;
        char* dst = meta;
        for (uint32_t i = 0; i < block_num; ++i) {
            memcpy(dst, src, META_BLOCK_SIZE);
            if (i < block_num - 1) {
                src += B;
                dst += META_BLOCK_SIZE;
            }
        }
        
        //compact
        src = buffer + META_BLOCK_SIZE;
        dst = buffer;
        uint32_t tot_item_size = sort_sch->item_per_blk * sort_sch->item_size;
        for (uint32_t i = 0; i < block_num; ++i) {
            for (uint32_t j = 0; j < tot_item_size; ++j) {
                char* pSrc = src + j;
                char* pDst = dst + j;
                *pDst = *pSrc;
            }
            if (i < block_num - 1) {
                src += B;
                dst += tot_item_size;
            }
        }
        
        //sort
        qsort(buffer, item_num, sort_sch->item_size, sort_compare);
        
        //recover items
        char* tmp = src;
        src = dst;
        dst = tmp;
        for (uint32_t i = 0; i < block_num; ++i) {
            for (int32_t j = tot_item_size - 1; j >= 0; --j) {
                char* pSrc = src + j;
                char* pDst = dst + j;
                *pDst = *pSrc;
            }
            if (i < block_num - 1) {
                src -= tot_item_size;
                dst -= B;
            }
        }
        
        //recover block meta-data
        src = meta;
        dst = buffer;
        for (uint32_t i = 0; i < block_num; ++i) {
            memcpy(dst, src, META_BLOCK_SIZE);
            if (i < block_num - 1) {
                src += META_BLOCK_SIZE;
                dst += B;
            }
        }
        delete[] meta;
        
        //update real_item_num
        src = buffer;
        for (uint32_t i = 0; i < block_num; ++i) {
            dst = src + META_BLOCK_SIZE;
            uint32_t real_count = 0;
            for (uint32_t j = 0; j < sort_sch->item_per_blk; ++j) {
                char flag = *((char *)dst);
                if (flag == 'r')
                    ++real_count;
                
                if (j < sort_sch->item_per_blk - 1)
                    dst += sort_sch->item_size;
            }
            memcpy(&(src[sizeof(uint32_t)]), &real_count, sizeof(uint32_t));
            if (i < block_num - 1)
                src += B;
        }
        
        // write block
        updateBlock(conn, left, right, buffer);
        /*********/
        comm_size += (size_t)(right - left) * (size_t)B;
        /*********/
    }
    
    void merge(const uint32_t first, const uint32_t second) {
        uint32_t fLeft = first * m_block;
        uint32_t fRight = fLeft + m_block;
        uint32_t sLeft = second * m_block;
        uint32_t sRight = std::min(sLeft + m_block, n_block);
        assert(sLeft < sRight);
        
        // read block
        char* fp = buffer;
        char* sp = buffer + m_block * B;
        readBlock(conn, fLeft, fRight, fp);
        readBlock(conn, sLeft, sRight, sp);
        /*********/
        comm_size += (size_t)((fRight - fLeft) + (sRight - sLeft)) * (size_t)B;
        /*********/
        
        uint32_t fLen = m_block;
        uint32_t sLen = sRight - sLeft;
        uint32_t totLen = fLen + sLen;
        uint32_t ftotItem = fLen * sort_sch->item_per_blk;
        uint32_t totItem = totLen * sort_sch->item_per_blk;        
        
        fLeft = 0;
        fRight = m_block;
        sLeft = m_block;
        sRight = m_block + sLen;
        
        // keep the final order: the item in the i-th slot should be put into the itemIndex[i]-th position after permutation
        uint32_t* itemIndex = new uint32_t[totItem];
        
        uint32_t sorted = 0;
        uint32_t fItem = 0, sItem = 0;
        uint32_t fItemInd = 0, sItemInd = ftotItem;
        fp += META_BLOCK_SIZE;
        sp += META_BLOCK_SIZE;
        char* fpCur = fp;
        char* spCur = sp;
        while (fLeft < fRight && sLeft < sRight) {
            // compare
            int32_t sort_cmp_res = sort_compare(fpCur, spCur);
            if (sort_cmp_res < 0) {
                itemIndex[fItemInd] = sorted;
                ++sorted;
                
                ++fItem;
                ++fItemInd;
                if (fItem >= sort_sch->item_per_blk) {
                    ++fLeft;
                    fItem = 0;
                    if (fLeft < fRight) {
                        fp += B;
                        fpCur = fp;
                    }
                }
                else fpCur += sort_sch->item_size;
            }
            else if (sort_cmp_res > 0) {
                itemIndex[sItemInd] = sorted;
                ++sorted;
                
                ++sItem;
                ++sItemInd;
                if (sItem >= sort_sch->item_per_blk) {
                    ++sLeft;
                    sItem = 0;
                    if (sLeft < sRight) {
                        sp += B;
                        spCur = sp;
                    }
                }
                else spCur += sort_sch->item_size;
            }
            else {
                itemIndex[fItemInd] = sorted;
                ++sorted;
                itemIndex[sItemInd] = sorted;
                ++sorted;
                
                ++fItem;
                ++fItemInd;
                if (fItem >= sort_sch->item_per_blk) {
                    ++fLeft;
                    fItem = 0;
                    if (fLeft < fRight) {
                        fp += B;
                        fpCur = fp;
                    }
                }
                else fpCur += sort_sch->item_size;
                ++sItem;
                ++sItemInd;
                if (sItem >= sort_sch->item_per_blk) {
                    ++sLeft;
                    sItem = 0;
                    if (sLeft < sRight) {
                        sp += B;
                        spCur = sp;
                    }
                }
                else spCur += sort_sch->item_size;
            }
        }
        
        while (fLeft < fRight) {
            while (fItem < sort_sch->item_per_blk) {
                itemIndex[fItemInd] = sorted;
                ++sorted;
                ++fItem;
                ++fItemInd;
            }
            fItem = 0;
            ++fLeft;
        }
        while (sLeft < sRight) {
            while (sItem < sort_sch->item_per_blk) {
                itemIndex[sItemInd] = sorted;
                ++sorted;
                ++sItem;
                ++sItemInd;
            }
            sItem = 0;
            ++sLeft;
        }
        
        assert(sorted == totItem);
        
        // permute based on the final order
        uint32_t item_cnt = 0;
        char* dst = buffer + META_BLOCK_SIZE;
        for (uint32_t i = 0; i < totLen; ++i) {
            char* dstCur = dst;
            for (uint32_t j = 0; j < sort_sch->item_per_blk; ++j) {
                uint32_t item_pos = itemIndex[item_cnt];
                while (item_cnt != item_pos) {
                    assert(item_cnt < item_pos);
                    uint32_t block_index = item_pos / sort_sch->item_per_blk;
                    uint32_t item_index = item_pos % sort_sch->item_per_blk;
                    char* srcCur = buffer + (block_index * B + META_BLOCK_SIZE + item_index * sort_sch->item_size);
                    
                    // swap
                    char tmp_item[sort_sch->item_size];
                    memcpy(tmp_item, srcCur, sort_sch->item_size);
                    memcpy(srcCur, dstCur, sort_sch->item_size);
                    memcpy(dstCur, tmp_item, sort_sch->item_size);
                    
                    itemIndex[item_cnt] = itemIndex[item_pos];
                    itemIndex[item_pos] = item_pos;
                    item_pos = itemIndex[item_cnt];
                }
                ++item_cnt;
                
                if (j < sort_sch->item_per_blk - 1)
                    dstCur += sort_sch->item_size;
            }
            if (i < totLen - 1)
                dst += B;
        }
        delete[] itemIndex;
        
        //update real_item_num
        char* src = buffer;
        for (uint32_t i = 0; i < totLen; ++i) {
            dst = src + META_BLOCK_SIZE;
            uint32_t real_count = 0;
            for (uint32_t j = 0; j < sort_sch->item_per_blk; ++j) {
                char flag = *((char *)dst);
                if (flag == 'r')
                    ++real_count;
                
                if (j < sort_sch->item_per_blk - 1)
                    dst += sort_sch->item_size;
            }
            memcpy(&(src[sizeof(uint32_t)]), &real_count, sizeof(uint32_t));
            if (i < totLen - 1)
                src += B;
        }
        
        // write block
        fLeft = first * m_block;
        fRight = fLeft + m_block;
        sLeft = second * m_block;
        sRight = std::min(sLeft + m_block, n_block);
        fp = buffer;
        sp = buffer + (m_block * B);
        updateBlock(conn, fLeft, fRight, fp);
        updateBlock(conn, sLeft, sRight, sp);
        /*********/
        comm_size += (size_t)((fRight - fLeft) + (sRight - sLeft)) * (size_t)B;
        /*********/
    }
    
    void oddeven_merge(const uint32_t lo, const uint32_t hi, const uint32_t r, const uint32_t len) {
        uint32_t step = r * 2;
        if (step < hi - lo) {
            oddeven_merge(lo, hi, step, len);
            oddeven_merge(lo + r, hi, step, len);
            for (uint32_t i = lo + r; i < hi - r; i += step)
                if (i + r < len) merge(i, i + r);
        }
        else if (lo < len) {
            // first operation for each "m_block": sort items in "m_block"s
            if (hi - lo == 1 && r == 1) internal_sort(lo, 2);
            // subsequent operations: each "m_block" has already been sorted; only need to merge
            else if (lo + r < len) merge(lo, lo + r);
        }
    }
    
    void oddeven_merge_sort_range(const uint32_t lo, const uint32_t hi, const uint32_t len) {
        if (hi - lo >= 1) {
            uint32_t mid = lo + (hi - lo) / 2;
            oddeven_merge_sort_range(lo, mid, len);
            oddeven_merge_sort_range(mid + 1, hi, len);
            oddeven_merge(lo, hi, 1, len);
        }
    }
    
    ServerConnector* conn;
    const uint32_t n_block;
    size_t comm_size;
};

#endif //__OBLIVIOUS_SORT_H__
