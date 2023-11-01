#ifndef __OBLIVIOUS_COMPACT_H__
#define __OBLIVIOUS_COMPACT_H__

#include "Schema.h"
#include "Util.h"
#include "ServerConnector.h"
#include "Basic.h"

#include <algorithm>
#include <memory>
#include <cstring>

using namespace std;

const Schema* compact_sch = NULL;

class ObliviousCompact {
public:
    ObliviousCompact(ServerConnector* _conn, const uint32_t _n_block, const Schema* _compact_sch) : n_block(_n_block) {
        conn = _conn;
        compact_sch = _compact_sch;
    }
    
    ~ObliviousCompact() {}
    
    uint32_t Compact() {
        /*********/
        access_num = 0;
        /*********/
        if (n_block >= 1) {
            data_consolidation();
            
            if (n_block > 1) {
                uint32_t max_level = (uint32_t)ceil(log2(n_block));
                uint32_t diff_level = (uint32_t)(log2(two_m_block));
                
                uint32_t begin_level = 0;
                uint32_t end_level;
                while (begin_level < max_level) {
                    end_level = begin_level + diff_level;
                    if (end_level > max_level) end_level = max_level;
                    data_compaction(begin_level, end_level, end_level >= max_level);
                    begin_level = end_level;
                }
            }
        }
        /*********/
        printf("finish oblivious compaction successfully!\n");
        return access_num;
        /*********/
    }
    
private:
    void data_consolidation() {
        char oItem[compact_sch->item_size];
        oItem[0] = 'd';
        uint32_t oItemContentSize = compact_sch->item_size - sizeof(char);
        std::string rnd_str = generate_random_block(oItemContentSize);
        memcpy(oItem + sizeof(char), rnd_str.c_str(), oItemContentSize);
        
        char out_block[B];
        char tmp_block[B];
        memset(out_block, 0, META_BLOCK_SIZE);
        
        int32_t empty_pos = 0;
        uint32_t round = (uint32_t)ceil((double)n_block / two_m_block);
        for (uint32_t h = 0; h < round; ++h) {
            //read block
            uint32_t left = h * two_m_block;
            uint32_t right = std::min((h + 1) * two_m_block, n_block);
            assert(left < right);
            uint32_t block_num = right - left;
            readBlock(conn, left, right, buffer);
            
            char* pCur = buffer;
            for (uint32_t i = 0; i < block_num; ++i) {
                int32_t out_real_count;
                memcpy(&out_real_count, out_block + sizeof(int32_t), sizeof(int32_t));
                int32_t real_count;
                memcpy(&real_count, pCur + sizeof(uint32_t), sizeof(uint32_t));
                
                bool full = false;
                char* srcCur = pCur + META_BLOCK_SIZE;
                char* dstCur = NULL;
                if (out_real_count >= compact_sch->item_per_blk) {
                    full = true;
                    dstCur = pCur + META_BLOCK_SIZE;
                }
                else {
                    full = false;
                    dstCur = out_block + META_BLOCK_SIZE + out_real_count * compact_sch->item_size;
                }
                
                for (uint32_t j = 0; j < compact_sch->item_per_blk; ++j) {
                    char flag = *srcCur;
                    if (flag == 'r') {
                        if(dstCur != srcCur) {
                            memcpy(dstCur, srcCur, compact_sch->item_size);
                            memcpy(srcCur, oItem, compact_sch->item_size);
                        }
                        dstCur += compact_sch->item_size;
                        if (!full) {
                            ++out_real_count;
                            --real_count;
                            if (out_real_count >= compact_sch->item_per_blk) {
                                full = true;
                                dstCur = pCur + META_BLOCK_SIZE;
                            }
                        }
                    }
                    srcCur += compact_sch->item_size;
                }
                
                if (full) {
                    memcpy(pCur + sizeof(uint32_t), &real_count, sizeof(uint32_t));
                    memcpy(out_block + sizeof(uint32_t), &out_real_count, sizeof(uint32_t));
                    if (!(h == 0 && i == 0)) {
                        // swap
                        memcpy(tmp_block, out_block, B);
                        memcpy(out_block, pCur, B);
                        memcpy(pCur, tmp_block, B);
                        
                        memcpy(pCur, &empty_pos, sizeof(int32_t));
                    }
                }
                else {
                    assert(real_count == 0);
                    memcpy(pCur + sizeof(uint32_t), &real_count, sizeof(uint32_t));
                    memcpy(out_block + sizeof(uint32_t), &out_real_count, sizeof(uint32_t));
                    if (!(h == 0 && i == 0)) {
                        int32_t dummy_pos = -1;
                        memcpy(pCur, &dummy_pos, sizeof(int32_t));
                        ++empty_pos;
                    }
                }
                
                if (i < block_num - 1)
                    pCur += B;
            }
            
            // write block
            if (h == 0) {
                if (left < right - 1)
                    updateBlock(conn, left, right - 1, buffer + B);
            }
            else updateBlock(conn, left - 1, right - 1, buffer);
            if (h == round - 1) {
                memcpy(out_block, &empty_pos, sizeof(int32_t));
                updateBlock(conn, right - 1, out_block);
            }
        }
        /*********/
        access_num += 2 * n_block;
        /*********/
    }
    
    void data_compaction(uint32_t begin_level, uint32_t end_level, bool last) {
        char oItem[compact_sch->item_size];
        oItem[0] = 'd';
        uint32_t oItemContentSize = compact_sch->item_size - sizeof(char);
        std::string rnd_str = generate_random_block(oItemContentSize);
        memcpy(oItem + sizeof(char), rnd_str.c_str(), oItemContentSize);
        
        uint32_t min_skip = 1 << begin_level;
        uint32_t max_skip = 1 << end_level;
        for (uint32_t h = 0; h < min_skip; ++h) {
            uint32_t begin_pos = 0;
            int32_t block_cnt = 0;
            uint32_t begin_id = h;
            uint32_t end_id = begin_id;
            
            char* pEnd = buffer;
            uint32_t bound_id = std::min(h + max_skip, n_block);
            for (; end_id < bound_id; end_id += min_skip) {
                readBlock(conn, end_id, pEnd);
                ++block_cnt;
                pEnd += B;
            }
            assert(block_cnt <= two_m_block);
            if (begin_pos + block_cnt == two_m_block)
                pEnd = buffer;
            
            while (begin_id < n_block) {
                bool find = false;
                char* pCur = buffer + begin_pos * B;
                for (uint32_t j = 0; j < block_cnt; ++j) {
                    if (begin_pos + j == two_m_block)
                        pCur = buffer;
                    
                    int32_t empty_pos;
                    memcpy(&empty_pos, pCur, sizeof(int32_t));
                    if (empty_pos >= 0) {
                        assert(empty_pos % min_skip == 0);
                        uint32_t move = empty_pos % max_skip;
                        if (move == j * min_skip) {
                        // TODO:
                        // if (move == j) {
                            find = true;
                            if (last) empty_pos = begin_id;
                            else empty_pos -= move;
                            memcpy(pCur, &empty_pos, sizeof(int32_t));
                            updateBlock(conn, begin_id, pCur);
                            
                            // mark the moved block
                            empty_pos = -1;
                            memcpy(pCur, &empty_pos, sizeof(int32_t));
                            break;
                        }
                    }
                    pCur += B;
                }
                if (!find) {
                    char dummy_block[B];
                    int32_t dummy_pos = -1;
                    if (last) dummy_pos = begin_id;
                    memcpy(dummy_block, &dummy_pos, sizeof(int32_t));
                    int32_t real_count = 0;
                    memcpy(dummy_block + sizeof(int32_t), &real_count, sizeof(int32_t));
                    char* pItem = dummy_block + META_BLOCK_SIZE;
                    for (uint32_t j = 0; j < compact_sch->item_per_blk; ++j) {
                        memcpy(pItem, oItem, compact_sch->item_size);
                        pItem += compact_sch->item_size;
                    }
                    updateBlock(conn, begin_id, dummy_block);
                }
                begin_id += min_skip;
                --block_cnt;
                ++begin_pos;
                if (begin_pos >= two_m_block)
                    begin_pos = 0;
                
                if (end_id < n_block) {
                    readBlock(conn, end_id, pEnd);
                    end_id += min_skip;
                    ++block_cnt;
                    pEnd += B;
                    if (begin_pos + block_cnt == two_m_block)
                        pEnd = buffer;
                }
            }
            assert(block_cnt == 0);
        }
        /*********/
        access_num += 2 * n_block;
        /*********/
    }
    
    ServerConnector* conn;
    const uint32_t n_block;
    uint32_t access_num;
};

#endif //__OBLIVIOUS_COMPACT_H__
