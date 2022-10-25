#include <cmath>
#include <vector>
#include <chrono>
#include <time.h>
#include "OneORAMIndexJoin.h"
#include "PathORAM.h"

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("Usage: ./one_type1_inlj_band_join <mode> <scale> <outsourced_height> <query_file>\n");
        return -1;
    }
    
    const uint32_t mode = atoi(argv[1]);
    const std::string method = "1bj";
    std::string scale = std::string(argv[2]);
    const uint32_t outsourced_height = atoi(argv[3]);
    std::string prefix = getCollectionPrefix(mode, method, scale, outsourced_height);
    std::string meta_prefix = "metas/" + prefix;
    OneORAMIndexJoin <PathORAM>* oij = new OneORAMIndexJoin <PathORAM> (mode, meta_prefix);
    
    FILE* fp = fopen(argv[4], "r");
    uint32_t n_tables;
    fscanf(fp, "%d", &n_tables);
    uint32_t table_id[n_tables];
    int32_t parent_id[n_tables];
    int32_t attr_id[n_tables];
    int32_t parent_attr_id[n_tables];
    for (uint32_t i = 0; i < n_tables; ++i)
        fscanf(fp, "%d", table_id + i);
    for (uint32_t i = 0; i < n_tables; ++i)
        fscanf(fp, "%d", parent_id + i);
    for (uint32_t i = 0; i < n_tables; ++i)
        fscanf(fp, "%d", attr_id + i);
    for (uint32_t i = 0; i < n_tables; ++i)
        fscanf(fp, "%d", parent_attr_id + i);
    
    double band_range[2];
    for (uint32_t i = 0; i < 2; ++i)
        fscanf(fp, "%lf", band_range + i);
    
    uint32_t n_proj_cols;
    fscanf(fp, "%d", &n_proj_cols);
    uint32_t proj_id[2][MAX_COLS];
    for (uint32_t i = 0; i < n_proj_cols; ++i)
        fscanf(fp, "%d %d", proj_id[0] + i, proj_id[1] + i);
    fclose(fp);
    
    auto start = std::chrono::high_resolution_clock::now();
    oij->ObliINLBJ(n_tables, table_id, parent_id, attr_id, parent_attr_id, band_range, n_proj_cols, proj_id);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    auto time = diff.count();
    printf("\n----------------------------------------\n");
    printf("Running Time: %.6lf\n", time);
    size_t server_size = oij->getServerSize();
    printf("Server Size is %lf MB\n", (double)server_size / (1 << 20));
    size_t client_size = oij->getClientSize();
    printf("Client Size is %lf MB\n", (double)client_size / (1 << 20));
    size_t comm_size = oij->getCommSize();
    printf("Comm Size is %lf MB\n", (double)comm_size / (1 << 20));
    printf("Access Count: %lu\n", oij->getAccessCount());
    printf("Read Time: %.6lf\n", oij->getReadTime());
    printf("Write Time: %.6lf\n", oij->getWriteTime());
    printf("EncDec Time: %.6lf\n", oij->getEncDecTime());
    printf("ORAM Time: %.6lf\n", oij->getORAMTime());
    printf("BTree Time: %.6lf\n", oij->getBTreeTime());
    printf("Res Time: %.6lf\n", oij->getResTime());
    printf("\n----------------------------------------\n");
    
    delete oij;
    return 0;
}
