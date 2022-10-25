#include <cmath>
#include <vector>
#include <chrono>
#include <time.h>
#include "ObliDatabaseJoin.h"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        printf("Usage: ./obli_database_equi_join <mode> <scale> <query_file>\n");
        return -1;
    }
    
    const uint32_t mode = atoi(argv[1]);
    const std::string method = "odbj";
    std::string scale = std::string(argv[2]);
    std::string prefix = getCollectionPrefix(mode, method, scale);
    std::string meta_prefix = "metas/" + prefix;
    ObliDatabaseJoin* odbj = new ObliDatabaseJoin(mode, meta_prefix);
    
    FILE* fp = fopen(argv[3], "r");
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
    
    uint32_t n_proj_cols;
    fscanf(fp, "%d", &n_proj_cols);
    uint32_t proj_id[2][MAX_COLS];
    for (uint32_t i = 0; i < n_proj_cols; ++i)
        fscanf(fp, "%d %d", proj_id[0] + i, proj_id[1] + i);
    fclose(fp);
    
    auto start = std::chrono::high_resolution_clock::now();
    odbj->ObliEquiJoin(n_tables, table_id, parent_id, attr_id, parent_attr_id, n_proj_cols, proj_id);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    auto time = diff.count();
    printf("\n----------------------------------------\n");
    printf("Running Time: %.6lf\n", time);
    size_t server_size = odbj->getServerSize();
    printf("Server Size is %lf MB\n", (double)server_size / (1 << 20));
    size_t client_size = odbj->getClientSize();
    printf("Client Size is %lf MB\n", (double)client_size / (1 << 20));
    size_t comm_size = odbj->getCommSize();
    printf("Comm Size is %lf MB\n", (double)comm_size / (1 << 20));
    printf("\n----------------------------------------\n");
    
    delete odbj;
    return 0;
}
