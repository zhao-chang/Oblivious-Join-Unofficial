#include <cmath>
#include <vector>
#include <chrono>
#include <time.h>
#include "ObliDatabaseJoin.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: ./obli_database_load <mode> <scale> \n");
        return -1;
    }
    
    const uint32_t mode = atoi(argv[1]);
    const std::string method = "odbj";
    std::string scale = std::string(argv[2]);
    std::string prefix = getCollectionPrefix(mode, method, scale);
    std::string dst = "oram." + prefix;
    std::string meta_prefix = "metas/" + prefix;
    
    auto start = std::chrono::high_resolution_clock::now();
    ObliDatabaseJoin* odbj = new ObliDatabaseJoin(mode, scale, dst, meta_prefix);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    auto time = diff.count();
    
    printf("\n----------------------------------------\n");
    printf("Loading Time: %.6lf\n", time);
    printf("Data Block Number: %u\n", odbj->getDataBlockNum());
    size_t server_size = odbj->getServerSize();
    printf("Server Size is %lf MB\n", (double)server_size / (1 << 20));
    size_t client_size = odbj->getClientSize();
    printf("Client Size is %lf MB\n", (double)client_size / (1 << 20));
    printf("\n----------------------------------------\n");
    
    delete odbj;
    return 0;
}
