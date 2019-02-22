#include <iostream>
#include <cmath>
int get_process_rank(int grid_len, int x, int y, int z) {
	// If one process needs to identify the rank of another process, based on grid coordinates
	return int(x * pow(grid_len, 2) + y * grid_len + z + 0.01);
}
int main(){
    std::cout<<get_process_rank(2, 0, 0, 0)<<std::endl;
    std::cout<<get_process_rank(2, 0, 0, 1)<<std::endl;
    std::cout<<get_process_rank(2, 1, 0, 0)<<std::endl;
    std::cout<<get_process_rank(2, 1, 0, 1)<<std::endl;
    std::cout<<get_process_rank(2, 0, 1, 0)<<std::endl;
    std::cout<<get_process_rank(2, 0, 1, 1)<<std::endl;
    std::cout<<get_process_rank(2, 1, 1, 0)<<std::endl;
    std::cout<<get_process_rank(2, 1, 1, 1)<<std::endl;
    return 0;
    }
