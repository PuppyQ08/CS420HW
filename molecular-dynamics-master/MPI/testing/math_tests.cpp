#include <iostream>
#include <cmath>
#include <string>

int get_grid_size(int num_processes) {
	double x = pow( double(num_processes), (1.0 / 3.0) );
	x += 0.1; 	// Make sure rounding error doesn't push it below the integer cube root
	return int(x);
}

void get_grid_coord(int rank, int num_processes, int* output) {
	int grid_len = get_grid_size(num_processes);
	output[2] = rank % grid_len;
	output[1] = (rank / grid_len) % grid_len;
	output[0] = rank / (num_processes / grid_len);
}

int main(int argc, char* argv[]) {
	int n = std::stoi(argv[1]);
	int m = std::stoi(argv[2]);
	std::cout << "Number of processes: " << n << std::endl;
	std::cout << "Grid size: " << get_grid_size(n) << std::endl;

	int grid_coord[3];
	get_grid_coord(m, n, grid_coord);
	std::cout << "Rank " << m << " has coordinates ";
	printf("(%d, %d, %d)\n", grid_coord[0], grid_coord[1], grid_coord[2]);
	
	return 0;
}



