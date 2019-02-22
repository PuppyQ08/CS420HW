#include <cmath>

inline double distance(double** positions, int a, int b) { 	// Calculates the straight line distance 
	double squared_dist = 0.0;								// between particles a and b
	for (int i = 0; i < 3; i++) { squared_dist += pow(positions[a][i] - positions[b][i], 2); }
	return sqrt(squared_dist);
}

inline double external_distance(double x1, double y1, double z1, double x2, double y2, double z2) {
	double dx = x2 - x1;
	double dy = y2 - y1;
	double dz = z2 - z1;
	return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}

inline double speed(double** velocity, int a) { // Calculates the speed of particle a
	double square_speed = 0.0;
	for (int i = 0; i < 3; i++) { square_speed += pow(velocity[a][i], 2); }
	return sqrt(square_speed);
}

inline double kinetic_energy(double** velocity, double m, int a) {
	double square_speed = 0.0;
	for (int i = 0; i < 3; i++) { square_speed += pow(velocity[a][i], 2); }
	return 0.5 * m * square_speed;
}

inline double LJ_potential(double r, double epsilon, double sigma) {
	/*
	Calculates the Leonard-Jones potential between two particles at a distance r, 
	epsilon is the depth of the well, and sigma is the distance at which the potential is 0.
	*/

	return 4.0 * epsilon * ( pow(sigma / r, 12) - pow(sigma / r, 6) );
}

inline double LJ_force(double r, double epsilon, double sigma ) { 
	/*
	Calculates the force between two particles at a distance r, resulting from a Leonard-Jones potential
	epsilon is the depth of the well, and sigma is the distance at which the potential is 0.

	x is just a temporary variable so the logic is clearer.
	*/
	double x = 2 * pow(sigma / r, 13) - pow(sigma / r, 7);  
	return 24 * epsilon * x / sigma;
}

inline int get_grid_dim(int num_processes) {
	double x = pow( double(num_processes), (1.0 / 3.0) );
	x += 0.1; 	// Make sure rounding error doesn't push it below the integer cube root
	return int(x);
}

void get_grid_coord(int rank, int num_processes, int* output) {
	int grid_len = get_grid_dim(num_processes);
	output[2] = rank % grid_len;
	output[1] = (rank / grid_len) % grid_len;
	output[0] = rank / (num_processes / grid_len);
}

int get_process_rank(int grid_len, int x, int y, int z) {
	// If one process needs to identify the rank of another process, based on grid coordinates
	return int(x * pow(grid_len, 2) + y * grid_len + z + 0.01);
}







