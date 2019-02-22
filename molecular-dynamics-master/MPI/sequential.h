#include <cmath>

inline double distance(double** positions, int a, int b) { 	// Calculates the straight line distance 
	double squared_dist = 0.0;								// between particles a and b
	for (int i = 0; i < 3; i++) { squared_dist += pow(positions[a][i] - positions[b][i], 2); }
	return sqrt(squared_dist);
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