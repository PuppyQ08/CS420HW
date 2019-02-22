#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <cmath>

#include "sequential.h"

using namespace std;

int main() {

	double t0 = double( clock() );	// start the clock (real time)
	int i, j, k;
	const int num_particles = 8192;
	double universe_size[3] = { 64.0, 64.0, 64.0 };	// Position bounds in x, y, z

	const double mass = 4.0;		// Helium
	const double kB = 8314.0;		// amu * nm^2 ns^-2 K^-1 
	const double V = pow(universe_size[0], 3);		// cubic nanometers

	double a = 0.0;					// Van-der-Waals parameters
	double b = 0.0;
	
	double speed_sum;
	double total_kinetic_energy;

	double t_initial = 0.0;				// (simulation time)
	double t_final = 0.01;	
	double t_step = 0.001;

	// to do: get inputs from a file


    //randomly generated coordinates of particles

	random_device rd;
	mt19937_64 rand(rd());	// 64 bit Mersenne Twister
	uniform_real_distribution<> generate_position(0.1, universe_size[0]-0.1);
	normal_distribution<> generate_velocity(0.0, 2000.0);
	
	double** positions = new double*[num_particles];
	double** velocities = new double*[num_particles];
	double** accelerations = new double*[num_particles];
	
	for (i = 0; i < num_particles; i++ ) {
		positions[i] = new double[3];
		velocities[i] = new double[3];
		accelerations[i] = new double[3];
	}

	for (i = 0; i < num_particles; i++) {
		for (j = 0; j < 3; j++) {
			positions[i][j] = generate_position(rand);
			velocities[i][j] = generate_velocity(rand);
			if (velocities[i][j] > 10000.0) { velocities[i][j] = 1000.0; }		// Make sure velocity isn't too high
			accelerations[i][j] = 0.0;
		}
	}

	double t1 = double(clock());
	cout << "Generated initial conditions in " << (t1 - t0) / CLOCKS_PER_SEC << " seconds.\n";

    /* writes the coordinates to output file at every time step
     * now is off. you can set it on if you like
     */
	bool output_everything = false;	


	// output files 
	ofstream f, v;
	f.open("MD_output.txt");
	v.open("average_speed.txt");

	// to do: format this more nicely
	f << "Initial positions:\n";
	for (i = 0; i < num_particles; i++) {
		f << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2];
		f << endl;
	}
	f << endl;

	double t = t_initial;

	while (t < t_final) {

        /* main part needed parallel*/
        
        speed_sum = 0.0;
		total_kinetic_energy = 0.0;

		for (i = 0; i < num_particles; i++) {

			double force[3] = { 0.0, 0.0, 0.0 };
			double prev_acceleration[3];

			// update positions
			for (j = 0; j < 3; j++) {
				positions[i][j] = positions[i][j] + velocities[i][j] * t_step + 0.5 * accelerations[i][j] * pow(t_step, 2);
				prev_acceleration[j] = accelerations[i][j];

				// Elastic collisions
				if (positions[i][j] < 0.0) { 

					if (positions[i][j] < -32.0) {
						positions[i][j] = 32.0;
					}

					else {
						positions[i][j] *= -1.0;
						velocities[i][j] *= -1.0;
					}
				}
				if (positions[i][j] > universe_size[j]) { 
					positions[i][j] = 2 * universe_size[j] - positions[i][j];
					velocities[i][j] *= -1.0;
				}
			}

			// Calculate forces
			for (j = 0; j < num_particles; j++) {
				if (i != j) {
					double r = distance(positions, i, j);

					if (r < 8.0) {	// cutoff radius

						double Fnorm = LJ_force(r, 0.4, 0.125);

						// multiply the norm of the force by each component of the unit vector
						for (k = 0; k < 3; k++) { force[k] += Fnorm * (positions[i][k] - positions[j][k]) / r; }
					}
				}
			}

			// calculate acceleration and velocity of particle i
			for (j = 0; j < 3; j++) { 
				accelerations[i][j] = force[j] / mass; 
				velocities[i][j] = velocities[i][j] + 0.5 * (prev_acceleration[j] + accelerations[i][j]) * t_step;
			}

			if (speed(velocities, i) > 20000.0) {
				// If the speed is too fast, scale down the speed while preserving the direction
				double scale_factor = speed(velocities, i) / 1000.0;
				velocities[i][0] = velocities[i][0] / scale_factor;
				velocities[i][1] = velocities[i][1] / scale_factor;
				velocities[i][2] = velocities[i][2] / scale_factor;
			}

			speed_sum += speed(velocities, i);
			total_kinetic_energy += kinetic_energy(velocities, mass, i);
		}

		if (output_everything) {
			f << "Step " << t << endl;
			for (j = 0; j < num_particles; j++) {
				f << positions[j][0] << ", " << positions[j][1] << ", " << positions[j][2];
				f << endl;
			}
			f << endl;
		}
		t += t_step;

		v << speed_sum / num_particles << endl;
	}
	
	double temperature = 2.0 * total_kinetic_energy / (3.0 * kB * num_particles);
	
	cout << "avg speed: " << speed_sum / num_particles << " m/s" << endl;
	cout << "T = " << temperature << " K" << endl;
	// VERY APPROXIMATE temperature based on kinetic energy

	f << "Final positions:\n";
	for (i = 0; i < num_particles; i++) {
		f << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2];
		f << endl;
	}
	f.close();
	v.close();
	
	double t2 = double( clock() );
	cout << "Run time: " << (t2 - t0) / CLOCKS_PER_SEC << " seconds\n";

	return 0;
}


