#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <ctime>
#include "particles.h"

using namespace std;

int main() {

	double t0 = double( clock() );
	int i, j;



	int t = 0;
	int t_final = 50;
	double t_step = 1.0;

    /* writes the coordinates to output file at every time step
     * now is off. you can set it on if you like
     */
	bool output_everything = false;	

	ofstream f;
	f.open("outputv2.txt");

    particles mdrun;
	while (t < t_final) {

        /* main part needed parallel*/

		for (i = 0; i < mdrun.getnum(); i++) {
            mdrun.calculate_force(i);
            mdrun.update_coordinates(i, t_step);
		}

		if (output_everything) {
            mdrun.write_file(f);
		}

		t += t_step;
	}
	
	f << "Final positions:\n";
    mdrun.write_file(f);
	f.close();
	
	double t1 = double( clock() );
	cout << "Time: " << (t1 - t0) / CLOCKS_PER_SEC << " seconds\n";

	return 0;
}



