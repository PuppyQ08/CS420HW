#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <cmath>
#include "particle.h"
#include "mdmpi.h"
#include <mpi.h>
#include <vector>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double t0 = MPI_Wtime();	// start the clock (real time)
	int i, j, k;
	const int num_particles = 2000;
	double cell_size[3] = { 64.0, 64.0, 64.0 };	// Position bounds in x, y, z
        double cutoff_radius = 8.0;

	const double mass = 4.0;		// Helium
	const double kB = 8314.0;		// amu * nm^2 ns^-2 K^-1 
	const double V = pow(cell_size[0], 3);		// cubic nanometers

	double a = 0.0;					// Van-der-Waals parameters
	double b = 0.0;
	
	double speed_sum;
	double total_kinetic_energy;

	double t_initial = 0.0;				// (simulation time)
	double t_final = 0.1;	
	double t_step = 0.001;
//___________________________________________________________________________________________
	int grid_dim = get_grid_dim(size);
	int grid_coordinates[3];
	get_grid_coord(rank, size, grid_coordinates); 	// Initialize grid_coordinates to the appropriate values

	bool has_edge[6] = { 0, 0, 0, 0, 0, 0 };
	for (i = 0; i < 3; i++) {
		if (grid_coordinates[2-i] == 0) { has_edge[2 * i] = true; }
		if (grid_coordinates[2-i] == grid_dim - 1) { has_edge[2 * i + 1] = true; }
	}


    //randomly generated coordinates of particles

	random_device rd;
	mt19937_64 rand(rd());	// 64 bit Mersenne Twister
	uniform_real_distribution<> generate_position(0.1, cell_size[0]-0.1);
	normal_distribution<> generate_velocity(0.0, 25.0);
	
	
        std::vector<double*> positions(num_particles);
        std::vector<double*> velocities(num_particles);
        std::vector<double*> accelerations(num_particles);
        //std::vector<bool> isexist(num_particles,1);
        
	//std::vector<double*> ghostshell;
	for (i = 0; i < num_particles; i++ ) {
		positions[i] = new double[3];
		velocities[i] = new double[3];
		accelerations[i] = new double[3];
	}

	for (i = 0; i < num_particles; i++) {
		for (j = 0; j < 3; j++) {
			positions[i][j] = generate_position(rand) + cell_size[j]*grid_coordinates[j];
			velocities[i][j] = generate_velocity(rand);
			if (velocities[i][j] > 1000.0) { velocities[i][j] = 10.0; }		// Make sure velocity isn't too high
			accelerations[i][j] = 0.0;
		}
	}
        //MPI Request of send and recv
        /*
	const int num_req = 3 * pow(grid_dim, 2) * (grid_dim - 1);
	MPI_Request* sreq = new MPI_Request[num_req];
	MPI_Request* rreq = new MPI_Request[num_req];
        */
	int adjacent_ranks[6];
        
        int summ = 0;
        for(i = 0;i<6;i++){if(!has_edge[i]){summ++;}}
        MPI_Request sreq[summ],rreq[summ];//TODO


	// figure out process ranks to send / receive from
	adjacent_ranks[0] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] - 1);
	adjacent_ranks[1] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] + 1);
	adjacent_ranks[2] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] - 1, grid_coordinates[2]);
	adjacent_ranks[3] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] + 1, grid_coordinates[2]);
	adjacent_ranks[4] = get_process_rank(grid_dim, grid_coordinates[0] - 1, grid_coordinates[1], grid_coordinates[2]);
	adjacent_ranks[5] = get_process_rank(grid_dim, grid_coordinates[0] + 1, grid_coordinates[1], grid_coordinates[2]);

	double t = t_initial;

	while (t < t_final) {

        /* main part needed parallel*/
        
        speed_sum = 0.0;
        total_kinetic_energy = 0.0;
        //recvbuf and sendbuf
        double**sendbuf = new double*[6];
        for(j = 0; j < 6;j++){
                sendbuf[j] = new double[num_particles];
        }
         
        int *sendidx = new int[6];
        for(i = 0;i<6;i++){sendidx[i] = 0;}

        int *recvidx = new int[6];
        double**recvbuf = new double*[6];
        for(j = 0; j < 6;j++){
                recvbuf[j] = new double[num_particles];
        }

                
        for (i = 0; i < num_particles; i++) {

                // update positions
                for (j = 0; j < 3; j++) {
                        positions[i][j] = positions[i][j] + velocities[i][j] * t_step + 0.5 * accelerations[i][j] * pow(t_step, 2);

                        // Elastic collisions
                        if (positions[i][j] < 0.0) {
                                positions[i][j] *= -1.0;
                                velocities[i][j] *= -1.0;
                        }
                        if (positions[i][j] > cell_size[j]*grid_dim) {
                                positions[i][j] = 2 * cell_size[j]*grid_dim - positions[i][j];
                                velocities[i][j] *= -1.0;
                        }
                        //lower face
                        if(positions[i][j]< grid_coordinates[j]*cell_size[j] + cutoff_radius){
                            if(!has_edge[2*j]){
                                  sendbuf[2*j][sendidx[2*j]*3] = positions[i][0]; 
                                  sendbuf[2*j][sendidx[2*j]*3+1] = positions[i][1]; 
                                  sendbuf[2*j][sendidx[2*j]*3+2] = positions[i][2]; 
                                  sendidx[2*j]++;
                            }
                        }
                        //higher face
                        if(positions[i][j]> (grid_coordinates[j]+1)*cell_size[j] - cutoff_radius){
                            if(!has_edge[2*j+1]){
                                  sendbuf[2*j+1][sendidx[2*j+1]*3] = positions[i][0]; 
                                  sendbuf[2*j+1][sendidx[2*j+1]*3+1] = positions[i][1]; 
                                  sendbuf[2*j+1][sendidx[2*j+1]*3+2] = positions[i][2]; 
                                  sendidx[2*j+1]++;
                            }
                        }

                }	
        } 
        //send and recv
        int sidx =0;
        int ridx =0;
        for(i =0; i <6; i++){
                if(!has_edge[i]){
                        MPI_Isend(&sendidx[i], 1, MPI_INT, adjacent_ranks[i], 0, MPI_COMM_WORLD, &sreq[sidx]);
                        MPI_Irecv(&recvidx[i], 1, MPI_INT, adjacent_ranks[i], 0, MPI_COMM_WORLD, &rreq[ridx]);
                        sidx++;
                        ridx++;
                }
        }
        MPI_Waitall(summ, sreq, MPI_STATUSES_IGNORE);
        MPI_Waitall(summ, rreq, MPI_STATUSES_IGNORE);
        sidx = 0;
        ridx = 0;
        for(i =0; i <6; i++){
                if(!has_edge[i]){
                        MPI_Isend(&sendbuf[i][0],sendidx[i], MPI_DOUBLE, adjacent_ranks[i], 1, MPI_COMM_WORLD, &sreq[2*sidx]);
                        MPI_Irecv(&recvbuf[i][0], recvidx[i], MPI_DOUBLE, adjacent_ranks[i], 1, MPI_COMM_WORLD,&rreq[2*ridx]);
                        sidx++;
                        ridx++;
                }
        }


        MPI_Waitall(summ, sreq, MPI_STATUSES_IGNORE);
        MPI_Waitall(summ, rreq, MPI_STATUSES_IGNORE);
        
        for(i = 0; i < num_particles;i++){
        double force[3] = { 0.0, 0.0, 0.0 };
        double prev_acceleration[3];
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
                //Calculate force from ghostshell
                for(j = 0; j < 6; j++){
                        if(!has_edge[j]){
                                for(int l = 0; l < recvidx[j];l++){
                                        double r = pow(recvbuf[j][3*l] - positions[i][0],2)+pow(recvbuf[j][3*l+1] - positions[i][1],2)+pow(recvbuf[j][3*l+2] - positions[i][2],2); 
                                        r = sqrt(r); 
                                        if(r < 8.0){
                                                double Fnorm = LJ_force(r, 0.4, 0.125);
                                                for (k = 0; k < 3; k++) { force[k] += Fnorm * (positions[i][k] - recvbuf[j][3*l+k]) / r; }
                                        }
                                }
                        }
                }

                // calculate acceleration and velocity of particle i
                for (j = 0; j < 3; j++) { 
                        accelerations[i][j] = force[j] / mass; 
                        velocities[i][j] = velocities[i][j] + 0.5 * (prev_acceleration[j] + accelerations[i][j]) * t_step;
                }

                if (speed(velocities, i) > 2000.0) {
                        // If the speed is too fast, scale down the speed while preserving the direction
                        double scale_factor = speed(velocities, i) / 100.0;
                        velocities[i][0] = velocities[i][0] / scale_factor;
                        velocities[i][1] = velocities[i][1] / scale_factor;
                        velocities[i][2] = velocities[i][2] / scale_factor;
                }

                speed_sum += speed(velocities, i);
                total_kinetic_energy += kinetic_energy(velocities, mass, i);
        }

		t += t_step;
        for(j = 0; j < 6;j++){
               delete sendbuf[j];
               delete recvbuf[j];
        }
                delete [] recvbuf;
                delete [] sendbuf;
                delete sendidx;
                delete recvidx;
                

	}
        double global_spdsum;	
        double global_temp;	
	double temperature = 2.0 * total_kinetic_energy / (3.0 * kB);
        MPI_Reduce(&speed_sum,&global_spdsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&temperature,&global_temp,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	//cout << "avg speed: " << speed_sum / num_particles << " m/s" << endl;
	//cout << "T = " << temperature << " K" << endl;
	// VERY APPROXIMATE temperature based on kinetic energy

	
	double t2 = MPI_Wtime();
        if(rank == 0){
	cout << "Avg Speed: " << global_spdsum/size << "m/s \n";
	cout << "T : " << global_temp << " K\n";
        }
	cout << "Run time: " << (t2 - t0) << " seconds\n";
	for (i = 0; i < num_particles; i++ ) {
		delete positions[i] ;
	      delete velocities[i] ;
		delete accelerations[i];
	}

	return 0;
}


