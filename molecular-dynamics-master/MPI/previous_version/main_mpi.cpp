#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <ctime>
#include <cmath>
#include <mpi.h>
#include "mdmpi.h"
#include "particle.h"

int main(int argc, char* argv[]) {
	
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// Start the timer
	if (rank == 0) { double t0 = MPI_Wtime(); }

	/*
	 cell_size refers to the dimensions of an individual cell. Each cell is computed by one process.
	 The processes and their corresponding volumes are arranged in a 3D grid.
	 grid_dim is the number of cells along one dimension of the grid.
	 Each process has integer coordinates (x,y,z) within the larger grid. These are stored in grid_coordinates.
	 */
    /* Coordinates :
     * z      y
     * ^     ^
     * |    /
     * |   /
     * |  /
     * | /___
     * |/   / -> one cell size
     * |--------------------> x
     */

	const int num_particles = 64000;
	double cell_size[3] = { 128.0, 128.0, 128.0 };	// Position bounds in x, y, z
	double cutoff_radius = 8.0;
	
	const double mass = 4.0;		// Helium
	const double kB = 8314.0;		// amu * nm^2 ns^-2 K^-1 
	const double V = pow(cell_size[0], 3);		// cubic nanometers
	int i, j, k;

	int grid_dim = get_grid_dim(size);
	int grid_coordinates[3];
	get_grid_coord(rank, size, grid_coordinates); 	// Initialize grid_coordinates to the appropriate values

	bool has_edge[6] = { 0, 0, 0, 0, 0, 0 };
	for (i = 0; i < 3; i++) {
		if (grid_coordinates[2-i] == 0) { has_edge[2 * i] = true; }
		if (grid_coordinates[2-i] == grid_dim - 1) { has_edge[2 * i + 1] = true; }
	}
	
	double a = 0.0;					// Van-der-Waals parameters
	double b = 0.0;
	double speed_sum;
	double total_kinetic_energy;
	
	double t_initial = 0.0;				// (simulation time)
	double t_final = 0.1;
	double t_step = 0.001;
	
	// to do: get inputs from a file
	
	
	//randomly generated coordinates of particles
	std::random_device rd;
	std::mt19937_64 rand(rd());	// 64 bit Mersenne Twister
	std::uniform_real_distribution<> generate_position(0.1, cell_size[0]-0.1);
	std::normal_distribution<> generate_velocity(0.0, 1000.0); // normal distribution centered at 0, std deviation of 1000
	/*change array to vector
	double** positions = new double*[num_particles];
	double** velocities = new double*[num_particles];
	double** accelerations = new double*[num_particles];
        */
	
	std::vector<double*> positions(num_particles);
	std::vector<double*> velocities(num_particles);
	std::vector<double*> accelerations(num_particles);
        
	std::vector<double*> ghostshell;

	// Boundary regions of the cell to SEND to nearby cells
	particle* boundary_face_xy1 = new particle[5000];
	particle* boundary_face_xy2 = new particle[5000];
	particle* boundary_face_xz1 = new particle[5000];
	particle* boundary_face_xz2 = new particle[5000];
	particle* boundary_face_yz1 = new particle[5000];
	particle* boundary_face_yz2 = new particle[5000];

	particle* sendbufs[6] = { boundary_face_xy1, boundary_face_xy2, boundary_face_xz1, 
		boundary_face_xz2, boundary_face_yz1, boundary_face_yz2 };

	// arrays to RECEIVE position data from nearby cells.
	particle* boundary_ghost_face_xy1 = new particle[5000];
	particle* boundary_ghost_face_xy2 = new particle[5000];
	particle* boundary_ghost_face_xz1 = new particle[5000];
	particle* boundary_ghost_face_xz2 = new particle[5000];
	particle* boundary_ghost_face_yz1 = new particle[5000];
	particle* boundary_ghost_face_yz2 = new particle[5000];

	particle* recvbufs[6] = { boundary_ghost_face_xy1, boundary_ghost_face_xy2, boundary_ghost_face_xz1, 
		boundary_ghost_face_xz2, boundary_ghost_face_yz1, boundary_ghost_face_yz2 };
	
	// MPI version of particle class
	MPI_Datatype MPI_particle;
	MPI_Type_contiguous(9, MPI_DOUBLE, &MPI_particle);
	MPI_Type_commit(&MPI_particle);
	
	int boundary_index[6] = { 0, 0, 0, 0, 0, 0 };
	int boundary_index_recvbuf[6];
	int boundary_array_alloc[6] = {5000, 5000, 5000, 5000, 5000, 5000};
	
	for (i = 0; i < num_particles; i++ ) {
		positions[i] = new double[3];
		velocities[i] = new double[3];
		accelerations[i] = new double[3];
	}
	//initial atom positions in each grid (processes).
	for (i = 0; i < num_particles; i++) {
		for (j = 0; j < 3; j++) {
                        //TODO:check this part if it is correct
                        //for each cell we shift the postion to designated number.
			positions[i][j] = generate_position(rand) + cell_size[j]*grid_coordinates[j] ;
			velocities[i][j] = generate_velocity(rand);
			if (velocities[i][j] > 50000.0) { velocities[i][j] = 1000.0; }		// Make sure velocity isn't too high
			accelerations[i][j] = 0.0;
		}
		// determine if position is within cutoff radius of walls
		// then put particles whithin cutoffs to correspoding boundary array.
		//int grid_dim = [x,y,z];
		
		// copy the position and velocity to a new particle object
		particle p_copy(positions[i][0], positions[i][1], positions[i][2],velocities[i][0],velocities[i][1],velocities[i][2],accelerations[i][0],accelerations[i][1],accelerations[i][2]);
		
		for (j = 0; j< 3; j++) {
			//for each direction(xyz) we check position is if whithin two cutoff away from two boundary.
            //Also if the particles are outside of boundaries we put them in correspoding boundary array to send to neighbors
			//the lower face
            //This one includes the situation that particles are outside of boundaries
			if(positions[i][j]< grid_coordinates[j]*cell_size[j] + cutoff_radius){
				if(j==0){boundary_face_yz1[boundary_index[2*j]] = p_copy; boundary_index[2 * j]++;}
				else if(j==1){boundary_face_xz1[boundary_index[2 * j]] = p_copy; boundary_index[2 * j]++;}
				else if(j==2){boundary_face_xy1[boundary_index[2 * j]] = p_copy; boundary_index[2 * j]++;}
			}
			
            //higher face
            //This one includes the situation that particles are outside of boundaries
			else if(positions[i][j]> grid_coordinates[j]*cell_size[j] - cutoff_radius){
				if(j==0){boundary_face_yz2[boundary_index[2 * j+1]] = p_copy; boundary_index[2 * j + 1]++;}
				else if(j==1){boundary_face_xz2[boundary_index[2 * j+1]] = p_copy; boundary_index[2 * j+1]++;}
				else if(j==2){boundary_face_xy2[boundary_index[2 * j+1]] = p_copy; boundary_index[2 * j+1]++;}
			}
		
			if ( boundary_index[2*j] >= boundary_array_alloc[2 * j] ) {
				if (j == 0) { resize_array(boundary_face_xy1, boundary_array_alloc[0]); boundary_array_alloc[0] += 1000; }
				if (j == 1) { resize_array(boundary_face_xz1, boundary_array_alloc[2]); boundary_array_alloc[2] += 1000; }
				if (j == 2) { resize_array(boundary_face_yz1, boundary_array_alloc[4]); boundary_array_alloc[4] += 1000; }
			}
		
			if ( boundary_index[2*j + 1] >= boundary_array_alloc[2 * j + 1]) {
				if (j == 0) { resize_array(boundary_face_xy2, boundary_array_alloc[1]); boundary_array_alloc[1] += 1000; }
				if (j == 1) { resize_array(boundary_face_xz2, boundary_array_alloc[3]); boundary_array_alloc[3] += 1000; }
				if (j == 2) { resize_array(boundary_face_yz2, boundary_array_alloc[5]); boundary_array_alloc[5] += 1000; }
			}
		}
		
	}

	MPI_Request req[6];
	int adjacent_ranks[6];

	// figure out process ranks to send / receive from
	adjacent_ranks[0] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] - 1);
	adjacent_ranks[1] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] + 1);
	adjacent_ranks[2] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] - 1, grid_coordinates[2]);
	adjacent_ranks[3] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] + 1, grid_coordinates[2]);
	adjacent_ranks[4] = get_process_rank(grid_dim, grid_coordinates[0] - 1, grid_coordinates[1], grid_coordinates[2]);
	adjacent_ranks[5] = get_process_rank(grid_dim, grid_coordinates[0] + 1, grid_coordinates[1], grid_coordinates[2]);
	
	// Do the initial copying
	for (i = 0; i < 6; i++) {
		if (!has_edge[i]) {
			int recv_index = i - ( 2 * (i % 2) - 1);
			MPI_Isend(&sendbufs[i], boundary_index[i], MPI_particle, adjacent_ranks[i], 0, MPI_COMM_WORLD, &req[4 * i]);
			MPI_Irecv(&recvbufs[recv_index], 5000, MPI_particle, adjacent_ranks[i], 0, MPI_COMM_WORLD, &req[4 * i + 1]);

			MPI_Isend(&boundary_index[i], 1, MPI_INT, adjacent_ranks[0], 1, MPI_COMM_WORLD, &req[4 * i + 2]);
			MPI_Irecv(&boundary_index[recv_index], 1, MPI_INT, adjacent_ranks[0], 1, MPI_COMM_WORLD, &req[4 * i + 3]);
                        //for each edge we gather from the neighbor and append those particles that run into this cell to the bulk array.
			for(int k = 0; k < boundary_index[revc_index];k++){
             //For each edge, we check if particles in ghost shell already run into this cell
             	bool incell = true;
             	for(j = 0; j < 3; j++){
                    //if the particle crossed to this cell.
                    if(grid_coordinates[j]*cell_size[j] > recvbufs[k].coord[j] && recvbufs[k].coord[j] > grid_coordinates[j] * (cell_size[j]+1)){
                                        	       //TODO:finish rest positions is vector so it is easier to append and erase.
                         incell = false;
                         }
             		}
             		if(incell){
                      	double* new_particlecoord = new double[3];
                        double* new_particlevel = new double[3];
                        double* new_particleaccel = new double[3];
                        for(j = 0;j< 3;j++){
                                                new_particlecoord[j] = recvbufs[k].coord[j];
                                                new_particlevel[j] = recvbufs[k].vel[j];
                                                new_particleaccel[j] = recvbufs[k].acce[j];
                                        }		
                                        positions.push_back[new_particlecoord];        
                                }     
                                else{ 
                                        //If the particle is not crossing into different cell
                                        //We just push it into the ghostshell vector then later we will iterate it to calculate force.
                                        double* rest_particle = new double[3];
                                        for(j = 0; j<3; j++){rest_particle[j] = recvbufs[k].coord[j]; }
                                        ghostshell.push_back(rest_particle);
                                }
                                 
                               
                        }
		}
	}

	MPI_Waitall(6, req, MPI_STATUSES_IGNORE);
	
	if (rank == 0) {
		double t1 = MPI_Wtime();
		std::cout << "Generated initial conditions in " << (t1 - t0) / CLOCKS_PER_SEC << " seconds.\n";
	}
	
	
	/* writes the coordinates to output file at every time step
	 * now is off. you can set it on if you like
	 */
	bool output_everything = false;
	
	// output files 
	std::ofstream f, v;
	f.open("MD_output.txt");
	v.open("average_speed.txt");
	
	// TODO: format this more nicely
	f << "Initial positions:\n";
	for (i = 0; i < num_particles; i++) {
		f << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2];
		f << std::endl;
	}
	f << std::endl;
	
	double t = t_initial;
	
	while (t < t_final) {
		
		speed_sum = 0.0;
		total_kinetic_energy = 0.0;
		
		for ( i = 0; i < num_particles; i++) {
			
			// update positions
            // TODO: update for both bulk and 6 boundaries.
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
			}
		}

        //TODO: Add Recv to boundary and bulk particles.
        //then update particle i in bulk and boundary from j in Recv,boundary and bulk.
		for (i = 0; i < num_particles; i++) {
			double force[3] = { 0.0, 0.0, 0.0 };
			double prev_acceleration[3];
			
			// Calculate forces from particles within the cell
			for (j = 0; j < num_particles; j++) {

				if (i != j) {
					double r = distance(positions, i, j);
					
					if (r < cutoff_radius) {
						
						double Fnorm = LJ_force(r, 0.4, 0.125);
						
						// multiply the norm of the force by each component of the unit vector
						for (k = 0; k < 3; k++) { force[k] += Fnorm * (positions[i][k] - positions[j][k]) / r; }
					}
				}
			}
			
			// TODO: calculate forces from boundaries outside the cell
			for (j = 0; j < 6; j++) {
				if (!has_edge[j] && ( positions[i][j / 3] < cutoff_radius ) ) {
					for (k = 0; k < boundary_index_recvbuf[j]; k++) {
						double cell_position[3];

						// convert the coordinates of the external cell to the internal coordinate frame
						cell_position[0] = recvbufs[j][k].x + ( ((j / 3) == 0) * 64 ) * ( 2 * (j % 2) - 1 ) ;
						cell_position[1] = recvbufs[j][k].y + (((j / 3) == 1) * 64) * (2 * (j % 2) - 1);
						cell_position[2] = recvbufs[j][k].z + (((j / 3) == 2) * 64) * (2 * (j % 2) - 1);

						double r = external_distance(positions[i][0], positions[i][1], positions[i][2],
							cell_position[0], cell_position[1], cell_position[2]);
						if (r < cutoff_radius) {
							double Fnorm = LJ_force(r, 0.4, 0.125);

							// multiply the norm of the force by each component of the unit vector
							for (k = 0; k < 3; k++) { force[k] += Fnorm * (positions[i][k] - positions[j][k]) / r; }
						}
					}
				}
			}

			// calculate acceleration and velocity of particle i
			for (j = 0; j < 3; j++) {
				prev_acceleration[j] = accelerations[i][j];
				accelerations[i][j] = force[j] / mass;
				velocities[i][j] = velocities[i][j] + 0.5 * (prev_acceleration[j] + accelerations[i][j]) * t_step;
			}
			
			if (speed(velocities, i) > 50000.0) {
				// If the speed is too fast, scale down the speed while preserving the direction
				double scale_factor = speed(velocities, i) / 10000.0;
				velocities[i][0] = velocities[i][0] / scale_factor;
				velocities[i][1] = velocities[i][1] / scale_factor;
				velocities[i][2] = velocities[i][2] / scale_factor;
			}
			
			speed_sum += speed(velocities, i);
			total_kinetic_energy += kinetic_energy(velocities, mass, i);
		}
		
		if (output_everything) {
			f << "Step " << t << std::endl;
			for (j = 0; j < num_particles; j++) {
				f << positions[j][0] << ", " << positions[j][1] << ", " << positions[j][2];
				f << std::endl;
			}
			f << std::endl;
		}
		t += t_step;
		
		v << speed_sum / num_particles << std::endl;
	}
	
	MPI_Type_free(&MPI_particle);
	
	double temperature = 2.0 * total_kinetic_energy / (3.0 * kB * num_particles);
	
	std::cout << "avg speed: " << speed_sum / num_particles << " m/s" << std::endl;
	std::cout << "T = " << temperature << " K" << std::endl;
	// VERY APPROXIMATE temperature based on kinetic energy
	
	f << "Final positions:\n";
	for (i = 0; i < num_particles; i++) {
		f << positions[i][0] << ", " << positions[i][1] << ", " << positions[i][2];
		f << std::endl;
	}
	f.close();
	v.close();
	
	if (rank == 0) {
		double t2 = MPI_Wtime();
		cout << "Run time: " << (t2 - t0) / CLOCKS_PER_SEC << " seconds\n";
	}
	
	
	MPI_Finalize();
	
	return 0;
}
