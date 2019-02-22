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
    double t0 =  MPI_Wtime() ; 

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

	int num_particles = 2000;
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
	double t_step = 0.1;
	
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
        std::vector<bool> isexist(num_particles,1);
        
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
	MPI_Type_contiguous(6, MPI_DOUBLE, &MPI_particle);
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
                
		
	}


	const int num_req = 3 * pow(grid_dim, 2) * (grid_dim - 1);
    MPI_Request* req = new MPI_Request[num_req];
	int adjacent_ranks[6];

	// figure out process ranks to send / receive from
	adjacent_ranks[0] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] - 1);
	adjacent_ranks[1] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] + 1);
	adjacent_ranks[2] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] - 1, grid_coordinates[2]);
	adjacent_ranks[3] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] + 1, grid_coordinates[2]);
	adjacent_ranks[4] = get_process_rank(grid_dim, grid_coordinates[0] - 1, grid_coordinates[1], grid_coordinates[2]);
	adjacent_ranks[5] = get_process_rank(grid_dim, grid_coordinates[0] + 1, grid_coordinates[1], grid_coordinates[2]);
	
	double t = t_initial;
	
	while (t < t_final) {
		
		speed_sum = 0.0;
		total_kinetic_energy = 0.0;
	        for(i = 0 ; i < 6;i++){boundary_index[i] = 0;}
                
                //-----------------------------------------------------------------------------------------------------------------------------------
		for ( i = 0; i < num_particles; i++) {
                            
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
			}



                        // determine if position is within cutoff radius of walls
                        // then put particles whithin cutoffs to correspoding boundary array.
                        
                        // copy the position and velocity to a new particle object
                        particle p_copy(positions[i][0], positions[i][1], positions[i][2],velocities[i][0],velocities[i][1],velocities[i][2],accelerations[i][0],accelerations[i][1],accelerations[i][2]);
                        
                        for (j = 0; j< 3; j++) {
                                //for each direction(xyz) we check position is if whithin two cutoff away from two boundary.
                                //Also if the particles are outside of boundaries we put them in correspoding boundary array to send to neighbors
                                //the lower face
                                if(positions[i][j]< grid_coordinates[j]*cell_size[j] + cutoff_radius){
                                        if(j==0){boundary_face_xy1[boundary_index[2*j]] = p_copy; boundary_index[2 * j]++;}
                                        else if(j==1){boundary_face_xz1[boundary_index[2 * j]] = p_copy; boundary_index[2 * j]++;}
                                        else if(j==2){boundary_face_yz1[boundary_index[2 * j]] = p_copy; boundary_index[2 * j]++;}
                                        //if the particle is outside of cell then we mark it as non exist
                                        if(positions[i][j] < grid_coordinates[j]*cell_size[j]){
                                        	isexist[i] = false;
                                        }
                                        	
                                }
                                
                                
                                //higher face
                                else if(positions[i][j]> (grid_coordinates[j]+1)*cell_size[j] - cutoff_radius){
                                        if(j==0){boundary_face_xy2[boundary_index[2 * j+1]] = p_copy; boundary_index[2 * j + 1]++;}
                                        else if(j==1){boundary_face_xz2[boundary_index[2 * j+1]] = p_copy; boundary_index[2 * j+1]++;}
                                        else if(j==2){boundary_face_yz2[boundary_index[2 * j+1]] = p_copy; boundary_index[2 * j+1]++;}
                                        if(positions[i][j] > (grid_coordinates[j] +1)*cell_size[j]){
                                        	isexist[i] = false;
                                        }
                                }
                                //TODO:somehow it doesnot resize well. Cause segfault. 
                                if ( boundary_index[2*j] >= boundary_array_alloc[2*j]-1 ) {
                                        if (j == 0) { resize_array(boundary_face_xy1, boundary_array_alloc[0]); boundary_array_alloc[0] += 1000; }
                                        if (j == 1) { resize_array(boundary_face_xz1, boundary_array_alloc[2]); boundary_array_alloc[2] += 1000; }
                                        if (j == 2) { resize_array(boundary_face_yz1, boundary_array_alloc[4]); boundary_array_alloc[4] += 1000; }
                                }
                        
                                if ( boundary_index[2*j + 1] >= boundary_array_alloc[2*j+1]-1 ) {
                                        if (j == 0) { resize_array(boundary_face_xy2, boundary_array_alloc[1]); boundary_array_alloc[1] += 1000; }
                                        if (j == 1) { resize_array(boundary_face_xz2, boundary_array_alloc[3]); boundary_array_alloc[3] += 1000; }
                                        if (j == 2) { resize_array(boundary_face_yz2, boundary_array_alloc[5]); boundary_array_alloc[5] += 1000; }
                                }
                                
                        }
		}
                if(rank==0){std::cout<<"first one"<<boundary_index[1]<<std::endl;}
/*
                for(i = 0; i < num_particles;i++){
                	if(isexist[i]==0){ 
                                positions.erase(positions.begin()+i);
                                velocities.erase(velocities.begin()+i);
                                accelerations.erase(accelerations.begin()+i);
                                isexist.erase(isexist.begin()+i);
                        }
                }
                */
                //Copying from 6 directions
                //-----------------------------------------------------------------------------------------------------------------------------------

                MPI_Barrier(MPI_COMM_WORLD);
                int idx = 0;
                for (i = 0; i < 6; i++) {
                        if (!has_edge[i]) {
                                /* 
                                MPI_Isend(&sendbufs[i], boundary_index[i], MPI_particle, adjacent_ranks[i], 0, MPI_COMM_WORLD, &req[4 * idx]);
                                MPI_Irecv(&recvbufs[i], 5000, MPI_particle, adjacent_ranks[i], 0, MPI_COMM_WORLD, &req[4 * idx + 1]);
*/
                                MPI_Isend(&boundary_index[i], 1, MPI_INT, adjacent_ranks[i], 1, MPI_COMM_WORLD, &req[2 * idx ]);
                                MPI_Irecv(&boundary_index_recvbuf[i], 1, MPI_INT, adjacent_ranks[i], 1, MPI_COMM_WORLD, &req[2 * idx + 1]);
                                idx++;
                                //for each edge we gather from the neighbor and append those particles that run into this cell to the bulk array.
                                //std::cout<<boundary_index_recvbuf[i]<<std::endl;
/*
                                for(int k = 0; k < boundary_index_recvbuf[i];k++){
                                        //For each edge, we check if particles in ghost shell already run into this cell
                                        bool incell = true;
                                        for(j = 0; j < 3; j++){
                                                //if the particle crossed to this cell.
                                                if(grid_coordinates[j]*cell_size[j] > recvbufs[i][k].coord[j] && recvbufs[i][k].coord[j] > grid_coordinates[j] * (cell_size[j]+1)){
                                                                 incell = false;
                                                }
                                        }                 
                                        if(incell){
                                                positions.push_back(new double[3]);        
                                                velocities.push_back(new double[3]);        
                                                accelerations.push_back(new double[3]);        
                                                for(j = 0;j< 3;j++){
                                                        positions.back()[j] = recvbufs[i][k].coord[j];
                                                        velocities.back()[j] = recvbufs[i][k].vel[j];
                                                        accelerations.back()[j] = recvbufs[i][k].acce[j];
                                                        isexist.push_back(true);
                                                        num_particles++;
                                                }		
                                        }     
                                        else{ 
                                                //If the particle is not crossing into different cell
                                                //We just push it into the ghostshell vector then later we will iterate it to calculate force.
                                                ghostshell.push_back(new double[3]);
                                                for(j = 0; j<3; j++){ghostshell.back()[j]= recvbufs[i][k].coord[j]; }
                                        }
                                       
                                       
                                }
                                */

                        }
                }
                //std::cout<<positions.size()<<"hdjsfkh"<<std::endl;
                //std::cout<<num_particles<<"qqy"<<std::endl;
                MPI_Waitall(num_req, req, MPI_STATUSES_IGNORE);
                std::cout<<rank<<std::endl;
                for(i = 0; i<6; i++){
                        if(!has_edge[i]){
                        if(rank==0){
                        std::cout<<"secondone"<<i<<boundary_index_recvbuf[i]<<std::endl;}
                        }
                        }
                        //if(rank==0){ std::cout<<i<<"hhhh"<<std::endl;}
                //-----------------------------------------------------------------------------------------------------------------------------------
                //then update particle i forces and accelarations in bulk and boundary from j in Recv,boundary and bulk.
		for (i = 0; i < num_particles; i++) {
                if(isexist[i] !=0){
			double force[3] = { 0.0, 0.0, 0.0 };
			double prev_acceleration[3];
                       // if(rank==0){ std::cout<<i<<"hhhh"<<std::endl;}
			
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
                        //Calculate forces from ghost shell
                        for(j = 0; j < ghostshell.size();j++){
					double r = distance(positions, i, j);
					
					if (r < cutoff_radius) {
						
						double Fnorm = LJ_force(r, 0.4, 0.125);
						
						// multiply the norm of the force by each component of the unit vector
						for (k = 0; k < 3; k++) { force[k] += Fnorm * (positions[i][k] - positions[j][k]) / r; }
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
            }
                
	        t += t_step;
                //-----------------------------------------------------------------------------------------------------------------------------------
                std::cout<<positions.size()<<"dsfjhhhhhhh"<<std::endl;
		
	}
                std::cout<<positions.size()<<"hhhhhhh"<<std::endl;
        //MPI_Gather(&speed_sum,1,MPI_DOUBLE,&speed_sum,1,MPI_DOUBLE,0,MPI_COMM_WORLD);	
        //MPI_Gather(&total_kinetic_energy,1,MPI_DOUBLE,&total_kinetic_energy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);	
	MPI_Type_free(&MPI_particle);
	delete[] req;
	
	double temperature = 2.0 * total_kinetic_energy / (3.0 * kB * num_particles);
	double t1 = MPI_Wtime();

	if (rank == 0) {		
		std::cout << "total time is " << (t1 - t0) / CLOCKS_PER_SEC << " seconds.\n";
	}
	
	MPI_Finalize();
    return 0;
    }
