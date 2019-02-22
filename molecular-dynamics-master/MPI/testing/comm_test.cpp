#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <mpi.h>

#include "mdmpi.h"
#include "particle.h"

int main(int argc, char* argv[]) {
	
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	// Start the timer
	if (rank == 0) { double t0 = double( clock() ); }

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

	const int num_particles = 10000;
	double cell_size[3] = { 100.0, 100.0, 100.0 };	// Position bounds in x, y, z
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
	double t_final = 0.002;
	double t_step = 0.001;
	
	// to do: get inputs from a file
	
	

	//randomly generated coordinates of particles
	std::random_device rd;
	std::mt19937_64 rand(rd());	// 64 bit Mersenne Twister
	std::uniform_real_distribution<> generate_position(0.1, cell_size[0]-0.1);
	std::normal_distribution<> generate_velocity(0.0, 1000.0); // normal distribution centered at 0, std deviation of 1000
	
	double** positions = new double*[num_particles];
	double** velocities = new double*[num_particles];
	double** accelerations = new double*[num_particles];
	
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

	std::vector<std::string> face_names;
	face_names.push_back("lower xy");
	face_names.push_back("upper xy");
	face_names.push_back("lower xz");
	face_names.push_back("upper xz");
	face_names.push_back("lower yz");
	face_names.push_back("upper yz");
	
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
			positions[i][j] = generate_position(rand);
			velocities[i][j] = generate_velocity(rand);
			if (velocities[i][j] > 50000.0) { velocities[i][j] = 1000.0; }		// Make sure velocity isn't too high
			accelerations[i][j] = 0.0;
		}
		// determine if position is within cutoff radius of walls
		// then put particles whithin cutoffs to correspoding boundary array.
		//int grid_dim = [x,y,z];
		
		// copy the positions to a new particle object
		particle p_copy(positions[i][0], positions[i][1], positions[i][2],velocities[i][0],velocities[i][1],velocities[i][2]);
		
		for (j = 0; j< 3; j++) {
			//for each direction(xyz) we check position is if whithin two cutoff away from two boundary.
            //Also if the particles are outside of boundaries we put them in correspoding boundary array to send to neighbors
			//the lower face
            //This one includes the situation that particles are outside of boundaries
			if ( positions[i][j] < cutoff_radius ) {
				sendbufs[2 * j][boundary_index[2 * j]] = p_copy;
				boundary_index[2 * j]++;
			}
			
            //higher face
            //This one includes the situation that particles are outside of boundaries
			else if (positions[i][j] > (cell_size[j] - cutoff_radius) ){
				sendbufs[2 * j + 1][boundary_index[2 * j + 1]] = p_copy;
				boundary_index[2 * j + 1]++;
			}
		
			if ( boundary_index[2*j] >= boundary_array_alloc[2*j] - 1 ) {
				if (j == 0) { resize_array(boundary_face_xy1, boundary_array_alloc[0]); boundary_array_alloc[0] += 1000; }
				if (j == 1) { resize_array(boundary_face_xz1, boundary_array_alloc[2]); boundary_array_alloc[2] += 1000; }
				if (j == 2) { resize_array(boundary_face_yz1, boundary_array_alloc[4]); boundary_array_alloc[4] += 1000; }
			}
		
			if ( boundary_index[2*j+1] >= boundary_array_alloc[2*j+1] - 1) {
				if (j == 0) { resize_array(boundary_face_xy2, boundary_array_alloc[1]); boundary_array_alloc[1] += 1000; }
				if (j == 1) { resize_array(boundary_face_xz2, boundary_array_alloc[3]); boundary_array_alloc[3] += 1000; }
				if (j == 2) { resize_array(boundary_face_yz2, boundary_array_alloc[5]); boundary_array_alloc[5] += 1000; }
			}
		}
		
	}

	//MPI_Request req[12];
	int adjacent_ranks[6];

	// figure out process ranks to send / receive from
	// doesn't check bounds. i.e. 0 <= rank < size
	adjacent_ranks[0] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] - 1);
	adjacent_ranks[1] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] + 1);
	adjacent_ranks[2] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] - 1, grid_coordinates[2]);
	adjacent_ranks[3] = get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1] + 1, grid_coordinates[2]);
	adjacent_ranks[4] = get_process_rank(grid_dim, grid_coordinates[0] - 1, grid_coordinates[1], grid_coordinates[2]);
	adjacent_ranks[5] = get_process_rank(grid_dim, grid_coordinates[0] + 1, grid_coordinates[1], grid_coordinates[2]);

	/*
	std::cout << "I am process number " << rank << std::endl;
	std::cout << "The process number above me in z direction is ";
	std::cout << get_process_rank(grid_dim, grid_coordinates[0], grid_coordinates[1], grid_coordinates[2] + 1);
	std::cout << "\n-----\n";
	*/
	
	for (i = 0; i < 6; i++) {
		//printf("(%d, %d, %d, %d)\t", rank, adjacent_ranks[i], boundary_index[i], boundary_index_recvbuf[i]);
	}
	printf("\n");
	MPI_Barrier;
	MPI_Request req_index[6];
	MPI_Request req_data[6];
	
	// test code for size = 8
	if (rank == 0) {
		MPI_Isend(&boundary_index[1], 1, MPI_INT, adjacent_ranks[1], 0, MPI_COMM_WORLD, &req_index[0]);
		MPI_Isend(&boundary_index[3], 1, MPI_INT, adjacent_ranks[3], 0, MPI_COMM_WORLD, &req_index[1]);
		MPI_Isend(&boundary_index[5], 1, MPI_INT, adjacent_ranks[5], 0, MPI_COMM_WORLD, &req_index[2]);
		MPI_Irecv(&boundary_index_recvbuf[1], 1, MPI_INT, adjacent_ranks[1], 0, MPI_COMM_WORLD, &req_index[3]);
		MPI_Irecv(&boundary_index_recvbuf[3], 1, MPI_INT, adjacent_ranks[3], 0, MPI_COMM_WORLD, &req_index[4]);
		MPI_Irecv(&boundary_index_recvbuf[5], 1, MPI_INT, adjacent_ranks[5], 0, MPI_COMM_WORLD, &req_index[5]);
		MPI_Waitall(6, req_index, MPI_STATUSES_IGNORE);
	}
	
	if (rank == 1) {
		MPI_Isend(&boundary_index[0], 1, MPI_INT, adjacent_ranks[0], 0, MPI_COMM_WORLD, &req_index[3]);
		MPI_Isend(&boundary_index[3], 1, MPI_INT, adjacent_ranks[3], 0, MPI_COMM_WORLD, &req_index[1]);
		MPI_Irecv(&boundary_index_recvbuf[0], 1, MPI_INT, adjacent_ranks[0], 0, MPI_COMM_WORLD, &req_index[0]);
		MPI_Irecv(&boundary_index_recvbuf[3], 1, MPI_INT, adjacent_ranks[3], 0, MPI_COMM_WORLD, &req_index[4]);
	}
	
	if (rank == 2) {
		MPI_Isend(&boundary_index[2], 1, MPI_INT, adjacent_ranks[2], 0, MPI_COMM_WORLD, &req_index[4]);
		MPI_Irecv(&boundary_index_recvbuf[2], 1, MPI_INT, adjacent_ranks[2], 0, MPI_COMM_WORLD, &req_index[1]);
	}
	
	if (rank == 3) {
		MPI_Isend(&boundary_index[2], 1, MPI_INT, adjacent_ranks[2], 0, MPI_COMM_WORLD, &req_index[4]);
		MPI_Irecv(&boundary_index_recvbuf[2], 1, MPI_INT, adjacent_ranks[2], 0, MPI_COMM_WORLD, &req_index[1]);
	}
	
	if (rank == 4) {
		MPI_Isend(&boundary_index[4], 1, MPI_INT, adjacent_ranks[4], 0, MPI_COMM_WORLD, &req_index[5]);
		MPI_Irecv(&boundary_index_recvbuf[4], 1, MPI_INT, adjacent_ranks[4], 0, MPI_COMM_WORLD, &req_index[2]);
	}
	
	
	if (rank == 0) {
		std::cout << "process 0 has boundary_index_recvbuf of " << boundary_index_recvbuf[1] << std::endl;
		//MPI_Isend(&sendbufs[1], boundary_index_recvbuf[1], MPI_particle, adjacent_ranks[1], 1, MPI_COMM_WORLD, &req_data[0]);
		//MPI_Irecv(&recvbufs[1], 5000, MPI_particle, adjacent_ranks[1], 1, MPI_COMM_WORLD, &req_data[1]);
		//MPI_Waitall(2, req_data, MPI_STATUSES_IGNORE);
	}
	
	if (rank == 1) {
		
		//MPI_Isend(&sendbufs[0], boundary_index_recvbuf[0], MPI_particle, adjacent_ranks[0], 1, MPI_COMM_WORLD, &req_data[1]);
		//MPI_Irecv(&recvbufs[0], 5000, MPI_particle, adjacent_ranks[0], 1, MPI_COMM_WORLD, &req_data[0]);
	}
	
	
	/*
	for (i = 0; i < 6; i++) {
		int recv_index = i - (2 * (i % 2) - 1);
		if (!has_edge[i]) {
			MPI_Isend(&boundary_index[i], 1, MPI_INT, adjacent_ranks[0], 0, MPI_COMM_WORLD, &req[4 * i]);
			MPI_Irecv(&boundary_index[recv_index], 1, MPI_INT, adjacent_ranks[0], 0, MPI_COMM_WORLD, &req[4 * i + 1]);

			//MPI_Isend(&sendbufs[i], boundary_index[i], MPI_particle, adjacent_ranks[i], 0, MPI_COMM_WORLD, &req[4 * i + 2]);
			//MPI_Irecv(&recvbufs[recv_index], 5000, MPI_particle, adjacent_ranks[i], 0, MPI_COMM_WORLD, &req[4 * i + 3]);
		}
	}
	
	MPI_Waitall(24, req, MPI_STATUSES_IGNORE);
	*/

	if (rank == 0) {
		double t1 = double(clock());
		//std::cout << "Generated initial conditions in " << (t1 - t0) / CLOCKS_PER_SEC << " seconds.\n";
		// This doesn't recognize t0 as a variable. I'm not sure why --Aidan
	}

	/*
	std::string fname = "proc_" + std::to_string(rank);
	fname += ".out";
	std::ofstream test;
	test.open(fname);

	test << "Boundary layer data for this cell: ";
	test << '(' << grid_coordinates[0] << ", " << grid_coordinates[1] << ", " << grid_coordinates[2] << ")\n\n";
	for (i = 0; i < 6; i++) {
		test << face_names[i] << " face\n";
		for (j = 0; j < boundary_index[i]; j++) { test << sendbufs[i][j]; }
	}

	test << "\n\n\n";
	test << "Boundary layer data from the adjacent cells:\n";
	for (i = 0; i < 6; i++) {
		if (has_edge[i]) {
			test << face_names[i] << "face\n";
			for (j = 0; j < 1000; j++) { test << recvbufs[i][j]; }
		}
	}
	
	test.close();
	*/
	
	if (rank == 0) {
		double t2 = double(clock());
		//std::cout << "Run time: " << (t2 - t0) / CLOCKS_PER_SEC << " seconds\n";
	}

	
	MPI_Finalize();
	
	return 0;
}
