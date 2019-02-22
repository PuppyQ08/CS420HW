#pragma once
#include <cmath>
#include <iostream>
using namespace std;
class particles {

public:

// Constructors, copy, destructor
	particles() { this->create(); }
	particles(const particles& p) { copy(p); }
	~particles() { this->destroy(); }

// Methods

	particles& operator= (particles& p);
	void reset_force(int i);
	void update_coordinates(int i , double dt);
	void calculate_force(int i);
	double getnum() { return num_particles; }
	double get_distance(int i, int j );
    std::ostream& write_file(std::ostream& out);
private:

	void create();
	void copy(const particles& p);
	void destroy();

// Representation:
    double bounds;
    int num_particles;
	double *mass;
	double *charge;

	double** x;
	double** x_prev;
	double** F;

};
/*write to file funtion to write calculation output to files */
std::ostream& particles::write_file(std::ostream& out) {
    for (int i = 0; i < num_particles; i++) {
        out << x[0][i] <<"    " << x[1][i] <<"    " << x[2][i] << endl;
    }
	return out;
}
/* default constructor*/
void particles::create() {
    //read the coordinates data
    ifstream iptfile("initial_coordinates.txt");
    iptfile >> num_particles >> bounds;

    //initial the memory to store coordinates and mass charge
	x = new double*[3];
	x_prev = new double*[3];
	F = new double*[3];
    mass = new double[num_particles];
    charge = new double[num_particles];
    for (int j = 0; j < 3; j++){
        x[j] =  new double[num_particles];
        x_prev[j] =  new double[num_particles];
        F[j] =  new double[num_particles];
    }
    for (int i = 0; i < num_particles; i++) {
        iptfile >> x[0][i] >> x[1][i] >> x[2][i] >> mass[i] >> charge[i] ;
    }
}




/*copy constructor*/
void particles::copy(const particles& p) {
	this->mass = new double [p.num_particles]; 
	this->charge = new double [p.num_particles]; 
	this->x = new double*[3];
	this->x_prev = new double*[3];
	this->F = new double*[3];
	for (int i = 0; i < 3; i++) {
        x[i] =  new double[p.num_particles];
        x_prev[i] =  new double[p.num_particles];
        F[i] =  new double[p.num_particles];
        for (int j = 0; j < p.num_particles; j++) {
		this->x[i][j] = p.x[i][j];
		this->x_prev[i][j] = p.x_prev[i][j];
		this->F[i][j] = p.F[i][j];
        }
	}
}


/*destroy helper function for destructor*/
void particles::destroy() {
    for (int i = 0; i < 3; i++) {
        delete x[i];
        delete x_prev[i];
        delete F[i];
    }

	delete[] x;
	delete[] x_prev;
	delete[] F;
	delete mass;
	delete charge;
}

/*overload particle equal operator */
particles& particles::operator= (particles& p) {
    destroy();
	this->copy(p);
	return *this;
}


/*For each particle, pass in time peroid and use calculated force to update the coordinates*/
void particles::update_coordinates(int j, double dt) {

	double temp[3];

	for (int i = 0; i < 3; i++) {
		temp[i] = 2 * x[i][j] - x_prev[i][j] + F[i][j] / mass[j] * dt * dt;
		if (temp[i] < 0) { temp[i] = -1 * temp[i]; }
        /*
         * If the new calculated coordinates are outside the 100*100 square
         * then bounce back
         */
		if (temp[i] > bounds) {
			double displacement = abs(bounds - temp[i]);
			temp[i] = bounds - displacement;
		}
		x_prev[i][j]= x[i][j];
		x[i][j] = temp[i];
	}

}

/*First get distance between this particle and passed in particles.
 * If distance is closed enough do calculation
 */
void particles::calculate_force(int i) {
    /*For particle i, we calculate force at k dimension between i and j,
     */
        reset_force(i);
        for (int j = 0; j < num_particles; j++) {
            if( i != j){
	            double distance = this->get_distance(i,j);

	            if (distance < 50.0) {
                    double direction[3];
                    double magnitude = charge[i] * charge[j] / pow(distance, 2);
                    double norm = 0.0;

                    for (int k = 0; k < 3; k++) {
                        direction[k] = x[k][i] - x[k][j];
                        norm += pow(direction[k], 2);
                    }

	}

        }
        }
}

/*To reset the force value*/
void particles::reset_force(int i){
    for (int k = 0; k < 3; k++) {
        F[k][i] = 0.0;
    }
}

/*To get distance between particle i and j*/
inline double particles::get_distance(int i, int j) {
	double output = 0.0;
	for (int k = 0; k < 3; k++) {
		output += pow(x[k][i] - x[k][j], 2);
	}
	return sqrt(output);
}
