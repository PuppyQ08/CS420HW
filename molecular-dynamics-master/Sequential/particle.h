#pragma once
#include <cmath>

class particle {

public:

// Constructors, copy, destructor
	particle() { this->create(); }		
	particle(double* x0, double m, double q, int ID) { this->create(x0, m, q, ID); }
	particle(const particle& p) { copy(p); }
	~particle() { this->destroy(); }

// Methods

	particle& operator= (particle& p);
	void print() { printf("particle %d: (%f, %f, %f)\n", ID_, x[0], x[1], x[2]); }
	void set_coordinates(double* coord);
	void set_prev_coordinates(double* coord);
	void reset_force() { F[0] = 0.0; F[1] = 0.0; F[2] = 0.0; }

	void update_coordinates(double* bounds, double dt);
	void calculate_force(particle& p);

	double* get_coordinates() { return x; }
	double* get_prev_coordinates() { return x_prev; }
	double* get_force() { return F; }

	double get_mass() { return mass; }
	double get_charge() { return charge; }
	int get_ID() { return ID_; }
	double get_distance(particle& p1);

private:

	void create();
	void create(double* x0, double m, double q, int ID);
	void copy(const particle& p);
	void destroy();

// Representation:

	int ID_;
	double mass;
	double charge;

	double* x;
	double* x_prev;
	double* F;

};

/* 
 * overload operator << to print out :
 * Particle "ID_number" : (position x y z)
 * pass into the output file
 */ 

std::ostream& operator<< (std::ostream& out, particle& p) {
	out << "Particle " << p.get_ID() << ": (";
	double* position = p.get_coordinates();
	out << position[0] << ", " << position[1] << ", " << position[2] << ")\n";
	return out;
}

/* default constructor*/
void particle::create() {
	mass = 1.0;
	charge = 0.0;
	x = new double[3];
	x_prev = new double[3];
	F = new double[3];

	for (int i = 0; i < 3; i++) {
		x[i] = 0.0;
		x_prev[i] = 0.0;
		F[i] = 0.0;
	}
}

/*constructor for initial particles position 
 * (randomly generated)
 */
void particle::create(double* x0, double m, double q, int ID) {
	mass = m;
	charge = q;
	ID_ = ID;
	x = new double[3];
	x_prev = new double[3];
	F = new double[3];
	for (int i = 0; i < 3; i++) {
		x[i] = x0[i];
		x_prev[i] = x0[i];
		F[i] = 0.0;
	}
}

/*copy constructor*/
void particle::copy(const particle& p) {
	this->mass = p.mass;
	this->charge = p.charge;
	this->ID_ = p.ID_;
	this->x = new double[3];
	this->x_prev = new double[3];
	this->F = new double[3];
	for (int i = 0; i < 3; i++) {
		this->x[i] = p.x[i];
		this->x_prev[i] = p.x_prev[i];
		this->F[i] = p.F[i];
	}
}

/*destroy helper function for destructor*/
void particle::destroy() {
	delete[] x;
	delete[] x_prev;
	delete[] F;
	mass = 0;
	charge = 0;
	ID_ = 0;
}

/*overload particle equal operator*/
particle& particle::operator= (particle& p) {
	delete[] x;
	delete[] x_prev;
	delete[] F;
	this->copy(p);
	return *this;
}

/*Pass in new coordinates from calculation*/
inline void particle::set_coordinates(double * new_coord) {
	for (int i = 0; i < 3; i++) { x[i] = new_coord[i]; }
}

/*Update previous coordinates for next calculation*/
inline void particle::set_prev_coordinates(double * new_coord) {
	for (int i = 0; i < 3; i++) { x_prev[i] = new_coord[i]; }
}


/*For each particle, pass in time peroid and use calculated force to update the coordinates*/
void particle::update_coordinates(double* bounds, double dt) {

	double temp[3];

	//for (int i = 0; i < 3; i++) { std::cout << F[i] << '\t'; }
	//std::cout << std::endl;

	for (int i = 0; i < 3; i++) {
		temp[i] = 2 * x[i] - x_prev[i] + F[i] / mass * dt * dt;
		if (temp[i] < 0) { temp[i] = -1 * temp[i]; }
        /*
         * If the new calculated coordinates are outside the 100*100 square 
         * then bounce back
         */
		if (temp[i] > bounds[i]) {
			double displacement = abs(bounds[i] - temp[i]);
			temp[i] = bounds[i] - displacement;
		}
		x_prev[i] = x[i];
		x[i] = temp[i];
	}

}

/*First get distance between this particle and passed in particles.
 * If distance is closed enough do calculation
 */
void particle::calculate_force(particle& p) {
	double distance = this->get_distance(p);
	
	if (distance < 50.0) {
		double direction[3];
		double magnitude = charge * p.get_charge() / pow(distance, 2);
		double norm = 0.0;

		for (int i = 0; i < 3; i++) {
			direction[i] = p.get_coordinates()[i] - x[i];
			norm += pow(direction[i], 2);
		}
		for (int i = 0; i < 3; i++) { F[i] += direction[i] * magnitude / sqrt(norm); }

	}
}

inline double particle::get_distance(particle& p1) {
	double output = 0.0;
	for (int i = 0; i < 3; i++) {
		output += pow(x[i] - p1.get_coordinates()[i], 2);
	}
	return sqrt(output);
}

