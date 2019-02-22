#include <cstdio>
#include <iostream>
#include <fstream>

class particle {
public:
	
	particle() { x = 0.0; y = 0.0; z = 0.0; vx = 0.0; vy = 0.0; vz = 0.0; }
	particle(double x_, double y_, double z_,double vx_, double vy_, double vz_) { x = x_; y = y_; z = z_; vx = vx_; vy = vy_; vz = vz_;}
	void print() { printf("x = (%f, %f, %f)\tv = (%f, %f, %f)\n", x, y, z,vx,vy,vz);  }
	
	double x, y, z, vx,vy,vz;
};

void resize_array(particle* array, int len) {
	particle* new_array = new particle[len + 1000];
	for (int i = 0; i < len; i++) { new_array[i] = array[i]; }
	delete[] array;
	array = new_array;
}

std::ostream& operator<< (std::ostream& out, particle& p) {
	out << '(' << p.x << ", " << p.y << ", " << p.z << ')';
	out << "\tv = (" << p.vx << ", " << p.vy << ", " << p.vz << ")\n";
}
