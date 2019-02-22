#include <cstdio>

class particle {
public:
	
	particle() { 
        //x = 0.0; y = 0.0; z = 0.0; vx = 0.0; vy = 0.0; vz = 0.0;
        for(int i = 0; i<3; i++){ coord[i] = 0.0;vel[i] = 0.0; acce[i]=0.0;}
        }
	particle(double x_, double y_, double z_,double vx_, double vy_, double vz_, double ax_,double ay_, double az_) { coord[0] = x_; coord[1] = y_; coord[2] = z_; vel[0]  = vx_; vel[1] = vy_; vel[2] = vz_;acce[0] = ax_;acce[1] = ay_;acce[2]= az_;}
	//void print() { printf("x = (%f, %f, %f)\tv = (%f, %f, %f)\n", x, y, z,vx,vy,vz);  }
        double coord[3];	
        double vel[3];
        double acce[3];
	//double x, y, z, vx,vy,vz;
};

void resize_array(particle* array, int len) {
	particle* new_array = new particle[len + 1000];
	for (int i = 0; i < len; i++) { new_array[i] = array[i]; }
	delete[] array;
	array = new_array;
}
