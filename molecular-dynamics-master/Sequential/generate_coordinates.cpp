#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <ctime>
#include <cstdlib>
using namespace std;

int main(int argc, char* argv[]) {
    if(argc == 3){

        int num_particles = atoi(argv[1]); 
        double universe_size[3] = {atof(argv[2]),atof(argv[2]),atof(argv[2])}; 
        /*randomly generated coordinates of particles*/

        random_device rd;
        mt19937 rand(rd());
        uniform_real_distribution<> generate(0.1, universe_size[0]-0.1);
        
        double mass = 1.0;
        double charge = -1.0;
        ofstream output;
        output.open("initial_coordinates.txt");
        output << num_particles <<"    " <<  universe_size[0] <<endl;
        for (int i = 0; i < num_particles; i++) {
            output << generate(rand) << "   " << generate(rand) << "    " << generate(rand) <<"    " <<  mass <<"    " << charge<<endl;
        }
        output.close();
    }
    else{
        cout<<"Not enough input"<<endl;
    }
    return 0;
}
