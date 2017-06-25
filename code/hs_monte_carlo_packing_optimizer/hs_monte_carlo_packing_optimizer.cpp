#include <map>
#include <set>
#include <list>
#include <cmath>
#include <ctime>
#include <deque>
#include <queue>
#include <stack>
#include <string>
#include <bitset>
#include <cstdio>
#include <limits>
#include <vector>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <array>
#include <random>
#include <regex>

using namespace std;

const double pi = 3.1415926535897;

//convenient object to describe a particle
struct Particle {
	double rx, ry, rz; 
	int type;
	int index;
};

//convenient object to store the details of the simulation
struct State {
	//explicitly loaded in
	int N; //number of particles
	double L; //box edge length
	//unordered_map<int, double> type_to_diameter; //maps type to diameter
	//unordered_map<int, string> type_to_atom; //maps type to atom for using VMD
	vector< double > type_to_diameter;
	vector< string > type_to_atom;
	vector< Particle > particles; //stores the particles
	//created via initialization function
	int num_types;
	double eta = 0;
	vector< int > type_to_N_type; //maps particle type to the number of them

	//does any initialization needed after loading in relevant details
	void Initialize(){
		//count the number of each type present
		Particle particle;
		num_types = type_to_diameter.size();
		type_to_N_type.clear();
		type_to_N_type.resize(num_types, 0);
		for (auto it = particles.begin(); it < particles.end(); it++){
			particle = *it;
			eta = eta + pi*pow(type_to_diameter[particle.type], 3) / 6.0;
			type_to_N_type[particle.type]++;
		}
		eta = eta / pow(L, 3);
	}
	
	//eventually store the cell list here
};

//function definitions
bool CheckParticleOverlap(State &state, Particle &particle);
void BinaryRandomSequentialAddition(State &state, double eta, int N, double d0, double d1, int max_attempts = 100000);
void WriteConfig(string file, State &state);
void WriteState(string file, State &state);
void ReadState(string file, State &state);
bool AttemptParticleTranslation(State &state, Particle &particle, double dx, double dy, double dz);
void MonteCarlo(State &state, long int steps, double dr_max);

//this will drive everything based on command line input values
int main(){
	//the simulation state
	State state;
	State state_new;

	//generate an initial state via binary random addition
	BinaryRandomSequentialAddition(state, 0.35, 100, 1.0, 1.0);
	state.Initialize();
	cout << state.eta << endl;

	//serialize a state
	WriteState("brsa_conf.state", state);

	//read a serialized state
	ReadState("brsa_conf.state", state_new);
	state_new.Initialize();
	cout << state.eta << endl;

	for (auto it = state.type_to_N_type.begin(); it < state.type_to_N_type.end(); it++){
		cout << distance(state.type_to_N_type.begin(), it) << "," << *it << endl;
	}

	//write an xyz file
	WriteConfig("brsa_conf.xyz", state_new);

	//move particles
	MonteCarlo(state_new, 100000, 0.2);
	
	return 0;
}

//checks if a specified particle overlaps with another via brute force methods
bool CheckParticleOverlap(State &state, Particle &particle){
	Particle existing_particle;
	double squared_distance;
	double min_squared_distance;
	double dx, dy, dz;

	//iterate over existing particles to see if overlap
	bool overlap_free = true;
	for (auto it = state.particles.begin(); it < state.particles.end(); it++){
		existing_particle = *it;
		dx = particle.rx - existing_particle.rx;
		dy = particle.ry - existing_particle.ry;
		dz = particle.rz - existing_particle.rz;
		dx = dx - round(dx / state.L) * state.L;
		dy = dy - round(dy / state.L) * state.L;
		dz = dz - round(dz / state.L) * state.L;

		squared_distance = pow(dx, 2) + pow(dy, 2) + pow(dz, 2);

		//check for overlap and if found push iterator to end
		min_squared_distance = pow(state.type_to_diameter[particle.type] + state.type_to_diameter[existing_particle.type], 2) / 4.0;
		if (squared_distance < min_squared_distance){
			if (particle.index != existing_particle.index){
				it = state.particles.end() - 1;
				overlap_free = false;
			}
		}
	}
	return overlap_free;
}

//using random sequential addition it attempts to generate a packing at some fixed volume fraction (no cell lists used yet)
void BinaryRandomSequentialAddition(State &state, double eta, int N, double d0, double d1, int max_attempts){
	double N0 = int(double(N)/2.0);
	double N1 = N - N0;
	state.N = N;
	state.L = cbrt((N0 / eta)*pow(d0, 3)*pi / 6.0 + (N1 / eta)*pow(d1, 3)*pi / 6.0);
	double d_small = (int)(d0 <= d1)*d0 + (int)(d1 < d0)*d1;
	double d_large = (int)(d0 >= d1)*d0 + (int)(d1 > d0)*d1;
	state.type_to_diameter.push_back(d_small); state.type_to_atom.push_back("A");
	state.type_to_diameter.push_back(d_large); state.type_to_atom.push_back("B");

	//random number generator
	mt19937_64 rng_x, rng_y, rng_z;
	rng_x.seed(12); rng_y.seed(100); rng_z.seed(3);
	uniform_real_distribution<double> unif(0, state.L);

	//items for attempting insertions
	int attempt, n;
	bool overlap_free;
	int particle_type;
	Particle particle;
	//array<double, 4> new_particle, existing_particle;

	//try to insert particles
	attempt = 0; n = 0;
	while (attempt < max_attempts && n < N){
		particle_type = (int)(n < N0) * 0 + (int)(n >= N0) * 1;
		particle.rx = unif(rng_x); particle.ry = unif(rng_y); particle.rz = unif(rng_z);
		particle.type = particle_type;
		particle.index = n;

		//check for an overlap
		overlap_free = CheckParticleOverlap(state, particle);

		//add new particle if possible
		if (overlap_free){
			attempt = 0; n++; state.particles.push_back(particle);
		}
		else{
			attempt++; 
		}

		//provide some feedback
		if (n % 100 == 0)
			cout << "Completed insertions: " << n << endl;

	}
	//provide some feedback
	if (attempt < max_attempts)
		cout << "BRSA completed!" << endl;

	getchar();
}

//take the state and write out an xyz file of it
void WriteConfig(string file, State &state){
	ofstream output(file, ios::app);
	output << state.N << endl;
	output << " " << "written by ryan" << endl;

	Particle particle;
	string atom;
	for (auto it = state.particles.begin(); it < state.particles.end(); it++){
		particle = *it;
		atom = state.type_to_atom[particle.type] + "       ";
		output << "  " << atom << particle.rx << "       " << particle.ry << "       " << particle.rz << endl;
	}
	output.close();
}

//serialize the state object
void WriteState(string file, State &state){
	ofstream output(file, ios::app);
	output << state.N << " " << state.L << endl; //save the N and L
	for (auto it = state.type_to_atom.begin(); it < state.type_to_atom.end(); it++){ //save the map from type to atom
		output << "(" << distance(state.type_to_atom.begin(), it) << "," << *it << ")" << " ";
	}
	output << endl;
	for (auto it = state.type_to_diameter.begin(); it < state.type_to_diameter.end(); it++){ //save the map from type to atom
		output << "(" << distance(state.type_to_diameter.begin(), it) << "," << *it << ")" << " ";
	}
	output << endl;
	Particle particle;
	//string atom_type;
	for (auto it = state.particles.begin(); it < state.particles.end(); it++){
		particle = *it;
		output << "(" << particle.rx << "," << particle.ry << "," << particle.rz << "," << particle.type << "," << particle.index << ")" << " ";
	}
	output << endl;
	output.close();
}

//read in the serialized state object
void ReadState(string file, State &state){
	ifstream in_file;
	string line;
	in_file.open(file);
	if (!in_file) {cerr << "Unable to open state file.txt"; exit(1);}

	//prepare the regex's and match
	regex N_and_L_re("([0-9]+)\\s+([0-9\\.]+)");
	regex type_and_atom_re("\\(([0-9]+),([a-zA-Z]+)\\)"); 
	regex type_and_diameter_re("\\(([0-9]+),([0-9\\.]+)\\)");
	regex particles_re("\\(([0-9\\.]+),([0-9\\.]+),([0-9\\.]+),([0-9]+),([0-9]+)\\)");
	sregex_iterator next, end;
	smatch match;

	//read in N and L
	getline(in_file, line);
	regex_search(line, match, N_and_L_re);
	state.N = stoi(match.str(1)); state.L = stod(match.str(2));

	//read in the type to atom
	getline(in_file, line);
	//regex_search(line, match, type_and_atom_re);
	next = sregex_iterator(line.begin(), line.end(), type_and_atom_re);
	end = sregex_iterator();
	state.type_to_atom.clear();
	while (next != end) {
		match = *next;
		cout << match.str(1) << "," << match.str(2) << "\n";
		//state.type_to_atom[stoi(match.str(1))] = match.str(2);
		state.type_to_atom.push_back(match.str(2));
		next++;
	}

	//read in the type to diameter
	getline(in_file, line);
	next = sregex_iterator(line.begin(), line.end(), type_and_diameter_re);
	end = sregex_iterator();
	state.type_to_diameter.clear();
	while (next != end) {
		match = *next;
		cout << match.str(1) << "," << match.str(2) << "\n";
		//state.type_to_diameter[stoi(match.str(1))] = stod(match.str(2));
		state.type_to_diameter.push_back(stod(match.str(2)));
		next++;
	}

	//read in particles
	getline(in_file, line);
	next = sregex_iterator(line.begin(), line.end(), particles_re);
	end = sregex_iterator();
	Particle particle;
	state.particles.clear();
	while (next != end) {
		match = *next;
		particle.rx = stod(match.str(1)); particle.ry = stod(match.str(2)); particle.rz = stod(match.str(3));
		particle.type = stoi(match.str(4));
		particle.index = stoi(match.str(5));
		state.particles.push_back(particle);
		next++;
	}
}

//attempts to translate a particle with periodic wrapping
bool AttemptParticleTranslation(State &state,  int index, double dx, double dy, double dz){
	Particle particle_translated = state.particles[index];
	particle_translated.rx = particle_translated.rx + dx;
	particle_translated.rx = particle_translated.rx - floor(particle_translated.rx / state.L)*state.L;
	particle_translated.ry = particle_translated.ry + dy;
	particle_translated.ry = particle_translated.ry - floor(particle_translated.ry / state.L)*state.L;
	particle_translated.rz = particle_translated.rz + dz;
	particle_translated.rz = particle_translated.rz - floor(particle_translated.rz / state.L)*state.L;
	//check if overlap and make move if possible
	if (CheckParticleOverlap(state, particle_translated)){
		state.particles[particle_translated.index] = particle_translated;
		return true;
	}
	else
		return false;
}

//attempts to change the particle type and adjust box volume as need be
/*bool AttemptParticleTypeChange(State &state, int index, int r_type_change){
	State particle_translated = state;


	unordered_map<int, int> type_to_N_type; //maps type to number of this type present



}*/

//perform monte carlo steps
void MonteCarlo(State &state, long int steps, double dr_max)
{
	//random number generator
	mt19937_64 rng_move_type, rng_index, rng_x, rng_y, rng_z, rng_type;
	rng_move_type.seed(1); rng_index.seed(991); rng_x.seed(34); rng_y.seed(103); rng_z.seed(333); rng_type.seed(6);
	uniform_int_distribution<int> r_index(0, state.N - 1);
	uniform_real_distribution<double> r_dr(-dr_max, dr_max);
	uniform_int_distribution<int> r_type_change(0, 1); //0=shrink, 1=grow
	bool translation_status;

	//loop and move particles randomly
	Particle particle_translated;
	bool overlap_free = true;
	for (long int i = 0; i < steps; i++){
		//random translation move
		translation_status = AttemptParticleTranslation(state, r_index(rng_index), r_dr(rng_x), r_dr(rng_y), r_dr(rng_z));

		//random particle type change
		//type_change_status = AttemptParticleTypeChange(state, r_index(rng_index), r_type_change);

		if (i % 1000 == 0){
			cout << i << endl;
			WriteConfig("trajectory.xyz", state);
		}
	}
}