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
#include <set>

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
	double V;
	double eta;
	vector< int > type_to_N_type; //maps particle type to the number of them
	vector< double > type_to_diameter_ratio; //maps particle type diameter ratio with respect to smallest
	//created via call to cell list builder
	//int N_cells;
	//double L_cell;
	//vector< vector< vector<double> > > cells;

	//does any initialization needed after loading in relevant details
	void Initialize(){
		//calculate the volume
		V = pow(L, 3);

		//count the number of each type present
		Particle particle;
		num_types = type_to_diameter.size();
		type_to_N_type.clear();
		type_to_N_type.resize(num_types, 0);
		eta = 0.0;
		for (auto it = particles.begin(); it < particles.end(); it++){
			particle = *it;
			eta = eta + pi*pow(type_to_diameter[particle.type], 3) / 6.0;
			type_to_N_type[particle.type]++;
		}

		//calculate the volume fraction
		eta = eta / V;

		//build a vector or particle diameter ratios with respect to smallest (0'th)
		type_to_diameter_ratio.clear();
		//type_to_diameter_ratio.resize(num_types, 0.0);
		for (auto it = type_to_diameter.begin(); it < type_to_diameter.end(); it++)
			type_to_diameter_ratio.push_back((double)(*it) / type_to_diameter[0]);
	}
};

// this is a cell list structure
struct CellList{
	int N_cells = -1;
	double L_cell;
	bool active = false;
	vector< vector< vector< list<int> > > > cells;

	//structure the cell list
	void PrepareCellList(State &state){
		int N_cells_prev = N_cells;
		bool active_prev = active;
		N_cells = (int)(state.L / state.type_to_diameter.back());
		active = N_cells > 3;
		if (active){
			//if not the same number of cells or if they have not been maintained do not rebuild
			if (N_cells != N_cells_prev || !active_prev){
				//build the cells
				int cell_x, cell_y, cell_z;
				cells.clear();
				L_cell = state.L / (float)N_cells;
				cells.resize(N_cells);
				for (cell_x = 0; cell_x < N_cells; cell_x++)
				{
					cells[cell_x].resize(N_cells);
					for (cell_y = 0; cell_y < N_cells; cell_y++)
					{
						cells[cell_x][cell_y].resize(N_cells);
					}
				}
				//load particles into the cells
				for (auto it = state.particles.begin(); it < state.particles.end(); it++){
					tie(cell_x, cell_y, cell_z) = FindCell((*it));
					cells[cell_x][cell_y][cell_z].push_back(distance(state.particles.begin(), it));
				}
			}
		}
	}

	//finds the cell a particle belongs to (never gets used outside of struct so no "active" check)
	tuple<int, int, int> FindCell(Particle &particle){
		return make_tuple((int)(particle.rx / L_cell), (int)(particle.ry / L_cell), (int)(particle.rz / L_cell));
	}

	//swaps a particle from one cell to another
	void SingleParticleUpdate(Particle &particle_new, Particle &particle_current){
		int cell_x_new, cell_y_new, cell_z_new;
		int cell_x_current, cell_y_current, cell_z_current;
		tie(cell_x_current, cell_y_current, cell_z_current) = FindCell(particle_new);
		//find the particle in it
		auto it = find(cells[cell_x_current][cell_y_current][cell_z_current].begin(),
			cells[cell_x_current][cell_y_current][cell_z_current].end(), particle_current.index);
		cells[cell_x_current][cell_y_current][cell_z_current].erase(it);
		cells[cell_x_new][cell_y_new][cell_z_new].push_back(particle_current.index);
	}
};

//function definitions
bool CheckParticleOverlapRSA(State &state, Particle &particle);
bool CheckParticleOverlap(State &state, CellList &cell_list, Particle &particle, bool advance_search = false);
bool CheckAllParticleOverlaps(State &state, CellList &cell_list);
void BinaryRandomSequentialAddition(State &state, double eta, int N, double d0, double d1, int max_attempts = 100000);
void WriteConfig(string file, State &state);
void WriteState(string file, State &state);
void ReadState(string file, State &state);
bool AttemptParticleTranslation(State &state, CellList &cell_list, Particle &particle, double dx, double dy, double dz);
bool AttemptParticleTypeChange(State &state, CellList &cell_list, int index, int r_type_change, double r_type_change_accept);
void WriteTypeStats(ofstream &type_stats_output, State &state);
void MonteCarlo(State &state, CellList &cell_list, long int equil_steps, long int prod_steps, long int skip_steps,
	string simulation_name, double dr_max, double frac_trans);

//this will drive everything based on command line input values
int main(){
	//the simulation state
	State state;
	State state_new;
	CellList cell_list;
	CellList cell_list_new;

	//generate an initial state via binary random addition
	BinaryRandomSequentialAddition(state, 0.25, 100, 1.0, 1.4);
	state.Initialize();
	cell_list.PrepareCellList(state);

	//serialize a state
	WriteState("brsa_conf.state", state);

	//read a serialized state
	ReadState("brsa_conf.state", state_new);
	state_new.Initialize();
	cout << "Volume fraction: " << state.eta << endl;
	cell_list.PrepareCellList(state);
	cout << "Number of cells along edge: " << cell_list.N_cells << endl;
	cout << "Cell list state: " << cell_list.active << endl;

	for (auto it = state.type_to_N_type.begin(); it < state.type_to_N_type.end(); it++){
		cout << distance(state.type_to_N_type.begin(), it) << "," << *it << endl;
	}

	//write an xyz file
	WriteConfig("brsa_conf.xyz", state_new);

	//move particles
	state_new.Initialize();
	cout << "t1: " << state.type_to_diameter_ratio[0] << endl;
	cout << "t2: " << state.type_to_diameter_ratio[1] << endl;
	getchar();
	MonteCarlo(state_new, cell_list_new, 0/*200000*/, 1000000, 100, "1_1.4", 0.2, 0.5);

	WriteState("final_conf.state", state_new);
	WriteConfig("final_conf.xyz", state_new);
	
	return 0;
}

//checks if a specified particle overlaps with another via brute force methods for RSA only!
bool CheckParticleOverlapRSA(State &state, Particle &particle){
	Particle existing_particle;
	double squared_distance;
	double min_squared_distance;
	double dx, dy, dz;

	//iterate over existing particles to see if overlap
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
				return false;
			}
		}
	}
	return true;
}

//checks if a specified particle overlaps with another via brute force methods
bool CheckParticleOverlap(State &state, CellList &cell_list, Particle &particle, bool advance_search){
	Particle existing_particle;
	double squared_distance;
	double min_squared_distance;
	double dx, dy, dz;

	//iterate over existing particles to see if overlap
	int advance = (int)advance_search*(particle.index + 1);
	for (auto it = state.particles.begin() + advance; it < state.particles.end(); it++){
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
				return false;
			}
		}
	}
	return true;
}

//checks for any possible overlap
bool CheckAllParticleOverlaps(State &state, CellList &cell_list){
	for (auto it = state.particles.begin(); it < state.particles.end(); it++){
		if (!CheckParticleOverlap(state, cell_list, *it, true))
			return false;
	}
	return true;
}

//using random sequential addition it attempts to generate a packing at some fixed volume fraction (no cell lists used yet)
void BinaryRandomSequentialAddition(State &state, double eta, int N, double d0, double d1, int max_attempts){
	double N0 = int(double(N) / 2.0);
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
		overlap_free = CheckParticleOverlapRSA(state, particle);

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
bool AttemptParticleTranslation(State &state, CellList &cell_list, int index, double dx, double dy, double dz){
	Particle particle_translated = state.particles[index];
	particle_translated.rx = particle_translated.rx + dx;
	particle_translated.rx = particle_translated.rx - floor(particle_translated.rx / state.L)*state.L;
	particle_translated.ry = particle_translated.ry + dy;
	particle_translated.ry = particle_translated.ry - floor(particle_translated.ry / state.L)*state.L;
	particle_translated.rz = particle_translated.rz + dz;
	particle_translated.rz = particle_translated.rz - floor(particle_translated.rz / state.L)*state.L;
	//check if overlap and make move if possible
	if (CheckParticleOverlap(state, cell_list, particle_translated)){
		state.particles[particle_translated.index] = particle_translated;
		return true;
	}
	else
		return false;
}

//attempts to change the particle type and adjust all diameters
bool AttemptParticleTypeChange(State &state, CellList &cell_list, int index, int r_type_change, double r_type_change_accept){
	//store the old info
	int type_current = state.particles[index].type;

	//get new possible type
	int type_new;
	if (r_type_change == 0)
		type_new = max(0, type_current - 1);
	else if (r_type_change == 1)
		type_new = min((int)state.type_to_N_type.size() - 1, type_current + 1);

	//calculate the acceptance probability for the mixing entropy correction
	double prob_mix_ent_accept = min(1.0, (double)(state.type_to_N_type[type_new] + 1) / (double)state.type_to_N_type[type_current]);
	//double prob_mix_ent_accept = 1.0; //REMOVE///////////////

	//accept or reject now so as to not waste time seeing if an overlap
	if (r_type_change_accept > prob_mix_ent_accept)
		return false;

	//make sure the change is not trivial as we just always accept (really just do nothing)
	if (type_current == type_new)
		return true;

	//register this change
	state.type_to_N_type[type_current]--;
	state.type_to_N_type[type_new]++;
	state.particles[index].type = type_new;

	//calculate needed rescale factor for diameters to maintain volume fraction
	double diameter_0 = 0.0;
	int type, Ni;
	for (auto it = state.type_to_N_type.begin(); it < state.type_to_N_type.end(); it++){
		type = distance(state.type_to_N_type.begin(), it); Ni = *it;
		diameter_0 = diameter_0 + (double)Ni*pow(state.type_to_diameter_ratio[type], 3.0);
	}
	diameter_0 = (pi / (6.0*state.V*state.eta))*diameter_0;
	diameter_0 = 1.0 / cbrt(diameter_0);

	//perform system wide rescale of particle diameters
	auto type_to_diameter_current = state.type_to_diameter;
	for (auto it = state.type_to_diameter_ratio.begin(); it < state.type_to_diameter_ratio.end(); it++){
		type = distance(state.type_to_diameter_ratio.begin(), it);
		state.type_to_diameter[type] = diameter_0*(*it);
	}

	//check for overlaps
	bool status = true;
	if (r_type_change == 0)
		status = CheckAllParticleOverlaps(state, cell_list);
	else if (r_type_change == 1)
		status = CheckParticleOverlap(state, cell_list, state.particles[index]);

	if (!status){
		state.type_to_N_type[type_current]++;
		state.type_to_N_type[type_new]--;
		state.particles[index].type = type_current;
		state.type_to_diameter = type_to_diameter_current;
		return false;
	}

	//if made it here then a valid type change was made
	return true;
}

//writes out particle type stats to file
void WriteTypeStats(ofstream &type_stats_output, State &state){
	type_stats_output << (double)state.type_to_N_type[0] / (double)state.N << endl;
}

//perform monte carlo steps
void MonteCarlo(State &state, CellList &cell_list, long int equil_steps, long int prod_steps, long int skip_steps, 
				string simulation_name, double dr_max, double frac_trans)
{
	//random number generator
	mt19937_64 rng_move_type, rng_index, rng_x, rng_y, rng_z, rng_type_change, rng_type_change_accept;
	rng_move_type.seed(1); rng_index.seed(991); rng_x.seed(34); rng_y.seed(103); rng_z.seed(333); rng_type_change.seed(6); rng_type_change_accept.seed(600);
	uniform_int_distribution<int> r_index(0, state.N - 1);
	uniform_real_distribution<double> r_dr(-dr_max, dr_max);
	uniform_int_distribution<int> r_type_change(0, 1); //0=shrink, 1=grow
	uniform_real_distribution<double> r_accept(0, 1);
	bool translation_status, type_change_status;
	double move_type;
	string filename;

	//equilibration cycle
	cout << "Equilibration cycle..." << endl << endl;
	for (long int i = 1; i <= equil_steps; i++){
		//choose a random type of move
		move_type = r_accept(rng_move_type);

		//random translation move
		if (move_type <= frac_trans)
			translation_status = AttemptParticleTranslation(state, cell_list, r_index(rng_index), r_dr(rng_x), r_dr(rng_y), r_dr(rng_z));
		//random particle type change
		else if (move_type > frac_trans)
			type_change_status = AttemptParticleTypeChange(state, cell_list, r_index(rng_index), r_type_change(rng_type_change), r_accept(rng_type_change_accept));

		//message user
		if (i % 10000 == 0){
			cout << "Completed " << i << " steps" << endl;
			WriteConfig("trajectory.xyz", state);
		}
	}

	//production cycle
	ofstream type_stats_output(simulation_name + "__type_stats.dat", ios::app);
	cout << endl << "Production cycle..." << endl << endl;
	for (long int i = 1; i <= prod_steps; i++){
		//choose a random type of move
		move_type = r_accept(rng_move_type);

		//random translation move
		if (move_type <= frac_trans)
			translation_status = AttemptParticleTranslation(state, cell_list, r_index(rng_index), r_dr(rng_x), r_dr(rng_y), r_dr(rng_z));
		//random particle type change
		else if (move_type > frac_trans)
			type_change_status = AttemptParticleTypeChange(state, cell_list, r_index(rng_index), r_type_change(rng_type_change), r_accept(rng_type_change_accept));

		//message user
		if (i % 10000 == 0){
			cout << "Completed " << i << " steps" << endl;
			WriteConfig("trajectory.xyz", state);
		}

		//write out trajectory stats
		if (i % skip_steps == 0){
			filename = simulation_name + "__type_stats.dat";
			WriteTypeStats(type_stats_output, state);
		}


	}
	type_stats_output.close();


}


//move particle and see if belongs to closest neighbors