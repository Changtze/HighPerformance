#include "simulation.h"
#include "physics.h"
#include "particle.h"
#include <iostream>
#include <boost/program_options.hpp>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cmath>

using namespace std;
namespace po = boost::program_options;


int main(int argc, char *argv[]){
	// SPECIFY AVAILABLE USER OPTIONS
	
	po::options_description opts(
	"Allowed options");
	double value;
	bool is_set = false;
	
	// AVAILABLE COMMAND-LINE ARGUMENTS
	opts.add_options()
	("help", "Print available options.")
	("Lx", po::value<double>()->default_value(20),
	"x length (Angstroms)")
	("Ly", po::value<double>()->default_value(20),
	"y length (Angstroms)")
	("Lz", po::value<double>()->default_value(20),
	"z length (Angstroms)")
	("dt", po::value<double>()->default_value(0.001),
	"Time-step")
	("T", po::value<double>(), "Final time")
	("ic-one", "Initial condition 1: one stationary particle")
	("ic-one-vel", "Initial condition 2: one moving particle")
	("ic-two", "Initial condition 3: two bouncing particles")
	("ic-two-pass1", "Initial condition 4: two passing particles")
	("ic-two-pass2", "Initial condition 5: two bouncing particles, close")
	("ic-two-pass3", "Initial condition 6: two passing particles, close, large")
	("ic-random", "Initial condition: N random particles")
	("percentage-type1", po::value<double>(&value)->notifier([&](double) { is_set = true; })->default_value(10.0),
	 "Percentage of type 1 particles with random IC")
	("N", po::value<int>(),
	 "Number of particels to spawn")
	("temp", po::value<double>(),
	 "Temperature (degrees Kelvin)");
	
	// PARSE CLI ARGUMENTS
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, opts), vm);
	po::notify(vm);
	
	// CHECK FOR --HELP OPTION
	if (vm.count("help")) {
		cout << opts << endl;
		return 1;
	}
	
	
	//// CLI ERROR CHECKING ////
	
	// validity of initial conditions
	vector<string> ic_options = {"ic-one", "ic-one-vel",
								 "ic-two-pass1", "ic-two-pass2",
								 "ic-two-pass3", "ic-random",
								 "ic-two"
								};
								 
	int ic_count = 0;  // counter to check for --ic-x inputs
	string ic_opt;  // initialising container for selected initial condition
	
	
	for (const auto& option : ic_options) {
		if (vm.count(option)) {
			ic_count++;  // increment for every -ic cli found
			//cout << option << endl;  // debugging purposes
			ic_opt = option;
		}
	}
	// ic_count check
	if (ic_count == 0) {
		cerr << "error: no initial condition provided." << endl;
		return 1;
	} else if (ic_count > 1) {
		cerr << "error: expected only 1 initial condition." << endl;
		return 1;
	}
	
	// Checking initial conditions to assign N
	string substr_1 = "one";
	string substr_2 = "two";	
	string substr_3 = "random";
	size_t found_1 = ic_opt.find(substr_1);
	size_t found_2 = ic_opt.find(substr_2);
	size_t found_3 = ic_opt.find(substr_3);
	
	int N;
	// check for N and percentange-type1 if --ic-random is provided
	if (ic_opt == "ic-random") {
		if (!vm.count("N")) {
			cerr << "error: N expected for initial condition. not parsed" << endl;
			return 1;
		} else {
			if (found_1 != string::npos) {
				N = 1;
			} else if (found_2 != string::npos) {
				N = 2;
			}  else {
				N = vm["N"].as<int>();
			}
		}
		if (!is_set){
			cerr << "error: percentage-type1 expected for random initial condition. not parsed." << endl;
			return 1;
		}
//		if (!vm.count("percentage-type1")) {
//			
//			cerr << "error: percentage-type1 expected for random initial condition. not parsed." << endl;
//			return 1;
//		} 
	}
	
	
	// Other options are mandatorydicate if temperature is fixed by user
	vector<string> mandatory_opts = {"Lx", "Ly", "Lz", "dt", "T"};
	for (const auto& option : mandatory_opts) {
		if (!vm.count(option)) {
			cerr << "error: " << option << " expected. not parsed." << endl;
			return 1;
		}
	}
	//// CLI ERROR CHECKING COMPLETE ////
	
	// Assigning simulation variables
	const auto Lx = vm["Lx"].as<double>();
	const auto Ly = vm["Ly"].as<double>();
	const auto Lz = vm["Lz"].as<double>();
	const double delta_t = vm["dt"].as<double>();
	const double T = vm["T"].as<double>();
	cout << "Final time: " << T << endl;
	cout << "X domain size: " << Lx << endl;
	cout << "Y domain size: " << Ly << endl;
	cout << "Z domain size: " << Lz << endl;
	cout << "Time step: " << delta_t << endl;

	// allow hard-coding of particle types for test cases
	// N = 1;
	// N_type0 = 1;
	// N_type1 = 1;
	//	cout << "Type 1 particle count: " << N_type1 << endl;
	//	cout << "Type 0 particle count: " << 

	
	// Beginning simulation
	cout << "Command-line inputs parsed successfully. Beginning simulation... " << endl;
	Simulation sim = Simulation(Lx, Ly, Lz, T, delta_t);
	
	if (vm.count("temp")) {
		sim.setFixTemp();
		sim.setTargetTemp(vm["temp"].as<double>());
	}
	
	// Initialise the correct intial condition
	// go through ic_options vector
	
	// Initialising test case 1
	if (ic_opt == "ic-one"){
		sim.TestCase1();
		cout << "Test case 1 initialised. " << endl;
	} else if (ic_opt == "ic-one-vel") {
		sim.TestCase2();
		cout << "Test case 2 initialised. " << endl;			
	} else if (ic_opt == "ic-two") {
		sim.TestCase3();
		cout << "Test case 3 initialised. " << endl;			
	} else if (ic_opt == "ic-two-pass1") {
		sim.TestCase4();
		cout << "Test case 4 initialised. " << endl;			
	} else if (ic_opt == "ic-two-pass2") {
		sim.TestCase5();
		cout << " Test case 5 initialised. " << endl;
	} else if (ic_opt == "ic-two-pass3") {
		sim.TestCase6();
		cout << "Test case 6 initialised. " << endl;			
	} else {
		cout << "Random test case initialised. " << endl;
		sim.Random(N, vm["percentage-type1"].as<double>());
	}

	// Begin simulation
	sim.simulate();
	
	cout << "Simulation successful." << endl;
	
	return 0;
	
}