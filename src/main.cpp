#include <iostream>
#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>

#include "population.h"

#define LOGD clog->debug()

using namespace std;
namespace po = boost::program_options;
namespace spd = spdlog;


int main(int argc, char** argv) {
	po::variables_map vm;
	int popsize;
	double innovrate;
	int numloci;
	int simlength;
	int inittraits;
	std::string logfile;
	int debug;

	auto clog = spd::stdout_logger_mt("console");

	clog->info() << "Neutral Cultural Transmission in C++ Framework V.0.1";

	try {
		po::options_description desc("Allowed options");
		desc.add_options()
		    ("help", "produce help message")
		    ("popsize", po::value<int>(), "Population size")
		    ("innovrate", po::value<double>(), "Innovation rate per time step per individual (e.g., 0.1 equals 10 percent chance of a mutation")
		    ("numloci", po::value<int>(), "Number of independent dimensions or loci to evolve within the population")
		    ("simlength", po::value<int>(), "Length of the simulation in elemental copying steps")
		    ("inittraits", po::value<int>(), "Number of initial traits present at each dimension/locus")
		    ("logfile", po::value<std::string>()->default_value("logfile.txt"), "Path for log file and filename (e.g., /tmp/test.log")
		    ("debug", po::value<int>(), "set debugging level");


		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);    

		if(vm.count("help")) {
	    	cout << desc << endl;
	    	return 1;
		}

		if(vm.count("debug")) {
			debug = vm["debug"].as<int>();
			if(debug < 1) { spd::set_level(spd::level::info); }
			else { spd::set_level(spd::level::debug); }

		}

		popsize = vm["popsize"].as<int>();
		innovrate = vm["innovrate"].as<double>();
		numloci = vm["numloci"].as<int>();
		simlength = vm["simlength"].as<int>();
		inittraits = vm["inittraits"].as<int>();

	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}

	clog->info("Constructing population");
	Population* pop = new Population();
	pop->set_population_size(popsize);
	pop->set_numloci(numloci);
	pop->set_inittraits(inittraits);

	SPDLOG_DEBUG(clog, pop->dbg_print())






    return 0;
}