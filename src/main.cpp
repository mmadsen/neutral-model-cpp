#include <iostream>
#include <boost/program_options.hpp>
#include <spdlog/spdlog.h>

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
	std::string logfile;
	int debug;

	auto clog = spd::stdout_logger_mt("console");

	clog->info() << "TESTPROGRAM V.0.1";

	try {
		po::options_description desc("Allowed options");
		desc.add_options()
		    ("help", "produce help message")
		    ("popsize", po::value<int>(), "Population size")
		    ("innovrate", po::value<double>(), "Innovation rate per time step per individual (e.g., 0.1 equals 10 percent chance of a mutation")
		    ("numloci", po::value<int>(), "Number of independent traits or loci to evolve within the population")
		    ("simlength", po::value<int>(), "Length of the simulation in elemental copying steps")
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

	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}

	SPDLOG_DEBUG(clog, "popsize: {} innovrate: {} numloci: {} simlength: {}", popsize, innovrate, numloci, simlength)

    return 0;
}