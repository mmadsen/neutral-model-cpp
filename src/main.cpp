#include <iostream>
#include <random>	
#include <chrono>
#include <ratio>
#include <sstream>
#include <memory>

#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include "population.h"
#include "statistics.h"
#include "defines.h"
#include "timer.h"
#include "globals.h"

using namespace std;
using namespace CTModels;
namespace spd = spdlog;



enum ruletype { BASICWF, WFIA };

int main(int argc, char** argv) {
	std::string VERSION = "0.0.1";
	int popsize;
	double innovrate = 10000.0;
	int numloci;
	int simlength;
	int inittraits;
	std::string logfile;
	int debug;
	std::random_device rd;
	std::mt19937_64 mt(rd());
	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	ruletype rt;
	spdlog::level::level_enum debug_level;



	CTModels::clog->info() << "Neutral Cultural Transmission in C++ Framework Version: " <<  VERSION;


	try {
		vector<std::string> allowed_types;
		allowed_types.push_back("basicwf");
		allowed_types.push_back("wfia");
		TCLAP::ValuesConstraint<std::string> allowedVals( allowed_types );

		TCLAP::CmdLine cmd("Neutral Cultural Transmission in C++ Framework", ' ', VERSION);

		TCLAP::ValueArg<int> p("p","popsize","Population size",true,100,"integer");
		TCLAP::ValueArg<int> nl("l", "numloci","Number of independent dimensions/loci to evolve within the population",true,4,"integer");
		TCLAP::ValueArg<double> ir("i", "innovrate","Innovation rate per time step per individual (e.g., 0.1 equals 10 percent change of an innovation per time step",true,0.001,"double");
		TCLAP::ValueArg<int> len("s", "simlength","Length of the simulation in generations of popsize individuals",true,1000,"integer");
		TCLAP::ValueArg<int> it("t","inittraits","Number of initial traits present at each dimension/locus",true,4,"integer");
		TCLAP::ValueArg<int> d("d", "debug", "Set debugging level, with 0 or absence indicating debug output is off, 1 indicating debug, >1 indicating TRACE",false,0,"integer");
		TCLAP::ValueArg<std::string> t("r","ruletype", "Copying rule to use",true,"basicwf",&allowedVals);
		TCLAP::ValueArg<std::string> f("f","logfile","Path to log file and filename (e.g., /tmp/test.log",false,"","string");


		cmd.add(p);
		cmd.add(nl);
		cmd.add(ir);
		cmd.add(len);
		cmd.add(it);
		cmd.add(d);
		cmd.add(f);
		cmd.add(t);
		cmd.parse( argc, argv );

		popsize = p.getValue();
		numloci = nl.getValue();
		inittraits = it.getValue();
		innovrate = ir.getValue();
		simlength = len.getValue();
		debug = d.getValue();
		logfile = f.getValue();
		std::string rule = t.getValue();
		if(rule == "basicwf") {
			CTModels::clog->debug("Using basicwf ruletype");
			rt = BASICWF;
		}
		else if( rule == "wfia") {
			rt = WFIA;
			CTModels::clog->debug("Using wfia ruletype");
		}
		else {
			std::cerr << "ERROR: ruletype " << rt << std::endl; 
			return 1;
		}

	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{ std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl; return 1; }

	if(debug == 0) { debug_level = spd::level::info; }
	else if(debug == 1) { debug_level = spd::level::debug; }
	else { debug_level = spd::level::trace; }

	spd::set_level(debug_level);

	timer.start("main");

	Population* pop = new Population(popsize, numloci, inittraits, innovrate);
	SPDLOG_TRACE(CTModels::clog, "Constructed population: {}", pop->dbg_params());
	pop->initialize();


	auto tf = pop->tabulate_trait_counts();
	print_trait_counts(tf);

	SPDLOG_DEBUG(CTModels::clog,"Evolving population for {} steps", simlength);

	switch(rt) {
		case BASICWF :
			for(int i = 0; i < simlength; i++)
				pop->step_basicwf();
			break;
		case WFIA :
			for(int i = 0; i < simlength; i++) 
				pop->step_wfia();
			break;
	}


	auto tf2 = pop->tabulate_trait_counts();
	auto ts = calculate_trait_statistics(tf2);

	print_trait_statistics(ts);
	print_trait_counts(tf2);

	timer.end("main");

	print_event_timing();



    return 0;
}