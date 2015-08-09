#include <iostream>
#include <random>	
#include <sstream>
#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include "population.h"
#include "defines.h"


using namespace std;
namespace spd = spdlog;


#define LOGD clog->debug()

void print_trait_counts(TraitFrequencies* tf,std::shared_ptr<spdlog::logger>& log) {
	// print with loci as rows, traits as columns
	int* locus_counts = tf->trait_counts;
	for(int locus = 0; locus < tf->numloci; locus++) {
		std::stringstream s;
		s << "locus " << locus << ": ";
		for(int trait = 0; trait < tf->max_num_traits; trait++) {
			s << setw(4) << locus_counts[locus * tf->max_num_traits + trait] << " ";
		}
		SPDLOG_DEBUG(log,"{}",s.str());
	}
}



int main(int argc, char** argv) {
	std::string VERSION = "0.0.1";
	int popsize;
	double innovrate;
	int numloci;
	int simlength;
	int inittraits;
	std::string logfile;
	int debug;
	std::random_device rd;
	std::mt19937_64 mt(rd());
	std::uniform_real_distribution<double> uniform(0.0, 1.0);
	TraitFrequencies* tf;

	auto clog = spdlog::stdout_logger_mt("console");

	clog->info() << "Neutral Cultural Transmission in C++ Framework Version: " <<  VERSION;


	try {
		TCLAP::CmdLine cmd("Neutral Cultural Transmission in C++ Framework", ' ', VERSION);

		TCLAP::ValueArg<int> p("p","popsize","Population size",true,100,"integer");
		TCLAP::ValueArg<int> nl("l", "numloci","Number of independent dimensions/loci to evolve within the population",true,4,"integer");
		TCLAP::ValueArg<double> ir("i", "innovrate","Innovation rate per time step per individual (e.g., 0.1 equals 10 percent change of an innovation per time step",true,0.001,"double");
		TCLAP::ValueArg<int> len("s", "simlength","Length of the simulation in generations of popsize individuals",true,1000,"integer");
		TCLAP::ValueArg<int> it("t","inittraits","Number of initial traits present at each dimension/locus",true,4,"integer");
		TCLAP::ValueArg<int> d("d", "debug", "Set debugging level, with 0 or absence indicating debug output is off",false,0,"integer");
		TCLAP::ValueArg<std::string> f("f","logfile","Path to log file and filename (e.g., /tmp/test.log",false,"","string");


		cmd.add(p);
		cmd.add(nl);
		cmd.add(ir);
		cmd.add(len);
		cmd.add(it);
		cmd.add(d);
		cmd.add(f);
		cmd.parse( argc, argv );

		popsize = p.getValue();
		numloci = nl.getValue();
		inittraits = it.getValue();
		innovrate = ir.getValue();
		simlength = len.getValue();
		debug = d.getValue();
		logfile = f.getValue();

	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; return 1; }


	if(debug < 1) { spd::set_level(spd::level::info); }
	else { spd::set_level(spd::level::debug); }




	Population* pop = new Population(popsize, numloci, inittraits, innovrate, mt, clog);
	SPDLOG_DEBUG(clog, "Constructed population: {}", pop->dbg_params());
	pop->initialize();


	tf = pop->tabulate_trait_counts();
	print_trait_counts(tf,clog);

	SPDLOG_DEBUG(clog,"Evolving population for {} steps", simlength);

	for(int i = 0; i < simlength; i++)
		pop->step();

	tf = pop->tabulate_trait_counts();
	print_trait_counts(tf,clog);








    return 0;
}