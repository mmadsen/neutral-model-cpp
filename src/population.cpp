#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <boost/format.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>
#include <memory>

#include "population.h"
#include "statistics.h"
#include "defines.h"
#include "timer.h"
#include "globals.h"
#include "parallel_random.h"

using namespace CTModels;

extern CTModels::Timer timer;

namespace CTModels {


/**
* Population destructor.  Normally not reached if we run a single population and then exit from main(),
* but in situations where we might run simulations in a loop, we don't want to leak memory for whole 
* populations. 
*/
Population::~Population() {
	SPDLOG_TRACE(clog,"deallocating block prev_population_traits {:p}", (void*)prev_population_traits); 
	SPDLOG_TRACE(clog,"deallocating block population_traits {:p}", (void*)population_traits); 
	SPDLOG_TRACE(clog,"deallocating block indiv_to_copy {:p}", (void*)indiv_to_copy);
	FREE(prev_population_traits);
	FREE(population_traits);
	FREE(indiv_to_copy);
}


void Population::initialize() {
	timer.start("population::initialize");
	// Initialize needed random number generators
	std::random_device rd_mt;
	std::mt19937_64 t_mt(rd_mt());

	this->mt = t_mt;


	// in debug printing, we want fixed columns, with the number of digits appropriate given the 
	// max number of traits or population size
	char buffer[32];
	pop_digits_printing = sprintf(buffer, "%ld", (long)popsize);
	//SPDLOG_DEBUG(clog,"Using {} digits to print population individual ID's", pop_digits_printing);


	// Construct a uniform integer distribution which will draw individuals from the population
	// std::random is not necessarily thread safe, so if we use this with OpenMP tasks or loops,
	// the strategy will be to draw as many as we need from the loop in a #pragma openmp single s
	// section, and then block them to the worker team to use.  
	std::uniform_int_distribution<int> u{0, this->popsize - 1};
	this->uniform_pop = u;

	// Construct a uniform distribution for choosing a random locus
	std::uniform_int_distribution<int> l{0, this->numloci - 1};
	this->uniform_locus = l;

	// Construct a Poisson distribution for innovation rates, with mean popsize * innovrate
	double mutation_rate = static_cast<double>(this->popsize) * this->innovation_rate; 
	//SPDLOG_DEBUG(log,"Constructing Poisson for innovation with mean {:0.4f}", mutation_rate);
	std::poisson_distribution<int> p{mutation_rate};
	this->poisson_dist = p;

	// next_trait stores the next new mutation/innovation for each locus/dimension, and each slot
	// is incremented when a new trait is handed out.  Since there are "inittraits" in the initial 
	// population, we initialize this array of values to inittraits + 1.  
	for(int i = 0; i < this->numloci; i++) {
		next_trait.push_back(this->inittraits + 1);
	}

	auto trait_bufsize = (numloci * popsize) * sizeof(int);
	population_traits = (int*) ALIGNED_MALLOC(trait_bufsize);
	prev_population_traits = (int*) ALIGNED_MALLOC(trait_bufsize);
	//SPDLOG_DEBUG(clog, "Pop initializing pop traits array {:p} as {}x{} block with size {}", (void*)population_traits, popsize, numloci, trait_bufsize);
	//SPDLOG_DEBUG(clog, "Pop initializing prev pop array {:p} as {}x{} block with size {}", (void*)prev_population_traits, popsize, numloci, trait_bufsize);

	int num_variates = popsize * numloci;
	generate_uniform_int(0, inittraits - 1, num_variates, population_traits);

	// Initialize a buffer to hold random numbers indicating which individuals are copied
	// in each time step
	auto indiv_bufsize = popsize * sizeof(int);
	indiv_to_copy = (int*) ALIGNED_MALLOC(indiv_bufsize);
	//SPDLOG_DEBUG(clog, "Pop initializing indiv_to_copy array {:p} as {} block with size {}", (void*)indiv_to_copy, popsize, indiv_bufsize);



	// For the first generation only, the previous population is the same as the initial population
	memcpy(prev_population_traits, population_traits, trait_bufsize);
	timer.end("population::initialize");
}



std::shared_ptr<TraitFrequencies> Population::tabulate_trait_counts() {
	timer.start("population::tabulate_trait_counts");
	// allocate space for the largest value in any locus
	// array of counts will be a rectangular array numloci * largest_locus_value
	// technically, the largest_locus_value is 1 greater than any trait value seen at any locus, 
	// so we trimmed by one
	// MEM:  dynamically allocated locus_counts is freed in the destructor of TraitFrequencies

	auto result = std::max_element(next_trait.begin(), next_trait.end());
	int largest_locus_value = *result - 1;

	// declared as a std::unique_ptr, because we want the TF object from the last time tabulate was called
	// to clean itsetf up once there isn't a reference anymore.
	std::shared_ptr<TraitFrequencies> tf;
	tf.reset(new TraitFrequencies(numloci,largest_locus_value));
	int* locus_counts = tf->trait_counts;

#pragma vector always
	for(int indiv = 0; indiv < popsize; indiv++) {
		for(int locus = 0; locus < numloci; locus++) {
			int trait_at_locus = population_traits[indiv * numloci + locus];
			++locus_counts[locus * largest_locus_value + trait_at_locus];
		}
	}
	timer.end("population::tabulate_trait_counts");
	return tf;
}






void Population::step_basicwf() {
	// Prepare by copying current state to previous state, before doing transmission 
	// algorithm
	swap_population_arrays();

	generate_uniform_int(0, popsize, popsize, indiv_to_copy);


	// Basic Wright-Fisher dynamics without innovation
#pragma omp parallel 
{
	int indiv;
	#pragma omp for private(indiv)
	for(indiv = 0; indiv < popsize; indiv++) {
		int tocopy = indiv_to_copy[indiv];

		for(int locus = 0; locus < numloci; locus++) {
			population_traits[indiv * numloci + locus] = prev_population_traits[tocopy * numloci + locus];
		}

	}
}


}


void Population::step_wfia() {
	// Prepare by copying current state to previous state, before doing transmission 
	// algorithm
	swap_population_arrays();

	generate_uniform_int(0, popsize, popsize, indiv_to_copy);

	// Wright-Fisher dynamics - in order to optimize performance, first we do all the copying so that
	// most of the copies will be vectorized and broken into work units.  Then, we come back for a small
	// number of individuals and give them mutated traits.  

#pragma omp parallel 
{
	int indiv;
	#pragma omp for private(indiv)
	for(indiv = 0; indiv < popsize; indiv++) {
		int tocopy = indiv_to_copy[indiv];

		for(int locus = 0; locus < numloci; locus++) {
			population_traits[indiv * numloci + locus] = prev_population_traits[tocopy * numloci + locus];
		}

	}
}

	// Now, we create innovations given the innovation rate, randomly throughout the population
	int num_mutations = poisson_dist(this->mt);
	//SPDLOG_TRACE(clog,"WFIA: num mutations this step: {}", num_mutations);

	for(int j = 0; j < num_mutations; j++) {
		int indiv_to_mutate = uniform_pop(this->mt);
		int locus_to_mutate = uniform_locus(this->mt);
		int new_trait = next_trait[locus_to_mutate];
		++next_trait[locus_to_mutate];
		population_traits[indiv_to_mutate * numloci + locus_to_mutate] = new_trait;
	}



}




void Population::swap_population_arrays() {
	// Start by swapping the population trait arrays, so that we capture the previous state for use

	int* tmp = prev_population_traits;
	prev_population_traits = population_traits;
	population_traits = tmp;

	// Clear the population_traits to zero so we can fill it by transmission from prev_population_traits
	int size_pop_array = numloci * popsize;
	memset(population_traits, 0, (size_pop_array * sizeof(int))); 
}



std::string Population::dbg_params() {
	boost::format fmt("[Population %4% | popsize: %1% numloci: %2% inittraits: %3%  innovation_rate: %5%]");
	fmt % this->popsize % this->numloci % this->inittraits % this % this->innovation_rate;
	return fmt.str();
}

void Population::dbg_log_population() {
	SPDLOG_TRACE(clog, "population state: (rows are individuals, columns are loci)");

	//print with indiv as rows, columns as loci
	for(int indiv = 0; indiv < popsize; indiv++) {
		std::stringstream s;
		s << "indiv: " <<  std::setw(pop_digits_printing) << indiv << ": ";
		for(int locus = 0; locus < numloci; locus++) {
			s << population_traits[indiv * numloci + locus] << " ";
		}
		SPDLOG_TRACE(clog,"{}",s.str());
	}
}

};
