#include <iostream>
#include <random>
#include <algorithm>
#include <sstream>
#include <boost/format.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>

#include "population.h"
#include "defines.h"

/**
* Population destructor.  Normally not reached if we run a single population and then exit from main(),
* but in situations where we might run simulations in a loop, we don't want to leak memory for whole 
* populations. 
*/
Population::~Population() {
	SPDLOG_DEBUG(log,"deallocating block prev_population_traits {:p}", (void*)prev_population_traits); 
	SPDLOG_DEBUG(log,"deallocating block population_traits {:p}", (void*)population_traits); 
	FREE(prev_population_traits);
	FREE(population_traits);
}


void Population::initialize() {

	// in debug printing, we want fixed columns, with the number of digits appropriate given the 
	// max number of traits or population size
	char buffer[32];
	pop_digits_printing = sprintf(buffer, "%ld", (long)popsize);
	SPDLOG_DEBUG(log,"Using {} digits to print population individual ID's", pop_digits_printing);


	// Construct a uniform integer distribution which will draw individuals from the population
	// std::random is not necessarily thread safe, so if we use this with OpenMP tasks or loops,
	// the strategy will be to draw as many as we need from the loop in a #pragma openmp single s
	// section, and then block them to the worker team to use.  
	std::uniform_int_distribution<int> u(0, this->popsize - 1);
	this->uniform_pop = u;

	// next_trait stores the next new mutation/innovation for each locus/dimension, and each slot
	// is incremented when a new trait is handed out.  Since there are "inittraits" in the initial 
	// population, we initialize this array of values to inittraits + 1.  
	for(int i = 0; i < this->numloci; i++) {
		next_trait.push_back(this->inittraits + 1);
	}

	// initialize population_traits.  We could use the new C++11 syntax
	// int* population_traits = new int[X][Y] but we also want the block aligned
	// and the access is slow because technically it's returning a pointer from the 
	// first array access before it does the second.  Access this array in one of 
	// two ways:
	//		population_traits[y * width + x] 
	//
	// which since I want to parallelize over individuals, I will interpret as:
	//
	// 		{individual X trait at locus Y} = population_traits[Y * popsize + X]
#if defined(__INTEL_COMPILER)
	population_traits = (int*) _mm_malloc((numloci * popsize) * sizeof(int), 64);
	prev_population_traits = (int*) _mm_malloc((numloci * popsize) * sizeof(int),64);
#else
	population_traits = (int*) malloc((numloci * popsize) * sizeof(int));
	prev_population_traits = (int*) malloc((numloci * popsize) * sizeof(int));
#endif

	std::uniform_int_distribution<int> initial_trait_dist(0, inittraits - 1);

	for(int indiv = 0; indiv < popsize; indiv++) {
		for(int locus = 0; locus < numloci; locus++) {
			population_traits[locus * popsize + indiv] = initial_trait_dist(this->mt);
		}
	}

	// For the first generation only, the previous population is the same as the initial population
	memcpy(prev_population_traits, population_traits, ((numloci * popsize) * sizeof(int)));

}



TraitFrequencies* Population::tabulate_trait_counts() {

	// allocate space for the largest value in any locus
	// array of counts will be a rectangular array numloci * largest_locus_value
	// technically, the largest_locus_value is 1 greater than any trait value seen at any locus, 
	// so we trimmed by one
	// MEM:  dynamically allocated locus_counts is freed in the destructor of TraitFrequencies

	auto result = std::max_element(next_trait.begin(), next_trait.end());
	int largest_locus_value = *result - 1;


	TraitFrequencies* tf = new TraitFrequencies(numloci,largest_locus_value,log);
	int* locus_counts = tf->trait_counts;

	for(int indiv = 0; indiv < popsize; indiv++) {
		for(int locus = 0; locus < numloci; locus++) {
			int trait_at_locus = population_traits[locus * popsize + indiv];
			++locus_counts[locus * largest_locus_value + trait_at_locus];
		}
	}

	return tf;


}


void Population::step() {
	SPDLOG_DEBUG(log, "Entering Population::step, implementing simple Wright-Fisher copying");
	// Prepare by copying current state to previous state, before doing transmission 
	// algorithm
	swap_population_arrays();

	// the RNG is not safe to use in parallel, and it's a hassle to give every thread an RNG,
	// so we generate a list of n=popsize individuals which will serve as the targets of copying
	// this list can then be broken into blocks for threads in an OpenMP team
	int* indiv_to_copy;

#if defined(__INTEL_COMPILER)
	indiv_to_copy = (int*) _mm_malloc(popsize * sizeof(int), 64);
#else
	indiv_to_copy = (int*) malloc(popsize * sizeof(int));
#endif

	for(int i = 0; i < popsize; i++) {
		indiv_to_copy[i] = uniform_pop(this->mt);
	}

	//DEBUG
	std::stringstream s;
	s << "random indiv: ";
	for(int i = 0; i < popsize; i++) {
		s << indiv_to_copy[i] << " ";
	}
	SPDLOG_DEBUG(log,"{}",s.str());

	//TODO: it really might be better for whole-individual copying algorithms to have rows be individuals
	// and columns be loci, since otherwise copying an individual is non-unit-stride across multiple rows
	// of the prev_population_traits directory and insertion is non-unit-stride.

	// Basic Wright-Fisher dynamics without innovation
	for(int i = 0; i < popsize; i++) {

	}


	// clean up 
	FREE(indiv_to_copy);

}


void Population::swap_population_arrays() {
	// Start by swapping the population trait arrays, so that we capture the previous state for use
	//SPDLOG_DEBUG(log, "preswap - prev_population_traits: {:p}", (void*)prev_population_traits);
	//SPDLOG_DEBUG(log, "preswap - population_traits: {:p}", (void*)population_traits);
	int* tmp = prev_population_traits;
	prev_population_traits = population_traits;
	population_traits = tmp;

	//SPDLOG_DEBUG(log, "postswap - prev_population_traits: {:p}", (void*)prev_population_traits);
	//SPDLOG_DEBUG(log, "postswap - population_traits: {:p}", (void*)population_traits);

	// Clear the population_traits to zero so we can fill it by transmission from prev_population_traits
	int size_pop_array = numloci * popsize;
	memset(population_traits, 0, (size_pop_array * sizeof(int))); 
}



std::string Population::dbg_params() {
	boost::format fmt("[Population %4% | popsize: %1% numloci: %2% inittraits: %3%]");
	fmt % this->popsize % this->numloci % this->inittraits % this;
	return fmt.str();
}

void Population::dbg_log_population() {
	SPDLOG_DEBUG(log, "population state: (rows are loci, columns are indivuals)");

	//print with loci as rows, columns as indiv
	for(int locus = 0; locus < numloci; locus++) {
		std::stringstream s;
		s << "locus " << locus << ": ";
		for(int indiv = 0; indiv < popsize; indiv++) {
			s << population_traits[locus * popsize + indiv] << " ";
		}
		SPDLOG_DEBUG(log,"{}",s.str());
	}
}

