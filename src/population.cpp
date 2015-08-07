#include <iostream>
#include <random>
#include <algorithm>
#include <sstream>
#include <boost/format.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>

#include "population.h"


Population::~Population() {
	free(population_traits);
}


void Population::initialize() {
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

	population_traits = (int*) malloc((numloci * popsize) * sizeof(int));

	std::uniform_int_distribution<int> initial_trait_dist(0, inittraits - 1);

	for(int indiv = 0; indiv < popsize; indiv++) {
		for(int locus = 0; locus < numloci; locus++) {
			population_traits[locus * popsize + indiv] = initial_trait_dist(this->mt);
		}
	}

}


TraitFrequencies* Population::tabulate_trait_freq() {

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







// void Population::test_random() {
// 	for(int i = 0; i < 20; i++) SPDLOG_DEBUG(log, "testing generator: {}", this->uniform_pop(this->mt));
// }

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

