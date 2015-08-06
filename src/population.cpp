#include <iostream>
#include <random>
#include <boost/format.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>



#include "population.h"


Population::~Population() {

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


	std::cout << "debug: initialized population (rows are loci, columns are indiv)" << std::endl;

	//print with loci as rows, columns as indiv
	for(int locus = 0; locus < numloci; locus++) {
		std::cout << "locus " << locus << ": ";
		for(int indiv = 0; indiv < popsize; indiv++) {
			std::cout << population_traits[locus * popsize + indiv] << " ";
		}
		std::cout << std::endl;
	}



}


void Population::test_random() {
	for(int i = 0; i < 20; i++) SPDLOG_DEBUG(log, "testing generator: {}", this->uniform_pop(this->mt));
}

std::string Population::dbg_print() {
	boost::format fmt("[Population %4% | popsize: %1% numloci: %2% inittraits: %3%]");
	fmt % this->popsize % this->numloci % this->inittraits % this;
	return fmt.str();
}


