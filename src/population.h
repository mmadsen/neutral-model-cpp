#pragma once

#include <random>
#include <spdlog/spdlog.h>
#include "defines.h"
#include "statistics.h"





namespace CTModels {

/** \class Population 
*
* Represents a population of individuals, which carry cultural traits along one or more dimensions
* or loci.  Traits can be limited or unlimited (practically speaking) in number.  A population object
* stores the current state of the population, and the state of the population in the immediately 
* previous time step.  A step() method implements the transmission algorithm for a given simulation.  
*
*/

class Population {
private:
	int popsize;
	int numloci;
	int inittraits;
	double innovation_rate;
	std::uniform_int_distribution<int> uniform_pop;
	std::uniform_int_distribution<int> uniform_locus;
	std::poisson_distribution<int> poisson_dist;
	std::mt19937_64 mt;
	std::shared_ptr<spdlog::logger> log;
	std::vector<int> next_trait;
	int* population_traits;
	int* prev_population_traits;
	int* locus_counts;
	int* indiv_to_copy;
	int trait_digits_printing = 0;
	int pop_digits_printing = 0;
	TraitFrequencies* current_trait_counts;

	void swap_population_arrays();


public:
	Population(int p,
		int n,
		int i,
		double r,
		std::shared_ptr<spdlog::logger>& l) : popsize(p), numloci(n), inittraits(i), innovation_rate(r), log(l)
	{}

	~Population();

	/**
	* Initializes a population given the population size, number of loci, and other values given at construction.
	* In this initial implementation, each individual gets uniform random integer values at numloci dimensions where
	* traits are constrained to be between [0, inittraits).  
	*/
	void initialize();

	/**
	* Tabulates frequencies of traits in the current population of individuals, separately for each locus/dimension.
	* Stores the trait counts as a TraitFrequencies object for later use.
	*/
	void tabulate_trait_counts();


	/**
	* Return the current trait counts.  The return value
	* is a TraitFrequencies object, which simply exposes the array as a public data member along with the dimensions
	* of the array (along with managing the count array memory for proper cleanup).
	*
	*/
	TraitFrequencies* get_current_trait_counts();

	/**
	* Calculate derived statistics from the current set of trait counts.  
	*
	*/
	TraitStatistics* calculate_trait_statistics();


	/**
	* Advances the simulation by one time step, implementing cultural transmission within the population.  
	* Implements a cultural transmission algorithm within the population.  The basic
	* WF model has no innovation.  
	*/
	void step_basicwf();

	/**
	* Advances the simultion by one time step, implementing cultural transmission in the population
	* along with innovation, using the infinite-alleles model of Kimura to simulate trait dimensions 
	* where there are no constraints on new innovations and every innovation is a new variant which has not
	* previously existed.  
	*
	*/
	void step_wfia();



	std::string dbg_params();
	void dbg_log_population();
};

};
