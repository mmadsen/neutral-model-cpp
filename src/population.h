#include <random>
#include <spdlog/spdlog.h>
#include "defines.h"

/** \class TraitFrequencies
*
* TraitFrequencies bundles a rectangular array of integers, formatted as a simple 1-D 
* array for performance, with the dimensions of the rectangular array.  The array 
* has loci as rows and traits within a locus as columns, so the array is addressed
* as:  trait_counts[locus * max_num_traits + trait]
*
*/

class TraitFrequencies {
public:
	int* trait_counts; 
	int numloci;
	int max_num_traits;
	std::shared_ptr<spdlog::logger> log;

	TraitFrequencies(int n, int m, std::shared_ptr<spdlog::logger>& l) : numloci(n), max_num_traits(m), log(l) {
		int size_count_array = numloci * max_num_traits;

#if defined(__INTEL_COMPILER)
		trait_counts = (int*) _mm_malloc(size_count_array * sizeof(int),64);
#else
		trait_counts = (int*) malloc(size_count_array * sizeof(int));
#endif

		//SPDLOG_TRACE(log, "initializing count array {:p} as {}x{} block with size {}", (void*)trait_counts, numloci, max_num_traits,size_count_array);
		memset(trait_counts, 0, (size_count_array * sizeof(int)));
	}
	~TraitFrequencies() { 
		//SPDLOG_TRACE(log,"deallocating block trait_counts {:p}", (void*)trait_counts); 
		FREE(trait_counts);
	}
};



/** \class TraitStatistics
*
* TraitStatistics bundle derived statistics from raw trait counts. 
*/

class TraitStatistics {
public:
	int* trait_richness_by_locus;
	int numloci;

	TraitStatistics(int numloci) : numloci(numloci) {
#if defined(__INTEL_COMPILER)
		trait_richness_by_locus = (int*) _mm_malloc(numloci * sizeof(int), 64);
#else
		trait_richness_by_locus = (int*) malloc(numloci * sizeof(int));
#endif
	}

	~TraitStatistics() {
		FREE(trait_richness_by_locus);
	}
};




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
	std::mt19937_64 mt_pois;
	std::mt19937_64 mt_locus;
	std::shared_ptr<spdlog::logger> log;
	std::vector<int> next_trait;
	int* population_traits;
	int* prev_population_traits;
	int* locus_counts;
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