#include <random>
#include <spdlog/spdlog.h>

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

		SPDLOG_DEBUG(log, "initializing count array {:p} as {}x{} block with size {}", (void*)trait_counts, numloci, max_num_traits,size_count_array);
		memset(trait_counts, 0, (size_count_array * sizeof(int)));
	}
	~TraitFrequencies() { 
		SPDLOG_DEBUG(log,"deallocating block trait_counts {:p}", (void*)trait_counts); 
#if defined(__INTEL_COMPILER)
		_mm_free(trait_counts);
#else
		free(trait_counts);
#endif
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
	int popsize;
	int numloci;
	int inittraits;
	std::uniform_int_distribution<int> uniform_pop;
	std::mt19937_64 mt;
	std::shared_ptr<spdlog::logger> log;
	std::vector<int> next_trait;
	int* population_traits;
	int* prev_population_traits;
	int* locus_counts;

	void swap_population_arrays();


public:
	Population(int p,
		int n,
		int i,
		std::mt19937_64 m,
		std::shared_ptr<spdlog::logger>& l) : popsize(p), numloci(n), inittraits(i), mt(m), log(l)
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
	* Counts the occurrence of traits in the current population of individuals, for each locus/dimension, storing
	* the results in a rectangular array of integers, with loci as rows and columns for each trait.  The return value
	* is a TraitFrequencies object, which simply exposes the array as a public data member along with the dimensions
	* of the array (along with managing the count array memory for proper cleanup).
	*/
	TraitFrequencies* tabulate_trait_counts();

	/**
	* Advances the simulation by one time step, implementing cultural transmission within the population.  
	* Implements a cultural transmission algorithm within the population.  First, the method copies the state 
	* of the population for reference as the "previous" population traits.  Then, transmission and innovation 
	* occurs.  
	*/
	void step();



	std::string dbg_params();
	void dbg_log_population();
};