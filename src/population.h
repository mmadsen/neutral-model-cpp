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
		trait_counts = (int*) malloc(size_count_array * sizeof(int));
		SPDLOG_DEBUG(log, "initializing count array {:p} as {}x{} block with size {}", (void*)trait_counts, numloci, max_num_traits,size_count_array);
		memset(trait_counts, 0, (size_count_array * sizeof(int)));
	}
	~TraitFrequencies() { 
		SPDLOG_DEBUG(log,"deallocating block trait_counts {:p}", (void*)trait_counts); 
		free(trait_counts);
	}
};


class Population {
	int popsize;
	int numloci;
	int inittraits;
	std::uniform_int_distribution<int> uniform_pop;
	std::mt19937_64 mt;
	std::shared_ptr<spdlog::logger> log;
	std::vector<int> next_trait;
	int* population_traits;
	int* locus_counts;

public:
	Population(int p,
		int n,
		int i,
		std::mt19937_64 m,
		std::shared_ptr<spdlog::logger>& l) : popsize(p), numloci(n), inittraits(i), mt(m), log(l)
	{}
	~Population();
	void initialize();

	TraitFrequencies* tabulate_trait_freq();

	void step();



	std::string dbg_params();
	void dbg_log_population();
};