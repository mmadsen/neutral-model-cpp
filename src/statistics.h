#pragma once
#include <spdlog/spdlog.h>
#include "defines.h"

namespace CTModels {

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


/* Defined in statistics.cpp */
void print_trait_counts(TraitFrequencies* tf,std::shared_ptr<spdlog::logger>& log);
void print_trait_statistics(TraitStatistics* ts, std::shared_ptr<spdlog::logger>& log);






};
