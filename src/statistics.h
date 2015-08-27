#pragma once
#include <spdlog/spdlog.h>
#include "defines.h"
#include "globals.h"



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

	TraitFrequencies(int n, int m) : numloci(n), max_num_traits(m) {
		int count_bufsize = (numloci * max_num_traits) * sizeof(int);
		trait_counts = (int*) ALIGNED_MALLOC(count_bufsize);

		//SPDLOG_DEBUG(clog, "TF initializing count array {:p} as {}x{} block with size {}", (void*)trait_counts, numloci, max_num_traits,count_bufsize);
		memset(trait_counts, 0, count_bufsize);
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
		int bufsize = numloci * sizeof(int);
		trait_richness_by_locus = (int*) ALIGNED_MALLOC(bufsize); 
		//SPDLOG_DEBUG(clog, "TF initializing richness array {:p} as {} block with size {}", (void*)trait_richness_by_locus, numloci, bufsize);

	}

	~TraitStatistics() {
		FREE(trait_richness_by_locus);
	}
};


/* Defined in statistics.cpp */
void print_trait_counts(TraitFrequencies* tf);
void print_trait_statistics(TraitStatistics* ts);
void print_event_timing();





};
