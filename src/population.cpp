#include <iostream>
#include <iomanip>
#include <random>
#include <algorithm>
#include <sstream>
#include <math.h>
#include <boost/format.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>
#include <omp.h>

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
	// Initialize needed random number generators
	std::random_device rd_mt;
	std::random_device rd_pois;
	std::random_device rd_locus;
	std::mt19937_64 t_mt(rd_mt());
	std::mt19937_64 t_mt_pois(rd_pois());
	std::mt19937_64 t_mt_locus(rd_locus());

	this->mt = t_mt;
	this->mt_pois = t_mt_pois;
	this->mt_locus = t_mt_locus;


	// in debug printing, we want fixed columns, with the number of digits appropriate given the 
	// max number of traits or population size
	char buffer[32];
	pop_digits_printing = sprintf(buffer, "%ld", (long)popsize);
	//SPDLOG_DEBUG(log,"Using {} digits to print population individual ID's", pop_digits_printing);


	// Construct a uniform integer distribution which will draw individuals from the population
	// std::random is not necessarily thread safe, so if we use this with OpenMP tasks or loops,
	// the strategy will be to draw as many as we need from the loop in a #pragma openmp single s
	// section, and then block them to the worker team to use.  
	std::uniform_int_distribution<int> u(0, this->popsize - 1);
	this->uniform_pop = u;

	// Construct a uniform distribution for choosing a random locus
	std::uniform_int_distribution<int> l(0, this->numloci - 1);
	this->uniform_locus = l;

	// Construct a Poisson distribution for innovation rates, with mean popsize * innovrate
	double mutation_rate = static_cast<double>(this->popsize) * this->innovation_rate; 
	//SPDLOG_DEBUG(log,"Constructing Poisson for innovation with mean {:0.4f}", mutation_rate);
	std::poisson_distribution<int> p(mutation_rate);
	this->poisson_dist = p;

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


	// for(int locus = 0; locus < numloci; locus++) {
	// 	for(int indiv = 0; indiv < popsize; indiv++) {
	// 		population_traits[indiv * numloci + locus] = initial_trait_dist(this->mt);
	// 	}
	// }


	// For the first generation only, the previous population is the same as the initial population
	memcpy(prev_population_traits, population_traits, ((popsize * numloci) * sizeof(int)));

}



void Population::tabulate_trait_counts() {

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
			int trait_at_locus = population_traits[indiv * numloci + locus];
			++locus_counts[locus * largest_locus_value + trait_at_locus];
		}
	}

	this->current_trait_counts = tf;
}

TraitFrequencies* Population::get_current_trait_counts() {
	return this->current_trait_counts;
}


TraitStatistics* Population::calculate_trait_statistics() {
	TraitStatistics* ts = new TraitStatistics(this->numloci);
	TraitFrequencies* tf = this->current_trait_counts;

	int* locus_counts = tf->trait_counts;
	for(int locus = 0; locus < tf->numloci; locus++) {
		int richness = 0;

		for(int trait = 0; trait < tf->max_num_traits; trait++) {
			if(locus_counts[locus * tf->max_num_traits + trait] > 0) 
				++richness;
		}
		ts->trait_richness_by_locus[locus] = richness;
	}


	return ts;
}



void Population::step_basicwf() {
	//SPDLOG_DEBUG(log, "Entering Population::step, implementing simple Wright-Fisher copying");
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

	// //DEBUG
	// std::stringstream s;
	// s << "random indiv: ";
	// for(int i = 0; i < popsize; i++) {
	// 	s << indiv_to_copy[i] << " ";
	// }
	// SPDLOG_DEBUG(log,"{}",s.str());


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


	// clean up 
	FREE(indiv_to_copy);

}


void Population::step_wfia() {
	//SPDLOG_DEBUG(log, "Entering Population::step_wfia, implementing infinite-alleles Wright-Fisher copying");
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

	// //DEBUG
	// std::stringstream s;
	// s << "random indiv: ";
	// for(int i = 0; i < popsize; i++) {
	// 	s << indiv_to_copy[i] << " ";
	// }
	// SPDLOG_DEBUG(log,"{}",s.str());


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
	int num_mutations = poisson_dist(this->mt_pois);
	SPDLOG_TRACE(log,"WFIA: num mutations this step: {}", num_mutations);

	for(int j = 0; j < num_mutations; j++) {
		int indiv_to_mutate = uniform_pop(this->mt);
		int locus_to_mutate = uniform_locus(this->mt_locus);
		int new_trait = next_trait[locus_to_mutate];
		++next_trait[locus_to_mutate];
		population_traits[indiv_to_mutate * numloci + locus_to_mutate] = new_trait;
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
	boost::format fmt("[Population %4% | popsize: %1% numloci: %2% inittraits: %3%  innovation_rate: %5%]");
	fmt % this->popsize % this->numloci % this->inittraits % this % this->innovation_rate;
	return fmt.str();
}

void Population::dbg_log_population() {
	SPDLOG_DEBUG(log, "population state: (rows are individuals, columns are loci)");

	//print with indiv as rows, columns as loci
	for(int indiv = 0; indiv < popsize; indiv++) {
		std::stringstream s;
		s << "indiv: " <<  std::setw(pop_digits_printing) << indiv << ": ";
		for(int locus = 0; locus < numloci; locus++) {
			s << population_traits[indiv * numloci + locus] << " ";
		}
		SPDLOG_DEBUG(log,"{}",s.str());
	}
}

