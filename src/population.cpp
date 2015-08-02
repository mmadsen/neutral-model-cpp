#include <iostream>
#include <random>
#include <boost/format.hpp>

#include "population.h"




void Population::initialize() {

}

void Population::set_population_size(int popsize) {
	this->popsize = popsize;
}
void Population::set_numloci(int numloci) {
	this->numloci = numloci;
}
void Population::set_inittraits(int inittraits) {
	this->inittraits = inittraits;
}

std::string Population::dbg_print() {
	boost::format fmt("[Population %4% | popsize: %1% numloci: %2% inittraits: %3%]");
	fmt % this->popsize % this->numloci % this->inittraits % this;
	return fmt.str();
}


