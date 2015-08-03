#include <iostream>
#include <random>
#include <boost/format.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>


#include "population.h"



void Population::initialize() {
	std::uniform_int_distribution<int> u(0, this->popsize - 1);
	this->uniform_pop = u;
}


void Population::test_random() {
	for(int i = 0; i < 20; i++) SPDLOG_DEBUG(log, "testing generator: {}", this->uniform_pop(this->mt));
}

std::string Population::dbg_print() {
	boost::format fmt("[Population %4% | popsize: %1% numloci: %2% inittraits: %3%]");
	fmt % this->popsize % this->numloci % this->inittraits % this;
	return fmt.str();
}


