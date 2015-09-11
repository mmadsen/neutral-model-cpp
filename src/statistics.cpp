#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <spdlog/spdlog.h>

#include "statistics.h"
#include "globals.h"
#include "timer.h"

namespace spd = spdlog;
using namespace CTModels;



namespace CTModels {

void print_trait_counts(std::shared_ptr<TraitFrequencies> tf) {
	// skip if logging level isn't high enough
	if(clog->level() == spd::level::trace) {

		// print with loci as rows, traits as columns
		int* locus_counts = tf->trait_counts;
		for(int locus = 0; locus < tf->numloci; locus++) {
			std::stringstream s;
			s << "locus " << locus << ": ";
			for(int trait = 0; trait < tf->max_num_traits; trait++) {
				s << std::setw(4) << locus_counts[locus * tf->max_num_traits + trait] << " ";
			}
			SPDLOG_TRACE(clog,"{}",s.str());
		}
	}
}


void print_trait_statistics(std::shared_ptr<TraitStatistics> ts) {
	// skip if logging level isn't high enough
	if(clog->level() == spd::level::trace || clog->level() == spd::level::debug) {
	

		int* locus_richness = ts->trait_richness_by_locus;
		for(int locus = 0; locus < ts->numloci; locus++) {
			std::stringstream s;
			s << "richness @ locus: " << locus << ": " << locus_richness[locus];
			SPDLOG_DEBUG(clog,"{}",s.str());
		}
	}
}


void print_event_timing() {
	// skip if logging level isn't high enough
	if(clog->level() == spd::level::trace || clog->level() == spd::level::debug) {
		std::vector<std::string> events = timer.get_timed_events();
		for(auto it = events.begin(); it != events.end(); ++it) {
			std::stringstream s;
			s << "event " << *it << " time: " << timer.interval_ms(*it);
			SPDLOG_DEBUG(clog, "{}", s.str());
		}
	}
}


std::shared_ptr<TraitStatistics> calculate_trait_statistics(std::shared_ptr<TraitFrequencies> tf) {
	timer.start("statistics::calculate_trait_statistics");
	std::shared_ptr<TraitStatistics> ts(new TraitStatistics(tf->numloci));

	int* locus_counts = tf->trait_counts;
	for(int locus = 0; locus < tf->numloci; locus++) {
		int richness = 0;

		for(int trait = 0; trait < tf->max_num_traits; trait++) {
			if(locus_counts[locus * tf->max_num_traits + trait] > 0) 
				++richness;
		}
		ts->trait_richness_by_locus[locus] = richness;
	}

	timer.end("statistics::calculate_trait_statistics");
	return ts;
}



};

