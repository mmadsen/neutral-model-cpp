#include <iostream>
#include <sstream>
#include <iomanip>
#include <spdlog/spdlog.h>
#include "statistics.h"

namespace spd = spdlog;


namespace CTModels {

void print_trait_counts(TraitFrequencies* tf,std::shared_ptr<spdlog::logger>& log) {
	// skip if logging level isn't high enough
	if(log->level() == spd::level::trace) {

		// print with loci as rows, traits as columns
		int* locus_counts = tf->trait_counts;
		for(int locus = 0; locus < tf->numloci; locus++) {
			std::stringstream s;
			s << "locus " << locus << ": ";
			for(int trait = 0; trait < tf->max_num_traits; trait++) {
				s << std::setw(4) << locus_counts[locus * tf->max_num_traits + trait] << " ";
			}
			SPDLOG_TRACE(log,"{}",s.str());
		}
	}
}


void print_trait_statistics(TraitStatistics* ts, std::shared_ptr<spdlog::logger>& log) {
	// skip if logging level isn't high enough
	if(log->level() == spd::level::trace || log->level() == spd::level::debug) {
	

		int* locus_richness = ts->trait_richness_by_locus;
		for(int locus = 0; locus < ts->numloci; locus++) {
			std::stringstream s;
			s << "richness @ locus: " << locus << ": " << locus_richness[locus];
			SPDLOG_DEBUG(log,"{}",s.str());
		}
	}
}






};

