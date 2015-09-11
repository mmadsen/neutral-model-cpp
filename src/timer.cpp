#include <chrono>
#include <ratio>
#include <vector>
#include <string>
#include <spdlog/spdlog.h>
#include <spdlog/logger.h>
#include "timer.h"
#include "globals.h"

using namespace CTModels;

namespace CTModels {

void Timer::start(std::string label) {
	auto start_point = std::chrono::high_resolution_clock::now();
	start_times[label] = start_point;
}

void Timer::end(std::string label) {
	auto end_point = std::chrono::high_resolution_clock::now();
	auto start = start_times[label];
	completed_times[label] = std::chrono::duration <double, std::milli>(end_point - start).count();
}

double Timer::interval_ms(std::string label) {
	return completed_times[label];
}

std::vector<std::string> Timer::get_timed_events() {
	std::vector<std::string> labels;
	for(auto it = completed_times.cbegin(); it != completed_times.cend(); ++it) {
		labels.push_back(it->first);
	}
	return labels;
}


};