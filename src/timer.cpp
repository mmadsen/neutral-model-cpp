#include <chrono>
#include <ratio>
#include "timer.h"

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


};