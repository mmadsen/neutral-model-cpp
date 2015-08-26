#pragma once

#include <chrono>
#include <ratio>
#include <unordered_map>
#include <string>
#include <vector>


namespace CTModels {


/** \class Timer 
*
* Simple timer class for timing sections of code using the C++11 chrono library, and the system
* "steady clock".  Allows multiple events to be timed and the results cached for later use under an interval label.    
*
*/

class Timer {
private:
	std::unordered_map<std::string, double> completed_times;
	std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point>	start_times;  
	// end time cache not needed

public:
	/**
	* Begins recording a time interval, under a string label.  Call end() to stop the time interval,
	* and calculate and cache the time interval for later retrieval.
	*
	*/
	void start(std::string label);

	/**
	* Ends recording a time interval for a given string label.  
	*
	*/
	void end(std::string label);

	/**
	* Reports the time interval for a string label, using milliseconds as the time scale.
	*
	*/
	double interval_ms(std::string label);

	/**
	* Returns a vector of string labels for events with timing information.
	*
	*/
	std::vector<std::string> get_timed_events();

};


};