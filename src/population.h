#include <random>
#include <spdlog/spdlog.h>


class Population {
	int popsize;
	int numloci;
	int inittraits;
	std::uniform_int_distribution<int> uniform_pop;
	std::mt19937_64 mt;
	std::shared_ptr<spdlog::logger> log;

public:
	Population(int p,
		int n,
		int i,
		std::mt19937_64 m,
		std::shared_ptr<spdlog::logger>& l) : popsize(p), numloci(n), inittraits(i), mt(m), log(l)
	{}
	void initialize();

	void test_random();

	std::string dbg_print();

};