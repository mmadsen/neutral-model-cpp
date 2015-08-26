#include <random>
#include <omp.h>

namespace CTModels {


void generate_uniform_int(int begin, int end, int num_variates, int* variates) {
#pragma omp parallel shared(variates)
{
	std::random_device rd;
	std::mt19937_64 eng(rd());
	std::uniform_int_distribution<int> uniform_int{begin, end - 1};

	#pragma omp for
	for(int i = 0; i < num_variates; i++) {
		variates[i] = uniform_int(eng);
	}
} // end pragma omp parallel

} // end function

}; // end namespace