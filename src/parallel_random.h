#pragma once


namespace CTModels {

/** \function generate_uniform_int 
*
* Simple wrapper function that encapsulates producing a large number of random variates from 
* the C++11 standard library using OpenMP parallelization to improve wall-clock performance.      
*
*/

void generate_uniform_int(int begin, int end, int num_variates, int* variates);


};