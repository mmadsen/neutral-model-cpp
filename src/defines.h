#pragma once

#if defined(__INTEL_COMPILER)
	#define FREE(x) _mm_free(x) 
#else
	#define FREE(x) free(x)
#endif

#if defined(__INTEL_COMPILER)
		#define ALIGNED_MALLOC(x) _mm_malloc(x, 64)
#else
		#define ALIGNED_MALLOC(x) malloc(x)
#endif

