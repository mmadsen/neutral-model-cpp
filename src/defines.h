

#if defined(__INTEL_COMPILER)
	#define FREE(x) _mm_free(x) 
#else
	#define FREE(x) free(x)
#endif