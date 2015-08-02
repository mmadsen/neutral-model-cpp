
class Population {
	int popsize = 0;
	int numloci = 0;
	int inittraits = 0;

public:
	void initialize();
	void set_population_size(int popsize);
	void set_numloci(int numloci);
	void set_inittraits(int inittraits);
	std::string dbg_print();

};