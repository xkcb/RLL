// codebook.h

#ifndef CODEBOOK_H_
#define CODEBOOK_H_

#include <array>
#include <exception>
#include <iostream>
#include <memory>
#include <vector>

#include "CodewordMap.h"
#include "constants.h"


class Exception : public std::exception
{
public:
	explicit Exception(const char * message)
	{
		std::strncpy(msg_, message, 1024);
	}

	virtual ~Exception() throw (){}

	virtual const char * what() const throw ()
	{
		return msg_;
	}
protected:
	char msg_[1024];
};

struct SUPER_STATE
{
	int32_t      rds      		    = INVALID_RDS;
	EfmplusState state    	        = NO_STATE;
	uint16_t	 ntz	  		    = 0;
	int8_t		 sign	  		    = 0;

	uint32_t     sum_var_actual     = INVALID_SUM_VAR;
	uint32_t     sum_var_prediction = INVALID_SUM_VAR;
	uint16_t     codeword           = INVALID_SEQUENCE;

	SUPER_STATE(uint16_t codeword) : codeword(codeword)
	{

	}
} ;

template<typename T>
struct STATS
{
	std::vector<T> probability;
	std::vector<T> variance;

	STATS (int size)
		: probability (size, 0.0),
		  variance (size, 0.0)
	{

	}

	STATS(int size1, int size2)
	{
		std::vector<double> v1(size2, 0.0);
		std::vector<std::vector<double> > probability1 (size1, v1);

		std::vector<double> v2(size2, 0.0);
		std::vector<std::vector<double> >variance1 (size1, v2);
		probability = probability1;
		variance = variance1;

	}
} ;

struct STATS_2D
{
	std::vector<std::vector<double> > probability;
	std::vector<std::vector<double> > variance;

	STATS_2D(int size1, int size2)
		: probability(size1, std::vector<double>(size2, 0.0)),
		  variance(size1, std::vector<double>(size2, 0.0))
	{
	}
} ;

struct CodebookState
{
	uint16_t     codeword ;
	EfmplusState next_state ;
};


using CodebookArray = CodebookState[CDB_LINES][NUM_EFMPLUS_STATES];
//using CodebookArray = std::array<std::array<CodebookState, NUM_EFMPLUS_STATES>, CDB_LINES>;

class Codebook
{
public:
	Codebook(CodebookArray &codebook, std::shared_ptr<CodewordMap> codeword_map);

	void print_codebook_as_header(const char * name, int flag) const;
	bool calculate_pdf_with_transition_table();

	CodebookArray & line;

private:
	SUPER_STATE get_best_between_main_and_alternate(int byte, const SUPER_STATE & st);
	void get_best_between_states_0_3(int byte, const SUPER_STATE & st, SUPER_STATE & min_st);

	bool states_have_valid_codeword() const;
	bool is_consistent() const;
	bool is_valid() const;

	void print_line(FILE * f, int line_index) const;
	void print_transition_table();

	double collect_codebook_stats(std::vector<double> & p);
	bool initialise_transition_table();
	void get_probability_vector(std::vector<double> &p);


	SUPER_STATE	choose_next_super_state(int byte, const SUPER_STATE & st, int * switch_cnt);

	std::array< std::array<uint32_t, TRANSITION_TABLE_STATES >, TRANSITION_TABLE_STATES> transition;

	
	std::shared_ptr<CodewordMap> codeword_map;
	std::shared_ptr<STATS<double> > global_full_state_stats;
	std::vector<double> rds_probability;
};

#endif /* CODEBOOK_H_ */
