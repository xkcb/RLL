// codebook.h

#ifndef CODEBOOK_H_
#define CODEBOOK_H_

#include <array>
#include <bitset>
#include <exception>
#include <iostream>
#include <memory>
#include <vector>

#include "CodewordMap.h"
#include "constants.h"

#include <blaze/math/StaticVector.h>
#include <blaze/math/StaticMatrix.h>

class Exception : public std::exception
{
public:
    explicit Exception(const char * message)
    {
        strncpy_s(msg_, message, 1024);
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

struct CodebookState
{
    uint16_t     codeword ;
    EfmplusState next_state ;
};


using CodebookArray = CodebookState[CDB_LINES][NUM_EFMPLUS_STATES];
//using CodebookArray = std::array<std::array<CodebookState, NUM_EFMPLUS_STATES>, CDB_LINES>;

template <size_t N>
struct Stats
{
    blaze::StaticVector<double, N, blaze::columnVector> prob;
    blaze::StaticVector<double, N, blaze::columnVector> var;
    void clear() { reset(prob); reset(var); }
};

class TransitionTable
{
public:
    TransitionTable();
    blaze::StaticMatrix<uint32_t, TRANSITION_TABLE_STATES, TRANSITION_TABLE_STATES> matrix;

    Stats<RDS_MAX + 1>             rds_stats;
    Stats<NUM_EFMPLUS_STATES>      efmplus_state_stats;
    Stats<NUM_SIGNS>               sign_bit_stats;
    Stats<TRANSITION_TABLE_STATES> state_stats;

    size_t size() { return TRANSITION_TABLE_STATES; }
    void solve(double eps);
    double get_codebook_variance() { return variance; }
private:
    blaze::StaticVector<double, TRANSITION_TABLE_STATES, blaze::columnVector> p; // this will end up being the same as probability
    void calculate_independent_probabilities();
    double variance;
    bool sanity_check();
};


class Codebook
{
public:
    Codebook(CodebookArray &codebook, std::shared_ptr<CodewordMap> codeword_map);

    void print_codebook_as_header(const char * name, int flag) const;
    bool calculate_pdf_with_transition_table();
    void optimise();
    void apply_strategy(uint32_t strategy, const std::vector<uint32_t>& indices);
    CodebookArray & line;

private:
    SUPER_STATE get_best_between_main_and_alternate(int byte, const SUPER_STATE & st);
    void get_best_between_states_0_3(int byte, const SUPER_STATE & st, SUPER_STATE & min_st);

    bool states_have_valid_codeword() const;
    bool is_consistent() const;
    bool is_valid() const;

    void print_line(FILE * f, int line_index) const;
    void print_transition_table();

    double collect_codebook_stats();
    bool initialise_transition_table();

    SUPER_STATE	choose_next_super_state(int byte, const SUPER_STATE & st, int * switch_cnt);

    TransitionTable transition;

    std::shared_ptr<CodewordMap> codeword_map;
    std::vector<double> rds_probability;

    std::bitset<15> valid_strategies_for_index(uint32_t index) const;
    bool can_be_applied(size_t strategy_, uint32_t index) const;
};



#endif /* CODEBOOK_H_ */
