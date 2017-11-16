// CodewordMap.h

#ifndef CODEWORDMAP_H_
#define CODEWORDMAP_H_

#include <list>
#include <map>
#include "constants.h"

//A (d,k) rll sequence is a binary sequence that satisfies the constraint
//that between two '1's there are at least 'd' and at most 'k' '0's

const int MAX_ONES_IN_A_SEQUENCE = 32; // as much as an integer can hold

struct SequenceDetails
{
	static const uint32_t length;
	uint32_t leading              = 0;
	uint32_t trailing             = 0;
	int32_t  rds                  = INVALID_RDS; //running digital sum at the end of the sequence
	bool     has_odd_num_of_ones  = false;
	uint32_t sum_of_nrzi_squares  = 0;
	int32_t  sum_per_bit_rds      = 0;

	SequenceDetails() = default;
	SequenceDetails(int rll_sequence);
};

typedef std::map<unsigned int, SequenceDetails> CODEWORD_MAP;

class CodewordMap
{
public:
	static const uint32_t rll_sequence_length = SIZE_OF_ENCODED_OUTPUT_IN_BITS; // rll sequence's length that satisfies the (d,k) constraint

	CodewordMap(int d, int k);
	CodewordMap();
	~CodewordMap();
	void print_codes();
	int get_d() const;
	int get_k() const;

    SequenceDetails & operator[](std::size_t idx);

private:
	void 			 create_codeword_map();
	void 			 create_sequences();

	std::list<int> * rll_coding(int length, int n, int left_restr, int right_restr, int start);
	void 			 merge(std::list<int> * root, std::list<int> & left, std::list<int> & right, unsigned int str);

	int min_constraint_d    = 2;
	int max_constraint_k 	= 10;
	int Amax 				= (rll_sequence_length + min_constraint_d) / (min_constraint_d + 1); //maximum number of '1's that a sequence may have
	int Amin 				= min_constraint_d / (max_constraint_k + 1); // minimum number of '1's that a sequence may have
	int rll_sequences_found = 0;
	std::list<int> * rll_list[MAX_ONES_IN_A_SEQUENCE];
	CODEWORD_MAP     codeword_map;
};


#endif /* CODEWORDMAP_H_ */
