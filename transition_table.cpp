//transition_table.cpp

#include <array>
#include <cassert>
#include <numeric>

#include "codebook.h"
#include "CodewordMap.h"


//--------------------------------------------------------------------------

#define BIT_VARIANCE 1
#define CLIP_RDS(x) ((abs(x) > (RDS_MAX)) ? (RDS_MAX) : (x))
#define BIT_TO_SIGN(x) (((x) == (0)) ? (1) : (-1));

static uint32_t SIGN_TO_BIT(int x)
{
	if (x >= 0)
		return 0;
	return 1;
}

static uint32_t super_state_with_rds(uint32_t rds, uint32_t state8, uint32_t sign_bit)
{
	return (rds * NUM_OF_STATE_NTZ_COMBINATIONS + state8) * NUM_SIGNS + sign_bit;
}

//--------------------------------------------------------------------------

bool switch_is_allowed(int leading, int ntz, int d, int k)
{
	return (leading + ntz >= d && leading + ntz <= k);
}

//--------------------------------------------------------------------------

bool new_wins_tie_breaker(const SUPER_STATE & new_st, const SUPER_STATE & min_st, const Stats<TRANSITION_TABLE_STATES> & full_state_stats)
{
	int new_sign = (new_st.rds >= 0) ? new_st.sign : -new_st.sign;
	int min_sign = (min_st.rds >= 0) ? min_st.sign : -min_st.sign;

	uint32_t new_full_state = (CLIP_RDS(abs(new_st.rds)) * NUM_OF_STATE_NTZ_COMBINATIONS + (old2new[new_st.state][new_st.ntz])) * NUM_SIGNS + SIGN_TO_BIT(new_sign);
	uint32_t min_full_state = (CLIP_RDS(abs(min_st.rds)) * NUM_OF_STATE_NTZ_COMBINATIONS + (old2new[min_st.state][min_st.ntz])) * NUM_SIGNS + SIGN_TO_BIT(min_sign);

	return  full_state_stats.var[new_full_state] < full_state_stats.var[min_full_state];

	return (full_state_stats.var[new_full_state]/full_state_stats.prob[new_full_state]
	        < full_state_stats.var[min_full_state]/full_state_stats.prob[min_full_state]);

}

//--------------------------------------------------------------------------

void calculate_new_rds_and_sum_var(int & new_rds_, uint32_t & new_sum_var_actual, uint32_t & new_sum_var_prediction,  const SequenceDetails & new_sequence, int rds, int8_t sign)
{
	int32_t new_rds = rds + sign * new_sequence.rds;
#ifdef BIT_VARIANCE
		new_sum_var_actual =  new_sequence.sum_of_nrzi_squares + 2 * new_sequence.sum_per_bit_rds * rds * SIGN_TO_BIT(sign) + SIZE_OF_ENCODED_OUTPUT_IN_BITS * rds * rds;

		new_sum_var_prediction =  new_sequence.sum_of_nrzi_squares + 2 * new_sequence.sum_per_bit_rds * rds * SIGN_TO_BIT(sign) + SIZE_OF_ENCODED_OUTPUT_IN_BITS * new_rds * new_rds;
#else
	//int new_sum_var = SIZE_OF_ENCODED_OUTPUT_IN_BITS * new_rds * new_rds;
#endif

	new_rds_ = new_rds;
}

//--------------------------------------------------------------------------

SUPER_STATE Codebook::get_best_between_main_and_alternate(int byte, const SUPER_STATE & st)
{
	//fprintf(stderr, "get_best_between_main_and_alternate for byte %d\n",byte);
	SUPER_STATE min_st(INVALID_SEQUENCE);
	for (uint32_t table_offset = 0; table_offset < CDB_LINES; table_offset += NUM_INDICES)
	{
		uint16_t codeword = line[byte + table_offset][st.state].codeword;

		if (codeword != INVALID_SEQUENCE)
		{
			SequenceDetails & p = (*codeword_map)[codeword];

			SUPER_STATE new_st(codeword);
			calculate_new_rds_and_sum_var(new_st.rds, new_st.sum_var_actual, new_st.sum_var_prediction, p, st.rds, st.sign);
			//fprintf(stderr, "calculated new rds and sum var\n");
			new_st.sign = p.has_odd_num_of_ones ? -st.sign : st.sign;
			new_st.ntz = p.trailing;
			new_st.state = line[byte + table_offset][st.state].next_state;
			//fprintf(stderr, "checking for draw\n");
			if ( (new_st.sum_var_prediction < min_st.sum_var_prediction) ||
			     ( (new_st.sum_var_prediction == min_st.sum_var_prediction) && new_wins_tie_breaker(new_st, min_st, transition.state_stats) ) )
			{
				min_st = new_st;
			}
		}
	}
	return std::move(min_st);
}

//--------------------------------------------------------------------------

void Codebook::get_best_between_states_0_3(int byte, const SUPER_STATE & st, SUPER_STATE & min_st)
{
	int d = codeword_map->get_d();
	int k = codeword_map->get_k();
	int switch_flag = 0;

	for (uint32_t table_offset = 0; table_offset < CDB_LINES; table_offset += NUM_INDICES)
	{
		uint16_t codeword = line[byte + table_offset][STATE_3 - st.state].codeword;

		if (codeword != INVALID_SEQUENCE)
		{
			SequenceDetails & p = (*codeword_map)[codeword];

			if (switch_is_allowed(p.leading, st.ntz, d, k))
			{
				SUPER_STATE new_st(codeword);

				calculate_new_rds_and_sum_var(new_st.rds, new_st.sum_var_actual, new_st.sum_var_prediction, p, st.rds, st.sign);

				new_st.sign = p.has_odd_num_of_ones ? -st.sign : st.sign;
				new_st.ntz = p.trailing;
				new_st.state = line[byte + table_offset][STATE_3 - st.state].next_state;

				if ( (new_st.sum_var_prediction < min_st.sum_var_prediction) ||
								 ( (new_st.sum_var_prediction == min_st.sum_var_prediction) && new_wins_tie_breaker(new_st, min_st, transition.state_stats) ) )
				{
					min_st = new_st;
					min_st.state = line[byte + table_offset][STATE_3 - st.state].next_state;
					min_st.codeword = codeword;
				}
				switch_flag = 1;
			}
		}
	}
}

//--------------------------------------------------------------------------

SUPER_STATE	Codebook::choose_next_super_state(int byte, const SUPER_STATE & st, int * switch_cnt)
{
	STATE_0_3_SWITCH state_0_3_switching = STATE_0_3_COMPLETE_SWITCH;

	SUPER_STATE min_st = get_best_between_main_and_alternate(byte, st);
	
	if (state_0_3_switching == STATE_0_3_COMPLETE_SWITCH || (state_0_3_switching == STATE_0_3_HALF_SWITCH && (byte & 0xFF) >= (NUM_INDICES-1)))
	{
		if (st.state == STATE_0 || st.state == STATE_3)
		{
			get_best_between_states_0_3(byte, st, min_st);

			//switch_cnt[st.state]++;
		}
	}
	return std::move(min_st);
}

//--------------------------------------------------------------------------

bool Codebook::initialise_transition_table()
{
	fprintf(stderr, "initialise_transition_table\n");
	int state_cnt[4];
	int switch_cnt[4] = { 0 };

	for (uint32_t rds = 0; rds < RDS_MAX + 1; rds+=2)
	{
		for (uint32_t state8 = 0 ; state8 < NUM_OF_STATE_NTZ_COMBINATIONS; state8++)
		{
			for (uint32_t sign_bit = 0; sign_bit < NUM_SIGNS; sign_bit++)
			{
				int8_t sign = BIT_TO_SIGN(sign_bit);
				uint32_t start_from_state = super_state_with_rds(rds, state8, sign_bit);

				if (start_from_state < 0 || start_from_state >= transition.size())
				{
					fprintf(stderr, "PANIC\n");
					return false;
				}
				// state should be in [0,1039], when rds_max=64
				double local_sum = 0.0;
				for (uint16_t input = 0; input < NUM_INDICES; input++)
				{
					SUPER_STATE current_st(INVALID_SEQUENCE);
					current_st.rds = rds;
					current_st.ntz = new2old[state8].ntz;
					current_st.state = new2old[state8].state;
					current_st.sign = sign;

					SUPER_STATE new_st = choose_next_super_state(input, current_st, switch_cnt);
					local_sum += new_st.sum_var_actual;

					//fprintf(stderr, "rds=%d, state8=%d, sign=%d, input=%d, new_sum_var=%d, min_rds=%d next_state=%d, 0x%04x\n",
					//	rds, state8, sign, input, new_st.sum_var_actual, new_st.rds, new_st.state, new_st.codeword);
					state_cnt[current_st.state]++;
					if (new_st.rds < 0)
						new_st.sign = -new_st.sign;

					uint32_t new_state = (CLIP_RDS(abs(new_st.rds)) * NUM_OF_STATE_NTZ_COMBINATIONS + old2new[new_st.state][new_st.ntz]) * NUM_SIGNS + SIGN_TO_BIT(new_st.sign);
					transition.matrix(new_state, start_from_state)++;
				}

				transition.state_stats.var[start_from_state] = local_sum / double(NUM_INDICES * SIZE_OF_ENCODED_OUTPUT_IN_BITS) ;


			}
		}
	}

	return true;
}

//--------------------------------------------------------------------------

void Codebook::print_transition_table()
{
	for (size_t i = 0UL; i < transition.matrix.rows(); ++i)
	{
		for (size_t j = 0UL; j<transition.matrix.columns(); ++j)
		{
			fprintf(stderr, "%2u ", transition.matrix(i,j));
		}
		fprintf(stderr, "\n");
	}
}

//--------------------------------------------------------------------------

bool Codebook::calculate_pdf_with_transition_table()
{
   /* Calculation of the super-state (rds,ntz,sign,state) probabilities using transition table*/
   /* Calculation of the variance that the codebook would give if it was applied to a uniform file */
   /* where each 8-bits have the same probability (1/256) */

	// New attempt - 2 iterations
	// in the second attempt the full_state_variances have been initialised and are used as tie breakers
	for (int pdf_iter = 0; pdf_iter < 2; pdf_iter++)
	{
		double bit_var = collect_codebook_stats();
		fprintf(stderr, "!!!!!!Estimated Bit variance: %lf\n",bit_var);
	}
	return true;
}

//--------------------------------------------------------------------------

double Codebook::collect_codebook_stats()
{
	fprintf(stderr, "collect_codebook_stats\n");
	initialise_transition_table();
	transition.solve(1e-14);
	double variance = transition.get_codebook_variance();

	return variance;
}

//--------------------------------------------------------------------------

TransitionTable::TransitionTable()
{
	printf("constructor called\n");
	for (uint32_t i = 0; i <= RDS_MAX; i++)
	{
		for (uint32_t j = 0; j < NUM_OF_SUPER_STATES; j++)
		{
			state_stats.prob[i * NUM_OF_SUPER_STATES + j] = int(i == 0) / double(NUM_OF_SUPER_STATES);
		}
	}
}

bool TransitionTable::sanity_check()
{
	double sum = 0.0;
	for (size_t i = 0; i < state_stats.prob.size(); i++)
	{
		sum += state_stats.prob[i];
	}
	return (sum - 1 < 1e-14);
}

template<typename T>
void normalise_vector(T & v)
{
	double sum = 0;
	for (size_t i = 0; i < v.size(); i++)
	{
		sum += v[i];
	}

	if (sum != 0)
	{
		for (size_t i = 0; i < v.size(); i++)
		{
			v[i] /= sum;
		}
	}
}

void TransitionTable::calculate_independent_probabilities()
{
	sign_bit_stats.clear();
	rds_stats.clear();
	efmplus_state_stats.clear();

	// calculate probability of rds={0,...,RDS_MAX), state8={0, ...8}, sign = {-1,1}
	for (uint32_t rds = 0; rds < RDS_MAX + 1; rds += 2)
	{
		for (uint32_t state8 = 0; state8 < NUM_OF_STATE_NTZ_COMBINATIONS; state8++)
		{
			double prob_sum = 0;
			double var_sum = 0;
			for (uint32_t sign_bit = 0; sign_bit < NUM_SIGNS; sign_bit++)
			{
				uint32_t state = super_state_with_rds(rds, state8, sign_bit);
				auto & probability = state_stats.prob;
				prob_sum += probability[state];
				sign_bit_stats.prob[sign_bit] += probability[state];

				var_sum += state_stats.var[state] * probability[state];
				sign_bit_stats.var[sign_bit] += state_stats.var[state] * probability[state];
			}
			auto efmplus_state = new2old[state8].state;
			efmplus_state_stats.prob[efmplus_state] += prob_sum;
			efmplus_state_stats.var[efmplus_state] += var_sum;
			rds_stats.prob[rds] += prob_sum;
			rds_stats.var[rds] += var_sum;
		}
	}

	normalise_vector(rds_stats.prob);
	normalise_vector(efmplus_state_stats.prob);
	normalise_vector(sign_bit_stats.prob);

	double sum = 0;
	for (uint32_t i = 0; i < efmplus_state_stats.prob.size(); i++)
	{
		sum += efmplus_state_stats.var[i];
		fprintf(stderr, "!!!!state %d, prob %lf, var %lf\n", i, efmplus_state_stats.prob[i], efmplus_state_stats.var[i]);
	}
	fprintf(stderr, "estimated variance from efmplus_state %lf\n", sum);

	sum = 0;
	for (uint32_t i = 0; i < rds_stats.prob.size(); i++)
	{
		sum += rds_stats.var[i];
		//fprintf(stderr, "!!!!state %d, prob %lf, var %lf\n", i, rds_stats.prob[i], rds_stats.var[i]);
	}
	fprintf(stderr, "estimated variance from rds %lf\n", sum);

	sum = 0;
	for (uint32_t i = 0; i < sign_bit_stats.prob.size(); i++)
	{
		sum += sign_bit_stats.var[i];
		fprintf(stderr, "!!!!state %d, prob %lf, var %lf\n", i, sign_bit_stats.prob[i], sign_bit_stats.var[i]);
	}
	fprintf(stderr, "estimated variance from sign %lf\n", sum);

	variance = sum;
}

void TransitionTable::solve(double eps)
{
	fprintf(stderr, "solve\n");

	for (uint32_t i = 0; i < p.size(); i++)
	{
		p[i] = state_stats.prob[i];
	}

	double sum, sum2, temp, dist;
	auto & ptr1 = p;
	auto & ptr2 = state_stats.prob;

	dist = sum2 = 1.0;
	for (uint32_t iter = 100; dist / sum2 >= eps && iter > 0; iter--)
	{
		ptr2 = matrix * ptr1;
		sum = 0.0;
		for (size_t i = 0; i < state_stats.prob.size(); i += 2 * NUM_OF_SUPER_STATES)
		{
			for (uint32_t j = i; j < i + NUM_OF_SUPER_STATES; j++)
				sum += ptr2[j];
		}
		sum2 = 0.0;
		dist = 0.0;
		for (size_t i = 0; i < state_stats.prob.size(); i += 2 * NUM_OF_SUPER_STATES)
		{
			for (uint32_t j = i; j < i + NUM_OF_SUPER_STATES; j++)
			{
				ptr2[j] /= sum;
				sum2 += ptr2[j] * ptr2[j];
				temp = ptr2[j] - ptr1[j];
				dist += temp*temp;
			}
		}
		std::swap(ptr1, ptr2);
	}
	ptr1 = ptr2;

	if (!sanity_check())
		fprintf(stderr, "error in sanity check\n");

	calculate_independent_probabilities();
}

