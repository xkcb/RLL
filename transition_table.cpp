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

static double get_total_codebook_variance(const STATS<double> & stats);
static void normalise_variance(STATS<double> & stats);

//--------------------------------------------------------------------------


template<typename T, size_t N>
void matrix_vector_multiplication(const std::array< std::array<T, N>, N> & A, const std::vector<double> & x, std::vector<double> & dest)
{
	assert(N == x.size());
	for (size_t i = 0; i < A.size(); i++)
	{
		auto row = A[i];
		dest[i] = std::inner_product(row.begin(), row.end(), x.begin(), 0.0);
	}
}

//--------------------------------------------------------------------------

template<typename T, size_t N>
void solve_iterative(std::vector<double> & p, std::vector<double> & b, const std::array< std::array<T, N>, N> & transition, double eps)
{
	fprintf(stderr, "solve iterative\n");

	double sum, sum2, temp, dist;
	std::vector<double> & ptr1 = p;
	std::vector<double> & ptr2 = b;
	std::vector<double> & ptr3 = p;
	
	dist = sum2 = 1.0;  
	for (uint32_t iter = 100; dist/sum2 >= eps && iter > 0; iter--)
	{
		ptr3 = ptr1;
		ptr1 = ptr2;
		ptr2 = ptr3;
		
		matrix_vector_multiplication(transition, ptr1, ptr2);
		sum = 0.0;
		for (size_t i = 0; i < b.size(); i += 2 * NUM_OF_SUPER_STATES)
		{
			for (uint32_t j = i; j < i + NUM_OF_SUPER_STATES; j++)
				sum += ptr2[j];
		}
		sum2 = 0.0;
		dist = 0.0;
		for (size_t i = 0; i < b.size(); i += 2 * NUM_OF_SUPER_STATES)
		{
			for (uint32_t j = i; j < i + NUM_OF_SUPER_STATES; j++)
			{
				ptr2[j] /= sum;
				sum2 += ptr2[j] * ptr2[j];
				temp = ptr2[j] - ptr1[j];
				dist += temp*temp;
			}
		}
	}
	ptr1 = ptr2;
}

//--------------------------------------------------------------------------

bool switch_is_allowed(int leading, int ntz, int d, int k)
{
	return (leading + ntz >= d && leading + ntz <= k);
}

//--------------------------------------------------------------------------

bool new_wins_tie_breaker(const SUPER_STATE & new_st, const SUPER_STATE & min_st, const std::shared_ptr<STATS<double>> full_state_stats) 
{
	int new_sign = (new_st.rds >= 0) ? new_st.sign : -new_st.sign;
	int min_sign = (min_st.rds >= 0) ? min_st.sign : -min_st.sign;

	uint32_t new_full_state = (CLIP_RDS(abs(new_st.rds)) * NUM_OF_STATE_NTZ_COMBINATIONS + (old2new[new_st.state][new_st.ntz])) * NUM_SIGNS + SIGN_TO_BIT(new_sign);
	uint32_t min_full_state = (CLIP_RDS(abs(min_st.rds)) * NUM_OF_STATE_NTZ_COMBINATIONS + (old2new[min_st.state][min_st.ntz])) * NUM_SIGNS + SIGN_TO_BIT(min_sign);

	return  full_state_stats->variance[new_full_state] < full_state_stats->variance[min_full_state];

	return (full_state_stats->variance[new_full_state]/full_state_stats->probability[new_full_state]
	        < full_state_stats->variance[min_full_state]/full_state_stats->probability[min_full_state]);

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
			     ( (new_st.sum_var_prediction == min_st.sum_var_prediction) && new_wins_tie_breaker(new_st, min_st, this->global_full_state_stats) ) )
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
								 ( (new_st.sum_var_prediction == min_st.sum_var_prediction) && new_wins_tie_breaker(new_st, min_st, this->global_full_state_stats) ) )
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
				for (uint16_t input = 0; input < NUM_INDICES; input++)
				{
					SUPER_STATE current_st(INVALID_SEQUENCE);
					current_st.rds = rds;
					current_st.ntz = new2old[state8].ntz;
					current_st.state = new2old[state8].state;
					current_st.sign = sign;

					SUPER_STATE new_st = choose_next_super_state(input, current_st, switch_cnt);
					//fprintf(stderr, "rds=%d, state8=%d, sign=%d, input=%d, new_sum_var=%d, min_rds=%d next_state=%d, 0x%04x\n",
					//	rds, state8, sign, input, new_st.sum_var_actual, new_st.rds, new_st.state, new_st.codeword);
					state_cnt[current_st.state]++;
					if (new_st.rds < 0)
						new_st.sign = -new_st.sign;

					uint32_t new_state = (CLIP_RDS(abs(new_st.rds)) * NUM_OF_STATE_NTZ_COMBINATIONS + old2new[new_st.state][new_st.ntz]) * NUM_SIGNS + SIGN_TO_BIT(new_st.sign);
					transition[new_state][start_from_state]++;
				}
			}
		}
	}
	return true;
}

//--------------------------------------------------------------------------

void Codebook::print_transition_table()
{
	for (auto &row : transition)
	{
		for (auto &elem : row)
		{
			fprintf(stderr, "%2u ",elem);
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
	global_full_state_stats = std::make_shared<STATS<double>>(TRANSITION_TABLE_STATES);

	// New attempt - 2 iterations
	// in the second attempt the full_state_variances have been initialised and are used as tie breakers
	for (int pdf_iter = 0; pdf_iter < 2; pdf_iter++)
	{
		double bit_var = collect_codebook_stats(global_full_state_stats->probability);
		fprintf(stderr, "!!!!!!Estimated Bit variance: %lf\n",bit_var);

		//codeword_probability(RDS_MAX, p, d, k, state_stats.probability);
	}
	return true;
}

//--------------------------------------------------------------------------

void Codebook::get_probability_vector(std::vector<double> &p)
{
	initialise_transition_table();

	/* p,b - probability vectors*/
	std::vector<double> b(TRANSITION_TABLE_STATES);
	fprintf(stderr, "create vector b (%d)\n", TRANSITION_TABLE_STATES);

	for (uint32_t i = 0; i <= RDS_MAX; i++)
	{
		for (uint32_t j = 0; j < NUM_OF_SUPER_STATES; j++)
		{
			b[i * NUM_OF_SUPER_STATES + j] = int(i == 0) / double(NUM_OF_SUPER_STATES);
		}
	}
	solve_iterative<uint32_t, TRANSITION_TABLE_STATES>(p, b, transition, 1e-14);
}

double Codebook::collect_codebook_stats(std::vector<double> &p)
{
	fprintf(stderr, "collect_codebook_stats\n");

	get_probability_vector(p);

	//check that the probabilities sum to 1
	double prob_sum = std::accumulate(p.begin(), p.end(), 0.0);
	std::cout << "SUM " << prob_sum << " should be 1" << std::endl;

	double frequency[RDS_MAX + 1][NUM_SIGNS][NUM_EFMPLUS_STATES];
	STATS <double> state_stats(NUM_EFMPLUS_STATES);
	STATS_2D state_sign_stats(NUM_EFMPLUS_STATES, NUM_SIGNS);
	STATS <double> sign_stats(NUM_SIGNS);
	STATS_2D index_stats(NUM_INDICES, TRANSITION_TABLE_STATES);
	STATS <double> rds_stats(RDS_MAX + 1);
	int switch_cnt[4] = { 0 }; //0: 1->4,


	for (uint32_t rds = 0; rds < RDS_MAX + 1; rds++)
		for (int8_t sign_bit = 0; sign_bit < NUM_SIGNS; sign_bit++)
			for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
				frequency[rds][sign_bit][state] = 0.0;

	for (uint32_t i = 0;  i < RDS_MAX + 1; i++)
	{
		double sum = 0.0;
		for (uint32_t j = 0; j < NUM_OF_SUPER_STATES; j++)
		{
			sum += p[(i) * NUM_OF_SUPER_STATES + j];
		}
		rds_stats.probability[i] = (i == 0) ? sum : sum / 2.0;
	}

	std::vector<uint32_t> rds_values(RDS_MAX + 1);

	for (uint32_t rds = 0; rds < rds_values.size(); rds++)
		rds_values[rds] = rds;


	


	// update frequency tables

 	for (uint32_t rds = 0; rds <= RDS_MAX; rds+=2)
	{
		for (uint32_t state8 = 0 ; state8 < NUM_OF_STATE_NTZ_COMBINATIONS; state8++)
		{
			for (uint32_t sign_bit = 0 ; sign_bit < NUM_SIGNS; sign_bit++)
			{
				int8_t sign = BIT_TO_SIGN(sign_bit);
				EfmplusState efmplus_state = new2old[state8].state;
				uint32_t state = super_state_with_rds(rds, state8, sign_bit);
				frequency[rds][sign_bit][efmplus_state] += p[state];
			
				double local_sum = 0.0;

				for (uint16_t input = 0; input < NUM_INDICES; input++)
				{
					SUPER_STATE st(INVALID_SEQUENCE);
					st.rds = rds;
					st.ntz = new2old[state8].ntz;
					st.state = efmplus_state;
					st.sign = sign;

					SUPER_STATE new_st = choose_next_super_state(input, st, switch_cnt);
					index_stats.variance[input][state] = new_st.sum_var_actual;
					local_sum += new_st.sum_var_actual;
					//fprintf(stderr, "rds=%u, state8=%u, sign=%d, input=%u, new_sum_var=%u, min_rds=%d, next_state=%u local_sum %d, codeword 0x%04x\n",
						//	rds, state8, sign, input, new_st.sum_var_actual, new_st.rds, new_st.state, local_sum, new_st.codeword);
				}
				global_full_state_stats->variance[state] = local_sum * global_full_state_stats->probability[state];
				//fprintf(stderr, "local_sum = %lf\n", local_sum);
				state_stats.variance[efmplus_state] += local_sum * p[state];
				sign_stats.variance[sign_bit] += local_sum * p[state];
				state_sign_stats.variance[efmplus_state][sign_bit] += local_sum * p[state];
			}
		}
	}

 	normalise_variance(state_stats);
 	normalise_variance(sign_stats);
 	//normalise_variance(state_sign_stats);
 	normalise_variance(*global_full_state_stats.get());

	for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
	{
		double sum = 0.0;
		for (uint32_t rds = 0; rds < RDS_MAX + 1; rds++)
		{
			for (int8_t sign = 0; sign < NUM_SIGNS; sign++)
			{
				sum += frequency[rds][sign][state];
			}
		}
		state_stats.probability[state] = sum;
	}

	for (int8_t sign = 0; sign < NUM_SIGNS; sign++)
	{
		double sum = 0.0;
		for (auto state = 0; state < NUM_EFMPLUS_STATES; state++)
		{
			for (uint32_t rds = 0; rds < RDS_MAX + 1; rds++)
			{
				sum+=frequency[rds][sign][state];
			}
		}
		sign_stats.probability[sign]=sum;
	}

	for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
	{
		for (int8_t sign = 0; sign < NUM_SIGNS; sign++)
		{
			double sum = 0.0;
			for (uint32_t rds = 0; rds <= RDS_MAX; rds++) {
				sum+=frequency[rds][sign][state];
			}
			state_sign_stats.probability[state][sign] =sum;
		}
	}


	double sum = get_total_codebook_variance(state_stats);
	fprintf(stderr, "VARIANCE ESTIMATION (from states) %f\n",sum);

	sum = get_total_codebook_variance(sign_stats);
	fprintf(stderr, "VARIANCE ESTIMATION (from sign) %f\n",sum);

	sum=0.0;
	for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
	{
		for (int8_t sign = 0; sign < NUM_SIGNS; sign++)
		{
			sum += state_sign_stats.variance[state][sign] * state_sign_stats.probability[state][sign];
		}
	}

	fprintf(stderr, "VARIANCE ESTIMATION (from state/sign) %f\n",sum);

	sum = get_total_codebook_variance(*global_full_state_stats.get());

	fprintf(stderr, "VARIANCE ESTIMATION (from full state) %f\n",sum);

	for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
	{
		if (state_stats.probability[state] != 0.0)
			state_stats.variance[state] /= state_stats.probability[state];
		fprintf(stderr, "State = %2d prob = %f variance = %f\n", state, state_stats.probability[state], state_stats.variance[state]);
	}

	for (int8_t sign = 0; sign < NUM_SIGNS; sign++)
	{
		if (sign_stats.probability[sign] != 0.0)
			sign_stats.variance[sign] /= sign_stats.probability[sign];
		fprintf(stderr, "Sign  = %2d prob = %f variance = %f\n", 1-2*sign, sign_stats.probability[sign], sign_stats.variance[sign]);
	}

	for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
	{
		for (int8_t sign = 0; sign < NUM_SIGNS; sign++)
		{
			if (state_sign_stats.probability[state][sign] != 0.0)
				state_sign_stats.variance[state][sign] /= state_sign_stats.probability[state][sign];
			fprintf(stderr, "State,sign = (%2d,%2d) prob = %f variance = %f\n",state+1, 1-2*sign, state_sign_stats.probability[state][sign],state_sign_stats.variance[state][sign]);
		}
	}

	return sum;
}

//--------------------------------------------------------------------------

double get_total_codebook_variance(const STATS<double> & stats)
{
	return std::accumulate(stats.variance.begin(), stats.variance.end(), 0.0);
}

//--------------------------------------------------------------------------

void normalise_variance(STATS<double> & stats)
{
	for (auto & var : stats.variance)
	{
		var /=  double(NUM_INDICES * SIZE_OF_ENCODED_OUTPUT_IN_BITS);
	}
}

//--------------------------------------------------------------------------
