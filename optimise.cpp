/*
* optimisation.cpp
*/

#include "codebook.h"
#include "constants.h"
#include <bitset>
#include <vector>

/*
void Codebook::optimise()
{
    std::array<uint32_t, 2> indices;
    std::array<std::bitset<16>, CDB_LINES> valid_strategies = calculate_valid_strategies();
    for (indices[0] = 0; indices[0] < 348 - 1; indices[0]++)
    {
        for (indices[1] = indices[0] + 1; indices[1] < 348; indices[1]++)
        {
            //fprintf(stderr, "index1 %u\n", indices[1]);
            std::bitset<16> common_strategies = valid_strategies[indices[0]] & valid_strategies[indices[1]];

            auto filtered_strategies = filter_pair(indices, common_strategies);
            for (size_t i = filtered_strategies.size() - 1; i > 0; i--)
            {
                if (filtered_strategies[i])
                {
                    if (improve_cdb(std::bitset<4>(i), indices))
                    {
                        valid_strategies[indices[0]] = valid_strategies_for_index(indices[0]);
                        valid_strategies[indices[1]] = valid_strategies_for_index(indices[1]);
                    }
                }
            }
        }
        std::string name = "optimised_" + std::to_string(indices[0]) + std::to_string(min_variance)+".h";
        print_codebook_as_header(name.c_str(), 0);

    }
}
*/
void Codebook::optimise()
{
    std::array<uint32_t, 2> indices;
    std::array<std::bitset<16>, CDB_LINES> valid_strategies = calculate_valid_strategies();

    for (auto strategy = 15; strategy > 0; strategy--)
    {
        for (indices[0] = 0; indices[0] < 348; indices[0]++)
        {
            for (indices[1] = indices[0] + 1; indices[1] < 348; indices[1]++)
            {
                std::bitset<16> common_strategies = valid_strategies[indices[0]] & valid_strategies[indices[1]];

                auto filtered_strategies = filter_pair(indices, common_strategies);
                if (filtered_strategies[strategy])
                {
                    //fprintf(stderr, "trying strategy %u between %u,%u\n", strategy, indices[0], indices[1]);
                    if (improve_cdb(std::bitset<4>(strategy), indices))
                    {
                        valid_strategies[indices[0]] = valid_strategies_for_index(indices[0]);
                        valid_strategies[indices[1]] = valid_strategies_for_index(indices[1]);
                    }
                }
           }
        }
        std::string name = "optimised_" + std::to_string(indices[0]) + std::to_string(min_variance) + ".h";
        print_codebook_as_header(name.c_str(), 0);
    }
}

void Codebook::apply_strategy(std::bitset<4> strategy, const std::array<uint32_t, 2> & indices)
{
    for (auto state = 0; state < NUM_EFMPLUS_STATES; state++)
    {
        if (strategy[(NUM_EFMPLUS_STATES - 1) - state])
            std::swap(line[indices[0]][state], line[indices[1]][state]);
    }
}

std::bitset<16> Codebook::filter_pair(const std::array<uint32_t, 2> & indices, const std::bitset<16> & common_strategies)
{
    std::bitset<16> filtered_pairs(0);
    for (size_t strategy = 1; strategy < common_strategies.size(); strategy++)
    {
        if (common_strategies[strategy])
        {
            apply_strategy(std::bitset<4>(strategy), indices);
            if (states_have_valid_codeword(indices))
                filtered_pairs.set(strategy);
            apply_strategy(std::bitset<4>(strategy), indices);
        }
    }
    return filtered_pairs;
}


bool Codebook::improve_cdb(std::bitset<4>strategy, std::array<uint32_t, 2> & indices)
{
    if (min_variance == 0)
        min_variance = transition.get_codebook_variance();

    apply_strategy(strategy, indices);
    double bit_var = 0;
    for (int pdf_iter = 0; pdf_iter < 1; pdf_iter++)
    {
        bit_var = collect_codebook_stats();
        //fprintf(stderr, "!!!!!!Estimated Bit variance: %lf, min %lf\n", bit_var, min_variance);
    }
    if (bit_var < min_variance && (( min_variance / bit_var) > 1.0001))
    {
        fprintf(stderr, "%lf keep strategy %u between indices, improvement %.3lf %u,%u\n", bit_var, strategy, min_variance/ bit_var, indices[0], indices[1]);
        min_variance = bit_var;
        return true;
    }
    else
    {
        apply_strategy(strategy, indices);
    }

    return false;
}

bool Codebook::can_be_applied(size_t strategy_, uint32_t index) const
{
	std::bitset<4> strategy(strategy_);

	for (auto state = 0; state < NUM_EFMPLUS_STATES; state++)
	{
		if (strategy[(NUM_EFMPLUS_STATES - 1) - state]) // strategy includes state
		{
			auto next_state = line[index][state].next_state;
			if (next_state != NO_STATE)
			{
				for (auto state2 = 0; state2 < NUM_EFMPLUS_STATES; state2++)
				{
					if ((!strategy[(NUM_EFMPLUS_STATES - 1) - state2]) && line[index][state2].next_state == next_state)
						return false;
				}
			}
		}
	}
	return true;
}

std::bitset<16> Codebook::valid_strategies_for_index(uint32_t index) const
{
	std::bitset<16> valid_strategies(0);

	for (size_t strategy = 1; strategy <= 15; strategy++)
	{
		if (can_be_applied(strategy, index))
		{
			valid_strategies.set(strategy);
		}
	}
	return valid_strategies;
}

std::array<std::bitset<16>, CDB_LINES> Codebook::calculate_valid_strategies()
{
    std::array<std::bitset<16>, CDB_LINES> valid_strategies;

    for (uint32_t i = 0; i < CDB_LINES; i++)
    {
        valid_strategies[i] = valid_strategies_for_index(i);
    }
    return valid_strategies;
}
