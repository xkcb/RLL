/*
 * constants.h
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_


template <typename Underlying, int UnitType, int low, int high>
class Unit
{
public:
    Unit(int i) : val(i)  //constructor only accepts args of type Unit and not any integer
	{
		if (i < low || i > high)
		{
			//std::cout << "value " << i << "low " << low << " high " << high << std::endl;
			throw std::invalid_argument("");
		}
	}
	operator int() const {return val;}

private:
	Underlying val;
};

enum {EfmplusStateType};

typedef enum
{
	NO_STATE = -1,
	STATE_0 = 0,
	STATE_1,
	STATE_2,
	STATE_3,
	NUM_EFMPLUS_STATES

} EFMPLUS_STATE;
using EfmplusState = Unit<int, EfmplusStateType, -1, 3>;

const uint32_t NUM_SIGNS = 2;
const uint32_t NUM_INDICES = 256;
const uint32_t CDB_LINE_ELEMENTS = NUM_EFMPLUS_STATES * 2;
const uint32_t CDB_LINES = NUM_INDICES * 2;
const uint32_t RDS_MAX = 32;
const uint32_t INVALID_RDS = 1000;
const uint32_t INVALID_SEQUENCE = 0XFFFF;
const uint32_t INVALID_SUM_VAR = 1000000;
const uint32_t NUM_OF_STATE_NTZ_COMBINATIONS = 8;
// NUM_OF_SUPER_STATES is 8(S0_0, S0_1, S1, S2, S3_6, S3_7, S3_8, S3_9)x 2 signs
const uint32_t NUM_OF_SUPER_STATES = NUM_OF_STATE_NTZ_COMBINATIONS * NUM_SIGNS;
const uint32_t TRANSITION_TABLE_STATES = (RDS_MAX + 1) * NUM_OF_SUPER_STATES;	
// every state of the transition table is a combination of (rds,ntz,sign,efmplus_state) namely (rds,0..9,+,state) (rds,0..9,-,state)
//+1 for the edge state (|rds|>rds_max)

const uint32_t SIZE_OF_INPUT_IN_BITS = 8;
const uint32_t SIZE_OF_ENCODED_OUTPUT_IN_BITS = 16;
struct STATE_NTZ_COMBINATION
{
	EfmplusState state;
	uint32_t     ntz;
} ;

const STATE_NTZ_COMBINATION new2old[NUM_OF_STATE_NTZ_COMBINATIONS] =
		{ {STATE_0, 0},
		  {STATE_0, 1},
		  {STATE_1, 4},
		  {STATE_2, 4},
		  {STATE_3, 6},
		  {STATE_3, 7},
		  {STATE_3, 8},
		  {STATE_3, 9}
		}; //combination of efmplus_state and ntz

const int old2new[4][10] = {
	{ 0, 1,-1,-1,-1,-1,-1,-1,-1,-1},
	{ 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
	{ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
	{-1,-1,-1,-1,-1,-1, 4, 5, 6, 7}
};

typedef enum
{
	STATE_0_3_NO_SWITCH = 0,
	STATE_0_3_HALF_SWITCH,
	STATE_0_3_COMPLETE_SWITCH
} STATE_0_3_SWITCH;


#endif /* CONSTANTS_H_ */
