// CodewordMap.cpp

#include "CodewordMap.h"
#include <bitset>

//--------------------------------------------------------------------------

static int get_ntz(int number, int number_length);
static int get_nlz(int number, int number_length);
static void get_rds_and_bit_var(int number, int number_length, int & rds, uint32_t & bit_var, int32_t & sum_per_bit_rds);


//--------------------------------------------------------------------------
const uint32_t SequenceDetails::length = CodewordMap::rll_sequence_length;

SequenceDetails::SequenceDetails(int rll_sequence)
{
    std::bitset<SequenceDetails::length> rll(rll_sequence);
    trailing = get_ntz(rll_sequence, SequenceDetails::length);
    leading = get_nlz(rll_sequence, SequenceDetails::length);
    has_odd_num_of_ones = rll.count() & 1;
    get_rds_and_bit_var(rll_sequence, SequenceDetails::length, rds, sum_of_nrzi_squares, sum_per_bit_rds);
    //fprintf(stderr, "sequence details: rds: %d, sum_var %u, sum %d\n", rds, sum_of_nrzi_squares, sum_per_bit_rds);
}

//--------------------------------------------------------------------------

CodewordMap::CodewordMap()
{
    create_sequences();
    create_codeword_map();
}

//--------------------------------------------------------------------------

CodewordMap::CodewordMap(int min_constraint_d_, int max_constraint_k_)
: min_constraint_d(min_constraint_d_),
  max_constraint_k(max_constraint_k_)
{

}

//--------------------------------------------------------------------------

CodewordMap::~CodewordMap()
{
    // TODO Auto-generated destructor stub
}

//--------------------------------------------------------------------------

int CodewordMap::get_d() const
{
    return min_constraint_d;
}

//--------------------------------------------------------------------------

int CodewordMap::get_k() const
{
    return max_constraint_k;
}

//--------------------------------------------------------------------------

SequenceDetails & CodewordMap::operator[](std::size_t idx)
{
    return codeword_map[idx];
};

//--------------------------------------------------------------------------

void CodewordMap::print_codes()
{
    FILE * file;

    fopen_s(&file, "codes.txt", "w");
    if (!file)
    {
        fprintf(stderr, "could not open file codes.txt for write");
        return ;
    }
    fprintf(file, "run-length limited codes (d, k, N) = (%d, %d, %d)\n",min_constraint_d, max_constraint_k, rll_sequence_length);

    for (int num_of_ones = Amin; num_of_ones <= Amax; num_of_ones++)
    {
        std::list<int> * list =  rll_list[num_of_ones - Amin];

        fprintf(file, "codes with %d '1'\n", num_of_ones);

        for (const auto & code : *list)
        {
            std::bitset<rll_sequence_length> codeword(code);
            fprintf(file, "0x%04x \t", code);

            fprintf(file, "%s ", codeword.to_string().c_str());

            fprintf(file, " %2d %2d rds: %2d, sum_var %d\n", codeword_map[code].leading, codeword_map[code].trailing, codeword_map[code].rds, codeword_map[code].sum_of_nrzi_squares);
        }
        fprintf(file, "(%lu sequences)\n", list->size());
        fprintf(file, "\n");
    }
    fprintf(file, "Total sequences %d\n", rll_sequences_found);
    fclose(file);
}

//--------------------------------------------------------------------------

void CodewordMap::create_codeword_map()
{
    for (int number_of_ones = Amin; number_of_ones <= Amax; number_of_ones++)
    {
        std::list<int> * list = rll_list[number_of_ones - Amin];
        for (const auto & val : *list)
        {
            codeword_map[val] = SequenceDetails{val};
        }
    }
    codeword_map[INVALID_SEQUENCE] = SequenceDetails{};
}

//--------------------------------------------------------------------------

void CodewordMap::create_sequences()
{
    rll_sequences_found = 0;

    for (int number_of_ones = Amin; number_of_ones <= Amax; number_of_ones++)
    {
        rll_list[number_of_ones - Amin] = rll_coding(rll_sequence_length, number_of_ones, 0, 0, 32 - rll_sequence_length);
        rll_sequences_found += rll_list[number_of_ones - Amin]->size();
    }
}

//--------------------------------------------------------------------------

std::list<int> * CodewordMap::rll_coding(int length, int number_of_ones, int left_restr, int right_restr, int start)
{
    unsigned int str;
    int jmin, jmax;
    int c;
    std::list<int> * root = new (std::list<int>);

    if (number_of_ones <= 0) return root;

    if (number_of_ones > 1)
    {
        int l = (number_of_ones - 1) >> 1;    // number of '1's that should fit left and right of the middle '1'
        int r = number_of_ones >> 1;
        jmin = start + (min_constraint_d + 1) * l + min_constraint_d * left_restr;
        if ((c = (start + length - 1) - (max_constraint_k + 1) * r - max_constraint_k) > jmin)
            jmin = c;
        jmax = start + (max_constraint_k + 1) * l + max_constraint_k;
        if ((c = (start + length - 1) - (min_constraint_d + 1) * r - min_constraint_d * right_restr) < jmax)
            jmax = c;
        str = 1 << (32 - jmin - 1);
        int l_left = jmin - start;					// available length left and right
        int l_right = (start + length - 1) - jmin;

        for (int i = jmin; i <= jmax; i++)
        {
            std::list<int> * left = rll_coding(l_left, l, left_restr, 1, start);
            std::list<int> * right = rll_coding(l_right, r, 1, right_restr, i + 1);
            merge(root, *left, *right, str);
            delete left;
            delete right;

            l_left++;
            l_right--;
            str >>= 1;
        }
    }
    else
    { // create sequences with one '1'
        if (left_restr & right_restr)
        {
            jmin = start + min_constraint_d;
            if ((c = start + length - max_constraint_k - 1) > jmin)
                jmin = c;
            jmax = start + max_constraint_k;
            if ((c = start + length - min_constraint_d - 1) < jmax)
                jmax = c;
        }
        else if (left_restr)
        {
            jmin = start + min_constraint_d;
            if ((c = start + length - max_constraint_k - 1) > jmin)
                jmin = c;
            jmax = start + max_constraint_k;
            if ((c = start + length - 1) < jmax)
                jmax = c;
        }
        else if (right_restr)
        {
            jmin = start + length - max_constraint_k - 1;
            if (jmin < start)
                jmin = start;
            jmax = start + max_constraint_k;
            if ((c = start + length - min_constraint_d - 1) < jmax )
                jmax = c;
        }
        else
        {
            jmax = start + max_constraint_k;
            jmin = start + length - max_constraint_k - 1;
        }

        str = 1 << (32 - jmin - 1);
        for (int i = jmin; i <= jmax; i++)
        {
            root->push_back(str);

            str >>= 1;
        }
    }
    return root;
}

//--------------------------------------------------------------------------

void CodewordMap::merge(std::list<int> * root, std::list<int> & left, std::list<int> & right, unsigned int str)
{
    if (left.empty())
    {
        if (right.empty())
        {
            root->push_back(str);
        }
        else
        {
            for (const auto & it : right)
            {
                root->push_back(it | str);
            }
        }
    }
    else if (right.empty())
    {
        for (const auto & it : left)
        {
            root->push_back(it | str);
        }
    }
    else
    {
        for (const auto &it_left : left)
        {
            for (const auto & it_right : right)
            {
                root->push_back(it_left | it_right | str);
            }
        }
    }
}

//--------------------------------------------------------------------------

static int get_ntz(int number, int number_length)
{
    int ntz = 0;
    for (ntz = 0; (ntz < number_length) && (!((1 << ntz) & number)); ntz++);
    return ntz;
}

//--------------------------------------------------------------------------

static int get_nlz(int number, int number_length)
{
    int nlz;
    for (nlz = 0; (nlz < number_length) && !((1<<(number_length-nlz-1)) & number); nlz++);
    return nlz;
}

//--------------------------------------------------------------------------

static void get_rds_and_bit_var(int number, int number_length, int & rds, uint32_t & bit_var, int32_t & sum_per_bit_rds)
{
    int sign = 1;

    rds = 0;
    bit_var = 0;
    sum_per_bit_rds = 0;

    for (int j = 0; j < number_length; j++)
    {
        int mask = 1 << (number_length - j - 1);
        if (mask & number)
        {
            sign = -sign;
        }
        rds += sign;
        sum_per_bit_rds += rds;
        bit_var += rds * rds;
    }
}

//--------------------------------------------------------------------------
