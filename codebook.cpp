//codebook.cpp

#include <fstream>

#include "codebook.h"


//--------------------------------------------------------------------------

Codebook::Codebook(CodebookArray &codebook,
                  std::shared_ptr<CodewordMap> codeword_map_)
    :line(codebook),
     codeword_map(codeword_map_)
{
    is_valid();
}

//--------------------------------------------------------------------------


bool Codebook::is_valid() const
{
    return states_have_valid_codeword() && is_consistent();
}

//--------------------------------------------------------------------------

bool Codebook::states_have_valid_codeword() const
{
    for (uint32_t i = 0; i < NUM_INDICES ; i++)
    {
        for (uint32_t state = 0; state < NUM_EFMPLUS_STATES ; state++)
        {
            uint16_t code1 = line[i][state].codeword;
            uint16_t code2 = line[(i + NUM_INDICES)][state].codeword;
            if (code1 == INVALID_SEQUENCE && code2 == INVALID_SEQUENCE)
            {
                char exception[1024];
                snprintf(exception, sizeof(exception), "Error in index %d state %d", i, state);
                throw Exception(&exception[0]);
            }
        }
    }
    return true;
}

//--------------------------------------------------------------------------

bool Codebook::states_have_valid_codeword(const std::array<uint32_t, 2> & indices) const
{
    for (auto & i : indices)
    {
        for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
        {
            uint16_t code1 = line[i][state].codeword;
            uint16_t code2 = line[(i + NUM_INDICES) % NUM_INDICES][state].codeword;
            if (code1 == INVALID_SEQUENCE && code2 == INVALID_SEQUENCE)
            {
                return false;
            }
        }
    }
    return true;
}

//--------------------------------------------------------------------------

//codewords with the same next state should be in the same line
bool Codebook::is_consistent() const
{
    for (uint32_t i = 0; i < CDB_LINES - 1;  i++)
    {
        for (uint32_t state1 = 0; state1 < NUM_EFMPLUS_STATES; state1++)
        {
            uint16_t codeword1 = line[i][state1].codeword;
            EfmplusState next_state1 = line[i][state1].next_state;

            if (codeword1 != INVALID_SEQUENCE)
            {
                for (uint32_t j = i + 1; j < CDB_LINES; j++)
                {
                    for (uint32_t state2 = 0; state2 < NUM_EFMPLUS_STATES; state2++)
                    {
                        uint16_t codeword2 = line[j][state2].codeword;
                        EfmplusState next_state2 = line[j][state2].next_state;
                        if (codeword1 == codeword2 && next_state1 == next_state2)
                        {
                            print_codebook_as_header("improved_EFMplus_table.h",0);
                            char exception[1024];
                            snprintf (exception, sizeof(exception), "indices %d %d inconsistent codeword 0x%x", i, j, codeword1);
                            throw Exception(exception);
                        }
                    }
                }
            }
        }
    }
    return true;
}

//--------------------------------------------------------------------------

void Codebook::print_line(FILE * f, int line_index) const
{
    fprintf(f, "/*%3d*/", line_index);
    fprintf(f, "{ ");
    for (uint32_t state = 0; state < NUM_EFMPLUS_STATES; state++)
    {
        fprintf(f, "{0x%04x, %d}, ", line[line_index][state].codeword, int(line[line_index][state].next_state));
    }
    fprintf(f, "},\n");
}

//--------------------------------------------------------------------------

void Codebook::print_codebook_as_header(const char * name, int flag) const
{
    FILE * file;

    fopen_s(&file, name, "w");
    if (!file)
        return;


    fprintf(file, "CodebookArray EFMplus_encoding_table = \n");
    fprintf(file, "{\n");
    int index;

    if (flag == 0)
    {
        for (int i = 0; i < 2; i++) // print MAIN and ALTERNATE lines
            for (index = 0; index < NUM_INDICES; index++)
                print_line(file, index + (i<<8));
    }

    else if (flag == 1)
    {
        for (int i = 1; i >= 0; i--)
            for (index = 0; index < NUM_INDICES; index++)
                print_line(file, index + (i<<8));
    }

    else if (flag==2)
    {
        int i=1;
        for(index = 0; index<123; index++)
            print_line(file, index + (i<<8));

        i=0;
        for(index = 123; index < NUM_INDICES; index++)
            print_line(file, index + (i<<8));

        for(index = 0; index<123; index++)
            print_line(file, index + (i<<8));

        i=1;
        for(index = 123; index < NUM_INDICES; index++)
            print_line(file, index + (i<<8));
    }

    else if (flag == 3)
    {
        int i=0;
        for(index = 0; index<64; index++)
                print_line(file, index + (i<<8));

        i=1;
        for(index = 64; index<192; index++)
            print_line(file, index + (i<<8));

        i=0;
        for(index = 192; index<NUM_INDICES; index++)
            print_line(file, index + (i<<8));

        i=1;
        for(index = 0; index<64; index++)
            print_line(file, index + (i<<8));

        i=0;
        for(index = 64; index<192; index++)
            print_line(file, index + (i<<8));

        i=1;
        for(index = 192; index<NUM_INDICES; index++)
            print_line(file, index + (i<<8));
    }

    fprintf(file, "};");
    fclose(file);
}

//--------------------------------------------------------------------------

void random_file()
{
    std::ofstream random_file("random_file1", std::ios::out | std::ios::binary);
    char buffer[1024];

    if (!random_file.is_open())
        return;

    for (int j = 0; j < 3000; j++)
    {
        for (int i = 0; i < (1<<10); i++)
        {
            buffer[i] = rand();
        }
        random_file.write(buffer, 1024);
    }
    random_file.close();
}

//--------------------------------------------------------------------------
