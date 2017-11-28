// test_main.cpp

#include "codebook.h"
#include "CodewordMap.h"
//#include "sony.h"
#include "optimised_34.h"

int main(int argc, char **argv)
{

    auto codeword_map = std::make_shared<CodewordMap> ( );
    fprintf(stderr, "hello from main\n");
    Codebook * codebook = nullptr;
    try
    {
        codebook = new Codebook(EFMplus_encoding_table, codeword_map);
    }
    catch (const std::exception & e)
    {
        std::cout << e.what();
        exit(-1);
    }
    fprintf(stderr, "codebook done\n");

    codebook->calculate_pdf_with_transition_table();
    for (auto i=0; i<10; i++)
        codebook->optimise();

    //fprintf(stderr, "calculated\n");

}
