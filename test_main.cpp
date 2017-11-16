// test_main.cpp

#include "codebook.h"
#include "CodewordMap.h"
#include "sony.h"

int main(int argc, char **argv)
{
    auto codeword_map = std::make_shared<CodewordMap> ( );
    fprintf(stderr, "hello from main\n");
    Codebook * codebook = new Codebook(EFMplus_encoding_table, codeword_map);
    fprintf(stderr, "codebook done\n");

    codebook->calculate_pdf_with_transition_table();
    //fprintf(stderr, "calculated\n");

}
