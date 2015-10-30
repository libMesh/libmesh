#include <cstdio>
#include <iostream>
#include "cpp_compress.hh"

int main(int argc, char* argv[])
{
    std::string out;
    for(;;)
    {
        char Buf[32768];
        if(!std::fgets(Buf, sizeof(Buf), stdin)) break;
        Buf[(sizeof Buf)-1] = '\0';
        out += Buf;
    }
    CPPcompressor Compressor;

    std::string prefix;
    if(argc == 2 && argv[1][0] != '-') prefix = argv[1];
    std::cout << Compressor.Compress(out, prefix);
}
