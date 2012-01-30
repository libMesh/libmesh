#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>

namespace
{
    const std::string versionString = "Function Parser for C++ v";
}

void writeVersion(std::fstream& is, const std::string& version)
{
    std::string line;
    while(std::getline(is, line))
    {
        const size_t ind = line.find(versionString);
        if(ind != line.npos)
        {
            is.seekp(size_t(is.tellg()) - line.length() +
                     ind + versionString.length() - 1);
            is.write(version.c_str(), version.length());
            is.seekg(is.tellp());
        }
    }
}

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <version string> <files>\n";
        return 1;
    }

    std::string version(argv[1]);
    if(version.length() < 6) version.resize(6, ' ');

    for(int i = 2; i < argc; ++i)
    {
        std::fstream is(argv[i]);
        if(!is)
        {
            std::perror(argv[i]);
            return 1;
        }

        writeVersion(is, version);
    }

    return 0;
}
