#define DEBUG 0
#include "matrix-conv.hpp"
#include "matrix-io.hpp"
#include "ioparser.hpp"
#include "clparser.hpp"
using namespace toolbox;
#include <fstream>
void banner() 
{
    std::cerr
            << " USAGE: matconv [-in-type auto|crs|full|coord]                                  \n"
            << "                 -out-type crs|full|coord                                       \n"
            << "                [-in filename] [-out filename] [-h] [ < in ] [ > out ]          \n"
            << "                                                                                \n"
            << " converts a matrix between two different input formats. the format for reading  \n"
            << " can be either guessed from the file, of provided via the -in-type option.      \n"
            << " output format must be provided via the -out-type option.                       \n"
            << " input and output can be performed either via stdin and stdout, or can be       \n"
            << " performed over the files given in the -in and -out options.                    \n"
            << "                                                                                \n";
}

int main(int argc, char **argv)
{
    std::string in_fmt, out_fmt, in_fname, out_fname;
    CLParser clp(argc, argv);
    bool fhelp;
    bool fok=
    clp.getoption(in_fmt,"in-type",std::string("auto")) &&
    clp.getoption(out_fmt,"out-type",std::string("full")) &&
    clp.getoption(in_fname,"in",std::string("")) &&
    clp.getoption(out_fname,"out",std::string("")) &&
    clp.getoption(fhelp,"h",false) ;
    if ( fhelp || ! fok) { banner(); return 0; }
    
    std::ifstream sin; bool finfile=false;
    if (in_fname !="")
    {   finfile=true; sin.open(in_fname.c_str());  if (!sin.good()) ERROR("Unable to read input file.\n");  }
    
    std::ofstream sout; bool foutfile=false;
    if (out_fname !="")
    {   foutfile=true; sout.open(out_fname.c_str());  if (!sout.good()) ERROR("Unable to write to output file.\n");  }
    
    std::cout.precision(10);
    sout.precision(10);
    if (in_fmt=="auto")
    {
        if (!io_label_seek(finfile?sin:std::cin,"mode"))
            ERROR("Unable to read mode tag in matrix file. Unable to determine type.");
        (finfile?sin:std::cin) >> in_fmt;
    }

    if (in_fmt=="crs")
    {
        CrsMatrix<double> min;
        (finfile?sin:std::cin) >> min;
        if (out_fmt == "crs") 
        { (foutfile?sout:std::cout)<<min; }
        else if (out_fmt == "coord") 
        {  CoordMatrix<double> mout(min); (foutfile?sout:std::cout)<<mout; }
        else if (out_fmt == "full") 
        {  FMatrix<double> mout(min); (foutfile?sout:std::cout)<<mout; }
        else ERROR("Unsupported matrix output format.");
    }
    else if (in_fmt=="coord")
    {
        CoordMatrix<double> min;
        (finfile?sin:std::cin) >> min;
        if (out_fmt == "crs") 
        { CrsMatrix<double> mout(min); (foutfile?sout:std::cout)<<mout; }
        else if (out_fmt == "coord") 
        { (foutfile?sout:std::cout)<<min; }
        else if (out_fmt == "full") 
        {  FMatrix<double> mout(min); (foutfile?sout:std::cout)<<mout; }
        else ERROR("Unsupported matrix output format.");
    }
    else if (in_fmt=="full")
    {
        FMatrix<double> min;
        (finfile?sin:std::cin) >> min;
        if (out_fmt == "crs") 
        { CrsMatrix<double> mout(min); (foutfile?sout:std::cout)<<mout; }
        else if (out_fmt == "coord") 
        {  CoordMatrix<double> mout(min); (foutfile?sout:std::cout)<<mout; }
        else if (out_fmt == "full") 
        { (foutfile?sout:std::cout)<<min; }
        else ERROR("Unsupported matrix output format.");
    }
    else ERROR("Unsupported matrix input format.");
    return 1;
}
