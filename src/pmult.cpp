#define TB_MPI 1
#include "tbdefs.hpp"
#include "matrix-pcrs.hpp"
#include <cstring>
#include "matrix-io.hpp"

using namespace toolbox;
int main(int argc, char ** argv)
{
#ifdef TB_MPI
    MPI_Init( &argc, &argv );
#endif
    
    std::string ma(argv[1]);
    std::string mb(argv[2]);
    PCrsMatrix<double>A,B,C;
    
    std::pout<<"Reading matrix A."<<std::endl;
    A.read(ma);
    std::pout<<"Reading matrix B."<<std::endl;
    B.read(mb);
    
    std::pout<<trace(A)<< " TRA\n";
    std::pout<<trace(B)<< " TRB\n";
    
    std::pout<<"Adding"<<std::endl;
    add(A,B,C);
    std::pout<<"Increment"<<std::endl;
    C+=B;
    std::pout<<"Multiplying"<<std::endl;
    A.trunc(0.);
    B.trunc(0.);
    
    mult(A,B,C);
    std::pout<<"Printing out result:"<<std::endl;
    C.dump(std::string("pmult-out"));
    std::pout<<"Success."<<std::endl;
    
    TBBenchmarks.print(std::pout);
#ifdef TB_MPI
    MPI_Finalize();
#endif
    exit(0);
}

