#include "tbdefs.hpp"
#include "mpi.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    char mpistr[MPI_MAX_ERROR_STRING];
    int ec=toolbox::str2int(std::string(argv[1])), sl;
    MPI_Error_string(ec,mpistr,&sl);
    std::cout<<"Error code "<<ec<<" corresponds to\n"<<mpistr<<"\n";
    return 0;
}
