#include "tbdefs.hpp"
#include "tools-autocorr.hpp"
#include "clparser.hpp"

using namespace toolbox;
#include <fstream>
#include <iomanip>
                 void banner()
{
    std::cerr
            << " USAGE: crosscorr -maxlag max-lag [-timestep unit-time]                         \n"
            << "                [-drop ndrop] [-raw]                                            \n"
            << "                [-mean_a exact_mean] [-mean_b exact_mean] [-s]  [-h] [-q]       \n"
            << "                < in  > out                                                     \n"
            << " compute the cross-correlation function of a series of real data couples in     \n"
            << " standard input. the mean and sigma used to compute and normalize the ACF can   \n"
            << " be given as parameters, or computed from the data sets.                        \n"
            << "                                                                                \n";
}


int main(int argc, char **argv)
{
    CCOptions<CrossCorrelation<double> > acopts;

    CLParser clp(argc, argv);
    bool fhelp, fshort, fquiet, fraw;
    unsigned long ncorr, ndrop;
    bool fok=
            clp.getoption(ncorr,"maxlag") &&
            clp.getoption(ndrop,"drop",(unsigned long) 0) &&
            clp.getoption(acopts.timestep,"timestep",1.) &&
            clp.getoption(acopts.f_exact_mean,"mean_a",false) &&
            clp.getoption(acopts.xmean_a,"mean_a",0.) &&
            clp.getoption(acopts.xmean_b,"mean_a",0.) &&
            clp.getoption(fraw,"raw",false) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fquiet,"q",false) &&
            clp.getoption(fshort,"s",false);

    if ( fhelp || ! fok) { banner(); return 0; }

    CrossCorrelation<double> CC(ncorr, acopts);

    double va, vb; unsigned long idata=0;
    while (std::cin>>va>>vb) {if (++idata>ndrop) CC.add(va,vb); }

    double mean_a=CC.mean_a(), sigma_a=CC.sigma_a();
    double mean_b=CC.mean_b(), sigma_b=CC.sigma_b();
    double cross=CC.cross_ab();
    std::cout.precision(6);
    /*
    std::cout<<(fshort?"":"#")<<std::setw(12)<<mean<<" "
            <<std::setw(12)<<sigma<<" "
            <<std::setw(12)<<tau<<" "
            <<std::setw(12)<<AC.actime2()<<"\n";
    if (fshort) return 0;
    std::cout<<"# ^  mean  ^ . ^  sigma ^ . ^   tau  ^ . ^  tau2  ^ .                              \n";
    */
    std::cout.precision(5);
    std::cout<<"################################################################################\n";
    std::cout<<"# cross-correlation function computed from "<<std::setw(10)<<CC.samples()<<" samples                  #\n";
    std::cout<<"#  computed with time unit:    "<<std::setw(10)<<acopts.timestep<<"                                      #\n";
    if (acopts.f_exact_mean)
        std::cout<<"#  computed with exact mean:   "<<std::setw(10)<<acopts.xmean_a<<","<<std::setw(10)<<acopts.xmean_b<<"                           #\n";
    std::cout<<"# **************************************************************************** #\n";
    std::cout<<"#  means:     "<<std::setw(10)<<mean_a
            <<" , "<<std::setw(10)<<mean_b
            <<"                                          #\n";
    std::cout<<"#  sigmas:    "<<std::setw(10)<<sigma_a<<" , "<<std::setw(10)<<sigma_b<<"                                          #\n";
    if (fraw)
        std::cout<<"#  non-normalized autocorrelation output                                       #\n";
    std::cout<<"################################################################################\n";
    std::cout<<"#    time   .     acf    .                                       \n";

    std::cout.precision(6);
    std::valarray<double> t,v;
    CC.fullanalysis(t,v);
    if (fraw) v*=(sigma_a*sigma_b);
    for (unsigned long i=0; i<ncorr; ++i)
    {
        std::cout <<std::setw(12)<< t[i]<<" "
                <<std::setw(12)<<v[i]<<"\n";
    }
    return 0;
}
