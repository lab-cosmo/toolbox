#include "tbdefs.hpp"
#include "tools-autocorr.hpp"
#include "clparser.hpp"

using namespace toolbox;
#include <fstream>
#include <iomanip>

#define PROGRESSSTRIDE 100


void banner() 
{
    std::cerr
            << " USAGE: autocorr -maxlag max-lag [-timestep unit-time]                          \n"
            << "                [-drop ndrop] [-coll] [-runlength lr] [-stride std]             \n"
            << "                [-mean exact_mean] [-sigma exact-sigma] [-s] [-raw]             \n"
            << "                [-tau exact_tau] [-tau2 exact-tau2] [-h] [-v]                   \n"
            << "                [< in | -input file1 file2 ...] > out                           \n"
            << " compute the autocorrelation function of a series of real data points given in  \n"
            << " standard input. the mean and sigma used to compute and normalize the ACF can   \n"
            << " be given as parameters, or computed from the data sets.                        \n"
            << " one can provide several series of data, either from different files (-input)   \n"
            << " or concatenating the input and providing the length of each run (-runlength).  \n"
            << " if -coll is specified, series are collapsed when finished to save memory.      \n"
            << " the error is computed from the integral of the ACF^2, which can be given as    \n"
            << " input. if the -s option is selected, the code writes only mean, sigma, tau and \n"
            << " tau2 on a single line, then exit. -v enables verbose output on stderr,         \n"
            << " including progress output. -raw prints non-normalized ACF.                     \n"
            << "                                                                                \n";
}


int main(int argc, char **argv)
{
    ACOptions<AutoCorrelation<double> > acopts;
    
    CLParser clp(argc, argv);

    bool fhelp, fshort, fverbose, fcoll, fraw;
    std::vector<std::string> datafiles;
    unsigned long ncorr, ndrop, runl, stride;
    bool fok=
            clp.getoption(ncorr,"maxlag") &&
            clp.getoption(ndrop,"drop",(unsigned long) 0) &&
            clp.getoption(runl,"runlength",(unsigned long) 0) &&
            clp.getoption(stride,"stride",(unsigned long) 1) &&
            clp.getoption(fcoll,"coll",false) &&
            clp.getoption(fraw,"raw",false) && 
            clp.getoption(datafiles,"input",std::vector<std::string>(0)) &&
            clp.getoption(acopts.timestep,"timestep",1.) &&
            clp.getoption(acopts.f_exact_sigma,"sigma",false) &&
            clp.getoption(acopts.xsigma,"sigma",1.) &&
            clp.getoption(acopts.f_exact_tau,"tau",false) &&
            clp.getoption(acopts.xtau,"tau",1.) &&
            clp.getoption(acopts.f_exact_tau2,"tau2",false) &&
            clp.getoption(acopts.xtau2,"tau2",1.) &&
            clp.getoption(acopts.f_exact_mean,"mean",false) &&
            clp.getoption(acopts.xmean,"mean",0.) &&
            clp.getoption(fhelp,"h",false) &&
            clp.getoption(fverbose,"v",false) &&
            clp.getoption(fshort,"s",false);
    
    if ( fhelp || ! fok || ! clp.chkunknown(true)) { banner(); return 0; }
    if (stride==0) stride=1;
    AutoCorrelation<double> AC(ncorr, acopts);

    double val; unsigned long idata=0, cseries=0;
    if (datafiles.size()==0)
    {
        while (std::cin>>val) {
            if (idata==0 && cseries>0) AC.series_next();
            if (++idata>ndrop && (idata-ndrop-1)%stride==0) AC<<val; 
            if (fverbose && idata%PROGRESSSTRIDE==0) std::cerr << " Read: "<< std::setw(15)<< idata<<" data points on data series "<<std::setw(8)<<cseries+1<<".\r";
            if (runl>0 && idata%runl==0)
            {
                idata=0; cseries++; 
                if (fverbose) std::cerr << " Collapsing data from series...\n";
                if (fcoll && cseries>1) AC.series_collapse();
            }
        }
    }
    else
    {
        if (fverbose) std::cerr << " Input data from file(s)\n";
        std::vector<std::string>::iterator it=datafiles.begin();
        while (it!=datafiles.end())
        {
            if (it!=datafiles.begin()) AC.series_next();
            std::ifstream ifile((*it).c_str());
            idata=0;
            while (ifile>>val) if (++idata>ndrop && (idata-ndrop-1)%stride==0) AC<<val;
            if (fverbose) std::cerr << " Read " <<std::setw(10)<<(idata-ndrop)<<" items from file "<<(*it)<<".\r";
            ++it; ++cseries;
            ifile.close();
            if (fcoll && cseries>1) 
            {
                if (fverbose) std::cerr << " Collapsing data from series...\n";
                AC.series_collapse();
            }
        }
    }
    
    if (fverbose) std::cerr << " Finished reading. Now computing ACF and outputting ...\n";
    double mean=AC.mean(), sigma=AC.sigma(), tau=AC.actime();
    std::cout.precision(8); std::cout.setf(std::ios::scientific); 
    std::cout<<(fshort?"":"# ")
            <<std::setw(15)<<mean<<" "
            <<std::setw(15)<<sigma<<" "
            <<std::setw(15)<<tau<<" "
            <<std::setw(15)<<AC.actime2()
            <<std::setw(15)<<(sigma/sqrt(AC.samples()*acopts.timestep*0.5/tau))<<" "
            <<"\n";
    if (fshort) return 0;
    std::cout<<"# ^  mean  ^ . ^  sigma ^ . ^   tau  ^ . ^  tau2  ^ . ^  error  ^ .                       \n";
    std::cout.precision(7); std::cout.setf(std::ios::scientific);
    std::cout<<"################################################################################\n";
    std::cout<<"#  autocorrelation function computed from "<<std::setw(12)<<AC.samples()<<" samples                   #\n";
    std::cout<<"#  computed with time unit:    "<<std::setw(12)<<acopts.timestep<<"                                      #\n";
    if (acopts.f_exact_mean)
        std::cout<<"#  computed with exact mean:   "<<std::setw(12)<<acopts.xmean<<"                                      #\n";
    if (acopts.f_exact_sigma)
        std::cout<<"#  computed with exact sigma:  "<<std::setw(12)<<acopts.xsigma<<"                                      #\n";
    if (acopts.f_exact_tau2)
        std::cout<<"#  computed with exact tau2:   "<<std::setw(12)<<acopts.xtau2<<"                                      #\n";
    if (fraw) 
        std::cout<<"#  non-normalized autocorrelation output                                       #\n";
    std::cout<<"# **************************************************************************** #\n";
    std::cout<<"#  mean:      "<<std::setw(12)<<mean
            <<" ± "<<
                std::setw(12)<<(sigma/sqrt(AC.samples()*acopts.timestep*0.5/tau))
            <<"                                          #\n";
    std::cout<<"#  sigma:     "<<std::setw(12)<<sigma<<"                                                       #\n";
    std::cout<<"#  a.c. time: "<<std::setw(12)<<tau<<"                                                       #\n"; 
    std::cout<<"################################################################################\n";
    std::cout<<"#    time   .     acf    .    Dacf    .    block   .   Dblock   .\n";
    
    std::cout.precision(6);
    std::valarray<double> t,v,b,ev,eb;
    
    if (fverbose) std::cerr << " Writing ACF and blocking analysis data ...\n";    
    AC.fullanalysis(t,v,b,ev,eb);
    if (fraw) v*=(sigma*sigma);
    for (unsigned long i=0; i<ncorr; ++i)
    {
        std::cout <<std::setw(12)<< t[i]<<" "
                <<std::setw(12)<<v[i]<<" "
                <<std::setw(12)<<ev[i]<<" "
                <<std::setw(12)<<b[i]<<" "
                <<std::setw(12)<<eb[i]<<"\n";
    }
    return 0;
}
