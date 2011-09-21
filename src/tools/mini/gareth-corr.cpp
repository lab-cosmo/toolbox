#include "tbdefs.hpp"
#include "tools-autocorr.hpp"
#include "clparser.hpp"
#include "ioparser.hpp"


using namespace toolbox;
#include <fstream>
#include <iomanip>
         
void banner() 
{
    std::cerr
            << " USAGE: gareth-corr -maxlag max-lag [-timestep time-step] < GARETH_FILE         \n"
            << " file must contain:                                                             \n"
            << " N-STEPS  N-SERIES  N-ZEROES                                                    \n"
            << " S01 S02 ...........     S0(N-SERIES - N-ZEROES)    etc...                      \n";
}


int main(int argc, char **argv)
{
    ACOptions<AutoCorrelation<double> > acopts;
    
    CLParser clp(argc, argv);
    bool fhelp;
    unsigned long ncorr, ndrop;
    bool fok=
            clp.getoption(ncorr,"maxlag") &&
            clp.getoption(acopts.timestep,"timestep",1.) &&
            clp.getoption(fhelp,"h",false);
    
    if ( fhelp || ! fok) { banner(); return 0; }
    
    unsigned long nsteps, nseries, nzeroes, ndata;
    std::cin >>nsteps >> nseries >> nzeroes; ndata=nseries-nzeroes;
    std::cerr<<ncorr<<", "<<ndata<<"," <<nsteps<<"\n";
    //reads all data!
    std::valarray<double> ldata(nsteps*ndata); unsigned long k=0;
    std::cerr<<"array allocated\n";
    for (int ist=0; ist<nsteps; ++ist) for (int id=0; id<ndata; ++id) std::cin>>ldata[k++];
    std::cerr<<"array allocated and read\n";
    std::cerr<<ncorr<<", "<<ndata<<"," <<nsteps<<"\n";
    double mean=0., mean2=0.; k=0;
    for (k=0; k<ldata.size(); ++k)  { mean+=ldata[k]; mean2+=ldata[k]*ldata[k]; }
    mean=mean/(ndata*nsteps); mean2=mean2/(ndata*nsteps);
    acopts.xmean=mean; acopts.f_exact_mean=true;
    AutoCorrelation<double>  AC(ncorr, acopts);

    std::valarray<double> gv(nsteps), gb(nsteps), gt, lv, lb, le;
    for (int id=0; id<ndata; ++id) 
    { 
        std::cerr<<"including acf "<<id<<"\n";
        AC.reset(ncorr);
        for (int ist=0; ist<nsteps; ++ist)
        {  AC<<ldata[ist*ndata+id]; }
            
        AC.fullanalysis(gt,lv,lb,le,le);
        std::cerr<<"excerpt from lv: "<<lv[0]<<"\n";
        gv+=lv*AC.sigma()*AC.sigma(); gb+=lb;
    }
    gv+=mean*mean*nzeroes;  gv*=1./nseries;  gv*=1./gv[0];
    
    
    //for (unsigned long i=0; i<ncorr; ++i) gb[i]+=0.5*i*nzeroes*mean*mean; gb*=1./nseries; 
    std::cerr<<"computing block averages\n";
    gb=0.;
    for (unsigned long i=0; i<ncorr; ++i) 
        for (unsigned long j=0; j<i; ++j)
             gb[i]+=(1.-j*1./i)*gv[j];
    gb*=acopts.timestep;
    
    //std::cerr<<gv<<"\n";
    //gb*=1./gv[0]; 
    
    std::cout.precision(6);
    std::cout<<"#"<<std::setw(12)<<mean<<" "
            <<std::setw(12)<<sqrt(mean2-mean*mean)<<"\n";
    std::cout<<"# ^  mean  ^ . ^  sigma ^ .                                                        \n";
    std::cout.precision(5);
    std::cout<<"################################################################################\n";
    std::cout<<"#  autocorrelation function computed from "<<std::setw(10)<<(nsteps*ndata)<<" samples                   #\n";
    std::cout<<"#  computed with time unit:    "<<std::setw(10)<<acopts.timestep<<"                                      #\n";
    std::cout<<"# **************************************************************************** #\n";
    std::cout<<"#  mean:      "<<std::setw(10)<<mean
            <<"                                          #\n";
    std::cout<<"#  sigma:     "<<std::setw(10)<<sqrt(mean2-mean*mean)<<"                                                       #\n";
    std::cout<<"################################################################################\n";
    std::cout<<"#    time   .     acf    .    block   .\n";
    
    std::cout.precision(6);
    for (unsigned long i=0; i<ncorr; ++i)
    {
        std::cout <<std::setw(12)<< gt[i]<<" "
                <<std::setw(12)<<gv[i]<<" "
                <<std::setw(12)<<gb[i]<<"\n";
    }
    return 0;
}
