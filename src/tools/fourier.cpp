/* This does not really use much of the toolbox infrastructure, 
   but is a handy interface to do Fourier transform from the CLI
   such as in BASH scripts and the like.
   */
   
#include "clparser.hpp"
#include "fftw3.h"
using namespace toolbox;

void banner()
{
   std::cerr<< " USAGE: fourier [options] < INPUT [ > OUTPUT ]                                  \n"
            << " reads a time series from stdin and computes its Fourier transform by FFT.      \n"
            << " Input should be just the time series, without headers or time code.            \n"           
            << " Data is interpreted as a sequence of points in time, spaced by an interval dt  \n"
            << " and starting from t=0. The FT is performed considering the function to be an   \n"
            << " even function of time, so it yields sqrt(2/pi)int_0^infty f(t) dt            \n"
            << " -dt [dt]    sets the time interval between input samples. {def: 1.0}           \n"
            << " -pad [npad] appends npad zeroes before doing the FT. increases the resolution. \n"
            << " -win [wnd]  applies a windowing function to the data before FT.                \n"
            << "             possible values: triangle | cosine | hanning | gauss-[2,3,4,5,6] | \n"
            << "             blackman-harris                                                    \n"
            << "                                                                                \n";
}

int main(int argc, char ** argv)
{
   CLParser clp(argc, argv);
   
   double dt=1.0; 
   unsigned long npad=0;
   std::string wnd;
   bool fhelp;
   bool fok=clp.getoption(dt, "dt", 1.0) &
            clp.getoption(fhelp, "h", false) &
            clp.getoption(npad, "pad", (unsigned long) 0) &   
            clp.getoption(wnd, "win", std::string(""));
   
   if (!fok || fhelp) { banner(); return -1; }
   
   std::vector<double> data(0); double y;
   // fills up the time series array from stdin 
   while (std::cin.good()) { std::cin >> y; data.push_back(y); }
   
   // prepares for fft: transfers into valarray, pads, boxes etc
   unsigned long ndata=data.size(), nfft=ndata+npad;
   std::valarray<double> vvt(nfft), vvw(nfft);
   vvt=0.0;  for (unsigned long i=0 ; i<ndata; ++i) vvt[i]=data[i];

   if (wnd=="triangle")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=(1.-double(i)/double(ndata-1));
   else if (wnd=="cosine")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=cos(0.5*constant::pi*double(i)/double(ndata-1));   
   else if (wnd=="hanning")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=0.5*(1.0+cos(constant::pi*double(i)/double(ndata-1)));   
   else if (wnd=="gauss-2")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=exp(-pow(double(i)/(2.0*double(ndata)/2.0),2));   
   else if (wnd=="gauss-3")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=exp(-pow(double(i)/(2.0*double(ndata)/3.0),2));   
   else if (wnd=="gauss-4")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=exp(-pow(double(i)/(2.0*double(ndata)/4.0),2));         
   else if (wnd=="gauss-5")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=exp(-pow(double(i)/(2.0*double(ndata)/5.0),2));         
   else if (wnd=="gauss-6")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=exp(-pow(double(i)/(2.0*double(ndata)/6.0),2));         
   else if (wnd=="blackman-harris")
      for (unsigned long i=0 ; i<ndata; ++i) vvt[i]*=(35875 + 48829*cos((i*constant::pi)/(-1. + ndata)) + 14128*cos((2*i*constant::pi)/(-1. + ndata)) + 1168*cos((3*i*constant::pi)/(-1. + ndata)))/   100000.;   
   else if (wnd!="") ERROR("Window function "<< wnd << " mispelled or not implemented.");
   // FFTW calls
   fftw_plan r2rplan=fftw_plan_r2r_1d(nfft, &vvt[0], &vvw[0], FFTW_REDFT00, FFTW_ESTIMATE);
   fftw_execute(r2rplan);        
   
   // fixes the normalization
   vvw*=dt/sqrt(2.*constant::pi);
   double dw=constant::pi/(dt*(nfft-1));

   // very terse output
   std::cout.precision(8); std::cout.width(15); std::cout.setf(std::ios::scientific); 
   for (unsigned long i=0; i<nfft; ++i) std::cout<<i*dw<<"  "<<vvw[i]<<std::endl;
   
  return 0;               
}

