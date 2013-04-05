#include "tbdefs.hpp"
#include "tools-histogram.hpp"
#include "clparser.hpp"

using namespace toolbox;
#include <fstream>
void banner()
{
    std::cerr
            << " USAGE: histogram -xi xi -xf xf [-n nbins] [-whard ] [-w] [-avg]  \n"
            << "                  [(-b box-w|-t tri-w|-g1 g1-w|-g2 g2-w|-g3 g3-w|-g5 g5-w)]      \n"
            << " compute the histogram of a series of data, with nbins (default:100) bins       \n"
            << " distributed evenly between xi and xf. optionally, box (b) or triangular (t)  \n"
            << " smoothing functions can be used. If -w is selected, then a weight is read     \n"
            << " after each point (input must be value weigth value weight ...).              \n"
            << " [-whard] states that hard walls should be used (i.e. density won't spillout if\n"
            << " the data point is inside the interval.                                       \n"
            << " If [-avg] is selected, on each line it will be read x y [w] and the average   \n"
            << " of y will be evaluated as a function of x, in a binned fashion.              \n";
}

int main(int argc, char **argv)
{
    HGOptions<Histogram<double> > hgopts;

    CLParser clp(argc, argv);
    bool fhelp, fweighted, faverage, fhard;
    double a,b,wb,wt,wg1,wg2,wg3,wg5;
    unsigned long nbins;
    bool fok=
            clp.getoption(nbins,"n",(unsigned long) 100) &&
            clp.getoption(a,"xi") &&
            clp.getoption(b,"xf") &&
            clp.getoption(wb,"b",0.) &&
            clp.getoption(wt,"t",0.) &&
            clp.getoption(wg1,"g1",0.) &&
            clp.getoption(wg2,"g2",0.) &&
            clp.getoption(wg3,"g3",0.) &&
            clp.getoption(wg5,"g5",0.) &&
            clp.getoption(fweighted,"w",false) &&
            clp.getoption(faverage,"avg",false) &&
            clp.getoption(fhard,"whard",false) &&
            clp.getoption(fhelp,"h",false);

    if ( fhelp || ! fok) { banner(); return 0; }

    Histogram<double> HG(a,b,nbins), HGY(a,b,nbins);
    HGOptions<Histogram<double> > hgo;

    HG.get_options(hgo);

    if (wb>0.)       { hgo.window=HGWBox; hgo.window_width=wb; }
    else if (wt>0.)  { hgo.window=HGWTriangle; hgo.window_width=wt; }
    else if (wg1>0.) { hgo.window=HGWGauss1; hgo.window_width=wg1; }
    else if (wg2>0.) { hgo.window=HGWGauss2; hgo.window_width=wg2; }
    else if (wg3>0.) { hgo.window=HGWGauss3; hgo.window_width=wg3; }
    else if (wg5>0.) { hgo.window=HGWGauss5; hgo.window_width=wg5; }
    else             { hgo.window=HGWDelta; hgo.window_width=0.; }

    if (fhard) hgo.walls=HGBHard;
    HG.set_options(hgo); HGY.set_options(hgo);

    double val, weight, y, ty, ny;

    if (faverage)
    {
       ty=ny=0.0;
       if (fweighted)
           while (std::cin>>val) { std::cin>>y;  std::cin>>weight; HG.add(val,weight); HGY.add(val,y*weight); ty+=weight*y; ny+=weight; }
       else
           while (std::cin>>val) { std::cin>>y; HG<<val; HGY.add(val,y); ty+=y; ny+=1.0; }
    }
    else
    {
       if (fweighted)
           while (std::cin>>val) { std::cin>>weight; HG.add(val,weight); }
       else
           while (std::cin>>val) { HG<<val; }
   }

   double above, below;
   HG.get_outliers(above,below);
   std::cout<<"# Total weight: "<< HG.get_totweight()<<std::endl;
   std::cout<<"# Fraction below: "<<below<<std::endl;
   std::cout<<"# Fraction above: "<<above<<std::endl;

   std::cout.precision(12);
   std::cout.setf(std::ios::scientific);
   std::cout.width(14);
   if (faverage)
   {
        std::valarray<double> wx, ww, wf, wy;
        HG.get_bins(wx,ww,wf);
        HGY.get_bins(wx,ww,wy);
        double ay=ty/ny;
        for (unsigned int i=0; i<wx.size(); ++i)
            std::cout<<wx[i]<<" "<<wy[i]/wf[i]*ay<<" "<<wf[i]<<" "<<ww[i]<<"\n";
   }
   else std::cout <<HG;
}
