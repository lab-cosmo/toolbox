#include "tbdefs.hpp"
#include "tools-histogram.hpp"
#include "clparser.hpp"
#include "ioparser.hpp"

using namespace toolbox;
#include <fstream>
void banner()
{
    std::cerr
            << " USAGE: ndhistogram -d dims -xi 'xi1,...,xin' -xf 'xf1,...,xfn'                 \n"
            << "                    [-n 'n1,...,nn'] [-(t|b|g1|g2|g3|g5) 'w1,...,wn']           \n"
            << "                    [-w] [-g] [-avg] [-whard|-wperi] [-adaptive err]            \n"
            << "                                                                                \n"
            << " compute the histogram of a series of data, in ndim dimensions.                 \n"
            << " data must be formatted as                                                      \n"
            << " D1(1)   D2(1)  .......   DN(1) [VALUE(1)] [ WEIGHT(1) ]                        \n"
            << " on every dimension n (default:100) bins are distributed evenly.                \n"
            << " between xi and xf. optionally, box (b), triangular (t) or Gaussian (-g?)  \n"
            << " smoothing functions can be used. w? determines the window in each dimension \n"
            << " [-whard] states that hard walls should be used (i.e. density won't spillout if\n"
            << " the data point is inside the interval).  [-wperi] uses periodic boundaries.    \n"
            << " If -g is selected, points will be output in gnuplot format. \n"
            << " If -avg is selected, then a quantity is read after each point, and the output  \n"
            << "    is the conditional average of that quantity constrained to the value of x   \n"
            << " If -w is selected, then a weight is read after each point.                     \n"
            << " -adaptive uses a simple procedure to adaptively smoothen the histogram as \n"
            << " it is built, trying to have a relative error of err on the histogram value \n"
            << " at each point.\n"
            ;
}


double (* as_radius) (const std::vector<double>& op, double val);

//parameters are r(fmax), r(fmin), fmax, fmin
double asr_linear(const std::vector<double>& op, double v)
{  return op[1]+(op[0]-op[1])*(v-op[3])/(op[2]-op[3]); }

//parameters are r(fmax), rmin, dr/dlogf, fmax
double asr_log(const std::vector<double>& op, double v)
{
   if ( exp(-(op[1]-op[0])/op[2])>v/op[3]) return op[1]; else return op[0]-log(v/op[3])*op[2];
}

double (* as_kernel) (double val);

double ask_linear(double x)
{  return 1-x; }

double ask_square(double x)
{  return (1-x)*(1-x); }

double ask_exp(double x)
{  return exp(-4*x); }


int main(int argc, char **argv)
{

    CLParser clp(argc, argv);
    bool fhelp, fweighted, fgnu, faverage, fhard, fperi;
    std::string a,b,wb,wt,wg1,wg2,wg3,wg5,nbins;
    unsigned long ndim; double aeps;
    bool fok=
            clp.getoption(ndim,"d") &&
            clp.getoption(a,"xi") &&
            clp.getoption(b,"xf") &&
            clp.getoption(wt,"t",std::string("")) &&
            clp.getoption(wb,"b",std::string("")) &&
            clp.getoption(wg1,"g1",std::string("")) &&
            clp.getoption(wg2,"g2",std::string("")) &&
            clp.getoption(wg3,"g3",std::string("")) &&
            clp.getoption(wg5,"g5",std::string("")) &&
            clp.getoption(nbins,"n",std::string("")) &&
            clp.getoption(aeps,"adaptive",0.) &&
            clp.getoption(fweighted,"w",false) &&
            clp.getoption(faverage,"avg",false) &&
            clp.getoption(fhard,"whard",false) &&
            clp.getoption(fperi,"wperi",false) &&
            clp.getoption(fgnu,"g",false) &&
            clp.getoption(fhelp,"h",false);

    if ( fhelp || ! fok) { banner(); return 0; }

    std::valarray<HGOptions<Histogram<double> > >hgo(ndim);
    std::valarray<double> va, vb, vw, vn, vasop; std::vector<double> asop;
    csv2floats(a,va); csv2floats(b,vb); csv2floats(nbins,vn);
    HGWindowMode wmode;
    if (wb!="") { csv2floats(wb,vw); wmode=HGWBox; }
    else if (wt!="") { csv2floats(wt,vw); wmode=HGWTriangle; }
    else if (wg1!="") { csv2floats(wg1,vw); wmode=HGWGauss1; }
    else if (wg2!="") { csv2floats(wg2,vw); wmode=HGWGauss2; }
    else if (wg3!="") { csv2floats(wg3,vw); wmode=HGWGauss3; }
    else if (wg5!="") { csv2floats(wg5,vw); wmode=HGWGauss5; }
    else wmode=HGWDelta;
    std::valarray<double> pva(va), pvb(vb); //extended intervals for periodic
    std::valarray<int> pbin(ndim);

    if (va.size()!=ndim || vb.size()!=ndim ) ERROR("Boundaries must be given for all dimensions.");
    if ((vn.size()!=ndim && vn.size()!=0)||(vw.size()!=ndim && vw.size()!=0)) ERROR("N. bins and width must be given for all dimensions, if specified.");
    for (int i=0; i<ndim; ++i)
    {
        hgo[i].window=wmode;
        hgo[i].window_width=(vw.size()==0 || hgo[i].window==HGWDelta?0.:vw[i]);

        if (fhard) hgo[i].walls=HGBHard;  //!this should be done per direction.
        else if (fperi) hgo[i].walls=HGBPeriodic;
        else hgo[i].walls=HGBNormal;

        hgo[i].boundaries.resize((vn.size()==0?101:vn[i]+1));
        for (int k=0; k<hgo[i].boundaries.size();k++)
             hgo[i].boundaries[k]=va[i]+k*(vb[i]-va[i])/(hgo[i].boundaries.size()-1);

        hgo[i].adaptive_eps = aeps;
    }

    NDHistogram<double> HG(hgo);
    NDHistogram<double> HGY(hgo);

    std::valarray<double> val(ndim); double weight, y, ty, ny;

    ty=ny=0.0;
    while (std::cin.good()) {
         for (int i=0; i<ndim; ++i) std::cin>>val[i];

         if (faverage) std::cin>>y;
         if (fweighted) std::cin>>weight; else weight = 1.0;

         HG.add(val,weight);
         if (faverage) { HGY.add(val,weight*y); ty+=weight*y; ny+=weight; }
    }
    double ay=ty/ny;

    double outliers;
    HG.get_outliers(outliers);
    std::cout<<"# Total weight: "<< HG.get_totweight() <<"\n# Fraction outside: "<<outliers;
    if (vw.size()>0) //!TODO write out which kernel is being used
    { std::cout<< "  -- kernel binning:  widths  "; for (int i=0; i<ndim; i++) std::cout<<vw[i]<<"  "; }
    std::cout<<"\n";
    std::cout.setf(std::ios::scientific);

    if (!fgnu)
    { 
         if (faverage)
         {
             std::cout<<"# Average normalization "<<ay<<"\n"; 
             std::cout<<HGY;
         } 
         else
         {
             std::cout <<HG;
         }
    }
    else
    {
        std::valarray<long> ind(ndim); std::valarray<double> cen(ndim); double val, valy;

        ind = 0;
        while (ind[0]<hgo[0].boundaries.size()-1)
        {
           HG.get_bin(ind,cen,val);
           // prints out the center of the histogram bin
           for (int i=0; i<ndim; i++)
              std::cout<<cen[i]<<"\t";
           if (faverage)
           {
              HGY.get_bin(ind,cen,valy);
              std::cout<<valy/val*ay<<"\n";
           }
           else
           {
              std::cout<<val<<"\n";
           }
           ind[ndim-1]++;
           for (int i=(ndim-1); i>0; --i)
              if (ind[i]>=hgo[i].boundaries.size()-1)
              {
                 ind[i-1]++; ind[i]=0;
                 if (i==1) std::cout<<"\n"; // newline on outer loop (gnuplot compatibility)
              }
        }

    }
}
