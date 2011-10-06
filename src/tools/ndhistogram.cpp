#include "tbdefs.hpp"
#include "tools-histogram.hpp"
#include "clparser.hpp"

using namespace toolbox;
#include <fstream>
void banner() 
{
    std::cerr
            << " USAGE: ndhistogram -d dims -xi 'xi1,...,xin' -xf 'xf1,...,xfn'                 \n"
            << "                    [-n 'n1,...,nn'] [-b 'b1,...,bn'] [-w] [-g]                 \n"
            << "                    [-as mode(parameters)]                                      \n"
            << "                                                                                \n"
            << " compute the histogram of a series of data, in ndim dimensions.                 \n"
            << " data must be formatted as                                                      \n"
            << " D1(1)   D2(1)  .......   DN(1) [ WEIGHT(1) ]                                   \n"
            << " on every dimension n (default:100) bins are distributed evenly.                \n"
            << " between xi and xf. optionally, box (b) smoothing functions can be used.        \n"
            << " If dimension is 2 and -g is selected, points will be output in gnuplot format. \n"
            << " If -w is selected, then a weight is read after each point.                     \n"
            << " If -as is selected, the final histogram will be smoothed adaptively, with      \n"
            << "    a smoothing radius going which depends on the central value.                \n"
            << "    the radius is determined by mode: linear(min,max) interpolates linearly     \n"
            << "    log(min,max,drad) uses logarithm of the ratio and a sigmoid function        \n";
}

void ssplit(const std::string&  istr, std::vector<double>& vl)
{
    vl.clear(); std::string ls=istr;
    int pos=0;
    
    while( (pos = ls.find_first_of(',')) != ls.npos )
    {
        if(pos > 0)
        {
            vl.push_back(str2float(ls.substr(0,pos)));
        }
        ls=ls.substr(pos+1);
    }
    if(ls.length() > 0)
    {
        vl.push_back(str2float(ls));
    }
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
    bool fhelp, fweighted, fgnu;
    std::string a,b,wb,nbins,asmooth,asfun;
    unsigned long ndim;
    bool fok=
            clp.getoption(ndim,"d") &&
            clp.getoption(a,"xi") &&
            clp.getoption(b,"xf") &&
            clp.getoption(wb,"b",std::string("")) &&
            clp.getoption(nbins,"n",std::string("")) &&
            clp.getoption(asmooth,"as",std::string("")) &&
            clp.getoption(asfun,"asf",std::string("linear")) &&
            clp.getoption(fweighted,"w",false) &&
            clp.getoption(fgnu,"g",false) &&
            clp.getoption(fhelp,"h",false);
    
    if ( fhelp || ! fok) { banner(); return 0; }
    
    std::valarray<HGOptions<Histogram<double> > >hgo(ndim);
    std::vector<double> va, vb, vw, vn, asop;
    ssplit(a,va); ssplit(b,vb); ssplit(wb,vw); ssplit(nbins,vn);
    
    std::string asmode, aspars;
    if (asmooth.size()>0) {
      asmode=asmooth.substr(0,asmooth.find("(",0)); aspars=asmooth.substr(asmooth.find("(",0)+1,asmooth.find(")",0)-1-asmooth.find("(",0));
      ssplit(aspars,asop); 
      if (asmode==std::string("linear")) as_radius=&asr_linear; 
      else if (asmode==std::string("log")) as_radius=&asr_log;
      else ERROR("Undefined smoothing mode "<<asmode); 
      
      if (asfun==std::string("linear")) as_kernel=&ask_linear;
      else if (asfun==std::string("square")) as_kernel=&ask_square;
      else if (asfun==std::string("exp")) as_kernel=&ask_exp;
      else ERROR("Undefined smoothing function "<<asfun); 
      
          std::cerr<<asfun<<" aops "<<asmode<<" "<<aspars<<"   "<<asop[0]<<"::"<<asop[1]<<"\n";
    }
    
    if (va.size()!=ndim || vb.size()!=ndim ) ERROR("Boundaries must be given for all dimensions.");
    if ((vn.size()!=ndim && vn.size()!=0)||(vw.size()!=ndim && vw.size()!=0)) ERROR("N. bins and width must be given for all dimensions, if specified.");
    for (int i=0; i<ndim; ++i)
    {
        hgo[i].window=(vw.size()==0?HGWDelta:HGWBox);
        hgo[i].window_width=(vw.size()==0?0.:vw[i]);
        hgo[i].boundaries.resize((vn.size()==0?101:vn[i]+1));
        for (int k=0; k<hgo[i].boundaries.size();k++)
            hgo[i].boundaries[k]=va[i]+k*(vb[i]-va[i])/(hgo[i].boundaries.size()-1);
    }
    
    NDHistogram<double> HG(hgo);
    
    std::valarray<double> val(ndim); double weight;
    if (fweighted)
        while (std::cin.good()) { 
            for (int i=0; i<ndim; ++i) std::cin>>val[i]; std::cin>>weight;
            HG.add(val,weight); 
            }
    else
        while (std::cin.good()) { 
            for (int i=0; i<ndim; ++i) std::cin>>val[i];
            HG<<val; 
        }
    
    double outliers;
    HG.get_outliers(outliers);
    std::cout<<"# Fraction outside: "<<outliers<<std::endl;
    std::cout.setf(std::ios::scientific);
    if (asmooth=="") {
    if (!fgnu) std::cout <<HG; 
    else
    {
        
        if (ndim!=2) ERROR("GNUPLOT format works only for dimension 2\n");
        std::valarray<long> ind(2); std::valarray<double> cen(2); double val;
        for (int i=0; i<hgo[0].boundaries.size()-1; ++i) 
        {
            for (int j=0; j<hgo[1].boundaries.size()-1; ++j)
            {
                ind[0]=i; ind[1]=j;
                HG.get_bin(ind,cen,val);
                std::cout<<cen[0]<<"\t"<<cen[1]<<"\t"<<val<<"\n";
            }
            std::cout<<std::endl;
        }
    }
    } else {
      if (!fgnu) { ERROR("Multidimensional smoothing not implemented\n"); }
      else
      {
        
        if (ndim!=2) ERROR("GNUPLOT format works only for dimension 2\n");

        std::valarray<long> ind(2); std::valarray<double> cen(2), dcen; double val, asr, dr, tv, w, tw;
        double maxv=HG.max(), minv=HG.min();  bool breakout=false;
        asop.push_back(maxv); asop.push_back(minv);
            std::cerr<<"aops"<<asop[0]<<","<<asop[1]<<","<<asop[2]<<","<<asop[3]<<"\n";
        for (int i=0; i<hgo[0].boundaries.size()-1; ++i) 
        {
            for (int j=0; j<hgo[1].boundaries.size()-1; ++j)
            {
                ind[0]=i; ind[1]=j;
                HG.get_bin(ind,cen,val);
                asr=as_radius(asop,val);
                //std::cerr<<"Value at center is "<<val<<" & radius "<<asr<<"\n";
                tw=tv=0.0; breakout=false;
                for(int di=0; !breakout && di<=((i<hgo[0].boundaries.size()-i-2)?i:hgo[0].boundaries.size()-i-2); di++)
                   for(int dj=0; dj<=((j<hgo[1].boundaries.size()-j-2)?j:hgo[1].boundaries.size()-j-2); dj++)
                   {
                     ind[0]=i+di; ind[1]=j+dj;
                     HG.get_bin(ind,dcen,val);
                     dr=sqrt((dcen[0]-cen[0])*(dcen[0]-cen[0])+(dcen[1]-cen[1])*(dcen[1]-cen[1]));
                     if (dr>asr) { if (dj==0) breakout=true; break; }
                     w=as_kernel(dr/asr); tw+=w; tv+=val*w;
                     ind[0]=i-di; ind[1]=j+dj; HG.get_bin(ind,dcen,val); tw+=w; tv+=val*w;
                     ind[0]=i-di; ind[1]=j-dj; HG.get_bin(ind,dcen,val); tw+=w; tv+=val*w;
                     ind[0]=i+di; ind[1]=j-dj; HG.get_bin(ind,dcen,val); tw+=w; tv+=val*w;
                   }
                
                //std::cerr<<"Smoothened val is "<<tv/tw<<"\n";
                std::cout<<cen[0]<<"\t"<<cen[1]<<"\t"<<tv/tw<<"\n";
                
            }
            std::cout<<std::endl;
        }
      }
    }
}
