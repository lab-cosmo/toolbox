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
            << "                                                                                \n"
            << " compute the histogram of a series of data, in ndim dimensions.                 \n"
            << " data must be formatted as                                                      \n"
            << " D1(1)   D2(1)  .......   DN(1) [ WEIGHT(1) ]                                   \n"
            << " on every dimension n (default:100) bins are distributed evenly.                \n"
            << " between xi and xf. optionally, box (b) smoothing functions can be used.        \n"
            << " If dimension is 2 and -g is selected, points will be output in gnuplot format. \n"
            << " If -w is selected, then a weight is read after each point.                     \n";
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

int main(int argc, char **argv)
{
    
    CLParser clp(argc, argv);
    bool fhelp, fweighted, fgnu;
    std::string a,b,wb,nbins;
    unsigned long ndim;
    bool fok=
            clp.getoption(ndim,"d") &&
            clp.getoption(a,"xi") &&
            clp.getoption(b,"xf") &&
            clp.getoption(wb,"b",std::string("")) &&
            clp.getoption(nbins,"n",std::string("")) &&
            clp.getoption(fweighted,"w",false) &&
            clp.getoption(fgnu,"g",false) &&
            clp.getoption(fhelp,"h",false);
    
    if ( fhelp || ! fok) { banner(); return 0; }
    
    std::valarray<HGOptions<Histogram<double> > >hgo(ndim);
    std::vector<double> va, vb, vw, vn;
    ssplit(a,va); ssplit(b,vb); ssplit(wb,vw); ssplit(nbins,vn); 
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
}
