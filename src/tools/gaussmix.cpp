#include "tbdefs.hpp"
#include "matrix-full.hpp"
#include "matrix-io.hpp"
#include "rndgen.hpp"
#include "linalg.hpp"
#include "clparser.hpp"
#include <iostream>

using namespace toolbox;

class Gaussian {
private:
    int sz;
    std::valarray<double> m;
    FMatrix<double> C, iC; 
    double ldet;
    
public:
    Gaussian(int nsz=0): sz(nsz), C(nsz,nsz), m(nsz) {};
    
    void update(const std::valarray<double>& nm, const FMatrix<double>& nC)
    {
        if (nm.size()!=nC.rows() || nm.size()!=nC.cols()) ERROR("Size mismatch in given parameters");
        sz=nm.size();
        m.resize(sz); m=nm;
        C.resize(sz,sz); C=nC;
        //it would be much better to use cholesky instead. TODO 
        MatrixInverse(C,iC);
        FMatrix<double> Q; std::valarray<double> q;
        EigenSolverSym(C,Q,q);
        ldet=0.; for (int i=0; i<sz; ++i) ldet+=log(q[i]>0.?q[i]:1e-100);
    }
    
    void getpars(std::valarray<double>& nm, FMatrix<double>& nC) const
    {
        nm.resize(sz); nm=m; nC=C;
    }
    
    double lp(const std::valarray<double>& x) const
    {
        if (x.size()!=sz) ERROR("Size mismatch in given vector");
        std::valarray<double> y(x); double rv=0.;
        y-=m;
        for (int i=0; i<sz; ++i) for (int j=0; j<sz; ++j) rv+=iC(i,j)*y[i]*y[j];
        rv*=-0.5;
        rv+=-log(2.*constant::pi)*double(sz)/2.;
        rv+=-0.5*ldet;
        
        return rv;
    }
    
    double p(const std::valarray<double>& x) const
    { return exp(lp(x)); }
};

void banner()
{
    std::cerr
            << " USAGE: gaussmix -d dimension [-ng n-gauss]             \n"
            << "                  [-ns n-iter-opt] [-eps tolerance] [-s smoothing]             \n"
            << "                 -h < INPUT > OUTPUT                                      \n"
            << " performs a gaussian mixture clustering of a set of input data points.          \n";
}

int main(int argc, char** argv)
{
    CLParser clp(argc,argv);
    
    unsigned long kgm, mxstep, ne, nd, seed;
    //reads input presented on stdin in full-matrix form
    FMatrix<double> rdata;
    double eps, smooth;
    bool fhelp, fok=
            clp.getoption(ne,"d") &&
            clp.getoption(kgm,"ng",3ul) &&
            clp.getoption(mxstep,"ns",100ul) &&
            clp.getoption(eps,"e",1e-5) &&
            clp.getoption(smooth,"s",0.) &&
            clp.getoption(seed,"seed",12345ul) &&
            clp.getoption(fhelp,"h",false);

    if (fhelp) { banner(); return -1; }

    std::valarray<double> val(ne);
    std::vector<std::valarray<double> > vdata;
    while (std::cin.good()) {    
       for (int i=0; i<ne; ++i) std::cin>>val[i];
       vdata.push_back(val);
    }
    nd=vdata.size(); 
    
    std::valarray<Gaussian> gauss(kgm); std::valarray<double> lpg(kgm);
    //initializes "randomly" the gaussian data
    FMatrix<double> Ck(ne,ne); std::valarray<double> mk(ne);
    StdRndUniform rndgen(seed);
    for (int k=0; k<kgm; ++k)
    {
        Ck*=0.; 
        int icnt=rndgen()*nd;
        for (int i=0; i<ne; ++i) { Ck(i,i)=1e-5; mk[i]=vdata[icnt][i]; }
        gauss[k].update(mk,Ck);
        lpg[k]=-log(double(kgm));
    }
    
    double llike, lolike=-1e-100;
    FMatrix<double> lpnk(nd,kgm), pnk(nd,kgm);
    for (int is=0; is<mxstep; ++is)
    {
        //E step
        //first, gets unnormalized logs of pnk
        for (int i=0; i<nd; ++i)
        for (int k=0; k<kgm; ++k)
        { lpnk(i,k)=gauss[k].lp(vdata[i])+lpg[k]; }
        
        
        //then, normalizes logs
        llike=0.;
        for (int i=0; i<nd; ++i)
        { 
            double mxl, sume, lpnorm;
            mxl=lpnk(i,0); for (int k=1; k<kgm; ++k) if (lpnk(i,k)>mxl) mxl=lpnk(i,k);
            sume=0.; for (int k=0; k<kgm; ++k) sume+=exp(lpnk(i,k)-mxl);
            //std::cerr<<"mxk, sume"<<mxl<<":"<<sume<<"\n";
            lpnorm=mxl+log(sume);
            for (int k=0; k<kgm; ++k) { lpnk(i,k)-=lpnorm; pnk(i,k)=exp(lpnk(i,k)); }
            sume=0.; for (int k=0; k<kgm; ++k) sume+=pnk(i,k);
        //    std::cerr<<" tpn " <<sume<<"\n";
            llike+=lpnorm;
        }
       
        std::cerr<<"Log likelihood: "<<llike<<"\n";
        if (fabs((llike-lolike)/llike)<eps) break;
        //M step (updates gaussian parameters)
        for (int k=0; k<kgm; ++k)
        {
            double tpnk=0.; for (int i=0; i<nd; ++i) tpnk+=pnk(i,k);
            lpg[k]=log(tpnk/nd);
            mk=0.;
            for (int i=0; i<nd; ++i) mk+=vdata[i]*pnk(i,k);
            mk*=1./tpnk;
            for (int i=0; i<ne; ++i) for (int j=0; j<ne; ++j)
            { 
                Ck(i,j)=0.; 
                for (int h=0; h<nd; ++h) Ck(i,j)+=pnk(h,k)*(vdata[h][i]-mk[i])*(vdata[h][j]-mk[j]);
                
                Ck(i,j)*=1./tpnk;
            }
            double ctr=trace(Ck)/ne*smooth;
            Ck*=(1.-smooth);
            for (int i=0; i<ne; ++i) Ck(i,i)+=ctr;
            
            std::cerr<<k<<" "<<lpg[k]<<" "<<mk[0]<<" "<<Ck(0,0)<<"\n";
            gauss[k].update(mk,Ck);
        }
    }
    std::cout<<"log-likelihood: "<<llike<<"\n";
    for (int k=0; k<kgm; ++k)
    {
        gauss[k].getpars(mk,Ck);
        std::cout<<k<<" "<<lpg[k]<<" ";
        for (int i=0; i<ne; ++i) std::cout<<mk[i]<<" ";
        for (int i=0; i<ne; ++i) for (int j=0; j<ne; ++j) std::cout<<Ck(i,j)<<" ";
        std::cout<<"\n";
    }
    return 0;
}
