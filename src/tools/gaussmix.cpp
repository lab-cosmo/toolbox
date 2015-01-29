#include "tbdefs.hpp"
#include "matrix-full.hpp"
#include "matrix-io.hpp"
#include "rndgen.hpp"
#include "linalg.hpp"
#include "clparser.hpp"
#include <iostream>
#include <fstream>

using namespace toolbox;

/* A simple class to compute a multivariate Gaussian probability distribution */

class Gaussian {
private:
    int sz;                   // the dimensionality of the Gaussian
    std::valarray<double> m;  // the mean
    FMatrix<double> C, iC;    // the covariance matrix and its inverse
    double ldet;              // logarithm of the determinant of C (for the normalization)
    double lpi;               // extra constant entering the normalization
public:
    // constructor
    Gaussian(int nsz=0): sz(nsz), C(nsz,nsz), m(nsz) {};


    // sets the parameters for the evaluation of the Gaussian PDF
    void update(const std::valarray<double>& nm, const FMatrix<double>& nC)
    {
        // resizes the private arrays
        if (nm.size()!=nC.rows() || nm.size()!=nC.cols()) ERROR("Size mismatch in given parameters");
        sz=nm.size();
        m.resize(sz); m=nm;
        C.resize(sz,sz); C=nC;


        // computes intermediate values that are handy to compute the PDF
        // inverse covariance matrix
        MatrixInverse(C,iC);  // it would be much better to use cholesky instead. TODO
        FMatrix<double> Q; std::valarray<double> q; 
        EigenSolverSym(C,Q,q);  // computes the determinant the silly way
        ldet=0.; for (int i=0; i<sz; ++i) ldet+=log(q[i]>0.?q[i]:1e-100);
        lpi=-log(2.*constant::pi)*double(sz)*0.5;  // the other constant in the normalization
    }

    // gets the parameters of the Gaussian PDF
    void getpars(std::valarray<double>& nm, FMatrix<double>& nC) const
    {
        nm.resize(sz); nm=m; nC=C;
    }

    // computes the log of the PDF evaluated at x
    double lp(const std::valarray<double>& x) const
    {
        if (x.size()!=sz) ERROR("Size mismatch in given vector");
        std::valarray<double> y(x); double rv=0.;
        y-=m;
        for (int i=0; i<sz; ++i) for (int j=0; j<sz; ++j) rv+=iC(i,j)*y[i]*y[j];
        rv*=-0.5;
        rv+=lpi;
        rv+=-0.5*ldet;

        return rv;
    }

    // computes the PDF evaluated at x
    double p(const std::valarray<double>& x) const
    { return exp(lp(x)); }
};

void banner()
{
    std::cerr
            << " USAGE: gaussmix -d dimension [-ng n-gauss]             \n"
            << "                  [-ns n-iter-opt] [-eps tolerance] [-s smoothing] [-init file] \n"
            << "                  [ -post file ] [-v] [-h] < INPUT > OUTPUT                                    \n"
            << " performs a gaussian mixture clustering of a set of input data points,          \n"
            << " given on stdin as x11 x12 x13 ... x1d x21 x22 ... x2d                         \n"
            << " -d specifies the dimensionality of the input points                            \n"
            << " -ng specifies the number of Gaussian clusters in the mixture   [3]            \n"
            << " -ns specifies the maximum number of E-M optimization steps  [100]             \n"
            << " -eps gives a tolerance on the change in log-likelihood to stop the optimization [1e-5]\n"
            << " -s gives a regularization parameter for the covariance matrices [0] \n"
            << " -init initializes clustering from the given file (compact format):  \n"
            << "       m1 m2 ... md c11 c12 .. c1d ... cd1 .. cdd lprobk ...\n"
            << " -v prints a verbose, more readable report on clustering \n"
            << " -post analizes the points in input giving the fuzzy clustering prob. based on the \n"
            << "        cluster information contained in the specified file. outputs \n"
            << "        p11 p12 ... p1k \n p21 p22 ... p2k ....\n"
            ;
}

int main(int argc, char** argv)
{
    CLParser clp(argc,argv);

    unsigned long kgm, mxstep, ne, nd, seed;
    //reads input presented on stdin in full-matrix form
    FMatrix<double> rdata; std::string filein, filepost;
    double eps, smooth;
    bool fhelp = false, fverb;
    bool fok=
            clp.getoption(ne,"d") &&
            clp.getoption(kgm,"ng",3ul) &&
            clp.getoption(mxstep,"ns",100ul) &&
            clp.getoption(eps,"e",1e-5) &&
            clp.getoption(smooth,"s",0.) &&
            clp.getoption(seed,"seed",12345ul) &&
            clp.getoption(filein,"init",std::string("")) &&
            clp.getoption(filepost,"post",std::string("")) &&
            clp.getoption(fverb,"v",false) &&
            clp.getoption(fhelp,"h",false);

    if (fhelp or !fok)
    { banner(); return -1; }

    std::valarray<Gaussian> gauss(kgm);     // makes kgm Gaussian classes (to compute the cluster probabilities)
    std::valarray<double> lpg(kgm);         // these are the log-weight of
    FMatrix<double> Ck(ne,ne); std::valarray<double> mk(ne);
    std::valarray<double> val(ne);

    /* post-processing mode */
    // restarts from file
    if (filepost!="")
    {
        std::ifstream fpp(filepost.c_str());
        for (unsigned long i=0; i<kgm; i++)
        {
            for (unsigned long j=0; j<ne; j++) fpp>>mk[j];
            for (unsigned long j=0; j<ne; j++) for (unsigned long k=0; k<ne; k++) fpp>>Ck(j,k);
            fpp>>lpg[i];  // not used
            gauss[i].update(mk,Ck);
        }

        double tpk; std::valarray<double> pk(kgm);
        while (std::cin.good()) {
            for (int i=0; i<ne; ++i) std::cin>>val[i];
            for (int k=0; k<kgm; ++k) pk[k]=exp(gauss[k].lp(val)+lpg[k]);
            tpk = pk.sum();
            for (int k=0; k<kgm; ++k) std::cout<<pk[k]/tpk<<" ";
            std::cout << std::endl;
        }

        return 0;
    }


    /* GM analysis mode */
    // reads data points from stdin and stores them in vdata
    std::vector< std::valarray<double> > vdata;
    while (std::cin.good()) {
       for (int i=0; i<ne; ++i) std::cin>>val[i];
       vdata.push_back(val);
    }
    nd=vdata.size();  // number of data points


    //initializes "randomly" the gaussian data
    StdRndUniform rndgen(seed);
    for (int k=0; k<kgm; ++k)
    {
        Ck*=0.;
        int icnt=rndgen()*nd;
        for (int i=0; i<ne; ++i) { Ck(i,i)=1e-5; mk[i]=vdata[icnt][i]; }
        gauss[k].update(mk,Ck);
        lpg[k]=-log(double(kgm));   // each cluster has a weight 1/kgm
    }


    // restarts from file
    if (filein!="") //!TODO CLEAN UP MAKING CONSISTENT WITH THE OUTPUT
    {
        std::ifstream fip(filein.c_str());
        for (unsigned long i=0; i<kgm; i++)
        {
            for (unsigned long j=0; j<ne; j++) fip>>mk[j];
            for (unsigned long j=0; j<ne; j++) for (unsigned long k=0; k<ne; k++) fip>>Ck(j,k);
            fip>>lpg[i];
            gauss[i].update(mk,Ck);
        }
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
            // finds the largest Gaussian weight for point i
            mxl=lpnk(i,0); for (int k=1; k<kgm; ++k) if (lpnk(i,k)>mxl) mxl=lpnk(i,k);
            sume=0.; for (int k=0; k<kgm; ++k) sume+=exp(lpnk(i,k)-mxl);  // computes the total probability relative to the maximum weight, for stability
            lpnorm=mxl+log(sume);

            for (int k=0; k<kgm; ++k) { lpnk(i,k)-= lpnorm; pnk(i,k)=exp(lpnk(i,k)); }
            sume=0.; for (int k=0; k<kgm; ++k) sume+=pnk(i,k);

            llike+=lpnorm;
        }

        // here we have the total log-likelyhood
        std::cerr<<"Log likelihood: "<<llike<<"\n";
        if (fabs((llike-lolike)/llike)<eps) break;   // exits if the log-likelyhood has changed very little

        //M step (updates gaussian parameters)
        for (int k=0; k<kgm; ++k)
        {
            double tpnk=0.; for (int i=0; i<nd; ++i) tpnk+=pnk(i,k);   // normalization for cluster k
            lpg[k]=log(tpnk/nd);

            // compute mean for cluster k
            mk=0.;
            for (int i=0; i<nd; ++i) mk+=vdata[i]*pnk(i,k);
            mk*=1./tpnk;

            // compute variance for cluster k
            for (int i=0; i<ne; ++i) for (int j=0; j<ne; ++j)
            {
                Ck(i,j)=0.;
                for (int h=0; h<nd; ++h) Ck(i,j)+=pnk(h,k)*(vdata[h][i]-mk[i])*(vdata[h][j]-mk[j]);

                Ck(i,j)*=1./tpnk;
            }

            // stabilizes the covariance with a smoothing
            double ctr=trace(Ck)/ne*smooth;
            Ck*=(1.-smooth);
            for (int i=0; i<ne; ++i) Ck(i,i)+=ctr;

            std::cerr<<k<<" "<<lpg[k]<<" "<<mk[0]<<" "<<Ck(0,0)<<"\n";

            // updates the cluster
            gauss[k].update(mk,Ck);
        }
    }

    std::cout<<std::scientific; std::cout.precision(10);
    if (fverb)
    {
       std::cout<<"# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # \n";
       std::cout<<"# Gaussian mixture clustering report\n# Log-likelihood: "<<llike<<"\n";
       std::cout<<"# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # \n";
       for (int k=0; k<kgm; ++k)
       {
           std::cout<<"# * Gaussian cluster n. "<<k<<"  log-like: "<<lpg[k]<<std::endl;
           gauss[k].getpars(mk,Ck);
           std::cout<<"# Mean: \n";
           for (int i=0; i<ne; ++i) std::cout<<mk[i]<<" ";
           std::cout<<"\n# Covariance: \n";
           for (int i=0; i<ne; ++i) { for (int j=0; j<ne; ++j) std::cout<<Ck(i,j)<<" "; std::cout<<std::endl; }
       }
   }

   std::cout.precision(8); std::cout.width(15); std::cout.setf(std::ios::scientific);
   //compact output
   std::cout<<"# MEAN 1 2 ...  COVARIANCE 11 12 ... 21 22 ....   WEIGHT"<<std::endl;
   for (int k=0; k<kgm; ++k)
   {
      gauss[k].getpars(mk,Ck);
      for (int i=0; i<ne; ++i) std::cout<<mk[i]<<" ";
      for (int i=0; i<ne; ++i) for (int j=0; j<ne; ++j) std::cout<<Ck(i,j)<<" ";
      std::cout<<lpg[k]<<" "<<std::endl;
   }

   return 0;
}
