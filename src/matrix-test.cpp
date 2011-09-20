#define TB_MPI 1 
#include "tbdefs.hpp"
#include "matrix-pcrs.hpp"
#include "matrix-full.hpp"
#include "matrix-crs.hpp"
#include "matrix-coord.hpp"
#include "matrix-conv.hpp"
#include "matrix-fun.hpp"
#include "matrix-io.hpp"
#include "rndgen.hpp"
#include "ioparser.hpp"
#include <valarray> 
using namespace toolbox;

typedef FMatrix<double> MType;
RndGaussian<double, MTRndUniform> rgg(RGPars<double>(0,1));

void fill_ckr(unsigned long n, unsigned long m, unsigned long s, MType& a)
{
    a.resize(n,m);
    for (unsigned long i=0; i<n; i+=s)
    for (unsigned long j=0; j<m; j+=s)
    {
        a(i,j)=i+j;
    }
}

void fill_rnd(unsigned long n, unsigned long m, unsigned long s, MType& a)
{
    unsigned long wi, ts;
    std::valarray<std::valarray<unsigned long> > iind(n);
    std::valarray<std::valarray<double> > ival(n);
    for (unsigned long i=0; i<n; i++)
    {
        wi=(unsigned long) (rgg.RUGenerator()()*(s+1));
        iind[i].resize(wi); ival[i].resize(wi);
        for (unsigned long j=0; j<wi; j++)
        {
            iind[i][j]=(unsigned long) (rgg.RUGenerator()()*m);
            for (unsigned long k=0; k<j;++k) 
                if(iind[i][j]==iind[i][k]) { iind[i][j]=(unsigned long) (rgg.RUGenerator()()*s); k=-1; }
            ival[i][j]=rgg();
        }
        for (unsigned long j=0; j<wi; j++)
        {
            for (unsigned long k=j+1; k<wi;++k) 
                if (iind[i][j]>iind[i][k]) { 
                    ts=iind[i][j]; iind[i][j]=iind[i][k]; iind[i][k]=ts;}
        }
    }
//    a=MType(n,m,iind,ival);
}

void fill_diag(unsigned long n, unsigned long m,  MType& a)
{
    a.resize(n,m);
    for (unsigned long i=0; i<(n>m?m:n); i++)
    {
        a(i,i)=i;
        a(i,0)=i;
        a(0,i)=i;
    }
}

void fill_bnd(unsigned long n, unsigned long m,  MType& a)
{
    a.resize(n,m);
    for (unsigned long i=0; i<(n>m?m:n); i++)
        for (unsigned long j=0; i<(n>m?m:n); i++)
        {
            a(i,j)=rgg.RUGenerator()()*(exp(-(i-j)*0.2*(i-j))+exp(-(i-j-n/2)*0.2*(i-j-n/2)));
        }
    MType b;
    transpose(a,b);
    a+=b;
}


#define DOTEST \
    hrt.start(); poly_standard(A,cf,C1); hrt.stop();\
    E.setops(MatrixOptions<MType>(0.)); E.resize(n,n); \
    std::cout<< "NORM OF 'standard' C: "<< normfrob(C1)<<" \n"; \
    E*=0.; incr(E,C); incr(E,C1,-1.); \
    std::cout<< "ERROR FOR 'standard' C: "<< normfrob(E)<<" \n"; \
    std::cout<< "TIME FOR 'standard' C: "<< hrt*1e-6<<" MTicks\n"; \
    \
    hrt.start(); poly_liang(A,cf,C1); hrt.stop();\
    E.setops(MatrixOptions<MType>(0.)); E.resize(n,n); \
    std::cout<< "NORM OF 'liang' C: "<< normfrob(C1)<<" \n"; \
    E*=0.; incr(E,C); incr(E,C1,-1.); \
    std::cout<< "ERROR FOR 'liang' C: "<< normfrob(E)<<" \n"; \
    std::cout<< "TIME FOR 'liang' C: "<< hrt*1e-6<<" MTicks\n"; \
    \
    hrt.start(); poly_vanloan1(A,cf,C1); hrt.stop();\
    E.setops(MatrixOptions<MType>(0.)); E.resize(n,n); \
    std::cout<< "NORM OF 'vanloan1' C: "<< normfrob(C1)<<" \n"; \
    E*=0.; incr(E,C); incr(E,C1,-1.); \
    std::cout<< "ERROR FOR 'vanloan1' C: "<< normfrob(E)<<" \n"; \
    std::cout<< "TIME FOR 'vanloan1' C: "<< hrt*1e-6<<" MTicks\n"; 
#include <fstream>

 
 
void printc(const double& cc)
{
    IOMap iom;
    iom.insert(cc,"test");
    std::cout <<iom;
}
int main(int argc, char **argv)
{

    MPI_Init( &argc, &argv );
    int myrank, mysize;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    RndGaussian<double,MTRndUniform> rgg(RGPars<double>(0.,1.)); 
    //rgg.RUGenerator().seed(time(NULL));
    rgg.RUGenerator().seed(332);
    HRTimer hrt;
    
 /* FMatrix<double> x1, x2, x3;
    CrsMatrix<double> y1, y2, y3;
    fill_rnd(10,10,6,y1);
    fill_rnd(10,10,6,y2);
    x1=y1; x2=y2;
    
    x1+=1.; y1+=1.;
    x3=x1; x3+=x2;
    y3=y1; y3+=y2;
    mult(x1,x2,x3);
    mult(y1,y2,y3);
    
    
    
    std::cerr<<"FULL \n"<<std::string(x3);
    std::cerr<<"CRS   \n"<<std::string(y3);
    std::cerr<<"FULL  "<<norm1(x3)<<" "<<norminf(x3)<<" "<<normfrob(x3)<<"\n";
    std::cerr<<"CRS   "<<norm1(y3)<<" "<<norminf(y3)<<" "<<normfrob(y3)<<"\n";
    return 0;
   */
    //1500 1500 10
    unsigned long n=10, m=10, p=8;
    
    std::cout<<" ******** TESTING STANDARD OPER. *********** \n";
    MType M1, M2, M3;
    
    //fill_rnd(n,m,20,M1);
    //fill_rnd(n,m,20,M2);
    M1.resize(n,m);
    M2.resize(m,n);
    
    for (unsigned long k=0; k<n*p; ++k)
    {
        long i=rgg.RUGenerator()()*n;
        long j=i+(rgg.RUGenerator()()-0.5)*p;
        if (j<0) j=i; 
        if (j>=m) j=i;
        if (j>=m) j=0;
        M1(i,j)+=k;
        M2(j,i)+=1;
    }
    //M1.optimize(); M2.optimize();
    transpose(M2,M3);
    std::valarray<double> pp(6);
    pp[0]=1.; pp[1]=0.5; pp[2]=0.6; pp[3]=-0.6; pp[4]=-0.5; pp[5]=-1.;
    
    M1.resize(1,1);
    M1(0,0)=-1.2;
    std::cerr<<M1<<"<< WAS M1\n";
    chebyshev_poly(M1,pp,M2);
    std::cerr<<M2<<"<< WAS M2\n";
    chebyshev_fastpoly(M1,pp,M3);
    std::cerr<<M3<<"<< WAS M3\n";
    
    /*M1.resize(n,n);
    for (int i=0; i<n; ++i)
    {
        M1(i,i)=i;
        M1(i,0)=i;
        M1(0,i)=i;
    }*/
    std::cout.flush();  std::cerr.flush();
    exit(0);
    MPI_Barrier(MPI_COMM_WORLD);
    /*!!
    for (int mn=0; mn<mysize; mn++) 
    {
        if (myrank==mn) 
        {
            std::cerr<<" I am on node: "<<mn<<"\n";
            std::cerr<<" which holds rows "<<M2.nroots[mn]<<" - "<<M2.nroots[mn+1]-1<<"\n";
            for (unsigned int i=0; i<m; ++i)
            {
                for (unsigned int j=0; j<n; ++j)
                    std::cerr<<M2(i,j)<<" ";
                std::cerr<<"\n";
            }
            std::cerr<<"************************************8\n";
        }
        std::cerr.flush();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    //*/
    //M3+=1;
    MPI_Barrier(MPI_COMM_WORLD);
    /*!!
    for (int mn=0; mn<mysize; mn++) 
    {
        if (myrank==mn) 
        {
            std::cerr<<" I am on node: "<<mn<<"\n";
            std::cerr<<" which holds rows "<<M3.nroots[mn]<<" - "<<M3.nroots[mn+1]-1<<"\n";
            for (unsigned int i=0; i<n; ++i)
            {
                for (unsigned int j=0; j<m; ++j)
                    std::cerr<<M3(i,j)<<" ";
                std::cerr<<"\n";
            }
            std::cerr<<"************************************8\n";
        }
        std::cerr.flush();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    //*/
    
    M3=M1;
    
    /*
    if (myrank==0) {
    std::cerr<<"\nmatrix trace is: "<<ntr<<"\n";
    std::cerr<<"\nmatrix norm  is: "<<nnr<<"\n";
    }
    M1+=1.;
    
    ntr=trace(M1); nnr=norminf(M1);
    if (myrank==0) {
        std::cerr<<"\nmatrix trace is: "<<ntr<<"\n";
        std::cerr<<"\nmatrix norm  is: "<<nnr<<"\n";
    }
    */
    
    /*
    std::cerr<<M2.size()<<"\n";
    */
    MPI_Barrier(MPI_COMM_WORLD);
    
    M3.resize(0,0);
    //M1+=5;
    IterOptions<double,3> inwops(100,1e-3, 0., ichk_default);
    inwops.thresh[0]=1.; //MINIMUM we want the residual to be below 1!
    inwops.flags[1]=ichk_sufficient;
    inwops.flags[2]=ichk_sufficient | ichk_change;
    //inverse_newton(M1,M3,inwops);
    MatrixOptions<MType> TO(0.);
    M2=M1;
    M3.setops(TO);
    
    double tass;
    hrt.start();
    for (unsigned long i=0; i<100; ++i) { M3=M2; }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix assign: "<<(tass=(hrt/100.*1e-6))<<"  MTicks\n";

    M3=M2; M3.setops(TO);
    hrt.start();
    for (unsigned long i=0; i<100; ++i) { M3*=(double) i; }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix scale: "<<(hrt/100.*1e-6)<<"  MTicks\n";
    
    hrt.start();
    for (unsigned long i=0; i<100; ++i) { M3=M1; M3.setops(TO); incr(M3,M2); }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix increment: "<<(hrt/100.*1e-6-tass)<<"  MTicks\n";
    
    hrt.start();
    for (unsigned long i=0; i<100; ++i) { M3=M1; M3.setops(TO); incr(M3,M1); }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix increment, same pattern: "<<(hrt/100.*1e-6-tass)<<"  MTicks\n";
    
    M3.setops(TO);
    hrt.start();
    for (unsigned long i=0; i<100; ++i) { add(M1,M2,M3); }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix-matrix add: "<<(hrt/100.*1e-6)<<"  MTicks\n";
    
    M3.setops(TO);
    hrt.start();
    for (unsigned long i=0; i<100; ++i) { add(M1,M1,M3); }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix-matrix add, same pattern: "<<(hrt/100.*1e-6)<<"  MTicks\n";
    
    hrt.start();
    for (unsigned long i=0; i<1000; ++i) { M3=M1; M3.setops(TO); incr(M3,(double) i); }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix diagonal +: "<<(hrt/1000.*1e-6-tass)<<"  MTicks\n";

    M3=M1; M3.setops(TO); incr(M3,(double) 1);
    hrt.start();
    for (unsigned long i=0; i<1000; ++i) { incr(M3,(double) i);}
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix diagonal +, filled: "<<(hrt/1000.*1e-6)<<"  MTicks\n";
    
    hrt.start();
    for (unsigned long i=0; i<100; ++i) {  M3=M1; M3.setops(TO);
        for (unsigned long k=0; k<M3.rows();++k) M3(k,k)+=i;
    }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix diagonal + (k,k): "<<(hrt/100.*1e-6-tass)<<"  MTicks\n";
    
    M3=M1; M3.setops(TO); incr(M3,(double) 1);
    hrt.start();
    for (unsigned long i=0; i<1000; ++i) {
        for (unsigned long k=0; k<M3.rows();++k) M3(k,k)+=i;
    }
    hrt.stop();
    if (myrank==0) std::cout<<" Timing for matrix diagonal + (k,k), filled: "<<(hrt/1000.*1e-6)<<"  MTicks\n";
    
    mult:
    M3.setops(TO);
    hrt.start();
    //M2*=0.; M2+=2.;
    M1.trunc(0.1); M2.trunc(0.1);
    
    for (unsigned long i=0; i<1; ++i)
    { mult(M1,M2,M3); }
    hrt.stop();
    
    double dt1, dt2, dtr;
    double nt1, nt2, ntr;
    dt1=trace(M1); dt2=trace(M2); dtr=trace(M3);
    nt1=norminf(M1); nt2=norminf(M2); ntr=norminf(M3);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank==0)
    {
        std::cerr<<"\n\n******************************\n";
        std::cerr<<"trace 1: "<<dt1<<"\n"
                <<"trace 2: "<<dt2<<"\n" 
                <<"trace res: "<<dtr<<"\n";
        std::cerr<<"norm 1: "<<nt1<<"\n"
                <<"norm 2: "<<nt2<<"\n"
                <<"norm res: "<<ntr<<"\n";
        std::cerr<<"******************************\n";
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (myrank==0) std::cout<<" Timing for matrix-matrix multiply: "<<(hrt/1.*1e-6)<<" MTicks\n";
    
    MPI_Finalize();
    exit(0); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    
    
    MType A, AT, C, C1, C2, C3, E;
    fill_rnd(n,m,40,A); 
    //fill_ckr(n,m,1,A);
    transpose(A,AT);
    mult(A,AT,C);
    A=C*-0.001;
    A+=1;
    A.setops(MatrixOptions<MType>(0.));
    C.setops(A.getops());

    //std::cerr<<"MAT DIFFERENCE "<<normfrob(A)<<"\n";
    //std::cerr<<std::string(A)<<"\n";
    
    std::valarray<double> cf(p);
    double f=1;
    
    // COEFFICIENTS FOR EXPONENTIAL
    std::cout<<" ******** TESTING FOR EXPONENTIAL ********** \n";
    
    std::cout<<" computing result to machine precision...\n";
    C.resize(0,0);
    C.setops(MatrixOptions<MType>(0.));
    exp(A,C,1e-15);
    std::cout<< "NORM OF 'exact' C: "<< normfrob(C)<<" \n"; 
    
    for (unsigned long i=0; i<p; ++i) { cf[i]=f; f*=1./(i+1.); }
    
    std::cout<<" computing with 1e-7 absolute truncation\n";
    C1.setops(MatrixOptions<MType>(1e-7));
    C2=C1; C3=C1;
    
    DOTEST
    
    std::cout<<" computing with 1e-3 absolute truncation\n";
    C1.setops(MatrixOptions<MType>(1e-3));
    C2=C1; C3=C1;
    
    DOTEST
    
    
    // COEFFICIENTS FOR EXPONENTIAL
    std::cout<<" ********  TESTING FOR INVERSE  ********** \n";
    
    std::cout<<" computing result to machine precision...\n";
    C.resize(0,0);
    C.setops(MatrixOptions<MType>(0.));
    
    incr(A,1);
    inverse_newton(A,C,IterOptions<double,1>(100,1e-5,0.,ichk_change));
    std::cout<< "NORM OF 'exact' C: "<< normfrob(C)<<" \n"; 
    incr(A,-1);
    
    for (unsigned long i=0; i<p; ++i) { cf[i]=(i%2==0?1:-1); }
    
    std::cout<<" computing with 1e-7 absolute truncation\n";
    C1.setops(MatrixOptions<MType>(1e-7));
    C2=C1; C3=C1;
    
    DOTEST
    
    std::cout<<" computing with 1e-3 absolute truncation\n";
    C1.setops(MatrixOptions<MType>(1e-3));
    C2=C1; C3=C1;
    
    DOTEST
    
    std::cout<<" computing with 1e-5 relative truncation\n";
    C1.setops(MatrixOptions<MType>(1e-5, 0., true, at_norminf));
    C2=C1; C3=C1;
    
    DOTEST
    
    std::cout<<" computing with 1e-3 relative truncation\n";
    C1.setops(MatrixOptions<MType>(1e-3, 0., true, at_norminf));
    C2=C1; C3=C1;
    
    DOTEST;
    
#ifdef BENCHMARK
    TBBenchmarks.print(std::cerr);
#endif
    return -1;
}
