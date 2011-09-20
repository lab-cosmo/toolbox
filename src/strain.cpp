#include "tbdefs.hpp"
#include "ioparser.hpp"
#include "clparser.hpp"
#include "matrix-full.hpp"
#include <string>
#include <vector>
extern "C" {
    void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
};

class pos {
    public:
        double x, y, z;
    
        pos() : x(0.), y(0.), z(0.) {}
    
        pos& operator += (const pos& p2)
        { x+=p2.x; y+=p2.y; z+=p2.z; return *this;}
    
        pos& operator -= (const pos& p2)
        { x-=p2.x; y-=p2.y; z-=p2.z; return *this;}
        
        pos& operator *= (const double& s)
        { x*=s; y*=s; z*=s; return *this;}
};
pos operator+(const pos& pa, const pos& pb)
{ pos pc(pa); pc+=pb; return pc; }
pos operator-(const pos& pa, const pos& pb)
{ pos pc(pa); pc-=pb; return pc; }
double dot(const pos& a, const pos& b)
{ return a.x*b.x+a.y*b.y+a.z*b.z; }
pos cross(const pos& a, const pos& b)
{
    pos c;
    c.z=a.x*b.y-a.y*b.x;
    c.y=a.z*b.x-a.x*b.z;
    c.x=a.y*b.z-a.z*b.y;
    return c;
}

using namespace toolbox;
toolbox::HRTimer hrt;

typedef struct _atom {
    std::string s;
    pos r;
} atom;

class frame {
public:
    int nat;
    std::string comment;
    std::valarray<atom> ats;
    
    frame() {}
    frame(const frame& nf) : nat(nf.nat), comment(nf.comment), ats(nf.ats) {}
    frame& operator=(const frame& nf) 
    {
        if (&nf==this) return *this;
        nat=nf.nat; comment=nf.comment; 
        ats.resize(nf.ats.size()); ats=nf.ats;
        return *this;
    }
}; 

std::istream& operator >> (std::istream& is, atom& rat)
{
    is >> rat.s >> rat.r.x >> rat.r.y >>rat.r.z;
    return is;
}

std::istream& operator >> (std::istream& is, frame& fr)
{
    char line[1024];
    is >> fr.nat;
    is.getline(line,1024,'\n');
    is.getline(line,1024,'\n');
    fr.comment=line;
    fr.ats.resize(fr.nat);
    for (int i=0; i<fr.nat; ++i)
        is>>fr.ats[i];
    return is;
}

std::istream& operator >> (std::istream& is, std::vector<frame>& fv)
{
    frame nfr;
    while (is.good())
    {
        is>>nfr;
        fv.push_back(nfr);
    }
}

std::ostream& operator << (std::ostream& os, const frame& fr)
{ return os; } //do-nothing lazy stub


void rmsd(const std::valarray<pos>& r, const std::valarray<pos>& p, FMatrix<double>& RMAT)
{
    FMatrix<double> D(4,4); D*=0.;
    double q[4],P[6][6],U[6];
    
    
    for (long at=0; at<r.size(); ++at)
    {
        U[0]=p[at].x+r[at].x; U[1]=p[at].y+r[at].y; U[2]=p[at].z+r[at].z;
        U[3]=r[at].x-p[at].x; U[4]=r[at].y-p[at].y; U[5]=r[at].z-p[at].z;
    
        for (int i=0; i<6; ++i)
        {
            P[i][i]=U[i]*U[i];
            for (int j=0; j<i; ++j)
                P[i][j]=P[j][i]=U[i]*U[j];
        }
        
        double tp=P[0][0]+P[1][1]+P[2][2];
        double tm=P[3][3]+P[4][4]+P[5][5];
        D(0,0)+=tm; 
        D(1,1)+=tp-P[0][0]+P[3][3]; 
        D(2,2)+=tp-P[1][1]+P[4][4]; 
        D(3,3)+=tp-P[2][2]+P[5][5];
        D(0,1)+=P[1][5]-P[2][4];
        D(0,2)+=P[2][3]-P[0][5];
        D(0,3)+=P[0][4]-P[1][3];
        D(1,2)+=P[3][4]-P[0][1];
        D(1,3)+=P[3][5]-P[0][2];
        D(2,3)+=P[4][5]-P[1][2];
    }
    
    FMatrix<double>md(4,4), mq(4,4), rr;
    for (int i=0; i<4; ++i)
        for (int j=0; j<=i; ++j)
            md(i,j)=md(j,i)=D(j,i);
    
    //diagonalize D using LAPACK....
    char jobz='V', uplo='L';
    int N=4; int lda=N; 
    double l[4];
    int lwork=100; double work[100];
    int info;

    
    dsyev_(&jobz, &uplo, &N, &(D(0,0)), &lda, l, work, &lwork, &info);
    if (info!=0)
        ERROR("Error "<<info<<"in DSYEV");
    
    for (int i=0; i<4; ++i) q[i]=D(0,i);
    /*
    std::cerr<<"\nl= ";    
    for (unsigned long i=0; i<4; ++i)
        std::cerr<<l[i]<<" ";
    std::cerr<<"\nq= ";
    for (int i=0; i<4; ++i) std::cerr<<q[i]<<" ";
    std::cerr<<"\n";
    */
    
    RMAT.resize(3,3);
    RMAT(0,0)=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
    RMAT(1,1)=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
    RMAT(2,2)=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
    RMAT(0,1)=2*(q[1]*q[2]+q[0]*q[3]);
    RMAT(0,2)=2*(q[1]*q[3]-q[0]*q[2]);
    RMAT(1,2)=2*(q[2]*q[3]+q[0]*q[1]);
    RMAT(1,0)=2*(q[1]*q[2]-q[0]*q[3]);
    RMAT(2,0)=2*(q[1]*q[3]+q[0]*q[2]);
    RMAT(2,1)=2*(q[2]*q[3]-q[0]*q[1]);
    
    /*
    std::cerr.precision(5); std::cerr.setf(std::ios::scientific);
    std::cerr<<"REFERENCE: \n";
    for (int i=0; i<3; ++i)
        std::cerr<<r[i].x<<","<<r[i].y<<","<<r[i].z<<"\n";
    std::cerr<<"VECTORS: \n";
    for (int i=0; i<3; ++i)
        std::cerr<<p[i].x<<","<<p[i].y<<","<<p[i].z<<"\n";
            
    FMatrix<double> lr(3,3), r1;
    for (int i=0; i<3; ++i)
    { lr(0,i)=p[i].x; lr(1,i)=p[i].y; lr(2,i)=p[i].z;}
    
    mult(RMAT,lr,r1);
    std::cerr<<"NEW: \n";
    for (int i=0; i<3; ++i)
        std::cerr<<r1(0,i)<<","<<r1(1,i)<<","<<r1(2,i)<<"\n";
    std::cerr<<"OMAT: \n";
    for (int i=0; i<3; ++i)
        std::cerr<<RMAT(0,i)<<","<<RMAT(1,i)<<","<<RMAT(2,i)<<"\n";*/
/*    np.resize(p.size());
    for (long at=0; at<r.size(); ++at)
    {
        np[at].x=R[0][0]*p[at].x+R[0][1]*p[at].y+R[0][2]*p[at].z;
        np[at].y=R[1][0]*p[at].x+R[1][1]*p[at].y+R[1][2]*p[at].z;
        np[at].z=R[2][0]*p[at].x+R[2][1]*p[at].y+R[2][2]*p[at].z;
    }
*/
}

void strain(const frame& al, FMatrix<double>& rstrain)
{
    //id nij, nji;
    double a0=3.8254483, d0=0.612372436*a0, dc=d0*1.1, dw=dc*1.01, dc2=dc*dc;
    std::valarray<std::vector<int > > neighs;
    neighs.resize(al.nat);
    rstrain.resize(al.nat,6);
    
    hrt.start();
    /* prepares neighbour list */
    /* linked cells-like approach */
    pos l=al.ats[0].r, r=al.ats[0].r;
    for(unsigned long i=0; i<al.nat; ++i)
    { 
        if (l.x>al.ats[i].r.x) l.x=al.ats[i].r.x;
        if (l.y>al.ats[i].r.y) l.y=al.ats[i].r.y;
        if (l.z>al.ats[i].r.z) l.z=al.ats[i].r.z;
        
        if (r.x<al.ats[i].r.x) r.x=al.ats[i].r.x;
        if (r.y<al.ats[i].r.y) r.y=al.ats[i].r.y;
        if (r.z<al.ats[i].r.z) r.z=al.ats[i].r.z;
    }
    int nx=(r.x-l.x)/dw+1, ny=(r.y-l.y)/dw+1, nz=(r.z-l.z)/dw+1, nxy=nx*ny;
    std::valarray<std::vector<int> > cells(nx*ny*nz);
    std::cerr<<dc<<" size : "<<nx<<","<<ny<<","<<nz<<"\n";
    //places atoms in cells
    pos d;
    int px, py, pz;
    for(unsigned long i=0; i<al.nat; ++i)
    {
        d=al.ats[i].r-l;
        px=d.x/dw; py=d.y/dw; pz=d.z/dw;
        cells[px+nx*py+nxy*pz].push_back(i);
    }
    
    int ix, iy, iz;
    double d2ij;
    for (ix=0; ix< nx; ix++)
        for (iy=0; iy< ny; iy++)
            for (iz=0; iz< nz; iz++)
    {
        int i1=ix+nx*iy+nxy*iz;
        std::vector<int> &cc=cells[i1];
        //on-site
        for(unsigned long i=0; i<cc.size(); ++i)
        {
            for(unsigned long j=0; j<i; ++j)
            {   d=al.ats[cc[i]].r-al.ats[cc[j]].r; d2ij=d.x*d.x+d.y*d.y+d.z*d.z;
                // std::cerr<<cc[i]<<","<<cc[j]<<" : "<<d2ij<<"," <<dc2<<"\n";
                if (d2ij<dc2) 
                {
                    neighs[cc[i]].push_back(cc[j]);
                    neighs[cc[j]].push_back(cc[i]);
                }
            }
        }
        int i2, wx, wy, wz;
        std::vector<int> co(0);
        for (wx=ix-1; wx<=ix+1; wx++)
            for (wy=iy-1; wy<=iy+1; wy++)
                for (wz=iz-1; wz<=iz+1; wz++)
        {
            i2=wx+nx*wy+nxy*wz;
            if (i2>0 && i2<i1) co.push_back(i2);
        }
            
        for (int ic=0; ic<co.size(); ++ic)
        {
            std::vector<int> &c2=cells[co[ic]];
            for(unsigned long i=0; i<cc.size(); ++i)
                for(unsigned long j=0; j<c2.size(); ++j)
                {   
                    d=al.ats[cc[i]].r-al.ats[c2[j]].r; d2ij=d.x*d.x+d.y*d.y+d.z*d.z;
                    if (d2ij<dc2) 
                    {
                        neighs[cc[i]].push_back(c2[j]);
                        neighs[c2[j]].push_back(cc[i]);
                    }
                }
        }
    }
    /*
    pos d;
    double d2ij;
    for(unsigned long i=0; i<al.nat; ++i)
    {
        for(unsigned long j=0; j<i; ++j)
        {
            d=al.ats[i].r-al.ats[j].r;
            
            if (fabs(d.x)<dc && fabs(d.y)<dc && fabs(d.z)<dc)
            {
                d2ij=d.x*d.x+d.y*d.y+d.z*d.z;
    //            std::cerr<<i<<","<<j<<" :: "<<d.x<<".."<<d2ij<<"\n";
                if (d2ij<dc2) 
                {
                    nij.idx=j; nji.idx=i; nij.dst=nji.dst=sqrt(d2ij);
                    neighs[i].push_back(nij);
                    neighs[j].push_back(nji);
                }
            }
        }
    }
    */
    
    hrt.stop();
    std::cerr<<"NEIGHBOUR LIST DONE "<<hrt*1e-6<<"\n";
    //FMatrix<double> slist(al.nat,6);
    
    
    /* reference vectors and inverse matrix*/
    std::valarray<pos> R(3), R0(3), RR(3);
    R0[0].x=1.; R0[0].y=0; R0[0].z=0.;
    R0[1].x=0.5; R0[1].y=sqrt(3.)/2.; R0[1].z=0.;
    R0[2].x=0.5; R0[2].y=sqrt(3.)/6.; R0[2].z=sqrt(2./3.);
    (R0[0])*=a0; R0[1]*=a0; R0[2]*=a0;
    
    FMatrix<double> IR(3,3), R1(3,3), eps(3,3), O(3,3);
    IR(0,0)=1.; IR(0,1)=-1./sqrt(3.); IR(0,2)=-1./sqrt(6.);
    IR(1,0)=0.; IR(1,1)= 2./sqrt(3.); IR(1,2)=-1/sqrt(6.);
    IR(2,0)=0.; IR(2,1)=0.;           IR(2,2)= sqrt(3./2.);
    IR*=1./a0;
    
    rstrain*=0.;
    hrt.start();
    for(unsigned long i=0; i<al.nat; ++i)
    {
        //std::cerr<<i<<"::"<<neighs[i].size()<<"\n";
        if (neighs[i].size()!=4)  {continue;}
        //{  slist(i,0)=slist(i,1)=slist(i,2)=slist(i,3)=slist(i,4)=slist(i,5)=0.;}
        else
        {
            for (int w=0; w<3; ++w)
                R[w]=al.ats[neighs[i][1+w]].r-al.ats[neighs[i][0]].r;
            
            //makes sure that R is right-handed
            if(dot(R[2],cross(R[0],R[1]))<0.)
            { RR[0]=R[0];  R[0]=R[1]; R[1]=RR[0]; }
            
            //rotate so as to align as much as possible with reference triple
            rmsd(R0,R,O);
            //now it would be that OR=
            R1(0,0)=R[0].x; R1(1,0)=R[0].y; R1(2,0)=R[0].z;
            R1(0,1)=R[1].x; R1(1,1)=R[1].y; R1(2,1)=R[1].z; 
            R1(0,2)=R[2].x; R1(1,2)=R[2].y; R1(2,2)=R[2].z;  
            
            /* compute the strain in the frame of reference of the local tetrahedron, non the reference*/
            mult(R1,IR,eps);
            mult(eps,O,R1);
            //mult(O,R1,eps);
            //mult(eps,IR,R1);
            eps=R1;
            
            rstrain(i,0)=eps(0,0)-1.;
            rstrain(i,1)=eps(1,1)-1.;
            rstrain(i,2)=eps(2,2)-1.;
            rstrain(i,3)=(eps(0,1)+eps(1,0))/2;
            rstrain(i,4)=(eps(0,2)+eps(2,0))/2;
            rstrain(i,5)=(eps(1,2)+eps(2,1))/2;
            
            /*
            std::cerr<<" =================  \n";
            for (int i=0; i<3; ++i)
            {
                for (int j=0; j<3; ++j)
                    std::cerr<<eps(i,j)<< " ";
                std::cerr<<"\n";
            }
            std::cerr<<" =================  \n";*/
            
            /*
            std::cerr<<R[0].x<<" "<<RR[0].x<<" "<<R0[0].x<<"\n";
            std::cerr<<R[0].y<<" "<<RR[0].y<<" "<<R0[0].y<<"\n";
            std::cerr<<R[0].z<<" "<<RR[0].z<<" "<<R0[0].z<<"\n\n";
            std::cerr<<R[1].x<<" "<<RR[1].x<<" "<<R0[1].x<<"\n";
            std::cerr<<R[1].y<<" "<<RR[1].y<<" "<<R0[1].y<<"\n";
            std::cerr<<R[1].z<<" "<<RR[1].z<<" "<<R0[1].z<<"\n\n";
            std::cerr<<R[2].x<<" "<<RR[2].x<<" "<<R0[2].x<<"\n";
            std::cerr<<R[2].y<<" "<<RR[2].y<<" "<<R0[2].y<<"\n";
            std::cerr<<R[2].z<<" "<<RR[2].z<<" "<<R0[2].z<<"\n\n";
            */
            
            
        }
    }
    hrt.stop();
    std::cerr<<" STRAIN DONE "<<hrt*1e-6<<std::endl;
}

int main(int argc, char **argv)
{
    CLParser clp(argc, argv);
    frame cf; FMatrix<double> sval;
    while (std::cin.good()) 
    {
        hrt.start();
        std::cin>>cf;
        if (!std::cin.good()) break;
        hrt.stop();
        std::cerr<<cf.nat<<" READ FRAME "<<hrt*1e-6<<"\n";
        strain(cf,sval);
        
        std::cout.precision(5);
        std::cout.setf(std::ios::scientific);
        std::cout.width(13);
        
        std::cout<<cf.nat<<"\n"<<"properties\n";
        for (int i=0; i< cf.nat; ++i)
        {
            std::cout<<std::setw(4)<<cf.ats[i].s<<" "
                    <<std::setw(12)<<cf.ats[i].r.x<<" "
                    <<std::setw(12)<<cf.ats[i].r.y<<" "
                    <<std::setw(12)<<cf.ats[i].r.z<<" "
                    <<std::setw(12)<<sval(i,0)<<" "
                    <<std::setw(12)<<sval(i,1)<<" "
                    <<std::setw(12)<<sval(i,2)<<" "
                    <<std::setw(12)<<sval(i,3)<<" "
                    <<std::setw(12)<<sval(i,4)<<" "
                    <<std::setw(12)<<sval(i,5)<<" "
                    <<std::setw(12)<<(sval(i,0)+sval(i,1)+sval(i,2))/3.<<" "
                    <<"\n";
        }
    }
    return 0;
}
