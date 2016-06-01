#include <algorithm>
#include "rndgen.hpp"
#include "ioatoms.hpp"
#include "linalg.hpp"
#include "clparser.hpp"
#include "tools-histogram.hpp"
#include "tools-autocorr.hpp"
#include "fftw3.h"
#include<fstream>
#include <math.h>
//ooops! beware, this requires extensions to STL. should program with conditional compilation
//#define __USEREGEX 1
#ifdef __USEREGEX
//#include <regex>
#include <boost/regex.hpp>
using namespace boost;
bool chk_sel(const std::string& label, const std::string& regex)
{
    boost::regex pattern(regex);
    return regex_match(label,pattern);
}
#else
bool chk_sel(const std::string& label, const std::string& regex)
{ return (regex=="*")||(label==regex); }
#endif

#include "matrix-full.hpp"
#include "matrix-io.hpp"
#include "ioparser.hpp"
using namespace toolbox;
void inline micdo(const double& a, double &d)
{
    if (a!=0.)
    {
        d/=a;
        d-=round(d);
        d*=a;
    }
}
void micpbc(const double &ax, const double &ay, const double &az,
            double& dx, double& dy, double& dz)
{
    micdo(ax,dx); micdo(ay,dy); micdo(az,dz);
}

void micmat(const FMatrix<double>& m, double &x, double &y, double &z)
{
    double ox=x, oy=y, oz=z;
    x=m(0,0)*ox+m(0,1)*oy+m(0,2)*oz;
    y=m(1,0)*ox+m(1,1)*oy+m(1,2)*oz;
    z=m(2,0)*ox+m(2,1)*oy+m(2,2)*oz;
}

void micrmsd(const std::vector<AtomData>& r, const std::vector<AtomData>& p, FMatrix<double>& RMAT)
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

    D(1,0)=D(0,1);
    D(2,0)=D(0,2);
    D(3,0)=D(0,3);
    D(2,1)=D(1,2);
    D(3,1)=D(1,3);
    D(3,2)=D(2,3);

    FMatrix<double> E; std::valarray<double> v;
    EigenSolverSym(D, E, v);

    for (int i=0; i<4; ++i) q[i]=E(i,0);

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
}

//align frame P minimizing MSD from frame R wrt atoms with name alab
void micalign(const AtomFrame& pr, AtomFrame& p, const std::vector<unsigned long>& asel, bool dopbc)
{
    FMatrix<double> rmat;
    AtomFrame r=pr;

    //Metric matrices
    FMatrix<double> CM(3,3), ICM;
    if (dopbc)
    {
        CM(0,0)=r.nprops["axx"];
        CM(1,0)=r.nprops["axy"];
        CM(2,0)=r.nprops["axz"];
        CM(0,1)=r.nprops["ayx"];
        CM(1,1)=r.nprops["ayy"];
        CM(2,1)=r.nprops["ayz"];
        CM(0,2)=r.nprops["azx"];
        CM(1,2)=r.nprops["azy"];
        CM(2,2)=r.nprops["azz"];
        MatrixInverse(CM,ICM);
    }

    if (r.ats.size()!=p.ats.size()) ERROR("Atom number in frames mismatches");
    //first, align "COM"
    AtomData rcom, pcom; double nlab;
    std::vector<AtomData> ral, pal;

    //first center the reference structure.
    //computes COM
    rcom.x=rcom.y=rcom.z=pcom.x=pcom.y=pcom.z=nlab=0.;
    for (int k=0; k<asel.size(); ++k)
    {
        unsigned long i=asel[k];
        if (r.ats[i].name!=p.ats[i].name) ERROR("Atom label in frames mismatches");
        nlab+=1.;
        rcom.x+=r.ats[i].x; rcom.y+=r.ats[i].y; rcom.z+=r.ats[i].z;
        pcom.x+=p.ats[i].x; pcom.y+=p.ats[i].y; pcom.z+=p.ats[i].z;
    }

    rcom.x/=nlab; rcom.y/=nlab; rcom.z/=nlab;
    pcom.x/=nlab; pcom.y/=nlab; pcom.z/=nlab;
    //center and apply PBC to reference structure
    for (int i=0; i<r.ats.size(); ++i)
    {
        r.ats[i].x-=rcom.x;  r.ats[i].y-=rcom.y;  r.ats[i].z-=rcom.z;
        if (dopbc)
        {
            micmat(ICM,r.ats[i].x,r.ats[i].y,r.ats[i].z);
            micpbc(1.,1.,1.,r.ats[i].x,r.ats[i].y,r.ats[i].z);
            micmat(CM,r.ats[i].x,r.ats[i].y,r.ats[i].z);
        }
    }

    //removes COM, applies PBC and builds alignment data vectors
    for (int i=0; i<r.ats.size(); ++i)
    {
        p.ats[i].x-=pcom.x;  p.ats[i].y-=pcom.y;  p.ats[i].z-=pcom.z;
        if (dopbc)
        {
            micmat(ICM,p.ats[i].x,p.ats[i].y,p.ats[i].z);
            micpbc(1.,1.,1.,p.ats[i].x,p.ats[i].y,p.ats[i].z);
            micmat(CM,p.ats[i].x,p.ats[i].y,p.ats[i].z);
        }
    }

    for (int k=0; k<asel.size(); ++k)
    {
        unsigned long i=asel[k];
        ral.push_back(r.ats[i]);
        pal.push_back(p.ats[i]);
    }

    micrmsd(ral,pal,rmat);
    for (int i=0; i<r.ats.size(); ++i)
    {
        micmat(rmat, p.ats[i].x, p.ats[i].y, p.ats[i].z);
        p.ats[i].x+=rcom.x; p.ats[i].y+=rcom.y; p.ats[i].z+=rcom.z;
    }
}

void af2cm(AtomFrame& af, FMatrix<double>& CM)
{
    if (CM.rows()!=3 || CM.cols()!=3) ERROR("Cell matrix should be already initialized here...");
    CM(0,0)=af.nprops["axx"];
    CM(1,0)=af.nprops["axy"];
    CM(2,0)=af.nprops["axz"];
    CM(0,1)=af.nprops["ayx"];
    CM(1,1)=af.nprops["ayy"];
    CM(2,1)=af.nprops["ayz"];
    CM(0,2)=af.nprops["azx"];
    CM(1,2)=af.nprops["azy"];
    CM(2,2)=af.nprops["azz"];
}

void cub_compute_cr(double x, double r, double dr, double h, double&f)
{
    if (x<=r)
    {
        f=1.;
    }
    else if (x<(r+dr))
    {
        double y=(x-r)/dr;
        double cr, crh, dcr; cr=y-1.; cr*=cr; cr*=(1.+y+y); crh=pow(cr,h);
        f=crh;
    }
    else
    {
        f=0.;
    }
}

void cub_compute_ca(double rx, double ry, double rz, double r, double& f)
{
    double nx=rx/r, ny=ry/r, nz=rz/r;
    double nx2=nx*nx, ny2=ny*ny, nz2=nz*nz;
    double nx4=nx2*nx2, ny4=ny2*ny2, nz4=nz2*nz2;
    double nx2y2=nx2*ny2, ny2z2=ny2*nz2, nz2x2=nz2*nx2;
    double px2y2=nx2+ny2, py2z2=ny2+nz2, pz2x2=nz2+nx2;
    double px2y2_2=px2y2*px2y2, py2z2_2=py2z2*py2z2, pz2x2_2=pz2x2*pz2x2;
    double px2y2_4=px2y2_2*px2y2_2, py2z2_4=py2z2_2*py2z2_2, pz2x2_4=pz2x2_2*pz2x2_2;


    double rpp12=
                nx4*ny4*(nx4+ny4)+ny4*nz4*(ny4+nz4)+nz4*nx4*(nz4+nx4)+2.*(
                         nx4*nx2y2*(py2z2_2-ny2z2)+
                ny4*ny2z2*(pz2x2_2-nz2x2)+
                nz4*nz2x2*(px2y2_2-nx2y2));
    double tx=(py2z2_4-2*ny2z2*(ny4+nz4));
    double ty=(pz2x2_4-2*nz2x2*(nz4+nx4));
    double tz=(px2y2_4-2*nx2y2*(nx4+ny4));

    double xpp12= nx2*(tx-nx4*(ny4+nz4))- ny4*nz4*py2z2*(2+9*nx2);
    double ypp12= ny2*(ty-ny4*(nz4+nx4))- nz4*nx4*pz2x2*(2+9*ny2);
    double zpp12= nz2*(tz-nz4*(nx4+ny4))- nx4*ny4*px2y2*(2+9*nz2);


// weird numbers are there to normalize so that perfect fcc=1 and isotropic liquid=0
    f=2288./79.*(rpp12-4./143.);
}

double get_cv(std::vector<AtomData>& al, unsigned long iat, std::vector<double>& pars, const FMatrix<double>& CM, const FMatrix<double>& ICM, unsigned long cvtype)
{
    int nat=al.size();
    int j;
    double rx, ry, rz, rij, crij, caij, tca, tcr, tcv;
    double r1, dr, h, r02;
    switch (cvtype) {
        case 1: // FCC order parameter, aligned to cartesian axes
        r1 = pars[0]; dr = pars[2]-pars[0]; h = (2.+pow((pars[2]-r1)/(pars[1]-r1),2))/6.; r02=pars[2]*pars[2];
        tca=tcr=0.;
        for(j=0;j<nat;++j) {
            if (j==iat) continue;
            rx=al[j].x-al[iat].x; ry=al[j].y-al[iat].y; rz=al[j].z-al[iat].z;

            micmat(ICM,rx,ry,rz);
            micpbc(1.,1.,1.,rx,ry,rz);
            micmat(CM,rx,ry,rz);
            rij=rx*rx+ry*ry+rz*rz;

            if (rij>r02) continue; //Useless to compute angular part if radial=0
            rij=sqrt(rij);
            cub_compute_cr(rij, r1, dr, h, crij);
            cub_compute_ca(rx, ry, rz, rij, caij);
            tcr+=crij; tca+=caij*crij;
        }

        tcv=tca/tcr;
        break;
        case 2: // same-kind coordination number
            r1 = pars[0]; dr = pars[2]-pars[0]; h = (2.+pow((pars[2]-r1)/(pars[1]-r1),2))/6.; r02=pars[2]*pars[2];
            tcv=0.;
            for(j=0;j<nat;++j) {
                if (j==iat) continue;
                rx=al[j].x-al[iat].x; ry=al[j].y-al[iat].y; rz=al[j].z-al[iat].z;

                micmat(ICM,rx,ry,rz);
                micpbc(1.,1.,1.,rx,ry,rz);
                micmat(CM,rx,ry,rz);
                rij=rx*rx+ry*ry+rz*rz;

                if (rij>r02) continue; //Useless to compute angular part if radial=0
                rij=sqrt(rij);
                cub_compute_cr(rij, r1, dr, h, crij);
                tcv+=crij;
            }
        break;
    }
    return tcv;
} //End of subroutine

#define BOHR2ANG 0.529177

void banner()
{
    std::cerr
            << " USAGE: trajworks  [options] < INPUT [ > OUTPUT ]                               \n"
            << " reads atomic trajectories and computes a number of average structural and      \n"
            << " dynamical properties.                                                          \n"
            << " ## GENERAL OPTIONS:                                                            \n"
            << " -ts       timestep (time span between two frames in input file)                \n"
            << " -o        prefix for output file names (otherwise everything goes to stdout)   \n"
            << " -fstart   frame to start at                                                    \n"
            << " -fstop    frame to stop at                                                     \n"
            << " -fstep    use just a frame every fstep frames in input                         \n"
            << " -ixyz     input is an xyz trajectory. velocities might be read as 4,5,6 field  \n"
            << " -ipdb     input is a PDB format, a portable format for trajectories            \n"
            << " -idlp     input is a DLPOLY HISTORY file                                       \n"
            << " -ref      reads reference for alignment, etc. from the given file              \n"
            << " -box file reads cell parameters from file ( axx axy axz\n ayx ayy ... )        \n"
            << " -qat file reads atomic charges from file ( label q1 \n label q2  ... )         \n"
            << " -unwrap   unwraps PBCs, i.e. makes sure trajectories are continuous....        \n"
            << " -vbox     enables variable box (to be read from DLP or sequentially from box)  \n"
            << " -lab file reads atoms labels from file ( lab1 lab2 lab3....  )                 \n"
            << " -hwin     window for all histo (triangle|box|delta|gauss-[1,2,3,5]) [delta]    \n"
            << " -hwinfac  size of the window, as a function of the bin size [1.0]              \n"
            << " -weights file read statistical log-weights of frames from an external file     \n"
            << " ## g(r) OPTIONS:   activate by -gr                                             \n"
            << " -gr1      label of the first specie  [*]                                       \n"
            << " -gr2      label of the second specie [*]                                       \n"
            << " -grmax    maximum distance to compute g(r) [5]                                 \n"
            << " -grbins   number of bins in the histogram for g(r) [100]                       \n"
            << " ## order parameter OPTIONS: activate by -cv                                    \n"
            << " -cvat     atom types to compute cv for [*]                                     \n"
            << " -cvtype   index of cv type [1]                                                 \n"
            << " -cvpars   cv parameters values [ xxx yyy zzz ]                                 \n"
            << " ## 3d density OPTIONS:  activate by -dens                                      \n"
            << " -dat      label of the monitored specie [*]                                    \n"
            << " -dbins    number of bins (either -dbins n or -dbins nx,ny,nz)                  \n"
            << " -dbinw    use a box function to smear the density over a DISTANCE w            \n"
            << " -dfold    fold the density in a smaller cell (-dfold fx,fy,fz)                 \n"
            << " -drange   only accumulate density in the specified FOLDED range                \n"
            << "           (-drange ax,bx,ay,by,az,bz. defaults to -0.5,0.5,-0.5,0.5,-0.5,0.5)  \n"
            << " -dtraj    also prints out a trajectory file of the (aligned) geometries        \n"
            << " -dpov     prints also a df3 density file, to be used in POV-Ray                \n"
            << " -dproj    also computes projection of density on xy,yz,xz planes and x,y,z axis\n"
            << " -dalign   align the box to the atoms with the given label, then accumulates    \n"
            << " ## principal component analysis options: activate by -pca                      \n"
            << " -pcat     label of the monitored specied [*]                                   \n"
            << " -pcaxyz   output a .xyz of the average structure, together with Gaussian file  \n"
            << " -pcalign  align the frames wrt this atoms on a RMSD basis                      \n"
            << " -pcanocov only computes the mean position of atoms, not covariance of PCA      \n"
            << " ## mean square displacement calculation (-msd)                                 \n"
            << " -msdat    label of the monitored specie [*]                                    \n"
            << " -msdlag   maximum time-lag to compute msd for [1000]                           \n"
            << " ## velocity-velocity correlation options (-vvac)                               \n"
            << " -vvat     label of the monitored specie [*]                                    \n"
            << " -vvlag    maximum time-lag to compute vvac for [1000]                          \n"
            << " -vvftpad  pad the ACF with zeroes before transforming [0]                      \n"
            << " -vvftbox  smear the ACF before transforming using a triangle function [false]  \n"
            << " -vvindex  select the indices of the atoms rather than their label []           \n"
            << " ## atom momentum distribution options (-pd )                                   \n"
            << " -pdat     label of the monitored specie [*]                                    \n"
            << " -pdvec    compute distribution projected on a direction [ distr. |p| if abs.]  \n"
            << " -pdbins   number of bins [100]                                                 \n"
            << " -pdmax    maximum momentum to bin [1000]                                       \n"
            << " ## 3D atom momentum distribution options (-p3d )                               \n"
            << " -p3dat     label of the monitored specie [*]                                   \n"
            << " -p3dinvert enforces inversion simmetry [false]                                 \n"
            << " -p3dbins   number of bins [100]                                                \n"
            << " -p3dmax    maximum momentum to bin [1000]                                      \n"
            << " ## projected momentum ( -pproj )                                               \n"
            << " -ppat1     label of the monitored specie [*]                                   \n"
            << " -ppat2     label of specie to look for nearest bond [*]                        \n"
            << " ## thermal ellipsoid tensor ( -thermal )                                       \n"            //may want to include options for center of mass removal and Berry phase
            ;
}

int main(int argc, char **argv)
{
    CLParser clp(argc, argv);

    bool fgdr=false, fvvac=false, fpdb=false, fxyz=false, fdlp=false, fmsd=false, fdipole=false, fdens=false, fdtraj=false, fdproj=false, fdpov=false, fpca=false, fpcaxyz=false, fpcanocov=false, fvbox=false, fcv=false, fpd=false, fvvacbox=false, fp3d=false, fp3dinv=false, funwrap=false, fpproj=false, fthermal=false, fcharge=false, fhelp;
    std::string lgdr1, lgdr2, dummy, lvvac, lmsd, ldens, ldalign, lpcalign, lpcat, prefix, fbox, fqat, sdbins, sdfold, sdrange, fref, lpdat,shwin, lp3dat, flab, lppat1, lppat2, lcv, fweights;
    double cogdr, dt, densw, pdmax, p3dmax, hwinfac; unsigned long fstart,fstop,fstep,gdrbins, vvlag, msdlag, ftpad, dbinsx, dbinsy, dbinsz, dfoldx, dfoldy, dfoldz, cvtype, pdbins, p3dbins;
    double drangeax, drangebx,  drangeay, drangeby,  drangeaz, drangebz;
    std::vector<double> cvpars, pdvec; std::vector<unsigned long> vvindex;
    bool fok=
            //general options
            clp.getoption(prefix,"o",std::string("")) &&
            clp.getoption(dt,"ts",1.) &&
            clp.getoption(fstart,"fstart",(unsigned long) 0) &&
            clp.getoption(fstop,"fstop",(unsigned long) 0) &&
            clp.getoption(funwrap,"unwrap",false) &&
            clp.getoption(fstep,"fstep",(unsigned long) 0) &&
            clp.getoption(fxyz,"ixyz",true) &&
            clp.getoption(fpdb,"ipdb",false) &&
            clp.getoption(fdlp,"idlp",false) &&
            clp.getoption(fref,"ref",std::string("")) &&
            clp.getoption(fbox,"box",std::string("")) &&
            clp.getoption(fqat,"qat",std::string("")) &&
            clp.getoption(flab,"lab",std::string("")) &&
            clp.getoption(fvbox,"vbox",false) &&
            clp.getoption(shwin,"hwin",std::string("delta")) &&
            clp.getoption(hwinfac,"hwinfac",1.) &&
            clp.getoption(fweights,"weights",std::string("")) &&
            //g(r) options
            clp.getoption(fgdr,"gr",false) &&            
            clp.getoption(lgdr1,"gr1",std::string("*")) &&
            clp.getoption(lgdr2,"gr2",std::string("*")) &&
            clp.getoption(cogdr,"grmax",5.) &&
            clp.getoption(gdrbins,"grbins",(unsigned long)100) &&
            //vvacf options
            clp.getoption(fvvac,"vvac",false) &&
            clp.getoption(fvvacbox,"vvftbox",false) &&
            clp.getoption(vvindex,"vvindex",std::vector<unsigned long>(0)) &&
            clp.getoption(vvlag,"vvlag",(unsigned long) 100) &&
            clp.getoption(ftpad,"vvftpad",(unsigned long) 0) &&
            clp.getoption(lvvac,"vvat",std::string("*")) &&
            //msd  options
            clp.getoption(fmsd,"msd",false) &&
            clp.getoption(msdlag,"msdlag",(unsigned long) 1000) &&
            clp.getoption(lmsd,"msdat",std::string("*")) &&
            //density options
            clp.getoption(fcv,"cv",false) &&
            clp.getoption(lcv,"cvat",std::string("*")) &&
            clp.getoption(cvtype,"cvtype",(unsigned long) 1) &&
            clp.getoption(cvpars,"cvpars",std::vector<double>(0)) &&
            //density options
            clp.getoption(fdens,"dens",false) &&
            clp.getoption(sdbins,"dbins",std::string("50")) &&
            clp.getoption(sdfold,"dfold",std::string("")) &&
            clp.getoption(sdrange,"drange",std::string("")) &&
            clp.getoption(densw,"dbinw",0.) &&
            clp.getoption(ldalign,"dalign",std::string("")) &&
            clp.getoption(fdtraj,"dtraj",false) &&
            clp.getoption(fdpov,"dpov",false) &&
            clp.getoption(fdproj,"dproj",false) &&
            clp.getoption(ldens,"dat",std::string("*")) &&
            //PCA options
            clp.getoption(fpca,"pca",false) &&
            clp.getoption(fpcanocov,"pcanocov",false) &&
            clp.getoption(fpcaxyz,"pcaxyz",false) &&
            clp.getoption(lpcalign,"pcalign",std::string("")) &&
            clp.getoption(lpcat,"pcat",std::string("")) &&
            //momentum distrib. options
            clp.getoption(fpd,"pd",false) &&
            clp.getoption(pdvec,"pdvec",std::vector<double>(0)) &&
            clp.getoption(pdbins,"pdbins",(unsigned long) 100) &&
            clp.getoption(pdmax,"pdmax",(double) 100.) &&
            clp.getoption(lpdat,"pdat",std::string("*")) &&
            //3D momentum distrib. options
            clp.getoption(fp3d,"p3d",false) &&
            clp.getoption(fp3dinv,"p3dinvert",false) &&
            clp.getoption(p3dbins,"p3dbins",(unsigned long) 100) &&
            clp.getoption(p3dmax,"p3dmax",(double) 100.) &&
            clp.getoption(lp3dat,"p3dat",std::string("*")) &&
            //projected momentum options
            clp.getoption(fpproj,"pproj",false) &&
            clp.getoption(lppat1,"ppat1",std::string("*")) &&
            clp.getoption(lppat2,"ppat2",std::string("*")) &&
            //mollified charge options
            clp.getoption(fcharge,"charge",false) &&
           // clp.getoption(qsmear,"sigma",(double) 1.0) &&
            //dipole options
            clp.getoption(fdipole,"dpl",false) &&
            clp.getoption(fhelp,"h",false)&&
            //thermal ellipsoids options
            clp.getoption(fthermal,"thermal",false)
            ;


    if (fhelp || ! fok) { banner(); exit(-1); }
    //general-purpose stuff

    HGWindowMode hwin; HGOptions<Histogram<double> > hgwo;
    if (shwin=="triangle") hwin=HGWTriangle;
    else if (shwin=="box") hwin=HGWBox;
    else if (shwin=="delta") hwin=HGWDelta;
    else if (shwin=="gauss-1") hwin=HGWGauss1;
    else if (shwin=="gauss-2") hwin=HGWGauss2;
    else if (shwin=="gauss-3") hwin=HGWGauss3;
    else if (shwin=="gauss-5") hwin=HGWGauss5;
    else ERROR("Unsupported histogram windowing mode");


    AtomFrame af; std::vector<AtomData> al1, al2, ldip;
    AtomData dip_tx, dip_cur; double tcharge=0.;
    Histogram<double> hgdr(0.,cogdr, gdrbins);
    hgdr.get_options(hgwo); hgwo.window=hwin; hgwo.walls=HGBHard;
    hgwo.window_width=(hgwo.boundaries[1]-hgwo.boundaries[0])*hwinfac;   hgdr.set_options(hgwo);

    double d12, dx, dy, dz, cog2=cogdr*cogdr, gdrw, gdrwtot;
    int nfr=0, npfr=0;
    
    //velocity-velocity correlation stuff
    double vvnat=0; std::valarray<AutoCorrelation<double> > vvacf; std::valarray<bool> fvvac_inc;

    //MSD stuff
    unsigned long imsd=0 , msdnat=0; std::valarray<unsigned long> nmsd(msdlag); 
    std::valarray<AtomFrame> msdbuff(msdlag);
    std::valarray<double> dmsd(msdlag); nmsd=0; dmsd=0.; FMatrix<bool>  fmsd_inc(msdlag,msdlag);

    //density histograms
    std::valarray<HGOptions<Histogram<double> > > hgo(3);
    NDHistogram<double> ndh(hgo);
    AtomFrame cubefr;
    if (sdbins.find_first_of(',',0)==std::string::npos)
    {
        dbinsz=dbinsy=dbinsx=str2int(sdbins);
    }
    else
    {
        int npos=0;
        npos=sdbins.find_first_of(',');
        while (npos!=std::string::npos)
        { sdbins[npos]=' ';  npos=sdbins.find_first_of(',',npos+1);  }
        std::istringstream iss(sdbins);
        iss>>dbinsx>>dbinsy>>dbinsz;
        if (iss.fail() || iss.bad()) ERROR("Invalid format for dbins (use either -dbins n or -dbins nx,ny,nz).");
    }

    dfoldx=dfoldy=dfoldz=1;
    if (sdfold!="")
    {
        int npos=0;
        npos=sdfold.find_first_of(',');
        while (npos!=std::string::npos)
        { sdfold[npos]=' ';  npos=sdfold.find_first_of(',',npos+1);  }
        std::istringstream iss(sdfold);
        iss>>dfoldx>>dfoldy>>dfoldz;
        if (iss.fail() || iss.bad()) ERROR("Invalid format for dfold (use -dfold fx,fy,fz).");
    }

    drangeax=drangeay=drangeaz=-0.5;
    drangebx=drangeby=drangebz= 0.5;
    if (sdrange!="")
    {
        int npos=0;
        npos=sdrange.find_first_of(',');
        while (npos!=std::string::npos)
        { sdrange[npos]=' ';  npos=sdrange.find_first_of(',',npos+1);  }
        std::istringstream iss(sdrange);
        iss>>drangeax>>drangebx>>drangeay>>drangeby>>drangeaz>>drangebz;
        if (iss.fail() || iss.bad()) ERROR("Invalid format for drange (use -drange ax,bx,ay,by,az,bz).");
    }
    //principal component analysis stuff
    FMatrix<double> pcaxx; std::valarray<double> pcax; AtomFrame pcafr;
    unsigned long pcasz;

    //momentum distribution stuff
    double npdat=0; std::valarray<bool> fpd_inc;
    Histogram<double> hpd(0.,pdmax,pdbins);
    FMatrix<double> pdtens2(3,3,0.); FMatrix<FMatrix<double> > pdtens4(3,3,FMatrix<double>(3,3,0.));

    if (pdvec.size()!=0)
    {
        if (pdvec.size()!=3) ERROR("Invalid vector given for -pdvec");
        double pdvnorm=sqrt(pdvec[0]*pdvec[0]+pdvec[1]*pdvec[1]+pdvec[2]*pdvec[2]);
        for (int i=0; i<3; ++i) pdvec[i]/=pdvnorm;
    }
    hpd.get_options(hgwo); hgwo.window=hwin;
    hgwo.window_width=(hgwo.boundaries[1]-hgwo.boundaries[0])*hwinfac;   hpd.set_options(hgwo);

    //3d momentum distribution stuff
    double np3dat=0; std::valarray<bool> fp3d_inc;
    std::valarray<HGOptions<Histogram<double> > > p3dhgo(3);
    for (int i=0; i<3; ++i) {
      p3dhgo[i].boundaries.resize(p3dbins+1);
      p3dhgo[i].window_width=2.*p3dmax/p3dbins*hwinfac;
      p3dhgo[i].window=hwin;
      for (int k=0; k<p3dbins+1;k++) p3dhgo[i].boundaries[k]=-p3dmax+(2.*k*p3dmax/p3dbins);
    }
    NDHistogram<double> p3dh(p3dhgo);

    //thermal ellipsoids stuff
    FMatrix<double> te_u0, te_uiuj;

    //initializes output streams
    std::ostream *ogdr, *ocv, *ovvac, *omsd, *odipole, *odens, *odtraj, *odpov, *opca, *opcaxyz, *opd, *op3d, *opproj, *otherm, *ocharge;
    std::cout.precision(8); std::cout.width(15); std::cout.setf(std::ios::scientific);

    if (prefix=="")
    {
        //everything goes to stdout
        ocv=ogdr=ovvac=omsd=odipole=odens=odtraj=opca=opcaxyz=opd=op3d=opproj=otherm=ocharge=&std::cout;
    }
    else
    {
        if (fgdr)
        {
            ogdr=new std::ofstream((prefix+std::string(".gdr")).c_str());
            (*ogdr).precision(8); (*ogdr).width(15); (*ogdr).setf(std::ios::scientific);
        }
        if (fcv)
        {
            ocv=new std::ofstream((prefix+std::string(".cv")).c_str());
            (*ocv).precision(8); (*ocv).width(15); (*ocv).setf(std::ios::scientific);
        }
        if (fpd)
        {
            opd=new std::ofstream((prefix+std::string(".pd")).c_str());
            (*opd).precision(8); (*opd).width(15); (*opd).setf(std::ios::scientific);
        }
        if (fp3d)
        {
            op3d=new std::ofstream((prefix+std::string(".p3d")).c_str());
            (*op3d).precision(8); (*op3d).width(10); (*op3d).setf(std::ios::scientific);
        }
        if (fpproj)
        {
            opproj=new std::ofstream((prefix+std::string(".pproj")).c_str());
            (*opproj).precision(8); (*opproj).width(10); (*opproj).setf(std::ios::scientific);
        }
        if (fvvac)
        {
            ovvac=new std::ofstream((prefix+std::string(".vvac")).c_str());
            (*ovvac).precision(8); (*ovvac).width(15); (*ovvac).setf(std::ios::scientific);
        }
        if (fmsd)
        {
            omsd=new std::ofstream((prefix+std::string(".msd")).c_str());
            (*omsd).precision(8); (*omsd).width(15); (*omsd).setf(std::ios::scientific);
        }
        if (fdipole)
        {
            odipole=new std::ofstream((prefix+std::string(".dpl")).c_str());
            (*odipole).precision(8); (*odipole).width(15); (*odipole).setf(std::ios::scientific);
        }
        if (fdens)
        {
            odens=new std::ofstream((prefix+std::string(".dens")).c_str());
            (*odens).precision(8); (*odens).width(15); (*odens).setf(std::ios::scientific);
            if (fdtraj)
            {
                odtraj=new std::ofstream((prefix+std::string(".dtraj.xyz")).c_str());
                (*odtraj).precision(8); (*odtraj).width(15); (*odtraj).setf(std::ios::scientific);
            }
            if (fdpov)
                odpov=new std::ofstream((prefix+std::string(".df3")).c_str(), std::ios_base::binary);
        }
        if (fcharge)
        {
            ocharge=new std::ofstream((prefix+std::string(".charge")).c_str());
        }
        if (fpca)
        {
            opca=new std::ofstream((prefix+std::string(".log")).c_str());
            (*opca).precision(8); (*opca).width(15); (*opca).setf(std::ios::scientific);
            if (fpcaxyz)
            {
                opcaxyz=new std::ofstream((prefix+std::string(".pca.xyz")).c_str());
                (*opcaxyz).precision(8); (*opcaxyz).width(15); (*opcaxyz).setf(std::ios::scientific);
            }
        }
        if (fthermal)
        {
            otherm=new std::ofstream((prefix+std::string(".therm")).c_str());
            (*otherm).precision(8); (*otherm).width(15); (*otherm).setf(std::ios::scientific);
        }
    }

    //someone might need a gaussian prn
    toolbox::RndGaussian<double,toolbox::MTRndUniform> rng(RGPars<double>(0,1));
    rng.RUGenerator().seed(100);

    //cell matrix reading
    std::ifstream ifbox; 
    FMatrix<double> CBOX(3,3), CM(3,3), ICM(3,3); double cvolume=0.; bool fhavecell=false;
   
    //charge dictionary
    std::ifstream ifqat;
    std::map<std::string, double> qmap;
    if (fqat!="") 
    {
       std::string qlabel; double qq;
       ifqat.open(fqat.c_str());
       while (ifqat.good())
       {
          ifqat>>qlabel>>qq; 
          qmap[qlabel]=qq;
       }
    }
    //weights reading
    std::ifstream ifweights;

    //named selection
    std::vector<unsigned long> denssel, al_denssel, pcasel, al_pcasel, gdrsel,vvacsel;
    AtomFrame reffr;

    std::ifstream iflab; std::vector<std::string > vlab(0);
    if (flab!="")
    {
        std::string lslab;
        iflab.open(flab.c_str());

        while((iflab>>lslab).good()) vlab.push_back(lslab);
    }

    if(fdlp) { std::getline(std::cin,dummy); std::getline(std::cin,dummy); } //jumps over header
    bool ffirstcell;
    AtomFrame uwframe;
    double statweight=1.0, stattot=0.0;
    if (fbox!="") ifbox.open(fbox.c_str());
    if (fweights!="") ifweights.open(fweights.c_str());
    ffirstcell=true;    
    while ((fdlp && ReadDLPFrame(std::cin,af))||(fpdb && ReadPDBFrame(std::cin,af))||(fxyz && ReadXYZFrame(std::cin,af)))
    {
        ++nfr;
        if (fstop!=0 && nfr>fstop)  break;
        // reads anyway box and weights, as they are meant to span the whole trajectory
        if (fbox!="" and (fvbox || ffirstcell) )
        { ffirstcell=false;  for (int i=0;i<3; ++i) for (int j=0;j<3;++j) ifbox>>CBOX(j,i); if (!ifbox.good()) ERROR("Format error in box file."); }
        if (fweights!="")
        {
            ifweights >> statweight; statweight = exp(statweight); // weights are stored as logs 
            if (!ifweights.good()) ERROR("Format error in weights file.");
        }
        if (fqat!="")
        {   // infers charges from atomic labels -- note that undefined labels will be set to have charge zero!
           for (int i=0; i<af.ats.size(); ++i)
           {
              af.ats[i].nprops["charge"]=qmap[af.ats[i].name];
           }
        }        
        if ((fstart!=0 && nfr<fstart) || (fstep>0 && nfr%fstep!=0))
        { std::cerr<<"Skipping frame   "<<std::setw(10)<<std::setiosflags(std::ios::right)<<nfr<<"\r"; continue; }

        std::cerr<<  "Processing frame "<<std::setw(10)<<std::setiosflags(std::ios::right)<<nfr<<"\r";
        if (vlab.size()!=0)
        {
            if (af.ats.size()!=vlab.size()) ERROR("Atom number mismatch at frame"<<nfr);
            for( int i=0 ; i<af.ats.size(); ++i ) af.ats[i].name=vlab[i];
        }
        npfr++; 
        if (!fhavecell) //sets up unit cell for PBC  
        {            
            fhavecell=true;
            if (fbox!="") {  CM = CBOX;  }
            else if (fdlp) { af2cm(af,CM); }
            else if (fpdb) { af2cm(af,CM); }
            else fhavecell=false;
            if (fhavecell)
            {
                MatrixInverse(CM,ICM);

                cvolume=fabs(CM(0,0)*CM(1,1)*CM(2,2)+CM(0,1)*CM(1,2)*CM(2,0)+CM(1,0)*CM(2,1)*CM(0,2)-
                        CM(0,2)*CM(1,1)*CM(2,0)+CM(0,1)*CM(1,0)*CM(2,2)+CM(2,1)*CM(1,2)*CM(0,0));
            }

            if (fref!="")
            {
                //reads reference frame from external file
                std::ifstream ifref(fref.c_str()); bool succ;
                if (fdlp) succ=ReadDLPConf(ifref,reffr);
                else if (fxyz) succ=ReadXYZFrame(ifref,reffr);
                else if (fpdb) succ=ReadPDBFrame(ifref,reffr);
                else ERROR("I don't know how to read this...");

                if (!succ) ERROR("Format error in ref file.");
                if (reffr.ats.size()!=af.ats.size()) ERROR("Atom number mismatch between reference and trajectory frames!");
            }
        }
        if (fvbox)
        {
            if (fbox!="") { CM = CBOX; }
            else if (fdlp) { af2cm(af,CM); }
            else if (fpdb) { af2cm(af,CM); }
            MatrixInverse(CM,ICM);
            //updates the average volume
            cvolume=cvolume+(fabs(CM(0,0)*CM(1,1)*CM(2,2)+CM(0,1)*CM(1,2)*CM(2,0)+CM(1,0)*CM(2,1)*CM(0,2)-
                    CM(0,2)*CM(1,1)*CM(2,0)+CM(0,1)*CM(1,0)*CM(2,2)+CM(2,1)*CM(1,2)*CM(0,0)) - cvolume)/npfr ;
        }
        if (funwrap)
        {
            if (uwframe.ats.size()==0) uwframe=af;
            else
            {
                //unwraps PBCs, i.e. makes sure trajectories are continuous by removing jumps in atomic coordinates.
                double uwdx, uwdy, uwdz;
                for (int i=0; i<af.ats.size(); ++i)
                {
                    uwdx=af.ats[i].x-uwframe.ats[i].x;
                    uwdy=af.ats[i].y-uwframe.ats[i].y;
                    uwdz=af.ats[i].z-uwframe.ats[i].z;
                    micmat(ICM,uwdx,uwdy,uwdz);
                    micpbc(1.,1.,1.,uwdx,uwdy,uwdz);
                    micmat(CM,uwdx,uwdy,uwdz);
                    af.ats[i].x=uwframe.ats[i].x+uwdx;
                    af.ats[i].y=uwframe.ats[i].y+uwdy;
                    af.ats[i].z=uwframe.ats[i].z+uwdz;
                }
                uwframe=af;
            }
        }
        if (fthermal)
        {
            if (te_u0.rows()==0)
            {
               //initializes vectors for average positions and variance
               te_u0.resize(af.ats.size(),3);
               te_uiuj.resize(af.ats.size(),6);
               te_u0*=0.0; te_uiuj*=0.0;
            }
            for (int i=0; i<af.ats.size(); ++i)
            {
               te_u0(i,0)+=af.ats[i].x;
               te_u0(i,1)+=af.ats[i].y;
               te_u0(i,2)+=af.ats[i].z;
               te_uiuj(i,0)+=af.ats[i].x*af.ats[i].x;
               te_uiuj(i,1)+=af.ats[i].y*af.ats[i].y;
               te_uiuj(i,2)+=af.ats[i].z*af.ats[i].z;
               te_uiuj(i,3)+=af.ats[i].x*af.ats[i].y;
               te_uiuj(i,4)+=af.ats[i].x*af.ats[i].z;
               te_uiuj(i,5)+=af.ats[i].y*af.ats[i].z;
            }
        }
        if (fdens)
        {
            if (fvbox) ERROR("OOPS... Variable box is not implemented for this calculation...."); //PLACEHOLDER
            if (hgo[0].boundaries.size()==0)
            {
                if (!fhavecell) ERROR("You must provide cell size in order to compute density hystogram!");
                //forces the data read above
                af.nprops["axx"]=CM(0,0);
                af.nprops["axy"]=CM(1,0);
                af.nprops["axz"]=CM(2,0);
                af.nprops["ayx"]=CM(0,1);
                af.nprops["ayy"]=CM(1,1);
                af.nprops["ayz"]=CM(2,1);
                af.nprops["azx"]=CM(0,2);
                af.nprops["azy"]=CM(1,2);
                af.nprops["azz"]=CM(2,2);

                hgo[0].boundaries.resize(dbinsx+1);
                hgo[1].boundaries.resize(dbinsy+1);
                hgo[2].boundaries.resize(dbinsz+1);
                for (int i=0; i<3; ++i)
                    hgo[i].window_width=densw/sqrt(CM(0,i)*CM(0,i)+CM(1,i)*CM(1,i)+CM(2,i)*CM(2,i));
                double va,vb;
                for (int i=0; i<3; ++i)
                {
                    hgo[i].walls=HGBPeriodic;
                    if (densw==0.) hgo[i].window=HGWDelta; else hgo[i].window=(hwin==HGWDelta?HGWTriangle:hwin);
                    switch(i)
                    {
                        case 0: va=drangeax; vb=drangebx; break;
                        case 1: va=drangeay; vb=drangeby; break;
                        case 2: va=drangeaz; vb=drangebz; break;
                    }

                    for (int k=0; k<hgo[i].boundaries.size();k++)
                        hgo[i].boundaries[k]=va+k*(vb-va)/(hgo[i].boundaries.size()-1);
                }
                //extra PBC folding (if the cell represents several unit cells and want density accumulated in a single one
                hgo[0].boundaries*=(1./dfoldx);
                hgo[1].boundaries*=(1./dfoldy);
                hgo[2].boundaries*=(1./dfoldz);
                ndh.set_options(hgo);

                if (fref!="") cubefr=reffr; else cubefr=af;    //take reference frame as provided


                if (ldalign!="")
                {
                    //we center the "COM" of labelled atoms  for alignment
                    AtomData ldcom; double ldnlab;
                    ldcom.x=ldcom.y=ldcom.z=ldnlab=0.;

                    //picks out the atoms for alignment
                    for (int i=0; i<cubefr.ats.size();++i)  if(chk_sel(cubefr.ats[i].name,ldalign))
                    {
                        ++ldnlab; ldcom.x+=cubefr.ats[i].x; ldcom.y+=cubefr.ats[i].y; ldcom.z+=cubefr.ats[i].z;
                        al_denssel.push_back(i);
                    }
                    ldcom.x/=ldnlab; ldcom.y/=ldnlab; ldcom.z/=ldnlab;

                    for (int i=0; i<cubefr.ats.size();++i)
                    { cubefr.ats[i].x-=ldcom.x;  cubefr.ats[i].y-=ldcom.y;  cubefr.ats[i].z-=ldcom.z; }


                }
                for (int i=0; i<cubefr.ats.size();++i)
                    if(chk_sel(cubefr.ats[i].name,ldens)) denssel.push_back(i);

            }

            AtomFrame laf=af;
            if (ldalign!="")  micalign(cubefr, laf, al_denssel,true);
            if (fdtraj)
            {
                (*odtraj) << laf.ats.size()<<"\nOutput from trajworks density run\n";
                for (unsigned long i=0; i<laf.ats.size(); ++i)
                {
                    (*odtraj)<<laf.ats[i].name<<" "<<laf.ats[i].x<<"  "<<laf.ats[i].y<<"  "<<laf.ats[i].z<<std::endl;
                }
            }
            al1.clear();  for (unsigned long i=0; i<denssel.size(); ++i)   al1.push_back(laf.ats[denssel[i]]);
            std::valarray<double> didata(3);
            for (unsigned long i=0; i<al1.size(); ++i)
            {
                //go to scaled coordinates
                didata[0]=al1[i].x;  didata[1]=al1[i].y; didata[2]=al1[i].z;
                micmat(ICM,didata[0],didata[1],didata[2]);
                if (ldalign=="") micpbc(1./dfoldx,1./dfoldy,1./dfoldz,didata[0],didata[1],didata[2]);  //if we have aligned, we ALREADY HAVE APPLIED PBC AT THE GOOD MOMENT!
                ndh<<didata;
            }
        }
        if (fcharge)
        {
            if (fvbox) ERROR("OOPS... Variable box is not implemented for this calculation...."); //PLACEHOLDER
            if (hgo[0].boundaries.size()==0)
            {
                if (!fhavecell) ERROR("You must provide cell size in order to compute density hystogram!");
                //forces the data read above
                af.nprops["axx"]=CM(0,0);
                af.nprops["axy"]=CM(1,0);
                af.nprops["axz"]=CM(2,0);
                af.nprops["ayx"]=CM(0,1);
                af.nprops["ayy"]=CM(1,1);
                af.nprops["ayz"]=CM(2,1);
                af.nprops["azx"]=CM(0,2);
                af.nprops["azy"]=CM(1,2);
                af.nprops["azz"]=CM(2,2);

                hgo[0].boundaries.resize(dbinsx+1);
                hgo[1].boundaries.resize(dbinsy+1);
                hgo[2].boundaries.resize(dbinsz+1);
                for (int i=0; i<3; ++i)
                    hgo[i].window_width=densw/sqrt(CM(0,i)*CM(0,i)+CM(1,i)*CM(1,i)+CM(2,i)*CM(2,i));
                double va,vb;
                for (int i=0; i<3; ++i)
                {
                    if (densw==0.) hgo[i].window=HGWDelta; else hgo[i].window=(hwin==HGWDelta?HGWTriangle:hwin);
                    hgo[i].walls=HGBPeriodic;
                    switch(i)
                    {
                        case 0: va=drangeax; vb=drangebx; break;
                        case 1: va=drangeay; vb=drangeby; break;
                        case 2: va=drangeaz; vb=drangebz; break;
                    }

                    for (int k=0; k<hgo[i].boundaries.size();k++)
                        hgo[i].boundaries[k]=va+k*(vb-va)/(hgo[i].boundaries.size()-1);
                }
                //extra PBC folding (if the cell represents several unit cells and want density accumulated in a single one
                hgo[0].boundaries*=(1./dfoldx);
                hgo[1].boundaries*=(1./dfoldy);
                hgo[2].boundaries*=(1./dfoldz);
                ndh.set_options(hgo);

                if (fref!="") cubefr=reffr; else cubefr=af;    //take reference frame as provided


                if (ldalign!="")
                {
                    //we center the "COM" of labelled atoms  for alignment
                    AtomData ldcom; double ldnlab;
                    ldcom.x=ldcom.y=ldcom.z=ldnlab=0.;

                    //picks out the atoms for alignment
                    for (int i=0; i<cubefr.ats.size();++i)  if(chk_sel(cubefr.ats[i].name,ldalign))
                    {
                        ++ldnlab; ldcom.x+=cubefr.ats[i].x; ldcom.y+=cubefr.ats[i].y; ldcom.z+=cubefr.ats[i].z;
                        al_denssel.push_back(i);
                    }
                    ldcom.x/=ldnlab; ldcom.y/=ldnlab; ldcom.z/=ldnlab;

                    for (int i=0; i<cubefr.ats.size();++i)
                    { cubefr.ats[i].x-=ldcom.x;  cubefr.ats[i].y-=ldcom.y;  cubefr.ats[i].z-=ldcom.z; }


                }
                for (int i=0; i<cubefr.ats.size();++i)
                    if(chk_sel(cubefr.ats[i].name,ldens)) denssel.push_back(i);

            }

            AtomFrame laf=af;
            if (ldalign!="")  micalign(cubefr, laf, al_denssel,true);
            al1.clear();  for (unsigned long i=0; i<denssel.size(); ++i)   al1.push_back(laf.ats[denssel[i]]);
            std::valarray<double> didata(3);
            for (unsigned long i=0; i<al1.size(); ++i)
            {
                //go to scaled coordinates
                didata[0]=al1[i].x;  didata[1]=al1[i].y; didata[2]=al1[i].z;
                micmat(ICM,didata[0],didata[1],didata[2]);
                if (ldalign=="") micpbc(1./dfoldx,1./dfoldy,1./dfoldz,didata[0],didata[1],didata[2]);  //if we have aligned, we ALREADY HAVE APPLIED PBC AT THE GOOD MOMENT!
               if (al1[i].props.size()<1) ERROR("charge data not available for atom "<<i<<".");
               //std::cerr<<densw<<al1[i].name<<(al1[i].name==std::string("O")?-2.0:1.0)<<std::endl;
                ndh.add(didata, al1[i].props[0]);
            }
        }

        if (fpca)
        {
            if (fvbox) ERROR("OOPS... Variable box is not implemented for this calculation...."); //PLACEHOLDER
            if (pcax.size()==0)
            {
                if (fref!="") pcafr=reffr; else pcafr=af;    //take reference frame as provided
                if (lpcalign!="")
                {
                    //we center the "COM" of labelled atoms  for alignment
                    AtomData lpccom; double lpcnlab;
                    lpccom.x=lpccom.y=lpccom.z=lpcnlab=0.;
                    for (int i=0; i<pcafr.ats.size();++i) if (chk_sel(pcafr.ats[i].name,lpcalign))
                    {
                        lpcnlab+=1.; lpccom.x+=pcafr.ats[i].x; lpccom.y+=pcafr.ats[i].y; lpccom.z+=pcafr.ats[i].z;
                        al_pcasel.push_back(i);
                    }
                    lpccom.x/=lpcnlab; lpccom.y/=lpcnlab; lpccom.z/=lpcnlab;

                    for (int i=0; i<pcafr.ats.size();++i)
                    { pcafr.ats[i].x-=lpccom.x;  pcafr.ats[i].y-=lpccom.y;  pcafr.ats[i].z-=lpccom.z; }

                }
                for (int i=0; i<pcafr.ats.size(); ++i) if (chk_sel(pcafr.ats[i].name, lpcat)) pcasel.push_back(i);
                pcasz=pcasel.size();
                pcax.resize(pcasz*3); pcax=0.;
                if(!fpcanocov) pcaxx.resize(pcasz*3,pcasz*3); pcaxx*=0.;
            }

            AtomFrame laf=af;
            if (lpcalign!="") micalign(pcafr, laf, al_pcasel,false);
            al1.clear();
            for (unsigned long i=0; i<pcasel.size(); ++i) al1.push_back(laf.ats[pcasel[i]]);
            if (al1.size()!=pcasz) ERROR("Number of atoms tagged for PCA changed along the trajectory");

            for (unsigned long i=0; i<al1.size(); ++i)
            {
                pcax[3*i]+=al1[i].x;
                pcax[3*i+1]+=al1[i].y;
                pcax[3*i+2]+=al1[i].z;
                if(!fpcanocov) for (unsigned long j=i; j<al1.size(); ++j)
                {
                    pcaxx(3*i,3*j)+=al1[i].x*al1[j].x;
                    pcaxx(3*i,3*j+1)+=al1[i].x*al1[j].y;
                    pcaxx(3*i,3*j+2)+=al1[i].x*al1[j].z;
                    pcaxx(3*i+1,3*j+1)+=al1[i].y*al1[j].y;
                    pcaxx(3*i+1,3*j+2)+=al1[i].y*al1[j].z;
                    pcaxx(3*i+2,3*j+2)+=al1[i].z*al1[j].z;
                }
            }
        }

        if (fcv)
        {
            if (!fhavecell) ERROR("You must provide cell data in order to compute g(r)!");

            al1.clear();
            if (lcv=="*") al1=af.ats;
            else for (unsigned long i=0; i<af.ats.size(); ++i)
                if (chk_sel(af.ats[i].name, lcv)) al1.push_back(af.ats[i]);

            //CV stuff
            std::valarray<double> cvv(al1.size());
            for (unsigned long i=0; i<al1.size(); ++i)
            {
                cvv[i]=get_cv(al1,i,cvpars,CM,ICM,cvtype);
            }
            (*ocv) << al1.size()<<"\nOutput from trajworks colvar run\n";
            for (unsigned long i=0; i<al1.size(); ++i)
            {
                (*ocv)<<al1[i].name<<" "<<al1[i].x<<"  "<<al1[i].y<<"  "<<al1[i].z<<"  "<<cvv[i]<<std::endl;
            }
        }

        if (fpd)
        {
            if (npdat==0)
            {
                fpd_inc.resize(af.ats.size()); fpd_inc=false;
		for(unsigned long i=0; i<af.ats.size(); ++i)  if (lpdat=="*" || chk_sel(af.ats[i].name,lpdat) ) { fpd_inc[i]=true; npdat++; }
            }
            FMatrix<double> dpab(3,3,0.);
            for (unsigned long m=0; m<af.ats.size(); ++m) if (fpd_inc[m])
            {
                if (af.ats[m].props.size()<3) ERROR("velocity data not available for atom "<<m<<".");

                for (int i=0; i<3; i++) for (int j=0; j<=i; j++)  pdtens2(i,j)+=(dpab(i,j)=af.ats[m].props[i]*af.ats[m].props[j]);
                for (int i=0; i<3; i++) for (int j=0; j<=i; j++) for (int k=0; k<=j; k++) for (int l=0; l<=k; l++)
                    pdtens4(i,j)(k,l)+=dpab(i,j)*dpab(k,l);
                if (pdvec.size()==0)
                  hpd<<sqrt(af.ats[m].props[0]*af.ats[m].props[0]+af.ats[m].props[1]*af.ats[m].props[1]+af.ats[m].props[2]*af.ats[m].props[2]);
                else
                { hpd<<fabs(af.ats[m].props[0]*pdvec[0]+af.ats[m].props[1]*pdvec[1]+af.ats[m].props[2]*pdvec[2]); }
            }
        }
        if (fpproj)
        {
            //!quick-and-dirty hack. must be optimized sometimes....
            for (unsigned long m=0; m<af.ats.size(); ++m)
            {

                double ppdist, ppmin; unsigned long ipp=0;
                if (lppat1=="*" || chk_sel(af.ats[m].name,lppat1))
                {
                ipp=m; ppmin=1e100;
                for (unsigned long p=0; p<af.ats.size(); ++p)
                {
                    if (p==m) continue;
                    if (lppat2=="*" || chk_sel(af.ats[p].name,lppat2) )
                    {
                        dx=af.ats[m].x-af.ats[p].x;
                        dy=af.ats[m].y-af.ats[p].y;
                        dz=af.ats[m].z-af.ats[p].z;
                        micmat(ICM,dx,dy,dz);
                        micpbc(1.,1.,1.,dx,dy,dz);
                        micmat(CM,dx,dy,dz);
                        ppdist=dx*dx+dy*dy+dz*dz;
                        if (ppdist<ppmin) {ppmin=ppdist; ipp=p; }
                    }
                }
                if (ipp==m) continue;
                double px,py,pz,rx,ry,rz;
                px=af.ats[m].props[0]; py=af.ats[m].props[1]; pz=af.ats[m].props[2];
                dx=af.ats[m].x-af.ats[ipp].x;
                dy=af.ats[m].y-af.ats[ipp].y;
                dz=af.ats[m].z-af.ats[ipp].z;
                micmat(ICM,dx,dy,dz);
                micpbc(1.,1.,1.,dx,dy,dz);
                micmat(CM,dx,dy,dz);
                ppdist=sqrt(dx*dx+dy*dy+dz*dz);
                dx*=1./ppdist; dy*=1./ppdist; dz*=1./ppdist;
                if (af.ats[m].props.size()<3) ERROR("velocity data not available for atom "<<m<<".");

                ppmin=px*dx+py*dy+pz*dz;
                (*opproj)<<ppmin<<" ";
                px-=dx*ppmin; py-=dy*ppmin; pz-=dz*ppmin;

                //project the remaining on two random directions
                rx=rng(); ry=rng(); rz=rng();
                ppmin=rx*dx+ry*dy+rz*dz;
                rx-=dx*ppmin; ry-=dy*ppmin; rz-=dz*ppmin;
                ppmin=1./sqrt(rx*rx+ry*ry+rz*rz);
                rx*=ppmin; ry*=ppmin; rz*=ppmin;

                ppmin=px*rx+py*ry+pz*rz;
                (*opproj)<<ppmin<<" "<<sqrt(px*px+py*py+pz*pz-ppmin*ppmin)*(rng()>0.?1.:-1.)<<"\n";
                }
            }
        }
        if (fp3d)
        {
            if (np3dat==0)
            {
                fp3d_inc.resize(af.ats.size()); fp3d_inc=false;
                for(unsigned long i=0; i<af.ats.size(); ++i)  if (lp3dat=="*" || chk_sel(af.ats[i].name,lp3dat) ) { fp3d_inc[i]=true; np3dat++; }
                if (np3dat==0) ERROR("No atoms matching "<<lp3dat<<" have been found");
            }
            std::valarray<double> p3data(3);
            for (unsigned long i=0; i<af.ats.size(); ++i) if (fp3d_inc[i])
            {
                if (af.ats[i].props.size()<3) ERROR("velocity data not available for atom "<<i<<".");

                p3data[0]=af.ats[i].props[0]; p3data[1]=af.ats[i].props[1]; p3data[2]=af.ats[i].props[2];
                p3dh<<p3data;
                if (fp3dinv) {p3data*=-1.; p3dh<<p3data; } //clearly, it would be better to symmetrize the histogram AFTERWARDS.
            }
        }

        if (fgdr)
        {
            if (!fhavecell) ERROR("You must provide cell data in order to compute g(r)!");
            al1.clear(); al2.clear();
            if (lgdr1=="*") al1=af.ats;
            else for (unsigned long i=0; i<af.ats.size(); ++i)
                if (af.ats[i].name==lgdr1) al1.push_back(af.ats[i]);
            if (lgdr2=="*") al2=af.ats;
            else for (unsigned long i=0; i<af.ats.size(); ++i)
                if (af.ats[i].name==lgdr2) al2.push_back(af.ats[i]);
            
            gdrw=0.0;
            for (unsigned long i=0; i<al1.size(); ++i)
                for (unsigned long j=0; j<al2.size(); ++j)
            {
                if (i==j) continue;
                gdrw+=1.0;
                
                dx=al1[i].x-al2[j].x;
                dy=al1[i].y-al2[j].y;
                dz=al1[i].z-al2[j].z;
                micmat(ICM,dx,dy,dz);
                micpbc(1.,1.,1.,dx,dy,dz);
                micmat(CM,dx,dy,dz);
                
                if (dx>cogdr || dy> cogdr || dz>cogdr) continue;
                d12=dx*dx+dy*dy+dz*dz;
                if (d12>cog2||d12==0.) continue;
                hgdr.add(sqrt(d12),statweight);  // this has the possibility of being weighted
            }
            if (gdrw>0.0) gdrwtot+=statweight*gdrw;   // total weight to be considered when renormalizing g(r)
        }
        if (fvvac)
        {
            if (fvbox) ERROR("OOPS... Variable box is not implemented for this calculation...."); //PLACEHOLDER uhm, box shouldn't really matters here....
            if (vvnat==0)
            {
                fvvac_inc.resize(af.ats.size()); fvvac_inc=false;
                if (vvindex.size()==0)
                {  for(unsigned long i=0; i<af.ats.size(); ++i)  if (lvvac=="*" || chk_sel(af.ats[i].name,lvvac) ) { fvvac_inc[i]=true; vvnat++; } }
                else
                  for(unsigned long j=0; j<vvindex.size(); ++j)
                  {
                      if(vvindex[j]<0 || vvindex[j]>=af.ats.size()) ERROR("Provided atom index "<<vvindex[j]<<" is out of bounds.");
                      if (fvvac_inc[vvindex[j]]==false) {fvvac_inc[vvindex[j]]=true; vvnat++; }
                  }
                std::cerr<<vvindex.size()<<" indices, "<<vvnat<<" atoms selected for acf\n";
                //init ACF array
                vvacf.resize(vvnat*3);
                ACOptions<AutoCorrelation<double> > aco; vvacf[0].get_options(aco); vvacf[0].reset(vvlag);
                aco.f_exact_mean=true; aco.xmean=0.; vvacf[0].set_options(aco);
                for(unsigned long i=0; i<vvnat*3; ++i) vvacf[i]=vvacf[0];
            }
            //else if (vvnat!=af.ats.size()) ERROR("v-v autocorrelation cannot handle varying atomic number");
            unsigned long i=0;
            for (unsigned long j=0; j<af.ats.size(); ++j) if (fvvac_inc[j])
            {
                if (af.ats[j].props.size()<3) ERROR("velocity data not available for atom "<<i<<".");
                vvacf[3*i]<<af.ats[j].props[0];
                vvacf[3*i+1]<<af.ats[j].props[1];
                vvacf[3*i+2]<<af.ats[j].props[2];
                i++;
            }
        }
        if (fmsd)
        {
            
            if (npfr==1){ // mark the species to be used
                          // to compute the MSD
                fmsd_inc.resize(msdlag,af.ats.size()); fmsd_inc.all()=false;
            }
            
            fmsd_inc.row((npfr-1)%msdlag)=false;
            if (lmsd=="*") { fmsd_inc.row((npfr-1)%msdlag)=true; }
            else for(unsigned long i=0; i<af.ats.size(); ++i)
                if(af.ats[i].name==lmsd) { 
                    fmsd_inc((npfr-1)%msdlag,i)=true;
                }
                
            /*std::valarray<bool> tempor(af.ats.size()); 
            tempor=fmsd_inc.row((npfr-1)%msdlag);
            std::cerr<<tempor<<"\n";*/
            // fill the buffer
            msdbuff[(npfr-1)%msdlag] = af;
            if((npfr-1)<(msdlag-1)) continue;
            
            // compute the avg MSD
            imsd = (npfr-1)-(msdlag-1);
            for (unsigned long j=0; j<msdlag; ++j) {
                double cmsd=0.;
                for (unsigned long o=0; o<af.ats.size(); ++o)
                    if(fmsd_inc(imsd%msdlag,o)) {
                       dx=msdbuff[(imsd+j)%msdlag].ats[o].x-msdbuff[imsd%msdlag].ats[o].x;
                       dy=msdbuff[(imsd+j)%msdlag].ats[o].y-msdbuff[imsd%msdlag].ats[o].y;
                       dz=msdbuff[(imsd+j)%msdlag].ats[o].z-msdbuff[imsd%msdlag].ats[o].z;
                       //cmsd+=(dx*dx+dy*dy+dz*dz)/msdnat[imsd%msdlag];
                       dmsd[j]+=dx*dx+dy*dy+dz*dz;
                       nmsd[j]++;
                    }                
                }
            
        }
        if (fdipole)
        {
            //if (fvbox) ERROR("OOPS... Variable box is not implemented for this calculation...."); //PLACEHOLDER
            //if (!fhavecell) ERROR("You must provide cell data in order to compute cell dipole!");
            tcharge=dip_tx.x=dip_tx.y=dip_tx.z=dip_cur.x=dip_cur.y=dip_cur.z=0.;
            for (unsigned long i=0; i<af.ats.size(); ++i)
            {
                dx=af.ats[i].x;
                dy=af.ats[i].y;
                dz=af.ats[i].z;
                /* MUST NOT APPLY PBCs!!!! au contraire, this breaks it all!
                micmat(ICM,dx,dy,dz);
                micpbc(1.,1.,1.,dx,dy,dz);
                micmat(CM,dx,dy,dz);
                */
                // hard-coded atomic charges
                if (af.ats[i].nprops.count("charge") ==0)
                {
                   if (af.ats[i].name == "H") af.ats[i].nprops["charge"]=1.0;
                   if (af.ats[i].name == "O") af.ats[i].nprops["charge"]=-2.0;
                } 
                dip_cur.x+=dx*af.ats[i].nprops["charge"];
                dip_cur.y+=dy*af.ats[i].nprops["charge"];
                dip_cur.z+=dz*af.ats[i].nprops["charge"];
                dip_tx.x+=dx;
                dip_tx.y+=dy;
                dip_tx.z+=dz;
                tcharge+=af.ats[i].nprops["charge"];
            }
            dip_cur.x-=tcharge*dip_tx.x;
            dip_cur.y-=tcharge*dip_tx.y;
            dip_cur.z-=tcharge*dip_tx.z;
            ldip.push_back(dip_cur);
        }

    }

    //normalize the histogram according to g(r) definition

    if (fgdr)
    {
        
        std::cerr<<"# PRINTING OUT g(r) FOR "<<lgdr1<<" - "<<lgdr2<<" PAIR\n";
        std::valarray<double> r, w, gr;
        hgdr.get_bins(r, w, gr);
        for (unsigned long i=0; i<r.size(); ++i)
        { gr[i]=gr[i]/(4*constant::pi*r[i]*r[i]); }
        
        //bins = bins/ndata -> Int bins = 1
        //!!CHECK NORMALIZATION IN CASE SAME SPECIES ARE USED!!
        if (cvolume == 0.) gr*=4./3.*constant::pi*cogdr*cogdr*cogdr;
        else gr*=hgdr.samples()/(gdrwtot/cvolume);
        for (unsigned long i=0; i<r.size(); ++i)
            (*ogdr)<<r[i]<<" "<<gr[i]<<std::endl;
    }
    if (fpd)
    {
        std::cerr<<"# PRINTING OUT PDF(p) FOR "<<lpdat<<" ATOM\n";
        std::valarray<double> p, w, pd;
        //prints out momenta of the distribution
        pdtens2*=(1./nfr)/npdat; pdtens4*=(1./npfr)/npdat;
        (*opd)<<" # SECOND-ORDER MOMENTS <p_i p_j>\n";
        for (int i=0; i<3; i++) { (*opd)<<" # "; for (int j=0; j<3; j++) (*opd)<<pdtens2(i>j?i:j,i>j?j:i)<<" "; (*opd)<<std::endl; }
        (*opd)<<" # FOURTH-ORDER MOMENTS <p_i p_j p_k p_l>-[<p_i p_j><p_k p_l+<p_i p_k><p_j p_l>+<p_i p_l><p_j p_k>]\n";
        std::vector<int> iii(4);
        for (int i=0; i<3; i++) for (int j=0; j<=i; j++) {
            (*opd)<<" # i= "<<i<<" j= "<<j<<std::endl;
            for (int k=0; k<3; k++) {
                (*opd)<<" # ";
                //ack! awful sort!
                for (int l=0; l<3; l++) {
                    iii[0]=i; iii[1]=j; iii[2]=(k>l?k:l); iii[3]=(k>l?l:k);
                    std::sort(iii.begin(),iii.end());
                    (*opd)<<pdtens4(iii[3],iii[2])(iii[1],iii[0])-
                            (pdtens2(iii[3],iii[2])*pdtens2(iii[1],iii[0])+
                             pdtens2(iii[3],iii[1])*pdtens2(iii[2],iii[0])+
                             pdtens2(iii[3],iii[0])*pdtens2(iii[2],iii[1]))<<" ";
                } (*opd)<<std::endl;  }}


        hpd.get_bins(p, w, pd);
        for (unsigned long i=0; i<p.size(); ++i)
            (*opd)<<p[i]<<" "<<pd[i]<<std::endl;
    }
    if (fp3d)
    {
        std::cerr<<"# PRINTING OUT P3DF(p) FOR "<<lp3dat<<" ATOM\n";


        //prints out in CUBE format
        (*op3d)<<"# PRELIMINARY P3DF(p) OUTPUT"
                <<"\n# Gaussian CUBE density format\n";
        //input is in A, but cube is always in bohr
        double cx,dx;
        cx=p3dmax*(-1.0+0.5/p3dbins); dx=p3dmax*2.0/p3dbins;
        (*op3d)<<1<<" "<<cx<<" "<<cx<<" "<<cx<<"\n";
        //remember CM and ICM have x=2 z=0
        (*op3d)<<p3dbins<<"  "<<dx<<" 0.0 0.0 \n";
        (*op3d)<<p3dbins<<" 0.0 "<<dx<<" 0.0 \n";
        (*op3d)<<p3dbins<<" 0.0 0.0 "<<dx<<" \n";
        (*op3d)<<"1 0.0 0.0 0.0 0.0 \n";
        std::valarray<double> nb;
        p3dh.get_bins(nb);
        int dozstep=p3dbins*p3dbins, doshift;
        //cube must be written out with z first
        for (int ix=0; ix<p3dbins; ++ix)  for (int iy=0; iy<p3dbins; ++iy)
        {
            doshift=ix+p3dbins*iy;
            for (int iz=0; iz<p3dbins; ++iz)
            {
                (*op3d)<<nb[doshift+iz*dozstep]<<"  ";
                if (iz%6==5) { (*op3d)<<std::endl; }
            }
            if (dbinsz%6!=0) (*op3d)<<std::endl;
        }
        (*op3d)<<std::endl;
    }
    if (fvvac)
    {

        std::valarray<double> vvt(vvlag+ftpad), vvw(vvlag+ftpad), vbt(vvlag+ftpad);
        vvt=0.; vvw=0.;
        for (unsigned long it=0; it<vvlag; ++it)
        {
            for (unsigned long i=0; i<vvnat; ++i) if (vvacf[3*i].sigma()>0.0) 
                vvt[it]+=vvacf[3*i][it]*vvacf[3*i].sigma()*vvacf[3*i].sigma()
                        +vvacf[3*i+1][it]*vvacf[3*i+1].sigma()*vvacf[3*i+1].sigma()
                        +vvacf[3*i+2][it]*vvacf[3*i+2].sigma()*vvacf[3*i+2].sigma();
        }
        vvt/=vvnat;
        vbt=vvt;
        if (fvvacbox) for (unsigned long it=0; it<vvlag; ++it) vbt[it]=vvt[it]*(1.-double(it)/double(vvlag));
        std::cerr<<"# PRINTING OUT v-v autocorrelation"<<(fvvacbox?" (triangle windowed) ":"")<<"\n";
        //we are goood guys, so we also compute the FT of the vvac straight away
        fftw_plan r2rplan=fftw_plan_r2r_1d(vvlag+ftpad, &vbt[0], &vvw[0],
                                           FFTW_REDFT00, FFTW_ESTIMATE);
        fftw_execute(r2rplan);
        vvw*=dt/(2.*constant::pi);
        for (unsigned long it=0; it<vvlag+ftpad; ++it)
            (*ovvac) <<it*dt<<"  "<<vvt[it]<<"  "<<constant::pi*it/(dt*(vvlag+ftpad-1))<<"  "<<vvw[it]<<std::endl;
    }
    if (fmsd)
    {
        // finishing to average what I can still average
        for (unsigned long ts=0; ts<msdlag; ++ts) {
            ++imsd;
            std::cerr<<  "Finishing to average the MSD "<<std::setw(10)<<std::setiosflags(std::ios::right)<<imsd<<"\r";
            for (unsigned long j=0; j<msdlag; ++j) {
                    if ((imsd+j)>(npfr-1)) break;
                    for (unsigned long o=0; o<af.ats.size(); ++o)
                        if(fmsd_inc(imsd%msdlag,o)) {
                           dx=msdbuff[(imsd+j)%msdlag].ats[o].x-msdbuff[imsd%msdlag].ats[o].x;
                           dy=msdbuff[(imsd+j)%msdlag].ats[o].y-msdbuff[imsd%msdlag].ats[o].y;
                           dz=msdbuff[(imsd+j)%msdlag].ats[o].z-msdbuff[imsd%msdlag].ats[o].z;
                           //cmsd+=(dx*dx+dy*dy+dz*dz)/msdnat[imsd%msdlag];
                           dmsd[j]+=dx*dx+dy*dy+dz*dz;
                           nmsd[j]++;
                        }                
                    
                }
        }
        
        std::cerr<<"# PRINTING OUT msd\n";
        for (unsigned long it=0; it<msdlag; ++it)
        {
            (*omsd) <<it*dt<<"  "<<dmsd[it]/nmsd[it]<<std::endl;
        }
        
    }
    if (fdipole)
    {
        std::cerr<<"# PRINTING OUT dipole\n";
        for (unsigned long it=0; it<ldip.size(); ++it)
        {
            (*odipole) <<it*dt<<"  "<<ldip[it].x<<"  "<<ldip[it].y<<"  "<<ldip[it].z<<std::endl;
        }
    }
     if (fthermal)
     {
         (*otherm) <<"Thermal ellipsoids: Uxx Uyy Uzz Uxy Uxz Uyz\n";
         for (int i=0; i<te_u0.rows(); i++)
         {
              (*otherm)<< te_uiuj(i,0)/npfr-te_u0(i,0)*te_u0(i,0)/(npfr*npfr)<< "   "
                       << te_uiuj(i,1)/npfr-te_u0(i,1)*te_u0(i,1)/(npfr*npfr)<< "   "
                       << te_uiuj(i,2)/npfr-te_u0(i,2)*te_u0(i,2)/(npfr*npfr)<< "   "
                       << te_uiuj(i,3)/npfr-te_u0(i,0)*te_u0(i,1)/(npfr*npfr)<< "   "
                       << te_uiuj(i,4)/npfr-te_u0(i,0)*te_u0(i,2)/(npfr*npfr)<< "   "
                       << te_uiuj(i,5)/npfr-te_u0(i,1)*te_u0(i,2)/(npfr*npfr)<< std::endl;
         }
         otherm->flush();
     }
    if (fcharge)
    {
        std::cerr<<"# PRINTING OUT CHARGE density\n";
        std::valarray<double> nb;
        double a2b1=(1./BOHR2ANG), a2b2=a2b1*a2b1, a2b3=a2b2*a2b1;
        double svolume=cvolume*(drangebx-drangeax)*(drangeby-drangeay)*(drangebz-drangeaz);

        double xvoxel=(drangebx-drangeax)*(drangeby-drangeay)*(drangebz-drangeaz)/(dfoldx*dfoldy*dfoldz)/(dbinsx*dbinsy*dbinsz);
        double vvoxel=svolume/(dfoldx*dfoldy*dfoldz)/(dbinsx*dbinsy*dbinsz);

        //volume of a voxel!
        ndh.get_bins(nb);
        double dfout; ndh.get_outliers(dfout);
        std::cerr<<ndh.samples()*(1.-dfout)/npfr<<" smpl inside\n"<<dfout<<" outliers\n";
        std::cerr<<svolume/(ndh.samples()*(1.-dfout)/npfr)<<" vol per at\n";
        std::cerr<<nb.sum()/nb.size()<<" integral\n"<<nb.size()<<" che oh\n";

        nb*=ndh.samples()/npfr*xvoxel/(vvoxel*a2b3)* (svolume*a2b3)/(svolume*a2b3);
        //nb*=(drangebx-drangeax)*(drangeby-drangeay)*(drangebz-drangeaz)/(dfoldx*dfoldy*dfoldz)*ndh.samples()/npfr /(vvoxel*nb.size()); //normalization so that int dx dy dz v(x,y,z)=N
        double dmax=nb.max(), dmin=nb.min();

        //prints out in CUBE format
        (*ocharge)<<"# VOXEL VOLUME: (bohr^3) " <<vvoxel*a2b3
                << " DENSITY (bohr^-3) RANGE : ["<<dmin<< ".."<<dmax
                << "] AVERAGE: "<<nb.sum()/nb.size()
                <<"\n# Gaussian CUBE density format\n";
        //input is in A, but cube is always in bohr
        FMatrix<double>SCM(CM); SCM*=(1./BOHR2ANG);
        double cx,cy,cz;
        cx=(drangeax+(drangebx-drangeax)*0.5/dbinsx)/dfoldx;
        cy=(drangeay+(drangeby-drangeay)*0.5/dbinsy)/dfoldy;
        cz=(drangeaz+(drangebz-drangeaz)*0.5/dbinsz)/dfoldz;
        micmat(SCM,cx,cy,cz);

        (*ocharge)<<cubefr.ats.size()<<" "<<cx<<" "<<cy<<" "<<cz<<"\n";
        //remember CM and ICM have x=2 z=0
        (*ocharge)<<dbinsx<<"  "
                <<SCM(0,0)*(drangebx-drangeax)/dbinsx/dfoldx<<" "
                <<SCM(1,0)*(drangebx-drangeax)/dbinsx/dfoldx<<" "
                <<SCM(2,0)*(drangebx-drangeax)/dbinsx/dfoldx<<"\n";
        (*ocharge)<<dbinsy<<"  "
                <<SCM(0,1)*(drangeby-drangeay)/dbinsy/dfoldy<<" "
                <<SCM(1,1)*(drangeby-drangeay)/dbinsy/dfoldy<<" "
                <<SCM(2,1)*(drangeby-drangeay)/dbinsy/dfoldy<<"\n";
        (*ocharge)<<dbinsz<<"  "
                <<SCM(0,2)*(drangebz-drangeaz)/dbinsz/dfoldz<<" "
                <<SCM(1,2)*(drangebz-drangeaz)/dbinsz/dfoldz<<" "
                <<SCM(2,2)*(drangebz-drangeaz)/dbinsz/dfoldz<<"\n";
        std::map<std::string, int> cubelab; int curlab=1;
        for (int i=0; i<cubefr.ats.size(); ++i)
        {
            micmat(ICM,cubefr.ats[i].x,cubefr.ats[i].y,cubefr.ats[i].z);
            micpbc(1.,1.,1.,cubefr.ats[i].x,cubefr.ats[i].y,cubefr.ats[i].z);
            micmat(SCM,cubefr.ats[i].x,cubefr.ats[i].y,cubefr.ats[i].z);
            //power of stl::map to number sequentially the species
            if (cubelab.count(cubefr.ats[i].name)==0) cubelab[cubefr.ats[i].name]=(curlab++);
            (*ocharge)<<cubelab[cubefr.ats[i].name] << "  0.0  " << cubefr.ats[i].x<<" "<<cubefr.ats[i].y<<" "<<cubefr.ats[i].z<<"\n";
        }

        int dozstep=dbinsx*dbinsy, doshift;
        //cube must be written out with z first
        for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy)
        {
            doshift=ix+dbinsx*iy;
            for (int iz=0; iz<dbinsz; ++iz)
            {
                (*ocharge)<<nb[doshift+iz*dozstep]<<"  ";
                if (iz%6==5) { (*ocharge)<<std::endl; }
            }
            if (dbinsz%6!=0) (*ocharge)<<std::endl;
        }
        (*ocharge)<<std::endl;

   }
    if (fdens)
    {
        std::cerr<<"# PRINTING OUT 3d density\n";
        std::valarray<double> nb;
        double a2b1=(1./BOHR2ANG), a2b2=a2b1*a2b1, a2b3=a2b2*a2b1;
        double svolume=cvolume*(drangebx-drangeax)*(drangeby-drangeay)*(drangebz-drangeaz);

        double xvoxel=(drangebx-drangeax)*(drangeby-drangeay)*(drangebz-drangeaz)/(dfoldx*dfoldy*dfoldz)/(dbinsx*dbinsy*dbinsz);
        double vvoxel=svolume/(dfoldx*dfoldy*dfoldz)/(dbinsx*dbinsy*dbinsz);

        //volume of a voxel!
        ndh.get_bins(nb);
        double dfout; ndh.get_outliers(dfout);
        std::cerr<<ndh.samples()*(1.-dfout)/npfr<<" smpl inside\n"<<dfout<<" outliers\n";
        std::cerr<<svolume/(ndh.samples()*(1.-dfout)/npfr)<<" vol per at\n";
        std::cerr<<nb.sum()/nb.size()<<" integral\n"<<nb.size()<<" che oh\n";

        nb*=ndh.samples()/npfr*xvoxel/(vvoxel*a2b3)* (svolume*a2b3)/(svolume*a2b3);
        //nb*=(drangebx-drangeax)*(drangeby-drangeay)*(drangebz-drangeaz)/(dfoldx*dfoldy*dfoldz)*ndh.samples()/npfr /(vvoxel*nb.size()); //normalization so that int dx dy dz v(x,y,z)=N
        double dmax=nb.max(), dmin=nb.min();

        //prints out in CUBE format
        (*odens)<<"# VOXEL VOLUME: (bohr^3) " <<vvoxel*a2b3
                << " DENSITY (bohr^-3) RANGE : ["<<dmin<< ".."<<dmax
                << "] AVERAGE: "<<nb.sum()/nb.size()
                <<"\n# Gaussian CUBE density format\n";
        //input is in A, but cube is always in bohr
        FMatrix<double>SCM(CM); SCM*=(1./BOHR2ANG);
        double cx,cy,cz;
        cx=(drangeax+(drangebx-drangeax)*0.5/dbinsx)/dfoldx;
        cy=(drangeay+(drangeby-drangeay)*0.5/dbinsy)/dfoldy;
        cz=(drangeaz+(drangebz-drangeaz)*0.5/dbinsz)/dfoldz;
        micmat(SCM,cx,cy,cz);

        (*odens)<<cubefr.ats.size()<<" "<<cx<<" "<<cy<<" "<<cz<<"\n";
        //remember CM and ICM have x=2 z=0
        (*odens)<<dbinsx<<"  "
                <<SCM(0,0)*(drangebx-drangeax)/dbinsx/dfoldx<<" "
                <<SCM(1,0)*(drangebx-drangeax)/dbinsx/dfoldx<<" "
                <<SCM(2,0)*(drangebx-drangeax)/dbinsx/dfoldx<<"\n";
        (*odens)<<dbinsy<<"  "
                <<SCM(0,1)*(drangeby-drangeay)/dbinsy/dfoldy<<" "
                <<SCM(1,1)*(drangeby-drangeay)/dbinsy/dfoldy<<" "
                <<SCM(2,1)*(drangeby-drangeay)/dbinsy/dfoldy<<"\n";
        (*odens)<<dbinsz<<"  "
                <<SCM(0,2)*(drangebz-drangeaz)/dbinsz/dfoldz<<" "
                <<SCM(1,2)*(drangebz-drangeaz)/dbinsz/dfoldz<<" "
                <<SCM(2,2)*(drangebz-drangeaz)/dbinsz/dfoldz<<"\n";
        std::map<std::string, int> cubelab; int curlab=1;
        for (int i=0; i<cubefr.ats.size(); ++i)
        {
            micmat(ICM,cubefr.ats[i].x,cubefr.ats[i].y,cubefr.ats[i].z);
            micpbc(1.,1.,1.,cubefr.ats[i].x,cubefr.ats[i].y,cubefr.ats[i].z);
            micmat(SCM,cubefr.ats[i].x,cubefr.ats[i].y,cubefr.ats[i].z);
            //power of stl::map to number sequentially the species
            if (cubelab.count(cubefr.ats[i].name)==0) cubelab[cubefr.ats[i].name]=(curlab++);
            (*odens)<<cubelab[cubefr.ats[i].name] << "  0.0  " << cubefr.ats[i].x<<" "<<cubefr.ats[i].y<<" "<<cubefr.ats[i].z<<"\n";
        }

        int dozstep=dbinsx*dbinsy, doshift;
        //cube must be written out with z first
        for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy)
        {
            doshift=ix+dbinsx*iy;
            for (int iz=0; iz<dbinsz; ++iz)
            {
                (*odens)<<nb[doshift+iz*dozstep]<<"  ";
                if (iz%6==5) { (*odens)<<std::endl; }
            }
            if (dbinsz%6!=0) (*odens)<<std::endl;
        }
        (*odens)<<std::endl;

        //if required, also outputs in DF3 format, for POVRAY

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | \
        (((x) << 8) & 0x00ff0000) | \
        (((x) >> 8) & 0x0000ff00) | \
        ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_INT(x)   (*(unsigned int *)&(x)   = SWAP_4(*(unsigned int *)&(x)))

        if (fdpov)
        {
            std::cerr<<"# PRINTING OUT POV-RAY DENSITY FORMAT\n";
            //binary output, I should care about sizes and endiandness, but I am lazy and running on my laptop
            unsigned short int tmp;
            tmp=((unsigned short int) dbinsx); FIX_SHORT(tmp); odpov->write((char *) &tmp,2);
            tmp=((unsigned short int) dbinsy); FIX_SHORT(tmp); odpov->write((char *) &tmp,2);
            tmp=((unsigned short int) dbinsz); FIX_SHORT(tmp); odpov->write((char *) &tmp,2);
            for (int i=0; i<nb.size(); ++i)
            {
                //tmp=((unsigned short int) (65535.*nb[i]/dmax)); odpov->write((char *) &tmp,2);
                //std::cerr<<((unsigned char) (65535.*(nb[i]-dmin)/(dmax-dmin)))<<"\n";
                tmp=((unsigned short int) (65535.*(nb[i]-dmin)/(dmax-dmin))); FIX_SHORT(tmp); odpov->write((char *) &tmp,1);
            }

            odpov->flush();
        }


        if (fdproj)
        {

            //also prints out projected density. since I don't want to have 1000 open files, I reuse odens.
            //in order to make this efficient I should be very careful, but I only want it to WORK, so fuck the cache coherence
#define H3D( ix, iy, iz)  (nb[ix+dbinsx*(iy+dbinsy*iz)])

            double clen, carea, laa, lba, lbb, b2a;
            std::cerr<<"# PRINTING OUT PROJECTED DENSITY\n";

            double pvolume=cvolume*(drangebx-drangeax)/dfoldx*(drangeby-drangeay)/dfoldy*(drangebz-drangeaz)/dfoldz;
            std::valarray<double> lb;

            //!DENS X
            lb.resize(dbinsx); lb=0.;
            for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy) for (int iz=0; iz<dbinsz; ++iz) lb[ix]+=H3D(ix,iy,iz);
            lb*=1./(dbinsy*dbinsz);

            if (odens!=&std::cout)
            {
                odens->flush(); delete odens; odens=new std::ofstream((prefix+std::string(".densx")).c_str());
                (*odens).precision(8); (*odens).width(15); (*odens).setf(std::ios::scientific);
            }
            clen=sqrt(cubefr.nprops["axx"]*cubefr.nprops["axx"]+cubefr.nprops["axy"]*cubefr.nprops["axy"]+cubefr.nprops["axz"]*cubefr.nprops["axz"]);
            vvoxel=clen*(drangebx-drangeax)/dfoldx/dbinsx;

            (*odens)<<"# VOXEL VOLUME: (bohr^1) " <<vvoxel*a2b1
                    << " DENSITY (bohr^-3) RANGE : ["<<lb.min()<< ".."<<lb.max()
                    << "] AVERAGE: "<<lb.sum()/lb.size()
                    <<"\n# Gaussian CUBE density format\n";
            (*odens)<<"0 "<<cx<<std::endl;
            (*odens)<<dbinsx<<"  "
                    <<clen*(drangebx-drangeax)/dbinsx/dfoldx*a2b1<<std::endl;
            for (int ix=0; ix<dbinsx; ++ix)
            {
                (*odens)<<lb[ix]<<"  ";
                if (ix%6==5) { (*odens)<<std::endl; }
            }
            if (dbinsx%6!=0) (*odens)<<std::endl;

            //!DENS Y
            lb.resize(dbinsy); lb=0.;
            for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy) for (int iz=0; iz<dbinsz; ++iz) lb[iy]+=H3D(ix,iy,iz);
            lb*=1./(dbinsx*dbinsz);

            if (odens!=&std::cout)
            {
                odens->flush(); delete odens; odens=new std::ofstream((prefix+std::string(".densy")).c_str());
                (*odens).precision(8); (*odens).width(15); (*odens).setf(std::ios::scientific);
            }
            clen=sqrt(cubefr.nprops["ayx"]*cubefr.nprops["ayx"]+cubefr.nprops["ayy"]*cubefr.nprops["ayy"]+cubefr.nprops["ayz"]*cubefr.nprops["ayz"]);
            vvoxel=clen*(drangeby-drangeay)/dfoldy/dbinsy;

            (*odens)<<"# VOXEL VOLUME: (bohr^1) " <<vvoxel*a2b1
                    << " DENSITY (bohr^-3) RANGE : ["<<lb.min()<< ".."<<lb.max()
                    << "] AVERAGE: "<<lb.sum()/lb.size()
                    <<"\n# Gaussian CUBE density format\n";
            (*odens)<<"0 "<<cy<<std::endl;
            (*odens)<<dbinsy<<"  "
                    <<clen*(drangeby-drangeay)/dbinsy/dfoldy*a2b1<<std::endl;
            for (int iy=0; iy<dbinsy; ++iy)
            {
                (*odens)<<lb[iy]<<"  ";
                if (iy%6==5) { (*odens)<<std::endl; }
            }
            if (dbinsy%6!=0) (*odens)<<std::endl;

            //!DENS Z
            lb.resize(dbinsz); lb=0.;
            for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy) for (int iz=0; iz<dbinsz; ++iz) lb[iz]+=H3D(ix,iy,iz);
            lb*=1./(dbinsx*dbinsy);

            if (odens!=&std::cout)
            {
                odens->flush(); delete odens; odens=new std::ofstream((prefix+std::string(".densz")).c_str());
                (*odens).precision(8); (*odens).width(15); (*odens).setf(std::ios::scientific);
            }
            clen=sqrt(cubefr.nprops["azx"]*cubefr.nprops["azx"]+cubefr.nprops["azy"]*cubefr.nprops["azy"]+cubefr.nprops["azz"]*cubefr.nprops["azz"]);
            vvoxel=clen*(drangebz-drangeaz)/dfoldz/dbinsz;

            (*odens)<<"# VOXEL VOLUME: (bohr^1) " <<vvoxel*a2b1
                    << " DENSITY (bohr^-3) RANGE : ["<<lb.min()<< ".."<<lb.max()
                    << "] AVERAGE: "<<lb.sum()/lb.size()
                    <<"\n# Gaussian CUBE density format\n";
            (*odens)<<"0 "<<cz<<std::endl;
            (*odens)<<dbinsz<<"  "
                    <<clen*(drangebz-drangeaz)/dbinsz/dfoldz*a2b1<<std::endl;
            for (int iz=0; iz<dbinsz; ++iz)
            {
                (*odens)<<lb[iz]<<"  ";
                if (iz%6==5) { (*odens)<<std::endl; }
            }
            if (dbinsz%6!=0) (*odens)<<std::endl;

            //!DENS2 XY
            lb.resize(dbinsx*dbinsy); lb=0.;
            for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy) for (int iz=0; iz<dbinsz; ++iz) lb[iy+dbinsy*ix]+=H3D(ix,iy,iz);
            lb*=1./(dbinsz);

            if (odens!=&std::cout)
            {
                odens->flush(); delete odens; odens=new std::ofstream((prefix+std::string(".densxy")).c_str());
                (*odens).precision(8); (*odens).width(15); (*odens).setf(std::ios::scientific);
            }
            laa=sqrt(cubefr.nprops["axx"]*cubefr.nprops["axx"]+cubefr.nprops["axy"]*cubefr.nprops["axy"]+cubefr.nprops["axz"]*cubefr.nprops["axz"]);
            lba=sqrt(cubefr.nprops["axx"]*cubefr.nprops["ayx"]+cubefr.nprops["axy"]*cubefr.nprops["ayy"]+cubefr.nprops["axz"]*cubefr.nprops["ayz"])/laa;
            lbb=sqrt(cubefr.nprops["ayx"]*cubefr.nprops["ayx"]+cubefr.nprops["ayy"]*cubefr.nprops["ayy"]+cubefr.nprops["ayz"]*cubefr.nprops["ayz"]-lba*lba);
            carea=laa*lbb;
            vvoxel=carea*(drangebx-drangeax)/dfoldx/dbinsx*(drangeby-drangeay)/dfoldy/dbinsy;


            (*odens)<<"# VOXEL VOLUME: (bohr^2) " <<vvoxel*a2b2
                    << " DENSITY (bohr^-3) RANGE : ["<<lb.min()<< ".."<<lb.max()
                    << "] AVERAGE: "<<lb.sum()/lb.size()
                    <<"\n# Gaussian CUBE density format\n";
            (*odens)<<"0 "<<cx<<" "<<cy<<std::endl;
            (*odens)<<dbinsx<<"  "<<laa*(drangebx-drangeax)/dbinsx/dfoldx*a2b1<<"  0."<<std::endl;
            (*odens)<<dbinsy<<"  "<<lba*(drangeby-drangeay)/dbinsy/dfoldy*a2b1<<"  "<<lbb*(drangeby-drangeay)/dbinsy/dfoldy*a2b1<<std::endl;
            for (int ix=0; ix<dbinsx; ++ix)
            {
                for (int iy=0; iy<dbinsy; ++iy)
                {
                    (*odens)<<lb[iy+dbinsy*ix]<<"  ";
                    if (iy%6==5) { (*odens)<<std::endl; }
                }
                if (dbinsy%6!=0) (*odens)<<std::endl;
            }

            //!DENS2 YZ
            lb.resize(dbinsy*dbinsz); lb=0.;
            for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy) for (int iz=0; iz<dbinsz; ++iz) lb[iz+dbinsz*iy]+=H3D(ix,iy,iz);
            lb*=1./(dbinsx);

            if (odens!=&std::cout)
            {
                odens->flush(); delete odens; odens=new std::ofstream((prefix+std::string(".densyz")).c_str());
                (*odens).precision(8); (*odens).width(15); (*odens).setf(std::ios::scientific);
            }
            laa=sqrt(cubefr.nprops["ayx"]*cubefr.nprops["ayx"]+cubefr.nprops["ayy"]*cubefr.nprops["ayy"]+cubefr.nprops["ayz"]*cubefr.nprops["ayz"]);
            lba=sqrt(cubefr.nprops["ayx"]*cubefr.nprops["azx"]+cubefr.nprops["ayy"]*cubefr.nprops["azy"]+cubefr.nprops["ayz"]*cubefr.nprops["azz"])/laa;
            lbb=sqrt(cubefr.nprops["azx"]*cubefr.nprops["azx"]+cubefr.nprops["azy"]*cubefr.nprops["azy"]+cubefr.nprops["azz"]*cubefr.nprops["azz"]-lba*lba);
            carea=laa*lbb;
            vvoxel=carea*(drangeby-drangeay)/dfoldy/dbinsy*(drangebz-drangeaz)/dfoldz/dbinsz;


            (*odens)<<"# VOXEL VOLUME: (bohr^2) " <<vvoxel*a2b2
                    << " DENSITY (bohr^-3) RANGE : ["<<lb.min()<< ".."<<lb.max()
                    << "] AVERAGE: "<<lb.sum()/lb.size()
                    <<"\n# Gaussian CUBE density format\n";
            (*odens)<<"0 "<<cy<<" "<<cz<<std::endl;
            (*odens)<<dbinsy<<"  "<<laa*(drangeby-drangeay)/dbinsy/dfoldy*a2b1<<"  0."<<std::endl;
            (*odens)<<dbinsz<<"  "<<lba*(drangebz-drangeaz)/dbinsz/dfoldz*a2b1<<"  "<<lbb*(drangebz-drangeaz)/dbinsz/dfoldz*a2b1<<std::endl;
            for (int iy=0; iy<dbinsy; ++iy)
            {
                for (int iz=0; iz<dbinsz; ++iz)
                {
                    (*odens)<<lb[iz+dbinsz*iy]<<"  ";
                    if (iz%6==5) { (*odens)<<std::endl; }
                }
                if (dbinsz%6!=0) (*odens)<<std::endl;
            }

            //!DENS2 ZX
            lb.resize(dbinsz*dbinsx); lb=0.;
            for (int ix=0; ix<dbinsx; ++ix)  for (int iy=0; iy<dbinsy; ++iy) for (int iz=0; iz<dbinsz; ++iz) lb[ix+dbinsx*iz]+=H3D(ix,iy,iz);
            lb*=1./(dbinsy);

            if (odens!=&std::cout)
            {
                odens->flush(); delete odens; odens=new std::ofstream((prefix+std::string(".denszx")).c_str());
                (*odens).precision(8); (*odens).width(15); (*odens).setf(std::ios::scientific);
            }
            laa=sqrt(cubefr.nprops["azx"]*cubefr.nprops["azx"]+cubefr.nprops["azy"]*cubefr.nprops["azy"]+cubefr.nprops["azz"]*cubefr.nprops["azz"]);
            lba=sqrt(cubefr.nprops["azx"]*cubefr.nprops["axx"]+cubefr.nprops["azy"]*cubefr.nprops["axy"]+cubefr.nprops["azz"]*cubefr.nprops["axz"])/laa;
            lbb=sqrt(cubefr.nprops["axx"]*cubefr.nprops["axx"]+cubefr.nprops["axy"]*cubefr.nprops["axy"]+cubefr.nprops["axz"]*cubefr.nprops["axz"]-lba*lba);
            carea=laa*lbb;
            vvoxel=carea*(drangebz-drangeaz)/dfoldz/dbinsz*(drangebx-drangeax)/dfoldx/dbinsx;


            (*odens)<<"# VOXEL VOLUME: (bohr^2) " <<vvoxel*a2b2
                    << " DENSITY (bohr^-3) RANGE : ["<<lb.min()<< ".."<<lb.max()
                    << "] AVERAGE: "<<lb.sum()/lb.size()
                    <<"\n# Gaussian CUBE density format\n";
            (*odens)<<"0 "<<cz<<" "<<cx<<std::endl;
            (*odens)<<dbinsz<<"  "<<laa*(drangebz-drangeaz)/dbinsz/dfoldz*a2b1<<"  0."<<std::endl;
            (*odens)<<dbinsx<<"  "<<lba*(drangebx-drangeax)/dbinsx/dfoldx*a2b1<<"  "<<lbb*(drangebx-drangeax)/dbinsx/dfoldx*a2b1<<std::endl;
            for (int iz=0; iz<dbinsz; ++iz)
            {
                for (int ix=0; ix<dbinsx; ++ix)
                {
                    (*odens)<<lb[ix+dbinsx*iz]<<"  ";
                    if (ix%6==5) { (*odens)<<std::endl; }
                }
                if (dbinsx%6!=0) (*odens)<<std::endl;
            }
        }
        std::cerr<<"# DENSITY OVER\n";
    }

    if (fpca)
    {
        //first, fills up covariance matrix
        pcax*=(1./npfr);
        if(!fpcanocov)
        {
            pcaxx*=(1./(npfr-1));
            for (unsigned long i=0; i<3*pcasz; ++i)
            {
                pcaxx(i,i)-=pcax[i]*pcax[i];
                for (unsigned long j=0; j<i; ++j)
                    pcaxx(i,j)=(pcaxx(j,i)-=pcax[i]*pcax[j]);
            }
        }
        if (fpcaxyz)
        {
            (*opcaxyz)<<pcasz<<"\nAVERAGE ATOM POSITIONS\n";
            for (unsigned long j=0; j<pcasel.size(); ++j)
            {
                unsigned int i=pcasel[j];
                (*opcaxyz)<<pcafr.ats[i].name<<"  "<<pcax[3*j]<<"  "<<pcax[3*j+1]<<"  "<<pcax[3*j+2];
                //also prints out the diagonal part of the covariance matrix....
                if (!fpcanocov) (*opcaxyz) <<"  "<< pcaxx(3*j,3*j)<<"  "<< pcaxx(3*j+1,3*j+1)<<"  "<< pcaxx(3*j+2,3*j+2);
                (*opcaxyz)<<"\n";
            }
            (*opcaxyz)<<std::endl;
            std::cerr<<"XYZ CENTER FILE WRITTEN\n";
        }

        if (!fpcanocov) {
            FMatrix<double> Q; std::valarray<double> ev;
            EigenSolverSym(pcaxx,Q,ev);

            //NOW WRITES F**ING GAUSSIAN LOG FILE
            (*opca)<<" Entering Gaussian System\nthis file is generated from the MLDGAU subroutine in the file secder.F\nStandard orientation:\n ---------------------------------------------------------------------\nCenter     Atomic     Atomic              Coordinates (Angstroms)\nNumber     Number      Type              X           Y           Z\n---------------------------------------------------------------------"<<std::endl;
            for (unsigned long i=0; i<pcasz; ++i)
               (*opca)<<i+1<<"  "<<1<<"   "<<0<<"    "<<pcax[3*i]<<"  "<<pcax[3*i+1]<<"  "<<pcax[3*i+2]<<std::endl;
            (*opca)<<" ---------------------------------------------------------------------\nbasis functions          primitive gaussians\nalpha electrons          beta electrons\n**********************************************************************"<<std::endl;
            (*opca)<<"\nHarmonic frequencies (cm**-1), IR intensities (KM/Mole),\nRaman scattering activities (A**4/AMU), Raman depolarization ratios,\nreduced masses (AMU), force constants (mDyne/A) and normal coordinates:\n";
            for (unsigned long i=0; i<pcasz; ++i)
            {
                (*opca)<< 3*i<< "  "<<3*i +1 <<"  "<< 3*i+2<<std::endl<<
                        " ?A                     ?A                     ?A\n Frequencies --   "<<
                        ev[3*i]<<"   "<<ev[3*i+1]<<"   "<< ev[3*i+2]<<std::endl<<
                        " Red. masses --     0.0000                 0.0000                 0.0000\n Frc consts  --     0.0000                 0.0000                 0.0000\n IR Inten    --     0.0000                 0.0000                 0.0000\n Raman Activ --     0.0000                 0.0000                 0.0000\nDepolar     --     0.0000                 0.0000                 0.0000\n Atom AN      X      Y      Z        X      Y      Z        X      Y      Z\n";
                for (unsigned long j=0; j<pcasz; ++j)
                {
                    (*opca)<<j<<"  1  "<< Q(3*j,3*i) << "  "<< Q(3*j+1,3*i) <<"  "<<Q(3*j+2,3*i) <<"  "
                            << Q(3*j,3*i+1) << "  "<< Q(3*j+1,3*i+1) <<"  "<<Q(3*j+2,3*i+1) <<"  "
                            << Q(3*j,3*i+2) << "  "<< Q(3*j+1,3*i+2) <<"  "<<Q(3*j+2,3*i+2) <<std::endl;
                }
            }
            (*opca)<< " Normal termination of Gaussian 98."<<std::endl;
        /*for (unsigned long i=0; i<3*pcasz; ++i) (*opca)<<ev[i]<<"\n";
        (*opca)<<"[FREQ]"<<std::endl;
        for (unsigned long i=0; i<3*pcasz; ++i) (*opca)<<ev[i]<<"\n";
        (*opca)<<"[FR-COORD]"<<std::endl;
        for (unsigned long i=0; i<pcasz; ++i)
            (*opca)<<lpcat<<"  "<<pcax[3*i]<<"  "<<pcax[3*i+1]<<"  "<<pcax[3*i+2]<<std::endl;
        (*opca)<<"[FR-NORM-COORD]"<<std::endl;
        for (unsigned long i=0; i<3*pcasz; ++i)
        {
            (*opca)<<"vibration "<<i<<std::endl;
            for (unsigned long j=0; j<pcasz; ++j)
                (*opca)<<Q(3*j,i)<<"  "<<Q(3*j+1,i)<<"  "<<Q(3*j+2,i)<<std::endl;
        }*/
        }
    }
    return 0;
}
