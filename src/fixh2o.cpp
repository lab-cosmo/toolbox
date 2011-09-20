#include "ioatoms.hpp"
#include "ioparser.hpp"
#include "clparser.hpp"
#include <valarray>
using namespace toolbox;

double norm(const std::valarray<double>& v)
{
    double rn=0.;
    for (int i=0; i<v.size(); ++i) rn+=v[i]*v[i];
    return sqrt(rn);
}

double sprod(const std::valarray<double>& u, const std::valarray<double>& v)
{
    double rn=0.;
    for (int i=0; i<v.size(); ++i) rn+=u[i]*v[i];
    return rn;
}

bool drude=true;
double l=0.9572, a=104.52/180.*constant::pi, dm=0.24034; //POLARIZABLE TIP4P (lamo+06cpl)

//bool drude=false;
//double l=0.9572, a=104.25/180.*constant::pi, dm=0.15;  //TIP4P
//double l=0.9572, a=104.25/180.*constant::pi, dm=0.1546;  //TIP4P/2005

int main(int argc, char **argv)
{
    AtomFrame af;
    std::string s;
    if (!ReadDLPConf(std::cin,af)) ERROR("No valid DLP file on stdin");

    std::vector<AtomData>& al=af.ats; 
    std::valarray<double> u(3), v(3), o(3), os(3), d(3), h1(3), h2(3);
    double n; double lu=l*cos(a/2.), lv=l*sin(a/2.);
//    std::cout<<al.size()<<"\n"<<"Converted to tip4p\n"; 
    std::cout.width(16); std::cout.precision(8); std::cout.setf(std::ios::scientific);
    std::cout<<"Water molecules fixed\n"
            <<"  0    3   "<<al.size()<<" 0.001\n";
//            <<"timestep   0  "<<al.size()<<"  0    3   "<<af.nprops["timestep"]<<"\n";
    std::cout<<af.nprops["axx"]<<" "<<af.nprops["axy"]<<" "<<af.nprops["axz"]<<"\n"
             <<af.nprops["ayx"]<<" "<<af.nprops["ayy"]<<" "<<af.nprops["ayz"]<<"\n"
             <<af.nprops["azx"]<<" "<<af.nprops["azy"]<<" "<<af.nprops["azz"]<<"\n";
    int ia=0;
    while(al[ia].name!=std::string("OW")) 
    {
        std::cout<<al[ia].name<<"  "<<1+ia<<"\n"; 
        std::cout<<al[ia].x <<" "<<al[ia].y<<" "<<al[ia].z<<"\n";
        ++ia;
    }
    for (; ia<al.size(); ia+=5)
    {
        o[0]=al[ia].x; o[1]=al[ia].y; o[2]=al[ia].z;
        os[0]=al[ia+1].x; os[1]=al[ia+1].y; os[2]=al[ia+1].z;  os-=o;
        d[0]=al[ia+2].x; d[1]=al[ia+2].y; d[2]=al[ia+2].z;
        h1[0]=al[ia+3].x; h1[1]=al[ia+3].y; h1[2]=al[ia+3].z;
        h2[0]=al[ia+4].x; h2[1]=al[ia+4].y; h2[2]=al[ia+4].z;
        
        u=(h1+h2)*0.5-o;
        v=h2-h1;
        u*=(1./norm(u));
        v-=(sprod(v,u)*u);
        v*=(1./norm(v));
        
        h1=o+lu*u+lv*v;
        h2=o+lu*u-lv*v;
        d=o+dm*u;
        std::cout<<"OW  " <<1+ia<< "\n"<<o[0]<<" "<<o[1]<<" "<<o[2]<<"\n";
        std::cout<<"OWS  "<<1+ ia+1<<"  \n"<<o[0]+os[0]<<" "<<o[1]+os[1]<<" "<<o[2]+os[2]<<"\n";
        std::cout<<"DW  "<<1+ia+2<<" \n"<<d[0]<<" "<<d[1]<<" "<<d[2]<<"\n";
        std::cout<<"HW  "<<1+ia+3<<" \n"<< h1[0]<<" "<<h1[1]<<" "<<h1[2]<<"\n";
        std::cout<<"HW  "<<1+ia+4<<" \n"<< h2[0]<<" "<<h2[1]<<" "<<h2[2]<<"\n";
    }
}
