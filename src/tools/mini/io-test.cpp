#include "tbdefs.hpp"
#include "ioparser.hpp"
#include <complex>
using namespace toolbox;


//this class define a read/write iostream operator. therefore, it can be used with the std method
class myclass {
public:
    double d;
    long l;
    std::string s;
};
std::ostream& operator << (std::ostream& ostr, const myclass& myc)
{
    ostr << myc.d<<" "<<myc.l <<" "<< myc.s;
    return ostr;
}
std::istream& operator >> (std::istream& istr, myclass& myc)
{
    istr >> myc.d >> myc.l >> myc.s;
}
    
//strategy two is to define, together with the class, a specialization of the IField template, 
//so we can do fancier things and, most important, avoid multiple copies of strings in nested fields
class myclass_over {
public:
        double d;
        long l;
        std::string s;
};
//in this case, we just avoid the label .... &end_label syntax
namespace toolbox{
    template<> inline bool IField<myclass_over>::operator<<(std::istream& istr) { istr >> value.d >> value.l >> value.s;  flags|=iff_set; return true; }
    template<> inline bool IField<myclass_over>::operator>>(std::ostream& ostr) const { ostr <<name<<" "<< value.d<<" "<<value.l <<" "<< value.s<<"\n"; 
        return true; }
}

//here we use an IOMap to build a structured input and to avoid multiple parsing issues
class myclass_complex {
    public:
        double d;
        long l;
        std::string s;
        myclass_over mco;
};

namespace toolbox{
    template<> inline bool IField<myclass_complex>::operator<<(std::istream& istr) 
    { 
        IOMap iom;
        iom.insert(value.d, "d"); iom.insert(value.l, "l"); iom.insert(value.s, "s"); iom.insert(value.mco, "mco");
        iom << istr;  flags|=iff_set;
    }
    template<> inline bool IField<myclass_complex>::operator>>(std::ostream& ostr) const { 
        IOMap iom;
        iom.insert(value.d, "d"); iom.insert(value.l, "l"); iom.insert(value.s, "s"); iom.insert(value.mco, "mco");
        ostr <<name<<" { \n";
        iom >> ostr;
        ostr<<" } \n";
    }
}

#include <fstream>
int main(int argc, char **argv)
{
    /*******************************************************************
                  test for IOMap as an input parser
    ********************************************************************/
     
    IOMap iom;
    iom.flags=if_warn_all;   //be overly verbose in input checking
    
    //defines a few variables which will hold the options
    std::string vcolore; double vtemperatura; long vnumero;
    std::valarray<double> vvar; std::complex<double> vcd;
    myclass vmyc; myclass_over vmyc_o; myclass_complex vmyc_c;
    
    //define the structure of the input file.
    //the linked variables will be automatically set as the input is parsed.
    iom.insert(vcolore,std::string("colore"));
    iom.insert(vtemperatura,"temperatura");
    iom.insert(vnumero,"numero");
    iom.insert(vvar,"array");
    //this field is using the default template field IField<std::complex<double>
    iom.insert(vcd,"complex"); 
    //this field is using specialized template
    iom.insert(vmyc_o,"myclass_over"); 
    iom.insert(vmyc_c,"cclass"); 
    
    //nested input field
    IOMap nested; nested.flags=if_warn_all;
    double vn1, vn2;
    std::valarray<std::string> vnsarr;
    nested.insert(vn1, "n1");
    nested.insert(vn2, "n2");
    nested.insert(vnsarr, "stringarray");
    nested.insert(vmyc,"myclass"); 
    
    iom.insert(nested, "nested");
    
    
    //the input is simply parsed this way
    std::ifstream iio("o");
    iom << iio;
    
    //check that variables have been set correctly
    std::cout << "**************************************************************\n" 
            <<   "*             VALUES OF SOME PARSED OPTIONS:                 *\n"
            <<   "**************************************************************\n" 
            << "colore : "<<vcolore << "\n"
            << "numero : "<<vnumero << "\n"
            << "complex : "<<vcd<< "\n";
    std::cout<< "array : ";
    for (int i=0; i<vvar.size(); ++i) std::cout<<vvar[i]<<" ";
    std::cout<< "\n";
    std::cout<<"nested.n1 : "<<vn1<<"\n";
    std::cout<<"nested.stringarray : ";
    for (int i=0; i<vnsarr.size(); ++i) std::cout<<vnsarr[i]<<" ";
    std::cout<< "\n";
    std::cout<<"nested.myclass : "<<vmyc.d<<" "<<vmyc.l<<" "<<vmyc.s<<"\n";
    std::cout << "**************************************************************\n";
     
    
    /*******************************************************************
                   test for IOMap as an output engine
    ********************************************************************/
    //this should be improved, but by now, since an initialized IOMap class 
    //can write out a re-parsable input, one coud build a strategy out of this,
    //to write out the options of various classes
    std::cout << "**************************************************************\n" 
            <<   "*               FULL DUMP OF IOMAP CONTENT:                  *\n"
            <<   "**************************************************************\n"; 
    std::cout <<iom;
    std::cout << "**************************************************************\n";
    
    return 0;
}
