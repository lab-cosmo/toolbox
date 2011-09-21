#include "rndgen.hpp"
using namespace toolbox;

int main(int argc, char **argv)
{
    //MTRndUniform RU;
        
    std::valarray<float> corg(4);
    corg[0]=1.; corg[1]=-0.6; corg[2]=0.2; corg[3]=0.0;
    RndCorrGaussian<float, RndGaussian<float, MTRndUniform> > prova(RCGPars<float>(0,1,corg));
    prova.RGGenerator().RUGenerator().seed(100);
    
    
    
    
    for (int i=0; i<100000; ++i)
    {
        std::cout << prova()<<"\t"<< prova.RGGenerator()()<<"\t"<< prova.RGGenerator().RUGenerator()() <<"\n";
    }
    return -1;
}
