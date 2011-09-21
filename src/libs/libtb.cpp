#define __EXTERNALS 1
#include "tbdefs.hpp"

//externals and very simple general-purpose functions

namespace toolbox {
#ifdef BENCHMARK
    CTBBD TBBenchmarks;
#endif
};

namespace std{
    toolbox::mpiostream pout(std::cout), perr(std::cerr);
};


namespace toolbox{
    double str2float(const std::string& str)
    {
        std::stringstream sstr(str);
        double rval; sstr >> rval;
        return rval;
    }

    int str2int(const std::string& str)
    {
        std::stringstream sstr(str);
        int rval; sstr >> rval;
        return rval;
    }
        
    std::string int2str(const long& ival)
    {
        std::stringstream sstr;
        sstr << ival;
        std::string rval; sstr >> rval;
        return rval;
    }
        
    std::string float2str(const double& ival)
    {
        std::stringstream sstr;
        sstr << ival;
        std::string rval; sstr >> rval;
        return rval;
    }
};
