/* A simple parsing library for command line options
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/
   
#include "clparser.hpp"

namespace toolbox{
CLParser::CLParser(int const& argc, char** const& argv, bool nfverbose)
{
    fverbose=nfverbose;
    ppars.resize(argc-1); fparsed.resize(argc-1); //argument zero is the program name. 
    for (long i=1; i<argc;++i) { ppars[i-1]=argv[i]; fparsed[i-1]=false; }
}

bool CLParser::chkunknown(bool print)
{
    bool allok=true;
    for (long i=0; i<fparsed.size();++i)
        if (!fparsed[i]) 
        {  
            allok=false;  
            if (print) std::cerr<<" toolbox::CLParser - Unkown argument "<< ppars[i]<< " on command line."<<std::endl; 
        }
    return allok; 
}
   
// A few specialized version of CLParser::getoption
template<>
bool CLParser::getoption(bool& value, const std::string& parname)
{ 
    //for a bool parameter we just check the existence of the key on the commandline
    std::string spar=std::string("-")+parname;
    unsigned long i=0; value=false;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!
    fparsed[i]=true; // signal the key has been checked
    value=true; return true;
}

template <>
bool CLParser::getoption(std::string& value, const std::string& parname)
{
    //for a string parameter we don't use the stringstream!
    std::string spar=std::string("-")+parname;
    unsigned long i=0;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!

    if (i+1==ppars.size() || ppars[i+1].size()==0) 
    {
        if (fverbose) WARNING((std::string("Key value -")+ parname  +std::string(" missing on command line\n")));
#ifdef DEBUG
        ERROR("Key value missing on command line\n");
#endif
        return false;
    }
    fparsed[i]=fparsed[i+1]=true;
    value=ppars[i+1]; return true;
}

//for double type, we must accept negative numbers!
template <>
bool CLParser::getoption(double& value, const std::string& parname)
{
    std::string spar=std::string("-")+parname;
    unsigned long i=0;
    while (i<ppars.size() && ppars[i]!=spar) ++i;
    if (i==ppars.size()) return false; //key not found!

    if (i+1==ppars.size() || ppars[i+1].size()==0) 
    {
        if (fverbose) WARNING((std::string("Key value -")+ parname  +std::string(" missing on command line\n")));
#ifdef DEBUG
        ERROR((std::string("Key value \"")+parname+std::string("\" missing on command line\n")));
#endif
        return false;
    }
    fparsed[i]=fparsed[i+1]=true;
    
    std::stringstream ls(ppars[i+1]); 
    ls >> value; 
    if (ls.good() || ls.eof()) return true; 
    else return false;
}

};
