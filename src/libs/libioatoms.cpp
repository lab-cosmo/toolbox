/* Read and write a few atomic data formats.
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License  
*/
   
#include "ioatoms.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "rpc/xdr.h"

#define XTC_MAGIC 1995


namespace toolbox {

bool ReadXYZFrame(std::istream& istr, AtomFrame& curfr)
{
    std::string line;
    std::stringstream ss;
    AtomData curat;
    unsigned long nat, i; double cprop; 
    if(!getline(istr, line)) return false;
    curfr.index=0; curfr.ats.resize(0);
    ss.clear(); ss<<line;
    ss>>nat;
    if (!getline(istr,curfr.comment)) ERROR("Read failed while getting frame "<<curfr.index<<".");
    for (i=0; i<nat && getline(istr, line); ++i)
    {
        ss.clear(); ss<<line;
        curat.props.resize(0);
        ss>>curat.name>>curat.x>>curat.y>>curat.z;
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        while ((ss>>cprop)) curat.props.push_back(cprop);
        curfr.ats.push_back(curat);
    }
    return true;
}

bool ReadPDBFrame(std::istream& istr, AtomFrame& curfr)
{
    std::string line, check, atlabel, none1, none2, start;
    std::stringstream ss;
    AtomData curat;
    unsigned long nat, ifr, ipbc, ikey, i; double cprop;
    if(!getline(istr, line)) return false;
    curfr.index=0; curfr.ats.resize(0);
    //ss.clear(); ss<<line;
    //ss>>start;
    //if (!getline(istr,curfr.comment)) ERROR("Read failed while getting frame "<<curfr.index<<".");
    for (i=0; getline(istr, line); ++i)
    {   
        ss.clear(); ss.str(line);
        curat.props.resize(0);
        ss>>check>>ikey>>atlabel>>curat.molname>>curat.nprops["mol"]>>curat.x>>curat.y>>curat.z>>none1>>none2>>curat.group>>curat.name;
        if(check == "END") continue;
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        while ((ss>>cprop)) curat.props.push_back(cprop);
        curfr.ats.push_back(curat);
    }
    return true;
}

void ReadXYZ(std::istream& istr, std::vector<AtomFrame>& frames)
{
    frames.resize(0);
    AtomFrame curfr;
    unsigned long nfr=0;
    while(ReadXYZFrame(istr,curfr)) {curfr.index=(++nfr); frames.push_back(curfr);}
}

void ReadPDB(std::istream& istr, std::vector<AtomFrame>& frames)
{
    frames.resize(0);
    AtomFrame curfr;
    unsigned long nfr=0;
    while(ReadPDBFrame(istr,curfr)) {curfr.index=(++nfr); frames.push_back(curfr);}
}
    
bool ReadDLPFrame(std::istream& istr, AtomFrame& curfr)
{
    std::string line, dummy;
    std::stringstream ss;
    AtomData curat;
    unsigned long nat, ifr, ipbc, ikey, i; double cprop; 
    if(!getline(istr, line)) return false;
    curfr.index=0; curfr.ats.resize(0);
    ss.clear(); ss<<line;
    ss>>dummy>>ifr>>nat>>ikey>>ipbc>>curfr.nprops["timestep"];
    if(ipbc>0)
    {
        istr>>curfr.nprops["axx"]>>curfr.nprops["axy"]>>curfr.nprops["axz"];
        istr>>curfr.nprops["ayx"]>>curfr.nprops["ayy"]>>curfr.nprops["ayz"];
        istr>>curfr.nprops["azx"]>>curfr.nprops["azy"]>>curfr.nprops["azz"];
        getline(istr,dummy);
    }
    for (i=0; i<nat && getline(istr, line); ++i)
    {
        ss.clear(); ss.str(line);
        curat.props.resize(0);
        ss>>curat.name>>dummy>>curat.nprops["mass"]>>curat.nprops["charge"];
        getline(istr, line); ss.clear(); ss.str(line);
        ss>>curat.x>>curat.y>>curat.z;
        
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>0)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>1)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        curfr.ats.push_back(curat);
    }
    return true;
}
    
bool ReadDLPConf(std::istream& istr, AtomFrame& curfr)
{
    std::string line, dummy;
    std::stringstream ss;
    AtomData curat;
    unsigned long nat, ifr, ipbc, ikey, i; double cprop; 
    if(!getline(istr, line)) return false;  //drops comment line
    if(!getline(istr, line)) return false;  //reads header
    curfr.index=0; curfr.ats.resize(0);
    ss.clear(); ss<<line;
    ss>>ikey>>ipbc>>nat>>dummy;
    if(ipbc>0)
    {
        istr>>curfr.nprops["axx"]>>curfr.nprops["axy"]>>curfr.nprops["axz"];
        istr>>curfr.nprops["ayx"]>>curfr.nprops["ayy"]>>curfr.nprops["ayz"];
        istr>>curfr.nprops["azx"]>>curfr.nprops["azy"]>>curfr.nprops["azz"];
        getline(istr,dummy);
    }
    nat=0;
    for (i=0; getline(istr, line); ++i)
    {
        nat++;
        ss.clear(); ss.str(line);
        curat.props.resize(0);
        ss>>curat.name>>dummy;
        getline(istr, line); ss.clear(); ss.str(line);
        ss>>curat.x>>curat.y>>curat.z;
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>0)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        if (ikey>1)
        {
            getline(istr, line); ss.clear(); ss.str(line);
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop); 
            ss>>cprop; curat.props.push_back(cprop);
        }
        if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        curfr.ats.push_back(curat);
    }
    return true;
}

bool WritePDBFrame(std::ostream& ostr, AtomFrame& frame)
{
    double lx, ly, lz, ca, cb, cc; 
    if (frame.nprops.count("axx")>0) 
    {
        lx=sqrt(frame.nprops["axx"]*frame.nprops["axx"]+frame.nprops["axy"]*frame.nprops["axy"]+frame.nprops["axz"]*frame.nprops["axz"]);
        ly=sqrt(frame.nprops["ayx"]*frame.nprops["ayx"]+frame.nprops["ayy"]*frame.nprops["ayy"]+frame.nprops["ayz"]*frame.nprops["ayz"]);
        lz=sqrt(frame.nprops["azz"]*frame.nprops["azx"]+frame.nprops["azy"]*frame.nprops["azy"]+frame.nprops["azz"]*frame.nprops["azz"]);
        cc=acos((frame.nprops["axx"]*frame.nprops["ayx"]+frame.nprops["axy"]*frame.nprops["ayy"]+frame.nprops["axz"]*frame.nprops["ayz"])/(lx*ly))*180/constant::pi;
        cb=acos((frame.nprops["axx"]*frame.nprops["azx"]+frame.nprops["axy"]*frame.nprops["azy"]+frame.nprops["axz"]*frame.nprops["azz"])/(lx*ly))*180/constant::pi;;
        ca=acos((frame.nprops["azx"]*frame.nprops["ayx"]+frame.nprops["azy"]*frame.nprops["ayy"]+frame.nprops["azz"]*frame.nprops["ayz"])/(lx*ly))*180/constant::pi;;
    }
    else
    {
        lx=ly=lz=0.0; 
        ca=cb=cc=90.0;
    }
    
    //outputs cell parameters (if present)
    ostr<<"CRYST1"
            <<std::setw(9)<<std::fixed<<std::setprecision(3)<<lx
            <<std::setw(9)<<std::fixed<<std::setprecision(3)<<ly
            <<std::setw(9)<<std::fixed<<std::setprecision(3)<<lz
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<ca
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<cb
            <<std::setw(7)<<std::fixed<<std::setprecision(2)<<cc
        <<" P 1       "
        <<"   1"
        <<std::endl;
    
    for (unsigned long i=0; i<frame.ats.size(); i++)
    {
        ostr<<"ATOM  "<<std::setw(5)<<i+1;
        ostr<<std::skipws<<std::setw(4)<<std::right<<frame.ats[i].name<<" ";
        ostr<<" "<<std::setw(3)<<std::right<<"X"<<" "<<" "<<"   1"<<" "<<"   ";
        ostr<<std::setw(8)<<std::fixed<<std::setprecision(3)<<frame.ats[i].x
            <<std::setw(8)<<std::fixed<<std::setprecision(3)<<frame.ats[i].y
            <<std::setw(8)<<std::fixed<<std::setprecision(3)<<frame.ats[i].z;
        ostr<<std::setw(6)<<std::fixed<<std::setprecision(2)
            <<(frame.ats[i].nprops.count("occupancy")>0?frame.ats[i].nprops["occupancy"]:0.0)
            <<std::setw(6)<<std::fixed<<std::setprecision(2)
            <<(frame.ats[i].nprops.count("beta")>0?frame.ats[i].nprops["beta"]:0.0);
        ostr<<"          "<<std::setw(2)<<frame.ats[i].name;
        ostr<<std::setw(2)<<(frame.ats[i].nprops.count("charge")>0?(int)frame.ats[i].nprops["charge"]:0);
        ostr<<std::endl;
    }
    ostr<<"END"<<std::endl;
}

 //ends namespace toolbox

/* I could not find a clean and portable way to use the RPC/XDR library from
 * C++. So I am directly reading bytes from the stream, even though this might 
 * be not very portable if one wants to move from system to system having 
 * binary-incompatible formats */

void xdr_dfloat(XDR* xdr, double& r)
{
    float rf;
    xdr_float(xdr, &rf);
    r=rf*10;   // XTC stores in nm, we typically want Angstrom.
}


int sizeofint(const int size) {
    unsigned int num = 1;
    int num_of_bits = 0;
   
    while (size >= num && num_of_bits < 32) {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}


int sizeofints( const int num_of_ints, unsigned int sizes[]) {
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i=0; i < num_of_ints; i++) {  
        tmp = 0;
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
            tmp = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0) {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;

}

static int receivebits(int buf[], int num_of_bits) {

    int cnt, num;
    unsigned int lastbits, lastbyte;
    unsigned char * cbuf;
    int mask = (1 << num_of_bits) -1;

    cbuf = ((unsigned char *)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int) buf[1];
    lastbyte = (unsigned int) buf[2];
   
    num = 0;
    while (num_of_bits >= 8) {
        lastbyte = ( lastbyte << 8 ) | cbuf[cnt++];
        num |=  (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -=8;
    }
    if (num_of_bits > 0) {
        if (lastbits < num_of_bits) {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | cbuf[cnt++];
        }
        lastbits -= num_of_bits;
        num |= (lastbyte >> lastbits) & ((1 << num_of_bits) -1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num;
}

static void receiveints(int buf[], const int num_of_ints, int num_of_bits,
        unsigned int sizes[], int nums[]) {
    int bytes[32];
    int i, j, num_of_bytes, p, num;
   
    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8) {
        bytes[num_of_bytes++] = receivebits(buf, 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
        bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
    }
    for (i = num_of_ints-1; i > 0; i--) {
        num = 0;
        for (j = num_of_bytes-1; j >=0; j--) {
            num = (num << 8) | bytes[j];
            p = num / sizes[i];
            bytes[j] = p;
            num = num - p * sizes[i];
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

static int magicints[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    8, 10, 12, 16, 20, 25, 32, 40, 50, 64,
    80, 101, 128, 161, 203, 256, 322, 406, 512, 645,
    812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501,
    8192, 10321, 13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536,
    82570, 104031, 131072, 165140, 208063, 262144, 330280, 416127, 524287, 660561,
    832255, 1048576, 1321122, 1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042,
    8388607, 10568983, 13316085, 16777216 };

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x):(y))
#endif
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x):(y))
#endif


bool ReadXTCFrame(FILE*  ifile, AtomFrame& curfr){
    
    AtomData curat; 
    curfr.index=0; curfr.ats.resize(0);  // flushes current frame data
    int nat, istep, dummy;  float ftime;
    XDR xdrif;
    xdrstdio_create(&xdrif, ifile, XDR_DECODE);
    xdr_int(&xdrif, &dummy);
    if (dummy != 1995) ERROR("Could not read magic number in XTC header");
    xdr_int(&xdrif, &nat);
    xdr_int(&xdrif, &istep);
    xdr_float(&xdrif, &ftime);
    
    xdr_dfloat(&xdrif, curfr.nprops["axx"]);
    xdr_dfloat(&xdrif, curfr.nprops["axy"]);
    xdr_dfloat(&xdrif, curfr.nprops["axz"]);
    xdr_dfloat(&xdrif, curfr.nprops["ayx"]);
    xdr_dfloat(&xdrif, curfr.nprops["ayy"]);
    xdr_dfloat(&xdrif, curfr.nprops["ayz"]);
    xdr_dfloat(&xdrif, curfr.nprops["azx"]);
    xdr_dfloat(&xdrif, curfr.nprops["azy"]);
    xdr_dfloat(&xdrif, curfr.nprops["azz"]);


    /* adapted from GROMACS library */
    float precision;
    int *size;
    int minint[3], maxint[3], mindiff;
    int lint1, lint2, lint3, oldlint1, oldlint2, oldlint3, smallidx;
    int minidx, maxidx, *lip, flag, k;
    int smallnum, smaller, larger, i, is_small, is_smaller, run, prevrun;
    int tmp, *thiscoord,  prevcoord[3];

    unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3, *luip;
    int bufsize, xdrid, lsize;
    unsigned int bitsize;
    float inv_precision, *lfp, lf;

    size=&nat;
    /* xdrs is open for reading */
   
    if (xdr_int(&xdrif, &lsize) == 0)
        return 0;
    if (*size != 0 && lsize != *size) {
        ERROR("Wrong number of coordinates in XTC read ")
    }
    *size = lsize;
    size3 = *size * 3;
    float *fp=new float[size3];
    int *ip=new int[size3];
    int *buf=new int[(int)(size3*1.2)];
    if (*size <= 9) {
        return (xdr_vector(&xdrif, (char *) fp, (unsigned int)size3,
                (unsigned int)sizeof(*fp), (xdrproc_t)xdr_float));
    }
    xdr_float(&xdrif, &precision);
    
    buf[0] = buf[1] = buf[2] = 0;
   
    xdr_int(&xdrif, &(minint[0]));
    xdr_int(&xdrif, &(minint[1]));
    xdr_int(&xdrif, &(minint[2]));

    xdr_int(&xdrif, &(maxint[0]));
    xdr_int(&xdrif, &(maxint[1]));
    xdr_int(&xdrif, &(maxint[2]));
           
    sizeint[0] = maxint[0] - minint[0]+1;
    sizeint[1] = maxint[1] - minint[1]+1;
    sizeint[2] = maxint[2] - minint[2]+1;
   
    /* check if one of the sizes is to big to be multiplied */
    if ((sizeint[0] | sizeint[1] | sizeint[2] ) > 0xffffff) {
        bitsizeint[0] = sizeofint(sizeint[0]);
        bitsizeint[1] = sizeofint(sizeint[1]);
        bitsizeint[2] = sizeofint(sizeint[2]);
        bitsize = 0; /* flag the use of large sizes */
    } else {
        bitsize = sizeofints(3, sizeint);
    }
   
    if (xdr_int(&xdrif, &smallidx) == 0)      
        return false;
    maxidx = MIN(LASTIDX, smallidx + 8) ;
    minidx = maxidx - 8; /* often this equal smallidx */
    smaller = magicints[MAX(FIRSTIDX, smallidx-1)] / 2;
    smallnum = magicints[smallidx] / 2;
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
    larger = magicints[maxidx];

    /* buf[0] holds the length in bytes */

    if (xdr_int(&xdrif, &(buf[0])) == 0)
        return false;
    if (xdr_opaque(&xdrif, (char *)&(buf[3]), (unsigned int)buf[0]) == 0)
        return false;
    buf[0] = buf[1] = buf[2] = 0;
   
    lfp = fp;
    inv_precision = 1.0 / precision;
    run = 0;
    i = 0;
    lip = ip;
    while ( i < lsize ) {
        thiscoord = (int *)(lip) + i * 3;

        if (bitsize == 0) {
            thiscoord[0] = receivebits(buf, bitsizeint[0]);
            thiscoord[1] = receivebits(buf, bitsizeint[1]);
            thiscoord[2] = receivebits(buf, bitsizeint[2]);
        } else {
            receiveints(buf, 3, bitsize, sizeint, thiscoord);
        }
       
        i++;
        thiscoord[0] += minint[0];
        thiscoord[1] += minint[1];
        thiscoord[2] += minint[2];
       
        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];
       
       
        flag = receivebits(buf, 1);
        is_smaller = 0;
        if (flag == 1) {
            run = receivebits(buf, 5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller--;
        }
        if (run > 0) {
            thiscoord += 3;
            for (k = 0; k < run; k+=3) {
                receiveints(buf, 3, smallidx, sizesmall, thiscoord);
                i++;
                thiscoord[0] += prevcoord[0] - smallnum;
                thiscoord[1] += prevcoord[1] - smallnum;
                thiscoord[2] += prevcoord[2] - smallnum;
                if (k == 0) {
                    /* interchange first with second atom for better
                     * compression of water molecules
                     */
                    tmp = thiscoord[0]; thiscoord[0] = prevcoord[0];
                            prevcoord[0] = tmp;
                    tmp = thiscoord[1]; thiscoord[1] = prevcoord[1];
                            prevcoord[1] = tmp;
                    tmp = thiscoord[2]; thiscoord[2] = prevcoord[2];
                            prevcoord[2] = tmp;
                    *lfp++ = prevcoord[0] * inv_precision;
                    *lfp++ = prevcoord[1] * inv_precision;
                    *lfp++ = prevcoord[2] * inv_precision;
                } else {
                    prevcoord[0] = thiscoord[0];
                    prevcoord[1] = thiscoord[1];
                    prevcoord[2] = thiscoord[2];
                }
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
        } else {
            *lfp++ = thiscoord[0] * inv_precision;
            *lfp++ = thiscoord[1] * inv_precision;
            *lfp++ = thiscoord[2] * inv_precision;          
        }
        smallidx += is_smaller;
        if (is_smaller < 0) {
            smallnum = smaller;
            if (smallidx > FIRSTIDX) {
                smaller = magicints[smallidx - 1] /2;
            } else {
                smaller = 0;
            }
        } else if (is_smaller > 0) {
            smaller = smallnum;
            smallnum = magicints[smallidx] / 2;
        }
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx] ;
    }

    
    curat.props.resize(0);
    lfp=fp;
    for (i=0; i<nat; i++)
    {
        curat.x=*(lfp++) * 10; 
        curat.y=*(lfp++) * 10; 
        curat.z=*(lfp++) * 10; 
        curfr.ats.push_back(curat);
    }
    delete [] fp; delete[] ip; delete[] buf;
    
    //~ if(!getline(istr, line)) return false;  //drops comment line
    //~ if(!getline(istr, line)) return false;  //reads header
    //~ curfr.index=0; curfr.ats.resize(0);
    //~ ss.clear(); ss<<line;
    //~ ss>>ikey>>ipbc>>nat>>dummy;
    //~ if(ipbc>0)
    //~ {
        //~ istr>>curfr.nprops["axx"]>>curfr.nprops["axy"]>>curfr.nprops["axz"];
        //~ istr>>curfr.nprops["ayx"]>>curfr.nprops["ayy"]>>curfr.nprops["ayz"];
        //~ istr>>curfr.nprops["azx"]>>curfr.nprops["azy"]>>curfr.nprops["azz"];
        //~ getline(istr,dummy);
    //~ }
    //~ nat=0;
    //~ for (i=0; getline(istr, line); ++i)
    //~ {
        //~ nat++;
        //~ ss.clear(); ss.str(line);
        //~ curat.props.resize(0);
        //~ ss>>curat.name>>dummy;
        //~ getline(istr, line); ss.clear(); ss.str(line);
        //~ ss>>curat.x>>curat.y>>curat.z;
        //~ if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        //~ if (ikey>0)
        //~ {
            //~ getline(istr, line); ss.clear(); ss.str(line);
            //~ ss>>cprop; curat.props.push_back(cprop); 
            //~ ss>>cprop; curat.props.push_back(cprop); 
            //~ ss>>cprop; curat.props.push_back(cprop);
        //~ }
        //~ if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        //~ if (ikey>1)
        //~ {
            //~ getline(istr, line); ss.clear(); ss.str(line);
            //~ ss>>cprop; curat.props.push_back(cprop); 
            //~ ss>>cprop; curat.props.push_back(cprop); 
            //~ ss>>cprop; curat.props.push_back(cprop);
        //~ }
        //~ if (ss.bad()) ERROR("Read failed in frame "<<curfr.index<<", on atom "<<i<<".");
        //~ curfr.ats.push_back(curat);
    //~ }
    return true;
    
 
}    

bool WriteXTCFrame(std::istream& istr, AtomFrame& frames){
    
    
}
};

// extern bool read_first_xtc(XDR *xd,char *filename, int *natoms,int *step,real *time, matrix box,rvec **x,real *prec);
// extern bool read_next_xtc(XDR *xd, int *natoms,int *step,real *time, matrix box,rvec *x,real *prec);
// extern bool write_xtc(XDR *xd, int natoms,int step, real time, matrix box, rvec *x, real prec);
// extern bool write_xtc(XDR *xd, int natoms,int step, real time, matrix box, rvec *x, real prec);
