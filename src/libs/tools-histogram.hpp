/* A library to compute histograms in one and many dimensions
   --------------------------------------------------
   Author: Michele Ceriotti, 2008
   Distributed under the GNU General Public License
*/

#ifndef __TOOLS_HISTOGRAM_H
#define __TOOLS_HISTOGRAM_H

#include "tbdefs.hpp"
#include "math.h"
#include <vector>
#include <limits>

namespace toolbox {
    template <class U> class Histogram;
    template <class U> class HGOptions;

    enum HGWindowMode {HGWDelta,  HGWBox, HGWTriangle,  HGWGauss1,  HGWGauss2,  HGWGauss3,   HGWGauss5 };
    enum HGBoundaryMode {HGBNormal,  HGBHard,  HGBHardNorm,  HGBPeriodic };


    template <class U> class HGOptions<Histogram<U> > {
    public:
        HGWindowMode window; HGBoundaryMode walls;

        U window_width;

        U adaptive_eps;

        double (*window_function)(const U& a, const U& b);
        std::valarray<U> boundaries;


        HGOptions(const HGWindowMode& nwindow=HGWDelta, const U nwindow_width=(U)0., const U nadaptive_eps=(U)0.,
                    const std::valarray<U>& nbnd=std::valarray<U>(0), const HGBoundaryMode& nwalls=HGBNormal) :
            window(nwindow), window_width(nwindow_width), adaptive_eps(nadaptive_eps), walls(nwalls)
            { boundaries.resize(nbnd.size()); boundaries=nbnd; }
        HGOptions(const HGOptions& hgo) : walls(hgo.walls), window(hgo.window), window_width(hgo.window_width),
                  adaptive_eps(hgo.adaptive_eps),
                  window_function(hgo.window_function), boundaries(hgo.boundaries) {}
        HGOptions& operator=(const HGOptions& hgo)
        {
            if (&hgo==this) return *this;
            window=hgo.window;
            walls=hgo.walls;
            window_width=hgo.window_width;
            window_function=hgo.window_function;
            adaptive_eps=hgo.adaptive_eps;
            boundaries.resize(hgo.boundaries.size()); boundaries=hgo.boundaries;
            return *this;
        }
    };

    template <class U> class Histogram {
    private:
        std::valarray<double> bins;
        double below, above;
        double ndata;
        HGOptions<Histogram<U> > opts;
        double range, center;

    public:
        void reset()
        {
            bins.resize(opts.boundaries.size()-1);
            bins=0.; ndata=0; above=below=0.;
            range=opts.boundaries[opts.boundaries.size()-1]-opts.boundaries[0];
            center=(opts.boundaries[opts.boundaries.size()-1]+opts.boundaries[0])*0.5;
        }

        double samples() { return ndata; }

        void get_options(HGOptions<Histogram<U> >& rop)
        { rop=opts; }

        void set_options(const HGOptions<Histogram<U> >& rop)
        { opts=rop; reset(); }

        Histogram(const HGOptions<Histogram<U> >& rop=HGOptions<Histogram<U> >())
        { set_options(rop); }

        Histogram(const U& a, const U&b, unsigned long n, HGWindowMode hgw=HGWDelta, double ww=0.)
        {
            opts.window=hgw; opts.window_width=ww;
            opts.boundaries.resize(n+1);
            for (unsigned long i=0; i<=n; ++i) opts.boundaries[i]=a+(b-a)*(1.*i)/n;
            reset();
        }

        Histogram(const Histogram<U>& ac) { set_options(ac.opts); }

        Histogram<U>& operator=(const Histogram<U>& ac)
        {
            if (&ac==this) return *this;
            set_options(ac.opts);
        }

    //copy & assign...
        template<class T> Histogram(const Histogram<T>& ac);
        template<class T> Histogram<U>& operator=(const Histogram<T>& ac);

        void add(const U& nel, double weight=1.0);

    //syntactic sugar to insert a new element into the series, to move to next or previous series
        inline void operator << (const U& nel) { add(nel); }
        inline void operator << (const std::valarray<U>& nseries) { for (unsigned long i=0; i<nseries.size(); ++i) add(nseries[i]); }

        void get_bins(std::valarray<double>& rbins) const
        {
            rbins.resize(bins.size());
            rbins=bins*(1./ndata);
        }

        void get_bins(std::valarray<U>& cent, std::valarray<U>& ws, std::valarray<double>& rbins) const
        {
            rbins.resize(bins.size()); cent.resize(bins.size()); ws.resize(bins.size());
            rbins=bins*(1./ndata);
            for (unsigned long i=0; i<bins.size();++i)
            {
                cent[i]=(opts.boundaries[i]+opts.boundaries[i+1])/2.;
                ws[i]=(opts.boundaries[i+1]-opts.boundaries[i]);
            }
        }

        void get_outliers(double& rabove, double& rbelow) const
        {
            rabove=above/ndata; rbelow=below/ndata;
        }

        double get_totweight() const
        {
            return ndata;
        }
    };

    template<class U>
    std::ostream& operator<<(std::ostream& os, const Histogram<U>& his)
    {
        std::valarray<U> wx, ww, wf;
        his.get_bins(wx,ww,wf);
        os.precision(12);
        os.setf(std::ios::scientific);
        os.width(14);
        for (unsigned int i=0; i<wx.size(); ++i)
            os<<wx[i]<<" "<<wf[i]<<" "<<ww[i]<<"\n";
        return os;
    }

    template<class U>
    double __hgwfbox(const U& a, const U& b)
    {
        if (b<=-0.5 || a >=0.5) return 0.;

        U ia, ib;
        ia=(a>-0.5?a:-0.5);
        ib=(b< 0.5?b:0.5);
        return ib-ia;
    }

    template<class U>
    double __hgwftri(const U& a, const U& b)
    {
        if (b<=-1. || a >=1.) return 0.;

        U ia, ib;
        ia=(a>-1.?a:-1.);
        ib=(b< 1.?b:1.);

        return (ib*(2.-fabs(ib))-ia*(2.-fabs(ia)))*0.5;
    }

#define G1ERF1 5.03149608121118568
#define G1DY   0.24197072451914334980
    template<class U>
    double __hgwfgauss1(const U& a, const U& b)
    {
        if (b<=-1. || a >=1.) return 0.;

        U ia, ib;
        ia=(a>-1.?a:-1.);
        ib=(b< 1.?b:1.);

        return ((erf(ib/constant::sqrt2)-erf(ia/constant::sqrt2)) *0.5 -(ib-ia)*G1DY )*G1ERF1;
    }

#define G2ERF1 1.354030373543121523
#define G2DY   0.10798193302637610390
    template<class U>
    double __hgwfgauss2(const U& a, const U& b)
    {
        if (b<=-1. || a >=1.) return 0.;

        U ia, ib;
        ia=(a>-1.?a:-1.);
        ib=(b< 1.?b:1.);

        return ((erf(ib*2/constant::sqrt2)-erf(ia*2/constant::sqrt2)) *0.5 -(ib-ia)*G2DY )*G2ERF1;
    }

#define G3ERF1 1.030174731161562310
#define G3DY   0.01329554523581402153
    template<class U>
    double __hgwfgauss3(const U& a, const U& b)
    {
        if (b<=-1. || a >=1.) return 0.;

        U ia, ib;
        ia=(a>-1.?a:-1.);
        ib=(b< 1.?b:1.);

        return ((erf(ib*3/constant::sqrt2)-erf(ia*3/constant::sqrt2)) *0.5 -(ib-ia)*G3DY )*G3ERF1;
    }

#define G5ERF1 1.00001544073670377005281
#define G5DY   7.43359757367148854e-6
    template<class U>
    double __hgwfgauss5(const U& a, const U& b)
    {
        if (b<=-1. || a >=1.) return 0.;

        U ia, ib;
        ia=(a>-1.?a:-1.);
        ib=(b< 1.?b:1.);

        return ((erf(ib*5/constant::sqrt2)-erf(ia*5/constant::sqrt2)) *0.5 -(ib-ia)*G5DY )*G5ERF1;
    }

    template<class U>
    void Histogram<U>::add(const U& pnel, double weight)
    {
        long bs=bins.size(),ia=0, ib=bs/2, ic=bs;
        U nel=pnel;

        if (opts.walls==HGBPeriodic)  //folds into the interval
        {   nel=nel/range;   nel=center+range*(nel-round(nel));   }

        //First, finds in which bin lies the center (this works also if bins have different sizes)
        while (ib>ia && ib<ic)
        {
            if (nel>opts.boundaries[ib+1])
            {   ia=ib;  ib=(ib+ic+1)/2; }
            else
            {
                if (nel>opts.boundaries[ib]) break;
                else { ic=ib; ib=(ia+ib)/2;  }
            }
        }
        //now, nel is between boundaries[ib] and boundaries[ib+1] OR below;
        double (*wf)(const U& a, const U& b);
        switch(opts.window)
        {
        case HGWDelta:
            if (ib==0)
            {
                if (nel>=opts.boundaries[0]) bins[0]+=weight/(opts.boundaries[1]-opts.boundaries[0]);
                else below+=weight;
            }
            else if(ib==bs) above+=weight;
            else bins[ib]+=weight/(opts.boundaries[ib+1]-opts.boundaries[ib]);
            break;
        case HGWBox:
            wf=__hgwfbox;
            break;
        case HGWTriangle:
            wf=__hgwftri;
            break;
        case HGWGauss1:
            wf=__hgwfgauss1;
            break;
        case HGWGauss2:
            wf=__hgwfgauss2;
            break;
        case HGWGauss3:
            wf=__hgwfgauss3;
            break;
        case HGWGauss5:
            wf=__hgwfgauss5;
            break;
        default:
            ERROR("Windowing mode not implemented yet!\n");
        }
        if (opts.window!=HGWDelta)
        {
            double nb=0.; double hw, ww;

            // use an adaptive width based on the current number of samples at the central bin.
            // note that this assumes that the point weight, if given, corresponds to a number of
            // (uncorrelated) samples otherwise the relation between spread and error is lost
            if (opts.adaptive_eps > 0.0)
            {
                if (bins[ib]==0.0)
                {   // spreads the sample over ALL the bins!
                    bins += weight/range;
                    ndata += weight;
                    return;
                }
                else
                {
					/* the idea here is as follows. the central bin already contains n=N P(x)*dx samples.
					 * if we have in total added N samples, assuming a binomial distribution, 
					 * this has a variance N P(x) dx (1-P(x)dx) ~ N P dx (the probability of being in a given bin
					 * is small) which is equal to n. so the relative error in P(x) is of the order of 
					 * 1/sqrt(n). If we use a kernel of width ww, we effectively average over
					 * O(ww/bin size) cells, which for simplicity we assume to be uncorrelated.
					 * Then, the error goes down to sqrt(binsize/(n ww)). if we want the relative error to be 
					 * aeps, then aeps^2=binsize/(n ww) so ww=binsize/(n aeps^2). since n/binsize = bin value ...
					*/
                    ww = 1.0/(bins[ib]*opts.adaptive_eps*opts.adaptive_eps);
                    if (ww>range) ww = range;
                    if (ww<opts.window_width) ww=opts.window_width;
                }
            }
            else ww = opts.window_width;

            switch(opts.walls)
            {
            case HGBPeriodic: // Periodic walls, apply minimum image convention
                //adds weights for bins smaller than ib (possibly going round)
                for (ia=ib-1; ia>=0; --ia)
                {
                    nb=wf((opts.boundaries[ia]-nel)/ww,(opts.boundaries[ia+1]-nel)/ww)
                            /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    if (nb==0) break; else bins[ia]+=nb*weight;
                }
                if (ia<0)  //spill over
                {
                    hw=nel+range;
                    //makes sure not to count points twice by stopping at ib
                    for (ia=bs-1; ia>=ib; --ia)
                    {
                        nb=wf((opts.boundaries[ia]-hw)/ww,(opts.boundaries[ia+1]-hw)/ww)
                                /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                        if (nb==0) break; else bins[ia]+=nb*weight;
                    }
                    if (ia<ib) break; // means the window is crazily broad and we went around! arguably we should actually raise an error
                }

                for (ia=ib; ia<bs; ++ia)
                {
                    nb=wf((opts.boundaries[ia]-nel)/ww,(opts.boundaries[ia+1]-nel)/ww)
                            /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    if (nb==0) break; else bins[ia]+=weight*nb;
                }

                if (ia==bs)
                {
                    hw=nel-range;
                    for (ia=0; ia<ib; ++ia)
                    {
                        nb=wf((opts.boundaries[ia]-hw)/ww,(opts.boundaries[ia+1]-hw)/ww)
                                /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                        if (nb==0) break; else bins[ia]+=nb*weight;
                    }
                    if (ia==ib) break; // means the window is crazily broad and we went around! arguably we should actually raise an error
                }
                break;

            case HGBNormal: // Normal walls, density will spill out
                //adds weights for bins smaller than ib
                for (ia=ib-1; ia>=0; --ia)
                {
                    nb=wf((opts.boundaries[ia]-nel)/ww,(opts.boundaries[ia+1]-nel)/ww)
                            /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    if (nb==0) break; else bins[ia]+=nb*weight;
                }
                if (ia<0)
                    below+=weight*wf(-std::numeric_limits<U>::max(),(opts.boundaries[0]-nel)/ww);
                for (ia=ib; ia<bs; ++ia)
                {
                    nb=wf((opts.boundaries[ia]-nel)/ww,(opts.boundaries[ia+1]-nel)/ww)
                            /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    if (nb==0) break; else bins[ia]+=weight*nb;
                }
                if (ia==bs)
                    above+=weight*wf((opts.boundaries[bs]-nel)/ww,std::numeric_limits<U>::max());
                break;

            case HGBHard:  // Hard walls, constraint the density to the available space -- with asymmetric kernel.
                // if the point is outside it is completely discarded
                if (nel<opts.boundaries[0]) {below+=weight; break; }
                if (nel>opts.boundaries[bs]) {above+=weight; break; }
                // integrate below nel
                hw=nel-opts.boundaries[0];
                if (hw>ww) hw=ww;

                //get the amount of density that would spill out
                if (hw==0) bins[0]+=0.5*weight/(opts.boundaries[1]-opts.boundaries[0]);
                else
                {
                    for (ia=ib-1; ia>=0; --ia)
                    {
                        nb=wf((opts.boundaries[ia]-nel)/hw,(opts.boundaries[ia+1]-nel)/hw)
                                /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                        if (nb==0) break; else bins[ia]+=nb*weight;
                    }
                    bins[ib]+=wf((opts.boundaries[ib]-nel)/hw,0)/(opts.boundaries[ib+1]-opts.boundaries[ib])*weight;
                }

                // integrate above nel
                hw=opts.boundaries[bs]-nel;
                if (hw>ww) hw=ww;

                if (hw==0) bins[bs-1]+=0.5*weight/(opts.boundaries[bs]-opts.boundaries[bs-1]);
                else
                {
                    for (ia=ib+1; ia<bs; ++ia)
                    {
                        nb=wf((opts.boundaries[ia]-nel)/hw,(opts.boundaries[ia+1]-nel)/hw)
                                /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                        if (nb==0) break; else bins[ia]+=weight*nb;
                    }
                    bins[ib]+=wf(0,(opts.boundaries[ib+1]-nel)/hw)/(opts.boundaries[ib+1]-opts.boundaries[ib])*weight;
                }

                break;
            case HGBHardNorm:  // Hard walls, constraint the density to the available space -- with renormalization.
                // if the point is outside it is completely discarded
                if (nel<opts.boundaries[0]) {below+=weight; break; }
                if (nel>opts.boundaries[bs]) {above+=weight; break; }

                hw=0.0;
                //get the amount of density that would spill out below
                if (nel-opts.boundaries[0]<ww)
                {  hw+=wf(-1.0,(opts.boundaries[0]-nel)/ww);  }
                if (opts.boundaries[bs]-nel<ww)
                {  hw+=wf((opts.boundaries[bs]-nel)/ww,1.0);  }
                //get the re-normalizing factor
                hw=1.0/(1.0-hw);

                above=1;
                for (ia=ib-1; ia>=0; --ia)
                {
                    nb=wf((opts.boundaries[ia]-nel)/ww,(opts.boundaries[ia+1]-nel)/ww)
                            /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    above-=nb*hw*weight*(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    if (nb==0) break; else bins[ia]+=nb*hw*weight;
                }
                for (ia=ib; ia<bs; ++ia)
                {
                    nb=wf((opts.boundaries[ia]-nel)/ww,(opts.boundaries[ia+1]-nel)/ww)
                            /(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    above-=nb*hw*weight*(opts.boundaries[ia+1]-opts.boundaries[ia]);
                    if (nb==0) break; else bins[ia]+=weight*hw*nb;
                }
            }
        }

        ndata+=weight;
    }

/***********************************************************
 *             N-DIMENSIONAL HISTOGRAM                     *
 ***********************************************************/

template <class U> class NDHistogram {
private:
    unsigned long dim;
    std::valarray<double> bins, vols, range, center;
    std::valarray<long> nbin;
    double outliers;
    double ndata;
    std::valarray<HGOptions<Histogram<U> > > opts;
    long c2b(const std::valarray<long>& vl)
    {
        long ts=1,rp=0;
        for (int i=0; i<dim;++i)
        { rp+=vl[i]*ts; ts*=nbin[i]; }
        return rp;
    }

public:
    void reset()
    {
        long tsz=1; nbin.resize(dim);
        for (int i=0; i<dim; ++i) tsz*=(nbin[i]=(opts[i].boundaries.size()-1));
        if (tsz<0) tsz=0;

        bins.resize(tsz);
        bins=0.; ndata=0; outliers=0.;

        std::valarray<long> cp(dim); cp=0;
        long k=0; vols.resize(tsz);
        while (cp[dim-1]<nbin[dim-1])
        {
            vols[k]=1.;
            for (int i=0; i<dim; ++i) vols[k]*=(opts[i].boundaries[cp[i]+1]-opts[i].boundaries[cp[i]]);
            k++; cp[0]++;
            for (int i=0; i<dim-1; ++i) if (cp[i]>=nbin[i]) {cp[i]=0; ++cp[i+1];}
        }

        range.resize(dim); center.resize(dim);
        for (int i=0; i<dim; ++i)
        {
            range[i]=opts[i].boundaries[opts[i].boundaries.size()-1]-opts[i].boundaries[0];
            center[i]=0.5*(opts[i].boundaries[opts[i].boundaries.size()-1]+opts[i].boundaries[0]);
        }
    }

    double samples() { return ndata; }

    void get_options(std::valarray<HGOptions<Histogram<U> > >& rop)
    { rop.resize(opts.size); rop=opts; }

    void set_options(const std::valarray<HGOptions<Histogram<U> > >& rop)
    { opts.resize(dim=rop.size()); opts=rop; reset(); }

    NDHistogram(const std::valarray<HGOptions<Histogram<U> > >& rop)
    { set_options(rop); }

//copy & assign...
    template<class T> NDHistogram(const NDHistogram<T>& ac);
    template<class T> NDHistogram<U>& operator=(const NDHistogram<T>& ac);

    void add(const std::valarray<U>& nel, double weight=1.0);

//syntactic sugar to insert a new element into the series, to move to next or previous series
    inline void operator << (const std::valarray<U>& nel) { add(nel); }

    double max() const { return bins.max()/ndata; }
    double min() const { return bins.min()/ndata; }


    void get_bins(std::valarray<double>& rbins) const
    {
        rbins.resize(bins.size());
        rbins=bins*(1./ndata);
    }

    void get_bins(std::valarray<std::valarray<U> >& cent, std::valarray<std::valarray<U> >& ws, std::valarray<double>& rbins) const
    {
        rbins.resize(bins.size()); cent.resize(dim); ws.resize(dim);
        rbins=bins*(1./ndata);
        for (unsigned long i=0; i<dim;++i)
        {
            cent[i].resize(nbin[i]); ws[i].resize(nbin[i]);
            for (unsigned long k=0; k<nbin[i];++k)
            {
                cent[i][k]=(opts[i].boundaries[k]+opts[i].boundaries[k+1])/2.;
                ws[i][k]=(opts[i].boundaries[k+1]-opts[i].boundaries[k]);
            }
        }
    }

    void get_bin(const std::valarray<long>& index, std::valarray<double>& center, double& val)
    {
        long ibin=c2b(index);
        center.resize(dim); val=bins[ibin]*(1./ndata);
        for (unsigned long i=0; i<dim;++i)
            center[i]=(opts[i].boundaries[index[i]]+opts[i].boundaries[index[i]+1])/2.;
    }

    void get_outliers(double& routliers) const
    {
        routliers=outliers/ndata;
    }

    double get_totweight() const
    {
        return ndata;
    }
};

template<class U>
std::ostream& operator<<(std::ostream& os, const NDHistogram<U>& his)
{
    std::valarray<U> wf;
    std::valarray<std::valarray<U> > wx, ww;
    his.get_bins(wx,ww,wf);
    os.precision(12);
    os.setf(std::ios::scientific);
    os.width(14);

    U outliers; his.get_outliers(outliers);
    os<<"# Fraction of outliers: "<<outliers<<"\n";
    for (unsigned int i=0; i<wx.size(); ++i)
    {
        os<<"# x("<<i<<"): ";
        for (unsigned int j=0; j<wx[i].size(); ++j) os<<wx[i][j]<<" ";
        os<<"\n# w("<<i<<"): ";
        for (unsigned int j=0; j<ww[i].size(); ++j) os<<ww[i][j]<<" ";
        os<<"\n";
    }

    for (unsigned int i=0; i<wf.size(); ++i)
        os<<wf[i]<<"\n";
    return os;
}

template<class U>
void NDHistogram<U>::add(const std::valarray<U>& pnel, double weight)
{
    // adds a data point to the n-dimensional histogram
    std::valarray<U> nel(pnel);

    for (int i=0; i<dim; ++i) if (opts[i].walls==HGBPeriodic)  //folds into the interval
    {   nel[i]=nel[i]/range[i];   nel[i]=center[i]+range[i]*(nel[i]-round(nel[i]));   }

    std::valarray<long> p(dim); //these are the coordinates of the center
    long bs, ia, ib, ic;

    //first, finds the "coordinates" of the center along all directions
    for (int i=0; i<dim; ++i)
    {
        bs=nbin[i]; ia=0; ic=bs; ib=(ic+1)/2;
        while (ib>ia && ib<ic)
        {
            if (nel[i]>opts[i].boundaries[ib+1])
            {   ia=ib;  ib=(ib+ic+1)/2; }
            else
            {
                if (nel[i]>opts[i].boundaries[ib]) break;
                else { ic=ib; ib=(ia+ib)/2;  }
            }
        }
        p[i]=ib;
    }
    ndata+=weight;

    //now, nel is between boundaries[p[i]] and boundaries[p[i]+1] OR below in each dim
    double (*wf)(const U& a, const U& b);
    std::valarray<std::valarray<double> > tbins(dim);
    int i;

    // so, a couple of words to explain how this works. we consider a D-dimensional
    // kernel which is a product of 1-d kernels K(x1,x2...)=K(x1)K(x2)....
    // so the integral of the kernel over each bin is the product of the integrals over
    // the "marginal bins" in each of the D dimensions.
    // so, we first compute these marginal integrals, and then make the products to get
    // the bin weights in D dimensions.
    for (i=0; i<dim; ++i)
    {
        tbins[i].resize(nbin[i]); tbins[i]=0.;

        std::valarray<long> ap(p);

        switch(opts[i].window)
        {
        case HGWDelta:
            if ((p[i]==0 && nel[i]< opts[i].boundaries[0]) || p[i]>=nbin[i])  break;
            tbins[i][p[i]]=1.0;
            break;
        case HGWBox:
            wf=__hgwfbox;
            break;
        case HGWTriangle:
            wf=__hgwftri;
            break;
        case HGWGauss1:
            wf=__hgwfgauss1;
            break;
        case HGWGauss2:
            wf=__hgwfgauss2;
            break;
        case HGWGauss3:
            wf=__hgwfgauss3;
            break;
        case HGWGauss5:
            wf=__hgwfgauss5;
            break;
        default:
            ERROR("Windowing mode not implemented yet!\n");
        }

        if (opts[i].window!=HGWDelta)
        {
            double nb=0.;   double hw, ww;

            // use an adaptive width based on the current number of samples at the central bin.
            // note that this assumes that the point weight, if given, corresponds to a number of
            // (uncorrelated) samples
            if (opts[i].adaptive_eps > 0.0)
            {
                if (bins[c2b(p)]==0.0)
                {   // spreads the sample over ALL the bins!
                    tbins[i] += 1.0/nbin[i];
                    ww = -1.0;
                }
                else
                {
					/* the idea here is analogous to the 1D case, only we must take care of dimensionality. 
					 * the central bin already contains n=N P(x)*dx samples.
					 * if we have in total added N samples, assuming a binomial distribution, 
					 * this has a variance N P(x) dx (1-P(x)dx) ~ N P dx (the probability of being in a given bin
					 * is small) which is equal to n. so the relative error in P(x) is of the order of 
					 * 1/sqrt(n). If we use a kernel of width ww, we effectively average over
					 * O(m=ww^D/binvolume) cells, which for simplicity we assume to be uncorrelated.
					 * Then, the error goes down to sqrt(1/(n m)). if we want the relative error to be 
					 * aeps, then aeps^2=1/(m ww) so ww=(binvolume/(n aeps^2))^1/d. since n/binvolume = bin value ...
					*/
                    ww = 1.0/pow(bins[c2b(p)]*opts[i].adaptive_eps*opts[i].adaptive_eps,1.0/dim);
                    
                    if (ww>range[i]) ww = range[i];
                    if (ww<opts[i].window_width) ww=opts[i].window_width;                    
                }
            }
            else ww = opts[i].window_width;            

            if (ww>=0) switch(opts[i].walls)
            {
            case HGBPeriodic: // Periodic walls, apply minimum image convention
                //adds weights for bins smaller than ib (possibly going round)
                for (ia=p[i]-1; ia>=0; --ia)
                {
                    nb=wf((opts[i].boundaries[ia]-nel[i])/ww,(opts[i].boundaries[ia+1]-nel[i])/ww);
                    if (nb==0) break; else tbins[i][ia]+=nb;
                }
                if (ia<0)  //spill over
                {
                    hw=nel[i]+range[i];
                    //makes sure not to count points twice
                    for (ia=nbin[i]-1; ia>=p[i]; --ia)
                    {
                        nb=wf((opts[i].boundaries[ia]-hw)/ww,(opts[i].boundaries[ia+1]-hw)/ww);
                        if (nb==0) break; else tbins[i][ia]+=nb;
                    }
                    if (ia<p[i]) break; // means the window is crazily broad and we went around! arguably we should actually raise an error
                }

                for (ia=p[i]; ia<nbin[i]; ++ia)
                {
                    nb=wf((opts[i].boundaries[ia]-nel[i])/ww,(opts[i].boundaries[ia+1]-nel[i])/ww);
                    if (nb==0) break; else tbins[i][ia]+=nb;
                }
                if (ia==nbin[i])
                {
                    hw=nel[i]-range[i];
                    for (ia=0; ia<p[i]; ++ia)
                    {
                        nb=wf((opts[i].boundaries[ia]-hw)/ww,(opts[i].boundaries[ia+1]-hw)/ww);
                        if (nb==0) break; else tbins[i][ia]+=nb;
                    }
                    if (ia==p[i]) break; // means the window is crazily broad and we went around! arguably we should actually raise an error
                }
                break;
            case HGBNormal: // Normal walls, density will spill out
                for (ia=p[i]-1; ia>=0; --ia)
                {
                    nb=wf((opts[i].boundaries[ia]-nel[i])/ww,(opts[i].boundaries[ia+1]-nel[i])/ww);
                    if (nb==0) break; else tbins[i][ia]+=nb;
                }
                for (ia=p[i]; ia<nbin[i]; ++ia)
                {
                    nb=wf((opts[i].boundaries[ia]-nel[i])/ww,(opts[i].boundaries[ia+1]-nel[i])/ww);
                    if (nb==0) break; else tbins[i][ia]+=nb;
                }
                break;

            case HGBHard:  // Hard walls, constraint the density to the available space -- unless a point is completely outside.
                // if point is outside, discards it completely
                if (nel[i]<opts[i].boundaries[0]) break;
                if (nel[i]>opts[i].boundaries[nbin[i]]) break;

                hw=nel[i]-opts[i].boundaries[0];
                if (hw>ww) hw=ww;

                if (hw==0) tbins[i][0]+=0.5;
                else
                {
                    for (ia=p[i]-1; ia>=0; --ia)
                    {
                        nb=wf((opts[i].boundaries[ia]-nel[i])/hw,(opts[i].boundaries[ia+1]-nel[i])/hw);
                        if (nb==0) break; else tbins[i][ia]+=nb;
                    }
                    tbins[i][p[i]]+=wf((opts[i].boundaries[p[i]]-nel[i])/hw,0);
                }

                hw=opts[i].boundaries[nbin[i]]-nel[i];
                if (hw>ww) hw=ww;

                if (hw==0) tbins[i][nbin[i]-1]+=0.5;
                else
                {
                    for (ia=p[i]+1; ia<nbin[i]; ++ia)
                    {
                        nb=wf((opts[i].boundaries[ia]-nel[i])/hw,(opts[i].boundaries[ia+1]-nel[i])/hw);
                        if (nb==0) break; else tbins[i][ia]+=nb;
                    }
                    tbins[i][p[i]]+=wf(0,(opts[i].boundaries[p[i]+1]-nel[i])/hw);
                }

                break;
            }
        }
    }
    
    //now we must make all the products and increment
    if (i<dim) { outliers+=weight; return; }
    std::valarray<long> minp(dim), maxp(dim); minp=0; maxp=0;
    double outs=1.;
    //boundaries of nonzero region (only accumulate over nonzero regions)
    for (i=0; i<dim; ++i)
    {
        int j;
        for (j=0; j<nbin[i] && tbins[i][j]==0.;) ++j;
        if (j==nbin[i]) { outliers+=weight; return; }
        minp[i]=j;
        if (minp[i]==0 && opts[i].walls==HGBPeriodic) j=nbin[i];  //Take care of nonzero kernels wrapping over the periodic boundary
        else for (;j<nbin[i]&&tbins[i][j]!=0.;) ++j;
        maxp[i]=j;
    }
    
    // we use a vector of indices to perform ndim nested loops
    std::valarray<long> cp(minp);
    while (cp[dim-1]<maxp[dim-1])
    {
        int j;

        long k=c2b(cp); // get the index of the bin corresponding to the current set of indices
        double tv=1.; for (i=0; i<dim; ++i) tv*=tbins[i][cp[i]];
        bins[k]+=tv/vols[k]*weight;
        outs-=tv;

        // increment the bin index
        cp[0]++;
        for (i=0; i<dim-1; ++i)
        {
            // in periodic binning, we can at least try to skip zero ranges here
            if (opts[i].walls==HGBPeriodic) while(cp[i]<maxp[i] && tbins[i][cp[i]] == 0.0) cp[i]++;
            if (cp[i]>=maxp[i]) {cp[i]=minp[i]; ++cp[i+1];} // wrap around
        }
    }
    outliers+=outs*weight;
}

}//ends namespace toolbox
#endif  //ends ifdef __TOOLS_HISTOGRAM_H
