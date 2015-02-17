#include "yocto/gfx/image.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/exception.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/gfx/ops/blob.hpp"
#include "yocto/gfx/ops/hist.hpp"
#include "yocto/gfx/ops/contrast.hpp"
#include "yocto/gfx/ops/blob.hpp"
#include "yocto/ptr/arc.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/wavelet/cwt.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;

static inline void put_rgb(void *addr, const rgba_t &c, const void *) throw()
{
    new (addr) rgb_t(c.r,c.g,c.b);
}

static inline rgba_t float2rgba(const void *addr,const void *) throw()
{
    const float   f = *(const float *)addr;
    const uint8_t u = conv::to_byte(1.0f-f);
    //const uint8_t u = conv::to_byte(f);
    return rgba_t(u,u,u,0xff);
}
/*
static inline rgba_t intensity2rgba(const void *addr,const void *) throw()
{
    const float   f = *(const float *)addr;
    const uint8_t u = conv::to_byte(f);
    //const uint8_t u = conv::to_byte(f);
    return rgba_t(u,u,u,0xff);
}
*/

namespace
{
    class slice
    {
    public:
        const unit_t lo;
        const unit_t hi;
        const unit_t count;
        inline slice( const unit_t __lo, const unit_t __hi ) throw() :
        lo(__lo),
        hi(__hi),
        count(hi-lo+1)
        {
        }

        inline ~slice() throw() {}

        inline slice(const slice &s ) throw() : lo(s.lo), hi(s.hi), count(s.count)
        {
        }


    private:
        YOCTO_DISABLE_ASSIGN(slice);
    };



    class slices : public vector<slice>, public counted
    {
    public:
        unit_t lo;
        unit_t hi;
        const size_t width;
        const size_t height;
        vector<float> x;
        vector<float> y;
        explicit slices(size_t w,size_t h) : vector<slice>(w,as_capacity), lo(0),hi(0),
        width(w),
        height(h),
        x(w,0),
        y(w,0)
        {
        }

        virtual ~slices() throw()
        {
        }

        typedef arc_ptr<slices> ptr;



    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(slices);
    };

    static inline float Sech( float x )
    {
        return 1.0f/Square(coshf(x));
    }

    class ParaWork
    {
    public:
        const slices &S;
        lockable     &access;
        const size_t  job_id;
        const string &outdir;

        explicit ParaWork(const slices &user_slices,
                          lockable     &user_access,
                          const size_t  user_job_id,
                          const string &user_outdir) :
        S(user_slices),
        access(user_access),
        job_id(user_job_id),
        outdir(user_outdir)
        {
        }

        ParaWork( const ParaWork &other ) throw() :
        S(other.S),
        access(other.access),
        job_id(other.job_id),
        outdir(other.outdir)
        {
        }

        inline void operator()(void)
        {
            const size_t n = S.size();
            vector<float>    shifts(n,0);
            vector<float>    scales(n,0);
            matrix<float>    W;
            numeric<float>::function Psi( cfunctor(Sech) );
            {
                scoped_lock guard(access);
                std::cerr << ".[" << job_id << "]";
                std::cerr.flush();
            }
            wavelet<float>::cwt(S.x, S.y, Psi, shifts, scales, W);

            image & IMG = image::instance();
            pixmap<float> Wimg(n,n);

            float wlo = W[1][1];
            wlo *= wlo;
            float whi = wlo;
            for(size_t i=1;i<=n;++i)
            {
                for(size_t j=1;j<=n;++j)
                {
                    float tmp = W[i][j]; tmp *= tmp;
                    if(tmp>whi) whi = tmp;
                    if(tmp<wlo) wlo = tmp;
                }
            }

            const float factor = 1.0f/(whi-wlo);
            for(size_t j=1;j<=n;++j)
            {
                for(size_t i=1;i<=n;++i)
                {
                    float tmp = W[i][j]; tmp *= tmp;
                    Wimg[j-1][i-1] = clamp<float>(0,(tmp-wlo)*factor,1);
                }
            }


            {
                const string outname = outdir + vformat("w%08u.png",unsigned(job_id));
                IMG["PNG"].save(outname, Wimg,image::get_rampf,NULL,NULL);
            }


        }

        ~ParaWork() throw()
        {

        }

    private:
        YOCTO_DISABLE_ASSIGN(ParaWork);
    };

}

#include "yocto/string/conv.hpp"
#include "yocto/threading/server.hpp"


int main(int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=3)
            throw exception("usage: %s directory extensions output_directory [xmin xmax]",prog);

        const string savext = "png";
        string       inpdir = argv[1];
        string       imgext = argv[2];
        string       outdir = argv[3];
        unit_t       xmin   = 0;
        if(argc>4)
        {
            xmin = strconv::to_int(argv[4],"xmin");
        }
        unit_t       xmax = -1;
        if(argc>5)
        {
            xmax = strconv::to_int(argv[5],"xmax");
        }
        vfs::as_directory(inpdir);
        vfs::as_directory(outdir);
        std::cerr << "Input  Dir=" << inpdir << std::endl;
        std::cerr << "Output Dir=" << outdir << std::endl;
        imgext.to_lower();

        image &IMG = image::instance();
        IMG.declare( new png_format()  );
        IMG.declare( new jpeg_format() );


        vfs & fs = local_fs::instance();

        std::cerr << "xmin=" << xmin << std::endl;
        std::cerr << "xmax=" << xmax << std::endl;

        xmin = max_of<unit_t>(xmin,0);
        if(xmax>=0)
        {
            if(xmax<xmin)
            {
                throw exception("xmax<xmin");
            }
        }

        // create output dir
        fs.create_dir(outdir,true);
        fs.remove_files_with_extension_in(outdir, savext);

        // scan input directory
        auto_ptr<vfs::scanner> scan( fs.new_scanner( inpdir ) );
        vector<slices::ptr>    work(1024,as_capacity);

        for( const vfs::entry *ep = scan->next(); ep; ep = scan->next() )
        {
            if(ep->is_directory())
                continue;
            const char *ep_ext = vfs::get_extension(ep->base_name);
            if(!ep_ext)
                continue;
            string ext_str = ep_ext;
            ext_str.to_lower();
            if(ext_str!=imgext)
                continue;

            const string &path = ep->path;
            std::cerr << "loading " << path << std::endl;

            // load image and convert to greyscale
            const bitmap::pointer bmp(IMG.load(path,3,put_rgb,NULL,NULL));
            pixmap3               img(bmp,NULL);
            pixmapf               pgs(img,rgb2gsf<rgba_t>);
            pixmapf               mask(pgs.w,pgs.h);

            // enhance contrast
            maximum_contrast(pgs);

            // automatic thresholding...
            histogram H;
            H.compute_from(pgs);
            const size_t t = H.threshold();
            threshold::apply(mask,t,pgs, threshold::keep_black);

            //clustering
            clusters cls;

            //blob
            blob B(mask,cls,true);
            cls.sort();
            std::cerr << "\t#cluster=" << cls.size() << std::endl;



            string outname = outdir + ep->base_name;
            vfs::change_extension(outname,savext);
            std::cerr << "\tsaving to " << outname << std::endl;


            size_t level = 0;
            if(cls.size()>0)
            {
                level = cls.front()->uuid;
            }
            //IMG["PNG"].save(outname, B,get_rgba_from_blob,&level,NULL);

            const size_t w = mask.w;
            if(work.size()>0 && work.back()->size() != w )
            {
                throw exception("width mismatch for '%s'", ep->base_name);
            }

            // update mask
            const unit_t h = mask.h;
            for(size_t j=0;j<h;++j)
            {
                pixmapf::row    &mj = mask[j];
                const blob::row &bj = B[j];
                for(size_t i=0;i<w;++i)
                {
                    if(bj[i]!=level)
                    {
                        mj[i] = 0.0f;
                    }
                }
            }

            IMG["PNG"].save(outname, mask,float2rgba,NULL,NULL);


            // convert image to a set of slices
            slices::ptr pS( new slices(w,h) );

            // that we add to the current set
            work.push_back(pS);

            // scan every slice, using the computed mask
            for(unit_t x=0;x<w;++x)
            {
                unit_t lo=0;
                unit_t hi=h;
                for(;lo<h;++lo)
                {
                    if(mask[lo][x]>0)
                        break;
                }


                for(--hi;hi>=0;--hi)
                {
                    if(mask[hi][x]>0)
                        break;
                }


                const slice s(lo,hi);
                pS->push_back(s);
                pS->x[x+1] = x+1;
                pS->y[x+1] = s.count;
            }

            //if(work.size()>=100) break;


        }


        const size_t n = work.size();
        std::cerr << "Analyzing " << n << " slices" << std::endl;

        if(n>0)
        {
            unit_t width  = work[1]->size();
            if(xmax<0||xmax>=width-1) xmax = width-1;
            if(xmin>xmax)             xmin = xmax;

            assert(xmin>=0);
            assert(xmax<width);
            assert(xmin<=xmax);

            std::cerr << "\txmin=" << xmin <<", xmax=" << xmax << std::endl;

            std::cerr << "Computing Shapes..." << std::endl;
            const size_t num_pop_back  = width - (xmax+1);
            const size_t num_pop_front = xmin;
            const size_t length        = xmax+1-xmin;
            for(size_t I=1;I<=n;++I)
            {
                std::cerr.flush();
                slices::ptr &pS = work[I];
                for(size_t k=num_pop_back; k>0;--k) { pS->pop_back();  pS->x.pop_back();  pS->y.pop_back();  }
                for(size_t k=num_pop_front;k>0;--k) { pS->pop_front(); pS->y.pop_front(); pS->y.pop_front(); }
                assert(pS->size() == length );
                const size_t h = work[1]->height;

                pixmapf shape(length,h);
                const unit_t hm = h/2;
                for(size_t i=0;i<length;++i)
                {
                    const slice &s = (*pS)[i+1];
                    const unit_t delta = s.count;
                    const unit_t z_lo  = hm - delta/2;
                    const unit_t z_hi  = hm + delta/2-1;
                    if(z_lo>=0)
                        shape[z_lo][i] = 1.0f;
                    if(z_hi<h)
                        shape[z_hi][i] = 1.0f;

                }

                {
                    const string outname = outdir + vformat("shape%08u.png",unsigned(I));
                    IMG["PNG"].save(outname, shape,float2rgba,NULL,NULL);
                }

            }

            std::cerr << "Computing Wavelets..." << std::endl;
            threading::server srv;
            for(size_t I=1;I<=n;++I)
            {
                const ParaWork Job( * work[I], srv.access, I, outdir);
                srv.enqueue(Job);
            }
            {
                srv.flush();
            }
            std::cerr << std::endl << std::endl;

#if 0
            std::cerr << "Wavelet transforms...." << std::endl;
            vector<float> shift(length,0);
            vector<float> scale(length,0);
            matrix<float> W;
            numeric<float>::function Psi( cfunctor(Sech) );
            pixmap<float> Wimg(length,length);
            for(size_t I=1;I<=n;++I)
            {
                std::cerr << ".";
                std::cerr.flush();
                slices::ptr &pS = work[I];
                assert(pS->size() == length );


                //Wavelet matrix
                wavelet<float>::cwt(pS->x, pS->y, Psi, shift, scale, W);
                float wlo = W[1][1];
                wlo *= wlo;
                float whi = wlo;
                for(size_t i=1;i<=length;++i)
                {
                    for(size_t j=1;j<=length;++j)
                    {
                        float tmp = W[i][j]; tmp *= tmp;
                        if(tmp>whi) whi = tmp;
                        if(tmp<wlo) wlo = tmp;
                    }
                }

                const float factor = 1.0f/(whi-wlo);
                for(size_t j=1;j<=length;++j)
                {
                    for(size_t i=1;i<=length;++i)
                    {
                        float tmp = W[i][j]; tmp *= tmp;
                        Wimg[j-1][i-1] = clamp<float>(0,(tmp-wlo)*factor,1);
                    }
                }

                {
                    const string outname = outdir + vformat("w%08u.png",unsigned(I));
                    IMG["PNG"].save(outname, Wimg,image::get_gsf,NULL,NULL);
                }

            }
            std::cerr << std::endl;
#endif

        }


        return 0;
    }
    catch(const exception &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << e.when() << std::endl;
    }
    catch(...)
    {
        std::cerr << "unhandled exception in " << prog << std::endl;
    }
    return 1;
}
