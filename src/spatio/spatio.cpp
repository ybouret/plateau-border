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

using namespace yocto;
using namespace gfx;

static inline void put_rgb(void *addr, const rgba_t &c, const void *) throw()
{
    new (addr) rgb_t(c.r,c.g,c.b);
}

static inline rgba_t float2rgba(const void *addr,const void *) throw()
{
    const float   f = *(const float *)addr;
    const uint8_t u = conv::to_byte(f);
    return rgba_t(u,u,u,0xff);
}

static inline rgba_t get_rgba_from_blob(const void *addr, const void *args)
{
    const size_t  value = *(const size_t *)addr;
    if(value<=0)
    {
        return rgba_t(0,0,0);
    }
    else
    {
        const size_t  level = *(const size_t *)args;
        const uint8_t u     = (level == value) ? 0xff : 0x00;
        return rgba_t(u,u,u);
    }
}


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
        
        inline slice(const slice &s ) throw() : lo(s.lo), hi(s.hi), count(s.count) {}
        
        
    private:
        YOCTO_DISABLE_ASSIGN(slice);
    };
    
    
    
    class slices : public vector<slice>, public counted
    {
    public:
        unit_t lo;
        unit_t hi;
        explicit slices(size_t w) : vector<slice>(w,as_capacity), lo(0),hi(0)
        {
        }
        
        virtual ~slices() throw()
        {
        }
        
        typedef arc_ptr<slices> ptr;
        
        
        
        void bracket() throw()
        {
            lo = hi = 0;
            const size_t w = size();
            if(w>0)
            {
                const array<slice> &S = *this;
                lo = S[1].lo;
                hi = S[1].hi;
                for(size_t i=2;i<=w;++i)
                {
                    {
                        if( S[i].lo < lo ) lo = S[i].lo;
                            if( S[i].hi > hi ) hi = S[i].hi;
                                }
                }
            }
        }
        
    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(slices);
    };
    
}

int main(int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=3)
            throw exception("usage: %s directory extensions output_directory",prog);
        
        const string savext = "png";
        string       inpdir = argv[1];
        string       imgext = argv[2];
        string       outdir = argv[3];
        vfs::as_directory(inpdir);
        vfs::as_directory(outdir);
        imgext.to_lower();
        
        image &IMG = image::instance();
        IMG.declare( new png_format()  );
        IMG.declare( new jpeg_format() );
        
        
        vfs & fs = local_fs::instance();
        
        // create output dir
        fs.create_sub_dir(outdir);
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
            const bitmap::pointer bmp(IMG.load(path,3,put_rgb,NULL));
            pixmap3               img(bmp,NULL);
            pixmapf               pgs(img,rgb2gsf<rgba_t>);
            pixmapf               mask(pgs.w,pgs.h);
            
            // enhance contrast
            maximum_contrast(pgs);
            
            // threshold
            histogram H;
            H.compute_from(pgs);
            const size_t t = H.threshold();
            threshold::apply(mask,t,pgs, threshold::keep_black);
            
            //cluster
            clusters cls;
            
            //blob
            blob B(mask,cls,true);
            cls.sort();
            std::cerr << "\t#cluster=" << cls.size() << std::endl;
            
            
            
            string outname = outdir + ep->base_name;
            vfs::change_extension(outname,savext);
            std::cerr << "\tsaving to " << outname << std::endl;
            
            //IMG["PNG"].save(outname, mask,float2rgba,NULL,NULL);
            
            size_t level = 0;
            if(cls.size()>0)
            {
                level = cls.front()->uuid;
            }
            IMG["PNG"].save(outname, B,get_rgba_from_blob,&level,NULL);
            
            const size_t w = mask.w;
            if(work.size()>0 && work.back()->size() != w )
            {
                throw exception("width mismatch for '%s'", ep->base_name);
            }
            const unit_t h = mask.h;
            slices::ptr pS( new slices(w) );
            work.push_back(pS);
            for(unit_t x=0;x<w;++x)
            {
                unit_t lo=0;
                for(;lo<h;++lo)
                {
                    if(mask[lo][x]>0)
                        break;
                }
                
                unit_t hi=h;
                for(--hi;hi>=0;--hi)
                {
                    if(mask[hi][x]>0)
                        break;
                }
                
                const slice s(lo,hi);
                pS->push_back(s);
            }
            
            //if(work.size()>10) break;
        }
        
        
        const size_t n = work.size();
        std::cerr << "Analyzing " << n << " slices" << std::endl;
        if(n>0)
        {
            for(size_t i=1;i<=n;++i)
            {
                work[i]->bracket();
            }
            unit_t lo = work[1]->lo;
            unit_t hi = work[1]->hi;
            for(size_t i=2;i<=n;++i)
            {
                if(work[i]->lo<lo) lo = work[i]->lo;
                if(work[i]->hi>hi) hi = work[i]->hi;
            }
            std::cerr << "lower=" << lo << std::endl;
            std::cerr << "upper=" << hi << std::endl;
            
            if(lo>hi)
            {
                throw exception("no data was found...");
            }
            
            const unit_t scale  = hi-lo+1;
            const size_t width  = work[1]->size();
            pixmapf      spatio(width,n);
            for(size_t j=1;j<=n;++j)
            {
                pixmapf::row &sp = spatio[j-1];
                const slices &sl = *work[j];
                for(size_t i=1,im=0;i<=width;++i,++im)
                {
                    const slice &s = sl[i];
                    if(s.count>0)
                    {
                        sp[im] = float(s.count)/scale;
                    }
                }
            }
            
            {
                const string outname = "spation.png";
                IMG["PNG"].save(outname, spatio,float2rgba,NULL,NULL);
            }
            
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
