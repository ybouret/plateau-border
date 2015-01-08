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

int main(int argc, char *argv[] )
{
    const char *prog = vfs::get_base_name(argv[0]);
    try
    {
        if(argc<=3)
            throw exception("usage: %s directory extensions output_directory",prog);
        
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
        
        // scan input directory
        auto_ptr<vfs::scanner> scan( fs.new_scanner( inpdir ) );
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
            blob(mask,cls,true);
            
            const string outname = outdir + ep->base_name;
            std::cerr << "\tsaving to " << outname << std::endl;
            IMG["PNG"].save(outname, mask,float2rgba,NULL,NULL);
            
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
