#include "yocto/graphics/ops/blob.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/tiff.hpp"
#include "yocto/graphics/ops/hist.hpp"

#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"

using namespace yocto;
using namespace graphics;

YOCTO_PROGRAM_START()
{

    image            &IMG = image::instance();
    threading::engine server(true);

    IMG.declare( new png_format()  );
    IMG.declare( new jpeg_format() );
    IMG.declare( new tiff_format() );
    
    const image::format &PNG = IMG["PNG"];
    histogram          H;
    histogram::patches hp;

    for(int arg=1;arg<argc;++arg)
    {
        const string  filename = argv[arg];
        const pixmap3 source( IMG.load3(filename,NULL) );
        const unit_t  w = source.w;
        const unit_t  h = source.h;

        PNG.save("source.png", source, NULL);

        //! build the histogram
        H.reset();
        histogram::create(hp,source,&server);
        histogram::launch(hp,source,&server);
        H.finish(hp,&server);


        //! find the threshold
        const size_t thr = H.threshold();

        //! keep the background
        pixmap3 target(w,h);
        threshold::apply(target,thr,source,threshold::keep_background);
        PNG.save("bg.png",target,NULL);

        // find the blobs
        get_named_color<blob::type> bproc;
        vector<size_t>              blobs;
        blob B(w,h,&server);
        B.build(blobs,target, 8, &server,0);
        std::cerr << "blobs=" << blobs << std::endl;
        PNG.save("blobs.png",B,bproc,NULL);


        if(false&&blobs.size()!=1)
        {
            throw exception("not one blob...");
        }

        pixmap3 surf(target);
        //B.transfer(1, surf, source);

        //______________________________________________________________________
        //
        // scan..
        //______________________________________________________________________
        vector<unit_t> thick(w,as_capacity);

        for(unit_t i=0;i<w;++i)
        {
            unit_t jmin = 0;
            for(unit_t j=0;j<h;++j)
            {
                if(B[j][i])
                {
                    jmin = j;
                    break;
                }
            }

            surf[jmin][i] = named_color::get("red");

            unit_t jmax = h-1;
            for(unit_t j=h-1;j>=0;--j)
            {
                if(B[j][i])
                {
                    jmax = j;
                    break;
                }
            }
            surf[jmax][i] = named_color::get("blue");
            thick.push_back(jmax-jmin);
        }
        PNG.save("surf.png",surf,NULL);
        {
            ios::wcstream fp("thick.dat");
            for(unit_t i=1;i<=w;++i)
            {
                fp("%ld %ld\n", long(i-1), long(thick[i]));
            }
        }
    }
    
}
YOCTO_PROGRAM_END()
