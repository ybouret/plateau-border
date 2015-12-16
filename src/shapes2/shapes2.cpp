#include "yocto/graphics/ops/blob.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/tiff.hpp"
#include "yocto/graphics/ops/hist.hpp"
#include "yocto/math/core/tao.hpp"

#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/glsf.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/trigconv.hpp"

using namespace yocto;
using namespace graphics;
using namespace math;

namespace
{
    class Soliton
    {
    public:
        explicit Soliton()
        {
        }

        virtual ~Soliton() throw()
        {

        }

        inline double Eval( double x, const array<double> &a )
        {
            assert(a.size()>=5);
            const double args = a[4]*(x-a[5]);
            const double csha = cosh(args);
            return (a[1] + x*a[2]) - a[3]/(csha*csha);
        }


    private:
        YOCTO_DISABLE_COPY_AND_ASSIGN(Soliton);
    };
}

YOCTO_PROGRAM_START()
{

    image            &IMG = image::instance();
    threading::engine server(true);

    IMG.declare( new png_format()  );
    IMG.declare( new jpeg_format() );
    IMG.declare( new tiff_format() );

    const image::format &PNG = IMG["PNG"];
    histogram            H;
    histogram::patches   hp;

    //ios::ocstream::overwrite("scaling.dat");

    const RGB blue    = named_color::get("blue");
    const RGB red     = named_color::get("red");

    for(int arg=1;arg<argc;++arg)
    {
        const string  filename = argv[arg];

        //______________________________________________________________________
        //
        // Load image
        //______________________________________________________________________

        std::cerr << "Processing " << filename << std::endl;
        const pixmap3 source( IMG.load3(filename,NULL) );
        const unit_t  w = source.w;
        const unit_t  h = source.h;

        // save a copy
        PNG.save("source.png", source, NULL);

        //______________________________________________________________________
        //
        // build the histogram
        //______________________________________________________________________
        H.reset();
        histogram::create(hp,source,&server);
        histogram::launch(hp,source,&server);
        H.finish(hp,&server);


        //______________________________________________________________________
        //
        //! find the threshold
        //______________________________________________________________________
        const size_t thr = H.threshold();

        //______________________________________________________________________
        //
        //! keep the background
        //______________________________________________________________________
        pixmap3 target(w,h);
        threshold::apply(target,thr,source,threshold::keep_background);
        PNG.save("bg.png",target,NULL);

        //______________________________________________________________________
        //
        // find the blobs
        //______________________________________________________________________
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

        //______________________________________________________________________
        //
        // scan edges and amplitude
        //______________________________________________________________________
        vector<double> xx(w,as_capacity);
        vector<double> yu(w,as_capacity); //! up value
        vector<double> yd(w,as_capacity); //! down value
        vector<double> ya(w,as_capacity); //! amplitude

        pixmap3 surf(source);
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

            surf[jmin][i] = red;

            unit_t jmax = h-1;
            for(unit_t j=h-1;j>=0;--j)
            {
                if(B[j][i])
                {
                    jmax = j;
                    break;
                }
            }
            surf[jmax][i] = blue;

            yd.push_back(jmin);
            yu.push_back(jmax);
            ya.push_back(jmax-jmin);
            xx.push_back(i);

        }

        string outname = vfs::get_base_name(filename);
        vfs::change_extension(outname, "png");
        outname = "profile_" + outname;
        std::cerr << "Saving in " << outname << std::endl;

        PNG.save(outname,surf,NULL);

        //______________________________________________________________________
        //
        // save amplitude
        //______________________________________________________________________
        const size_t n = w;
        outname = vfs::get_base_name(filename);
        vfs::change_extension(outname, "dat");
        outname = "profile_" + outname;
        std::cerr << "Saving in " << outname << std::endl;
        {
            ios::wcstream fp(outname);
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g\n", xx[i], ya[i]);
            }
        }

        //______________________________________________________________________
        //
        // find the min
        //______________________________________________________________________
        size_t imin= 1;
        double vmin= ya[1];
        for(size_t i=2;i<=n;++i)
        {
            const double tmp = ya[i];
            if(tmp<vmin)
            {
                imin = i;
                vmin = tmp;
            }
        }


        //______________________________________________________________________
        //
        // Prepare Fit
        //______________________________________________________________________
        Soliton               soliton;
        GLS<double>::Function F(&soliton, & Soliton::Eval);
        GLS<double>::Samples  samples;
        const size_t          nvar = 5;
        vector<double>        fa(n);
        GLS<double>::Sample  &sample = samples.append(xx, ya, fa);

        samples.prepare(nvar);

        vector<double> aorg(nvar);
        vector<bool>   used(nvar,false);
        vector<double> aerr(nvar);


        // named variables
        double & start = aorg[1];
        double & slope = aorg[2];
        double & ampli = aorg[3];
        double & scale = aorg[4];
        double & xmin  = aorg[5];


        start = ya[1];
        slope = (ya[n]-ya[1])/(xx[n]-xx[1]);
        xmin  = xx[imin];
        ampli = (start+xmin*slope)-ya[imin];
        scale = 1;

        // compute scaling
        {
            const double ymin = ya[imin];
            const double yl   = ymin + 0.5*(ya[1]-ymin);
            size_t il = imin;
            while(il>1)
            {
                if(ya[il]>=yl)
                {
                    break;
                }
                --il;
            }

            assert(il>=1);

            const double yr = ymin + 0.5*(ya[n]-ymin);
            size_t ir = imin;
            while(ir<n)
            {
                if(ya[ir]>=yr)
                {
                    break;
                }
                ++ir;
            }
            assert(ir<=n);

            const double hw = xx[ir]-xx[il];
            std::cerr << "hw=" << hw << std::endl;
            scale = 2.0/max_of<double>(hw,1);
        }

        (void)sample.computeD2(F,aorg);
        {
            ios::wcstream fp("f0.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g\n", xx[i], fa[i]);
            }
        }

        used[1] = used[2] = used[3] = used[4] = used[5] = true;
        if( ! samples.fit_with(F, aorg, used, aerr ) )
        {
            std::cerr << "Couldn't fit scaling LEVEL-1" << std::endl;
        }

        GLS<double>::display(std::cerr, aorg, aerr);
        {
            ios::wcstream fp("f1.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g\n", xx[i], fa[i]);
            }
        }



#if 0
        //______________________________________________________________________
        //
        // analyze
        //______________________________________________________________________
        const size_t n      = w;


        size_t imin= 1;
        double vmin= ya[1];
        for(size_t i=2;i<=n;++i)
        {
            const double tmp = ya[i];
            if(tmp<vmin)
            {
                imin = i;
                vmin = tmp;
            }
        }
        const double xmin = xx[imin];

        {
            ios::wcstream fp("profile.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", xx[i], yd[i], yu[i]);
            }
        }
        
        string outname = vfs::get_base_name(filename);
        vfs::change_extension(outname, "png");
        outname = "fit_" + outname;
        std::cerr << "Saving in " << outname << std::endl;
        PNG.save(outname,surf,NULL);
#endif
        
        
    }
    
}
YOCTO_PROGRAM_END()
