#include "yocto/graphics/ops/blob.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/tiff.hpp"
#include "yocto/graphics/ops/hist.hpp"

#include "yocto/program.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/fit/lsf.hpp"


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
            return a[1] + x*a[2] - a[3]/(csha*csha);
        }


        inline bool ToDo(LeastSquares<double>::Function       &F,
                         const LeastSquares<double>::Samples  &S)
        {
            return true;
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
        vector<double> xx(w,as_capacity);
        vector<double> yy(w,as_capacity);
        vector<double> yu(w,as_capacity);
        vector<double> yd(w,as_capacity);

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
            yd.push_back(jmin);
            yu.push_back(jmax);
            yy.push_back(jmax-jmin);
            xx.push_back(i);
        }
        PNG.save("surf.png",surf,NULL);
        {
            ios::wcstream fp("thick.dat");
            for(unit_t i=1;i<=w;++i)
            {
                fp("%g %g\n", xx[i], (yy[i]));
            }
        }

        //______________________________________________________________________
        //
        // scan..
        //______________________________________________________________________
        
        vector<double> yf(w);

        Soliton              soliton;
        LeastSquares<double> Fit;
        LeastSquares<double>::Function F(  &soliton, & Soliton::Eval );
        LeastSquares<double>::Callback CB( &soliton, & Soliton::ToDo );
        LeastSquares<double>::Samples  samples;
        samples.append(xx, yy, yf);

        const size_t   nv = 5;
        vector<double> aorg(nv,0);
        vector<bool>   used(nv,false);
        vector<double> aerr(nv,0);

        samples.prepare(nv);

        unit_t imin = 1;
        double vmin = yy[1];
        for(unit_t i=2;i<=w;++i)
        {
            const double tmp = yy[i];
            if(tmp<vmin)
            {
                imin=i;
                vmin=tmp;
            }
        }

        aorg[1] = yy[1];
        aorg[2] = (yy[w]-yy[1])/double(w);
        aorg[3] = (aorg[1] + imin*aorg[2]) - vmin;
        aorg[4] = 1; //4.0/w;
        aorg[5] = imin;


        std::cerr << "level-1" << std::endl;
        used[4] = true;
        if(!Fit(samples,F,aorg,used,aerr,&CB))
        {
            std::cerr << "couldn't fit level-1" << std::endl;
            continue;
        }


        Fit.display(std::cerr, aorg, aerr);
        {
            ios::wcstream fp("f1.dat");
            for(unit_t i=1;i<=w;++i)
            {
                fp("%g %g\n", xx[i], yf[i]);
            }
        }

        std::cerr << "level-2" << std::endl;
        used[3] = used[5] = true;
        if(!Fit(samples,F,aorg,used,aerr,&CB))
        {
            std::cerr << "couldn't fit level-2" << std::endl;
            continue;
        }
        Fit.display(std::cerr, aorg, aerr);
        {
            ios::wcstream fp("f2.dat");
            for(unit_t i=1;i<=w;++i)
            {
                fp("%g %g\n", xx[i], yf[i]);
            }
        }

        std::cerr << "level-3" << std::endl;
        used[1]=used[2] = true;
        if(!Fit(samples,F,aorg,used,aerr,&CB))
        {
            std::cerr << "couldn't fit level-2" << std::endl;
            continue;
        }
        Fit.display(std::cerr, aorg, aerr);
        {
            ios::wcstream fp("f3.dat");
            for(unit_t i=1;i<=w;++i)
            {
                fp("%g %g\n", xx[i], yf[i]);
            }
        }




    }
    
}
YOCTO_PROGRAM_END()
