#include "yocto/graphics/ops/blob.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/tiff.hpp"
#include "yocto/graphics/ops/hist.hpp"
#include "yocto/math/core/tao.hpp"

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
            //std::cerr << "Fitting..." << std::endl;
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
        std::cerr << "Processing " << filename << std::endl;
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

        //______________________________________________________________________
        //
        // scan..
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
            ya.push_back(jmax-jmin);
            xx.push_back(i);

        }

        PNG.save("surf.png",surf,NULL);

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

        vector<double> fu(n), fd(n);

        LeastSquares<double>::Samples samples;
        LeastSquares<double>::Sample::Pointer su( new   LeastSquares<double>::Sample(xx,yu,fu) );
        LeastSquares<double>::Sample::Pointer sd( new   LeastSquares<double>::Sample(xx,yd,fd) );

        samples.push_back(su);
        samples.push_back(sd);

        Soliton soliton;
        LeastSquares<double>::Function F(  &soliton, & Soliton::Eval );
        LeastSquares<double>::Callback cb( &soliton, & Soliton::ToDo );
        LeastSquares<double> Fit;


        const size_t   gvar = 8;
        const size_t   nvar = 5;
        samples.prepare(nvar,gvar);

        vector<double> aorg(gvar);

        double &uStart = aorg[1];
        double &uSlope = aorg[2];
        double &uAmpli = aorg[3];
        double &uCoeff = aorg[4];
        double &uShift = aorg[5];

        uStart = yu[1];
        uSlope = (yu[n]-yu[1])/(xx[n]-xx[1]);
        uAmpli = (uStart+xmin*uSlope) - yu[xmin];
        uCoeff = 1;
        uShift = xmin;

        su->connect(1,1);
        su->connect(2,2);
        su->connect(3,3);
        su->connect(4,4);
        su->connect(5,5);


        double       &dStart = aorg[6];
        double       &dSlope = aorg[7];
        double       &dAmpli = aorg[8];
        //const double &dCoeff = uCoeff;
        //const double &dShift = uShift;

        dStart = yd[1];
        dSlope = (yd[n]-yd[1])/(xx[n]-xx[1]);
        dAmpli = (dStart+xmin*dSlope) - yd[xmin];

        sd->connect(1,6);
        sd->connect(2,7);
        sd->connect(3,8);
        sd->connect(4,4);
        sd->connect(5,5);


        vector<bool>   used(gvar,false);
        vector<double> aerr(gvar,0);

        // Level 1: adjust scaling
        std::cerr << "LEVEL 1" << std::endl;
        used[4] = true;
        if(!Fit(samples,F,aorg,used,aerr,&cb))
        {
            std::cerr << "Couldn't fit Level-1" << std::endl;
            continue;
        }
        Fit.display(std::cerr, aorg, aerr);

        {
            ios::wcstream fp("f1.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", xx[i], fd[i], fu[i]);
            }
        }

        // Level 2: adjust bumps
        std::cerr << "LEVEL 2" << std::endl;
        used[3] = true;
        used[5] = true;
        used[8] = true;

        if(!Fit(samples,F,aorg,used,aerr,&cb))
        {
            std::cerr << "Couldn't fit Level-2" << std::endl;
            continue;
        }
        Fit.display(std::cerr, aorg, aerr);
        {
            ios::wcstream fp("f2.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", xx[i], fd[i], fu[i]);
            }
        }

        // Level 3: full fit
        std::cerr << "LEVEL 3" << std::endl;
        tao::ld(used, true);
        if(!Fit(samples,F,aorg,used,aerr,&cb))
        {
            std::cerr << "Couldn't fit Level-3" << std::endl;
            continue;
        }
        Fit.display(std::cerr, aorg, aerr);
        {
            ios::wcstream fp("f3.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", xx[i], fd[i], fu[i]);
            }
        }

        surf.copy(source);
        for(size_t i=1;i<=n;++i)
        {
            const unit_t x(xx[i]);
            const unit_t d(Floor(fd[i]+0.5));
            const unit_t u(Floor(fu[i]+0.5));
            surf[d][x] = named_color::get("magenta");
            surf[u][x] = named_color::get("orange");
        }
        PNG.save("sfit.png",surf,NULL);

    }

}
YOCTO_PROGRAM_END()
