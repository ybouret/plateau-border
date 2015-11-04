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

        //______________________________________________________________________
        //
        // scan..
        //______________________________________________________________________
        vector<double> xx(w,as_capacity);
        vector<double> yy(w,as_capacity);
        vector<double> yu(w,as_capacity);
        vector<double> yd(w,as_capacity);
        vector<double> ym(w,as_capacity);

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
            yy.push_back(jmax-jmin);
            ym.push_back(0.5*(yu.back()+yd.back()));
            xx.push_back(i);

            surf[unit_t(ym.back())][i] = named_color::get("green");
        }

        PNG.save("surf.png",surf,NULL);

        //______________________________________________________________________
        //
        // analyze
        //______________________________________________________________________
        const size_t n      = w;
        const double width  = xx[n]-xx[1];
        const double aLeft  = yy[1];
        const double aRight = yy[n];
        const double aSlope = (aRight-aLeft)/width;
        const double mLeft  = ym[1];
        const double mRight = ym[n];
        const double mSlope = (mRight-mLeft)/width;

        std::cerr << aLeft << "+(" << aSlope << ")*x, " << mLeft << "+(" << mSlope << ")*x" << std::endl;

        size_t imin=1;
        double vmin=aLeft;
        for(size_t i=2;i<=n;++i)
        {
            const double tmp = yy[i];
            if(tmp<vmin)
            {
                imin = i;
                vmin = tmp;
            }
        }

        {
            ios::wcstream fp("profile.dat");
            for(size_t i=1;i<=n;++i)
            {
                fp("%g %g %g\n", xx[i], yy[i], ym[i]);
            }
        }

        vector<double> fy(n), fm(n);

        LeastSquares<double>::Samples samples;
        LeastSquares<double>::Sample::Pointer pY( new   LeastSquares<double>::Sample(xx,yy,fy) );
        LeastSquares<double>::Sample::Pointer pM( new   LeastSquares<double>::Sample(xx,ym,fm) );

        samples.push_back(pY);
        samples.push_back(pM);

        Soliton soliton;
        LeastSquares<double>::Function F(  &soliton, & Soliton::Eval );
        LeastSquares<double>::Callback cb( &soliton, & Soliton::ToDo );

        //______________________________________________________________________
        //
        // Preparing all variables
        //______________________________________________________________________
        const size_t nvar = 5; //!< variables for one curve
        const size_t gvar = 2*nvar - 1;

        samples.prepare(nvar,gvar);

        vector<double> aorg(gvar);
        double &yStart = aorg[1];
        double &ySlope = aorg[2];
        

    }
    
}
YOCTO_PROGRAM_END()
