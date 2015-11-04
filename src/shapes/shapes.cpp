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
        const double _yLeft  = yy[1];
        const double _yRight = yy[n];
        const double _ySlope = (_yRight-_yLeft)/width;
        const double _mLeft  = ym[1];
        const double _mRight = ym[n];
        const double _mSlope = (_mRight-_mLeft)/width;

        std::cerr << _yLeft << "+(" << _ySlope << ")*x, " << _mLeft << "+(" << _mSlope << ")*x" << std::endl;

        size_t imin=1;
        double vmin= _yLeft;
        for(size_t i=2;i<=n;++i)
        {
            const double tmp = yy[i];
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
        LeastSquares<double> Fit;
        
        //______________________________________________________________________
        //
        // Preparing all variables
        //______________________________________________________________________
        const size_t nvar = 5; //!< variables for one curve
        const size_t gvar = 8; //!< 2*5 -2 shared variables

        samples.prepare(nvar,gvar);

        vector<double> aorg(gvar);
        double &yStart = aorg[1];
        double &ySlope = aorg[2];
        double &yAmpli = aorg[3];
        double &yScale = aorg[4];
        double &yShift = aorg[5];
        
        //double &mScale = aorg[4]; //!< same scaling
        //double &mShift = aorg[5]; //!< same shift
        
        double &mStart = aorg[6]; //!< start
        double &mSlope = aorg[7];
        double &mAmpli = aorg[8];

        
        // link global to local
        pY->connect(1, 1);
        pY->connect(2, 2);
        pY->connect(3, 3);
        pY->connect(4, 4);
        pY->connect(5, 5);

        
        pM->connect(1,6);
        pM->connect(2,7);
        pM->connect(3,8);
        pM->connect(4,4);
        pM->connect(5,5);
        
        // initialize values
        yStart = _yLeft;
        ySlope = _ySlope;
        yAmpli = (yStart + xmin * ySlope) - yy[imin];
        yScale = 1;
        yShift = xmin;
        
        mStart = _mLeft;
        mSlope = _mSlope;
        mAmpli = (mStart+xmin*mSlope) - ym[imin];
        
        vector<bool>   used(gvar,false);
        vector<double> aerr(gvar);
        used[3] = true;
        used[4] = true;
        used[5] = true;
        
        
        if(!Fit(samples,F,aorg,used,aerr,&cb))
        {
            std::cerr << "Couldn't Fit Level-1" << std::endl;
        }
        
        Fit.display(std::cerr, aorg, aerr);
        {
            ios::wcstream fp("f1.dat");
            for(size_t i=1;i<=n;++i)
            {
                 fp("%g %g %g %g %g\n", xx[i], yy[i], ym[i], fy[i], fm[i]);
            }
        }
        
        
        
    }
    
}
YOCTO_PROGRAM_END()
