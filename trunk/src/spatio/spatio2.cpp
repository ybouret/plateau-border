#include "yocto/gfx/image.hpp"
#include "yocto/fs/local-fs.hpp"
#include "yocto/gfx/image/png.hpp"
#include "yocto/gfx/image/jpeg.hpp"
#include "yocto/exception.hpp"
#include "yocto/ptr/auto.hpp"
#include "yocto/gfx/ops/blob.hpp"
#include "yocto/gfx/ops/hist.hpp"
#include "yocto/gfx/ops/contrast.hpp"
#include "yocto/ptr/arc.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/ios/ocstream.hpp"
#include "yocto/math/wavelet/cwt.hpp"
#include "yocto/threading/server.hpp"

#include "yocto/program.hpp"

using namespace yocto;
using namespace gfx;
using namespace math;

static inline double Sech(const double x) throw()
{
    return 1.0/Square(cosh(x));
}

class Frame : public counted_object
{
public:
    typedef arc_ptr<Frame> Pointer;

    bitmap::pointer bmp;  //!< original bitmap, rgb
    pixmap3         img;  //!< using shared bmp
    pixmapf         pgs;  //!< pixmap greyscale
    pixmapf         mask; //!< thresholded data
    vector<double>  x;    //!< in pixels
    vector<double>  y;    //!< in pixels
    const string   &outdir; //!< common output fir
    const string    epname; //!< original entry pointer name
    const string   &savext;
    size_t          xmin;
    size_t          xmax;
    size_t          length;
    double          alpha;

    explicit Frame(bitmap       *user_bmp,
                   const string &user_outdir,
                   const char   *user_epname,
                   const string &user_savext) :
    bmp(user_bmp),
    img(bmp,NULL),
    pgs(img,rgb2gsf<rgb_t>),
    mask(pgs.w,pgs.h),
    x(mask.w,as_capacity),
    y(mask.w,as_capacity),
    outdir( user_outdir ),
    epname( user_epname ),
    savext(user_savext),
    xmin(0),
    xmax(img.w-1),
    length(img.w),
    alpha(0)
    {
        // processing

        // enhance gresycale contrast
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

        //build mask
        size_t level = 0;
        if(cls.size()>0)
        {
            level = cls.front()->uuid;
        }

        const size_t w = mask.w;
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

        //analyze mask to get shape
        for(size_t i=0;i<w;++i)
        {
            unit_t lo=0;
            unit_t hi=h;
            for(;lo<h;++lo)
            {
                if(mask[lo][i]>0)
                    break;
            }


            for(--hi;hi>=0;--hi)
            {
                if(mask[hi][i]>0)
                    break;
            }

            x.push_back(i);
            y.push_back((hi-lo)+1);
        }

        const image &IMG = image::instance();
        string outname = outdir + epname;
        vfs::change_extension(outname, savext);
        std::cerr << "\tsaving to " << outname << std::endl;
        IMG["PNG"].save(outname,mask, image::get_rampf, NULL, NULL);

    }

    virtual ~Frame() throw()
    {
    }

    static inline void PutRGB(void *addr, const rgba_t &c, const void *) throw()
    {
        new (addr) rgb_t(c.r,c.g,c.b);
    }

    inline void trim(size_t user_xmin, size_t user_xmax) throw()
    {
        assert(xmin<=xmax);
        assert(xmax<length);
        xmin = user_xmin;
        xmax = user_xmax;

        const size_t nback = length - (xmax+1);
        for(size_t i=0;i<nback;++i)
        {
            x.pop_back();
            y.pop_back();
        }

        for(size_t i=0;i<xmin;++i)
        {
            x.pop_front();
            y.pop_front();
        }
        length = xmax+1-xmin;
    }

    class Job
    {
    public:
        Frame    &frame;
        lockable &access;
        inline  Job(Frame    &user_frame,
                    lockable &user_access
                    ) :
        frame(user_frame),
        access(user_access)
        {}

        inline ~Job() throw() {}

        Job(const Job &J ) :
        frame(J.frame),
        access(J.access)
        {

        }

        inline void operator()()
        {
            {
                scoped_lock guard(access);
                std::cerr << ".";
            }
            const size_t n = frame.length;

            numeric<double>::function Psi( cfunctor(Sech) );
            vector<double> shifts(n,0);
            vector<double> scales(n,0);
            matrix<double> W;
            wavelet<double>::cwt(frame.x, frame.y, Psi, shifts, scales, W);
            wavelet<double>::rescale(W);

            string outname = frame.outdir + "w_" + frame.epname;
            vfs::change_extension(outname, frame.savext);

            //create a picture
            const size_t res_w = frame.bmp->w;
            const size_t res_h = frame.bmp->h+n;
            pixmap3      res(res_w,res_h);

            // colorize wavelets
            size_t offset = frame.bmp->h;
            for(size_t i=1;i<=n;++i)
            {
                for(size_t j=1;j<=n;++j)
                {
                    const float w(W[i][j]);
                    rgb_t C;
                    if(w>=0)
                    {
                        C.r = conv::to_byte(w);
                    }
                    else
                    {
                        C.g = conv::to_byte(-w);
                    }
                    //const rgb_t  C = rgb_t::make_ramp(w, 0, 1);
                    res[(j-1)+offset][(i-1)+frame.xmin] = C;
                }
            }

            //duplicate original drawing
            for(size_t j=0;j<frame.img.h;++j)
            {
                for(size_t i=0;i<frame.img.w;++i)
                {
                    res[j][i] = frame.img[j][i];
                }
            }

            // save the picture
            const image &IMG = image::instance();

            IMG["PNG"].save(outname, res, image::get_rgb_dup, NULL, NULL);

        }

    private:
        YOCTO_DISABLE_ASSIGN(Job);
    };


private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Frame);
};

class Frames : public vector<Frame::Pointer>
{
public:
    explicit Frames() : vector<Frame::Pointer>(1000,as_capacity)
    {
    }

    virtual ~Frames() throw()
    {
    }

private:
    YOCTO_DISABLE_COPY_AND_ASSIGN(Frames);
};

#include "yocto/string/conv.hpp"

YOCTO_PROGRAM_START()
{
    if(argc<=3)
        throw exception("usage: %s directory extensions output_directory [xmin xmax]",program);

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
    if(inpdir==outdir)
    {
        throw exception("same directories");
    }
    imgext.to_lower();

    image &IMG = image::instance();
    IMG.declare( new png_format()  );
    IMG.declare( new jpeg_format() );


    vfs & fs = local_fs::instance();

    // create output dir
    fs.create_dir(outdir,true);
    fs.remove_files_with_extension_in(outdir, savext);

    // scan input directory
    auto_ptr<vfs::scanner> scan( fs.new_scanner( inpdir ) );
    Frames frames;

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

        // load, preprocess and store frame
        Frame               *frame = new Frame(IMG.load(path,3,Frame::PutRGB,NULL,NULL),outdir,ep->base_name,savext);
        const Frame::Pointer fptr(frame);
        frames.push_back(fptr);


        // check consistency
        if(frames.size()>1)
        {
            if( frames.back()->length != frames.front()->length )
            {
                throw exception("inconsistent width");
            }
        }

        //if(frames.size()>=10) break;

    }

    const size_t num_frames = frames.size();
    if(num_frames>0)
    {
        const size_t length = frames[1]->length;
        if(xmax<0||xmax>=length-1) xmax = length-1;
        if(xmin>xmax)              xmin = xmax;

        assert(xmin>=0);
        assert(xmax<length);
        assert(xmin<=xmax);

        std::cerr << "Study on [" << xmin << "," << xmax << "]" << std::endl;
        for(size_t i=1;i<=num_frames;++i)
        {
            frames[i]->trim(xmin,xmax);
        }
        
        threading::server srv;
        for(size_t i=1;i<=num_frames;++i)
        {
            const Frame::Job J( *frames[i], srv.access);
            const threading::server::job todo(J);
            srv.enqueue(todo);
        }
        srv.flush();
        
    }
    
    
}
YOCTO_PROGRAM_END()
