#include "yocto/program.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/point3d.hpp"
#include "yocto/code/bzset.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/sequence/vector.hpp"
#include "yocto/container/matrix.hpp"
#include "yocto/spade/format/stl.hpp"

using namespace yocto;
using namespace math;
using namespace spade;

static const double beta  = 120.0;
static const double alpha = 180.0 - beta;
static size_t       n     = 20;


static inline
size_t points_per_arch()
{
    return 3*n;
}

static inline double compute_delta(double R)
{
    return R/sin( Deg2Rad(beta/2.0) );
}


typedef point3d<double> vertex;

static inline
void compute_centers(const double R, vertex &C1, vertex &C2, vertex &C3)
{
    bzset(C1);
    bzset(C2);
    bzset(C3);

    const double delta = compute_delta(R);
    C1.x = 0;
    C1.y = delta;

    const double angle = 90-beta/2;
    C2.x =  delta * cos( Deg2Rad(angle) );
    C2.y = -delta * sin( Deg2Rad(angle) );

    C3.x = -C2.x;
    C3.y =  C2.y;
}

static inline
void generate_arche(const double R, vector<vertex> &v, const double z)
{
    vertex C1,C2,C3;
    compute_centers(R,C1,C2,C3);
    v.free();

    const double lo = -alpha/2;

    for(size_t i=0;i<=n;++i)
    {
        const double angle = lo + (i*alpha)/n;
        const double aa    = Deg2Rad(angle);
        const double ca    = cos( aa );
        const double sa    = sin( aa );

        vertex r;
        r.z = z;
        r.x = C1.x + R * sa;
        r.y = C1.y - R * ca;
        v.push_back(r);
    }

    const double b = Deg2Rad(beta);
    double cb = cos( b );
    double sb = sin( b );

    for(size_t i=1;i<=n;++i)
    {
        const double angle = lo + (i*alpha)/n;
        const double aa    = Deg2Rad(angle);
        const double ca    = cos( aa );
        const double sa    = sin( aa );
        vertex r;
        r.z = z;
        r.x = C1.x + R * sa;
        r.y = C1.y - R * ca;

        vertex u;
        u.x =  r.x * cb + r.y*sb;
        u.y = -r.x * sb + r.y*cb;
        u.z = z;
        v.push_back(u);
    }

    const double bb = b+b;
    cb = cos( bb );
    sb = sin( bb );

    for(size_t i=1;i<n;++i)
    {
        const double angle = lo + (i*alpha)/n;
        const double aa    = Deg2Rad(angle);
        const double ca    = cos( aa );
        const double sa    = sin( aa );
        vertex r;
        r.z = z;
        r.x = C1.x + R * sa;
        r.y = C1.y - R * ca;

        vertex u;
        u.x =  r.x * cb + r.y*sb;
        u.y = -r.x * sb + r.y*cb;
        u.z = z;
        v.push_back(u);
    }




}




double am    = 0.5;
double scale = 1;

YOCTO_PROGRAM_START()
{
    vertex C1,C2,C3;

    compute_centers(1.0, C1, C2, C3);


    {
        ios::wcstream fp("centers.dat");
        fp("%g %g\n", C1.x, C1.y );
        fp("%g %g\n", C2.x, C2.y );
        fp("%g %g\n", C3.x, C3.y );
    }

    vector<vertex> v;
    generate_arche(1.0,v,0.0);
    {
        ios::wcstream fp("arche.dat");
        for(size_t i=1;i<=v.size();++i)
        {
            fp("%g %g\n", v[i].x, v[i].y);
        }
        fp("%g %g\n", v[1].x, v[1].y);
    }
    std::cerr << "#v=" << v.size() << std::endl;


    const size_t nr = points_per_arch();
    const size_t nz = 200;

    const double Lz = 3.0;
    const vertex inside(0,0,0);
    const vertex bot(0,0,-Lz);
    const vertex top(0,0,Lz);

    matrix<vertex> shape(nz,nr);
    for(size_t i=1;i<=nz;++i)
    {
        const double z = -Lz + ((2.0*Lz)*(i-1))/double(nz-1);
        const double R = 1.0 - (1.0-am) / Square(cosh(scale*z));
        v.free();
        generate_arche(R,v,z);
        assert(v.size()==nr);
        for(size_t j=1;j<=nr;++j)
        {
            shape[i][j] = v[j];
        }
    }

    typedef stl::facet<double> facet_t;

    vector<facet_t> facets;

    stl::close_contour(facets, shape[1], bot, inside);

    for(size_t i=1;i<nz;++i)
    {
        stl::make_ribbon(facets, shape[i], shape[i+1], inside);
    }

    stl::close_contour(facets, shape[nz], top, inside);
    {
        ios::wcstream fp("shape.stl");
        stl::save_binary(fp,facets);
    }

}
YOCTO_PROGRAM_END()


