#include "yocto/program.hpp"
#include "yocto/math/types.hpp"
#include "yocto/math/trigconv.hpp"
#include "yocto/math/point3d.hpp"
#include "yocto/code/bzset.hpp"

#include "yocto/ios/ocstream.hpp"
#include "yocto/sequence/vector.hpp"

using namespace yocto;
using namespace math;

static const double beta  = 120.0;
static const double alpha = 180.0 - beta;
static size_t       n     = 20;

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
        v.push_back(u);
    }




}






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


}
YOCTO_PROGRAM_END()

