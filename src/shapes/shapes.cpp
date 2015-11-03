#include "yocto/graphics/ops/blob.hpp"
#include "yocto/graphics/image/jpeg.hpp"
#include "yocto/graphics/image/png.hpp"
#include "yocto/graphics/image/tiff.hpp"

#include "yocto/program.hpp"

using namespace yocto;
using namespace graphics;

YOCTO_PROGRAM_START()
{

    image &IMG = image::instance();
    
    IMG.declare( new png_format()  );
    IMG.declare( new jpeg_format() );
    IMG.declare( new tiff_format() );
    
    //const image::format &PNG = IMG["PNG"];
    
    for(int arg=1;arg<argc;++arg)
    {
        const string  filename = argv[arg];
        const pixmap3 source( IMG.load3(filename,NULL) );
    }
    
}
YOCTO_PROGRAM_END()
