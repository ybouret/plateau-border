PROJECT(spatio)

ADD_EXECUTABLE(shapes shapes.cpp)
TARGET_LINK_YOCTO(shapes mk graphics png z jpeg)

INSTALL( TARGETS shapes RUNTIME DESTINATION bin )


