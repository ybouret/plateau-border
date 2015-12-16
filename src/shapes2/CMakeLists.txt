
ADD_EXECUTABLE(shapes2 shapes2.cpp)
TARGET_LINK_YOCTO(shapes2 mk graphics png z jpeg)

INSTALL( TARGETS shapes2 RUNTIME DESTINATION bin )


