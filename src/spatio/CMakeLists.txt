PROJECT(spatio)

ADD_EXECUTABLE(spatio spatio.cpp)
TARGET_LINK_YOCTO(spatio mk gfx png z jpeg)

ADD_EXECUTABLE(spatio2 spatio2.cpp)
TARGET_LINK_YOCTO(spatio2 mk gfx png z jpeg)

INSTALL( TARGETS spatio spatio2 RUNTIME DESTINATION bin )

