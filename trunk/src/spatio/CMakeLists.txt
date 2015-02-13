PROJECT(spatio)

ADD_EXECUTABLE(spatio spatio.cpp)
TARGET_LINK_YOCTO(spatio mk gfx png z jpeg)

