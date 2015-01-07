PROJECT(spatio)

ADD_EXECUTABLE(spatio spatio.cpp)
TARGET_LINK_YOCTO(spatio gfx png z jpeg)

