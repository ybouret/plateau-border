ADD_EXECUTABLE(genshape genshape.cpp)
TARGET_LINK_YOCTO(genshape mk graphics png z jpeg)

INSTALL( TARGETS genshape RUNTIME DESTINATION bin )


