project(VHFDummy)

set(INCLUDE_DIRS
    ${NASSP_SOURCE_ROOT}/src_sys
    ${ORBITERSDK_DIR}/include
    ${ORBITERSDK_DIR}/XRSound
)
set(SOURCES
    ${NASSP_SOURCE_ROOT}/src_sys/connector.cpp
    ${NASSP_SOURCE_ROOT}/src_sys/vhfdummy.cpp
    ${NASSP_SOURCE_ROOT}/src_sys/VHFDummyConnector.cpp
)

link_directories(${ORBITERSDK_DIR}/lib)
link_directories(${ORBITERSDK_DIR}/XRSound)
link_directories(${CMAKE_BINARY_DIR}/Projects/PanelSDK)

add_library(${PROJECT_NAME} SHARED ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE
    Orbiter
    Orbitersdk
    PanelSDK
    XRSound
)
