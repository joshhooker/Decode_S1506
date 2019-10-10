cmake_minimum_required(VERSION 3.9)
project(12be_dp_analysis)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

#Find ROOT
find_package(ROOT REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

set(SCRIPTS
    BeCutIC.root
    Sd1r_AlphaCalibration.txt
    Sd1r_InBeamCalibration.txt
    Sd1s_Calibration.txt
    Sd2r_Calibration.txt
    Sd2s_Calibration.txt
    Yd_Calibration.txt
    Yu_Calibration.txt
    )

foreach(_script ${SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${_script}
        COPYONLY
    )
endforeach()

add_executable(analysis12Be Analysis12Be.cpp)
target_link_libraries(analysis12Be ${ROOT_LIBRARIES})

add_executable(reduce12Be Reduce12Be.cpp)
target_link_libraries(reduce12Be ${ROOT_LIBRARIES})