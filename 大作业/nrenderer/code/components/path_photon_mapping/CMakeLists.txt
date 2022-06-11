cmake_minimum_required(VERSION 3.18)

# 设置名称， 会在"/components”文件夹下生成  名称.dll
set(MY_COMPONENT_NAME "PathPhotonMapping")

file(GLOB_RECURSE COMP_HEADER_FILES "./include/*.h" "./include/*.hpp")
source_group("Header Files" FILES ${COMP_HEADER_FILES})
file(GLOB_RECURSE COMP_SOURCE_FILES "./src/*.cpp")
add_library(${MY_COMPONENT_NAME} SHARED "${COMP_SOURCE_FILES}" "${COMP_HEADER_FILES}" "./src/KDtree.cpp" "./src/shaders/Glass.cpp")
target_link_libraries(${MY_COMPONENT_NAME} NRServer)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${MY_COMPONENT_NAME} OpenMP::OpenMP_CXX)
endif()

include_directories("./include")
