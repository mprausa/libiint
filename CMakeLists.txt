cmake_minimum_required(VERSION 3.0)
project(libiint)
set(CMAKE_CXX_STANDARD 11)

option(BUILD_MMA_INTERFACE "Build Mathematica interface" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

foreach(p LIB INCLUDE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

find_package(yaml-cpp REQUIRED)
find_package(FLINT REQUIRED)
find_package(ARB REQUIRED)
find_package(GMP REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${YAML_CPP_INCLUDE_DIR})
include_directories(${FLINT_INCLUDE_DIR})
include_directories(${ARB_INCLUDE_DIR})

file(GLOB SOURCES src/arb/*.cpp src/iint/*.cpp)
add_library(iint STATIC ${SOURCES})

install(DIRECTORY include/ DESTINATION "${INSTALL_INCLUDE_DIR}" PATTERN "iint/*.h" PATTERN "arb/*.h")
install(TARGETS iint DESTINATION "${INSTALL_LIB_DIR}")

if(BUILD_MMA_INTERFACE)
    find_package(Mathematica)
    if(NOT INSTALL_MMA_DIR)
        set(INSTALL_MMA_DIR  "${Mathematica_USERBASE_DIR}/Applications/IInt")
    endif()

    file(GLOB MMA_SOURCES mma/interface/*.cpp)

    include_directories(${Mathematica_MathLink_INCLUDE_DIR})

    Mathematica_MathLink_ADD_EXECUTABLE(IInt mma/interface/MathLink.tmpp ${MMA_SOURCES})

    target_link_libraries(IInt iint ${FLINT_LIBRARIES} ${ARB_LIBRARIES} ${YAML_CPP_LIBRARIES} ${GMP_LIBRARIES})

    install(TARGETS IInt DESTINATION "${INSTALL_MMA_DIR}/interface")
    install(FILES mma/IInt.m DESTINATION "${INSTALL_MMA_DIR}")
endif(BUILD_MMA_INTERFACE)

