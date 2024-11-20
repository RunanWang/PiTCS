cmake_minimum_required(VERSION 3.22)
project(PiTCS)

set(CMAKE_CXX_STANDARD 17)

# Unify using clang as compiler
if (APPLE)
    SET(CMAKE_C_COMPILER /opt/homebrew/opt/llvm/bin/clang++)
    SET(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)
else ()
    # SET(CMAKE_C_COMPILER /usr/bin/clang++)
    # SET(CMAKE_CXX_COMPILER /usr/bin/clang++)
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -stdlib=libstdc++")
endif ()

# Debug/Release
#set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Find gflags
find_package(gflags REQUIRED)
set(FLAG_LIBRARIES gflags)

# Find glog
find_package(glog REQUIRED)
set(LOG_LIBRARIES glog::glog)

find_package(fmt REQUIRED)
set(FMT_LIBRARIES fmt::fmt)
#set(CMAKE_CXX_FLAGS "-fsanitize=address")
# Find openmp
if (APPLE)
    include_directories("/opt/homebrew/opt/llvm/include" "/opt/homebrew/opt/libomp/include")
    link_directories("/opt/homebrew/opt/llvm/lib" "/opt/homebrew/opt/libomp/lib")
    set(OpenMP_CXX_FLAGS "-fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else ()
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(${CMAKE_CXX_FLAGS})
endif ()

add_executable(PiTCS
        main.cpp constant.h
        Utils/Timer.cpp Utils/Timer.h Utils/Memory.h Utils/logging.h

        DataStructure/Graph.h DataStructure/Graph.cpp
        DataStructure/MLGWithSchema.cpp DataStructure/MLGWithSchema.h

        DataStructure/Bucket.cpp DataStructure/Bucket.h
        DataStructure/EqualTree.cpp DataStructure/EqualTree.h DataStructure/UnionFind.cpp DataStructure/UnionFind.h
        DataStructure/PivotalCombo.h DataStructure/PivotalCombo.cpp

        Algorithm/Algo.cpp Algorithm/Algo.h
        Algorithm/PivotalTruss.cpp Algorithm/PivotalTruss.h
        test/testHeader.h
        Algorithm/TriangleConnTruss.cpp
        Algorithm/TriangleConnTruss.h
        test/SOTATestPi.cpp
)

target_link_libraries(PiTCS
        ${FLAG_LIBRARIES}
        ${LOG_LIBRARIES}
        ${FMT_LIBRARIES}
        )