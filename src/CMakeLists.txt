cmake_minimum_required (VERSION 2.6)

project(TPSC)

add_executable(TPSC generique.cpp hamiltonien.cpp green0.cpp chi0.cpp chi.cpp green.cpp cond_opt.cpp tpsc_main.cpp)

target_compile_options(TPSC PRIVATE -std=c++11 -fopenmp -O0)

target_include_directories(TPSC PRIVATE /usr/local/opt/llvm/include /usr/local/Cellar/llvm/7.0.0/lib/clang/7.0.0/include)

target_link_libraries(TPSC fftw3 blas lapack omp)

install(TARGETS TPSC DESTINATION bin)

