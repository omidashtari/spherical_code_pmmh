cmake_minimum_required(VERSION 3.0)

project(spherical_code_pmmh_project Fortran)

# Set gfortran as default compiler
if(NOT DEFINED CMAKE_Fortran_COMPILER)
    message(STATUS "No Fortran compiler specified, defaulting to gfortran")
    set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "Fortran compiler" FORCE)
endif()

# Set release as default compilation
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (default Release)" FORCE)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "Using GNU Fortran compiler")
    set(CMAKE_Fortran_FLAGS_DEBUG    "-O0 -fdefault-real-8 -fdefault-double-8 -g -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_RELEASE  "-O3 -fdefault-real-8 -fdefault-double-8")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    message(STATUS "Using Intel Fortran compiler")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -real-size 64 -debug all -check all -traceback -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -real-size 64")

# elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
#     message(STATUS "Using NVIDIA HPC (nvfortran) compiler")
#     set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -r8 -g -Mbounds -Mchkptr -Mchkstk -traceback")
    # set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -r8 -Minfo=all")

endif()

# --- MKL / BLAS / LAPACK selection ---

# User override to disable MKL
option(USE_MKL "Use Intel MKL if available" ON)

set(HAVE_MKL FALSE)

# if(USE_MKL AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
if(USE_MKL)
    message(STATUS "Attempting to find MKL...")
    find_package(MKL QUIET)
    if(MKL_FOUND)
        set (MKLROOT "$ENV{MKLROOT}")
        message(STATUS "MKL found in $ENV{MKLROOT} and will be used")
        set(HAVE_MKL TRUE)
        find_library(MKLCORE mkl_core HINTS ${MKL_ROOT}/lib/intel64)
        find_library(MKLSEQ mkl_sequential HINTS ${MKL_ROOT}/lib/intel64)
        if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
            find_library(MKL64 mkl_intel_lp64 HINTS ${MKL_ROOT}/lib/intel64)
        endif()
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            find_library(MKL64 mkl_gf_lp64 HINTS ${MKL_ROOT}/lib/intel64)
        endif()
    else()
        message(WARNING "MKL requested but not found; falling back to LAPACK/BLAS")
    endif()
endif()

if(NOT HAVE_MKL)

    message(STATUS "Looking for BLAS and LAPACK")

    # --- BLAS ---
    if(USER_BLAS_PATH)
        message(STATUS "Using user-specified BLAS path: ${USER_BLAS_PATH}")
        find_library(BLAS_LIB NAMES blas PATHS ${USER_BLAS_PATH} REQUIRED)
    else()
        find_package(BLAS QUIET)
        if(NOT BLAS_FOUND)
            message(FATAL_ERROR "BLAS library not found. Please specify -DUSER_BLAS_PATH=/path/to/blas")
        endif()
        set(BLAS_LIB ${BLAS_LIBRARIES})
    endif()

    # --- LAPACK ---
    if(USER_LAPACK_PATH)
        message(STATUS "Using user-specified LAPACK path: ${USER_LAPACK_PATH}")
        find_library(LAPACK_LIB NAMES lapack PATHS ${USER_LAPACK_PATH} REQUIRED)
    else()
        find_package(LAPACK QUIET)
        if(NOT LAPACK_FOUND)
            message(FATAL_ERROR "LAPACK library not found. Please specify -DUSER_LAPACK_PATH=/path/to/lapack")
        endif()
        set(LAPACK_LIB ${LAPACK_LIBRARIES})
    endif()

    # Add gfortran runtime if Intel compiler linking against GNU-built BLAS/LAPACK
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        list(APPEND EXTRA_LINK_LIBS gfortran)
    endif()
endif()

# Optional user inputs
set(USER_FFTW_PATH "" CACHE PATH "Path to FFTW library")
set(USER_SHTNS_PATH "" CACHE PATH "Path to SHTNS library")

# FFTW
if(USER_FFTW_PATH)
    message(STATUS "Using user-specified FFTW path: ${USER_FFTW_PATH}")
    find_library(FFTW_LIB NAMES fftw3 PATHS ${USER_FFTW_PATH} REQUIRED)
else()
    find_library(FFTW_LIB NAMES fftw3)
endif()

if(NOT FFTW_LIB)
    message(FATAL_ERROR "FFTW library not found. Please specify -DUSER_FFTW_PATH=/path/to/fftw")
endif()

# SHTNS
if(USER_SHTNS_PATH)
    message(STATUS "Using user-specified SHTNS path: ${USER_SHTNS_PATH}")
    find_library(SHTNS_LIB NAMES shtns PATHS ${USER_SHTNS_PATH} REQUIRED)
else()
    find_library(SHTNS_LIB NAMES shtns)
endif()

if(NOT SHTNS_LIB)
    message(FATAL_ERROR "SHTNS library not found. Please specify -DUSER_SHTNS_PATH=/path/to/shtns")
endif()

# Collect module files
file(GLOB MOD_FILES src/mods/*.f90)

# Add executable
add_executable(spherical_code_pmmh
    src/spherical_code_pmmh.f90
    ${MOD_FILES}
)

# Include directories
target_include_directories(spherical_code_pmmh PRIVATE
    src/mods
    src/include
)

# Link libraries
if(HAVE_MKL)
    target_link_libraries(spherical_code_pmmh
        ${MKL64}
        ${MKLSEQ}
        ${MKLCORE}
        pthread
        dl
        ${SHTNS_LIB}
        ${FFTW_LIB}
        m
        c
    )
else()
    target_link_libraries(spherical_code_pmmh
        ${LAPACK_LIB}
        ${BLAS_LIB}
        ${EXTRA_LINK_LIBS}
        ${SHTNS_LIB}
        ${FFTW_LIB}
        m
        c
    )
endif()

