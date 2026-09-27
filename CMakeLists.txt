cmake_minimum_required(VERSION 3.18)

project(RESIDUUM
	VERSION 0.1.0
	DESCRIPTION "Residuum - A Modular Finite Element Solver"
	LANGUAGES CXX
)

# Set global output dirs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Cpp Standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

# build type and basic compiler summary
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

# Build Options
option(BUILD_TESTS "Build unit and integration tests" ON)
option(BUILD_APPLICATIONS "Build CLI application executables" ON)
option(BUILD_BENCHMARKS "Build benchmark executables" ON)

# Default compiler options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
	add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# AddressSanitizer + UndefinedBehaviorSanitizer
option(ENABLE_SANITIZERS "Build with AddressSanitizer + UndefinedBehaviorSanitizer" OFF)

if(ENABLE_SANITIZERS)
	if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
		add_compile_options(-fsanitize=address,undefined -fno-omit-frame-pointer)
		add_link_options(-fsanitize=address,undefined)
	else()
		message(WARNING "ENABLE_SANITIZERS requested but not supported for compiler ${CMAKE_CXX_COMPILER_ID}")
	endif()
endif()

include(GNUInstallDirs)

# Header-only fem/linalg/equation layer, and the umbrella include/ path every
# compiled library also exposes publicly instead of setting it individually.
add_library(residuum_headers INTERFACE)
add_library(Residuum::residuum_headers ALIAS residuum_headers)

target_include_directories(residuum_headers INTERFACE
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(residuum_headers INTERFACE cxx_std_20)

# subdirectories
add_subdirectory(src/mesh)
add_subdirectory(src/io)
add_subdirectory(src/solver)
add_subdirectory(src/utils)

# build CLI applications
if(BUILD_APPLICATIONS)
	add_subdirectory(src/application/mesh)
	add_subdirectory(src/application/heateq)
endif()

# build tests
if(BUILD_TESTS)
	enable_testing()
	find_package(GTest REQUIRED)
	include(GoogleTest)
	add_subdirectory(tests)
endif()

# build benchmark executables
if(BUILD_BENCHMARKS)
	add_subdirectory(benchmark)
endif()

# ---- Installation ----

set(RESIDUUM_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/Residuum")

install(TARGETS residuum_headers residuum_mesh residuum_io residuum_solver residuum_expression
	EXPORT ResiduumTargets
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT ResiduumTargets
	FILE ResiduumTargets.cmake
	NAMESPACE Residuum::
	DESTINATION ${RESIDUUM_INSTALL_CMAKEDIR}
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
	${CMAKE_CURRENT_BINARY_DIR}/ResiduumConfigVersion.cmake
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/ResiduumConfig.cmake.in
	${CMAKE_CURRENT_BINARY_DIR}/ResiduumConfig.cmake
	INSTALL_DESTINATION ${RESIDUUM_INSTALL_CMAKEDIR}
)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/ResiduumConfig.cmake
	${CMAKE_CURRENT_BINARY_DIR}/ResiduumConfigVersion.cmake
	DESTINATION ${RESIDUUM_INSTALL_CMAKEDIR}
)

# Summary
message(STATUS "")
message(STATUS "========================================")
message(STATUS "  Residuum Solver Configuration Summary")
message(STATUS "========================================")
message(STATUS "  Version:      ${PROJECT_VERSION}")
message(STATUS "  Build type:   ${CMAKE_BUILD_TYPE}")
message(STATUS "  Options:")
message(STATUS "    Tests:        ${BUILD_TESTS}")
message(STATUS "    Applications: ${BUILD_APPLICATIONS}")
message(STATUS "    Benchmarks:   ${BUILD_BENCHMARKS}")
message(STATUS "========================================")
message(STATUS "")
