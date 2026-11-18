cmake_minimum_required(VERSION 3.18)
project(ResiduumConsumerDemo CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# A standalone CMake project
find_package(Residuum REQUIRED)

add_executable(consumer_demo main.cpp)

target_link_libraries(consumer_demo
	PRIVATE
		Residuum::residuum_mesh
		Residuum::residuum_expression
)
