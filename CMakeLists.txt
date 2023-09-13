cmake_minimum_required(VERSION 3.10.2)
project(DITHERNUMBERSPLAYGROUND)
set(CMAKE_CXX_STANDARD 17)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/imgui ${CMAKE_CURRENT_SOURCE_DIR}/system ${CMAKE_CURRENT_SOURCE_DIR}/system/utilities)
add_executable(${PROJECT_NAME}
		main.cpp
		string_math1.cpp
		string_math1.hpp
		dithernumbers.cpp
		dithernumbers.hpp
		)
target_link_libraries(${PROJECT_NAME}
)