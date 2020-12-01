cmake_minimum_required(VERSION 3.5)

project(FC VERSION 0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

### OPTIONS

option(USE_DEBUG "Enter debug mode" OFF)

if (USE_DEBUG)
  add_definitions(-DDEBUG)
endif()

### LIBRARIES



### EXECUTABLE 

add_executable(reject src/reject.cpp)
add_executable(lnv    src/lnv.cpp)

target_include_directories(reject PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(lnv    PUBLIC ${PROJECT_SOURCE_DIR}/include)


### INSTALL

# Binaries
install (TARGETS reject lnv
	DESTINATION ${CMAKE_BINARY_DIR}/bin)

# Library
#install (TARGETS library_target
#LIBRARY DESTINATION lib)

# Header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
	DESTINATION ${CMAKE_BINARY_DIR}/include)

# Config
#install (FILES cmake-examples.conf
#DESTINATION etc)
