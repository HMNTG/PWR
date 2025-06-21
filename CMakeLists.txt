cmake_minimum_required(VERSION 3.10)
project(MPI_C_Project C)

# MPI suchen
find_package(MPI REQUIRED)

# Spracheinstellungen
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

include_directories(mpi_program "src/")

# src files
file(GLOB_RECURSE all_SRC
	"src/*.c"
)
set(src)
list(APPEND src
	${all_SRC}
	main.c
)

# Executable bauen
add_executable(mpi_program ${src})

# MPI einbinden
target_include_directories(mpi_program PUBLIC ${MPI_C_INCLUDE_PATH})
target_link_libraries(mpi_program ${MPI_C_LIBRARIES})
