# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mbrewster/Desktop/FEM_Prime

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mbrewster/Desktop/FEM_Prime

# Include any dependencies generated for this target.
include CMakeFiles/FEM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/FEM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FEM.dir/flags.make

CMakeFiles/FEM.dir/src/matrix_utils.c.o: CMakeFiles/FEM.dir/flags.make
CMakeFiles/FEM.dir/src/matrix_utils.c.o: src/matrix_utils.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mbrewster/Desktop/FEM_Prime/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/FEM.dir/src/matrix_utils.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/FEM.dir/src/matrix_utils.c.o   -c /home/mbrewster/Desktop/FEM_Prime/src/matrix_utils.c

CMakeFiles/FEM.dir/src/matrix_utils.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/FEM.dir/src/matrix_utils.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/mbrewster/Desktop/FEM_Prime/src/matrix_utils.c > CMakeFiles/FEM.dir/src/matrix_utils.c.i

CMakeFiles/FEM.dir/src/matrix_utils.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/FEM.dir/src/matrix_utils.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/mbrewster/Desktop/FEM_Prime/src/matrix_utils.c -o CMakeFiles/FEM.dir/src/matrix_utils.c.s

CMakeFiles/FEM.dir/src/matrix_utils.c.o.requires:

.PHONY : CMakeFiles/FEM.dir/src/matrix_utils.c.o.requires

CMakeFiles/FEM.dir/src/matrix_utils.c.o.provides: CMakeFiles/FEM.dir/src/matrix_utils.c.o.requires
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/matrix_utils.c.o.provides.build
.PHONY : CMakeFiles/FEM.dir/src/matrix_utils.c.o.provides

CMakeFiles/FEM.dir/src/matrix_utils.c.o.provides.build: CMakeFiles/FEM.dir/src/matrix_utils.c.o


CMakeFiles/FEM.dir/src/fem.c.o: CMakeFiles/FEM.dir/flags.make
CMakeFiles/FEM.dir/src/fem.c.o: src/fem.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mbrewster/Desktop/FEM_Prime/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/FEM.dir/src/fem.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/FEM.dir/src/fem.c.o   -c /home/mbrewster/Desktop/FEM_Prime/src/fem.c

CMakeFiles/FEM.dir/src/fem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/FEM.dir/src/fem.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/mbrewster/Desktop/FEM_Prime/src/fem.c > CMakeFiles/FEM.dir/src/fem.c.i

CMakeFiles/FEM.dir/src/fem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/FEM.dir/src/fem.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/mbrewster/Desktop/FEM_Prime/src/fem.c -o CMakeFiles/FEM.dir/src/fem.c.s

CMakeFiles/FEM.dir/src/fem.c.o.requires:

.PHONY : CMakeFiles/FEM.dir/src/fem.c.o.requires

CMakeFiles/FEM.dir/src/fem.c.o.provides: CMakeFiles/FEM.dir/src/fem.c.o.requires
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/fem.c.o.provides.build
.PHONY : CMakeFiles/FEM.dir/src/fem.c.o.provides

CMakeFiles/FEM.dir/src/fem.c.o.provides.build: CMakeFiles/FEM.dir/src/fem.c.o


# Object files for target FEM
FEM_OBJECTS = \
"CMakeFiles/FEM.dir/src/matrix_utils.c.o" \
"CMakeFiles/FEM.dir/src/fem.c.o"

# External object files for target FEM
FEM_EXTERNAL_OBJECTS =

bin/FEM: CMakeFiles/FEM.dir/src/matrix_utils.c.o
bin/FEM: CMakeFiles/FEM.dir/src/fem.c.o
bin/FEM: CMakeFiles/FEM.dir/build.make
bin/FEM: CMakeFiles/FEM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mbrewster/Desktop/FEM_Prime/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable bin/FEM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FEM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FEM.dir/build: bin/FEM

.PHONY : CMakeFiles/FEM.dir/build

CMakeFiles/FEM.dir/requires: CMakeFiles/FEM.dir/src/matrix_utils.c.o.requires
CMakeFiles/FEM.dir/requires: CMakeFiles/FEM.dir/src/fem.c.o.requires

.PHONY : CMakeFiles/FEM.dir/requires

CMakeFiles/FEM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/FEM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/FEM.dir/clean

CMakeFiles/FEM.dir/depend:
	cd /home/mbrewster/Desktop/FEM_Prime && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mbrewster/Desktop/FEM_Prime /home/mbrewster/Desktop/FEM_Prime /home/mbrewster/Desktop/FEM_Prime /home/mbrewster/Desktop/FEM_Prime /home/mbrewster/Desktop/FEM_Prime/CMakeFiles/FEM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FEM.dir/depend

