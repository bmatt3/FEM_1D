# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/mbrewster/Desktop/FEM_Prime/CMakeFiles /home/mbrewster/Desktop/FEM_Prime/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/mbrewster/Desktop/FEM_Prime/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named FEM

# Build rule for target.
FEM: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 FEM
.PHONY : FEM

# fast build rule for target.
FEM/fast:
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/build
.PHONY : FEM/fast

src/fem.o: src/fem.c.o

.PHONY : src/fem.o

# target to build an object file
src/fem.c.o:
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/fem.c.o
.PHONY : src/fem.c.o

src/fem.i: src/fem.c.i

.PHONY : src/fem.i

# target to preprocess a source file
src/fem.c.i:
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/fem.c.i
.PHONY : src/fem.c.i

src/fem.s: src/fem.c.s

.PHONY : src/fem.s

# target to generate assembly for a file
src/fem.c.s:
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/fem.c.s
.PHONY : src/fem.c.s

src/matrix_utils.o: src/matrix_utils.c.o

.PHONY : src/matrix_utils.o

# target to build an object file
src/matrix_utils.c.o:
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/matrix_utils.c.o
.PHONY : src/matrix_utils.c.o

src/matrix_utils.i: src/matrix_utils.c.i

.PHONY : src/matrix_utils.i

# target to preprocess a source file
src/matrix_utils.c.i:
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/matrix_utils.c.i
.PHONY : src/matrix_utils.c.i

src/matrix_utils.s: src/matrix_utils.c.s

.PHONY : src/matrix_utils.s

# target to generate assembly for a file
src/matrix_utils.c.s:
	$(MAKE) -f CMakeFiles/FEM.dir/build.make CMakeFiles/FEM.dir/src/matrix_utils.c.s
.PHONY : src/matrix_utils.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... FEM"
	@echo "... edit_cache"
	@echo "... src/fem.o"
	@echo "... src/fem.i"
	@echo "... src/fem.s"
	@echo "... src/matrix_utils.o"
	@echo "... src/matrix_utils.i"
	@echo "... src/matrix_utils.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

