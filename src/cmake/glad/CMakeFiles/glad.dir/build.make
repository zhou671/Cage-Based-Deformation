# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake

# Include any dependencies generated for this target.
include glad/CMakeFiles/glad.dir/depend.make

# Include the progress variables for this target.
include glad/CMakeFiles/glad.dir/progress.make

# Include the compile flags for this target's objects.
include glad/CMakeFiles/glad.dir/flags.make

glad/CMakeFiles/glad.dir/src/glad.c.o: glad/CMakeFiles/glad.dir/flags.make
glad/CMakeFiles/glad.dir/src/glad.c.o: /home/yzhou162/preceputal-smi/libigl/external/glad/src/glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object glad/CMakeFiles/glad.dir/src/glad.c.o"
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/glad.dir/src/glad.c.o   -c /home/yzhou162/preceputal-smi/libigl/external/glad/src/glad.c

glad/CMakeFiles/glad.dir/src/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/glad.dir/src/glad.c.i"
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/yzhou162/preceputal-smi/libigl/external/glad/src/glad.c > CMakeFiles/glad.dir/src/glad.c.i

glad/CMakeFiles/glad.dir/src/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/glad.dir/src/glad.c.s"
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad && /usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/yzhou162/preceputal-smi/libigl/external/glad/src/glad.c -o CMakeFiles/glad.dir/src/glad.c.s

# Object files for target glad
glad_OBJECTS = \
"CMakeFiles/glad.dir/src/glad.c.o"

# External object files for target glad
glad_EXTERNAL_OBJECTS =

glad/libglad.a: glad/CMakeFiles/glad.dir/src/glad.c.o
glad/libglad.a: glad/CMakeFiles/glad.dir/build.make
glad/libglad.a: glad/CMakeFiles/glad.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libglad.a"
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad && $(CMAKE_COMMAND) -P CMakeFiles/glad.dir/cmake_clean_target.cmake
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/glad.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
glad/CMakeFiles/glad.dir/build: glad/libglad.a

.PHONY : glad/CMakeFiles/glad.dir/build

glad/CMakeFiles/glad.dir/clean:
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad && $(CMAKE_COMMAND) -P CMakeFiles/glad.dir/cmake_clean.cmake
.PHONY : glad/CMakeFiles/glad.dir/clean

glad/CMakeFiles/glad.dir/depend:
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src /home/yzhou162/preceputal-smi/libigl/external/glad /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/glad/CMakeFiles/glad.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : glad/CMakeFiles/glad.dir/depend

