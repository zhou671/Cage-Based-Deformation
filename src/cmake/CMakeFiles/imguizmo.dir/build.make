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
include CMakeFiles/imguizmo.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/imguizmo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/imguizmo.dir/flags.make

CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.o: CMakeFiles/imguizmo.dir/flags.make
CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.o: /home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.o -c /home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp

CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp > CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.i

CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp -o CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.s

# Object files for target imguizmo
imguizmo_OBJECTS = \
"CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.o"

# External object files for target imguizmo
imguizmo_EXTERNAL_OBJECTS =

libimguizmo.a: CMakeFiles/imguizmo.dir/home/yzhou162/preceputal-smi/libigl/external/imguizmo/ImGuizmo.cpp.o
libimguizmo.a: CMakeFiles/imguizmo.dir/build.make
libimguizmo.a: CMakeFiles/imguizmo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libimguizmo.a"
	$(CMAKE_COMMAND) -P CMakeFiles/imguizmo.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/imguizmo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/imguizmo.dir/build: libimguizmo.a

.PHONY : CMakeFiles/imguizmo.dir/build

CMakeFiles/imguizmo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/imguizmo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/imguizmo.dir/clean

CMakeFiles/imguizmo.dir/depend:
	cd /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake /home/yzhou162/preceputal-smi/Cage-Based-Deformation-MVC/src/cmake/CMakeFiles/imguizmo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/imguizmo.dir/depend

