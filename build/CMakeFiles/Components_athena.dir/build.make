# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/software/linux-debian-x86_64/gcc-11.2.0/cmake-3.21.4-zkiq5ma7qlm65gws4hyn4vnfml2ssskd/bin/cmake

# The command to remove a file.
RM = /opt/software/linux-debian-x86_64/gcc-11.2.0/cmake-3.21.4-zkiq5ma7qlm65gws4hyn4vnfml2ssskd/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ryan/eic/eic_endcap_insert

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ryan/eic/eic_endcap_insert/build

# Utility rule file for Components_athena.

# Include any custom commands dependencies for this target.
include CMakeFiles/Components_athena.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Components_athena.dir/progress.make

CMakeFiles/Components_athena: libathena.components

libathena.components: lib/libathena.so
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/ryan/eic/eic_endcap_insert/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating libathena.components"
	cd /home/ryan/eic/eic_endcap_insert/build/lib && LD_LIBRARY_PATH=/home/ryan/eic/eic_endcap_insert/build/lib:/usr/local/lib:/home/ryan/eic/local/lib:/lib/x86_64-linux-gnu:/usr/local/lib:/usr/local/lib64:/opt/detector/athena-nightly/lib:/home/ryan/eic/local/lib /usr/local/bin/listcomponents_dd4hep -o libathena.components libathena.so

Components_athena: CMakeFiles/Components_athena
Components_athena: libathena.components
Components_athena: CMakeFiles/Components_athena.dir/build.make
.PHONY : Components_athena

# Rule to build all files generated by this target.
CMakeFiles/Components_athena.dir/build: Components_athena
.PHONY : CMakeFiles/Components_athena.dir/build

CMakeFiles/Components_athena.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Components_athena.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Components_athena.dir/clean

CMakeFiles/Components_athena.dir/depend:
	cd /home/ryan/eic/eic_endcap_insert/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ryan/eic/eic_endcap_insert /home/ryan/eic/eic_endcap_insert /home/ryan/eic/eic_endcap_insert/build /home/ryan/eic/eic_endcap_insert/build /home/ryan/eic/eic_endcap_insert/build/CMakeFiles/Components_athena.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Components_athena.dir/depend
