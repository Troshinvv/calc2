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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Users/Misha/CLionProjects/calc2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug

# Utility rule file for G__Calculator.

# Include the progress variables for this target.
include CMakeFiles/G__Calculator.dir/progress.make

CMakeFiles/G__Calculator: G__Calculator.cxx
CMakeFiles/G__Calculator: libCalculator_rdict.pcm
CMakeFiles/G__Calculator: libCalculator.rootmap


G__Calculator.cxx: ../src/CalculatorLinkDef.h
G__Calculator.cxx: ../src/correlation.h
G__Calculator.cxx: ../src/functions.h
G__Calculator.cxx: ../src/correlation.h
G__Calculator.cxx: ../src/functions.h
G__Calculator.cxx: ../src/CalculatorLinkDef.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating G__Calculator.cxx, libCalculator_rdict.pcm, libCalculator.rootmap"
	/usr/bin/cmake -E env LD_LIBRARY_PATH=/home/mikhail/root-6.24.06/install/lib: /home/mikhail/root-6.24.06/install/bin/rootcling -v2 -f G__Calculator.cxx -s /mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug/libCalculator.so -rml libCalculator.so -rmf /mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug/libCalculator.rootmap -I/mnt/c/Users/Misha/CLionProjects/calc2 -I/home/mikhail/boost_1_80_0 -I/home/mikhail/root-6.24.06/install/include -I/home/mikhail/QnTools/install/lib/cmake/QnTools/../../../include/QnTools src/correlation.h src/functions.h /mnt/c/Users/Misha/CLionProjects/calc2/src/CalculatorLinkDef.h

libCalculator_rdict.pcm: G__Calculator.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libCalculator_rdict.pcm

libCalculator.rootmap: G__Calculator.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libCalculator.rootmap

G__Calculator: CMakeFiles/G__Calculator
G__Calculator: G__Calculator.cxx
G__Calculator: libCalculator_rdict.pcm
G__Calculator: libCalculator.rootmap
G__Calculator: CMakeFiles/G__Calculator.dir/build.make

.PHONY : G__Calculator

# Rule to build all files generated by this target.
CMakeFiles/G__Calculator.dir/build: G__Calculator

.PHONY : CMakeFiles/G__Calculator.dir/build

CMakeFiles/G__Calculator.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/G__Calculator.dir/cmake_clean.cmake
.PHONY : CMakeFiles/G__Calculator.dir/clean

CMakeFiles/G__Calculator.dir/depend:
	cd /mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/Misha/CLionProjects/calc2 /mnt/c/Users/Misha/CLionProjects/calc2 /mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug /mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug /mnt/c/Users/Misha/CLionProjects/calc2/cmake-build-debug/CMakeFiles/G__Calculator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/G__Calculator.dir/depend

