HZZWorkspace
============

Workspace maker for analysis involving shape information


Description
-----------

This package provides a means to build a workspace that contains 
necessary information to perform statistical interpretations on observed data, 
and tools to check the validity of inputs, calculate statistical significance, 
calculate upper limits and other things.

Several examples have been placed in directory _examples_. 
Detailed twiki page can be found [HiggsZZRunIIWorkspaces](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsZZRunIIWorkspaces).

Updated spring-summer 2019 to use CMake for compilation, packaged with ASG AnalysisBase. 


Setup with CMake
----------------

```bash
cd <workDir>
mkdir source build run

cd build
asetup AnalysisBase,21.2.68,here

cd ../source
git clone https://:@gitlab.cern.ch:8443/HZZ/HZZSoftware/HZZWorkspace.git
cd HZZWorkspace
# At the moment, the CMake version is not merged to master and lives in branch
'2-make-compatible-with-asg-cmake', so a bit more stuff is needed:
git checkout -b 2-make-compatible-with-asg-cmake origin/2-make-compatible-with-asg-cmake
cd RooFitExtensions
git submodule init
git submodule update
# This can be simplified once this is in master using
# git clone --recursive https://:@gitlab.cern.ch:8443/HZZ/HZZSoftware/HZZWorkspace.git

cd ../..
# Should now be in <workDir>/source
# Link the top level CMakeLists
ln -s HZZWorkspace/CMakeLists.topLevel.txt CMakeLists.txt

cd ../build
cmake ../source
make

cd ../run
source ../build/x86_64-slc6-gcc62-opt/setup.sh
```

