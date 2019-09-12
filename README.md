HZZWorkspace
============

Workspace maker for analysis involving shape information


Description
-----------

This package provides a means to build a workspace that contains
necessary information to perform statistical interpretations on observed data,
and tools to check the validity of inputs, calculate statistical significance,
calculate upper limits and other things.

Several examples have been placed in directory `examples`.
<!--Detailed twiki page can be found [HiggsZZRunIIWorkspaces](https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsZZRunIIWorkspaces).-->

Updated spring-summer 2019 to use CMake for compilation, packaged with ASG AnalysisBase.
WikiPages are available [here](https://gitlab.cern.ch/HZZ/HZZSoftware/HZZWorkspace/wikis/home)
and a list of related presentations is available [here](https://gitlab.cern.ch/HZZ/HZZSoftware/HZZWorkspace/wikis/Resources/Related-Talks).

Setup with CMake
----------------

```bash
cd <workDir>
mkdir source build run

cd build
asetup AnalysisBase,{latest release here},here
# Confirmed numpy import working with 21.2.89 (does not work in 21.2.68)

cd ../source
git clone https://:@gitlab.cern.ch:8443/HZZ/HZZSoftware/HZZWorkspace.git
cd HZZWorkspace
# CMake version is default on master
git clone --recursive https://:@gitlab.cern.ch:8443/HZZ/HZZSoftware/HZZWorkspace.git


# To check out an alternate branch directly, with submodules, a bit more stuff is needed:
git checkout -b {branchname} origin/{branchname}
cd RooFitExtensions
git submodule init
git submodule update
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
