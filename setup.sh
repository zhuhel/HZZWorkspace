#!/bin/bash
script_name=${BASH_SOURCE[0]}
if [[ $SHELL == *"zsh"* ]]; then
    script_name=${(%):-%N}
fi
currentDir=$PWD

##setup gcc and python
/bin/grep ' release [2345]\.' /etc/redhat-release >/dev/null 2>&1 && \
  echo "WARNING: This version of the HSG7 ROOT build will not work on an SLC5 machine. Please use a different machine (eg. lxplus.cern.ch)" >&2

# GCC 4.9.3
PATH="/afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin:$PATH"
LD_LIBRARY_PATH="/afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/lib64:$LD_LIBRARY_PATH"

# Python 2.7.4
#PYTHONDIR="/afs/cern.ch/sw/lcg/external/Python/2.7.4/x86_64-slc6-gcc48-opt"
#PATH="$PYTHONDIR/bin:$PATH"
#LD_LIBRARY_PATH="$PYTHONDIR/lib:$LD_LIBRARY_PATH"
#PYTHONPATH=$PYTHONPATH:"$PYTHONDIR/lib/python2.7"
#export PATH LD_LIBRARY_PATH PYTHONDIR PYTHONPATH
# additional lib.
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/user/s/shsun/public/python/numpy/1.10.4/install/lib/python2.7/site-packages/:/afs/cern.ch/user/s/shsun/public/python/scipy/install/lib/python2.7/site-packages/:/afs/cern.ch/user/x/xju/public/python/lib/python2.7/site-packages/

#ROOT
#RootDir=/afs/cern.ch/atlas/project/HSG7/root/current/x86_64-slc6-gcc48/
#RootDir=/afs/cern.ch/atlas/project/HSG7/root/root_v6-04-02/x86_64-slc6-gcc49/
#cd $RootDir
#source bin/thisroot.sh
source $ATLAS_LOCAL_ROOT_BASE/user/atlasLocalSetup.sh
lsetup "root 6.08.06-HiggsComb-x86_64-slc6-gcc49-opt"

# increase stack size - needed for large workspaces
ulimit -S -s unlimited

cd $currentDir
#set ws code env
export HZZWSCODEDIR="$( cd "$( dirname "$script_name" )" && pwd )"
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HZZWSCODEDIR}/lib:/afs/cern.ch/work/k/kecker/public/HZZTensorWS/RooLagrangianMorphFunc/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HZZWSCODEDIR}/lib
export PATH=$PATH:${HZZWSCODEDIR}/bin:${HZZWSCODEDIR}/test-bin
export HZZWSDIR="/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/Workspaces/"

if [ ! -f ./Hzzws_Dict_rdict.pcm ]; then
    ln -s $HZZWSCODEDIR/Root/Hzzws_Dict_rdict.pcm .
fi

D1=$PWD

# set up environment for running Python code
export PATH=$HZZWSCODEDIR/Python/bin:$PATH
export PYTHONPATH=$HZZWSCODEDIR/Python/modules:$PYTHONPATH
