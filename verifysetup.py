import platform, sys
import numpy
import ROOT

if platform.system() == "Darwin":
    print "macOS", platform.mac_ver()[0]
print "OS:", platform.platform()
print "python:", sys.version
print "numpy: ", numpy.__version__
print "ROOT:", ROOT.gROOT.GetVersion()
