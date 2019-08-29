//
// ATLAS Style, based on a style file from BaBar
//
#ifndef __atlas_styple_c__
#define __atlas_styple_c__

#include <iostream>

// #include "AtlasStyle.h"

#include "TROOT.h"

/// Put these methods, which used to be in a .c macro,
/// into a class as static methods, to make them accessible 
/// from within pyROOT. 
class AtlasStyleHelper{
public:
    static TStyle* AtlasStyle();
    static void SetAtlasStyle ();
};

#endif
