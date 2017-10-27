import ROOT
import os
ROOT.gROOT.LoadMacro(os.getenv("HZZWSCODEDIR")+"/scripts/AtlasStyle.C") 
ROOT.SetAtlasStyle()
