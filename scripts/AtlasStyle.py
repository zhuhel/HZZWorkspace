import ROOT
import os

SCRIPT_DIR = os.getenv("HZZWSCODEDIR")+"/scripts"
if not hasattr(ROOT, "SetAtlasStyle"):
	ROOT.gROOT.LoadMacro(SCRIPT_DIR+"/AtlasStyle.C")
ROOT.SetAtlasStyle()
