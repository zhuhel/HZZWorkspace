# Python code for HZZ workspaces

*Developed by Stefan von Buddenbrock (stef.von.b@cern.ch)*

Contents:
1. Motivation
2. Setting up the code
3. Structure of the code
4. Modules
5. Programs
6. Examples

## Motivation

Many of the programs and scripts in the HZZWorkspace package are currently written in C++.
They work well and achieve their goals, but tend to be quite long and complicated.
This is because working with systematics often has to do with iterating over big lists of strings.
It is therefore possible to improve the readibility of the code using Python, since Python's features for iterators and string handling make it ideal.
In addition to this, Python is often a better choice for making plots and outputs (using ROOT or otherwise), and so the visualisation of systematic effects could be improved.

## Setting up the code

The Python programs are put into the `Python/` directory in the main HZZWorkspace package.
The most up to date version will always be in a Git branch, currently called `Python_Development`.
The environment is set up in the same way as always, by doing `source setup.sh` in the top directory.
This will add the Python programs to the `$PATH` variable and the appropriate Python modules to the `$PYTHONPATH` varaible.

## Structure of the code

In the `Python/` directory, the following can be found:
* `bin/` -- The directory which contains programs intended to run with an input file and produce an output.
* `modules` -- All of the helper modules and common tools that can be used by the programs (to parse config files, make plots and tables, etc.)
* `configuration/` -- Other files needed to set up the features used by the Python modules and programs.

## Modules

Below is a brief description of the different modules in the code.

#### Utilities.py

Currently consists of functions used to parse config files and return them as Python dictionaries.
There are several functions which can check the input config files and return information about any faults.
In addition, some functions exist for conditioning outputs to be plotted more easily.

#### Plotter.py

Uses ROOT to make the various plots needed for visualising results.
Note that a function in this module should necessarily write an output PDF file (instead of returning a TCanvas, for instance), since ROOT usually deletes objects outside of the function scope.
For each function, a `publicity` argument can be included, containing the text which goes after the **_ATLAS_** logo on the plots (default is "Simulation Internal").

#### Writer.py

Contains code for making tables of results (preferrably in LaTeX).
These sorts of functions are much more focused towards producing verbose outputs for inspection and debugging.

## Programs

Below is the description of each program found in the `bin/` directory.

#### SysProd.py

This program mimics the primary function of the SysProd.cxx file in the C++ version of the code.
It takes an input configuration file, containing a list of samples to be run over and categorisations.
The steps of the calculation are written to the screen, and debug information will be written to `log.txt`.
The config file should also point to a list of nuisance parameters over which the scan is made.
Subject to the specified cuts in the config file, SysProd.py will integrate the M4l distribution and return a systematic variation of the integral compared to the nominal distribution.
The result will be output into a text file, and some simple plots will be made for each category.
Optionally, all of the numbers can be written to the log file in JSON format.

**Usage**: `SysProd.py config_file [-o output_dir] [--json]`

## Examples

#### Comparing nuisance parameters for two different versions of MiniTree production

Suppose you want to identify the difference in the nuisance parameter effects in two different versions of MiniTrees.
This is how I would do it.

1. Assuming you have the code checked out from GitLab (remember, this might require checking out a specific branch), set up the environment by doing `source setup.sh`.
2. Set up two different directories for each version of the MiniTrees, and put the appropriate config files in each directory. Note that you might want to use the same list of nuisance parameters in each run, so that you can compare them all.
3. For each config file in each directory, run the SysProd.py code as shown above. You should get output text files called `norm_[sample_name].txt` for each sample included in the config file.
4. Now that you have two separate output files, you can use the Plotter.py module to make a comparison plot. This can be done interactively using a Python prompt, or in a Python script. Try the following:

```
python
from Plotter import plot_compare_NPs
import ROOT

ROOT.gROOT.SetBatch(True) # This is just to prevent canvases from popping up

plot_compare_NPs(["versionA/norm_sample.txt", "versionB/norm_sample.txt"],
                 ["Version_A", "Version_B"], "Sample")
```

Congratulations! This will make a plot (for each category) comparing version A to version B.
The first argument of `plot_compare_NPs()` contains a list of the different text files to input for the plotting.
The second argument gives each input file a readable name, in the same order as the first argument.
This is for the purpose of the legend, so that you know which plot is which.
The third argument is an identifier, so that you can remember on the plot what you were trying to demonstrate.
Usually it will be what is common between the two versions (in this case, the sample was the same, so we can remind ourselves of that while the legend will tell us what is different).

Note: there are other arguments you can explore in the code to help you get exactly what you want. But this is a good starting point.
