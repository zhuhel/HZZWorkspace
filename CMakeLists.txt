## automatically generated CMakeLists.txt file

# Declare the package
atlas_subdir( HZZWorkspace )

# Declare ASG dependencies 
atlas_depends_on_subdirs( 
  PRIVATE 
  PUBLIC 
  RooFitExtensions 
  ${extra_dep} ) 

# Declare external dependencies ... default here is to include ROOT
find_package( ROOT COMPONENTS MathCore RIO Core Tree Hist 
			      RooFitCore RooFit RooStats HistFactory
			      Graf Graf2d Gpad 
			      Minuit MathMore Matrix REQUIRE )

# Declare dictionary 
atlas_add_root_dictionary( HZZWorkspaceLib HZZWorkspaceLibCintDict 
  ROOT_HEADERS HZZWorkspace/*.h Root/LinkDef.h
  EXTERNAL_PACKAGES ROOT ) 

# Declare package as a library
# Note the convention that library names get "Lib" suffix
# Any package you depend on you should add
# to LINK_LIBRARIES line below (see the example)
atlas_add_library( HZZWorkspaceLib Root/*.cxx HZZWorkspace/*.h ${HZZWorkspaceLibCintDict}
                   PUBLIC_HEADERS HZZWorkspace
                   PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                   PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES}
                   LINK_LIBRARIES RooFitExtensions
)


# if you add an application (exe) to this package
# declare it like this (note convention that apps go in the util dir)
# atlas_add_executable( MyApp util/myApp.cxx
#                       LINK_LIBRARIES HZZWorkspaceLib
# )
file(GLOB_RECURSE files "utils/*.cxx")
foreach(_exeFile ${files})
    get_filename_component(_theExec ${_exeFile} NAME_WE)
    get_filename_component(_theLoc ${_exeFile} DIRECTORY)
    # we specify a folder for programs we do not want to compile. Useful during r21 transition...
    atlas_add_executable( ${_theExec}
    ${_exeFile}
    INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} 
    LINK_LIBRARIES ${ROOT_LIBRARIES} HZZWorkspaceLib RooFitExtensions )
endforeach()


# Install python modules, joboptions, and share content
#atlas_install_python_modules( python/*.py )
#atlas_install_joboptions( share/*.py )
#atlas_install_data( data/* )
# You can access your data from code using path resolver, e.g.
# PathResolverFindCalibFile("HZZWorkspace/file.txt")



