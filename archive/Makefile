
MAKEFLAGS = --no-print-directory -r -s -j2

#
# Include the architecture definitions from the ROOT sources
#
#  Makefile.arch can be in two different locations depending on the system
#  you're compiling on. The Fink installed version of ROOT has this file
#  in a different location than the "normally installed" ROOT versions...
#
ARCH_LOC_1 := $(wildcard $(shell root-config --prefix)/etc/Makefile.arch)
# ARCH_LOC_2 := $(wildcard $(shell root-config --prefix)/share/root/test/Makefile.arch)
ARCH_LOC_2 := $(wildcard $(shell root-config --prefix)/test/Makefile.arch)
ARCH_LOC_3 := $(wildcard $(shell root-config --prefix)/share/doc/root/test/Makefile.arch)
ifneq ($(strip $(ARCH_LOC_1)),)
  $(info Using $(ARCH_LOC_1))
  include $(ARCH_LOC_1)
else
  ifneq ($(strip $(ARCH_LOC_2)),)
    $(info Using $(ARCH_LOC_2))
    include $(ARCH_LOC_2)
  else
    ifneq ($(strip $(ARCH_LOC_3)),)
      $(info Using $(ARCH_LOC_3))
      include $(ARCH_LOC_3)
    else
      $(error Could not find Makefile.arch!)
    endif
  endif
endif

CXX           = g++
ROOTCXXFLAGS  = $(shell root-config --cflags)

CXXFLAGS      = -std=c++11 -fPIC -g -O0 -Wall $(ROOTCXXFLAGS)

ROOTLIBS   = $(shell root-config --libs) -lRooFitCore -lRooFit -lRooStats -lHistFactory -lMinuit -lMathMore

LIBFLAGS  =  -O2 -shared -m64 $(ROOTLIBS) -rdynamic -dynamiclib -L/usr/local/Cellar/boost/1.55.0_2/lib
# INCLUDES = -I. -I/usr/local/Cellar/boost/1.55.0_2/include -I/afs/cern.ch/work/k/kecker/public/HZZTensorWS/RooLagrangianMorphFunc/inc $(CXXFLAGS)
INCLUDES = -I. -I/usr/local/Cellar/boost/1.55.0_2/include $(CXXFLAGS)


# Set the locations of some files
SRCDIR = Root
LIBRARY = Hzzws
OBJDIR = obj
INCDIR = Hzzws
DEPDIR  = $(OBJDIR)/dep
LIB_PATH = lib

# BINFLAGS = -L./$(LIB_PATH) -lHzzws -L/afs/cern.ch/work/k/kecker/public/HZZTensorWS/RooLagrangianMorphFunc/lib -lRooLagrangianMorphFuncMy $(ROOTLIBS) -m64 -rdynamic
BINFLAGS = -L./$(LIB_PATH) -lHzzws $(ROOTLIBS) -m64 -rdynamic

DICTHEAD  = $(SRCDIR)/$(LIBRARY)_Dict.h
DICTFILE  = $(SRCDIR)/$(LIBRARY)_Dict.$(SrcSuf)
DICTOBJ   = $(OBJDIR)/$(LIBRARY)_Dict.$(ObjSuf)
DICTLDEF  = $(INCDIR)/$(LIBRARY)_LinkDef.h
SKIPCPPLIST = $(DICTFILE)
SKIPHLIST   = $(DICTHEAD) $(DICTLDEF)
LIBFILE   = $(LIB_PATH)/lib$(LIBRARY).a
SHLIBFILE = $(LIB_PATH)/lib$(LIBRARY).$(DllSuf)

# List of all header and source files to build
HLIST= $(shell grep -l ClassDef $(INCDIR)/*h)
# Add a file
HLISTLAST = $(HLIST) Hzzws/RooStatsHelper.h
CPPLIST = $(filter-out $(SKIPCPPLIST),$(wildcard $(SRCDIR)/*.$(SrcSuf)))

### List of all object files ###
OBJECTSORG=$(patsubst %.$(SrcSuf),%.o,$(CPPLIST))
OBJECTS=$(subst Root,obj,$(OBJECTSORG))

### executes ###
TARGET1=$(wildcard utils/*.cxx)
TARGET2=$(subst utils,bin,$(TARGET1))
TARGET_BIN = $(patsubst %.cxx,%,$(TARGET2))

### programs for test ###
TEST_TARGET1 = $(wildcard test/*.cxx)
TEST_TARGET2 = $(subst test/,test-bin/,$(TEST_TARGET1))
TARGET_TESTBIN = $(patsubst %.cxx,%,$(TEST_TARGET2))

###################################################################################
SILENT=
# all: build ./lib/libHzzws.so $(TARGET_BIN) $(TARGET_TESTBIN)
default: build shlib  shlib $(TARGET_BIN) $(TARGET_TESTBIN)

build:
	@mkdir -p obj
	@mkdir -p lib
	@mkdir -p bin
	@mkdir -p test-bin

$(DICTFILE): $(HLISTLAST) $(DICTLDEF)
	@echo "Generating dictionary $@ $(HLISTLAST)"
	#@$(shell root-config --exec-prefix)/bin/rootcint -f $(DICTFILE) -c -p $(INCLUDES) $^
	@$(shell root-config --exec-prefix)/bin/rootcling -f $(DICTFILE) -c -p $(INCLUDES) $^

# Rule to comile the dictionary
$(DICTOBJ): $(DICTFILE)
	@echo "Compiling dictionary $<"
	@mkdir -p $(OBJDIR)
	$(CXX) $(C++FLAGS) -O2 -c $(INCLUDES) -o $@ $<

#

./test-bin/% :  obj/%.o
	$(SILENT)echo Linking `basename $@`
	@echo $(CXX)  $< $(BINFLAGS) -o $@
	$(SILENT)$(CXX)  $< $(BINFLAGS) -o $@

./bin/% : ./obj/%.o
	$(SILENT)echo Linking `basename $@`
	@echo $(CXX) -o $@ $< $(BINFLAGS)
	$(SILENT)$(CXX) -o $@ $< $(BINFLAGS)

./obj/%.o : ./Root/%.cxx
	echo "Compiling $<"
	$(CXX) $(INCLUDES) -g -c $<  -o $@

./obj/%.o : ./test/%.cxx
	echo "Compiling $<"
	$(CXX) $(INCLUDES) -g -c $<  -o $@

./obj/%.o : ./utils/%.cxx
	echo "Compiling $<"
	$(CXX) $(INCLUDES) -g -c $<  -o $@


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include $(foreach var,$(notdir $(CPPLIST:.$(SrcSuf)=.d)),$(DEPDIR)/$(var))
endif
endif

$(DEPDIR)/%.d: %.$(SrcSuf)
	@mkdir -p $(DEPDIR)
	if test -f $< ; then \
		echo "Making $(@F)"; \
		$(SHELL) -ec '$(CPP) -MM $(C++FLAGS) $(INCLUDES) $< | sed '\''/Cstd\/rw/d'\'' > $@'; \
	fi

# Rule to combine objects into a unix shared library
$(SHLIBFILE): $(OBJECTS) $(DICTOBJ)
	@echo "Making shared library: $(SHLIBFILE)"
	@rm -f $(SHLIBFILE)
ifneq (,$(findstring macosx,$(ARCH)))
	@$(LD) $(LDFLAGS)  -dynamiclib -single_module -undefined dynamic_lookup -O2 $(addprefix $(OBJDIR)/,$(OBJECTS)) $(DICTOBJ) -o $(SHLIBFILE)
else
	@echo $(LD) $(LDFLAGS) $(SOFLAGS)  -O2 $(OBJECTS) $(DICTOBJ) -o $(SHLIBFILE)
	@$(LD) $(LDFLAGS) $(SOFLAGS)  -O2 $(OBJECTS) $(DICTOBJ) -o $(SHLIBFILE)
endif

# Useful build targets
shlib: $(SHLIBFILE)

clean:
	@rm -rf ./obj
	@rm -rf ./lib
	@rm -rf ./bin
	@rm -rf ./test-bin
