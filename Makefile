# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Mon Mar 28, 2022 at 01:02 AM -0400

BINPATH := bin
VPATH := include:src
CPP_FILES	:=	$(wildcard src/*.cpp)
EXE_FILES	:=	$(patsubst src/%.cpp,$(BINPATH)/%.exe,$(CPP_FILES))

TIME_STAMP	:=	$(shell date +"%y_%m_%d_%H_%M")

COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags)
LINKFLAGS	:=	$(shell root-config --libs)
ADDCXXFLAGS	:=	-O2 -march=native -mtune=native
ADDLINKFLAGS	:=	-lyaml-cpp -lRooFitCore -lRooFit -lRooStats -lRooUnfold


###########
# General #
###########
.PHONY: exe clean

exe: $(EXE_FILES)

clean:
	@rm -rf ./gen/*


#######
# RDX #
#######
.PHONY: build-tagged-histo build-rdx-true-to-tag-2016

build-rdx-tag-2016:
	$(eval OUT_DIR	:=	./gen/rdx-$(TIME_STAMP)-tag-2016)
	./scripts/build_histo_tagged.py -c ./spec/rdx-2016.yml -o $(OUT_DIR)

build-rdx-true-to-tag-2016:
	$(eval OUT_DIR	:=	./gen/rdx-$(TIME_STAMP)-true_to_tag-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-2016.yml -o $(OUT_DIR) -y 16

build-rdx-merged-2016:
	$(eval OUT_DIR	:=	./gen/rdx-$(TIME_STAMP)-merged-2016)
	./scripts/merge_histo.py -c ./spec/rdx-2016.yml -o $(OUT_DIR)


########
# Test #
########
.PHONY: test-pidcalib2-wrapper test-unfold

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-2016.yml -o ./gen --dry-run

test-unfold: $(BINPATH)/UnfoldMisID.exe
	./bin/UnfoldMisID.exe -c ./spec/rdx-2016.yml --dryRun


###############
# Compile C++ #
###############

$(BINPATH)/%.exe: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) $(ADDLINKFLAGS)
