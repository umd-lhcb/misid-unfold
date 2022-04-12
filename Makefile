# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Tue Apr 12, 2022 at 02:27 PM -0400

BINPATH := ./bin
GENPATH := ./gen
VPATH := include:src:docs
CPP_FILES	:=	$(wildcard src/*.cpp)
EXE_FILES	:=	$(patsubst src/%.cpp,$(BINPATH)/%.exe,$(CPP_FILES))
TEX_FILES	:=	$(wildcard docs/*.tex)
PDF_FILES	:=	$(patsubst docs/%.tex,$(GENPATH)/%.pdf,$(TEX_FILES))

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
	@rm -rf $(GENPATH)/*


#######
# RDX #
#######
.PHONY: build-tagged-histo build-rdx-true-to-tag-2016 plot-rdx-2016

build-rdx-tag-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-tag-2016)
	./scripts/build_histo_tagged.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016

build-rdx-true-to-tag-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016

build-rdx-merged-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-merged-2016)
	./scripts/merge_histo.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016

build-rdx-unfolded-2016: $(BINPATH)/UnfoldMisID.exe
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-unfolded-2016)
	@mkdir -p $(OUT_DIR)
	$< --debug --iteration 20 \
		-e ./histos/rdx-22_04_12_14_03-merged-2016/merged.root \
		-y ./histos/rdx-22_04_12_14_03-tag-2016/tagged.root \
		-o $(OUT_DIR) \
		-c ./spec/rdx-run2.yml | tee $(OUT_DIR)/stdout.log

plot-rdx-2016:
	./scripts/plot_histo.py -o $(GENPATH) \
		-i ./histos/rdx-22_03_30_12_40-unfolded-2016/unfolded.root
	./scripts/plot_histo.py -o $(GENPATH) -s Tag \
		-i ./histos/rdx-22_03_27_18_05-tag-2016/tagged.root


########
# Test #
########
.PHONY: test-pidcalib2-wrapper test-unfold

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2.yml -o $(GENPATH) --dry-run

test-unfold: $(BINPATH)/UnfoldMisID.exe
	./bin/UnfoldMisID.exe -c ./spec/rdx-run2.yml --dryRun


###############
# Compile C++ #
###############

$(BINPATH)/%.exe: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) $(ADDLINKFLAGS)


#######
# Doc #
#######
.PHONY: doc

doc: $(PDF_FILES)

$(GENPATH)/%.pdf: %.tex
	latexmk $<
