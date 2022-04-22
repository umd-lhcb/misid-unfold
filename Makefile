# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Fri Apr 22, 2022 at 04:46 AM -0400

BINPATH := ./bin
GENPATH := ./gen
VPATH := include:src:docs
CPP_FILES	:=	$(wildcard src/*.cpp)
EXE_FILES	:=	$(patsubst src/%.cpp,$(BINPATH)/%,$(CPP_FILES))
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
	@mkdir -p $(OUT_DIR)
	./scripts/merge_histo.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016 | tee $(OUT_DIR)/stdout.log

build-rdx-unfolded-2016: $(BINPATH)/UnfoldMisID
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-unfolded-2016)
	@mkdir -p $(OUT_DIR)
	$< --debug --iteration 20 \
		-e ./histos/rdx-22_04_13_23_16-merged-2016/merged.root \
		-y ./histos/rdx-22_04_12_14_03-tag-2016/tagged.root \
		-o $(OUT_DIR) \
		-c ./spec/rdx-run2.yml | tee $(OUT_DIR)/stdout.log

build-rdx-weights-2016: \
	$(BINPATH)/ApplyMisIDWeight \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_02--mu_misid--data--2016--md.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-weights-2016)
	$(eval AUX_NTP	:=	$(basename $(notdir $(word 2, $^)))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	$< -Y 2016 -i $(word 2, $^) -o $(OUT_DIR)/$(AUX_NTP) -c ./spec/rdx-run2.yml

build-generic-dif-smearing:
	$(eval OUT_DIR	:=	$(GENPATH)/generic-$(TIME_STAMP)-dif_smearing)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_dif.py --plot -o $(OUT_DIR) \
		-k ./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		-p ./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root


############
# RDX plot #
############
.PHONY: plot-rdx-bin_vars-2016 plot-rdx-fit_vars-2016

plot-rdx-bin_vars-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-bin_vars-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_histo.py -o $(OUT_DIR) \
		-i ./histos/rdx-22_04_15_01_04-unfolded-2016/unfolded.root
	./scripts/plot_histo.py -o $(OUT_DIR) -s Tag \
		-i ./histos/rdx-22_04_12_14_03-tag-2016/tagged.root

plot-rdx-bin_vars-ana-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-bin_vars-ana-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_histo.py -o $(OUT_DIR) -s Tag \
		-i ./histos/rdx-22_04_12_14_03-tag-2016/tagged.root \
		--show-title 0 1 0 --show-legend 1 0 0
	./scripts/plot_histo.py -o $(OUT_DIR) \
		-i ./histos/rdx-22_04_15_01_04-unfolded-2016/unfolded.root \
		--show-title 0 1 0 --show-legend 1 0 0


plot-rdx-fit_vars-2016: \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_02--mu_misid--data--2016--md.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-fit_vars-2016)
	$(eval AUX_NTP	:=	$(basename $<)--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_fit_vars.py -o $(OUT_DIR) -i $< -a $(AUX_NTP)


########
# Test #
########
.PHONY: test-pidcalib2-wrapper test-unfold

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2.yml -o $(GENPATH) --dry-run

test-unfold: $(BINPATH)/UnfoldMisID
	$< -c ./spec/rdx-run2.yml --dryRun \
        -y ./histos/rdx-22_04_22_00_45-tag-2016/tagged.root \
        -e ./histos/rdx-22_04_12_14_03-merged-2016/merged.root

test-gen-aux-filename: ./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_02--mu_misid--data--2016--md.root
	$(eval AUX_NTP	:=	$(basename $(notdir $^))--aux_misid.root)
	@echo $(AUX_NTP)

test-get-particle-name: ./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_02--mu_misid--data--2016--md.root
	$(eval NTP_NAME	:=	$(basename $(notdir $^)))
	$(eval PARTICLE	:=	$(word 1, $(subst --, ,$(NTP_NAME))))
	@echo $(NTP_NAME)
	@echo $(PARTICLE)


###############
# Compile C++ #
###############

$(BINPATH)/%: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) $(ADDLINKFLAGS)


#######
# Doc #
#######
.PHONY: doc

doc: $(PDF_FILES)

$(GENPATH)/%.pdf: %.tex
	latexmk $<
