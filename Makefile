# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Mon Sep 12, 2022 at 05:06 AM -0400

BINPATH := ./bin
GENPATH := ./gen
VPATH := src:docs
CPP_FILES	:=	$(wildcard src/*.cpp)
EXE_FILES	:=	$(patsubst src/%.cpp,$(BINPATH)/%,$(CPP_FILES))
TEX_FILES	:=	$(wildcard docs/*.tex)
PDF_FILES	:=	$(patsubst docs/%.tex,$(GENPATH)/%.pdf,$(TEX_FILES))

TIME_STAMP	:=	$(shell date +"%y_%m_%d_%H_%M")

COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags) -Iinclude
LINKFLAGS	:=	$(shell root-config --libs)
ADDCXXFLAGS	:=	-O2 -march=native -mtune=native
ADDLINKFLAGS	:=	-lyaml-cpp -lRooFitCore -lRooFit -lRooStats -lRooUnfold


###########
# General #
###########
.PHONY: exe applyer clean build-test build-test-lxplus build-test-nix plot-test

exe: $(EXE_FILES)

applyer: $(BINPATH)/ApplyMisIDWeight

clean:
	@rm -rf $(GENPATH)/*

test-nix:
	nix build ".#misid-unfold-applyer"


#######
# RDX #
#######

# Generation of true to tag misID efficiency
.PHONY: build-rdx-true-to-tag-2016-glacier build-rdx-true-to-tag-2016-lxplus
build-rdx-true-to-tag-2016-glacier:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_glacier-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016 -m glacier

build-rdx-true-to-tag-2016-lxplus:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016 -m lxplus


.PHONY: test-pidcalib2-wrapper-glacier test-pidcalib2-wrapper-lxplus
test-pidcalib2-wrapper-glacier:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2.yml -o $(GENPATH) --dry-run -m glacier

test-pidcalib2-wrapper-lxplus:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2.yml -o $(GENPATH) --dry-run -m lxplus


.PHONY: build-rdx-true-to-tag-2016-local
build-rdx-true-to-tag-2016-local:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local-2016)
	./scripts/build_histo_eff.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016


.PHONY: build-rdx-merged-2016
build-rdx-merged-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-merged-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/merge_histo.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016 | tee $(OUT_DIR)/stdout.log


.PHONY: build-generic-dif-smearing
build-generic-dif-smearing:
	$(eval OUT_DIR	:=	$(GENPATH)/generic-$(TIME_STAMP)-dif_smearing)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_dif.py --plot -o $(OUT_DIR) \
		-k ./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		-p ./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root


# Unfold the misID weights
.PHONY: build-rdx-tag-2016
build-rdx-tag-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-tag-2016)
	./scripts/build_histo_tagged.py -c ./spec/rdx-run2.yml -o $(OUT_DIR) -y 2016

build-rdx-unfolded-2016: $(BINPATH)/UnfoldMisID
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-unfolded-2016)
	@mkdir -p $(OUT_DIR)
	$< --debug --iteration 20 \
		-e ./histos/rdx-22_09_12_05_03-merged-2016/merged.root \
		-y ./histos/rdx-22_06_23_12_07-tag-2016/tagged.root \
		-o $(OUT_DIR) \
		-c ./spec/rdx-run2.yml | tee $(OUT_DIR)/stdout.log


.PHONY: test-unfold
test-unfold: $(BINPATH)/UnfoldMisID
	$< -c ./spec/rdx-run2.yml --dryRun \
		-y ./histos/rdx-22_06_23_12_07-tag-2016/tagged.root \
		-e ./histos/rdx-22_04_13_23_16-merged-2016/merged.root


# Application of misID weights
.PHONY: build-rdx-weights-2016 build-rdx-weights-2016-md build-rdx-weights-2016-mu
build-rdx-weights-2016: build-rdx-weights-2016-md build-rdx-weights-2016-mu

build-rdx-weights-2016-md: \
	$(BINPATH)/ApplyMisIDWeight \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--md.root \
	./histos/generic-22_04_21_20_12-dif_smearing/dif.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-weights-2016)
	$(eval AUX_NTP	:=	$(basename $(notdir $(word 2, $^)))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	$< -Y 2016 -i $(word 2, $^) -x $(word 3, $^) -o $(OUT_DIR)/$(AUX_NTP) -c ./spec/rdx-run2.yml | tee $(OUT_DIR)/stdout.log

build-rdx-weights-2016-mu: \
	$(BINPATH)/ApplyMisIDWeight \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--mu.root \
	./histos/generic-22_04_21_20_12-dif_smearing/dif.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-weights-2016)
	$(eval AUX_NTP	:=	$(basename $(notdir $(word 2, $^)))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	$< -Y 2016 -i $(word 2, $^) -x $(word 3, $^) -o $(OUT_DIR)/$(AUX_NTP) -c ./spec/rdx-run2.yml | tee $(OUT_DIR)/stdout.log


###############################
# Ghost efficiency comparison #
###############################
.PHONY: ghost-eff-gen

# TODO: The input files are removed. Update needed
# you probably want to refer back to 0.2.4
ghost-eff-gen:
	$(eval OUT_DIR	:=	$(GENPATH)/ghost_eff-$(TIME_STAMP)-true_to_tag-2016)
	./scripts/build_histo_eff.py -c ./spec/rdx-ghost_eff_comp.yml -o $(OUT_DIR) -y 2016


############
# RDX plot #
############
.PHONY: plot-rdx-bin_vars-2016 plot-rdx-fit_vars-2016 \
	plot-rdx-bin_vars-ana-2016 plot-rdx-fit_vars-ana-2016 \
	plot-rdx-fit_vars_dsb-ana-2016

plot-rdx-bin_vars-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-bin_vars-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_histo.py -o $(OUT_DIR) -s Tag \
		-i ./histos/rdx-22_06_23_12_07-tag-2016/tagged.root \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb
	./scripts/plot_histo.py -o $(OUT_DIR) \
		-i ./histos/rdx-22_09_12_05_03-unfolded-2016/unfolded.root \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb

plot-rdx-bin_vars-ana-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-bin_vars-ana-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_histo.py -o $(OUT_DIR) -s Tag --extension pdf \
		-i ./histos/rdx-22_06_23_12_07-tag-2016/tagged.root \
		--show-title 0 1 0 --show-legend 1 0 0 \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb
	./scripts/plot_histo.py -o $(OUT_DIR) --extension pdf \
		-i ./histos/rdx-22_09_12_05_03-unfolded-2016/unfolded.root \
		--show-title 0 1 0 --show-legend 1 0 0 \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb


plot-rdx-fit_vars-2016: \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--md.root \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--mu.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-fit_vars-2016)
	$(eval AUX_NTP_MD	:=	$(basename $<)--aux_misid.root)
	$(eval AUX_NTP_MU	:=	$(basename $(word 2, $^))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_fit_vars.py -o $(OUT_DIR) -i $< $(word 2, $^) -a $(AUX_NTP_MD) $(AUX_NTP_MU)

plot-rdx-fit_vars-ana-2016: \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--md.root \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--mu.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-fit_vars-ana-2016)
	$(eval AUX_NTP_MD	:=	$(basename $<)--aux_misid.root)
	$(eval AUX_NTP_MU	:=	$(basename $(word 2, $^))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_fit_vars.py -o $(OUT_DIR) -i $< $(word 2, $^) -a $(AUX_NTP_MD) $(AUX_NTP_MU) \
		--show-title 1 0 0 --show-legend 1 0 0


plot-rdx-fit_vars_dsb-ana-2016: \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--md.root \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid-study-step2/D0--22_04_24--mu_misid--data--2016--mu.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-fit_vars_dsb-ana-2016)
	$(eval AUX_NTP_MD	:=	$(basename $<)--aux_misid.root)
	$(eval AUX_NTP_MU	:=	$(basename $(word 2, $^))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_fit_vars.py -o $(OUT_DIR) -i $< $(word 2, $^) -a $(AUX_NTP_MD) $(AUX_NTP_MU) \
		--show-title 1 0 0 --show-legend 1 0 0 \
		--cuts d0_m_ok b_m_ok in_fit_range


########
# Test #
########
.PHONY: test-gen-aux-filename test-get-particle-name

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

$(BINPATH)/ApplyMisIDWeight: ApplyMisIDWeight.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) -lyaml-cpp

$(BINPATH)/%: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) $(ADDLINKFLAGS)


#######
# Doc #
#######
.PHONY: doc

doc: $(PDF_FILES)

$(GENPATH)/%.pdf: %.tex
	latexmk $<
