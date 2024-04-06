# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sat Oct 15, 2022 at 12:45 AM -0400

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

OS := $(shell uname)
ifeq ($(OS),Darwin)
  $(info OS is $(OS) (macOS), adding -lc++fs to LINKFLAGS)
  LINKFLAGS := $(shell root-config --libs) -lc++fs
endif

#################
# Configuration #
#################

EFFICIENCIES	:=	./histos/rdx-22_09_12_05_03-merged-2016/merged.root
TAGGED	:=	./histos/rdx-22_06_23_12_07-tag-2016/tagged.root
UNFOLDED	:=	./histos/rdx-22_10_15_00_44-unfolded-2016/unfolded.root


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
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2-misid_ctrl.yml -o $(OUT_DIR) -y 2016 -m glacier

build-rdx-true-to-tag-2016-lxplus:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2-misid_ctrl.yml -o $(OUT_DIR) -y 2016 -m lxplus


.PHONY: test-pidcalib2-wrapper-glacier test-pidcalib2-wrapper-lxplus
test-pidcalib2-wrapper-glacier:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2-misid_ctrl.yml -o $(GENPATH) --dry-run -m glacier

test-pidcalib2-wrapper-lxplus:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-run2-misid_ctrl.yml -o $(GENPATH) --dry-run -m lxplus


.PHONY: build-rdx-true-to-tag-2016-local
build-rdx-true-to-tag-2016-local:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local-2016)
	./scripts/build_histo_eff.py -c ./spec/rdx-run2-misid_ctrl.yml -o $(OUT_DIR) -y 2016


.PHONY: build-rdx-merged-2016
build-rdx-merged-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-merged-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/merge_histo.py -c ./spec/rdx-run2-misid_ctrl.yml -o $(OUT_DIR) -y 2016 | tee $(OUT_DIR)/stdout.log


.PHONY: build-generic-dif-smearing
build-generic-dif-smearing:
	$(eval OUT_DIR	:=	$(GENPATH)/generic-$(TIME_STAMP)-dif_smearing)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_dif.py --plot -o $(OUT_DIR) \
		./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root


# Unfold the misID weights
.PHONY: build-rdx-tag-2016
build-rdx-tag-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-tag-2016)
	./scripts/build_histo_tagged.py -c ./spec/rdx-run2-misid_ctrl.yml -o $(OUT_DIR) -y 2016

build-rdx-unfolded-2016: $(BINPATH)/UnfoldMisID
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-unfolded-2016)
	@mkdir -p $(OUT_DIR)
	$< --debug --iteration 5 \
		-e $(EFFICIENCIES) -y $(TAGGED) -o $(OUT_DIR) \
		-c ./spec/rdx-run2-misid_ctrl.yml | tee $(OUT_DIR)/stdout.log


.PHONY: test-unfold
test-unfold: $(BINPATH)/UnfoldMisID
	$< -c ./spec/rdx-run2-misid_ctrl.yml --dryRun \
		-y $(TAGGED) -e $(EFFICIENCIES)


# Test of application of misID weights
.PHONY: test-apply-rdx-weights-2016
test-apply-rdx-weights-2016: \
	$(BINPATH)/ApplyMisIDWeight \
	./ntuples/0.9.6-2016_production/Dst_D0-mu_misid/Dst_D0--22_03_01--mu_misid--LHCb_Collision16_Beam6500GeV-VeloClosed-MagDown_Real_Data_Reco16_Stripping28r2_90000000_SEMILEPTONIC.DST.root \
	./histos/generic-22_09_23_04_48-dif_smearing/dif.root
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-weights-2016)
	$(eval AUX_NTP	:=	$(basename $(notdir $(word 2, $^)))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	$< -a -Y 2016 -i $(word 2, $^) -x $(word 3, $^) \
		--kSmrBrName k_smr --piSmrBrName pi_smr \
		-o $(OUT_DIR)/$(AUX_NTP) -c ./spec/rdx-run2-misid_ctrl.yml | tee $(OUT_DIR)/stdout.log


# Aux. efficiencies for ProbNNk > 2 on true ghost
.PHONY: build-rdx-aux-probnnk-on-ghost
build-rdx-aux-probnnk-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_oldcut_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_probnnk.yml -y 2016 -o $(OUT_DIR)


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
.PHONY: plot-rdx-bin_vars-2016 plot-rdx-bin_vars-ana-2016

plot-rdx-bin_vars-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-bin_vars-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_histo.py -o $(OUT_DIR) -s Tag \
		-i $(TAGGED) \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb
	./scripts/plot_histo.py -o $(OUT_DIR) \
		-i $(UNFOLDED) \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb

plot-rdx-bin_vars-ana-2016:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-bin_vars-ana-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/plot_histo.py -o $(OUT_DIR) -s Tag --extension pdf \
		-i $(TAGGED) \
		--show-title 0 1 0 --show-legend 1 0 0 \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb
	./scripts/plot_histo.py -o $(OUT_DIR) --extension pdf \
		-i $(UNFOLDED) \
		--show-title 0 1 0 --show-legend 1 0 0 \
		-p D0 D0_bsb Dst Dst_bsb Dst_dsb Dst_dsb_bsb Dst_ws_Mu_dsb_bsb


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
