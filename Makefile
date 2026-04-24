# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sat Oct 15, 2022 at 12:45 AM -0400

BINPATH := ./bin
GENPATH := ./gen
VPATH   := src:docs
CPP_FILES := $(wildcard src/*.cpp)
EXE_FILES := $(patsubst src/%.cpp,$(BINPATH)/%,$(CPP_FILES))
TEX_FILES := $(wildcard docs/*.tex)
PDF_FILES := $(patsubst docs/%.tex,$(GENPATH)/%.pdf,$(TEX_FILES))

YML_FILE  := ./spec/rdx-run2.yml
YML_FILE_RUN2ANG  := ./spec/rdx-run2-run2ang.yml

YML_PIDCALIB_EFFS := ./spec/pidcalib_sel_effs.yml
YML_PIDCALIB_EFFS_RUN2ANG := ./spec/pidcalib_sel_effs_run2ang.yml

YML_D0_BKG := ./spec/d0_bkg_constraints.yml
YML_D0_BKG_RUN2ANG := ./spec/d0_bkg_constraints_run2ang.yml

TIME_STAMP	:=	$(shell date +"%y_%m_%d_%H_%M")

COMPILER     := $(shell root-config --cxx)
CXXFLAGS     := $(shell root-config --cflags) -Iinclude
LINKFLAGS    := $(shell root-config --libs)
ADDCXXFLAGS  := -O3 -march=native -mtune=native -Wall
ADDLINKFLAGS := -lyaml-cpp -lRooFitCore -lRooFit

OS := $(shell uname)
ifeq ($(OS),Darwin)
$(info OS is $(OS) (macOS), adding -lc++fs to LINKFLAGS)
	LINKFLAGS := $(shell root-config --libs) -lc++fs
endif


#################
# Configuration #
#################

EFFICIENCIES := ./histos/rdx-26_04_24_06_52-merged/merged.root
UNFOLDED := ./histos/rdx-26_04_24_07_02-unfolded/unfolded.root
TAGGED := ./histos/rdx-26_04_14_07_02-tag/tagged.root
DIF    := ./histos/generic-26_04_14_06_57-dif_smearing/dif.root


###########
# General #
###########
.PHONY: clean build-test

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

# K, pi and p efficiencies with PIDCalib using Run 1 RDx PID cuts
.PHONY: build-rdx-true-to-tag-2016-glacier build-rdx-true-to-tag-2017-glacier build-rdx-true-to-tag-2018-glacier
build-rdx-true-to-tag-2016-glacier:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_glacier-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(OUT_DIR) -y 2016 -m glacier -u run1rdx 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2017-glacier:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_glacier-2017)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(OUT_DIR) -y 2017 -m glacier -u run1rdx 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2018-glacier:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_glacier-2018)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(OUT_DIR) -y 2018 -m glacier -u run1rdx 2>&1 | tee $(OUT_DIR)/stdout.log

# K, pi and p efficiencies with PIDCalib using Run 2 Ang PID cuts
.PHONY: build-rdx-true-to-tag-2016-glacier-run2ang build-rdx-true-to-tag-2017-glacier-run2ang build-rdx-true-to-tag-2018-glacier-run2ang
build-rdx-true-to-tag-2016-glacier-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_glacier_run2ang-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2016 -m glacier -u run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2017-glacier-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_glacier_run2ang-2017)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2017 -m glacier -u run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2018-glacier-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_glacier_run2ang-2018)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2018 -m glacier -u run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

# e efficiencies with PIDCalib using Run 1 RDx PID cuts (not accounting for mu_uBDT cut)
.PHONY: build-rdx-true-to-tag-2016-lxplus build-rdx-true-to-tag-2017-lxplus build-rdx-true-to-tag-2018-lxplus
build-rdx-true-to-tag-2016-lxplus:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(OUT_DIR) -y 2016 -m lxplus 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2017-lxplus:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus-2017)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(OUT_DIR) -y 2017 -m lxplus 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2018-lxplus:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus-2018)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(OUT_DIR) -y 2018 -m lxplus 2>&1 | tee $(OUT_DIR)/stdout.log

# e efficiencies with PIDCalib using Run 1 RDx PID cuts (not accounting for mu_uBDT cut)
.PHONY: build-rdx-true-to-tag-2016-lxplus-run2ang build-rdx-true-to-tag-2017-lxplus-run2ang build-rdx-true-to-tag-2018-lxplus-run2ang
build-rdx-true-to-tag-2016-lxplus-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus_run2ang-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2016 -m lxplus 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2017-lxplus-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus_run2ang-2017)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2017 -m lxplus 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2018-lxplus-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_lxplus_run2ang-2018)
	@mkdir -p $(OUT_DIR)
	./scripts/pidcalib_wrapper.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2018 -m lxplus 2>&1 | tee $(OUT_DIR)/stdout.log

.PHONY: test-pidcalib2-wrapper-glacier test-pidcalib2-wrapper-glacier-run2ang test-pidcalib2-wrapper-lxplus
test-pidcalib2-wrapper-glacier:
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(GENPATH) --dry-run -m glacier -u run1rdx

test-pidcalib2-wrapper-glacier-run2ang:
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(GENPATH) --dry-run -m glacier -u run2ang

test-pidcalib2-wrapper-lxplus:
	./scripts/pidcalib_wrapper.py -c $(YML_FILE) -o $(GENPATH) --dry-run -m lxplus

# Ghost efficiencies and e conditional efficiencies using Run 1 RDx PID cuts
.PHONY: build-rdx-true-to-tag-2016-local build-rdx-true-to-tag-2017-local build-rdx-true-to-tag-2018-local
build-rdx-true-to-tag-2016-local:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_eff.py -c $(YML_FILE) -o $(OUT_DIR) -y 2016 -u run1rdx 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2017-local:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local-2017)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_eff.py -c $(YML_FILE) -o $(OUT_DIR) -y 2017 -u run1rdx 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2018-local:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local-2018)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_eff.py -c $(YML_FILE) -o $(OUT_DIR) -y 2018 -u run1rdx 2>&1 | tee $(OUT_DIR)/stdout.log

# Ghost efficiencies and e conditional efficiencies using Run 2 Angx PID cuts
.PHONY: build-rdx-true-to-tag-2016-local-run2ang build-rdx-true-to-tag-2017-local-run2ang build-rdx-true-to-tag-2018-local-run2ang
build-rdx-true-to-tag-2016-local-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local_run2ang-2016)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_eff.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2016 -u run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2017-local-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local_run2ang-2017)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_eff.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2017 -u run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-true-to-tag-2018-local-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-true_to_tag_local_run2ang-2018)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_eff.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) -y 2018 -u run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

# Starting here are the steps needed to perform the fits to the pidcalib calibration samples

# Check effect of WM cut on pidcalib comb. bkg and determine fit ranges
# Seems unaffected by muon PID cuts
build-check-comb: $(BINPATH)/CheckCombBkg
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-check_comb)
	@mkdir -p $(OUT_DIR)/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) 2>&1 | tee $(OUT_DIR)/stdout.log

# Produce bkg constraints needed to calculate K/pi -> mu misid efficiencies
build-d0-bkg-decays: $(BINPATH)/d0BkgDecays
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-d0_bkg_decays)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE) -f d0_bkg_constraints.yml 2>&1 | tee $(OUT_DIR)/stdout.log

build-d0-bkg-decays-run2ang: $(BINPATH)/d0BkgDecays
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-d0_bkg_decays_run2ang)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -f d0_bkg_constraints_run2ang.yml --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

# Produce pidcalib selection efficiencies
build-pidcalib-sel-effs: $(BINPATH)/PIDCalibSelEffs
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-pidcalib_sel_effs)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE) -f pidcalib_sel_effs.yml 2>&1 | tee $(OUT_DIR)/stdout.log

build-pidcalib-sel-effs-run2ang: $(BINPATH)/PIDCalibSelEffs
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-pidcalib_sel_effs_run2ang)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -f pidcalib_sel_effs_run2ang.yml --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

# K/pi -> mu misid ISO+CTRL efficiencies
build-rdx-misid-mc-corrections: \
	build-rdx-misid-mc-corrections-2016 \
	build-rdx-misid-mc-corrections-2017 \
	build-rdx-misid-mc-corrections-2018

build-rdx-misid-mc-corrections-2016: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2016)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2016 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2017: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2017)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2017 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2018: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2018)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2018 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-run2ang: \
	build-rdx-misid-mc-corrections-2016-run2ang \
	build-rdx-misid-mc-corrections-2017-run2ang \
	build-rdx-misid-mc-corrections-2018-run2ang

build-rdx-misid-mc-corrections-2016-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2016)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2016 --run2ang --float_rich 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2017-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2017)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2017 --run2ang --float_rich 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2018-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2018)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2018 --run2ang --float_rich 2>&1 | tee $(OUT_DIR)/stdout.log

# K/pi -> mu misid VMU efficiencies
build-rdx-misid-mc-corrections-vmu: \
	build-rdx-misid-mc-corrections-2016-vmu \
	build-rdx-misid-mc-corrections-2017-vmu \
	build-rdx-misid-mc-corrections-2018-vmu

build-rdx-misid-mc-corrections-2016-vmu: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2016-vmu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2016 --vmu 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2017-vmu: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2017-vmu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2017 --vmu 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2018-vmu: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2018-vmu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2018 --vmu 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-vmu-run2ang: \
	build-rdx-misid-mc-corrections-2016-vmu-run2ang \
	build-rdx-misid-mc-corrections-2017-vmu-run2ang \
	build-rdx-misid-mc-corrections-2018-vmu-run2ang

build-rdx-misid-mc-corrections-2016-vmu-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2016-vmu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2016 --vmu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2017-vmu-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2017-vmu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2017 --vmu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2018-vmu-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2018-vmu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2018 --vmu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

# K/pi -> mu misid FAKE_MU efficiencies
build-rdx-misid-mc-corrections-fake_mu: \
	build-rdx-misid-mc-corrections-2016-fake_mu \
	build-rdx-misid-mc-corrections-2017-fake_mu \
	build-rdx-misid-mc-corrections-2018-fake_mu

build-rdx-misid-mc-corrections-2016-fake_mu: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2016-fake_mu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2016 --fake_mu 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2017-fake_mu: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2017-fake_mu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2017 --fake_mu 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2018-fake_mu: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-2018-fake_mu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE) -b $(YML_D0_BKG) -s $(YML_PIDCALIB_EFFS) -y 2018 --fake_mu 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-fake_mu-run2ang: \
	build-rdx-misid-mc-corrections-2016-fake_mu-run2ang \
	build-rdx-misid-mc-corrections-2017-fake_mu-run2ang \
	build-rdx-misid-mc-corrections-2018-fake_mu-run2ang

build-rdx-misid-mc-corrections-2016-fake_mu-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2016-fake_mu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2016 --fake_mu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2017-fake_mu-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2017-fake_mu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2017 --fake_mu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-misid-mc-corrections-2018-fake_mu-run2ang: $(BINPATH)/GetMisIDCorrections
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-mc-corrections-run2ang-2018-fake_mu)
	@mkdir -p $(OUT_DIR)/figs/params
	@mkdir -p $(OUT_DIR)/figs/fits
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) -b $(YML_D0_BKG_RUN2ANG) -s $(YML_PIDCALIB_EFFS_RUN2ANG) -y 2018 --fake_mu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

# Merge efficiencies into single file
.PHONY: build-rdx-merged build-rdx-merged-run2ang
build-rdx-merged:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-merged)
	@mkdir -p $(OUT_DIR)
	./scripts/merge_histo.py -c $(YML_FILE) -o $(OUT_DIR) 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-merged-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-merged-run2ang)
	@mkdir -p $(OUT_DIR)
	./scripts/merge_histo.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) 2>&1 | tee $(OUT_DIR)/stdout.log

# Produce DiF smearing sample
.PHONY: build-generic-dif-smearing build-generic-dif-smearing-run2ang
build-generic-dif-smearing:
	$(eval OUT_DIR	:=	$(GENPATH)/generic-$(TIME_STAMP)-dif_smearing)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_dif.py --plot -o $(OUT_DIR) \
		./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root \
		2>&1 | tee $(OUT_DIR)/stdout.log

build-generic-dif-smearing-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/generic-$(TIME_STAMP)-dif_smearing-run2ang)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_dif.py --run2ang --plot -o $(OUT_DIR) \
		./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/K-mix/K--17_06_28--mix--2011-2012--md-mu--greg.root \
		./ntuples/ref-rdx-run1/Pi-mix/Pi--17_06_28--mix--2011-2012--md-mu--greg.root \
		2>&1 | tee $(OUT_DIR)/stdout.log

# Create histograms with tagged yields
.PHONY: build-rdx-tag build-rdx-tag-run2ang
build-rdx-tag:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-tag)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_tagged.py -c $(YML_FILE) -o $(OUT_DIR) 2>&1 | tee $(OUT_DIR)/stdout.log

build-rdx-tag-run2ang:
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-tag-run2ang)
	@mkdir -p $(OUT_DIR)
	./scripts/build_histo_tagged.py -c $(YML_FILE_RUN2ANG) -o $(OUT_DIR) --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

# Unfold the misID weights
build-rdx-unfolded: $(BINPATH)/UnfoldMisID
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-unfolded)
	@mkdir -p $(OUT_DIR)
	$< --iteration 5 -e $(EFFICIENCIES) -y $(TAGGED) -o $(OUT_DIR) \
		-c $(YML_FILE) 2>&1 | tee $(OUT_DIR)/stdout.log

test-unfold: $(BINPATH)/UnfoldMisID
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-test-unfold)
	@mkdir -p $(OUT_DIR)
	$< -c $(YML_FILE) --dryRun --debug -y $(TAGGED) -e $(EFFICIENCIES) -o $(OUT_DIR)  | tee $(OUT_DIR)/stdout.log


# Test of application of misID weights
test-apply-rdx-weights: test-apply-rdx-weights-2016 test-apply-rdx-weights-2017 test-apply-rdx-weights-2018

test-apply-rdx-weights-2016: $(BINPATH)/ApplyMisIDWeight \
	./ntuples/0.9.18-misid_pid_kept/data/Dst_D0--25_09_15--mu_misid--LHCb_Collision16_Beam6500GeV-VeloClosed-MagDown_Real_Data_Reco16_Stripping28r2_90000000_SEMILEPTONIC.DST--000-dv.root \
	$(DIF)
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-weights-2016)
	$(eval AUX_NTP	:=	$(basename $(notdir $(word 2, $^)))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	$< --debug -a -Y 2016 -i $(word 2, $^) -x $(word 3, $^) \
		--kSmrBrName k_smr --piSmrBrName pi_smr \
		-o $(OUT_DIR)/$(AUX_NTP) -c $(YML_FILE) | tee $(OUT_DIR)/stdout.log

test-apply-rdx-weights-2017: $(BINPATH)/ApplyMisIDWeight \
	./ntuples/0.9.18-misid_pid_kept/data/Dst_D0--25_09_15--mu_misid--LHCb_Collision17_Beam6500GeV-VeloClosed-MagDown_Real_Data_Reco17_Stripping29r2_90000000_SEMILEPTONIC.DST--000-dv.root \
	$(DIF)
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-weights-2017)
	$(eval AUX_NTP	:=	$(basename $(notdir $(word 2, $^)))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	$< --debug -a -Y 2017 -i $(word 2, $^) -x $(word 3, $^) \
		--kSmrBrName k_smr --piSmrBrName pi_smr \
		-o $(OUT_DIR)/$(AUX_NTP) -c $(YML_FILE) | tee $(OUT_DIR)/stdout.log

test-apply-rdx-weights-2018: $(BINPATH)/ApplyMisIDWeight \
	./ntuples/0.9.18-misid_pid_kept/data/Dst_D0--25_09_15--mu_misid--LHCb_Collision18_Beam6500GeV-VeloClosed-MagDown_Real_Data_Reco18_Stripping34_90000000_SEMILEPTONIC.DST--000-dv.root \
	$(DIF)
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-weights-2018)
	$(eval AUX_NTP	:=	$(basename $(notdir $(word 2, $^)))--aux_misid.root)
	@mkdir -p $(OUT_DIR)
	$< --debug -a -Y 2018 -i $(word 2, $^) -x $(word 3, $^) \
		--kSmrBrName k_smr --piSmrBrName pi_smr \
		-o $(OUT_DIR)/$(AUX_NTP) -c $(YML_FILE) | tee $(OUT_DIR)/stdout.log


# Aux. efficiencies for iso track PID cuts on true ghost
.PHONY: build-rdx-aux-isoPID-on-ghost build-rdx-aux-probnnk-on-ghost build-rdx-aux-probnnp-on-ghost build-rdx-aux-probnnpnng_gtlt-on-ghost build-rdx-aux-probnnknng_ltlt-on-ghost build-rdx-aux-probnnknng_lt02lt02-on-ghost build-rdx-aux-probnnknng_gtlt-on-ghost build-rdx-aux-probnnknng_gt02lt02-on-ghost
build-rdx-aux-isoPID-on-ghost: build-rdx-aux-probnnk-on-ghost build-rdx-aux-probnnp-on-ghost build-rdx-aux-probnnpnng_gtlt-on-ghost build-rdx-aux-probnnknng_ltlt-on-ghost build-rdx-aux-probnnknng_lt02lt02-on-ghost build-rdx-aux-probnnknng_gtlt-on-ghost build-rdx-aux-probnnknng_gt02lt02-on-ghost
build-rdx-aux-probnnk-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_nnk_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_probnnk.yml -y 2016 -o $(OUT_DIR)_2016
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_probnnk.yml -y 2017 -o $(OUT_DIR)_2017
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_probnnk.yml -y 2018 -o $(OUT_DIR)_2018
build-rdx-aux-probnnp-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_nnp_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_probnnp.yml -y 2016 -o $(OUT_DIR)_2016
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_probnnp.yml -y 2017 -o $(OUT_DIR)_2017
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_probnnp.yml -y 2018 -o $(OUT_DIR)_2018
build-rdx-aux-probnnpnng_gtlt-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_nnpnng_gtlt_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnpnng_gtlt.yml -y 2016 -o $(OUT_DIR)_2016
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnpnng_gtlt.yml -y 2017 -o $(OUT_DIR)_2017
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnpnng_gtlt.yml -y 2018 -o $(OUT_DIR)_2018
build-rdx-aux-probnnknng_ltlt-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_nnknng_ltlt_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_ltlt.yml -y 2016 -o $(OUT_DIR)_2016
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_ltlt.yml -y 2017 -o $(OUT_DIR)_2017
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_ltlt.yml -y 2018 -o $(OUT_DIR)_2018
build-rdx-aux-probnnknng_lt02lt02-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_nnknng_lt02lt02_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_lt02lt02.yml -y 2016 -o $(OUT_DIR)_2016
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_lt02lt02.yml -y 2017 -o $(OUT_DIR)_2017
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_lt02lt02.yml -y 2018 -o $(OUT_DIR)_2018
build-rdx-aux-probnnknng_gtlt-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_nnknng_gtlt_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_gtlt.yml -y 2016 -o $(OUT_DIR)_2016
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_gtlt.yml -y 2017 -o $(OUT_DIR)_2017
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_gtlt.yml -y 2018 -o $(OUT_DIR)_2018
build-rdx-aux-probnnknng_gt02lt02-on-ghost:
	$(eval OUT_DIR	:=	$(GENPATH)/root-run2-rdx_iso_nnknng_gt02lt02_ghost)
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_gt02lt02.yml -y 2016 -o $(OUT_DIR)_2016
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_gt02lt02.yml -y 2017 -o $(OUT_DIR)_2017
	@./scripts/build_histo_eff.py -c ./spec/rdx-run2-iso_nnknng_gt02lt02.yml -y 2018 -o $(OUT_DIR)_2018


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

misid-plots: $(BINPATH)/MisIDPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-plots)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE) 2>&1 | tee $(OUT_DIR)/stdout.log

misid-plots-vmu: $(BINPATH)/MisIDPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-plots-vmu)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE) --vmu 2>&1 | tee $(OUT_DIR)/stdout.log

misid-plots-fake_mu: $(BINPATH)/MisIDPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-plots-fake_mu)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE) --fake_mu 2>&1 | tee $(OUT_DIR)/stdout.log

misid-plots-run2ang: $(BINPATH)/MisIDPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-plots-run2ang)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

misid-plots-vmu-run2ang: $(BINPATH)/MisIDPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-plots-vmu-run2ang)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) --vmu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

misid-plots-fake_mu-run2ang: $(BINPATH)/MisIDPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-misid-plots-fake_mu-run2ang)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) --fake_mu --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

pidcalib-plots: $(BINPATH)/MakePIDCalibPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-pidcalib-plots)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE) 2>&1 | tee $(OUT_DIR)/stdout.log

pidcalib-plots-run2ang: $(BINPATH)/MakePIDCalibPlots
	$(eval OUT_DIR	:=	$(GENPATH)/rdx-$(TIME_STAMP)-pidcalib-plots-run2ang)
	@mkdir -p $(OUT_DIR)
	$< -o $(OUT_DIR) -c $(YML_FILE_RUN2ANG) --run2ang 2>&1 | tee $(OUT_DIR)/stdout.log

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

$(BINPATH)/compareEffs: compareEffs.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS)

$(BINPATH)/UnfoldMisID: UnfoldMisID.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) $(ADDLINKFLAGS) -lRooUnfold

$(BINPATH)/%: %.cpp
	$(COMPILER) $(CXXFLAGS) $(ADDCXXFLAGS) -o $@ $< $(LINKFLAGS) $(ADDLINKFLAGS)


#######
# Doc #
#######
.PHONY: doc

doc: $(PDF_FILES)

$(GENPATH)/%.pdf: %.tex
	latexmk $<
