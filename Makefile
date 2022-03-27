# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sun Mar 27, 2022 at 05:54 PM -0400

TIME_STAMP	:=	$(shell date +"%y_%m_%d_%H_%M")

###########
# General #
###########

.PHONY: clean

clean:
	@rm -rf ./gen/*


#######
# RDX #
#######

.PHONY: build-tagged-histo

build-rdx-tag-2016:
	$(eval OUT_DIR	:=	gen/rdx-$(TIME_STAMP)-tag-2016)
	./scripts/build_histo_tagged.py -c ./spec/run2-rdx.yml -o ./gen


.PHONY: lxplus-build-pidcalib-histo-2016

build-rdx-true-to-tag-2016:
	$(eval OUT_DIR	:=	gen/rdx-$(TIME_STAMP)-true_to_tag-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/run2-rdx.yml -o $(OUT_DIR) -y 16


########
# Test #
########

.PHONY: test-pidcalib2-wrapper

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/run2-rdx.yml -o ./gen --dry-run
