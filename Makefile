# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sun Mar 27, 2022 at 06:03 PM -0400

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

.PHONY: build-tagged-histo build-rdx-true-to-tag-2016

build-rdx-tag-2016:
	$(eval OUT_DIR	:=	./gen/rdx-$(TIME_STAMP)-tag-2016)
	./scripts/build_histo_tagged.py -c ./spec/rdx-2016.yml -o $(OUT_DIR)

build-rdx-true-to-tag-2016:
	$(eval OUT_DIR	:=	./gen/rdx-$(TIME_STAMP)-true_to_tag-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-2016.yml -o $(OUT_DIR) -y 16


########
# Test #
########

.PHONY: test-pidcalib2-wrapper

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/rdx-2016.yml -o ./gen --dry-run
