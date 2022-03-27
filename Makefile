# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sat Mar 26, 2022 at 08:35 PM -0400

TIME_STAMP	:=	$(shell date +"%y_%m_%d_%H_%M")


.PHONY: build-tagged-histo

build-tagged-histo:
	./scripts/build_histo_tagged.py -c ./spec/run2-rdx.yml -o ./gen


.PHONY: lxplus-build-pidcalib-histo

lxplus-build-pidcalib-histo:
	$(eval OUT_DIR	:=	gen/pidcalib-$(TIME_STAMP))
	./scripts/pidcalib_wrapper.py -c ./spec/run2-rdx.yml -o $(OUT_DIR)


.PHONY: test-pidcalib2-wrapper

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/run2-rdx.yml -o ./gen --dry-run
