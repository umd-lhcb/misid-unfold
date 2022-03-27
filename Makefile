# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sat Mar 26, 2022 at 10:54 PM -0400

TIME_STAMP	:=	$(shell date +"%y_%m_%d_%H_%M")


.PHONY: clean

clean:
	@rm -rf ./gen/*


.PHONY: build-tagged-histo

build-tagged-histo:
	./scripts/build_histo_tagged.py -c ./spec/run2-rdx.yml -o ./gen


.PHONY: lxplus-build-pidcalib-histo-2016

lxplus-build-pidcalib-histo-2016:
	$(eval OUT_DIR	:=	gen/pidcalib-$(TIME_STAMP)-2016)
	./scripts/pidcalib_wrapper.py -c ./spec/run2-rdx.yml -o $(OUT_DIR) -y 2016


.PHONY: test-pidcalib2-wrapper

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/run2-rdx.yml -o ./gen --dry-run
