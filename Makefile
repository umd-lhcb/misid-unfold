# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Sat Mar 26, 2022 at 04:58 AM -0400

.PHONY: build-tagged-histo

build-tagged-histo:
	./scripts/build_histo_tagged.py -c ./spec/run2-rdx.yml -o ./gen


.PHONY: test-pidcalib2-wrapper

test-pidcalib2-wrapper:
	./scripts/pidcalib_wrapper.py -c ./spec/run2-rdx.yml -o ./gen --dry-run
