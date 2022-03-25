# Author: Yipeng Sun
# License: BSD 2-clause
# Last Change: Fri Mar 25, 2022 at 12:24 AM -0400

.PHONY: build-tagged-histo

build-tagged-histo:
	./scripts/build_histo_tagged.py -c ./spec/run2-rdx.yml -o ./gen
