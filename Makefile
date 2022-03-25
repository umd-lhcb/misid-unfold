.PHONY: build-tagged-histo

build-tagged-histo:
	./scripts/build_histo_tagged.py -c ./spec/run2-rdx.yml -o ./gen
