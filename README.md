# misid-unfold

Unfolding Muon misID samples.


## `git-annex` setup

If you haven't cloned this project, do these steps first:

```shell
git clone git@github.com:umd-lhcb/misid-unfold
cd misid-unfold

# do this only once, right after cloning
git remote add glacier git@10.229.60.85:misid-unfold
git annex init --version=7
git annex sync
```


## Running the misID unfolding procedure

1. Read `spec/rdx-run2.yml` to see all required inputs!

2. **If you want to update the **true -> tag** efficiencies with the official `pidcalib` sample**:
    1. Clone this project on `lxplus`, without setting up `git annex`
    2. `make build-rdx-true-to-tag-2016`
    3. Copy the histograms in `gen` folder to `histos` folder, commit them then update the spec YAML.

3. **Now everything's local**. Make sure to download all required ntuples with `git annex get`.

4. `make build-rdx-tag-2016`, Copy and commit histograms then update YAML (this
    is implied for all subsequent steps).

    This is to build tagged histograms from misID control sample.

5. `make build-rdx-merged-2016`

    This is to merge all `pidcalib` histograms into a single file.

6. `make build-rdx-unfold-2016`

    This is to do actual misID unfolding.

7. `make build-generic-dif-smearing`

    This is to extract `K, pi` momentum smearing info into a small ntuple.

8. Finally, `make build-rdx-weights-2016`

    This apply misID weights and momentum smearing on misID control samples.
    It's typically not applied here. The rule here is for demo purpose only.
