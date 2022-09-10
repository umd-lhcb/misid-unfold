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


## The misID unfolding procedure

Read `spec/rdx-run2.yml` to see all required inputs!

### Generation of true to tag misID efficiencies

1. Update $K, \pi, p$ efficiencies w/ _modified_ `pidcalib` samples:

     _modified_ refers to the `pidcalib` ntuples w/ UBDT branches, which are
     only available at `glacier`

    1. Clone this project on `glacier`
    2. Run `nix develop` in the project root
    3. In the resulting shell, run `make build-rdx-true-to-tag-2016-glacier`

2. Update $e$ efficiency w/ _original_ `pidcalib` samples:

    Note that the _original_ samples lack UBDT branch, so we need to add them later

    1. Clone this project on `lxplus`
    2. In the project root, run `make build-rdx-true-to-tag-2016-lxplus`

3. Update ghost and $e$ conditional efficiencies w/ local $J/\psi$ samples


#### Remark on `cut` vs. `pid_cut` options in the YAML file

These names are consistent with `pidcalib2` naming scheme, that is:

- `cut`: The cuts in the denominator
- `pid_cut`: The _additional_ cuts in the nominator

Therefore, the efficiency is defined as:

$$
\epsilon = \frac{\text{event passing cut and PID cut}}{\text{event passing cut}}
$$

#### Remark on test run

One can test these rules with:

- `make test-pidcalib2-wrapper-glacier`
- `make test-pidcalib2-wrapper-lxplus`

locally, without actually running the efficiency generation program.


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


## Apply misID and DiF smearing weights

To apply these weight on some sample ntuple:

```
make test-rdx-weights
```

Some test plots including both the nominal region and the DSB region can be
generated with:

```
make plot-rdx-fit_vars_dsb-ana-2016
```
