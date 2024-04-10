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

Read `spec/rdx-run2.yml` for more info!
The output folders are located in `gen` folder.

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

    Note that you might need to use a specific version of `pidcalib`. Edit `cmd_prefix` in `scripts/pidcalib_wrapper.py` and replace it with:
    ```
    cmd_prefix = "lb-conda pidcalib/2022-09-02 " if mode == "lxplus" else ""
    ```
    `pidcalib/2022-09-02` corresponds to version `1.0.6`.

4. Update ghost and $e$ conditional efficiencies w/ local incl. $J/\psi$ MC
   samples w/ `make build-rdx-true-to-tag-2016-local`

    Starting from this step, all operations are done locally or on `glacier`.
    And `nix develop` is always assumed

5. Merge all efficiencies obtained above into a single file locally
    w/ `make build-rdx-merged-2016`

    Copy output folders of previous steps from `gen` to `histo` folder, and
    configure the YAML file (the `input_histos` section) before you run the
    `make` command!

6. Extract $K, \pi$ momentum smearing info into a small ntuple
    w/ `make build-generic-dif-smearing`

#### Remark on negative efficiencies

Sometimes the efficiencies are negative, due to sWeight. To correct for that,
we followed a procedure suggested by Phoebe. For more info, see p. 3 of [this slide](https://github.com/umd-lhcb/group-talks/blob/master/phys_group_meetings/22-04-13_yipeng_rdx_status.pdf).

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


### Generation of tagging efficiencies

This step relies on our misID ntuple with a special process which add
PIDCalib-like branches like `Brunel_ProbNNghost` (which is just an alias
to `mu_ProbNNghost`).

The special process is needed so that we can apply the same cuts on
these ntuples as if we are applying cuts on PIDCalib ntuples, eliminating
the need to translate branch names.

To build these special ntuples, go to `lhcb-ntuples-gen`, then type:

```
make rdx-ntuple-run2-misid_study
```

Once that is done, configure `spec/rdx-run2.yml` properly, then:

```
build-rdx-tag-2016
```


### Unfold the misID efficiencies

Make sure to download all required ntuples with `git annex get`!
Also, make sure the input files are set correct for unfolding
in the `Makefile`.

```
make build-rdx-unfolded-2016
```

One can also test-run the unfolding (w/o actually doing unfolding!)
by `make test-unfold`.

### Apply misID and DiF smearing weights

```
make build-rdx-weights-2016
```

This apply misID weights and momentum smearing on misID control samples.
It's typically not applied here. The rule here is for demo purpose only.

::info
**Mind the YML file**
By default, the ```make``` commands listed above use the ```rdx-run2.yml``` configuration file, which contains the signal uBDT cut.
To use the misid uBDT cut, the ```USE_CTRL_SAMPLE``` variable must be set appropriately:

```
make build-rdx-weights-2016 USE_CTRL_SAMPLE=true
```

Assigning unexpected values to this variable will generate a warning.
:::
