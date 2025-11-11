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

1. Generate $K$, $\pi$ and $p$ efficiencies with _modified_ `pidcalib` samples (i.e. ntuples with UBDT branches, which are only available at `glacier`).

    1. Clone this project on `glacier`
    2. Run `nix develop` in the project root
    3. In the resulting shell, run `make build-rdx-true-to-tag-2016-glacier`. The output files will be saved in `gen/rdx-<timestamp>-true_to_tag_glacier-2016`.
    4. Copy the output to `histos/default`.
    5. Also in `histos/default`, make a copy of the folder produced above with `-corrected` appended to its name.

2. Generate $K\to\mu$ and $\pi\to\mu$ misid efficiencies through fits to the modified `pidcalib` samples.

    1. Clone this project on `glacier`
    2. Run `nix develop` in the project root
    3. In the resulting shell, run `make build-d0-decays` to produce constraints for the fits. The output files will be saved in `gen/rdx-<timestamp>-d0_decays`.
    4. Copy `d0_decays.yml` from the output folder to `spec/`.
    5. Now run `make build-rdx-misid-mc-corrections` and `make build-rdx-misid-mc-corrections-fake_mu` to actually perform the fits. The output files will be saved in `gen/rdx-<timestamp>-misid-mc-corrections` and `gen/rdx-<timestamp>-misid-mc-corrections-fake_mu`, respectively.
    6. Copy the `*TrueToMuTag.root` files from the output folders into the `-corrected` folder produced in step $1.5$ along with the other $K$, $\pi$ and $p$ efficiencies. This will overwrite the $K\to\mu$ and $\pi\to\mu$ efficiencies from PIDCalib, which are still preserved in the non`-corrected` folder.

    > [!TIP]
    > The two `make` commands in step $2.5$ may be executed concurrently, and this is probably preferable since each takes about 12 hours!
    > One alternative is to use `screen` sessions:
    >
    > ```shell
    > screen -S misid-iso # Create named screen session
    > cd misid-unfold
    > nix develop
    > make build-rdx-misid-mc-corrections
    > # Detach session with CTRL+A+D
    > # Repeat for `make build-rdx-misid-mc-corrections-fake_mu`
    > ```
    >
    > The screen sessions may be reattached with e.g. `screen -r misid-iso`, and killed from inside with CTRL+A+K.

    > [!IMPORTANT]
    > To run step 2 with different PID cuts, changes must be made to two C++ files: `src/d0BkgDecays.cpp` and `src/GetMisIDCorrection.cpp`.
    >
    > In `src/d0BkgDecays.cpp`, the event loop must be modified around lines 248-259.
    >
    > In `src/GetMisIDCorrection.cpp`, the **three** event loops must be modified similarly, around lines 999-1010, 1336-1347 and 1514-1526.

3. Generate $e$ efficiencies with _original_ `pidcalib` samples (i.e. ntuples that lack UBDT branch, so we need to account for this later).

    1. Clone this project on `lxplus`
    2. In the project root, run `make build-rdx-true-to-tag-2016-lxplus`. The output files will be saved in `gen/rdx-<timestamp>-true_to_tag_lxplus-2016`.
    3. Copy the output to `histos/default`.

    Note that you might need to use a specific version of `pidcalib`. Edit `cmd_prefix` in `scripts/pidcalib_wrapper.py` and replace it with:

    ```shell
    cmd_prefix = "lb-conda pidcalib/2022-09-02 " if mode == "lxplus" else ""
    ```

    `pidcalib/2022-09-02` corresponds to version `1.0.6`.

4. Generate ghost efficiencies and $e$ conditional efficiencies with local incl. $J/\psi$ MC samples. These $e$ conditional efficiencies will be used to account for the missing UBDT cuts in step 2.

    1. Clone this project on `glacier` (in this step, running on `glacier` is not a requirement)
    2. Run `nix develop` in the project root
    3. In the resulting shell, run `make build-rdx-true-to-tag-2016-local`. The output files will be saved in `gen/rdx-<timestamp>-true_to_tag_local-2016`.
    4. Copy the output to `histos/default`.

Starting from the next step, all operations are done locally or on `glacier` and `nix develop` is always assumed.

5. Merge all efficiencies obtained above into a single file locally with `make build-rdx-merged-2016`.

    1. If not done already, copy the output folders from the previous steps from `gen` to `histos/default` folder.
    2. The inputs to this step are listed in the `input_histos/2016/true_to_tag/default` section of `rdx-run2.yml`, so
    make sure to update it before you run the `make` command!
    3. This stage will produce `gen/rdx-<timestamp>-merged-2016/merged.root` as output. Copy it to `histos/default` as well.
    4. Update the `EFFICIENCIES` variable in `Makefile`.

> [!IMPORTANT]
> The efficiencies produced in the steps above assume the nominal UBDT cut `mu_uBDT > 0.25`.
> For the misid validation sample (which requires `mu_uBDT < 0.25` instead) you should repeat the steps above (**except step 2!**) setting the ```USE_CTRL_SAMPLE``` variable appropriately, e.g.:
>
> ```shell
> make build-rdx-true-to-tag-2016-glacier USE_CTRL_SAMPLE=true
> ```
>
> Assigning unexpected values to this variable will generate a warning and will cause the nominal uBDT cut to be used.
>
> For step 2, instead of `make build-rdx-misid-mc-corrections` you should run `make build-rdx-misid-mc-corrections-vmu`.
> The `-fake_mu` command is not affected by the change in the uBDT cut, so you don't need to run it again.
>
> Furthermore:
>
> 1. The outputs of steps 1 through 5 should be placed in `histos/ctrl_sample`.
> 2. After step 4, you should update `input_histos/2016/true_to_tag/misid_ctrl` in `rdx-run2.yml`
> 3. After step 5, you should update the `EFFICIENCIES_VMU` variable in `Makefile`.
>

6. Extract $K$ and $\pi$ momentum smearing info into a small ntuple with `make build-generic-dif-smearing`.

    1. This stage will produce `gen/generic-<timestamp>-dif_smearing/dif.root` as output. Copy it to `histos` (not `histos/default` or `histos/ctrl_sample`! Those are reserved for nominal- or vmu-specific ntuples).
    2. Update the `DiF` variable in `Makefile`.

#### Remark on negative efficiencies

Sometimes the efficiencies are negative, due to sWeight. To correct for that,
we followed a procedure suggested by Phoebe. For more info, see p. 3 of [this slide](https://github.com/umd-lhcb/group-talks/blob/master/phys_group_meetings/22-04-13_yipeng_rdx_status.pdf).

#### Remark on `cut` vs. `pid_cut` options in the YAML file

These names are consistent with `pidcalib2` naming scheme, that is:

- `cut`: The cuts in both the denominator and numerator of the efficiency
- `pid_cut`: The _additional_ cuts applied only to the numerator

Therefore, the efficiency is defined as:

$$
\epsilon = \frac{\text{event passing cut and pid-cut}}{\text{event passing cut}}
$$

The cuts and pid_cuts are defined in `rdx-run2.yml` as described below.

- For steps 1 and 2:
  - In the case of non-muon tags (e.g. $K$ true to $\pi$ tag), the pid_cuts for each tag are defined in `tags`->`tag`, while the cuts are defined in `pidcalib_config`->`tags`->`cut`.
  - For the case of muon tags, both cuts and pid_cuts are defined in `pidcalib_config`->`tags_addon`->`nom`/`denom`.
    Here, `nom` refers to the muon tag efficiency that appear in the numerator of the transfer factor, and it amounts to the efficiency to pass the nominal mu PID requirements.
    On the other hand, `denom` refers to the muon tag efficiency that appears in the denominator of the transfer factor, which is the efficiency to pass the mu PID requirements of the fake muon sample.
    - Note that `nom` also has an `add_pid_cut` entry; The requirement listed here is appended to `pid_cut`.

- For step 3, the cuts are listed in `local_pid_config`, however the structure is similar:
  - For non-muon tags (in this step, only $g$), the pid_cuts are defined in `tags`->`tag`, while the cuts are defined in `local_pid_config`->`2016`->`g`->`tags`->`cut`.
  - For the muon tags, both cuts and pid_cuts are defined in `local_pid_config`->`2016`->`g`->`tags_addon`->`nom`/`denom` and `local_pid_config`->`2016`->`e`->`tags_addon`->`rel`.
    Recall that, in this step, for $e$ we only calculate corrections to account for the missing uBDT information in step 2, hence the single `rel` correction for the `nom` efficiency computed in step 2.

#### Remark on test run

One can test these rules with:

- `make test-pidcalib2-wrapper-glacier`
- `make test-pidcalib2-wrapper-lxplus`

locally, without actually running the efficiency generation program.

### Generation of tagged yields

This step relies on our misID ntuple with a special process which add
PIDCalib-like branches like `Brunel_ProbNNghost` (which is just an alias
to `mu_ProbNNghost`).

The special process is needed so that we can apply the same cuts on
these ntuples as if we are applying cuts on PIDCalib ntuples, eliminating
the need to translate branch names.

To build these special ntuples, go to `lhcb-ntuples-gen`, then type:

```shell
make rdx-ntuple-run2-misid_study
```

Once that is done, configure `spec/rdx-run2.yml` properly, then:

```shell
make build-rdx-tag-2016
```

This will produce `gen/generic-<timestamp>-tag/tagged.root` as output, which you should also copy to `histos`.
Finally, update the `TAGGED` variable in `Makefile`.

### Unfold the misID efficiencies

Make sure to download all required ntuples with `git annex get`!
Then, make sure to set the inputs to this step correctly. The unfolding procedure requires the `merged.root` and `tagged.root` produced above. Go to `Makefile` and set the `EFFICIENCIES`, `EFFICIENCIES_VMU` and `TAGGED` variables to point to your `merged.root` and `tagged.root` files, respectively, as explained above.

```shell
make build-rdx-unfolded-2016
```

This will produce `gen/generic-<timestamp>-unfoled/unfolded.root` as output, which you should (you guessed!) copy to `histos`.
Once copied, update the `UNFOLDED` variable in `Makefile`.

One can also test-run the unfolding (w/o actually doing unfolding!)
by `make test-unfold`.

> [!IMPORTANT]
> To proceed with the analysis workflow, `unfolded.root` and `dif.root` must be copied over to _lhcb-ntuples-gen_ (into `/run2-rdx/reweight/misid/histos/`).
> The unfolded efficiencies and momentum smearing histograms will be used to produce the misid step 2 ntuples for the fit.

### Apply misID and DiF smearing weights

```shell
make build-rdx-weights-2016
```

This compiles and runs `src/ApplyMisIDWeight.cpp` to apply misID weights and momentum smearing on misID control samples.

**Note**: This rule here is for demo purposes only. This script is actually used in _lhcb-ntuples-gen_ to calculate and apply the transfer factors for the fake mu sample.
Notably, the factor 10 due to the **prescale** on the fake mu HLT2 line is hardcoded in this C++ script.
