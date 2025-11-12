# scripts directory

Contains python scripts needed for the misid undolfing procedure.


## compare_effs.py

Usage: From the root directory, run `./scripts/compare_effs.py <path_to_old_effs> <path_to_new_effs>`, e.g.

```shell
~/misid-unfold$ ./scripts/compare_effs.py histos/default/rdx-24_09_10_10_16-true_to_tag_glacier-2016 histos/default/rdx-25_04_08_06_47-true_to_tag_glacier-2016
```

The script scans the "new" directory and compares each efficiency histogram with the file of same name store in the "old" directory.
Files not present in both directories are ignored. It also is expected that the efficiency histograms in the root files are named "eff".
