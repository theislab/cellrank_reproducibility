# How to run the memory benchmarks

From the current directory, simply run:
```bash
bash run_all.sh <METHOD> <NUMBER_OF_JOBS>
```

Where `<MEHOD>` is one of: gpcca, palantir, velocyto, fateid, stemnet.
Note that `<NUMBER_OF_JOBS>` only affects CellRank and Palantir.
It was set to either `32` for multi-core benchmark or `1` for single-core.

The respective commands for multi-core benchmarks were:
```bash
bash run_all.sh gpcca 32
bash run_all.sh palantir 32
bash run_all.sh velocyto 32
bash run_all.sh stemnet
bash run_all.sh fateid
```

The respective commands for single-core benchmarks on 100k cells were:
```bash
bash run_all.sh gpcca 1
bash run_all.sh palantir 1
```
