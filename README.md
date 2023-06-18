[![Rust Build](https://github.com/marschall-lab/panacus/actions/workflows/rust_build.yml/badge.svg)](https://github.com/marschall-lab/panacus/actions/workflows/rust_build.yml) [![Anaconda-Server Badge](https://anaconda.org/bioconda/panacus/badges/version.svg)](https://conda.anaconda.org/bioconda) [![Anaconda-Server Badge](https://anaconda.org/bioconda/panacus/badges/platforms.svg)](https://anaconda.org/bioconda/panacus) [![Anaconda-Server Badge](https://anaconda.org/bioconda/panacus/badges/license.svg)](https://anaconda.org/bioconda/panacus)

# A Counting Tool for Pangenome Graphs

![panacus is a counting tool for pangenome graphs](docs/panacus-illustration.png?raw=true "panacus is a counting tool for pangenome graphs")

`panacus` is a tool for calculating statistics for [GFA](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) files. It supports GFA files with `P` and
`W` lines, but requires that the graph is `blunt`, i.e., nodes do not overlap and consequently, each link (`L`) points from the end of one segment
(`S`) to the start of another.

`panacus` supports the following calculations:

- coverage histogram
- pangenome growth statistics
- path-/group-resolved coverage table

## Dependencies

`panacus` is written in [RUST](https://www.rust-lang.org/) and requires a working RUST build system for installation. See [here](https://www.rust-lang.org/tools/install) for more details.

- clap
- itertools
- quick-csv
- rand
- rayon
- regex
- rustc-hash
- strum, strum_macros

`panacus` provides a Python script for visualizing the calculated counting statistics and requires the following Python libraries:

- matplotlib
- numpy
- pandas
- scikit-learn
- scipy
- seaborn

## Installation

### From bioconda channel

Make sure you have [conda](https://conda.io)/[mamba](https://anaconda.org/conda-forge/mamba) installed!

```shell
mamba install -c conda-forge -c bioconda panacus
```

### From binary release 
#### Linux x86\_64
```shell
wget --no-check-certificate -c https://github.com/marschall-lab/panacus/releases/download/0.2.1/panacus-0.2.1_linux_x86_64.tar.gz
tar -xzvf panacus-0.2.1_linux_x86_64.tar.gz

# suggestion: add tool to path in your ~/.bashrc
export PATH="$(readlink -f panacus-0.2.1_linux_x86_64/bin)":$PATH

# you are ready to go! 
panacus --help
```

#### Mac OSX arm64
```shell
wget --no-check-certificate -c https://github.com/marschall-lab/panacus/releases/download/0.2.1/panacus-0.2.1_macos_arm64.tar.gz
tar -xzvf panacus-0.2.1_macos_arm64.tar.gz

# suggestion: add tool to path in your ~/.bashrc
export PATH="$(readlink -f panacus-0.2.1_macos_arm64/bin)":$PATH

# you are ready to go! 
panacus --help
```

### From repository
```shell
git clone git@github.com:marschall-lab/panacus.git

cd panacus
cargo build --release

mkdir bin
ln -s ../target/release/panacus bin/
ln -s ../scripts/panacus-visualize.py bin/panacus-visualize

# suggestion: add tool to path in your ~/.bashrc
export PATH="$(readlink -f bin)":$PATH

# you are ready to go! 
panacus --help

```

## Run

```console
$ panacus
Calculate count statistics for pangenomic data

Usage: panacus <COMMAND>

Commands:
  histgrowth          Run in default mode, i.e., run hist and growth successively and output the results of the latter
  hist                Calculate coverage histogram from GFA file
  growth              Construct growth table from coverage histogram
  ordered-histgrowth  Compute growth table for order specified in grouping file (or, if non specified, the order of paths in the GFA file)
  table               Compute coverage table for count items
  help                Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

## Pangenome Growth Statistics

Here's a quick example for computing pangenome growth statistics on the HPRC v.1.0 pggb, chr 22: 

1. Download and unpack the graph:
```shell
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr22.hprc-v1.0-pggb.gfa.gz
gunzip chr22.hprc-v1.0-pggb.gfa.gz
```
2. Prepare file to select subset of paths corresponding to haplotypes:
```shell
grep '^P' chr22.hprc-v1.0-pggb.gfa | cut -f2 | grep -ve 'grch38\|chm13' > chr22.hprc-v1.0-pggb.paths.haplotypes.txt
```
3. Run `panacus histgrowth` to calculate pangenome growth for nodes (default) with coverage/quorum thresholds 1/0, 2/0, 1/1, 1/0.5, and 1/0.1 using up to 4 threads:
```shell
RUST_LOG=info panacus histgrowth -t4 -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -s chr22.hprc-v1.0-pggb.paths.haplotypes.txt chr22.hprc-v1.0-pggb.gfa > chr22.hprc-v1.0-pggb.histgrowth.node.tsv
```
4. Visualize growth curve and estimate growth parameters:
```shell
panacus-visualize -e chr22.hprc-v1.0-pggb.histgrowth.node.tsv > chr22.hprc-v1.0-pggb.histgrowth.node.pdf
```

![pangenome growth of nodes in hprc-v1.0-pggb.gfa](docs/chr22.hprc-v1.0-pggb.histgrowth.node.png?raw=true "pangenome growth statistics on the HPRC v.1.0 pggb, chr 22")

## Ordered Pangenome Growth Statistics

Sometimes it is interesting to look at the pangenome growth when samples are processed in a specific order rather than considering all all possible
orders. `panacus`' capability to construct such plots is illustrated here by the example of the GRCh38-based HPRC v.1.0 minigraph-cactus graph (all
chromosomes). The example reproduces Figure 3g(left) from the publication [A draft human pangenome
reference](https://doi.org/10.1038/s41586-023-05896-x) that quantifies pangenome growth of the amount of non-reference (GRCh38) sequence of the
minigraph-cactus based human pangenome reference graph.

1. Download and unpack the graph:
```shell
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.gfa.gz
gunzip hprc-v1.0-mc-grch38.gfa.gz
```
2. Establish order of samples in the growth statistics:
```shell
echo 'HG03492 HG00438 HG00621 HG00673 HG02080 HG00733 HG00735 HG00741 HG01071 HG01106 HG01109 HG01123 HG01175 HG01243 HG01258 HG01358 HG01361 HG01928
HG01952 HG01978 HG02148 HG01891 HG02055 HG02109 HG02145 HG02257 HG02486 HG02559 HG02572 HG02622 HG02630 HG02717 HG02723 HG02818 HG02886 HG03098
HG03453 HG03486 HG03516 HG03540 HG03579 NA18906 NA20129 NA21309' | tr ' ' '\n' > hprc-v1.0-mc-grch38.order.samples.txt
```
3. Exclude paths from reference genome GRCh38
```shell
grep '^P' hprc-v1.0-mc-grch38.gfa | cut -f2 | grep -ie 'grch38' > hprc-v1.0-mc-grch38.paths.grch38.txt
```
4. Run `panacus ordered-histgrowth` to calculate pangenome growth for basepairs with coverage thresholds 1,2,3, and 42 using up to 4 threads:
```shell
RUST_LOG=info panacus ordered-histgrowth -c bp -t4 -l 1,2,3,42 -S -e hprc-v1.0-mc-grch38.paths.grch38.txt hprc-v1.0-mc-grch38.gfa > hprc-v1.0-mc-grch38.ordered-histgrowth.bp.tsv
```
(The log will report some errors regarding missing order information of CHM13 paths. These paths will be ignored in the plot, which is the intended
behavior of this command line call)

5. Visualize growth curve and estimate growth parameters :
```shell
panacus-visualize hprc-v1.0-mc-grch38.ordered-histgrowth.bp.tsv > hprc-v1.0-mc-grch38.ordered-histgrowth.bp.pdf
```

![ordered pangenome growth of bps in hprc-v1.0-mc-grch38.gfa](docs/hprc-v1.0-mc-grch38.ordered-histgrowth.bp.png?raw=true "pangenome growth of non-reference sequence in the HPRC v.1.0 MC GRCh38 graph")

