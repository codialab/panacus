[![Rust Build](https://github.com/codialab/panacus/actions/workflows/rust_build.yml/badge.svg)](https://github.com/codialab/panacus/actions/workflows/rust_build.yml) [![Anaconda-Server Badge](https://anaconda.org/bioconda/panacus/badges/version.svg)](https://anaconda.org/bioconda/panacus) [![Anaconda-Server Badge](https://anaconda.org/bioconda/panacus/badges/platforms.svg)](https://anaconda.org/bioconda/panacus) [![Anaconda-Server Badge](https://anaconda.org/bioconda/panacus/badges/license.svg)](https://anaconda.org/bioconda/panacus)

# A Counting Tool for Pangenome Graphs

![panacus is a counting tool for pangenome graphs](docs/panacus-illustration.png?raw=true "panacus is a counting tool for pangenome graphs")

`panacus` is a tool for calculating statistics for pangenome files. It supports GFA files with `P` and
`W` lines, but requires that the graph is `blunt`, i.e., nodes do not overlap and consequently, each link (`L`) points from the end of one segment
(`S`) to the start of another. These GFA files can be sqz-ified (see [sqz](https://github.com/codialab/sqz)). Additionally, it also supports
VCF files.

`panacus` supports the following calculations:

- coverage histogram
- pangenome growth statistics
- path similarity
- allele/non-reference features-plots
- node plots resolved by length and coverage
- ...

## Quickstart
1. Install panacus using conda/mamba:
```shell
mamba install -c conda-forge -c bioconda panacus
```
2. Create a file `report.yaml` with the following content:
```yaml
- !Gfa
  graph: ../graphs/test_graph.gfa    # Change this to a GFA file on your system
  count_type: Bp
  analyses:
    - !Hist
    - !Growth
      coverage: 1,1,2
      quorum: 0,0.9,0
```
3. Run panacus:
```shell
panacus report report.yaml > report.html
```
4. Take a look at the generated html file using your favorite browser!

For more info on what to write into `report.yaml` see the [documentation](https://github.com/codialab/panacus/wiki).

## Installation
### From bioconda channel

Make sure you have [conda](https://conda.io)/[mamba](https://anaconda.org/conda-forge/mamba) installed!

```shell
mamba install -c conda-forge -c bioconda panacus
```

### From binary release
#### Linux x86\_64
```shell
wget --no-check-certificate -c https://github.com/codialab/panacus/releases/download/v0.5.3/panacus-0.5.3_x86_64-unknown-linux-musl.tar.gz
tar -xzvf panacus-0.5.3_x86_64-unknown-linux-musl.tar.gz

# install the Python libraries necessary for panacus-visualize
pip install --user matplotlib numpy pandas scikit-learn scipy seaborn

# suggestion: add tool to path in your ~/.bashrc
export PATH="$(readlink -f panacus-0.5.3_x86_64-unknown-linux-musl/bin)":$PATH

# you are ready to go!
panacus --help
```

#### Mac OSX arm64
```shell
wget --no-check-certificate -c https://github.com/codialab/panacus/releases/download/v0.5.3/panacus-0.5.3_aarch64-apple-darwin.tar.gz
tar -xzvf panacus-0.5.3_aarch64-apple-darwin.tar.gz

# install the Python libraries necessary for panacus-visualize
pip install --user matplotlib numpy pandas scikit-learn scipy seaborn

# suggestion: add tool to path in your ~/.bashrc
export PATH="$(readlink -f panacus-0.5.3_aarch64-apple-darwin/bin)":$PATH

# you are ready to go!
panacus --help
```

### From repository
`panacus` requires a working RUST build system (version >= 1.74.1) to build from source. See [here](https://www.rust-lang.org/tools/install) for more details.
```shell
git clone git@github.com:codialab/panacus.git

cd panacus
cargo build --release

mkdir bin
ln -s ../target/release/panacus bin/
ln -s ../scripts/panacus-visualize.py bin/panacus-visualize

# install the Python libraries necessary for panacus-visualize
pip install --user matplotlib numpy pandas scikit-learn scipy seaborn

# suggestion: add tool to path in your ~/.bashrc
export PATH="$(readlink -f bin)":$PATH

# you are ready to go!
panacus --help
```

## Examples
Examples can be found in the [examples directory](/examples/). To get more information about how to use `panacus` check out the [documentation](https://github.com/codialab/panacus/wiki).

## Citation
Parmigiani, L., Garrison, E., Stoye, J., Marschall, T. & Doerr, D. Panacus: fast and exact pangenome growth and core size estimation. Bioinformatics, https://doi.org/10.1093/bioinformatics/btae720 (2024).
