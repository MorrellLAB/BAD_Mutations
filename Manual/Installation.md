# Installation Guide
To install `BAD_Mutations`, simply clone the repository:

```bash
git clone https://github.com/MorrellLAB/BAD_Mutations.git
```

## Dependencies
### With Conda
`BAD_Mutations` comes with a [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) `.yml` file that will install the dependencies. This file is located in the `Supporting/bad_mutations.yml` file under the repository. Build the environment as follows:

```bash
conda env create -f /path/to/BAD_Mutations/Supporting/bad_mutations.yml
```

Then, activate the environment:

```bash
conda activate bad_mutations
```

#### A Note on PASTA
The `PASTA` aligner requires that you set a path to the `MAFFT` aligner. This can be achieved by setting an environment variable. You can automatically set the correct path by running the following commands, once you have activated the `bad_mutations` Conda environment:

```bash
mkdir -p \
    ${CONDA_PREFIX}/etc/conda/activate.d \
    ${CONDA_PREFIX}/etc/conda/deactivate.d

touch \
    ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh \
    ${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh
```

Edit `${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh` to look like the following:

```bash
#!/bin/sh
export PASTA_TOOLS_DEVDIR="${CONDA_PREFIX}/bin"
```

Edit `${CONDA_PREFIX}/etc/conda/deactivate.d/env_vars.sh` to look like the following:

```bash
#!/bin/sh
unset PASTA_TOOLS_DEVDIR
```

Deactivate and re-activate the Conda environment set the environment variable. You should not be able to run `BAD_Mutations` with all its dependencies.

### Without Conda
Withou Conda, you must install the following, and make sure they are in your `$PATH`:

- Python 3.x
- [PASTA](https://github.com/smirarab/pasta)
- [HyPhy](http://hyphy.org/resources/)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [BioPython](https://biopython.org/)

If you are on MacOS, then the [Homebrew](https://brew.sh/) package manager should have build scripts for these tools.
