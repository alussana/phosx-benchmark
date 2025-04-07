# phosx-benchmark

Benchmark of the kinase activity inference method [PhosX](https://github.com/alussana/phosx).

> Research paper: [https://doi.org/10.1101/2024.03.22.586304](https://doi.org/10.1101/2024.03.22.586304)

## Build the container image

Requires [Docker](https://www.docker.com) and [Apptainer](https://apptainer.org).

```bash
docker build -t phosx-benchmark - < env/Dockerfile
docker save -o env/phosx-benchmark.tar.gz phosx-benchmark
apptainer build env/phosx-benchmark.sif docker-archive://env/phosx-benchmark.tar.gz
```

## Customise `nextflow.config`

Modify `process.executor`, `process.queue`, `workDir`, and `env.out_dir` according to the infrastructure where the workflow will be executed. Find the Nextflow [configuration file](https://www.nextflow.io/docs/latest/config.html) documentation.

Alternatively, as a minimal example to run the workflow locally, just replace the `nextflow.config` with `misc/nextflow-local.config` (a backup `nextflow.config~` will be created):

```bash
mv misc/nextflow-local.config nextflow.config -b
```

Finally set `params.n_cores` to the preferred number of cores that the workflow should use.

## Run the workflow

```bash
nextflow run main.nf -resume -with-dag misc/flowchart.svg
```

## Cite

### BioRxiv 

BibTeX:

```bibtex
@article{Lussana2024,
  title = {PhosX: data-driven kinase activity inference from phosphoproteomics experiments},
  url = {http://dx.doi.org/10.1101/2024.03.22.586304},
  DOI = {10.1101/2024.03.22.586304},
  publisher = {Cold Spring Harbor Laboratory},
  author = {Lussana,  Alessandro and Petsalaki,  Evangelia},
  year = {2024},
  month = mar 
}
```

### Bioinformatics

BibTeX:

```bibtex
@article{10.1093/bioinformatics/btae697,
    author = {Lussana, Alessandro and Müller-Dott, Sophia and Saez-Rodriguez, Julio and Petsalaki, Evangelia},
    title = {PhosX: data-driven kinase activity inference from phosphoproteomics experiments},
    journal = {Bioinformatics},
    volume = {40},
    number = {12},
    pages = {btae697},
    year = {2024},
    month = {11},
    abstract = {The inference of kinase activity from phosphoproteomics data can point to causal mechanisms driving signalling processes and potential drug targets. Identifying the kinases whose change in activity explains the observed phosphorylation profiles, however, remains challenging, and constrained by the manually curated knowledge of kinase–substrate associations. Recently, experimentally determined substrate sequence specificities of human kinases have become available, but robust methods to exploit this new data for kinase activity inference are still missing. We present PhosX, a method to estimate differential kinase activity from phosphoproteomics data that combines state-of-the-art statistics in enrichment analysis with kinases’ substrate sequence specificity information. Using a large phosphoproteomics dataset with known differentially regulated kinases we show that our method identifies upregulated and downregulated kinases by only relying on the input phosphopeptides’ sequences and intensity changes. We find that PhosX outperforms the currently available approach for the same task, and performs better or similarly to state-of-the-art methods that rely on previously known kinase–substrate associations. We therefore recommend its use for data-driven kinase activity inference.PhosX is implemented in Python, open-source under the Apache-2.0 licence, and distributed on the Python Package Index. The code is available on GitHub (https://github.com/alussana/phosx).},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btae697},
    url = {https://doi.org/10.1093/bioinformatics/btae697},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/40/12/btae697/60972735/btae697.pdf},
}
```
