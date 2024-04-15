# phosx-benchmark

Benchmark of the kinase activity inference method [PhosX](https://github.com/alussana/phosx).

> Research paper: [https://doi.org/10.1101/2024.03.22.586304](https://doi.org/10.1101/2024.03.22.586304)

## Build the container image

Requires [Docker](https://www.docker.com) and [Apptainer](https://apptainer.org).

```bash
docker build -t phosx-benchmark - < env/Dockerfile
docker save -o env/phosx-benchmark.tar.gz phosx-benchmark
singularity build env/phosx-benchmark.sif docker-archive://env/phosx-benchmark.tar.gz
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

![flowchart](misc/flowchart.svg)

## Cite

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