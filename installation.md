# Installation of SAM-harmonization pipeline

```bash
# setup conda env
conda env create -f environment.yml
```

To add the conda packages to an existing environment:

1. Comment out (`#`) the first line of the `environment.yml` containing the `name: sam-harmonization`
2. Activate the desiered conda env
3. Run `conda env update --file environment.yml`

Installing the package:

```bash
# sam-harmonization
git clone https://git.wur.nl/brank001/sam-harmonization
cd src
perl Makefile.PL
make
make install
```

Alternatively, you can use the following docker image:
[docker-registry.wur.nl/brank001/sam-harmonization](https://git.wur.nl/brank001/sam-harmonization/container_registry)

```bash
docker pull docker-registry.wur.nl/brank001/sam-harmonization
```
