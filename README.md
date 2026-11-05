# Monitoring Dynamics in Ecosystem Service Indicators with Remote Sensing

## Companion Resources

This repository contains companion resources for the book:

**Pedro Cabral, _Monitoring Dynamics in Ecosystem Service Indicators with Remote Sensing_, Elsevier.**

The resources are designed to support teaching, learning, and practical exploration of selected concepts and methodological frameworks presented in the book. They include teaching exercises, code, example datasets, notebooks, analytical workflows, and supporting documentation.

The exercises are intended to help readers move from environmental observations to interpretable ecosystem service indicators and decision-relevant information. Topics covered by the companion materials may include ecosystem dynamics, scenario-based projection, uncertainty and confidence assessment, spatial analysis, remote sensing, and environmental decision support.

## Exercises

| # | Exercise | What it covers | Book chapters |
|---|---|---|---|
| 1 | [Carbon storage change in mainland Portugal, 1990–2018](exercises/portugal_carbon_1990_2018/) | CORINE Land Cover preparation, InVEST carbon model, transition attribution, valuation ranges, ensemble confidence bands, MODIS NDVI trends, fire disturbance and recovery, SDG 15.3.1 productivity trajectory | 1–6, 8 |

[`docs/chapter_map.md`](docs/chapter_map.md) lists, chapter by chapter, where each part of the book is put into practice.

**Quick start**

```bash
conda env create -f environment.yml
conda activate es-rs-book
jupyter lab exercises/portugal_carbon_1990_2018
```

## Resources

Materials will be added and updated as they are developed and used in graduate teaching and professional training.

The repository will be organized into:

- `exercises/` — self-contained teaching exercises
- `notebooks/` — executable examples and demonstrations
- `scripts/` — reusable analytical code
- `data/` — example or derived datasets where redistribution is permitted
- `docs/` — supporting documentation

## Archived and citable version

Citable versions of the companion resources are archived on Figshare:

**Figshare DOI:** https://doi.org/10.6084/m9.figshare.33439819

Please cite the Figshare record together with the book when directly using or adapting these companion materials.

## Software and data

The resources are designed to use open-source software and openly accessible data wherever practicable. Individual exercises provide information about their software requirements, data sources, and licences.

Third-party datasets are not redistributed unless their licences explicitly permit redistribution. Where appropriate, instructions are provided for obtaining data from the original source.

## License

Code in this repository is released under the MIT License unless otherwise stated. Documentation, teaching materials, and datasets may be subject to separate licences specified within the relevant directories or files.

## Author

**Pedro Cabral**  
NOVA Information Management School (NOVA IMS)  
NOVA University Lisbon, Portugal
