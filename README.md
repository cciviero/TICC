# TICC

*The Inversion Code Collection*

**ATTENTION**:
This is a compilation of various codes written in the past (from Princton, Munich...).
All codes are in FFINVERSION on GitHub (https://github.com/kasra-hosseini/FFINVERSION).

This folder contains codes that one needs for the finite frequency inversion process and some visualisation codes for statistics or creating vtk files.

**Not** in this folder is **obspyDMT**, or a code of your choice to download/organise a database for processing.
**Not** in this folder are **Karins STF codes**.

Order of running the codes:

* [get database, eg. download data with obspyDMT]
* [genearte STF's from database]

* **pyffproc** (cross correlation with synthetic waveforms)

* **meshgrid** (creates grid)

* **raydata_raymatrix** (kernels, common corrections, etc)

* **assemble_matrix** (phases you wish to invert for)

* **solvetomo** (solve inversion)

Details about how to run the codes (e.g. input, output, etc) are in separate README files in the corresponding folder.