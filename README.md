# TICC

*The Tomography Inversion Code Collection*

**ATTENTION**:
This is a compilation of various codes written in the past (from Princeton, Munich...).
All codes are in FFINVERSION on GitHub (https://github.com/kasra-hosseini/FFINVERSION).

This folder contains codes that one needs for the finite frequency inversion process and some visualisation codes for statistics or creating vtk files.

**Not** in this folder is **obspyDMT**, or a code of your choice to download/organise a database for processing.
**Not** in this folder are **Karins STF codes**.

N.B.! While installing all the packages remember to:
1) create a new environment with python 3.8 installed
2) start with main package instaseis: conda install instaseis
3) start downloading all the other

Order of running the codes:

* [get database, eg. download data with obspyDMT]
    ### example command used to download data in Europe:
    obspyDMT --datapath /mnt/data/eu-data --mag_type=Mw --min_mag 5.0 --min_date 2010-01-01 --max_date 2011-01-01 --station_rect=-50/60/20/80 --data_source=all --cut_time_phase P --preset=600 --offset=900 --waveform_format=sac --cha=BHZ,HHZ,SHZ,EHZ --req_parallel —req_np=8 --event_catalog=NEIC_USGS --instrument_correction --parallel_process --process_np=6 --pre_process process_unit_default

* [genearte STF's from database]

* **pyffproc** (cross correlation with synthetic waveforms)

* **meshgrid** (creates grid)

* **raydata_raymatrix** (kernels, common corrections, etc)

* **assemble_matrix** (phases you wish to invert for)

* **solvetomo** (solve inversion)

Details about how to run the codes (e.g. input, output, etc) are in separate README files in the corresponding folder.
