# Median Shift Lunar Correction for VIIRS

This repository is supplementary to the manuscript of "Median Shift Lunar Correction for VIIRS", submitted to IEEE Geoscience and Remote Sensing Letters.

## Details

### Python scripts

In Python directory is contained  the *vnp46A1_to_geotiff.py* file which:

* converts hdf files to geotiffs.

* geotiffs are returned both in the initial radiation values and lunar corrected values, reprojected in Greek Grid and resampled in 500m. resolution.

* the required flags are applied to keep useful information only.



The execution is as follows:
*vnp46A1_to_geotiff.py* inDir outDir BRDF_Dir exportOriginalDNB

* inDir  = directory which contains VNP46A1 hdf files.

* outDir = directory in which geotiffs are saved.

* BRDF_Dir = directory in which BRDF VNP43MA4 hdf files are contained. (AOI must be covered by the extents of BRDF file for a successful lunar correction)

* exportOriginalDNB=True/False to export or not geotiffs with original radiance values.

In order to correct the lunar radiation are used the MT 2009 model<sup>[1](#myfootnote1)</sup> and the relevant process as it is described in bibliography<sup>[2](#myfootnote2)</sup>. The relevant variables/settings for the extent, the reprojection and resampling of data are hardcoded.

For the needs of the current research, all data that have been scanned before 20:00 UTC time (usually in the early hours of the current date) are exported to another geotiff file separately with a previous date flag.

The vnp46A1_to_geotiff.py file is executed first in order.

### R scripts

*geotiffs2grid.R*: generates daily composite raster stacks (for original and lunar corrected geotiffs) and save them as R grid filetype.
Important parameter is the variable  *DNB_geotiffs_DIR* which indicates the directory with the geotiffs resulting from the execution of *vnp46A1_to_geotiff.py* script.

The *MainScript.R* file 

*MainScript.R* is the main file which executed per AOI and includes the basic code for the proposed Median Shift method, the validation metrics and the accompanying plots. 

Important parameter is the variable "name" which defines the AOI  on which the script will be executed. The value of *name*  
must be registered in the column *NAME* inside  *data/extents.csv*.


### Ancillary files

*R/config.yml*, configuration file. Contains general settings and settings per AOI.

*data/extents.csv*, extents of each AOI.

*R/miller_plot.R*, generates Figure 1.

*R/map_with_extents.R*, generates Figure 2.

*R/myfunctions.R*, contains ancillary functions.


## Authors

* [Stathakis Dimitris](https://gr.linkedin.com/in/dstath)

* [Liakos Leonidas](https://gr.linkedin.com/in/leonidasliakos)


## Acknowledgements
This research is co-financed by Greece and the European Union (European Social Fund- ESF) through the Operational Programme «Human Resources Development, Education and Lifelong Learning» in the context of the project "Strengthening Human Resources Research Potential via Doctorate Research" (MIS-5000432), implemented by the State Scholarships Foundation (ΙΚΥ, grand #2018-050-0502-13146). 


## References

<a name="myfootnote1"><sup>[1]</sup></a> S. D. Miller and R. E. Turner, "A Dynamic Lunar Spectral Irradiance
Data Set for NPOESS/VIIRS Day/Night Band Nighttime Environmental
Applications," IEEE Transactions on Geoscience and Remote Sensing,
vol. 47, no. 7, pp. 2316-2329, July 2009.

<a name="myfootnote2"><sup>[2]</sup></a> M.O. Román, Z. Wang, Q. Sun, V. Kalb, S. Miller, A. Molthan, et al.,
"NASA’s Black Marble nighttime lights product suite", Remote
Sens. Environ., vol. 210, pp. 113–143, 1 June 2018.

