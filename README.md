# Copernicus Climate Change - Early Upper Air Service (CEUAS)

> Service Contract C3S 311c Lot2

<img src="https://climate.copernicus.eu/sites/default/files/custom-uploads/ESOTC2020/press_resources/Copernicus%20vecto%20def%20%20Europe's%20eyes%20on%20Earth%20with%20flag.png" width="300px"> <img src="https://climate.copernicus.eu/sites/default/files/custom-uploads/ESOTC2020/press_resources/C3S%E2%80%93POS%E2%80%93LINE.png" width="300px">
<img src="https://climate.copernicus.eu/sites/default/files/custom-uploads/ESOTC2020/press_resources/ECMWF_LOGO_POS_ImplementedBy.png" width="300px"> <img src="https://communications.univie.ac.at/fileadmin/user_upload/d_oeffentlichkeitsarbeit/Logos/2016-02/Uni_Logo_2016.png" width="300px"> <img src="https://img.univie.ac.at/fileadmin/_processed_/csm_logo_imgw_color_with_text_2100x660_8c1ff1816b.png" width="300px">

Collaborators:
* L. Haimberger
* M. Blaschek
* F. Ambrogi
* U. Voggenberger

Status: Development

Date: 06-2021


Host: [UNIVIE](https://www.univie.ac.at)

## Copernicus Climate Change Service Volume II

**Support for climate reanalysis including satellite data rescue**

*Lot 2: Historic upper-air data. The objective for this Lot is to develop and maintain a quality-controlled
global database containing all known digitised in-situ upper-air weather observations made prior to
1979, together with metadata and information needed for data assimilation such as bias adjustments
and uncertainty estimates.*

### Common Data Model (CDM)

More information can be found in [`CEUAS/cdm`](./CEUAS/cdm). Within this task the CDM tables are downloaded and converted.

### Metadata

More information can be found in [`CEUAS/meta`](./CEUAS/meta). Within this task the station inventory is generated from multiple sources and a consolidated list is produced to be used in processing and merging stations.

### CDS Backend

More information can be found in [`CEUAS/public/cds-backend`](./CEUAS/public/cds-backend). Within this task the backend server and handling of the requests are coded. When requesting data from the Climate Data Store (CDS) the request is forwarded to the backend server. The request is then converted into subtasks and the data files are processed and results are returned. An example on how to use can be found here: [Example.ipynb](./CEUAS/public/cds-backend/Example.ipynb)

### Harvest

More information can be found in [`CEUAS/public/harvest`](./CEUAS/public/harvest). Within this task the source datasets are downloaded from different providers and convert to a CDM compliant netCDF format.

### Merge

More information can be found in [`CEUAS/public/merge`](./CEUAS/public/merge). Within this task the CDM compliant netCDF files from the harvest task are merged and a consistent and combined dataset is produced.

### Adjust

More information can be found in [`CEUAS/public/adjust`](./CEUAS/public/adjust). Within this task the bias estimates are calculated based on background depatures from [ERA5](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview) ([Hersbach et al. (2018)](#references)). The calculated bias estimates, based on [Haimberger et al. (2008)](https://doi.org/10.1175/2008JCLI1929.1), are added to the CDM compliant dataset. An example can be found here: [Example Homogenization.ipynb](./CEUAS/public/adjust/Example_Homogenization.ipynb)

### Uncertainties

More information can be found in [`CEUAS/public/uncertainties`](./CEUAS/public/uncertainties). Within this task the observation errors are estimated based on a method from [Desrozier et al. (2005), Desrozier (2011)](#references), that uses first guess and analysis departures to estimate a combined observation error.

## References:
*Desroziers, G., Berre, L., Chapnik, B. and Poli, P. (2005), Diagnosis of observation, background and analysis-error statistics in observation space. Q.J.R. Meteorol. Soc., 131: 3385-3396. (Accessed on 23.06.2021), 10.1256/qj.05.108*

*Desroziers, G. (2011): Observation error specification, ECMWF Annual Seminar 2011 [online] Available from: https://www.ecmwf.int/node/14958*

*Haimberger, L., Tavolato, C., and Sperka, S. (2008). Toward Elimination of the Warm Bias in Historic Radiosonde Temperature Records—Some New Results from a Comprehensive Intercomparison of Upper-Air Data. Journal of Climate 21, 18, 4587-4606. (Accessed on 23.06.2021), 10.1175/2008JCLI1929.1*

*Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I., Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2018): ERA5 hourly data on pressure levels from 1979 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). (Accessed on 23.06.2021), 10.24381/cds.bd0915c6*

