# Integration into the CDS

[Confluence Page](https://confluence.ecmwf.int/display/COPSRV/How+to+integrate+a+dataset+into+the+CDS+Catalogue)



## Required files:

1. `abstract.md`, Summary of the data set
   1. NAME of the dataset ?
   2. Add Reference to gridded version and intercomparisons
2. `dataset.yaml` Data set name, license, ...
3. `overview.png` Main figure
4. `variables.yaml` Definition of the variables and description
5. `abstract.yaml` Definition of the data set (resolution)
6. `related_variables.yaml`Definition of additional variables available in the output
7. `documentation.yaml` Defining the Documentation Tab, e.g. PUG
8. `Product_User_Guide.pdf`, PUG



Check YAML files with [Yaml-Lint](http://www.yamllint.com/)

Convert between format with [convert table](https://tableconvert.com/)

