
"""
Utilities

Retrieve ERA-Interim fields via the WEBAPI
Convert level units to Pa (or whatever you want)

"""

# retrieve_data_ECMWF.py -s 01-2018 -e 12-2018 --lev --sfc -b
from .retrieve_data_ECMWF import default_request
from .era_level_to_pa import check_level_units