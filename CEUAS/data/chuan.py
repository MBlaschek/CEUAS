# -*- coding: utf-8 -*-


def chuan(ident, filename, variables=None, levels=None, **kwargs):
    """ Read CHUAN station

    Args:
        ident (str): WMO ID
        filename (str): filename to read from
        variables (list): select only these variables
        levels (list): interpolate to these pressure levels [Pa]
        **kwargs:

    Returns:
        Dataset : xarray Dataset
    """
    from .igra.read import to_std_levels
    import xarray as xr

    if '.nc' in filename:
        data = xr.open_dataset(filename, **kwargs)
    else:
        # same data format as uadb -> use function
        data = to_std_levels(ident, filename, levels=levels, uadb=True, **kwargs)

    if variables is not None:
        avail = list(data.data_vars.keys())
        if not isinstance(variables, list):
            variables = [variables]

        variables = [iv for iv in variables if iv in avail]
        if len(variables) > 0:
            data = data[variables]  # subset

    data.attrs.update({'dataset': 'CHUAN, ds352.0'})
    return data
