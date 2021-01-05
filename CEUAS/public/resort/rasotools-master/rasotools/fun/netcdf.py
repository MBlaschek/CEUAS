# -*- coding: utf-8 -*-

__all__ = ['view']


def view(filename, show_attributes=False):
    """ NCVIEW

    Args:
        filename:
        show_attributes:

    """
    import os
    import netCDF4 as nc

    if not os.path.isfile(filename):
        raise IOError("Unknonw file: %s" % filename)

    if 'gz' in filename:
        import gzip
        with gzip.open(filename) as g:
            with nc.Dataset("dummy", 'r', memory=g.read()) as f:
                # HEADER INFORMATION
                print("File: %s" % filename)
                _header(f)
                # other variables ?
                if show_attributes:
                    print()
                    _detailed_information(f)

                print()
    else:
        with nc.Dataset(filename, 'r', ) as f:
            # HEADER INFORMATION
            print("File: %s" % filename)
            _header(f)
            # other variables ?
            if show_attributes:
                print()
                _detailed_information(f)

            print()


def _header(f):
    print("Dimensions: ")
    ns = max([len(i) for i in f.dimensions.keys()])
    for idim in f.dimensions:
        print("   %*s : %-8d %s" % (ns,
                                    idim, f.dimensions[idim].size,
                                    '[unlimited]' if f.dimensions[idim].isunlimited() else ''))

    print("Variables: ")
    ns = max([len(i) for i in f.variables.keys()])
    for ivar in f.variables:
        if ivar in f.dimensions:
            continue
        print("   %*s : %-20s %-15s [%s]" % (
            ns, ivar, ",".join(f.variables[ivar].dimensions), str(f.variables[ivar].shape),
            str(getattr(f.variables[ivar], 'units', '1'))))


def _detailed_information(f):
    ns = max([len(i) for i in f.variables.keys()])
    for ivar in f.variables:
        if ivar in f.dimensions.keys():
            continue
        print(">>   %*s : %-20s %-15s [%s]" % (
            ns, ivar, ",".join(f.variables[ivar].dimensions), str(f.variables[ivar].shape),
            str(getattr(f.variables[ivar], 'units', '1'))))
        ms = max([len(i) for i in f.variables[ivar].ncattrs()])
        for iatt in f.variables[ivar].ncattrs():
            print("   %*s : %s" % (ms, iatt, f.variables[ivar].getncattr(iatt)))  # getattr(f.variables[ivar], iatt))

        print()

    print("#" * 80)
    for ivar in f.dimensions:
        print(">>   %*s : %s %s [%s]" % (
            ns, ivar, ",".join(f.variables[ivar].dimensions), str(f.variables[ivar].shape),
            str(getattr(f.variables[ivar], 'units', '1'))))
        ms = max([len(i) for i in f.variables[ivar].ncattrs()])
        for iatt in f.variables[ivar].ncattrs():
            print("     %*s : %s" % (ms, iatt, f.variables[ivar].getncattr(iatt)))  # getattr(f.variables[ivar], iatt))
