# -*- coding: utf-8 -*-


def now(timespec='auto'):
    """ Datetime string

    Returns:
        str : datetime now
    """
    import datetime
    return datetime.datetime.now().isoformat(timespec=timespec)


def find_files(directory, pattern, recursive=True):
    """ find files

    Args:
        directory (str): directory path
        pattern (str):  regex string: '*.nc'
        recursive (bool): recursive search?

    Returns:
        list: of files
    """
    import os
    import fnmatch
    matches = []
    if recursive:
        for root, dirnames, filenames in os.walk(directory):
            for filename in fnmatch.filter(filenames, pattern):
                matches.append(os.path.join(root, filename))
    else:
        matches.extend(fnmatch.filter([os.path.join(directory, ifile) for ifile in os.listdir(directory)], pattern))
    return matches


def message(*args, mname=None, verbose=0, level=0, logfile=None, **kwargs):
    if logfile is not None:
        # with open(kwargs['filename'], 'a' if not kwargs.get('force', False) else 'w') as f:
        with open(logfile, 'a') as f:
            f.write(_print_string(*args, **kwargs) + "\n")

    elif verbose > level:
        text = _print_string(*args, **kwargs)
        if mname is not None:
            text = "[%s] " % mname + text

        print(text)
    else:
        pass


def bool_verbose(verbose=0, levels=0, **kwargs):
    return verbose > levels


def _print_string(*args, adddate=False, **kwargs):
    if adddate:
        return "[" + now() + "] " + " ".join([str(i) for i in args])
    else:
        return " ".join([str(i) for i in args])


def dict2str(tmp, sep=', '):
    return sep.join("{!s}={!r}".format(k, v) for (k, v) in tmp.items())


def print_fixed(liste, sep, width, offset=0):
    offset = " " * offset
    out = offset + liste[0]
    n = len(out)
    for i in liste[1:]:
        if (n + len(i) + 1) > width:
            out += sep + "\n" + offset + i
            n = len(offset + i)
        else:
            out += sep + " " + i
            n += len(i) + 2

    return out


def check_kw(name, value=None, **kwargs):
    if value is None:
        return name in kwargs.keys()
    return kwargs.get(name, None) == value


def update_kw(name, value, **kwargs):
    kwargs.update({name: value})
    return kwargs


def levelup(**kwargs):
    kwargs.update({'level': kwargs.get('level', 0) + 1})
    return kwargs


def leveldown(**kwargs):
    kwargs.update({'level': kwargs.get('level', 0) - 1})
    return kwargs


def suche123(eins, zwei, drei, test=None):
    if eins is not test:
        return eins
    elif zwei is not test:
        return zwei
    else:
        return drei


def dict_add(d1, d2):
    d2 = d2.copy()
    for i, j in d1.items():
        if i in d2.keys():
            if ',' in j:
                j = j.split(',')
            else:
                j = [j]

            if ',' in d2[i]:
                k = d2[i].split(',')
            else:
                k = [d2[i]]

            for l in k:
                if l not in j:
                    j.append(l)

            d1[i] = ",".join(j)
            d2.pop(i)
    d1.update(d2)
    return d1


def list_in_list(jlist, ilist):
    """ compare_lists lists and use wildcards

    Args:
        jlist (list): list of search patterns
        ilist (list): list of available choices

    Returns:
        list : common elements
    """
    out = []
    for ivar in jlist:
        if '*' in ivar:
            new = [jvar for jvar in ilist if ivar.replace('*', '') in jvar]
            out.extend(new)
        elif ivar in ilist:
            out.append(ivar)
        else:
            pass
    return out


def get_from_list(ilist, pattern):
    for ivar in ilist:
        if pattern in ivar:
            return ivar
    return None


def dict_in_dict(a, b, method='difference'):
    if method == 'difference':
        return {k: v for k, v in set(a.items()) - set(b.items())}
    elif method == 'union':
        return {k: v for k, v in set(a.items()).union(set(b.items()))}
    elif method == 'intersection':
        return {k: v for k, v in set(a.items()).intersection(set(b.items()))}
    else:
        return {}


def get_data(path=None):
    """
    get a filename form the Module AX directory of the module

    Parameters
    ----------
    path : str
        filename

    Returns
    -------
    str
         path to the file
    """
    import os
    location = os.path.dirname(__file__).replace('/fun', '/ax')
    if path is None:
        print("Choose one: ")
        print("\n".join(os.listdir(os.path.abspath(location))))
    else:
        return os.path.join(os.path.abspath(location), path)


def load_rttov():
    """ Load RTTOV profile limits for quality checks
    """
    import os
    import pandas as pd
    from .. import config
    filename = get_data('rttov_54L_limits.csv')
    if os.path.isfile(filename):
        data = pd.read_csv(filename)
        names = {}
        units = {}
        for i in list(data.columns):
            j = i.split('(')
            iname = j[0].strip().replace(' ', '_')
            iunit = j[1].replace(')', '').strip()
            names[i] = iname
            units[iname] = iunit
        data = data.rename(names, axis=1)
        data = data.set_index('Pressure')
        data = data.to_xarray()
        for i in list(data.data_vars):
            data[i].attrs['units'] = units[i]
        data['Pressure'].attrs['units'] = units['Pressure']

        setattr(config, 'rttov_profile_limits', data)
