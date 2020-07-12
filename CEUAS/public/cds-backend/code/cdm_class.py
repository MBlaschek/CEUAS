"""

Create a class object based on Xarray that can read the CDM and use the information
or h5Py


tmp = cdm_dataset('..merged.nc')
tmp
 - groups
    - variables + metadata + shape
    - read information for datetime and z_coordinate
    - observed variables

read_group(variables, datetime_var) - convert to datetimes

write_group(name, data, attributes)

read_to_cube(group, variables)
 - group : observation_table
 - groupby:
    - observed_variable -> 85, ...
    - date_time
    - z_coordinate


"""
import numpy as np
import h5py
import xarray as xr

odb_codes = {'t': 85, 'rh': 38, 'td': 36, 'dpd': 34, 'z': 117, 'dd': 106, 'ff': 107, 'u': 104, 'v': 105,
             'q': 39}
cdm_codes = {'z': 1, 't': 2, 'u': 3, 'v': 4, 'dd': 111, 'ff': 112, 'td': 59, 'dpd': 299, 'rh': 29, 'p': 999}

std_plevs = np.asarray([10, 20, 30, 50, 70, 100, 150, 200, 250, 300, 400, 500, 700, 850, 925, 1000])


@xr.register_dataset_accessor('cdm')
class CDMAccessor(object):
    def __init__(self, index, **kwargs):
        self.index = index
        self._names = []
        for ikey, ival in kwargs.items():
            setattr(self, ikey, ival)
            self._names.append(ikey)

    def keys(self):
        return self._names


class CDMVariable:
    def __init__(self, origin, name: str, **kwargs):
        self._name = name
        self._origin = origin
        self._names = []
        for ikey, ival in kwargs.items():
            setattr(self, ikey, ival)
            self._names.append(ikey)

    def __setitem__(self, key, value):
        self._names.append(key)
        self.__setattr__(key, value)

    def __getitem__(self, item):
        return self._origin[item]

    def __repr__(self):
        return " ".join(["{}".format(str(getattr(self, i))) for i in self._names])

    def keys(self):
        return self._names

    def attrs(self, name=None):
        if name is not None:
            if name in self._origin.attrs.keys():
                return self._origin.attrs[name]
            else:
                return None
        else:
            return self._origin.attrs.keys()


class CDMGroup(CDMVariable):
    def __repr__(self):
        return self._name + ":\n" + "\n".join(["{:_<50} : {}".format(i, getattr(self, i)) for i in self.keys()])


class CDMDataset:
    # __slots__ = ['filename', 'file', 'groups', 'data']

    def __init__(self, filename: str):
        self.filename = filename
        self.file = h5py.File(filename, 'r')
        print("[OPEN] %s" % self.filename)
        self.groups = []
        # Check everything
        self.inquire()

    def __getitem__(self, item):
        return self.__getattribute__(item)

    def __repr__(self):
        text = "Filename: " + self.filename
        text += "\nGroups: \n"
        text += "\n".join(["- {} - {} : {}".format('G' if isinstance(self[igroup], CDMGroup) else 'V', igroup,
                                                   getattr(self[igroup], 'shape')) for igroup in self.groups])
        return text

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

    def inquire(self):
        try:
            for igroup in self.file.keys():
                self.groups.append(igroup)
                if isinstance(self.file[igroup], h5py.Group):
                    jgroup = CDMGroup(self.file[igroup], igroup)

                    for ivar in self.file[igroup].keys():
                        shape = self.file[igroup][ivar].shape
                        jgroup[ivar] = CDMVariable(self.file[igroup][ivar], ivar, shape=shape)

                    jgroup['shape'] = len(jgroup.keys())
                    setattr(self, igroup, jgroup)  # attach to class
                else:
                    setattr(self, igroup, CDMVariable(self.file[igroup], igroup, shape=self.file[igroup].shape))

        except Exception as e:
            print(repr(e))
            self.close()

    def close(self):
        self.file.close()
        print("[CLOSED] %s" % self.filename)

    def reopen(self, mode: str = 'r', **kwargs):
        self.file.close()
        print("[REOPEN] %s [%s]" % (self.filename, mode))
        self.file = h5py.File(self.filename, mode=mode, **kwargs)

    def read_group(self, group: str, variables: str or list, decode_byte_array: bool = True, **kwargs):
        """ Return data from variables and concat string-byte arrays


        Args:
            group:
            variables:
            decode_byte_array:
            **kwargs:

        Returns:

        """
        if not isinstance(variables, list):
            variables = [variables]

        data = {}
        if group in self.groups:
            jgroup = getattr(self, group)
            for ivar in variables:
                if ivar in jgroup.keys():
                    if len(jgroup[ivar].shape) > 1 and decode_byte_array:
                        data[ivar] = jgroup[ivar][()].astype(object).sum(axis=1).astype(str)  # concat byte-char-array
                    else:
                        data[ivar] = jgroup[ivar][()]  # get the full numpy array
        return data

    def write_group(self, group: str, data: dict, **kwargs):
        if group in self.groups:
            # WOWOWOWOW
            pass
        else:
            if self.file.mode == 'r':
                self.reopen(mode='a')
                # igroup = self.file.create_group()

    def read_data_to_cube(self, group: str, variables: str or list,
                          datetime: str = None,
                          z_coordinate: str = None,
                          observed_variable: str = None,
                          **kwargs):
        """

        Args:
            group:
            variables:
            datetime:
            z_coordinate:
            observed_variable:
            **kwargs:

        Returns:


        Notes:
            Converting from ragged arraay to Cube requires a index
            save the index for the group
            for later writing
            return a Xarray Dataset with cdm Accessor
        """
        if group not in self.groups:
            raise ValueError("Group not found", group)

        if isinstance(variables, str):
            variables = [variables]

        for ivar in variables:
            if ivar not in self[group].keys():
                raise ValueError("Variable not in Group", ivar, group)

        if datetime is not None:
            if datetime in self[group].keys():
                raise ValueError("Variable not in Group", datetime, group)

        if z_coordinate is not None:
            if z_coordinate in self[group].keys():
                raise ValueError("Variable not in Group", z_coordinate, group)

        if observed_variable is not None:
            read_ragged = True
            if observed_variable in self[group].keys():
                raise ValueError("Variable not in Group", observed_variable, group)
        else:
            read_ragged = False

        if datetime is None and z_coordinate is None:
            raise RuntimeError("Requires at least one varibale to make an index from [datetime, z_coordinate]")

        if datetime is not None:
            time = self[group][datetime][()]  # everything
            units = self[group][datetime].attrs('units')
            if units is None:
                units = 'seconds since 1900-01-01 00:00:00'
            else:
                units = units.decode()

        if z_coordinate is not None:
            plev = self[group][z_coordinate][()]  #  everything
            punits = self[group][z_coordinate].attrs('units')
            if punits is None:
                pass

            iplev = np.in1d(plev, std_plevs * 100)  # only std pressure (can have NaN)
            plev = plev.astype(np.int32) // 100  # needs to be integer for indexing, hPa as well

        # recordindex, recordtimestamp -> start end of profiles
        # profiles = np.split(times, recidx)  -> list of profiles []
        if read_ragged:
            # variable x time x plev
            # for every variable run makeindex_variable(observed_variable, ivar, date_time, -1, -1, z_ccordinate, 0, mask)
            pass

        else:
            # just time or just pressure or together [time x plev]
            pass


@njit(cache=True)
def makeindex_variable(obsvar, var, datevar, start, stop, pvar, plist, mask):
    for i in range(mask.shape[0]):
        # Variable numbers agree
        if obsvar[i] == var[0]:
            # There is date info
            found = False
            if start[0] == -1:
                if start[0] <= datevar[i] <= stop[0]:
                    found = True

            if plist[0] == 0:
                for j in range(plist.shape[0]):
                    if pvar[i] == plist[j]:
                        found = True
                        break
            mask[i] = found