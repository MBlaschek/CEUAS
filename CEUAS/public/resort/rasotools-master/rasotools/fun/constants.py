class constant(object):
    """ Add pint unit aware conversion to constants

    """
    def __init__(self, value, units, desc):
        self.value = value
        self.units = units
        self.description = desc

    def __repr__(self):
        return "%14e [%s] %s" % (self.value, self.units, self.description)

    def __call__(self):
        return self.value

    #  -- MATH RULES --
    def __mul__(self, b):
        if isinstance(b, constant):
            return self.value * b.value
        else:
            return self.value * b

    def __rmul__(self, b):
        return self.value * b

    def __truediv__(self, b):
        if isinstance(b, constant):
            return self.value / b.value
        else:
            return self.value / b

    def __rtruediv__(self, b):
        if isinstance(b, constant):
            return b.value / self.value
        else:
            return b / self.value

    def __pow__(self, b):
        if isinstance(b, constant):
            return self.value ** b.value
        else:
            return self.value ** b

    def __rpow__(self, b):
        if isinstance(b, constant):
            return b.value ** self.value
        else:
            return b ** self.value

    def __sub__(self, b):
        if isinstance(b, constant):
            return self.value - b.value
        else:
            return self.value - b

    def __rsub__(self, b):
        if isinstance(b, constant):
            return b.value - self.value
        else:
            return b - self.value

    def __add__(self, b):
        if isinstance(b, constant):
            return self.value + b.value
        else:
            return self.value + b

    def __radd__(self, b):
        return self.value + b


#  --- ...  Math constants

pi = constant(3.1415926535897931, "rad", 'pi')
sqrt2 = constant(1.414214e+0, '', "square root of 2")
sqrt3 = constant(1.732051e+0, '', "square root of 3")

# Physical constants
Cpd = constant(1005.7, 'J K-1 kg-1', 'Spec. heat dry air at const. press.')
Cvd = constant(718., 'J K-1 kg-1', 'Spec. heat dry air at const. volume')
Cpv = constant(1870., 'J K-1 kg-1', 'Specific heat of water vapour')
Cl = constant(4190., 'J K-1 kg-1', 'Spec heat liquid water')
Rv = constant(461.5, 'J K-1 kg-1', 'Gas constant for water vapour')
Rd = constant(287.04, 'J K-1 kg-1', 'Gas constant dry air')
Lv = constant(2.501e6, ' ', 'Latent heat vaporisation')
Lf = constant(3.336e5, ' ', 'Lat heat of fusion')
rowl = constant(1000., 'kg m-3', 'Density of liquid water')
stebol = constant(5.673e-8, 'W m-2 K-4 ', 'Stefan Boltzmann constant')
epsilon = constant(287.04 / 461.5, '-', 'Rd/Rv')

# Orbital and planetary parameters
eccen = constant(0.016715, '-', 'Eccentricity of orbit')
obliq = constant(23.4441, 'deg', 'Planetary obliquity')
prece = constant(102.7, 'deg', 'Precession of vernal equinox')
orb_year = constant(1995, '-', 'Earth year for which orbital params computed')
lod = constant(86400., 's', 'Length of day')
daysperyear = constant(365.25, '-', 'Number of days in a year')
calday = constant(80.5, 'days', 'Julian calendar day')
r = constant(6374000., 'm', 'Planetary radius')
g = constant(9.81, 'm/s2', 'Gravitational acceleration')
omega = constant(2. * pi / 86400., 'rad/s', 'Planetary rotation rate')
radius = constant(1., 'AU', 'Mean orbital radius')
scon = constant(1367., 'W/m2', 'Solar irradiance at mean orbital radius')

#  --- ...  Geophysics/Astronomy constants

rerth = constant(6.3712e+6, 'm', 'radius of earth')
g = constant(9.80665e+0, 'm/s2', 'gravity           ')
omega = constant(7.2921e-5, '1/s', 'ang vel of earth  ')
p0 = constant(1.01325e5, 'pa', 'std atms pressure ')
# solr  =_constant(1.36822e+3  ,'W/m2'  ,'solar constant    ')# -aer(2001)
# solr  =_constant(1.3660e+3  ,'W/m2'     ,'solar constant    ') #(W/m2)-liu(2002)
solr = constant(1.36742732e+3, 'W/m2', 'solar constant    ')  # (W/m2)-gfdl(1989) - OPR as of Jan 2006

#  --- ...  Thermodynamics constants

rgas = constant(8.314472, 'J/mol/K', 'molar gas constant  ')
rd = constant(2.8705e+2, 'J/kg/K', 'gas constant air   ')
rv = constant(4.6150e+2, 'J/kg/K', 'gas constant H2O   ')
cp = constant(1.0046e+3, 'J/kg/K', 'spec heat air @p   ')
cv = constant(7.1760e+2, 'J/kg/K', 'spec heat air @v    ')
cvap = constant(1.8460e+3, 'J/kg/K', ' spec heat H2O gas  ')
cliq = constant(4.1855e+3, 'J/kg/K', 'spec heat H2O liq   ')
csol = constant(2.1060e+3, 'J/kg/K', 'spec heat H2O ice   ')
hvap = constant(2.5000e+6, 'J/kg', 'lat heat H2O cond   ')
hfus = constant(3.3358e+5, 'J/kg', 'lat heat H2O fusion ')
psat = constant(6.1078e+2, 'Pa', 'pres at H2O 3pt    ')
t0c = constant(2.7315e+2, 'K', 'temp at 0C         ')
ttp = constant(2.7316e+2, 'K', 'temp at H2O 3pt    ')
tice = constant(2.7120e+2, 'K', 'temp freezing sea  ')
jcal = constant(4.1855E+0, '', 'joules per calorie  ')

#  Secondary constants

rocp = constant(rd / cp, '1', 'rd / cp')  # J/kg/k / J/Kg/K
cpor = constant(cp / rd, '1', ' cp / rd ')
rog = constant(rd / g, '', 'rd / g')  # J/kg/K / m/s2
fvirt = constant(rv / rd - 1., '1', 'rv/rd -1')
eps = constant(rd / rv, '1', ' rd/rv')
epsm1 = constant(rd / rv - 1., '1', 'rd/rv -1')
dldt = constant(cvap - cliq, '1', 'cvap - cliq')
xpona = constant(-1 * dldt / rv, '', '-(cvap-cliq)/rv')
xponb = constant(-1 * dldt / rv + hvap / (rv * ttp), '', '-(cvap-cliq)/rv + hvap/(rv*ttp)')

#  --- ...  Other Physics/Chemistry constants

plnk = constant(6.6260693e-34, 'J/s', 'planck constant')  # )  (J/s) -nist(2002)
sbc = constant(5.6730e-8, 'W/m2/K4', 'stefan-boltzmann  ')  # (W/m2/K4)
# sbc   =_constant(5.670400e-8    ,'stefan-boltzmann    (W/m2/K4) -nist(2002)
# avgd  =_constant(6.02214e23     ,'avogadro constant   (1/mol) -aer
avgd = constant(6.0221415e23, '1/mol', 'avogadro constant   ')  # -nist(2002)
gasv = constant(22413.996e-6, 'm3/mol', 'vol of ideal gas at 273.15k, 101.325kpa ')  # (m3/mol) -nist(2002)
# amd   =_constant(28.970    ,'g/mol'     ,'molecular wght of dry air ') #(g/mol)
amd = constant(28.9644, 'g/mol', 'molecular wght of dry air ')  # (g/mol)
amw = constant(18.0154, 'g/mol', 'molecular wght of water vapor ')  # (g/mol)
amo3 = constant(47.9982, 'g/mol', 'molecular wght of o3 ')  # (g/mol)
# amo3  =_constant(48.0      ,'g/mol '    ,'molecular wght of o3  ') #(g/mol)
amco2 = constant(44.011, 'g/mol', 'molecular wght of co2 ')  # (g/mol)
amo2 = constant(31.9999, 'g/mol', 'molecular wght of o2  ')  # (g/mol)
amch4 = constant(16.043, 'g/mol', 'molecular wght of ch4 ')  # (g/mol)
amn2o = constant(44.013, 'g/mol', 'molecular wght of n2o ')  # (g/mol)
