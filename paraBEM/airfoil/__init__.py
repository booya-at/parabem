from __future__ import division

from copy import copy
import requests
import numpy

from paraBEM import Panel2, PanelVector2


class Airfoil(object):
    data_profiles = None

    def __init__(self, coordinates, name="airfoil"):
        self.coordinates = numpy.array(coordinates)
        self.name = name
        self.noseindex = 0
        self.find_nose()
        self._panels = None
        self.vertices = None

    def __repr__(self):
        return(str(self.coordinates))

    @property
    def x_values(self):
        """Get XValues of airfoil. upper side neg, lower positive"""
        return [((i > self.noseindex) * 2 - 1) * val[0]
                for i, val in enumerate(self.coordinates)]

    @x_values.setter
    def x_values(self, xval):
        """Set X-Values of airfoil to defined points."""
        # assume that coordinates are allready ordered (-1, ..., 0, ... 1)
        # self.coordinates ... numpy array

        coords = copy(self.coordinates)
        self.coordinates = []

        def x(i):
            return coords[i, 0] * ((i > self.noseindex) * 2 - 1)

        def t(i, xval):
            # interpolation:
            # x1-----x1+ t*(x - x1)--------x
            return (xval - x(i-1)) / (x(i) - x(i-1))

        def xy(i, t):
            return coords[i - 1] * (1 - t) + coords[i] * t

        i = 0
        for x_new in xval:
            while x(i) <= x_new:
                if i >= len(coords) - 2:
                    break
                else:
                    i += 1
            ti = t(i, x_new)
            self.coordinates.append(xy(i, ti))
        self.coordinates = numpy.array(self.coordinates)
        self.find_nose()
        self.normalize()

    @property
    def numpoints(self):
        return len(self.coordinates)

    @numpoints.setter
    def numpoints(self, numpoints):
        self.x_values = self.cos_2_distribution(numpoints)

    @property
    def panels(self):
        if self._panels:
            return self._panels
        pans = []
        vertices = []
        coords = self.coordinates[:-1]
        vertices = [PanelVector2(*i) for i in coords]
        vertices[0].wake_vertex = True
        for i, coord in enumerate(coords):
            j = (i+1 if (i + 1) < len(coords) else 0)
            pan = Panel2([vertices[i], vertices[j]])
            pans.append(pan)
        self._vertices = vertices
        self._panels = pans
        return self._panels

    def find_nose(self):
        i = 0
        while (self.coordinates[i + 1][0] < self.coordinates[i][0]
               and i < len(self.coordinates)):
            i += 1
        self.noseindex = i

    def normalize(self, noseindex=None):
        """
        Normalize the airfoil.
        This routine does:
            *Put the nose back to (0,0)
            *De-rotate airfoil
            *Reset its length to 1
        """

        def norm_squared(vec_2):
            return vec_2[0] ** 2 + vec_2[1] ** 2

        p1 = self.coordinates[0]
        nose = self.coordinates[noseindex or self.noseindex]
        diff = p1 - nose  # put nose to (0,0)

        # Angle: a.b=|a|*|b|*sin(alpha)
        sin_sq = diff.dot([0, -1]) / norm_squared(diff)
        cos_sq = diff.dot([1, 0]) / norm_squared(diff)
        # de-rotate and scale
        matrix = numpy.array([[cos_sq, -sin_sq], [sin_sq, cos_sq]])
        self.coordinates = numpy.array([matrix.dot(i - nose) for i in self.coordinates])
        self.coordinates[-1] = self.coordinates[0]

    @classmethod
    def compute_naca(cls, naca=1234, numpoints=100):
        """Compute and return a four-digit naca-airfoil"""
        # See: http://people.clarkson.edu/~pmarzocc/AE429/The%20NACA%20airfoil%20series.pdf
        # and: http://airfoiltools.com/airfoil/naca4digit
        m = int(naca / 1000) * 0.01  # Maximum Camber Position
        # second digit: Maximum Thickness position
        p = int((naca % 1000) / 100) * 0.1
        t = (naca % 100) * 0.01  # last two digits: Maximum Thickness(%)
        x_values = [
            1 - numpy.sin((x * 1. / (numpoints - 1)) * numpy.pi / 2) for x in range(numpoints)]

        upper = []
        lower = []
        a0 = 0.2969
        a1 = -0.126
        a2 = -0.3516
        a3 = 0.2843
        a4 = -0.1036            # modified for closed profile

        for x in x_values:
            if x < p:
                mean_camber = (m / (p ** 2) * (2 * p * x - x ** 2))
                gradient = 2 * m / (p ** 2) * (p - x)
            else:
                mean_camber = (
                    m / ((1 - p) ** 2) * ((1 - 2 * p) + 2 * p * x - x ** 2))
                gradient = 2 * m / (1 - p ** 2) * (p - x)

            thickness_this = t / 0.2 * \
                (a0 * numpy.sqrt(x) + a1 * x + a2 * x **
                 2 + a3 * x ** 3 + a4 * x ** 4)
            #theta = math.atan(gradient)
            costheta = (1 + gradient ** 2) ** (-0.5)
            sintheta = gradient * costheta
            upper.append([x - thickness_this * sintheta,
                          mean_camber + thickness_this * costheta])
            lower.append([x + thickness_this * sintheta,
                          mean_camber - thickness_this * costheta])
        return cls(upper + lower[::-1][1:], name="NACA_" + str(naca))

    @classmethod
    def joukowsky(cls, m=-0.1+0.1j, numpoints=100):
        from paraBEM.airfoil.conformal_mapping import JoukowskyAirfoil
        airfoil = JoukowskyAirfoil(m)
        profile = [[c.real, c.imag] for c in airfoil.coordinates(numpoints)]

        # find the smallest xvalue to reset the nose
        x = numpy.array([i[0] for i in profile])
        profile = cls(profile, "joukowsky_" + str(m))
        profile.normalize(numpy.where(x == min(x))[0][0])
        profile.normalize()
        profile.numpoints = numpoints
        return profile

    @classmethod
    def vandevooren(cls, tau=0.05, epsilon=0.05, numpoints=100):
        from paraBEM.airfoil.conformal_mapping import VanDeVoorenAirfoil
        airfoil = VanDeVoorenAirfoil(tau=tau, epsilon=epsilon)
        profile = [[c.real, c.imag] for c in airfoil.coordinates(numpoints)]

        # find the smallest xvalue to reset the nose
        x = numpy.array([i[0] for i in profile])
        profile = cls(profile, "VanDeVooren_tau=" + str(tau) + "_epsilon=" + str(epsilon))
        profile.normalize(numpy.where(x == min(x))[0][0])
        profile.normalize()
        profile.numpoints = numpoints
        return profile

    @classmethod
    def trefftz_kutta(cls, m=-0.1+0.1j, tau=0.05, numpoints=100):
        from paraBEM.airfoil.conformal_mapping import TrefftzKuttaAirfoil
        airfoil = TrefftzKuttaAirfoil(midpoint=m, tau=tau)
        profile = [[c.real, c.imag] for c in airfoil.coordinates(numpoints)]

        # find the smallest xvalue to reset the nose
        x = numpy.array([i[0] for i in profile])
        profile = cls(profile, "TrefftzKuttaAirfoil_m=" + str(m) + "_tau=" + str(tau))
        profile.normalize(numpy.where(x == min(x))[0][0])
        profile.normalize()
        profile.numpoints = numpoints
        return profile

    @classmethod
    def import_from_database(cls, name="joukowsky"):
        resp = requests.get(
            "http://m-selig.ae.illinois.edu/ads/coord/" + name + ".dat")
        return Airfoil.import_from_text(resp.iter_lines(), name=name)

    @classmethod
    def list_from_database(cls, name=""):
        if not cls.data_profiles:
            text = requests.get("http://m-selig.ae.illinois.edu/ads/coord/")
            profiles = []
            for i in text.iter_lines():
                if b".dat" in i:
                    j = i.split(b'"')
                    for k in j:
                        if b".dat" in k and k.decode()[0] is not '>':
                            profiles.append(k.decode()[:-4])
            cls.data_profiles = profiles
        sim_profiles = []
        for i in cls.data_profiles:
            if name in i:
                sim_profiles.append(i)
        return sim_profiles

    @classmethod
    def import_from_text(cls, text, name=None):
        name = name
        profile = []
        for i, line in enumerate(text):
            if i > 0:
                split_line = line.split()
                if len(split_line) == 2:
                    profile.append([float(i) for i in split_line])
                else:
                    name = line
        return cls(profile, name)

    @classmethod
    def import_from_dat(cls, path):
        """
        Import an airfoil from a '.dat' file
        """
        name = 'imported from {}'.format(path)
        with open(path, "r") as p_file:
            return Airfoil.import_from_text(p_file, name)

    def export_dat(self, pfad):
        """
        Export airfoil to .dat Format
        """
        with open(pfad, "w") as out:
            if self.name:
                out.write(str(self.name))
            for i in self.coordinates:
                out.write("\n" + str(i[0]) + "\t" + str(i[1]))
        return pfad

    @staticmethod
    def cos_distribution(numpoints):
        """
        return cosinus distributed x-values
        """
        numpoints -= numpoints % 2
        xtemp = lambda x: ((x > 0.5) - (x < 0.5)) * \
            (1 - numpy.sin(numpy.pi * x))
        return [xtemp(i / numpoints) for i in range(numpoints + 1)]

    @staticmethod
    def cos_2_distribution(numpoints):
        """
        return cosinus distributed x-values
        double-cosinus -> neat distribution at nose and trailing edge
        """
        numpoints -= numpoints % 2
        xtemp = lambda x: ((x > 0.5) - (x < 0.5)) * \
            (1 + numpy.cos(2 * numpy.pi * x)) / 2
        return [xtemp(i / numpoints) for i in range(numpoints + 1)]
