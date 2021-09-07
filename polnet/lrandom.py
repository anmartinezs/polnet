"""
Classes with the implementation of the random models.
A networks is a combination of a polymer in a volume
"""

__author__ = 'Antonio Martinez-Sanchez'

import math
import random
import numpy as np
from abc import ABC, abstractmethod

MAX_TRIES_ELLIP = 1e6


class PGen(ABC):
    """
    Abstract class for random polymer models
    """

    @abstractmethod
    def gen_next_length(self):
        raise NotImplemented


class PGenUniformInRange(PGen):
    """
    Class to model polymer with uniform random length with a range
    """

    def __init__(self, min_l, max_l):
        """
        Constructor
        :param min_l: minimum length
        :param max_l: maximum length
        """
        assert (min_l > 0) and (max_l > min_l)
        self.__min_l, self.__max_l = min_l, max_l
        self.__step = self.__max_l - self.__min_l

    def gen_next_length(self):
        """
        Generate a new length
        :return: a random float value
        """
        return self.__step * random.random() + self.__min_l


class PGenHelixFiber(PGen):
    """
    Class to model random distribution for helix fiber parameters
    """

    def gen_persistence_length(self, min_p):
        """
        Generate a persistence length according an exponential distribution with lambda=1 and minimum value
        :param min_p: minimum persistence value
        :return:
        """
        return min_p + np.random.exponential()

    def gen_zf_length(self, min_zf=0, max_zf=1):
        """
        Generates a z-axis factor within a range
        :param min_zf: minimum value (default 0)
        :param max_zf: maximum value (default 1)
        """
        assert (min_zf >= 0) and (max_zf <= 1)
        return (max_zf - min_zf) * random.random() + min_zf


class PGenHelixFiberB(PGenHelixFiber):
    """
    Class to model random distribution for helix fiber parameters with branches
    """

    def gen_branch(self, b_prob=0.5):
        """
        Generates a boolean that is True (branching) with some input probability
        :param b_prob: branching probability [0, 1) (default 0.5)
        :return: a boolean
        """
        if random.random() <= b_prob:
            return True
        else:
            return False


class SurfGen(ABC):
    """
    Abstract class for modeling 3D parametric surface random parameter generators
    """

    @abstractmethod
    def gen_parameters(self):
        raise NotImplemented


class EllipGen(SurfGen):
    """
    Class for model the paramaters for modelling an Ellipsoid
    """

    def __init__(self, radius_rg, max_ecc=1, max_tries=int(1e6)):
        """
        Constructor
        :param radius_rg: ranges for semi-axis parameters
        :param max_ecc: maximum eccentricity for both planes (default 1), in range [0, 1]
        :param max_tries: raises an RuntimeError exception if there is no proper settings are found before trying
                 'max_tries' times
        """
        assert hasattr(radius_rg, '__len__') and (len(radius_rg) == 2) and (radius_rg[0] <= radius_rg[1])
        assert (max_ecc >= 0) and (max_ecc <= 1)
        assert isinstance(max_tries, int) and (max_tries > 0)
        self.__radius_rg, self.__max_ecc, self.__max_tries = radius_rg, max_ecc, max_tries

    def gen_parameters(self):
        """
        Generates randomly the three semi-axes parameters
        :return: an array with the three semi-axes sorted in descending order
        """
        for i in range(self.__max_tries):
            axes = np.sort(np.asarray((random.uniform(self.__radius_rg[0], self.__radius_rg[1]),
                                       random.uniform(self.__radius_rg[0], self.__radius_rg[1]),
                                       random.uniform(self.__radius_rg[0], self.__radius_rg[1]))))[::-1]
            ecc1, ecc2 = math.sqrt(1 - (axes[2] / axes[0])**2), math.sqrt(1 - (axes[2] / axes[1])**2)
            if (ecc1 <= self.__max_ecc) and (ecc2 <= self.__max_ecc):
                return axes
        raise RuntimeError


class SphGen(SurfGen):
    """
    Class for model the parameters for modelling an Sphere
    """

    def __init__(self, radius_rg):
        """
        Constructor
        :param radius_rg: ranges for radius
        """
        assert hasattr(radius_rg, '__len__') and (len(radius_rg) == 2) and (radius_rg[0] <= radius_rg[1])
        self.__radius_rg = radius_rg

    def gen_parameters(self):
        """
        Generates randomly sphere radii
        :return: a float with the radius value
        """
        return random.uniform(self.__radius_rg[0], self.__radius_rg[1])


class TorGen(SurfGen):
    """
    Class for model the paramaters for modelling a Torus
    """

    def __init__(self, radius_rg):
        """
        Constructor
        :param radius_rg: ranges for radii parameters
        """
        assert hasattr(radius_rg, '__len__') and (len(radius_rg) == 2) and (radius_rg[0] <= radius_rg[1])
        self.__radius_rg = radius_rg

    def gen_parameters(self):
        """
        Generates randomly the two Torus radii
        :return: a array with the two radii ('a' torus radius, and 'b' tube radius, where a > b)
        """
        return np.sort(np.asarray((random.uniform(self.__radius_rg[0], self.__radius_rg[1]),
                                   random.uniform(self.__radius_rg[0], self.__radius_rg[1]))))[::-1]