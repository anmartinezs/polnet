"""
Classes with the implementation of the random models.
A networks is a combination of a polymer in a volume
"""

__author__ = 'Antonio Martinez-Sanchez'

import random
import numpy as np
from abc import ABC, abstractmethod


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


class EllipGen:
    """
    Class for model the paramaters for modelling an Ellipsoid
    """

    def gen_parameters(self, radius_rg):
        """
        Generates randomly the three semi-axis parameters
        :param radius_rg: ranges for semi-axis parameters
        :return: a 3-tuple of floats
        """
        assert hasattr(radius_rg, '__len__') and (len(radius_rg) == 2) and (radius_rg[0] <= radius_rg[1])
        return random.uniform(radius_rg[0], radius_rg[1]), \
               random.uniform(radius_rg[0], radius_rg[1]), \
               random.uniform(radius_rg[0], radius_rg[1])