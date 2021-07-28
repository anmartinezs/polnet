"""
Classes with the implementation of the random models.
A networks is a combination of a polymer in a volume
"""

__author__ = 'Antonio Martinez-Sanchez'

import random
import numpy as np
from abc import ABC, abstractmethod


class PGen:
    """
    Abstract class for random models
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

