"""Random parameter generators for helical fiber networks.

:class:`FlmsParamGen` draws fiber lengths and persistence lengths
from analytical distributions, providing stochastic parameter
sampling for
:class:`~polnet.sample.filaments.flms_network.NetHelixFiber`.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import random
import numpy as np


class FlmsParamGen:
    """Stochastic parameter generator for helical fiber networks.

    Draws fiber lengths, persistence lengths, z-axis factors, and
    density correction factors from analytical distributions
    (uniform, exponential).
    """

    def gen_length(self, min_l, max_l):
        """Draw a fiber length from a uniform distribution.

        Args:
            min_l (float): Minimum length (must be >= 0 and
                <= max_l).
            max_l (float): Maximum length.

        Returns:
            float: A sampled length in [min_l, max_l].

        Raises:
            ValueError: If min_l < 0 or min_l > max_l.
        """
        if min_l < 0 or min_l > max_l:
            raise ValueError("min_l must be >= 0 and <= max_l.")
        return (max_l - min_l) * random.random() + min_l

    def gen_persistence_length(self, min_p):
        """Draw a persistence length from a shifted Exp(1) distribution.

        Samples from Exp(lambda=1) and adds the minimum value so
        the result is always >= min_p.

        Args:
            min_p (float): Minimum persistence length.

        Returns:
            float: Sampled persistence length >= min_p.
        """
        return min_p + np.random.exponential()

    def gen_zf_length(self, min_zf=0, max_zf=1):
        """Draw a z-axis growth factor from a uniform distribution.

        Controls the bias towards in-plane (xy) versus out-of-plane
        (z) fiber growth.

        Args:
            min_zf (float): Minimum z-factor (default 0).
            max_zf (float): Maximum z-factor (default 1).

        Returns:
            float: A sampled z-factor in [min_zf, max_zf].

        Raises:
            ValueError: If min_zf < 0 or max_zf > 1.
        """
        if min_zf < 0 or max_zf > 1:
            raise ValueError("min_zf must be >= 0 and " "max_zf must be <= 1.")
        return (max_zf - min_zf) * random.random() + min_zf

    def gen_den_cf(self, low, high):
        """Draw a density correction factor from a uniform distribution.

        Args:
            low (float): Minimum value (must be >= 0).
            high (float): Maximum value (must be >= low).

        Returns:
            float: A sampled correction factor in [low, high].

        Raises:
            ValueError: If low < 0 or high < low.
        """
        if low < 0 or high < low:
            raise ValueError("low must be >= 0 and high >= low.")
        return (high - low) * random.random() + low


class HxParamGenBranched(FlmsParamGen):
    """Extends :class:`FlmsParamGen` with stochastic branching.

    Adds a Bernoulli trial to decide at each growth step whether a
    new branch should be initiated.
    """

    def gen_branch(self, b_prob=0.5):
        """Decide stochastically whether to initiate a branch.

        Draws from a Bernoulli(b_prob) distribution.

        Args:
            b_prob (float): Branching probability in [0, 1)
                (default 0.5).

        Returns:
            bool: True if a branch should be created.
        """
        return random.random() <= b_prob
