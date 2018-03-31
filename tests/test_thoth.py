#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `thoth` package."""

import unittest

import numpy as np
import numpy.testing as npt
# https://docs.scipy.org/doc/numpy-dev/reference/generated/numpy.testing.assert_allclose.html

from thoth import thoth


class TestThoth(unittest.TestCase):
    """Tests for `thoth` package.

    .. warning:: The expected values are copied from the first test run. # TODO

    """

    def assertAllClose(self, *args, **kwargs):
        """Call numpy.testing.assert_allclose() with verbose=True"""
        kwargs['verbose'] = True
        # kwargs['rtol'] = 1e-07  # relative tolerance
        # kwargs['atol'] = 0      # absolute tolerance
        return npt.assert_allclose(*args, **kwargs)

    def test_100_entropy(self):
        a_prob = thoth.prob_from_array([2, 3, 12, 5, 7])
        result = thoth.norm_prob(a_prob)
        expected = []
        self.assertEqual(result, expected)

    def test_200_entropy(self):
        a_prob = thoth.prob_from_array([2, 3, 12, 5, 7])
        result = thoth.entropy(a_prob)
        expected = 2.063651255829042
        self.assertEqual(result, expected)

    def test_300_entropy_nsb(self):
        a_prob = thoth.prob_from_array([2, 3, 12, 5, 7])
        result = thoth.entropy_nsb(a_prob)
        expected = 2.112851131072417
        self.assertEqual(result, expected)

    def test_400_calc_entropy(self):
        result = thoth.calc_entropy([2, 3, 12, 5, 7], 100000)
        expected = [2.17115643, 2.00301591, 2.34095622, 1.89913747, 2.57355048]
        self.assertAllClose(result, expected)

    def test_500_jsd(self):
        result = thoth.calc_jsd([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], 0.5, 100000)
        expected = np.array(
            [0.0524735, -0.09003991,  0.19567123, -0.26433774,  0.28254133])
        self.assertAllClose(result, expected)

    def test_600_mi(self):
        arrays = np.array([[1, 2, 3], [3, 5, 2], [3, 2, 4]])
        result = thoth.calc_mi(arrays, 100000)
        expected = np.array(
            [-0.05951304, -0.18778695,  0.0685165, -0.37435192,  0.12983761])
        self.assertAllClose(result, expected)
