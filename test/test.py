#!/usr/bin/env python2

'''
(c) Fritjof Lammers

fritjoflammers@gmail.com

'''


from parse_TAREAN import csv_parser
import unittest
import sys


class TestParser(unittest.TestCase):

    def testparser(self):
        """
        Test if CSV parser runs through.Plain and simple.
        :return:
        """
        csv_parser("Test.csv", "Test.out.csv")