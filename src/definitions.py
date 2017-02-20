# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21

@author: stonneau
"""

class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class OptimError(Error):
	"""Raised when problem is not solved"""
	def __init__(self, Error):
		self.message = Error
		
	def __repr__(self):
		return Error

	def __str__(self):
		return Error

__EPS = 1e-6
