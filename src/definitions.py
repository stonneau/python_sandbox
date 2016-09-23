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
	def __init__(self, message):
		self.message = message
		
	def __repr__(self):
		return message

	def __str__(self):
		return message

__EPS = 1e-6
