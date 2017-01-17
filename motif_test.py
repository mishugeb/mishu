# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 22:29:36 2017

@author: Ashiqul
"""
from __future__ import division
#file is to test the import of the feature_module
import os 
import sys
print os.getcwd()
os.chdir("C:/Users/Ashiqul/Desktop")
import feature_module as m
reload(feature_module)
protein = "ASDG"
print m.extract_motifs_pos(protein)

