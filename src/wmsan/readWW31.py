#!/usr/bin/env python3
# -*-coding:Utf-8 -*

# Preamble
__author__ = "Lisa Tomasetto"
__copyright__ = "Copyright 2024, UGA"
__credits__ = ["Lisa Tomasetto"]
__version__ = "0.1"
__maintainer__ = "Lisa Tomasetto"
__email__ = "lisa.tomasetto@univ-grenoble-alpes.fr"
__status__ = "Development"

"""
This program aims at reading WW3 default bathymetry files (adapted from Lucia GUALTIERI and Pierre BOUÉ)
"""

##################################################################################
import math
import sys, os
import time
import datetime
import matplotlib.pyplot as plt
import numpy as np

def read_dpt(filename = '../../data/ww3.07121700.dpt'):
	"""Read default bathymetry dpt file"""
	f = open(filename, 'r')
	line1 = f.readline()
	DATA = f.read()
	f.close()

	DATA = DATA.split()
	DATA = [int(x) for x in DATA]

	line1 = line1.split()
	modelname = line1[0] + ' '+ line1[1]
	modelname2 = line1[0] + ' '+ line1[1]
	date = [line1[2], line1[3]]
	lonr = [ float(line1[4]), float(line1[5])]
	nx = int(line1[6])
	latr = [float(line1[7]), float(line1[8])]
	ny = int(line1[9])
	var = str(line1[10])
	scale = float(line1[11])
	tv1  = len(var)
	unit = str(line1[12])
	idla = int(line1[13])
	idlafm = int(line1[14])
	form = str(line1[15])
	nodata = int(line1[16])
	temp = np.empty((ny,nx))
	ind = 0
	for i in range(ny):
		for j in range(nx):
			temp[i][j] = int(DATA[ind])
			ind = ind+1

	mat1 = np.ones((ny,nx))
	mat1 = mat1 * np.NaN

	for i in range(len(temp)):
		for j in range(len(temp[0])):
			if temp[i][j] == nodata:
				mat1[i][j] = float('NaN')
			else:
				mat1[i][j] = scale * float(temp[i][j])

	mat1 = np.flipud(mat1)

	for i in range(len(mat1)):
		for j in range(len(mat1[0])):
			if mat1[i][j] == nodata * scale:
				mat1[i][j] = float('NaN')

	if ny > 1:
		lat = np.transpose(np.linspace(latr[0],latr[1],ny))
		lon = np.transpose(np.linspace(lonr[0],lonr[1],nx))
	else:
		lat = np.transpose(np.linspace(latr[0],latr[1],np.floor(np.sqrt(nx))))
		lon = np.transpose(np.linspace(lonr[0],lonr[1],np.floor(np.sqrt(nx))))

	return lat, lon, mat1, var