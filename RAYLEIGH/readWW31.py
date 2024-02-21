#!/usr/bin/env python3
# -*-coding:Utf-8 -*
# Traduction script Matlab corel

import math
import sys, os
import time
import datetime
import matplotlib.pyplot as plt
import numpy as np

#----------------------------------------------------------------------------------#
#launch = 1


#if launch ==1:
def read_dpt(filename):
	"""Reads WAVEWATCH IDLA=3, IDFM=1 output files"""
	f = open(filename, 'r')
        #f = open('ww3.07121700.dpt', 'r')
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

#	print "Converting temp..."
#	print "temp[0] : ", temp[0]
#	for i in range(len(temp)):
#		for j in range(len(temp[0])):
#			temp[i][j] = float(temp[i][j])
#			mat1[i][j] = temp[i][j] * scale
#
#	mat1 = np.flipud(mat1)
#
#	for i in range(len(mat1)):
#		for j in range(len(mat1[0])):
#			if mat1[i][j] == nodata * scale:
#				mat1[i][j] = float('NaN')


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

#	print "... temp converted"
#	print "type(temp[0]) : ", type(temp[0])

#	print "mat1[0] : ", mat1[0]
#	print "len(mat1) : ", len(mat1)
#	print "len(mat1[0]) : ", len(mat1[0])
#	print "nx = ", nx
#	print "ny = ", ny

	if ny > 1:
		lat = np.transpose(np.linspace(latr[0],latr[1],ny))
		lon = np.transpose(np.linspace(lonr[0],lonr[1],nx))
	else:
		lat = np.transpose(np.linspace(latr[0],latr[1],np.floor(np.sqrt(nx))))
		lon = np.transpose(np.linspace(lonr[0],lonr[1],np.floor(np.sqrt(nx))))

#	print "lat, lon, mat1, var : "
#	print "type(lat) : ", type(lat)
#	print "type(lon) : ", type(lon)
#	print "type(mat1) : ", type(mat1)
#	print "type(var) : ", type(var)

	return lat, lon, mat1, var

