#!/home/vatir/.conda/envs/vpy-2.7.3/bin/python
# -*- coding: utf-8 -*-

"""
Created on Wed Jul 16 11:09:07 2014

@author: Koan
"""
from subprocess import check_output, STDOUT
if __name__ == '__main__':
	#Code=["sleep 10","date","sleep 20"]
	
	import numpy as np

	N = 1
	
	WorkingName = "TriStackedRingODEStatic"
	#WorkingName = "TriStackedRingODESlowSynDeg"
	#WorkingName = "TriStackedRingODE1GenPerHour"
	#WorkingName = "TriStackedRingODE2GenPerHour"
	
	import h5py
	File = h5py.File(".\\Data\\{}{}.h5".format(N, WorkingName), "r")
	import importlib
	if N == 1:
		CD = importlib.import_module("Data.{}{}".format("",WorkingName))
		CD.CDImportName = "Data.{}{}".format("",WorkingName)
	else:
		CD = importlib.import_module("Data.{}{}".format(N-1,WorkingName))
		CD.CDImportName = "Data.{}{}".format(N-1,WorkingName)
	
	FullOutArray = File["Results"][:]
	Errors = File["Errors"][:]
	HasData = File["HasData"][:]
	Runtimes = File["Runtimes"][:]
	A0Range = File["A0Range"][:]
	KD1Range = File["KD1Range"][:]
	KD2Range = File["KD2Range"][:]
	TimePoints = CD.KeepTimePoints
	SpeciesCount = CD.SpeciesCount
	
	
	#KD1Indicies, KD2Indicies = np.where(Errors==True) # Get from File

	
	for Index in range(len(KD1Indicies)):
	#for Index in range(18,30):
		i = KD1Indicies[Index]
		j = KD2Indicies[Index]
		
		Code = ("python Worker.py {} {} {} {} {} 1".format(CD.CDImportName, i, i, j, j))
		JobName = "{}{}:{}:{}:{}:{}".format(CD.UniquePrefix, CD.BatchType, i, i, j, j)
		print "Index: {}".format(Index)
		print JobName
		print check_output(Code, stderr = STDOUT)
