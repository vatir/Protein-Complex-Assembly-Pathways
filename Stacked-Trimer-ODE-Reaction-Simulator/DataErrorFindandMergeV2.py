#!/home/vatir/.conda/envs/vpy-2.7.3/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import h5py
import DataManager as DM
N = 1

WorkingName = "TriStackedRingODEStatic"
WorkingName = "TriStackedRingODEStatic"

WorkingName = " TriStackedRingODELargeSynDeg28e-04T"
#WorkingName = "TriStackedRingODESlowSynDeg"
#WorkingName = "TriStackedRingODE1GenPerHour"
#WorkingName = "TriStackedRingODE2GenPerHour"

#File = h5py.File(".\\Data\\{}{}.h5".format(N, WorkingName), "r")
import importlib
if N == 1:
	CD = importlib.import_module("Data.{}{}".format("",WorkingName))
	CD.CDImportName = "Data.{}{}".format("",WorkingName)
else:
	CD = importlib.import_module("Data.{}{}".format(N-1,WorkingName))
	CD.CDImportName = "Data.{}{}".format(N-1,WorkingName)
	
def FindFiles(CD):
	Files = []
	for Filename in os.listdir(CD.LocalDataDir):
		Files.append(Filename)
	return Files

from ODERun import ODEFunc
Func = ODEFunc()
Func.GetSDataFromFiles()
Species = Func.Species

#CompleteFilename = DM.CreateNewComplete(CD.CompleteDataDir, CD)
#CompleteFilename = ".\\Data\\1TriStackedRingODEStatic.h5"

KD1IndiciesOut = []
KD2IndiciesOut = []
ResultsOut = []
ErrorsOut = []
RuntimesOut = []
from collections import OrderedDict
Max = OrderedDict()
Sum = OrderedDict()
#Max = {}
#Sum = {}

BinSpecies = [str(bin(int(Species.keys()[i]))) for i in range(CD.SpeciesCount)]
ProteinsInSpecies = [sum(x == "1" for x in BinSpecies[i]) for i in range(CD.SpeciesCount)]
NowGood = []
StillBad = []

for CurrentFile in FindFiles(CD):
	#LocalDir = "C:\Users\Koan\Dropbox\Proteasome Project\VSProjects\Cluster\Data\Local\\"
	LocalFilename = "{}{}".format(CD.LocalDataDir , CurrentFile)
#	print LocalFilename
	#DM.ImportLocalToComplete(CD, LocalFilename, CompleteFilename, ProteinsInSpecies)
	with h5py.File(LocalFilename, 'r') as Current:
		T = Current["Results"][:]
		CResults = np.array(Current["Results"][0])
		DM.ConvertToFraction(CD, CResults, ProteinsInSpecies)
		Index = "{},{}".format(Current["KD1Indices"][0], Current["KD2Indices"][0])
		#print Index
		try:
			Max[Index] = np.fromiter(flatten([Max[Index],np.max(CResults)]), np.float)
			#print "Appending"
		except:
			Max[Index] = [np.max(CResults),]
		try:
			Sum[Index] = np.fromiter(flatten([Sum[Index],np.sum(CResults)]), np.float)
		except:
			Sum[Index] = [np.sum(CResults),]

#for key in Max.keys():
#	if len(Max[key])>1:
#		print "----------------------------------------------------------"
#		print key
#		print "Max: {}".format(Max[key])
#		print "Sum: {}".format(Sum[key])
		
		#if not Current["Errors"][0,0]:
		NowGood.append((Current["KD1Indices"][0], Current["KD2Indices"][0]))
		print "Good: {}".format(NowGood[-1])
		
		try:
			if not ((KD1IndiciesOut[-1] == Current["KD1Indices"][0]) and (KD2IndiciesOut[-1] == Current["KD2Indices"][0])):
				KD1IndiciesOut.append(Current["KD1Indices"][0])
				KD2IndiciesOut.append(Current["KD2Indices"][0])
				ResultsOut.append(Current["Results"][:])
				ErrorsOut.append(Current["Errors"][:])
				RuntimesOut.append(Current["Runtimes"][:])
				CResults = np.array(Current["Results"][0])
				DM.ConvertToFraction(CD, CResults, ProteinsInSpecies)
				print "KD1: {} KD2 {}".format(KD1IndiciesOut[i],KD2IndiciesOut[i])
				print np.max(CResults)
				print np.sum(CResults[0])
				
			else:
				print "---------------------------------------------------------------"
				CResults = np.array(Current["Results"][0])
				DM.ConvertToFraction(CD, CResults, ProteinsInSpecies)
				print "KD1: {} KD2 {}".format(KD1IndiciesOut[i],KD2IndiciesOut[i])
				print np.max(CResults)
				print np.sum(CResults[0])
		except:
				KD1IndiciesOut.append(Current["KD1Indices"][0])
				KD2IndiciesOut.append(Current["KD2Indices"][0])
				ResultsOut.append(Current["Results"][:])
				ErrorsOut.append(Current["Errors"][:])
				RuntimesOut.append(Current["Runtimes"][:])

#	else:
#		StillBad.append((Current["KD1Indices"][0], Current["KD2Indices"][0]))
#		print "Bad: {}".format(StillBad[-1])

#NowGood = np.unique(NowGood)
#StillBad = np.unique(StillBad)

ToBeFixed = []
for SB in StillBad:
	Con = True
	for NG in NowGood:
		if (SB[0] == NG[0]) and (SB[1] == NG[1]):
			Con = False
	if Con:
		if len(ToBeFixed) == 0:
			ToBeFixed.append(SB)
		else:
			Temp = ToBeFixed[:]
			Unique = True
			for Entry in Temp:
				if ((SB[0] == Entry[0]) and (SB[1] == Entry[1])):
					Unique = False
			if Unique:
				ToBeFixed.append(SB)
					
#			if Unique:
#				ToBeFixed.append(SB)

print len(ToBeFixed)
print ToBeFixed

try:
	KD1Indicies = np.array(map(np.array,ToBeFixed))[:,0]
	KD2Indicies = np.array(map(np.array,ToBeFixed))[:,1]
except:
	pass

#[[KD1IndiciesOut[i],KD2IndiciesOut[i]] for i in range(len(KD1IndiciesOut))]
#ProteinsInSpecies = [sum(x == "1" for x in BinSpecies[i]) for i in range(CD.SpeciesCount)]

OutputShape = [
	len(KD1IndiciesOut),
	len(KD2IndiciesOut),
	len(CD.A0Range),
	len(CD.KeepTimePoints),
	CD.SpeciesCount,
	]

with h5py.File(".\\Data\\1TriStackedRingODEStatic.h5", "r+") as ErrorFixFile:
	for i in range(len(KD1IndiciesOut)):
		CResults = np.array(ResultsOut[i][0])
		DM.ConvertToFraction(CD, CResults, ProteinsInSpecies)
#		print "KD1: {} KD2 {}".format(KD1IndiciesOut[i],KD2IndiciesOut[i])
#		print np.max(CResults)
#		print np.sum(CResults[0])
		
		ErrorFixFile["Results"][KD1IndiciesOut[i], KD2IndiciesOut[i]] = CResults[0]
		ErrorFixFile["Errors"][KD1IndiciesOut[i], KD2IndiciesOut[i]] = np.array(ErrorsOut[i], np.bool)
		ErrorFixFile["Runtimes"][KD1IndiciesOut[i], KD2IndiciesOut[i]] = np.array(RuntimesOut[i], np.float32)
		ErrorFixFile["HasData"][KD1IndiciesOut[i], KD2IndiciesOut[i]] = np.array(1)

#Current.close()
#print File["Errors"][:]


"""
FullOutArray = File["Results"][:]
Errors = File["Errors"][:]
HasData = File["HasData"][:]
Runtimes = File["Runtimes"][:]
A0Range = File["A0Range"][:]
KD1Range = File["KD1Range"][:]
KD2Range = File["KD2Range"][:]
TimePoints = CD.KeepTimePoints
SpeciesCount = CD.SpeciesCount
"""
