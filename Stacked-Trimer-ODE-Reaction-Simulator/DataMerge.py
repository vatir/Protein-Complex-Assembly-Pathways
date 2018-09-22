#!/home/vatir/.conda/envs/vpy-2.7.3/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import h5py

#try:
#	DataModule = sys.argv[1]
#except:
#	print "No CLI Args found!"
#	DataModule = "CurrentJobData"
import importlib

N = 1
Name = "STV2SynDeg28e-04T"
Name = "STV2StaticT"
Name = "NoBondsInPosStaticT"
Name = "NoBondsInPosSynDeg28e-04T"
Name = "NoPosSynDeg28e-04T"
Name = "NoPosStaticT"


Name = "STV3StaticT"
#Name = "STV3SynDeg14e-04T"
#Name = "STV3SynDeg28e-04T"
#Name = "STV3SynDeg56e-04T"
#Name = "STV3SynDeg17e-05T"
#Name = "STV3SynDeg35e-05T"
#Name = "STV3SynDeg69e-05T"
#Name = "STV3SynDeg27e-07T"
#Name = "STV3SynDeg54e-07T"
#Name = "STV3SynDeg11e-06T"
#Name = "STV3SynDeg22e-06T"
#Name = "STV3SynDeg43e-06T"
#Name = "STV3SynDeg87e-06T"
#Name = "STV3SynDeg32e-08T"
#Name = "STV3SynDeg10e-09T"
#Name = "STV3SynDeg10e-10T"
#Name = "STV3SynDeg10e-11T"
#Name = "STV3SynDeg10e-13T"
Name = "STV3SynDeg10e-14T"
#Name = "STV3SynDeg10e-15T"
#Name = "STV3SynDeg10e-16T"


#Name = "NoDisSynDeg28e-04T"
#Name = "LargeStaticT"

#Name = "LargeSynDeg15e-12T"
#Name = "LargeSynDeg15e-16T"
#Name = "LargeSynDeg15e-50T"

#Name = "LargeSynDeg12e-05T"
#Name = "LargeSynDeg30e-06T"
#Name = "LargeSynDeg59e-06T"
#Name = "LargeSynDeg28e-05T"
#Name = "LargeSynDeg90e-06T"

CD = importlib.import_module("Data.{}".format(Name))
#CD = importlib.import_module("Data.TriStackedRingODEStatic")


if N == 1:
	DataName = "{}{}".format("", Name)
else:
	DataName = "{}{}".format(N-1, Name)
CD = importlib.import_module("Data.{}".format(DataName))

DataName = "{}{}".format(N, Name)
CD.DataFilenamePrefix = DataName
CD.MainDataFilename = DataName
CD.DataFilenamePrefix = DataName

import DataManager as DM

def FindFiles(CD):
	Files = []
	for Filename in os.listdir(CD.LocalDataDir):
		try:
			if Filename[:len(CD.DataFilenamePrefix)]==CD.DataFilenamePrefix:
				Files.append(Filename)
				print "Found Filename: {}".format(Filename)
		except:
			pass
	return Files

from ODERun import ODEFunc
Func = ODEFunc()
Func.GetSDataFromFiles()
Species = Func.Species

if __name__ == '__main__':
	CompleteFilename = DM.CreateNewComplete(CD.CompleteDataDir, CD)
	#CompleteFilename = ".\\Data\\1TriStackedRingODEStatic.h5"

	BinSpecies = [str(bin(int(Species.keys()[i]))) for i in range(CD.SpeciesCount)]
	ProteinsInSpecies = [sum(x == "1" for x in BinSpecies[i]) for i in range(CD.SpeciesCount)]

	print CompleteFilename
	for CurrentFile in FindFiles(CD):
		LocalFilename = "{}{}".format(CD.LocalDataDir, CurrentFile)

		print "{}{}".format(CD.LocalDataDir, CurrentFile)
		DM.ImportLocalToComplete(CD, LocalFilename, CompleteFilename, ProteinsInSpecies)
