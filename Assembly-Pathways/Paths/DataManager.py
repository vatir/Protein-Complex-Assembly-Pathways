#!/data/koanb/vpy-2.7.3/bin/python2.7
# -*- coding: utf-8 -*-

"""
Created on Wed Jul 16 11:09:07 2014

@author: Koan
"""

from CurrentJobData import SetLogRange
import numpy as np
try:
	import h5py # hdf5 not installed on login nodes but it is on the workers
except:
	pass
import os

def NewFile(DataDir, Filename = "Complete"):
	Ext = ".h5"
	CompleteDataFilename = DataDir+"{}{}".format(Filename, Ext)
	if not os.path.isfile(CompleteDataFilename):
		H5CompleteFile = h5py.File(CompleteDataFilename, "w")
		return H5CompleteFile
	for n in range(1,100):
		CompleteDataFilename = DataDir+"{}_{}{}".format(Filename,n, Ext)
		if not os.path.isfile(CompleteDataFilename):
			H5CompleteFile = h5py.File(CompleteDataFilename, "w")
			return H5CompleteFile 

def FindLastCompleteFile(DataDir, Filename = "Complete"):
	Ext = ".h5"
	CompleteDataFilename = DataDir+"{}{}".format(Filename, Ext)
	if not os.path.isfile(CompleteDataFilename):
		return CompleteDataFilename 
	for n in range(1,100):
		CompleteDataFilename = DataDir+"{}_{}{}".format(Filename,n, Ext)
		if not os.path.isfile(CompleteDataFilename):
			return DataDir+"{}_{}{}".format(Filename,n-1, Ext)


def CreateNewComplete(DataDir, CD): # CD Should be a CD module
	OutputShape = [
		len(CD.KD1Range),
		len(CD.KD2Range),
		len(CD.A0Range),
		len(CD.KeepTimePoints),
		CD.SpeciesCount,
		]

	File = NewFile(DataDir, CD.MainDataFilename)
	Filename = File.filename

	File['Results'] = np.zeros(OutputShape, dtype=np.float)
	File['HasData'] = np.zeros(OutputShape, dtype=np.bool)
	File['Runtimes'] = np.zeros(OutputShape[:2], dtype=np.float32)
	File['Errors'] = np.zeros(OutputShape[:2], dtype=np.bool)
	File['A0Range'] = CD.A0Range
	File['KD1Range'] = CD.KD1Range
	File['KD2Range'] = CD.KD2Range
	File['VarValuesIndicies'] = CD.VarValuesIndicies
	File['Delta'] = CD.Delta
	File['Kp'] = CD.Kp
	File['TimeLength'] = CD.TimeLength
	File['SpeciesCount'] = CD.SpeciesCount

	File.close()

	return Filename

def CreateNewLocal(DataDir, KD1Indices, KD2Indices, CD): # CD Should be a CD module
	OutputShape = [
		len(KD1Indices),
		len(KD2Indices),
		len(CD.A0Range),
		len(CD.KeepTimePoints),
		CD.SpeciesCount,
		]
	KD1Indices = np.array(KD1Indices)
	KD2Indices = np.array(KD2Indices)
	File = NewFile(DataDir, Filename = "{} {} {} {} {}".format(CD.DataFilenamePrefix, KD1Indices[0], KD1Indices[-1], KD2Indices[0], KD2Indices[-1]))
	Filename = File.filename
	
	File['KD1Indices'] = np.array(KD1Indices)
	File['KD2Indices'] = np.array(KD2Indices)
	File['Results'] = np.zeros(OutputShape, dtype=np.float)
	File['Runtimes'] = np.zeros(OutputShape[:2], dtype=np.float32)
	File['Errors'] = np.zeros(OutputShape[:2], dtype=np.bool)

	File.close()

	return Filename 

def ConvertToFraction(CD, Data, ProteinsInSpecies):
	for k in range(Data.shape[0]):
		for m in range(Data.shape[2]):
			for i in range(Data.shape[3]):
				try:
					Data[k,:,m,i] = Data[k,:,m,i]*ProteinsInSpecies[i]/CD.A0Range
				except:
					print "%i, %i, %i"%(k,m,i)

def ImportLocalToComplete(CD, Local, Complete, ProteinsInSpecies = False):
	"""
	Local: Filename
	Complete: Filename
	"""
	Local = h5py.File(Local, "r")
	Complete = h5py.File(Complete, "r+")

	KD1Indices = np.array(Local["KD1Indices"])
	KD2Indices = np.array(Local["KD2Indices"])

	for i, KD1 in enumerate(KD1Indices):
		Results = np.array(Local['Results'][i])
		if ProteinsInSpecies: ConvertToFraction(CD, Results, ProteinsInSpecies)
		Complete['Results'][KD1,KD2Indices] = Results
		Complete['HasData'][KD1,KD2Indices] = np.ones_like(Local['Results'][i],np.bool)
		Complete['Runtimes'][KD1,KD2Indices] = np.array(Local['Runtimes'][i])
		Complete['Errors'][KD1,KD2Indices] = np.array(Local['Errors'][i])

	Local.close()
	Complete.close()

if __name__ == '__main__':
	import CurrentJobData as CD
	#CompleteFilename = CreateNewComplete(CompleteDataDir)
	CompleteFilename = FindLastCompleteFile(CompleteDataDir)
	print CompleteFilename 
	CompleteFile = h5py.File(CompleteFilename)
	#print np.array(CompleteFile["A0Range"])
	CompleteFile.close()

	KD1Indices = [0,99]
	KD2Indices = np.arange(100)

	LocalFilename = CreateNewLocal(LocalDataDir, KD1Indices, KD2Indices)

	print LocalFilename 
	ImportLocalToComplete(LocalFilename, CompleteFilename)

	#CurrentRun = np.transpose(
	#						np.vstack(
	#							[
	#							KD1*np.ones(len(A0Range)),
	#							KD2*np.ones(len(A0Range)),
	#							A0Range,
	#							]
	#						)
	#					)

	#LongestIndex = max([
	#				len(A0Range),
	#				len(KD1Range),
	#				len(KD2Range)]
	#				)

	#np.save(f, VarValuesIndicies)
	#np.save(f, KD1Range)
	#np.save(f, KD2Range)
	#np.save(f, A0Range)
	#np.save(f, TimePoints)
	#np.save(f, Func.SpeciesConvert.keys())
	#np.save(f, Func.SpeciesConvert.values())
	#np.save(f, FullOutArray)
	#np.save(f, np.array([Failed,]))
