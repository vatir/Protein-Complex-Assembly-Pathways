#!/home/vatir/.conda/envs/vpy-2.7.3/bin/python
# -*- coding: utf-8 -*-

import sys
import numpy as np

try:
	args = map(np.int,sys.argv[2:])
	DataModule = sys.argv[1]
except:
	print "No CLI Args found!"
	args = []
	DataModule = "CurrentJobData"

#if len(args) == 4:
#	if args[0]==args[1]:
#		KD1Indices = np.array([args[0],])
#	else:
#		KD1Indices = np.arange(args[0],args[1])
#	if args[2]==args[3]:
#		KD2Indices = np.array([args[2],])
#	else:
#		KD2Indices = np.arange(args[2],args[3])
if len(args) == 4:
	KD1Indices = np.arange(args[0],args[1]+1)
	KD2Indices = np.arange(args[2],args[3]+1)
elif len(args) == 0:
	KD1Indices = np.array([80,])
	#KD2Indices = np.array([75,])
	#KD1Indices = np.array([80,])
	KD2Indices = np.arange(50,100)
elif len(args) == 5:
	if args[0]==args[1]:
		KD1Indices = np.array([args[0],])
	else:
		KD1Indices = np.arange(args[0],args[1])
	if args[2]==args[3]:
		KD2Indices = np.array([args[2],])
	else:
		KD2Indices = np.arange(args[2],args[3])
	PrintAll = args[4]
else:
	KD1Indices = args[0]
	KD2Indices = args[1]

import importlib
CD = importlib.import_module(DataModule)

import DataManager as DM
from ODERun import ODEFunc
from time import time
import h5py

Func = ODEFunc()
Func.GetSDataFromFiles()
Func.GetRDataFromFiles()
Func.SetTimeRange(CD.TimeStart, CD.TimeStop, CD.TimeLength)
FuncType = ODEFunc.SetInitialConc

Failed = False
def ChangeAndRun(Var):
	ODEFunc.GenCoef(
				 Func,
				 Kd1 = Var[0],
				 Kd2 = Var[1],
				 Kp= CD.Kp,
				 A0=Var[2],
				 Delta = Var[3]
				 )
	Out, Message = Func.RunODEINT(
		ReturnOutMessage=True,
		rtol = CD.rtol,
		atol = CD.atol,
		mxstep = CD.mxstep,
		mxordn = CD.mxordn,
		mxords = CD.mxords,
		hmin = CD.hmin,
		hmax = CD.hmax
		)
	#global PrintAll
	#if PrintAll:
	#	Out, AllOutput = Func.RunODEINT(ReturnAll=True)
	if Message != 'Integration successful.':
		global Failed
		Failed = True
		#if PrintAll: 
		#	print AllOutput.keys()
		#	print AllOutput
		print "KD1: %e KD2: %e A0: %e Out: \"%s\"" % (Var[0], Var[1], Var[2], Message)
	return Out # All

if __name__ == '__main__':

	#FullOutArray = np.zeros([
	#						len(CD.KD1Range),
	#						len(CD.KD2Range),
	#						len(CD.A0Range),
	#						len(CD.KeepTimePoints),
	#						CD.SpeciesCount,
	#						], np.float)

	Runtimes = list()

	LocalFilename = DM.CreateNewLocal(CD.LocalDataDir, KD1Indices, KD2Indices, CD)
	LocalFile = h5py.File(LocalFilename, "r+")

	#["KD1", "KD2", "A0", "Time", "SpeciesIndex"]
	for KD1LinearIndex, KD1Index in enumerate(KD1Indices):
		for KD2LinearIndex, KD2Index in enumerate(KD2Indices):
			StartTime = time()
			Failed = False

			import CurrentJobData as CDBase

			KD1Value = CD.KD1Range[KD1Index]
			KD2Value = CD.KD2Range[KD2Index]
			
			CurrentRun = np.transpose(
										np.vstack(
											[
											KD1Value*np.ones(len(CD.A0Range)),
											KD2Value*np.ones(len(CD.A0Range)),
											CD.A0Range,
											CD.Delta*np.ones(len(CD.A0Range)),
											]
										)
									)

			OutArray = np.array(map(ChangeAndRun, CurrentRun))
			
			Runtimes.append(np.float32(time() - StartTime))
			print "RunTime: %f" % (Runtimes[-1])

			KeepTimePointsIndicies = np.in1d(Func.Time.ravel(), CD.KeepTimePoints)
			
			#FullOutArray[KD1Index, KD2Index] = OutArray[:,KeepTimePointsIndicies, :]
			LocalFile['Results'][KD1LinearIndex,KD2LinearIndex] = OutArray[:,KeepTimePointsIndicies, :]
			LocalFile['Runtimes'][KD1LinearIndex,KD2LinearIndex] = Runtimes[-1]
			LocalFile['Errors'][KD1LinearIndex,KD2LinearIndex] = Failed

			del OutArray

			print "KD1: %.10e" %(KD1Value)
			print "KD2: %.10e" %(KD2Value)
			print "-----------------------------------"
	
	print "Total RunTime: %f" % (np.sum(Runtimes))
	print "KD1 Range: {}".format(KD1Indices)
	print "KD2 Range: {}".format(KD2Indices)

	LocalFile.close()
