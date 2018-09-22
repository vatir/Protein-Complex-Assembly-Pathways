#!/home/vatir/.conda/envs/vpy-2.7.3/bin/python
# -*- coding: utf-8 -*-

"""
Created on Wed Jul 16 11:09:07 2014

@author: Koan
"""

import numpy as np
from scipy.integrate import odeint
from collections import OrderedDict
import itertools

class ODEFunc(object):
	def __init__(self, SpeciesCount = 12):
		self.SpeciesCount = SpeciesCount
		self.Conc_0 = np.zeros(SpeciesCount, np.float)
		self.CreateTempData()
		self.Params = {}
		
	def SetConc(self, Conc):
		self.Conc_0 = np.array(Conc[:], np.float)

	def CreateTempData(self):
		self.dydat = np.zeros(self.SpeciesCount,np.float)
		self.ydatM = np.zeros((self.SpeciesCount,)*2,np.float)
		
	def Der(self, ydat, Time_Initial):
		ydat = np.array(ydat,np.float)
		self.ydatM = np.outer(ydat,ydat)
		return np.sum(self.Coef1D*ydat,axis=1)+np.sum(np.sum(self.ydatM*self.Coef2D,axis=1),axis=1)+self.Synth
	
	def SetCoef(self, Coef1D, Coef2D):
		self.Coef1D = Coef1D
		self.Coef2D = Coef2D
		self.SpeciesRange = np.arange(len(self.Coef1D))
	
	def GetSDataFromFiles(self, filename="species3"):
		try:
			self.Species = OrderedDict()
			with open(filename,"r") as SpeciesFile:
				SpeciesList = SpeciesFile.read().strip().split("\n")
				for S in SpeciesList:
					L = S.split(";;")
					self.Species[int(L[0])] = L[1]
		
			self.Species = OrderedDict(sorted(self.Species.items(), key=lambda x: x[0]))
			self.SpeciesConvert = OrderedDict()
			for i, S in enumerate(map(int,self.Species.keys())):
				self.SpeciesConvert[int(S)] = i

			self.SpeciesCount = len(self.Species.keys())
			self.CreateTempData()
		except:
			return False

	def GetRDataFromFiles(self, filename = "react3"):
		try:
			self.Reactions = []
			with open(filename,"r") as ReactFile:
				ReactList = ReactFile.read().strip().split("\n")
				ReactList = map(lambda x: x.split("\t"), ReactList)
				for R in ReactList:
					R = map(lambda x: x.split(";;"), R)
					R = map(lambda x: x[0], R)
					R = map(lambda x: x.split(","), R)
					# Drop Representation
					R = list(itertools.chain(*R))
					R = map(int, R)
					R[0]=self.SpeciesConvert[R[0]]
					R[1]=self.SpeciesConvert[R[1]]
					if R[0] > R[1]:
						R0 = R[0]
						R1 = R[1]
						R[0] = R1
						R[1] = R0
					R[3]=self.SpeciesConvert[R[3]]
					self.Reactions.append(R)

				# Sort Reactions by Product, R1, R2 (R1 is always less than R2)
				SortList = map(list,list(np.array(self.Reactions)[:,(3,0,1)]))
				L = range(len(SortList))
				for i in L:
					SortList[i].append(i)
				SortOrder = np.array(sorted(SortList))[:,3]
				self.Reactions = map(list,list(np.array(self.Reactions)[SortOrder]))
				
				return True
		except Exception as e:
			print e
			return False
	def GenCoef(
						self,
						Kd1,
						Kd2,
						Kp,
						Delta,
						BondCoef = (-9.0/0.6),
						A0 = False,
						):
		"""
		If SynDeg is being used make sure to update Conc_0 before updating GenCoef
		"""
		if A0:
			self.SetInitialConc(A0)

		Coef2D = np.zeros((self.SpeciesCount,)*3,np.float)
		Coef1D = np.zeros((self.SpeciesCount,)*2,np.float)
		for R in self.Reactions:
			Coef2D[R[0],R[0],R[1]] -= R[2]
			Coef2D[R[1],R[0],R[1]] -= R[2]
			Coef2D[R[3],R[0],R[1]] += R[2]
			
			Keff = (Kd1**R[4])*(Kd2**R[5])*np.exp((R[4]+R[5]-1)*BondCoef)
			
			Coef1D[R[0],R[3]] += R[6]*Keff
			Coef1D[R[1],R[3]] += R[6]*Keff
			Coef1D[R[3],R[3]] -= R[6]*Keff
		
		for j in range(self.SpeciesCount):
			for i in range(self.SpeciesCount):
				Coef2D[j,i,i] = Coef2D[j,i,i]/2.0
		
		Coef1D = Coef1D*Kp
		Coef2D = Coef2D*Kp

		self.Synth = np.zeros(self.SpeciesCount, np.float)
		self.Synth[0] = self.Conc_0[0]*Delta # Synthesis
		for Species in range(self.SpeciesCount):
			Coef1D[Species,Species] -= Delta # Dilution Effect

		self.Kp = Kp
		self.Kd1 = Kd1
		self.Kd2 = Kd2
		self.SetCoef(Coef1D, Coef2D)
	
	def UpdateKp(self, Kp):
		self.Coef1D = Kp*self.Coef1D/self.Kp
		self.Coef2D = Kp*self.Coef2D/self.Kp
		self.Kp = Kp
		
	def UpdateKd1(self, Kd1):
		self.GenCoef(Kd1, self.Kd2,self.Kp)
		self.Kd1 = Kd1
		
	def UpdateKd2(self, Kd2):
		self.GenCoef(self.Kd1, Kd2,self.Kp)
		self.Kd2 = Kd2
	
	def SetTimeRange(self, Start, End, N = 1e5, IncludeZero = True, IncludePowersof10 = True):
		# Overkill method for finding what powers of 10 are missing
		Time = np.logspace(Start,End,num=N)

		if (0.0 not in Time) and IncludeZero:
			Time = np.hstack((0.0,Time))

		if IncludePowersof10:
			Add = 10.0**np.arange(Start,End+1)[
				np.array(
					[(10.0**i not in Time) for i in range(Start,End+1)]
					)
				]
			if len(Add) > 0:
				Time = np.sort(np.hstack([Time,Add]))
		self.Time = Time

	def SetInitialConc(self, A0):
			self.Conc_0 = np.zeros(self.SpeciesCount, np.float)
			self.Conc_0[0] = A0

	def RunODEINT(
							self,
							#rtol = 1e10, # Used 1e-8 for 2Static (1e-10 Normally)
							#atol = 1e-10, # Used 1e-20 for 2Static (1e-25 Normally)
							rtol = 1e-16, # Used 1e-8 for 2Static (1e-10 Normally)
							atol = 1e-20, # Used 1e-20 for 2Static (1e-25 Normally)
							ResultsOnly=True,
							ReturnOutMessage=False,
							ReturnAll=False,
							mxstep=int(5e8),
							mxordn=int(5e7),
							mxords=int(5e7),
							hmin = 0,
							hmax = 0,
							):
		self.LastResults = odeint(
												self.Der,
												self.Conc_0,
												self.Time,
												full_output=True,
												rtol=rtol,
												atol=atol,
												mxstep = mxstep,
												mxordn  = mxordn ,
												mxords  = mxords ,
												hmin = hmin,
												hmax = hmax
												)
		if ReturnAll:
			return self.LastResults[0], self.LastResults[1]
		if ReturnOutMessage:
			return self.LastResults[0], self.LastResults[1]["message"]
		if ResultsOnly:
			return self.LastResults[0]
		return self.LastResults
		
if __name__ == "__main__":
	pass
