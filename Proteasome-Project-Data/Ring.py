import CurrentJobData as CD
import numpy as np
from scipy.integrate import odeint

class Ring3(object):
	"""description of class"""
	def __init__(self, KD = 1e-8):
		# (Monomer, Dimer, Trimer)
		self.Synth = 0.0
		self.Deg = 0.0
		self.Conc0 = np.zeros(3,dtype="f8")
		self.KD = KD
		self.KOn = CD.Kp
		self.KOff = self.KD*self.KOn
		self.Time = CD.KeepTimePoints
		self.A0Range = CD.A0Range

	def SetConc0(self, Conc0):
		self.Synth = Conc0*self.Deg
		self.Conc0[0] = Conc0

	def SetDeg(self, Deg):
		self.Deg = Deg
		if self.Deg != 0.0:
			self.Synth = self.Conc0[0]*self.Deg
		else:
			self.Synth = 0.0

	def SetKD(self, KD):
		self.KD = KD
		self.KOff = self.KD*self.KOn

	# For now NEVER change KOn
	#def SetKOn(self, KOn):
	#	self.KOn = KOn
	#	self.KOff = self.KD*self.KOn

	def Equations(self, Conc, Time0 = 0.0):
		Results = np.zeros(3,dtype="f8")
		C0 = 1.0
		Alpha = self.KOn
		Beta	= self.KOff
		Gamma = Alpha*C0*(self.KD**2.0)*np.exp(-9.0/0.6)
		Results[0] = 2.0*Beta*Conc[1]-2.0*Alpha*Conc[0]**2.0-Alpha*Conc[0]*Conc[1]+3.0*Gamma*Conc[2]+self.Synth-self.Deg*Conc[0]
		Results[1] = Alpha*Conc[0]**2.0-Beta*Conc[1]-Alpha*Conc[0]*Conc[1]+3.0*Gamma*Conc[2]-self.Deg*Conc[1]
		Results[2] = Alpha*Conc[0]*Conc[1]-3.0*Gamma*Conc[2]-self.Deg*Conc[2]
		return Results

	def RunODEINT(
							self,
							rtol = 1e-30,
							atol = 1e-15,
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
												self.Equations,
												self.Conc0,
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

	def GetA0Course(self, TimeIndex):
		#Index = np.argmax(self.Time == Time)
		try:
			return map(lambda x: self.GetTimeCourse(x)[TimeIndex], self.A0Range)
		except:
			print "Error in Ring A0 Data Plotting."
			return np.zeros_like(self.A0Range)

	def GetTimeCourse(self, A0, FinalFrac = True):
		self.SetConc0(A0)
		if FinalFrac:
			try:
				return self.RunODEINT()[:,-1]*3.0/self.Conc0[0]
			except:
				print "Error in Ring Time Data Plotting."
				return np.zeros_like(self.Time)
		else:
			return self.RunODEINT()

if __name__ == "__main__":
	Current = Ring3(KD = 10.0**-8.0)
	Current.SetConc0(10.0**-8.0)
	#Current.SetDeg(0.0)
	print Current.Conc0
	print Current.KOff
	import matplotlib.pyplot as plt
	fig = plt.figure(facecolor="w")
	axes = plt.subplot(111)
	axes.plot(Current.A0Range, Current.GetA0Course(1e8))
	axes.set_xscale('log')
	plt.show()
	#print Current.GetTimeCourse(1e-6)

	#for A0 in CD.A0Range:
	#	Current.SetConc0(A0)
		#print "A0 {} : {}".format(A0, Current.RunODEINT(ReturnOutMessage=True)[1])