#!/home/vatir/.conda/envs/vpy-2.7.3/bin/python
# -*- coding: utf-8 -*-

"""
Created on Wed Jul 16 11:09:07 2014

@author: Koan
"""

from platform import system
import numpy as np
from subprocess import check_output, STDOUT, call
from tempfile import NamedTemporaryFile


class JobHolder(object):
	def __init__(self, *args, **kwargs):
		self.JobID = False	# False until submitted then set to Job ID
		
		self.SetDefault()

		return super(JobHolder, self).__init__(*args, **kwargs)

	def SetDefault(self):
		self.NoEmail
		self.User = "vatir"
		self.JobName = "Default"
		self.JobType = "Base"
		self.NodeCount = 1
		self.ProcsPerNode = 1
		#self.MaxMemmory = str(4500) + "m"
		self.MaxMemmory = str(2) + "g"
		self.MaxWallTime = "6:00:00" # sixhour limit
		#self.MaxWallTime = "60:00:00:00" # compbio limit
		self.Hardware = ":intel" # Always ask for intel hardware
		self.QueueName = "sixhour" # Ask for the upto 6 hour walltime queue
		#self.QueueName = "compbio" # Default Queue
		# self.Hardware = "" # Any hardware will do
		self.Shell = "/bin/bash"
		#self.TempDir = "/data/koanb/temp/"
		self.LogDir = "/panfs/pfs.local/work/compbio/vatir/PANFS-DATA-KOANB/LogsPBS/"
		self.TempDir = self.LogDir
		self.TskArrayNum = 10 # Default Value
		self.qsub = "/opt/moab/bin/msub"
		self.echo = "/bin/cat"

	def SetPBSBlock(self):
		self.PBSBlock = []		# Text of the PBS Job
		self.PBSBlock.append("#MSUB -N {JobName}".format(JobName=self.JobName))
		self.PBSBlock.append("#MSUB -q {QueueName}".format(QueueName=self.QueueName))
		self.PBSBlock.append("#MSUB -l nodes={nodes}:ppn={ppn}{hw},mem={mem},walltime={wall}".format(
			nodes=self.NodeCount,
			ppn=self.ProcsPerNode,
			mem=self.MaxMemmory,
			wall=self.MaxWallTime,
			hw=self.Hardware
			))
		if self.email:
			self.PBSBlock.append("#MSUB -M {email}".format(email=self.email))
			self.PBSBlock.append("#MSUB -m abe")
		else:
			self.PBSBlock.append("#MSUB -m n")
		self.PBSBlock.append("#MSUB -S {}".format(self.Shell))
		self.PBSBlock.append("#MSUB -d {}".format(self.TempDir))

		FilenameBase = "JobName-{}-JobType-{}".format(self.JobName, self.JobType)
		self.OutFile = "{}.out".format(FilenameBase)
		self.ErrFile = "{}.err".format(FilenameBase)

		self.PBSBlock.append("#MSUB -o {}".format(self.LogDir+self.OutFile))
		self.PBSBlock.append("#MSUB -e {}".format(self.LogDir+self.ErrFile))
		self.PBSBlock.append("#MSUB -t {}".format(self.TskArrayNum))
		

	def SetCodeBlock(self, Code):
		self.PBSBlock.append("echo \"Script Started at:\"")
		self.PBSBlock.append("date")
		for Line in Code:
			self.PBSBlock.append(Line)
		self.PBSBlock.append("echo \"Script Ended at:\"")
		self.PBSBlock.append("date")

	def GenScipt(self, Code):
		self.SetPBSBlock()
		self.SetCodeBlock(Code)

	def SubmitScipt(self, Code, EchoOutput = False):
		self.GenScipt(Code)

		Script = NamedTemporaryFile() # Filelike Object to hold the PBS script
		for Line in self.PBSBlock:
			Script.write(Line+"\n")
		#Script.flush()
		Script.seek(0)
		if EchoOutput:
			Command = self.echo
		else:
			Command = self.qsub

		if system() == "Windows":
			Command = "cat"

		self.JobID = check_output([Command,Script.name], stderr = STDOUT)
		Script.close()

		return self.JobID

#	@property
#	def Email(self):
#		self.email = "kbriggs@trileo.net"

	@property
	def NoEmail(self):
		self.email = False



def FindNextArrayNum(Username):
	try:
		CurrentMax = 0
		CurrentJobs = check_output(['/usr/local/bin/qstat',"-u {}".format(Username), "-t"])
		if CurrentJobs != '':
			for Line in CurrentJobs.split("\n")[5:-1]:
				JobID = Line.split()[0]
				ArrayID = int(JobID.split("[")[1].split("]")[0])
				if (ArrayID != '') and (CurrentMax < ArrayID):
					CurrentMax = ArrayID
		return CurrentMax + 1
	except Exception as e:
		print "Current Job Query Failed!!! (Interactive job probably in queue)"
		print e
		return 1


if __name__ == '__main__':
	Worker = JobHolder()
	#Code=["sleep 10","date","sleep 20"]
	
	import DataManager as DM

	import CurrentJobData as CJB

	call(["cp", "CurrentJobData.py", CJB.CDFileName])

	import importlib
	CD = importlib.import_module(CJB.CDImportName)

	
	#Worker.JobName = "Grid"
	#Code.append("python Worker.py {} {} {} {}".format(KD1Indices[0], KD1Indices[-1], KD2Indices[0], KD2Indices[-1]))
	#

	#print Worker.SubmitScipt(Code, EchoOutput = True)
	#for i in range(len(CD.KD1Range)):
	KD1List = []
	KD2List = []

	KD1Len = len(CD.KD1Range)
	KD2Len = len(CD.KD2Range)

	len(CD.KD2Range)
	from Common import Chunks
	for Chunk in Chunks(range(KD1Len), CD.JobCount):
		KD1List.append([Chunk[0],Chunk[-1]])
	if CD.JobCount > KD1Len:
		for Chunk in Chunks(range(KD2Len), int(CD.JobCount/KD1Len)):
			KD2List.append([Chunk[0],Chunk[-1]])
	else:
		KD2List = [[0,len(CD.KD2Range)-1]]

	#for KD1 in KD1List:
	#	for KD2 in KD2List:
	#		print "{} {} {} {}".format(KD1[0], KD1[1], KD2[0], KD2[1])

	#for KD1 in KD1List[:1]:
	#	for KD2 in KD2List[:1]:
	Worker.TskArrayNum = FindNextArrayNum(Worker.User)
	Worker.MaxMemmory = CD.JobMem

	i = 0
	for KD1 in KD1List:
		for KD2 in KD2List:
	#for KD1 in KD1List[:2]:
	#	for KD2 in KD2List[:2]:
			Code = []
			Code.append("cd /panfs/pfs.local/work/compbio/vatir/PANFS-DATA-KOANB/Current")
			Code.append("python Worker.py {} {} {} {} {}".format(CJB.CDImportName, KD1[0], KD1[1], KD2[0], KD2[1]))
			Worker.JobName = "{}:{}:{}:{}:{}".format(CJB.CDName, KD1[0], KD1[1], KD2[0], KD2[1])
			print Worker.JobName
			print Worker.SubmitScipt(Code, EchoOutput = False)
			print i
			i = i + 1

	print len(KD1List)*len(KD2List)


	#LocalFilename = DM.CreateNewLocal(DM.LocalDataDir, KD1Indices, KD2Indices)
	#DM.ImportLocalToComplete(LocalFilename, CompleteFilename)
