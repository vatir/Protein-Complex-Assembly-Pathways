import numpy as np
def SetLogRange(Start, End, N = 1e5, IncludeZero = True, IncludePowersof10 = True):
	# Overkill method for finding what powers of 10 are missing
	R = np.logspace(Start,End,num=N)

	if (0.0 not in R) and IncludeZero:
		R = np.hstack((0.0,R))

	if IncludePowersof10:
		Add = 10.0**np.arange(Start,End+1)[
			np.array(
				[(10.0**i not in R) for i in range(Start,End+1)]
				)
			]
		if len(Add) > 0:
			R = np.sort(np.hstack([R,Add]))
	return R

import os.path
def CheckForImportFile(Dir, Filename, Ext):
	DataFilename = Dir+Filename+Ext
	if not os.path.isfile(DataFilename):
		return 0
	for n in range(1,100):
		DataFilename = Dir+"{}{}{}".format(n, Filename, Ext)
		if not os.path.isfile(DataFilename):
			return n

# --------------------------------------------------------
BatchName = "STV2"
#BatchName = "NoBondsInPos"
#BatchName = "NoDis"
#BatchName = "NoPos"
JobCount = 600
# --------------------------------------------------------
# odeint Parameters
rtol = 1e-10 # Used 1e-8 for 2Static (1e-10 Normally)
atol = 1e-20 # Used 1e-20 for 2Static (1e-25 Normally)
mxstep=int(5e8)
mxordn=int(5e7)
mxords=int(5e7)
hmin = 1e-50 # 0 means solver determined
hmax = 1e20 # 0 means solver determined

# --------------------------------------------------------
#Delta = 0.0 # No SynDeg
#Delta = 2.0/3600.0 # Half Hour Degredation Coeffecient (SynDeg if not equal to 0.0)
#Delta = 1.0/3600.0 # Once an Hour Degredation Coeffecient (SynDeg if not equal to 0.0)
#Delta = 0.01/3600.0 # Very Slow
#Delta = 1.0/600.0 # Once every 10 minutes

# Association Adjustment for Keff Kd Orig with no factor 2
Hours = 128
Delta = 1.0/(3600.0*float(Hours)) # Once a Hour

JobMem = str(16) + "g"
#JobMem = str(4) + "g"

Kp = 1e6
TimeLength = 1e2
TimeStart = -3
TimeStop = 12

SpeciesCount = 12
# --------------------------------------------------------
A0Range = SetLogRange(-12, -3, N = 1e2, IncludeZero=False, IncludePowersof10 = True)
KD1Range = SetLogRange(-12, -3, N = 1e2, IncludeZero=False, IncludePowersof10 = True)
KD2Range = SetLogRange(-12, -3, N = 1e2, IncludeZero=False, IncludePowersof10 = True)
#KeepTimePoints = np.array([10.0**i for i in range(TimeStart,TimeStop+1)])  # Only powers of 10
KeepTimePoints = SetLogRange(TimeStart, TimeStop, N = 1e2, IncludeZero=False, IncludePowersof10 = True)  # Keep same length as A0 for Tau calculations

# --------------------------------------------------------

VarValuesIndicies = np.array(["KD1", "KD2", "A0", "Time", "SpeciesIndex"])

CompleteDataDir = "/panfs/pfs.local/work/compbio/vatir/PANFS-DATA-KOANB/Current/Data/"
LocalDataDir = "/panfs/pfs.local/work/compbio/vatir/PANFS-DATA-KOANB/Current/Data/Local/"

if Delta == 0.0:
	BatchType = "StaticT"
else:
	BatchType = "SynDeg{:.1e}T".format(Delta).replace(".","")

from platform import system
if system() == "Windows":
	CompleteDataDir = ".\\Data\\"
	LocalDataDir = ".\\Data\\Local\\"

CDName = "{}{}".format(BatchName, BatchType)
UniquePrefix = ""
Index = CheckForImportFile(CompleteDataDir, CDName, ".py")
if Index != 0:
	UniquePrefix = str(Index)
CDName = "{}{}".format(UniquePrefix, CDName)

CDFileName = "{}{}.py".format(CompleteDataDir, CDName)
CDImportName = "Data.{}".format(CDName)
MainDataFilename = CDName
DataFilenamePrefix = CDName
