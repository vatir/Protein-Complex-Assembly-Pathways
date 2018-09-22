# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 16:41:09 2014

@author: Koan
"""

#import re
#Dir = "Assemblies"
#Filename = "6merHomoAll.xml"
#i = 0
#n = True
#with open(".\\Data\\{}\\{}".format(Dir, Filename), "r") as In:
#	with open(".\\Data\\{}\\Modded{}".format(Dir, Filename), "w") as Out:
#		for Line in In.readlines():
#			OutLine, Reps = re.subn("pisa_multimers", "pisa_multimers{}".format(i), Line)
#			if (Reps > 0) and (n):
#				n = False
#			elif (Reps > 0):
#				n = True
#				i += 1
#				print "Found: {}".format(i)
#			Out.write(OutLine)

#with open(".\\Data\\{}\\Modded{}".format(Dir, Filename), "r") as In:
#	for i in range(1):
#		print In.readline()


if __name__ == '__main__':
	import sys
	import numpy as np

	#Filename = "E:\Proteasome Project Data\PisaWebData\Assemblies\\6merHomoAll.xml"
	#Filename = "E:\Proteasome Project Data\PisaWebData\Interfaces\\6merHomoAll.dat"
	#TargetNum = 10

	try:
		Filename = str(sys.argv[1])
		TargetNum = int(sys.argv[2])
	except:
		print "No CLI Args found! Exiting ..."
		quit()

	R350ContainingCount = 0
	EntryData = []
	EntryWhiteSpaceCount = 0
	CurrentEntryNum = 0
	with open(Filename, 'r') as File:
		InData = False
		for Line in File:
			if ("pdb_entry" in Line) and not ("/pdb_entry" in Line):
				InData = True
				CurrentEntryNum += 1
				EntryWhiteSpaceCount = len(Line) - len(Line.lstrip())
			if ("/pdb_entry" in Line):
				InData = False
				if CurrentEntryNum == TargetNum:
					EntryData.append(Line.strip("\n")[EntryWhiteSpaceCount:])
					break

			if CurrentEntryNum == TargetNum:
				EntryData.append(Line.strip("\n")[EntryWhiteSpaceCount:])

	for Line in EntryData:
		print Line