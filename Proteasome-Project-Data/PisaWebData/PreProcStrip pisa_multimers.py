# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 16:41:09 2014

@author: Koan
"""
						
if __name__ == '__main__':
	import sys
	Dir = "E:\Proteasome Project Data\PisaWebData\\"
	InName = "6merHomoAll.Assemblies.xml"
	OutName = "6merHomoAll.Assemblies.Stripped.xml"
	
	StripList = (
							"<pisa_multimers>\n",
							"</pisa_multimers>\n",
							"  <status>Ok</status>\n",
							)
	with open(Dir+InName, 'r') as In:
		with open(Dir+OutName, 'w') as Out:
			Out.write("<data>\n")
			for Line in In.readlines():
				if not (Line in StripList):
					Out.write(Line)
			Out.write("</data>\n")
		