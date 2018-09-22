# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 16:41:09 2014

@author: Koan
"""
#import cStringIO.StringIO as StringIO
from lxml import etree
						
if __name__ == '__main__':
	import sys
	Dir = "E:\Proteasome Project Data\PisaWebData\\"
	Name = "6merHomoAll.Assemblies.Stripped.xml"
	Filename = Dir+Name
	TargetNum = 10
	
	try:
		Filename = str(sys.argv[1])
		TargetNum = int(sys.argv[2])
	except:
		print "No CLI Args found! Exiting ..."
#		quit()
	
	R350ContainingCount = 0
	EntryData = []
	EntryWhiteSpaceCount = 0
	CurrentEntryNum = 0

	depth = 0
	prefix_width = 8
	prefix_dots = '.' * prefix_width
	line_template = '{prefix:<0.{prefix_len}}{event:<8}{suffix:<{suffix_len}} {node.tag:<12} {node_id}'
	
	with open(Filename, 'r') as File:
		# get an iterable
		#context = etree.iterparse(File, events=("start", "end"))
		context = etree.iterparse(File, tag=("pdb_code"))
		
		# turn it into an iterator
		context = iter(context)
		
		# get the root element
		event, root = context.next()
	
		for (event, node) in context:
#			if event == 'end':
#				depth -= 1
			#print event
			print node.text
#			prefix_len = depth * 2
#			if depth < 2:
#				print line_template.format(prefix=prefix_dots,
#													    prefix_len=prefix_len,
#				                               suffix='',
#				                               suffix_len=(prefix_width - prefix_len),
#				                               node=node,
#				                               node_id=id(node),
#				                               event=event,
#				                               )
			    
#			if event == 'start':
#				depth += 1
			
			node.clear()
		root.clear()
#	XMLTemp = StringIO()
#	
#	#XMLTemp.write('<?xml version="1.0"?>'+"\n")
#	#XMLTemp.write('<data>'+"\n")
#	for Line in EntryData:
#		XMLTemp.write(Line.strip())
#	#XMLTemp.write('</data>'+"\n")
#	XMLTemp.seek(0)
#	Doc = ET.parse(XMLTemp)
##	for i in Doc.:
##		k = i
##		print k
#	XMLTemp.seek(0)
#	#print XMLTemp.read()
