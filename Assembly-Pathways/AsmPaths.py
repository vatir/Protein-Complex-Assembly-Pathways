from collections import OrderedDict, deque
import numpy as np
import itertools
import matplotlib.pylab as plt

class ImportedData(object):
	def __init__(self, SpeciesCount = 12):
		self.SpeciesCount = SpeciesCount
		self.GetSDataFromFiles()
		self.GetRDataFromFiles()

	def GetSDataFromFiles(self, filename="species3"):
		try:
			self.Species = OrderedDict()
			with open(filename,"r") as SpeciesFile:
				SpeciesList = SpeciesFile.read().strip().split("\n")
				for S in SpeciesList:
					L = S.split(";;")
					self.Species[L[0]] = L[1]
		
			self.SpeciesConvert = OrderedDict()
			for i, S in enumerate(np.sort(map(int,self.Species.keys()))):
				self.SpeciesConvert[S] = i

			self.SpeciesCount = len(self.Species.keys())
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
					R[3]=self.SpeciesConvert[R[3]]
					self.Reactions.append(R)
		except Exception as e:
			print e
			return False

def SpeciesToIndex(S, Species):
	return np.where(np.array(Species.values())==S)[0]

def SpeciesToName(S, Species):
	return Species.values()[S]

#def UniqueNodeName(Graph, NodeName):
#	Max = 1
#	for Node in Graph.nodes():
#		N = Node.split(" ")
#		try:
#			if len(N) in (2,3) and int(N[-1]) > Max:
#				Max = N[-1] + 1
#		except:
#			pass
#	return "{} {}".format(NodeName, Max)

def UniqueNodeName(Graph, NodeName):
	return "{} {}".format(NodeName,len(Graph.nodes()))

def FillOut(Graph, Reactions, Species, TODO, Dual=True):
	while len(TODO) > 0:
		C = TODO.pop()
		LastNodeName = C[-1]
		for R in Reactions:
			RV0 = (R[3] == C[0])
			if Dual:
				RV1 = (R[3] == C[1])
			else:
				RV1 = False
			if (RV0 or RV1):
				R1 = Species.values()[R[0]]
				if Dual:
					R2 = Species.values()[R[1]]
					NewNodeName = UniqueNodeName(Graph, "{} {}".format(R1, R2))
					Graph.add_node(NewNodeName, {"R1" : R[0],"R2" : R[1]})
				else:
					NewNodeName = UniqueNodeName(Graph, "{}".format(R1))
					Graph.add_node(NewNodeName, {"R1" : R[0]})
				Graph.add_edge(LastNodeName, NewNodeName)
				if not ((R[0] == 0) and (R[1] == 0)):
					TODO.append((R[0], R[1], NewNodeName))

def CheckReaction(R, P, Reactions, Species, R2=False):
	for Rec in Reactions:
		print SpeciesToIndex(R, Species)
		R11 = (SpeciesToIndex(R, Species) == Rec[0])
		R12 = (SpeciesToIndex(R, Species) == Rec[1])
		#R1 = (SpeciesToIndex(R, Species) == Rec[0]) and Rec[0] != "1,0,0;0,0,0"
		#R2 = (SpeciesToIndex(R, Species) == Rec[1]) and Rec[1] != "1,0,0;0,0,0"
		if R2:
			R21 = (SpeciesToIndex(R2, Species) == Rec[0])
			R22 = (SpeciesToIndex(R2, Species) == Rec[1])
			if (R11 or R12) and (R21 or R22) and Rec[3] == SpeciesToIndex(P, Species):
				return True
		else:
			if (R11 or R12) and Rec[3] == SpeciesToIndex(P, Species):
				return True
		return False

def BSplit(Nodes, Reactions, Species):
	print Nodes
	BGraph = nx.DiGraph()
	TODO = deque()
	#NewNodes = deque(maxlen=3)
	NewNodes = deque()
	BGraph.add_node(Nodes[0], data=True)
	NewNodes.append(Nodes[0])
	Nodes = Nodes[1:]
	Nodes.reverse()
	for Node in Nodes:
		N = Node.split(" ")
		if len(N) in (2,3):
			N1 = N[0]
			N2 = N[1]
			#if CheckReaction(N1, LastNode, Reactions, Species, R2=N2):
			#	BGraph.add_edge(LastNode, N1)
			#	BGraph.add_edge(LastNode, N2)
			try:
				LN = [NewNodes[-1]]
				LN.append(NewNodes[-2])
			except:
				LN = [NewNodes[-1]]
			
			for L in LN:
				if CheckReaction(N1, L.split(" ")[0], Reactions, Species, R2=N2):
					for Nx in (N1, N2):
						NewNodes.append(Nx)
						NewNodeName = UniqueNodeName(BGraph, "{}".format(Nx))
						BGraph.add_node(NewNodeName)
						BGraph.add_edge(L, NewNodeName)
				else:
					for Nx in (N1, N2):
						NewNodes.append(Nx)
						NewNodeName = UniqueNodeName(BGraph, "{}".format(Nx))
						BGraph.add_node(NewNodeName)
						for PN in BGraph.nodes():
							PN = PN.split(" ")[0]
							if CheckReaction(Nx, PN, Reactions, Species):
								BGraph.add_edge(PN, Nx)
				#
				#BGraph.add_node(NewNodeName, {"R1" : Nx})
				#TODO.append((Nx, NewNodeName))
				#if len(NewNodes) > 2:
				#	BGraph.add_edge(NewNodes[1],NewNodes[0])

			#BGraph.add_edge(NewNodes[2],NewNodes[0])
			#BGraph.add_edge(NewNodes[2],NewNodes[1])
	#NewNodes.reverse()

	#for Nx in NewNodes:
	#	NewNodeName = UniqueNodeName(BGraph, "{}".format(Nx))
	#	BGraph.add_node(NewNodeName, {"R1" : Nx})
	#	for Node in BGraph.nodes():
	#		if CheckReaction(Node, Nx, Reactions, Species):
	#			BGraph.add_edge(Nx,Node)
	print NewNodes
	#FillOut(BGraph, Reactions, Species, TODO, Dual=False)
	return BGraph

def Plot(self, with_labels=True, Show=True):
	Fig = plt.figure()
	Axes = Fig.add_subplot(111)
	#print self.nodes()
	nx.write_dot(self,'test.dot')
	#Fig.set_title("Tree Plot")
	pos=nx.graphviz_layout(self,prog='dot')
	Axes.plot(nx.draw(self,pos,with_labels=with_labels,arrows=False))
	if Show:
		plt.show()
	return Fig

if __name__ == "__main__":
	import networkx as nx
	Asm = nx.DiGraph()

	D = ImportedData()


	RootSpeciesNum = int(max(map(int,D.Species.keys())))
	RootSpeciesName = D.Species[str(RootSpeciesNum)]
	RootSpeciesCode = D.SpeciesConvert[RootSpeciesNum]
	
	TODO = deque()

	Asm.add_node(RootSpeciesName)
	#  0  1 2 3  4  5 6
	# R1 R2 F P B1 B2 R
	for R in D.Reactions:
		if R[3] == RootSpeciesCode:
			R1 = D.Species.values()[R[0]]
			R2 = D.Species.values()[R[1]]
			NewNodeName = "{} {}".format(R1, R2)
			#NewNodeName = UniqueNodeName(Asm, NewNodeName)
			Asm.add_node(NewNodeName, {"R1" : R[0],"R2" : R[1]})
			Asm.add_edge(RootSpeciesName, NewNodeName)
			if not ((R[0] == 0) and (R[1] == 0)):
				TODO.append((R[0], R[1], NewNodeName))

	FillOut(Asm, D.Reactions, D.Species, TODO)

	for Node in Asm.nodes():
		Entries = Node.split(" ")
		if len(Entries) == 3:
			N1, N2, N = Entries
			if N1 == D.Species.values()[0] and N2 == D.Species.values()[0]:
				TODO.append(Node)

	BTrees = []
	for Node in [TODO[0]]:
		Nodes = nx.shortest_path(Asm, RootSpeciesName, Node)
		BTrees.append(BSplit(Nodes, D.Reactions, D.Species))
	
	Plot(BTrees[0])
	#Plot(Asm)
	pass