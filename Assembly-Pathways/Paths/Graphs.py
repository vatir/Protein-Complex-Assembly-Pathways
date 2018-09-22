import networkx as nx
import networkx.algorithms.isomorphism as iso
from Helpers import GraphPlot, TimeStamp, TreePlot, GraphCompare
import itertools
from copy import copy
import numpy as np
from collections import deque

class GraphClass(nx.DiGraph, GraphPlot, object):
	"""
	Base Class for Networkx Graphs
	Note: DiGraphs are not fully tested with masking but seem to work at least for stacked rings.
	"""
	def __init__(self, Graph, data = None, **attr):
		super(GraphClass, self).__init__(data, **attr)

		self.CurrentlyMasked = False
		self.add_nodes_from(Graph.nodes(data=True))
		self.add_edges_from(Graph.edges(data=True))
		if issubclass(type(Graph), GraphClass):
			self.CurrentlyMasked = copy(Graph.CurrentlyMasked)
			if self.CurrentlyMasked:
				self._MaskedNodes = copy(Graph._MaskedNodes)

		self._OrigNodes = copy(Graph.nodes(data=True))
		self._OrigEdges = copy(Graph.edges(data=True))

	def TxtRep(self, Length = False, Canonical = True):
		if not Length:
			Length = (len(self._OrigNodes))
		S = list(Length*"0")
		for Node in self.nodes():
			T, N = Node.split("_")
			N = int(N)
			if "A" == T:
				N = N-1
			if "B" == T:
				N = N+(int(len(S)/2)-1)
			S[N] = T
		if Canonical:
			if "A" in S:
				Shift = np.argmax(np.array(S) != "0")
			else:
				Shift = np.argmax(np.array(S) != "0")-(int(len(S)/2))
			S1 = np.array(S)[:int(len(S)/2)]
			S2 = np.array(S)[int(len(S)/2):]
			S = np.array([np.roll(S1,-Shift), np.roll(S2,-Shift)]).flatten()
		return "".join(S)

	@property
	def CheckIfMasked(self):
		return self.CurrentlyMasked

	def NodeComplement(self, Nodes):
		return [Node for Node in self.nodes() if not (Node in Nodes)]

	@property
	def EdgeTypes(self):
		Types = dict()
		for Edge in self.edges(data=True):
			try:
				Types[Edge[2]["Type"]] += 1
			except:
				Types[Edge[2]["Type"]] = 1
		return Types
		
	def Similar(self, Graph, SubNodes = "All", InterfaceDict = False):
		Local = self
		if not (SubNodes == "All"):
			Local = self.SubCopy(SubNodes)
		return GraphCompare(Local, Graph)

	def Copy(self):
		New = GraphClass(self.copy())
		if self.CurrentlyMasked:
			New.CurrentlyMasked = True
			New._MaskedNodes = copy(self._MaskedNodes)
		return New

	def MaskNodes(self, Nodes):
		self.CurrentlyMasked = True
		self._MaskedNodes = Nodes
		self.remove_nodes_from(Nodes)
		pass

	def UnMaskNodes(self):
		if self.CurrentlyMasked:
			self.CurrentlyMasked = False
			self.add_nodes_from(self._OrigNodes)
			self.add_edges_from(self._OrigEdges)

	def FullNodes(self):
		# Return a copy with all the original nodes but not the masked edges
		Graph = self.Copy()
		Graph.add_nodes_from(self._MaskedNodes)
		self._MaskedNodes = []
		return Graph

	def SubCopy(self, IncludedNodes):
		# Returns a GraphClass object with only the nodes in IncludedNodes
		Graph = self.Copy()
		Graph.MaskNodes(self.NodeComplement(IncludedNodes))
		return Graph

	def DiGraphRep(self): # Plain nx.DiGraph version
		Graph = nx.DiGraph()
		for Node in self.nodes(data=True):
			Graph.add_node(*Node)
		for Edge in self.edges(data=True):
			Graph.add_edge(*Edge)
		return Graph

	@property
	def Undirected(self):
		Graph = nx.Graph()
		Graph.add_nodes_from(self.to_undirected().nodes())
		Graph.add_edges_from(self.to_undirected().edges())
		return Graph

	@property
	def Connected(self):
		Graphs = []
		for Graph in nx.connected_component_subgraphs(self.Undirected):
			Graphs.append(self.SubCopy(Graph.nodes()))
		return Graphs

	def ConnectedSubNodes(self, Nodes):
		self.MaskNodes(Nodes)
		Graphs = []
		for Graph in nx.connected_component_subgraphs(self.Undirected):
			Graphs.append(self.SubCopy(Graph.nodes()))
		self.UnMaskNodes()
		return Graphs

from StdGraphs import *
from collections import OrderedDict, deque, Counter

class SubGraphData(object):
	def __init__(self, Graph, data = None, **attr):
		#super(SubGraphData, self).__init__(data = None, **attr)
		self._BuildAll(Graph)

	def _is_binding_partners(self, Product, Reactant1, Reactant2):
		# Find the binding partner of Reactant that results in Product
		BP = self._OG
		#RemoveNodes = [n for n in BindingPartner.nodes() if n in Reactant1.nodes()]
		BP.MaskNodes(Reactant1.nodes())
		#BindingPartner.remove_nodes_from(RemoveNodes)

		if (len(BP.nodes())+len(Reactant2.nodes())) != len(Product.nodes()):
			BP.UnMaskNodes()
			return 0

		AllConnectedComplenents=sorted(BP.Connected, key = len, reverse=True)

		BPs = 0
		for Component in AllConnectedComplenents:
			if Reactant2.Similar(Component):
				BPs += 1
		BP.UnMaskNodes()
		return BPs
	
	def _node_permutations(self, Graph):
		# All possible node removals
		All = deque()
		for i in range(1, len(Graph.nodes())-1):
			for NodeList in itertools.combinations(Graph.nodes(), i):
				All.append(NodeList)
		return All

	def _canonical(self, Graph1, Graph2):
		# Return True if Graph1 should be replaced with Graph2
		# Graph1 and Graph2 must have the same number of nodes
		N1 = map(lambda x: x.split("_"), Graph1.nodes())
		N2 = map(lambda x: x.split("_"), Graph2.nodes())
		G1 = 0
		G2 = 0
		for i in range(len(N1)):
			if N1[i][0] <= N2[i][0]:
				G1 += 1000
			else:
				G2 += 1000
			if N1[i][1] <= N2[i][1]:
				G1 += 1
			else:
				G2 += 1

	def _add_subgraphs(self, Graph, RemoveNodes, Count):
		# Update the known Count with the induced connected subgraphs of the Graph
		GraphList = Graph.ConnectedSubNodes(RemoveNodes)
		for Graph in GraphList:
			for Key in Count.keys():
				if Graph.Similar(Key):
					if self._canonical(Key, Graph):
						Count.update((Key,))
					else:
						Count[Graph] = Count[Key] + 1
						Count.pop(Key)
					break
			else:
				Count.update((Graph,))

	def _BuildAll(self, Graph):
		# Build and Cache everything
		self._OG = Graph
		self._SubGraphs = self._BuildSubGraphs(self._OG)

		self._Graphs = [self._OG]
		self._Graphs.extend(self._SubGraphs.keys())

		self._GraphIndices = self._BuildReverseLookup(self._Graphs)

		self._NodeLength = self._BuildLengthIndex(self._Graphs)
		self._BP = self._BuildPartnerArray(self._Graphs)

	def _BuildReverseLookup(self, Graphs):
		GraphIndices = dict()
		i = 0
		for G in Graphs:
			GraphIndices[G] = i
			i += 1
		return GraphIndices

	def _BuildPartnerArray(self, Graphs):
		N = len(Graphs)
		# 2-D Array of binding partners (Indcies match _Graphs)
		# First Index is the Product, Second is the Reactant
		BP = np.zeros((N,N), dtype="object")
		for i in range(N):
			SG = self._BuildSubGraphs(Graphs[i]).keys()
			for R1 in SG:
				R2N = len(Graphs[0].nodes())-len(R1.nodes())
				for R2 in self._NodeLength[R2N]:
					BP[i,self._GraphIndices[R2]] += self._is_binding_partners(Graphs[i],R1,R2)
		return BP

	def _BuildSubGraphs(self, Graph):
		# Find all inducable subgraphs of Graph and return a Counter object (fast hash counter)
		SubGraphs = Counter()
		map(lambda x: self._add_subgraphs(Graph, x, SubGraphs), self._node_permutations(Graph))
		return SubGraphs

	def _BuildLengthIndex(self, Graphs):
		# Build an list of all known Graphs by length
		# List will always start with 0 as the index will be the node count
		try:
			SubGraphs = SubGraphs.keys()
		except:
			pass

		Index = dict()
		for Graph in Graphs:
			try:
				Index[len(Graph.nodes())].append(Graph)
			except KeyError:
				Index[len(Graph.nodes())] = [Graph]
		List = [[]]
		for NodeCount in range(1, max(Index.keys())+1):
			try:
				List.append(Index[NodeCount])
			except KeyError:
				List.append([])
		return List

	def is_known(self, Graph):
		# Check if Graph is known to local
		pass

	@property
	def Graphs(self):
		return self._Graphs

	@property
	def MostCommonGraphs(self):
		return self._SubGraphs.most_common()

class TreeLevel(object):
	def __init__(self, Graph, RootData = False, data = None, **attr):
		#super(TreeLevel, self).__init__(data = None, **attr)

		self.AllKnownGraphs = OrderedDict()
		self.OrigGraph = Graph.Copy()
		if RootData:
			self.AllKnownGraphs = copy(RootData)
			for Key in self.AllKnownGraphs.keys():
				self.AllKnownGraphs[Key] = 0
		self.AddAllInducibleSubgraphs()

	@property
	def Children(self):
		ChildKeys = set()
		PartnerKeys = set()

		for Key in self.AllKnownGraphs.keys():
			if self.AllKnownGraphs[Key] > 0:
				if len(Key.nodes()) > 0: # This determines the smallest subunit to be added 0 shows all dimerization reactions
					Node1 = Key
					Node2 = self.ChildBindingPartner(Key)
					if len(PartnerKeys) > 0:
						NoMatch = True
						for Partner in PartnerKeys:
							if GraphCompare(Node1, Partner):
								NoMatch = False
						if NoMatch:
							ChildKeys.add(Node1)
							PartnerKeys.add(Node2)
					else:
						ChildKeys.add(Node1)
						PartnerKeys.add(Node2)
		return ChildKeys

	def CheckKnown(self, Graph):
		# Check if the graph is in the known hash table and if not add it.
		New = True
		if len(self.AllKnownGraphs.keys()) == 0:
			self.AllKnownGraphs[Graph.Copy()] = 1
		for key in self.AllKnownGraphs.keys():
			if Graph.Similar(key):
				self.AllKnownGraphs[key] += 1
				New = False
		if New:
			self.AllKnownGraphs[Graph.Copy()] = 1
		return New

	def AddSubGraphs(self, NodeList):
		#self.OrigGraph.MaskNodes(NodeList)
		#for Graph in self.OrigGraph.Connected:
		for Graph in self.OrigGraph.SubCopy(NodeList).Connected:
			self.CheckKnown(Graph)
		#self.OrigGraph.UnMaskNodes()
	
	def InducedConnectedSubgraphs(self, Graphs, NodeList):
		# Returns a list of the connected components of the induced subgraph
		connected_subgraphs = []
		for Graph in Graph.SubCopy(NodeList).Connected:
			connected_subgraphs.append(Nodes.nodes())
		return connected_subgraphs # As a list of Nodes

	def AddAllInducibleSubgraphs(self):
		for i in range(1, len(self.OrigGraph.nodes())):
			for NodeList in itertools.combinations(self.OrigGraph.nodes(), i):
				self.AddSubGraphs(NodeList)

	def AllNodeSubsets(self, set):
		# At some point think about how to make this faster, (probably a decision tree)
		All = []
		for i in range(1+len(set)):
			for s in itertools.combinations(set, i):
				All.append(list(s))
		return All

	def ChildBindingPartner(self, Child):
		BindingPartner = self.OrigGraph.Copy()
		BindingPartner.remove_nodes_from(n for n in BindingPartner.nodes() if n in Child.nodes())
		AllConnectedComplenents=sorted(BindingPartner.Connected, key = len, reverse=True)

		#print "Connected Comps: {}".format(len(AllConnectedComplenents))
		#for G in AllConnectedComplenents:
		#	print "Node Count: {}".format(len(G.nodes()))
		BindingPartner = AllConnectedComplenents[0]
		return GraphClass(BindingPartner)


class AsmLandscape(nx.DiGraph, TreePlot, object):
	def __init__(self, RootGraph, data = None, **attr):
		super(AsmLandscape, self).__init__(data = None, **attr)

		self.PseudoRoot = dict()
		self.RootData = False
		self.Generate(RootGraph)

	def UniqueNodeName(self, Node1, Node2):
		MaxI = 0
		for Node in self.nodes():
			Node = Node.split(" ")
			#Node = map(lambda x: x.strip(), Node)
			if len(Node) == 3:
				CI = int(Node[2])
				if Node[0] == Node1:
					if MaxI < CI:
						MaxI = CI
				if Node[1] == Node1:
					if MaxI < CI:
						MaxI = CI
				if Node[0] == Node2:
					if MaxI < CI:
						MaxI = CI
				if Node[1] == Node2:
					if MaxI < CI:
						MaxI = CI
		NewNodeName = "{} {} {}".format(Node1, Node2, MaxI+1)
		return NewNodeName, MaxI+1

	def NotDone(self, Level):
		for Node in Level.AllKnownGraphs.keys():
			if 1 > (Node.nodes()):
				return False
		return True

	def AddLevel(self, NLevel, Graph, RootName, Root = False, PreviousLevel = False):
		if not Root:
			Graph2 = PreviousLevel.ChildBindingPartner(Graph)
			if len(Graph.nodes()) > len(Graph2.nodes()):
				Level = TreeLevel(Graph, self.RootData)
				Level2 = TreeLevel(Graph2, self.RootData)
			else:
				Level = TreeLevel(Graph2, self.RootData)
				Level2 = TreeLevel(Graph, self.RootData)
			Children = Level.Children
			Children.update(Level2.Children)
			ChildrenTemp = copy(Children)
			ToRemove = set()
			for Child in Children:
				ChildrenTemp.remove(Child)
				for Test in ChildrenTemp:
					if GraphCompare(Graph2, Child) or GraphCompare(Graph, Child):
						ToRemove.add(Child)
			Children -= ToRemove

		else:
			Level = TreeLevel(Graph)
			Children = Level.Children
		#RootName, i = self.UniqueNodeName(Level.OrigGraph.TxtRep(self.RootNodeCount))
		NewTargetNodes = []
		NewNode = ""
		for Node in Children:
			Protein1Name = Node.TxtRep(self.RootNodeCount)
			Protein2Name = Level.ChildBindingPartner(Node).TxtRep(self.RootNodeCount)
			NewNode, NLevel = self.UniqueNodeName(Protein1Name, Protein2Name)
			if len(Node.nodes()) > 0:
				# print len(Node.nodes())
				NewTargetNodes.append((NewNode, Node.Copy(), NLevel, Level))
			self.add_node(NewNode)
			self.add_edge(RootName, NewNode)
		return NewNode, NewTargetNodes, Level.AllKnownGraphs, Level

	def Generate(self, RootGraph, Show_Timing = True):

		if Show_Timing:
			Time = TimeStamp()
		#print "Root Level Started"

		self.RootNodeCount = len(RootGraph.nodes())
		self.RootGraph = RootGraph.Copy()
		self.RootTxt = RootGraph.TxtRep(self.RootNodeCount)

		NewNodes = self.AddLevel(1, RootGraph, self.RootTxt, Root=True)[1]
		ToDo = deque(NewNodes)
		while len(ToDo) > 0:
			Next = ToDo.pop()
			NewNodes = self.AddLevel(Next[2], Next[1], Next[0], Root=False, PreviousLevel=Next[3])[1]
			if len(NewNodes) > 0:
				ToDo.extend(NewNodes)
			#print "----------------"
			#print "Node:\n{}\nFinished: {} s".format(Next[0], Time.SinceLast)
		print "Total Time: {} s".format(Time.SinceStart)
		

if __name__ == "__main__":
	from time import time
	Starttime = time()
	Graph = StackedRing(3)
	#Graph.Plot(with_labels=True)
	#T = TreeLevel(Graph)
	#print "NodeSubset Count: {}".format(len(T.AllNodeSubsets(Graph.nodes())))
	#print "Completed in {:.2e} seconds".format(time()-Starttime)
	#print "Unique Graphs: {}".format(len(T.AllKnownGraphs.values()))
	#print T.AllKnownGraphs.values()
	#print "Total Counts of Unique Graphs: {}".format(sum(T.AllKnownGraphs.values()))

	#N = 1
	#T.AllKnownGraphs.keys()[N].Plot(with_labels=True)
	#T.ChildBindingPartner(T.AllKnownGraphs.keys()[N]).Plot(with_labels=True)

	#map(lambda x: x.DiPlot(with_labels=True), T.AllKnownGraphs.keys())
	#map(lambda x: T.BindingPartner(x).DiPlot(with_labels=True), T.AllKnownGraphs.keys())

	AL = AsmLandscape(Graph)
	#AL.Plot(with_labels=False)
	AL.Plot(with_labels=True)
	#for Node in AL.nodes():
	#	print Node

	SG = SubGraphData(Graph)
	#map(lambda x: x.DiPlot(with_labels=True), SG.Graphs)
	#map(lambda x: x[0].DiPlot(with_labels=True), SG.MostCommonGraphs)
	#SG.MostCommonGraph.DiPlot(with_labels=True)
	pass