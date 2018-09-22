import networkx as nx
from Helpers import GraphPlot
import itertools as it

class Star(nx.Graph,GraphPlot, object):
	"""
	Star Case
	"""
	def __init__(self, N = 6, Prefix1 = "A_", Prefix2 = "B_", data = None, **attr):
		super(Star, self).__init__(data, **attr)
		if int(N) == N:
			RingLength = int(N/2.0)

		# Add first ring nodes
		map(lambda x: self.add_node("{}{}".format(Prefix1, x+1)), range(RingLength))
		# Add second ring nodes
		map(lambda x: self.add_node("{}{}".format(Prefix2, x+1)), range(RingLength))

		# Add first ring edges
		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix1,(x+1)%RingLength+1)), range(RingLength))
		# Add second ring edges
		map(lambda x: self.add_edge("{}{}".format(Prefix2,x+1),"{}{}".format(Prefix2,(x+1)%RingLength+1)), range(RingLength))

		# Add first set of cross edges
		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix2,(x+1)%RingLength+1)), range(RingLength))
		# Add second set of cross edges
		map(lambda x: self.add_edge("{}{}".format(Prefix2,x+1),"{}{}".format(Prefix1,(x+1)%RingLength+1)), range(RingLength))

class Ring(nx.Graph,GraphPlot, object):
	"""
	Ring Case
	"""
	def __init__(self, N = 6, Prefix = "A_", data = None, **attr):
		super(Ring, self).__init__(data, **attr)
		# Add nodes
		map(lambda x: self.add_node("{}{}".format(Prefix,x+1)), range(N))
		# Add edges
		map(lambda x: self.add_edge("{}{}".format(Prefix,x+1),"{}{}".format(Prefix,(x+1)%N+1)), range(N))

class StackedTrimer(nx.Graph,GraphPlot, object):
	"""
	Stacked Trimer Case
	"""
	def __init__(self, N = 6, Prefix1 = "A_", Prefix2 = "B_", data = None, **attr):
		super(StackedTrimer, self).__init__(data, **attr)
		if int(N) == N:
			RingLength = int(N/2.0)
		else:
			raise TypeError
		# Add first ring nodes
		map(lambda x: self.add_node("{}{}".format(Prefix1,x+1)), range(RingLength))
		# Add second ring nodes
		map(lambda x: self.add_node("{}{}".format(Prefix2,x+1)), range(RingLength))
		# Add first ring edges
		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix1,(x+1)%RingLength+1)), range(RingLength))
		# Add second ring edges
		map(lambda x: self.add_edge("{}{}".format(Prefix2,x+1),"{}{}".format(Prefix2,(x+1)%RingLength+1)), range(RingLength))
		# Add links ring edges
		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix2,x+1)), range(RingLength))

class GraphClass(nx.Graph, GraphPlot, object):
	def __init__(self, Graph, data = None, **attr):
		super(GraphClass, self).__init__(data, **attr)
		map(lambda x: self.add_node(x), Graph.nodes())
		map(lambda x: self.add_edge(x[0], x[1]), Graph.edges())


class PathwayTree(nx.DiGraph, object):
	AllKnownGraphs = dict()
	def __init__(self, Graph, data = None, **attr):
		super(PathwayTree, self).__init__(data = None, **attr)
		

	def CheckKnown(self, Graph):
		# Check if the graph is in the known hash table and if not add it.
		New = True
		if len(self.AllKnownGraphs.keys()) == 0:
			self.AllKnownGraphs[GraphClass(Graph)] = 1
		for key in self.AllKnownGraphs.keys():
			if nx.is_isomorphic(key, Graph):
				self.AllKnownGraphs[key] += 1
				New = False
		if New:
			self.AllKnownGraphs[GraphClass(Graph)] = 1
		return New

	def CheckList(self, GraphList):
		for Graph in GraphList:
			self.CheckKnown(Graph)
	
	def induced_connected_subgraphs(self, Graph, NodeList):
		# Returns a list of the connected components of the induced subgraph
		graph = nx.copy.deepcopy(Graph)

		#Brute force method as new-style class super function seems to be causing an error
		#graph = nx.Graph()
		#map(lambda x: graph.add_node(x), Graph.nodes())
		#map(lambda x: graph.add_edge(x[0], x[1]), Graph.edges())
		graph.remove_nodes_from(NodeList)
		connected_subgraphs = nx.connected_component_subgraphs(graph)
		"""
		subgraphs = []
		for subgraph in subgraph_edges:
			output_graph = nx.Graph()
			edgelist = []
			for edge in edge_set:
				edgelist.append(edge)
				output_graph.add_edge(edge[0], edge[1])
			print edgelist
			subgraphs.append(GraphClass(output_graph))
		"""
		return connected_subgraphs

	def AddAllInducibleSubgraphs(self, Graph):
		for node_set in self.AllNodeSubsets(Graph.nodes()):
			self.CheckList(self.induced_connected_subgraphs(Graph, node_set))

	def AllNodeSubsets(self, set):
		All = []
		for i in range(len(set)+1):
			for s in it.combinations(set, i):
				All.append(list(s))
		return All
	

if __name__ == "__main__":
	Graph = StackedTrimer(6)
	#Graph.Plot(with_labels=True)
	T = PathwayTree("")
	print "NodeSubset Count: {}".format(len(T.AllNodeSubsets(Graph.nodes())))
	T.AddAllInducibleSubgraphs(Graph)
	print "Unique Graphs: {}".format(len(T.AllKnownGraphs.values()))
	print T.AllKnownGraphs.values()
	print "Total Counts of Unique Graphs: {}".format(sum(T.AllKnownGraphs.values()))
	map(lambda x: x.Plot(with_labels=True), T.AllKnownGraphs.keys())
	#for G in T.AllKnownGraphs.keys():
	#	print nx.is_connected(G)

