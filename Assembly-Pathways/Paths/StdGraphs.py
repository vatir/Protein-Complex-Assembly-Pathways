from Helpers import GraphPlot
from Graphs import GraphClass
import networkx as nx
from copy import copy

class Star(GraphClass, GraphPlot, object):
	"""
	Star Case

	Don't use for now needs edge symetry information clarified.
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
		EdgeAttr = dict()
		#EdgeAttr["Type"] = "{}Ring".format(Prefix1)
		EdgeAttr["Type"] = "Intra-Ring"
		EdgeAttr["Directed"] = True
		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix1,(x+1)%RingLength+1),EdgeAttr), range(RingLength))
		# Add second ring edges
		EdgeAttr = dict()
		#EdgeAttr["Type"] = "{}Ring".format(Prefix2)
		EdgeAttr["Type"] = "Intra-Ring"
		EdgeAttr["Directed"] = True
		map(lambda x: self.add_edge("{}{}".format(Prefix2,x+1),"{}{}".format(Prefix2,(x+1)%RingLength+1),EdgeAttr), range(RingLength))

		# Add first set of cross edges
		EdgeAttr = dict()
		EdgeAttr["Type"] = "Inter-Ring Type 1"
		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix2,(x+1)%RingLength+1),EdgeAttr), range(RingLength))
		# Add second set of cross edges
		EdgeAttr = dict()
		EdgeAttr["Type"] = "Inter-Ring Type 2"
		map(lambda x: self.add_edge("{}{}".format(Prefix2,x+1),"{}{}".format(Prefix1,(x+1)%RingLength+1),EdgeAttr), range(RingLength))

class Ring(GraphClass,GraphPlot, object):
	"""
	Ring Case
	"""
	def __init__(self, N = 6, Prefix = "A_", data = None, **attr):
		super(Ring, self).__init__(data, **attr)
		# Add nodes
		map(lambda x: self.add_node("{}{}".format(Prefix,x+1)), range(N))
		# Add edges
		EdgeAttr = dict()
		EdgeAttr["Type"] = "{}Ring".format(Prefix)
		EdgeAttr["Directed"] = True 
		# A directed edge is only equivalent to another if the order of the nodes is the same. 
		# (In these cases the number in the node name goes from lower to higher)
		map(lambda x: self.add_edge("{}{}".format(Prefix,x+1),"{}{}".format(Prefix,(x+1)%N+1), EdgeAttr), range(N))

#class StackedRing(GraphClass,GraphPlot, object):
#	"""
#	Stacked Ring Case
#	"""
#	def __init__(self, N = 6, Prefix1 = "A_", Prefix2 = "B_", data = None, **attr):
#		super(StackedRing, self).__init__(data, **attr)
#		if int(N) == N:
#			RingLength = int(N/2.0)
#		else:
#			raise TypeError
#		# Add first ring nodes
#		map(lambda x: self.add_node("{}{}".format(Prefix1,x+1)), range(RingLength))
#		# Add second ring nodes
#		map(lambda x: self.add_node("{}{}".format(Prefix2,x+1)), range(RingLength))
#		# Add first ring edges
#		EdgeAttr = dict()
#		#EdgeAttr["Type"] = "{}Ring".format(Prefix1)
#		EdgeAttr["Type"] = "Intra-Ring"
#		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix1,(x+1)%RingLength+1), EdgeAttr), range(RingLength))
#		# Add second ring edges
#		EdgeAttr = dict()
#		#EdgeAttr["Type"] = "{}Ring".format(Prefix2)
#		EdgeAttr["Type"] = "Intra-Ring"
#		map(lambda x: self.add_edge("{}{}".format(Prefix2,x+1),"{}{}".format(Prefix2,(x+1)%RingLength+1), EdgeAttr), range(RingLength))
#		# Add links inter-ring edges
#		EdgeAttr = dict()
#		EdgeAttr["Type"] = "Inter-Ring"
#		map(lambda x: self.add_edge("{}{}".format(Prefix1,x+1),"{}{}".format(Prefix2,x+1), EdgeAttr), range(RingLength))

class NStackedRing(GraphClass, GraphPlot, object):
	"""
	K rings stacked: Prefix should be a list with the names of the rings. N is the number of proteins per ring.
	The defaults generate a proteasome like structure.
	"""
	def __init__(self, N = 7, Prefix = ["B1_", "A1_", "A2_", "B2_"], data = None, **attr):
		super(NStackedRing, self).__init__(nx.DiGraph(), data, **attr)
		self.Initialize(N, Prefix)

	def Initialize(self, N, Prefix):
		# Generate Rings and their Intra-Ring Edges
		#for Name in Prefix:
		#	map(lambda x: self.add_node("{}{}".format(Name,x+1)), range(N))
		#	EdgeAttr = dict()
		#	EdgeAttr["Type"] = "Intra-Ring"
		#	# A directed edge is only equivalent to another if the order of the nodes is the same. 
		#	# (In these cases the number in the node name goes from lower to higher)
		#	map(lambda x: self.add_edge("{}{}".format(Name,x+1),"{}{}".format(Name,(x+1)%N+1), EdgeAttr), range(N))

		# Add links inter-ring edges
		EdgeAttr = dict()
		EdgeAttr["Type"] = "Inter-Ring"
		for i in range(len(Prefix)-1):
			map(lambda x: self.add_edge("{}{}".format(Prefix[i],x+1),"{}{}".format(Prefix[i+1],x+1), EdgeAttr), range(N))
			map(lambda x: self.add_edge("{}{}".format(Prefix[i+1],x+1),"{}{}".format(Prefix[i],x+1), EdgeAttr), range(N))

		# Create DiGraph Version for asymetric edges
		Forward = True # Alternate ring bond directions
		EdgeAttr = dict()
		EdgeAttr["Type"] = "Intra-Ring"
		for Name in Prefix:
			if Forward:
				map(lambda x: self.add_edge("{}{}".format(Name,x+1),"{}{}".format(Name,(x+1)%N+1), EdgeAttr), range(N))
				Forward = False
			else:
				map(lambda x: self.add_edge("{}{}".format(Name,(x+1)%N+1),"{}{}".format(Name,x+1), EdgeAttr), range(N))
				Forward = True

		self._OrigNodes = copy(self.nodes(data=True))
		self._OrigEdges = copy(self.edges(data=True))


class StackedRing(NStackedRing, object):
	"""
	Stacked Ring Case (Special case of the NStackedRing)
	"""
	def __init__(self, N = 3, Prefix1 = "A_", Prefix2 = "B_", data = None, **attr):
		super(StackedRing, self).__init__(N, [Prefix1, Prefix2], data, **attr)

if __name__ == "__main__":
	NStackedRing().DiPlot(with_labels=True)
	StackedRing().DiPlot(with_labels=True)
	#Ring().Plot(with_labels=True)
	#Star().Plot(with_labels=True)
	#print NStackedRing(7).nodes()
	#for Edge in NStackedRing(7).edges(data=True):
	#	print Edge
