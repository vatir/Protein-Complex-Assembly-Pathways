import numpy as np
import matplotlib.pylab as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso

def GraphCompare(Graph1, Graph2):
	if len(Graph1.nodes()) != len(Graph2.nodes()):
		return False
	if len(Graph1.edges()) != len(Graph2.edges()):
		return False
	em = iso.categorical_edge_match("Type", "")
	if not nx.is_isomorphic(Graph1, Graph2, edge_match=em):
		return False
	return True

def BindingPartner(Reactant, Product):
	pass

class GraphPlot(object):
	"""
	Holder for extra fuctions for all Graphs.
	"""
	def __init__(self, data = None, **attr):
		return super(GraphPlot, self).__init__(data, **attr)
	def Plot(self, Show = True, **attr):
		return self.PlotSpring(Show, **attr)
	def PlotSpring(self, Show = True, **attr):
		Fig = plt.figure()
		Axes = Fig.add_subplot(111)
		Axes.plot(nx.draw_spring(self.to_undirected(), **attr))
		if Show:
			plt.show(Fig)
		return Fig

	def PlotCircular(self, Show = True, **attr):
		Fig = plt.figure()
		Axes = Fig.add_subplot(111)
		Axes.plot(nx.draw_circular(self.to_undirected(), **attr))
		if Show:
			plt.show(Fig)
		return Fig

	def DiPlot(self, Show = True, **attr):
		Fig = plt.figure()
		Axes = Fig.add_subplot(111)
		Axes.plot(nx.draw_spring(self, **attr))
		if Show:
			plt.show(Fig)
		return Fig

from collections import OrderedDict
class TreePlot(object):
	"""
	Holder for extra fuctions for all Tree Style Graphs.
	"""
	def __init__(self, data = None, **attr):
		#return super(TreePlot, self).__init__(data, **attr)
		pass
	def _NodeLengths(self, N):
		try:
			N1, N2 = N.split(" ")[:2]
			SizeN1 = np.sum(np.array([i != "0" for i in N1]))
			SizeN2 = np.sum(np.array([i != "0" for i in N2]))
			print "---------------------------"
			print "{} : {}".format(N1, SizeN1)
			print "{} : {}".format(N2, SizeN2)
			return SizeN1, SizeN2
		except:
			SizeN1 = np.sum(np.array([i != "0" for i in N]))
			return SizeN1, 0


	def Plot(self, Show = True, with_labels = True, **attr):
		Fig = plt.figure(figsize=(12, 4), dpi=70)
		Axes = Fig.add_subplot(111)
		#print self.nodes()
		nx.write_dot(self,'test.dot')
		#Fig.set_title("Tree Plot")
		pos=OrderedDict(nx.graphviz_layout(self,prog='dot'))
		
		# NodeColor
		P = "#AD18AB" # Purple
		B = "#003CFF" # Blue
		O = "#F55714" # Orange
		R = "#FF0000" # Red

		NodeColor = []
		for p in self.nodes():
			SN1, SN2 = self._NodeLengths(p)
			if SN1 == 6:
				NodeColor.append(P)
			elif (SN1 == 1) and (SN2 == 1):
				NodeColor.append(R)
			elif (SN1 == 1) or (SN2 == 1):
				NodeColor.append(B)
			else:
				NodeColor.append(O)



		#Axes.plot(nx.draw(self,pos,with_labels=with_labels,arrows=False))
		nx.draw(self,pos,arrows=False, nodes = pos.keys(), ax=Axes, width=3.0,node_shape="o", node_size =350, node_color=NodeColor, with_labels=with_labels)
		#Axes.plot(TreeData,linewidth=10.0)
		sep = " "
		for i, T in enumerate(Axes.texts):	# Remove unique node labels.
			#if T._text.split(" ")[:-1] == "1,0,0;0,0,0":
			#	T._text = "1"
			#else:
			#T._text = "\n"*4+T._text.split(" ")[:-1].replace(" ","\n")

			#SN1, SN2 = self._NodeLengths(T._text)
			#if SN1 == 6:
			#	Axes.axes.collections[0]._facecolors[i] = np.array([ 0.67843137,  0.09411765,  0.67058824,  1.        ])
			#	Axes.axes.collections[0]._facecolors_original[i] = np.array([ 0.67843137,  0.09411765,  0.67058824,  1.        ])
			#elif (SN1 == 1) and (SN2 == 1):
			#	Axes.axes.collections[0]._facecolors[i] = np.array([ 0.        ,  0.0,  0.0        ,  1.        ])
			#	Axes.axes.collections[0]._facecolors_original[i] = np.array([ 0.        ,  0.0,  0.0        ,  1.        ])
			#elif (SN1 == 1) or (SN2 == 1):
			#	Axes.axes.collections[0]._facecolors[i] = np.array([ 0.96078431,  0.34117647,  0.07843137,  1.        ])
			#	Axes.axes.collections[0]._facecolors_original[i] = np.array([ 0.96078431,  0.34117647,  0.07843137,  1.        ])
			#else:
			#	Axes.axes.collections[0]._facecolors[i] = np.array([ 0.96078431,  0.34117647,  0.07843137,  1.        ])
			#	Axes.axes.collections[0]._facecolors_original[i] = np.array([ 0.96078431,  0.34117647,  0.07843137,  1.        ])

			if with_labels and self.degree()[T._text] > -1 and True:
				T._text = "\n"*4+sep.join(T._text.split(" ")[:-1]).replace(" ","\n")
			else:
				T._text = ""

		if Show:
			plt.show()
		return Fig

from time import time
class TimeStamp(object):
	def __init__(self, *args, **kwargs):
		self._StartTime = time()
		self._CurrentTime = time()
		return super(TimeStamp, self).__init__(*args, **kwargs)

	@property
	def SinceStart(self):
		self._CurrentTime = time()
		return "{:.2e}".format(self._CurrentTime-self._StartTime)

	@property
	def SinceLast(self):
		NewTime = time()
		OldTime = self._CurrentTime
		self._CurrentTime = NewTime
		return "{:.2e}".format(NewTime-OldTime)
