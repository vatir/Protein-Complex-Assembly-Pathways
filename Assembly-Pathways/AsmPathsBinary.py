from collections import OrderedDict, deque
import numpy as np
import itertools
import networkx as nx
from networkx.drawing.nx_pydot import write_dot
from networkx.drawing.nx_pydot import pydot_layout
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pylab as plt

import matplotlib.pylab as plt    

mpl.rcParams['figure.dpi'] = 1600
mpl.rcParams['figure.figsize'] = (
    20,
    20
    )
mpl.rcParams['savefig.transparent'] = True

    
def SpeciesToIndex(S, Species):
    return int(np.where(np.array(Species.values())==S)[0])

def SpeciesToName(S, Species):
    return Species.values()[S]

class Probabilities(object):
    def __init__(self, Data, Concentrations=False):
        self.Data = Data
        self.ProteinsInSpecies = [sum(x == "1" for x in Data.Species.values()[i]) for i in range(Data.SpeciesCount)]
        self._GenReactionsByProduct()
        self.UpdateNorms(Concentrations)

    def LookupReaction(self, Product, R1 = -1, R2 = -1):
        Out = deque()
        for R in self.Data.Reactions:
            if R[3] == Product:
                Out.append(R)
        if R1 == -1 and R2 == -1:
            return Out
        R1 = int(R1)
        R2 = int(R2)
        for R in Out:
            if R[0] == R1 and R[1] == R2:
                return R
            if R[1] == R1 and R[0] == R2:
                return R
        return False

    def _GenReactionsByProduct(self):
        ReactionsByProduct = OrderedDict()
        for Species in self.Data.SpeciesConvert.values():
            for R in self.LookupReaction(Species):
                try:
                    ReactionsByProduct[Species].append(R)
                except:
                    ReactionsByProduct[Species] = list()
                    ReactionsByProduct[Species].append(R)
        self.ReactionsByProduct = ReactionsByProduct

    def UpdateNorms(self, Concentrations = False): 
        # Generate the normalization factors (Return list in the same order as the Data.Species list)
        if hasattr(Concentrations, '__iter__'):
            NormsBySpecies = np.zeros(len(Concentrations))
        else:
            Concentrations = np.ones(self.Data.SpeciesCount)
            NormsBySpecies = np.zeros(self.Data.SpeciesCount)

        for i, Rs in self.ReactionsByProduct.items():
            #print "{} : {}".format(i,Rs)
            for R in Rs:
                # R[2]: Forward
                # R[0]: Reactant 1
                # R[1]: Reactant 2
                # R[3]: Product
                NormsBySpecies[i] += R[2]*Concentrations[R[0]]*Concentrations[R[1]]
                #print "{} : {}".format(R, NormsBySpecies[i])
        self.Concentrations = Concentrations 
        self.NormsBySpecies = NormsBySpecies

    def UpdateConcentrations(self, Fraction, Concentration): 
        #RealConc = np.zeros_like(Fraction)
        RealConc = Fraction*Concentration/self.ProteinsInSpecies
        self.UpdateNorms(RealConc)

    def ReactionFlux(self, Product, R1, R2):
        Forward = self.LookupReaction(Product, R1, R2)
        if hasattr(Forward, '__iter__'):
            Forward = Forward[2]
        else:
            print "Reaction for P: {} R1: {} R2: {} Not Found!!!!".format(Product, R1, R2)
        if self.Concentrations[R1]*self.Concentrations[R2] != 0.0:
            return Forward*self.Concentrations[R1]*self.Concentrations[R2]/(self.NormsBySpecies[Product])
        return 0.0

    def NodeFlux(self, Graph, Node = "Root"):
        if Node == "Root":
            Node = self.Data.Species.values()[-1] # Root Node Name
        Product = nx.get_node_attributes(Graph,"Product")[Node]
        R1 = nx.get_node_attributes(Graph,"Reactant1")[Node]
        R2 = nx.get_node_attributes(Graph,"Reactant2")[Node]
        return self.ReactionFlux(Product,R1,R2)

    def PathIndependentFlux(self, Graph, Probability = True):
        Product = nx.get_node_attributes(Graph,"Product").values()
        R1 = nx.get_node_attributes(Graph,"Reactant1").values()
        R2 = nx.get_node_attributes(Graph,"Reactant2").values()
        Flux = 1.0
        for i in range(len(Product)):
            Flux *= self.ReactionFlux(Product[i],R1[i],R2[i])
        if not Probability:
            # Remove normalization from last step
            Flux *= self.NormsBySpecies[-1]
        return Flux

    def GetChildren(self, Graph, Node):
        return Graph.successors(Node)

class ImportedData(object):
    def __init__(self):
        self.GetSDataFromFiles()
        self.GetRDataFromFiles()
        #         Sort Species
        self.MonomerCount = int(len(self.Species.values()[0])/2)

    def GetSDataFromFiles(self, filename="species3"):
        try:
            self.Species = OrderedDict()
            with open(filename,"r") as SpeciesFile:
                SpeciesList = SpeciesFile.read().strip().split("\n")
                for S in SpeciesList:
                    L = S.split(";;")
                    self.Species[int(L[0])] = L[1]
        
            self.Species = OrderedDict(sorted(self.Species.items(), key=lambda x: x[0]))
            self.SpeciesConvert = OrderedDict()
            for i, S in enumerate(map(int,self.Species.keys())):
                self.SpeciesConvert[int(S)] = i

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
                    if R[0] > R[1]:
                        R0 = R[0]
                        R1 = R[1]
                        R[0] = R1
                        R[1] = R0
                    R[3]=self.SpeciesConvert[R[3]]
                    self.Reactions.append(R)

                # Sort Reactions by Product, R1, R2 (R1 is always less than R2)
                SortList = map(list,list(np.array(self.Reactions)[:,(3,0,1)]))
                L = range(len(SortList))
                for i in L:
                    SortList[i].append(i)
                SortOrder = np.array(sorted(SortList))[:,3]
                self.Reactions = map(list,list(np.array(self.Reactions)[SortOrder]))
                
                return True
        except Exception as e:
            print e
            return False

    def ReactionValue(self, Node, Graph): # Used for cannonical definition
        def AddToNode(N):
            if N >= 0:                
                ##  0  1 2 3  4  5 6
                ## R1 R2 F P B1 B2 R
                Graph.add_node(Node, {
                                            "RIndex":N,
                                            "Forward":self.Reactions[N][2],
                                            "Reverse":self.Reactions[N][6],
                                            "BondType1":self.Reactions[N][4],
                                            "BondType2":self.Reactions[N][5],
                                            "Product":self.Reactions[N][3],
                                            "Reactant1":self.Reactions[N][0],
                                            "Reactant2":self.Reactions[N][1],
                                            })
            else:
                Graph.add_node(Node, RIndex = -1)
        P = SpeciesToIndex(Node.split(" ")[0], self.Species)
        if P == min(self.SpeciesConvert.values()):
            AddToNode(-1)
            return 0
        #if P == min(self.SpeciesConvert.values()) or P == max(self.SpeciesConvert.values()):
        #    AddToNode(0)
        #    return 0
        R1, R2 = Graph.succ[Node].keys()
        
        R1 = SpeciesToIndex(R1.split(" ")[0], self.Species)
        R2 = SpeciesToIndex(R2.split(" ")[0], self.Species)
        for i, R in enumerate(self.Reactions):
            if P == R[3]:
                if R1 in R[:2] and R2 in R[:2]:
                    AddToNode(i)
                    return i+1
        else:
            raise KeyError

def AddReactionDataToGraph(Graph, Data):
    map(lambda i: Data.ReactionValue(i, Graph), Graph.nodes())
    return list()

def LeftRight(N1, N2, Graph, Data):
    R1 = Graph.node[N1]["RIndex"]
    R2 = Graph.node[N2]["RIndex"]
    if R1 > R2:
        return N1, N2
    elif R1 < R2:
        return N2, N1
    elif R1 == R2:
        R1 = CannonicalGraph(Graph, Data, Start = N1).values()
        R2 = CannonicalGraph(Graph, Data, Start = N2).values()
        for i in range(max(len(R1), len(R2))):
            try:
                if R1[i] > R2[i]:
                    return N1, N2
                if R1[i] < R2[i]:
                    return N2, N1
            except Exception as e:
                print e
                if len(R1) > len(R2):
                    return N1, N2
                elif len(R2) > len(R1):
                    return N2, N1
            return N1, N2

def CannonicalGraph(Graph, Data, Start = "1,1,1;1,1,1"):
    AddReactionDataToGraph(Graph, Data)
    CForm = OrderedDict()
    Todo = deque()
    Todo.append(Start)
    while len(Todo) > 0:
        C = Todo.popleft()
        CForm[C] = Graph.node[C]["RIndex"]
        Children = Graph.succ[C].keys()
        if len(Children) != 2:
            pass
        else:
            C1, C2 = Children
            R1 = Graph.node[C1]["RIndex"]
            R2 = Graph.node[C2]["RIndex"]
            I1, I2 = LeftRight(C1, C2, Graph, Data)
            Todo.appendleft(I1)
            Todo.append(I2)
    return CForm

Max = 1
def UniqueNodeName(Graph, NodeName):
    global Max
    Max += 1
    return "{} {}".format(NodeName,Max)

# Import Images
Images = dict()
from glob import glob
for item in glob("./AllPatterns/*.png"):
    try:
        Filename = item.split("/")[2]
        print Filename
        Images[int(Filename[:2])] = Filename
    except:
        print Filename
        raise SystemExit

def ImageLookup(NodeName):
    """
    :string Nodename:
    :return string Filename:
    """
    ID = np.sum(np.reshape(np.array(map(lambda x: np.array(x.split(","), dtype=int), NodeName.split(" ")[0].split(";"))), 6) * [32, 16, 8, 4, 2, 1])

    return Images[ID]

from matplotlib.offsetbox import OffsetImage, AnnotationBbox
def imscatter(x, y, image, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread(image)
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(x, y)
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists


        
def Plot(Graph, with_labels=True, Show=False):
    Fig = plt.figure()
    Axes = Fig.add_subplot(111)
    write_dot(Graph,'test.dot')
    #Fig.set_title("Tree Plot")
    pos=pydot_layout(Graph, prog="dot")
    # nx.draw(Graph,pos,with_labels=with_labels,arrows=False, ax=Axes, hold=True)
    nx.draw(Graph, pos, with_labels=False, arrows=False, ax=Axes, hold=True, nodeshape = "o", width=10.0)

    #ZoomSize = 1.8*10**-0
    ZoomSize = 0.2
    for NodeName, Loc in pos.items():
        imscatter(Loc[0], Loc[1], "{}{}".format('./AllPatterns/', ImageLookup(NodeName)), zoom=ZoomSize)
    # Axes.plot(plot_data)
    for T in Axes.texts:    # Remove unique node labels.
        if T._text.split(" ")[0] == "1,0,0;0,0,0":
            T._text = "1"
        else:
            T._text = "\n"*4+T._text.split(" ")[0].replace(";","\n")
    if Show:
        #plt.show()
                pass
    return Fig

def PathwayFlux(SortedPathways, FluxFunc):
    Flux = deque()
    for Path in SortedPathways.values():
        # By multiplying by the final step normilazation the cancelation will weight the system by the final flux
        Flux.append(FluxFunc(Path))
    return Flux

def PathBarPlot(SortedPathways, FluxFunc, with_labels=True, Show=False, Axes = "None"):
    Index = range(1, len(SortedPathways.values())+1)
    Flux = PathwayFlux(SortedPathways, FluxFunc)

    if Axes == "None":
        Fig = plt.figure()
        Axes = Fig.add_subplot(111)
    Axes.bar(Index,
             Flux,
             align= "center",
             )
    Axes.set_xlim((0.25,len(Index)+1))
    
    Axes.set_xticks(range(1,len(Index)+1,5))

    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    Axes.set_xlabel("Total Flux: {:.4f}".format(np.sum(Flux)))

    YMax = np.sum(Flux)
    Axes.set_ylim((0.0,YMax))
    Breaks = np.unique(np.array(SortedPathways.keys())[:,0],return_counts=True,return_index=True)[2]
    Breaks = np.cumsum(Breaks[::-1])[:-1]
    Flux = np.array(Flux)
    Marks = map(np.sum,np.split(Flux, Breaks))
    Breaks = list(Breaks)
    Breaks.append(Index[-1])
    Breaks = np.array(Breaks)+1.0-0.475
    for i in range(len(Breaks)):
        Axes.axvline(Breaks[i], 0.0,1.0,color="grey")
        Axes.axvline(Breaks[i], 0.0,Marks[i]/YMax,color="red")
    if Show:
        #plt.show()
                pass
    if Axes == "None":
        return Fig
    else:
        return Axes

def Reactants(P, Reactions, Species):
    PR = deque()
    for R in Reactions:
        PT = (SpeciesToIndex(P, Species) == R[3])
        if PT:
            PR.append((R[0], R[1], R[3]))
    return PR

def AsmPathways():
    Asm = nx.DiGraph()
    
    D = ImportedData()

    RootSpeciesNum = int(max(map(int,D.Species.keys())))
    RootSpeciesName = D.Species[RootSpeciesNum]
    RootSpeciesCode = D.SpeciesConvert[RootSpeciesNum]
    
    ##  0  1 2 3  4  5 6
    ## R1 R2 F P B1 B2 R

    MCount = 0
    Tasks = OrderedDict()
    Tree = nx.DiGraph()
    Tree.add_node("1,1,1;1,1,1")
    Tasks[Tree]=deque()
    Tasks[Tree].append("1,1,1;1,1,1")
    while np.sum(map(len,Tasks.values())) > 0:
        for Key in Tasks.keys():
            while len(Tasks[Key]) > 0:
                LastNode = Tasks[Key].pop()
                G = Key

                Rs = Reactants(LastNode.split(" ")[0], D.Reactions,D.Species)

                if Rs > 1:
                    Tree = nx.DiGraph()
                    Tree.add_nodes_from(G.nodes(data=True))
                    for Edge in G.edges():
                        Tree.add_edge(Edge[0], Edge[1])
                for i, R in enumerate(Rs):
                    if i > 0:
                        NewTree = nx.DiGraph()
                        NewTree.add_nodes_from(Tree.nodes(data=True))
                        for Edge in Tree.edges():
                            NewTree.add_edge(Edge[0], Edge[1])
                        G = NewTree
                        Tasks[G] = deque()

                    SN1 = SpeciesToName(R[0], D.Species)
                    NodeName1 = UniqueNodeName(G, SN1)

                    G.add_node(NodeName1)
                    G.add_edge(LastNode, NodeName1)

                    SN2 = SpeciesToName(R[1], D.Species)
                    NodeName2 = UniqueNodeName(G, SN2)

                    G.add_node(NodeName2)

                    G.add_edge(LastNode, NodeName2)

                    if SN1 != "1,0,0;0,0,0":
                        Tasks[G].append(NodeName1)
                    else:
                        MCount += 1
                    if SN2 != "1,0,0;0,0,0":
                        Tasks[G].append(NodeName2)
                    else:
                        MCount += 1

                    Degree = G.degree()
                    for K in Degree.keys():
                        if K.split(" ")[0] == "1,1,1;1,1,1":
                            pass
                        elif K.split(" ")[0] == "1,0,0;0,0,0":
                            pass
                        elif Degree[K] != 3 and K not in Tasks[G]:
                            Tasks[G].append(K)

    # Unique Graphs
    Graphs = {}
    InvalidGraphs= {}

    def InList(L1, L2, Index = False):
        N1 = np.array(L1)
        for k, i in enumerate(L2):
            Ni = np.array(i)
            if (N1 == Ni).all():
                if Index:
                    return True, k
                else:
                    return True
        return False

    for G in Tasks:
        GC = CannonicalGraph(G, D).values()
        GC = np.array(GC)
        GC = tuple(GC)
        Graphs[GC]=G

    #print "Unique Graphs Found: {}".format(len(Graphs.keys()))
    
    return D, OrderedDict([x for x in sorted(Graphs.items(),reverse=True)])

    #for k in Graphs.keys():
    #    print InList(k, SortedPathways.keys())

if __name__ == "__main__":
    D, SortedPathways = AsmPathways()
    OutputGraphs = True
    # print mpl.rcParams
    if OutputGraphs:
        for i, Key in enumerate(SortedPathways.values()):
            F = Plot(
                Key,
                Show=False
                )
            F.savefig("./Images/Pathway{}.pdf".format(i+1), 
                                bbox_inches = 'tight',
                                transparent=True,
                                #dpi = 2400
                                #dpi = 'figure'
                                )
            plt.close(F)
            del F

    # Everything in this block is for later
    Later = True
    if not Later:
        # Probability Calculations
    
        #Concentrations = np.arange(-11,1)*-1.0*10.0**-7.0/float(len(D.Species))
        Concentrations = np.ones(len(D.Species))/float(len(D.Species))
        #Concentrations = np.zeros(len(D.Species)); Concentrations[0]
    
        P = Probabilities(D,Concentrations)
        #print P.NormsBySpecies
        #P.UpdateNorms(Concentrations)
        #print P.NormsBySpecies
        # List of Pathways: SortedPathways.values().nodes(data=True)
        try:
            P.UpdateNorms(FullOutArray[KD1==KD1Range,KD2==KD2Range, A0Target, TimeTarget,:][0])
        except:
            pass
    
        Total = 0.0
        for Path in SortedPathways.values():
            P.NodeFlux(Path)
            Total += P.NodeFlux(Path)
        #print "Total: {}".format(Total)
        Total = 0.0
        for Path in SortedPathways.values():
            P.PathIndependentFlux(Path)
            Total += P.PathIndependentFlux(Path)
        #print "Total: {}".format(Total)
        #PathBarPlot(SortedPathways, P.NodeFlux, Show=True)
    
        Concentrations = np.arange(-11,1)*-1.0*10.0**-7.0/float(len(D.Species))
        
        # Test Array Flux
        from time import time
        StartTime = time()
        CurrentTime = StartTime
        CurrentConcentration = np.load("ConcTestingData.npy")
        
    
        TimeArray = np.logspace(-3,6,num=CurrentConcentration.shape[0])
        TimeSpace = np.diff(TimeArray)
    
        FinalFlux = np.zeros(len(SortedPathways.values()))
        WFluxSum = np.zeros(len(SortedPathways.values()))
        Flux3D = np.zeros([len(TimeArray)-1,len(SortedPathways.values())])
        #print "--------------------"
        for i in range(CurrentConcentration.shape[0]-1):
            P.UpdateConcentrations(CurrentConcentration[i,:], 1.0)
            Total = 0.0
            WTotal = 0.0
            for k, Path in enumerate(SortedPathways.values()):
                T = P.PathIndependentFlux(Path)
                WT = T*TimeSpace[i]
                FinalFlux[k] = T
                WFluxSum[k] += WT
                Total += T
                WTotal += WT
            Flux3D[i] = FinalFlux[:]
            #print "Index: {}".format(i)
            #print Total
            #print WTotal
            #print "--------------------"
        #print Flux3D
        #WFluxTotal = np.sum(WFluxSum)/float(TimeArray[-1]-TimeArray[0])
        #print FluxTotal
        #print WFluxTotal
        CurrentTime = time()
        print "Pathways Calc Time: {:.5e}".format(CurrentTime-StartTime)
    
        
        #WFlux = WFluxSum/float(TimeArray[-1]-TimeArray[0])
        #PathBarPlot(SortedPathways, P.PathIndependentFlux, Show=True)
        """
        import matplotlib.pylab as plt
        Fig = plt.figure()
        Axes = Fig.add_subplot(111)
    
        Index = range(1,len(FinalFlux)+1)
    
        WFlux = WFluxSum/float(TimeArray[-1]-TimeArray[0])
        Axes.bar(Index, WFlux, color="blue", align= "center", alpha = 0.5)
    
        FinalFlux = np.array(FinalFlux)
        Axes.bar(Index, FinalFlux, color="red",align= "center", alpha = 0.5)
        Axes.set_xlim((0.5,len(Index)+1))
        Axes.set_ylim((0.0,1.0))
        Axes.set_xticks(range(1,len(Index)+1,5))
    
        left, width = .25, .5
        bottom, height = .25, .5
        right = left + width
        top = bottom + height
        Axes.set_xlabel("Total Weighted Flux: {:.4f}".format(np.sum(WFlux)))
        Breaks = np.unique(np.array(SortedPathways.keys())[:,0],return_counts=True,return_index=True)[2]
        Breaks = np.cumsum(Breaks[::-1])[:-1]
        Marks = map(np.sum,np.split(FinalFlux, Breaks))
        WMarks = map(np.sum,np.split(WFlux , Breaks))
        Breaks = list(Breaks)
        Breaks.append(Index[-1])
        Breaks = np.array(Breaks)+1.0-0.5
        for i in range(len(Breaks)):
            Axes.axvline(Breaks[i], 0.0,1.0,color="grey",alpha=0.5)
            Axes.axvline(Breaks[i], 0.0,WMarks[i],color="blue")
            Axes.axvline(Breaks[i], 0.0,Marks[i],color="red")
        plt.show()
        """
        
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        ExtFig = plt.figure(facecolor="w")
        
        ExtAxes = ExtFig.add_subplot(121, projection='3d')
        Index = range(1,len(FinalFlux)+1)
        for i, Time in enumerate(TimeArray[:-1]):
            color = ["b"]*len(Index)
            color[np.argmax(Flux3D[i,:])] = "r"
            ExtAxes.bar(Index, Flux3D[i,:], zs=np.log10(Time), zdir="y", color=color, alpha = 0.5)
    
        ExtAxes = ExtFig.add_subplot(122, projection='3d')
        Index = range(CurrentConcentration.shape[1])
        for i, Time in enumerate(TimeArray[:-1]):
            color = ["b"]*len(Index)
            color[np.argmax(CurrentConcentration[i,:])] = "r"
            ExtAxes.bar(Index, CurrentConcentration[i,:], zs=np.log10(Time), color=color, zdir="y", alpha = 0.5)
        plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
        #plt.show(ExtFig)
