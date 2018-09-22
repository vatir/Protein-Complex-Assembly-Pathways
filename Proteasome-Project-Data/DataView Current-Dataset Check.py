#!/usr/bin/env python
# --------------------------------------------
# Options

# GUI

try:
    get_ipython().config
    iPython = True
except:
    iPython = False
WebInterface = True
ShowSecondary3DScreen = True
ShowSecondaryPathwaysConcScreen = True
ShowPaths = True

ShowRingCurve = False
# --------------------------------------------

# Add Binary Asm Pathways
import AsmPathsBinary as AB
D, SortedPathways = AB.AsmPathways()
P = AB.Probabilities(D) 
def PathwayPlot(Fraction, Concentration, Axes):
    P.UpdateConcentrations(Fraction, Concentration)
    return AB.PathBarPlot(SortedPathways, P.PathIndependentFlux, Show=False, Axes=Axes)
# -----------------------------------------------------------------------------------------
def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx]

def ValidbyTotalConc(TotalConc, Ep = 2e-1):
    Valid = (TotalConc > (1.0 - Ep)) & (TotalConc < (1.0 + Ep))
    return Valid

from scipy.interpolate import InterpolatedUnivariateSpline
def InflectPoints(TP, Frac, TotalConc, DerEp = 1e-3, D1Cut = 1e-1):
    try:
        Valid = ValidbyTotalConc(TotalConc)
        TP = TP[Valid]
        Frac = Frac[Valid]
        LogTP = np.log10(TP)
        s5 = InterpolatedUnivariateSpline(LogTP, Frac, k=5)
        #s5.set_smoothing_factor(1e-3)
        D2 = s5.derivative(2)
        R2 = D2.roots()
    
        s4 = InterpolatedUnivariateSpline(LogTP, Frac, k=4)
        #s4.set_smoothing_factor(1e-3)
        D1 = s4.derivative(1)
        #R1 = D1.roots()
        
        ITest = R2[D1(R2) >= D1Cut]
    
        DevCheck1 = np.sign(D1(ITest-DerEp))==np.sign(D1(ITest+DerEp))
        DevCheck2 = np.sign(D2(ITest-DerEp))!=np.sign(D2(ITest+DerEp))
        DevCheck = np.logical_and(DevCheck1,DevCheck2)
        #DevCheck = DevCheck1
        #DevCheck = DevCheck2
        # print "Radius of Curvature: {}".format(np.abs((1+D1(R2)**1.5)/D2(R2)))
    
        return ITest[DevCheck]
    except:
        print "Inflection Point Error!!!!"
        return []

# -----------------------------------------------------------------------------------------
# Cleanup for use in iPython

# -----------------------------------------------------------------------------------------
# General Imports
import numpy as np
from datetime import datetime
from time import time

# -----------------------------------------------------------------------------------------
# Load Dataset
import h5py
import importlib

N = 1
#WorkingName = "STV2StaticT"
WorkingName = "STV3StaticT"
#WorkingName = "TriStackedRingODELargeSynDeg28e-04TAfterFix"

# New testing datasets
#N = 5
#WorkingName = "TriStackedRingODELargeSynDeg17e-03T"

#N = 1
#WorkingName = "TriStackedRingODELargeSynDeg28e-05T"
#WorkingName = "TriStackedRingODELargeSynDeg15e-06T"
#WorkingName = "TriStackedRingODELargeSynDeg90e-06T"

#WorkingName = "TriStackedRingODELargeStaticT"

#WorkingName = "TriStackedRingODELargeSynDeg28e-04T"
#WorkingName = "TriStackedRingODELargeSynDeg17e-03T"
# Small Testing Dataset
#WorkingName = "TriStackedRingODE1GenPerHour"

def ImportData():
    global N
    global WorkingName
    global FullOutArray
    global Errors
    global HasData
    global Runtimes
    global A0Range
    global KD1Range
    global KD2Range
    global Kp
    global TimePoints
    global SpeciesCount
    global Delta
    global Deg

    try:
        del FullOutArray
    except:
        pass

    File = h5py.File("./{}{}.h5".format(N, WorkingName), "r")
    if N == 1:
        CD = importlib.import_module("{}{}".format("",WorkingName))
    else:
        CD = importlib.import_module("{}{}".format(N-1,WorkingName))

    FullOutArray = File["Results"][:]
    FullOutArray = FullOutArray[::-1,::-1,:,:,:]
    
    
    Errors = File["Errors"][:]
    HasData = File["HasData"][:]
    Runtimes = File["Runtimes"][:]
    A0Range = File["A0Range"][:]
    Delta = CD.Delta
    Deg = Delta

    KD1Range = File["KD1Range"][:]
    KD1Range = KD1Range[::-1]
    
    KD2Range = File["KD2Range"][:]
    KD2Range = KD2Range[::-1]

    Kp = File["Kp"][()]
    TimePoints = CD.KeepTimePoints
    SpeciesCount = CD.SpeciesCount
    print "Working Name: {}".format(WorkingName)
    print "Shape: {}".format(Errors.shape)
    print "Errors: {}".format(np.sum(Errors==True))
    print "HasData: {}".format(np.sum(HasData==True))

ImportData()
# End: Load Dataset
# -----------------------------------------------------------------------------------------
# Start Plotting


# Choose Matplotlib Backend

if WebInterface and not iPython:
    import matplotlib
    matplotlib.use("webagg") # Working
    import sys
    try:
        port = sys.argv[1] 
    except:
        port = 9000
    matplotlib.rcParams['webagg.port'] = port
    print "Webserver on Port: {}".format(port)
#    matplotlib.use('module://mplh5canvas.backend_h5canvas')
    import matplotlib.pyplot as plt
elif iPython:
    import matplotlib
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt
    plt.ioff()
else:
    import matplotlib
    import matplotlib.pyplot as plt
    
import mplh5canvas
from matplotlib.widgets import Button
from matplotlib.widgets import Slider
import matplotlib.cm as cm # Color Map
    
#path = 'C:\\Anaconda2\\Lib\\site-packages\\matplotlib\\mpl-data\\fonts\\ttf\\pala.ttf' # Normal
#path = 'C:\\Anaconda2\\Lib\\site-packages\\matplotlib\\mpl-data\\fonts\\ttf\\palab.ttf' # Bold
path = './Fonts/pala.ttf' # Normal
path = './Fonts/palab.ttf' # Bold

prop = matplotlib.font_manager.FontProperties(fname=path)
prop = matplotlib.font_manager.FontProperties(fname=path)
matplotlib.rcParams['font.family'] = prop.get_name()
matplotlib.rcParams['font.weight'] = 'normal'
matplotlib.rcParams["font.size"] = 15


#matplotlib.rcParams["figure.dpi"] = 600
#matplotlib.rcParams["figure.figsize"] = list(np.array([8.0,6.0])*(72.0/600.0))
"""
MainFigure.savefig("E:\\test.png", format="png", dpi=1800)


"""


#matplotlib.rc('font', family='Palatino') 
#matplotlib.rc('text', usetex='false') 
#matplotlib.rcParams.update({'font.size': 22})



# Setup Plot Params
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

# Setup Plot Substucture
MainFigure = plt.figure("Main", figsize=(10*1.8,10))
MainFigure.patch.set_facecolor('w')

FvsTime = plt.subplot(221)
FvsConc = plt.subplot(223)
ContourPlt = plt.subplot(122, aspect='equal')


# Setup Initial Values
Species = 11

# A0 = 3.5e-8
# Time = 10800.0
# KD1 = 1e-7
# KD2 = 1e-7

A0 = 4.0e-6
# A0 = 1.0e-7
Time = 1.0e8
KD1 = 1e-12
KD2 = 1e-12

# Setup Initial Derived Values
Tau = Kp * A0 * Time
TimeTarget = np.argmax(TimePoints==find_nearest(TimePoints, Time))
A0Target = np.argmax(A0Range==find_nearest(A0Range, A0))
#TimeTarget = np.argmax(Time==TimePoints)
#A0Target = np.argmax(A0==A0Range)

# -----------------------------------------------------------------------------------------

# Testing CRN Data addition for flux weighting
# Initial Setup
from ODERun import ODEFunc
ODE = ODEFunc(SpeciesCount)
ODE.GetSDataFromFiles()
ODE.GetRDataFromFiles()
for i, entry in enumerate(ODE.Species.items()):
    print "{:2} : {:2} : {}".format(i, entry[0], entry[1])
# Update
#Fraction = FullOutArray[KD1Range==KD1,KD2Range==KD2,A0Target,TimeTarget,:][0]
#Fraction = A0Range[A0Target]*Fraction/float(P.ProteinsInSpecies[-1])
#ODE.SetConc(Fraction)
#ODE.GenCoef(
#                Kd1   = KD1,
#                Kd2   = KD2,
#                Kp    = Kp,
#                Delta = Delta
#                )
#FinalFlux = ODE.Der(ODE.Conc_0, None)[-1]

# -----------------------------------------------------------------------------------------
if ShowSecondaryPathwaysConcScreen:
    from mpl_toolkits.mplot3d import Axes3D
    def ExtConcPlots(
                        Frac,
                        ExtFig = None,
                        ExtAxes1 = None,
                        ExtAxes2 = None,
                        ExtFig3 = None,
                        Ext3Axes1 = None,
                        Redraw = True,
                     ):
        if len(A0Range) > 5:
            A0Start = np.argmax(A0Range > 1e-9)-1
            A0Stop = np.argmax(A0Range > 1e-5)
        else:
            A0Start = 0
            A0Stop = len(A0Range)

        ConcRange = A0Range[A0Start:A0Stop]

        FinalFlux = np.zeros(len(SortedPathways.values()))
        Flux3D = np.zeros([len(ConcRange),len(SortedPathways.values())])

        for i in range(A0Stop-A0Start):
            P.UpdateConcentrations(Frac[A0Start+i,:], A0Range[A0Start+i])

            for k, Path in enumerate(SortedPathways.values()):
                T = P.PathIndependentFlux(Path, Probability = True)
                FinalFlux[k] = T
            Flux3D[i] = FinalFlux
        if not Redraw:
            ExtFig = plt.figure("ExtConc",facecolor="w", figsize=(10*1.8,10))
            ExtAxes1 = ExtFig.add_subplot(121, projection='3d')
            ExtAxes2 = ExtFig.add_subplot(122, projection='3d')
            ExtFig3 = plt.figure("Pathway Average",facecolor="w", figsize=(10*1.8,10))
            Ext3Axes1 = ExtFig3.add_subplot(111)
        else:
            plt.figure("ExtConc")
            ExtAxes1.cla()
            ExtAxes2.cla()
            plt.figure("Pathway Average")
            Ext3Axes1.cla()
        plt.figure("Pathway Average")
        ExtFig.suptitle("Time At Redraw: {}".format(datetime.fromtimestamp(time())))
        plt.draw()
        FluxIndex = range(1, len(FinalFlux)+1)
        # -------------------------------------------------------------
        global SaveAverage
        global SaveStdDev
        Average = []
        StdDev  = []
        for i in FluxIndex:
            Average.append(np.mean(Flux3D[:,i-1]))
            StdDev.append(np.std(Flux3D[:,i-1]))
        Average = np.array(Average)
        StdDev  = np.array(StdDev)
        SaveAverage = Average
        SaveStdDev = StdDev
        np.save("./Working/PathwayData - Averages.npy",SaveAverage)
        np.save("./Working/PathwayData - StdDev.npy",SaveStdDev)
        Ext3Axes1.bar(FluxIndex,
             Average,
             yerr = StdDev,
             align= "center",
             )
        Ext3Axes1.set_xlim((0.25,len(FluxIndex)+1))
    
        Ext3Axes1.set_xticks(range(1,len(FluxIndex)+1,5))

        left, width = .25, .5
        bottom, height = .25, .5
        right = left + width
        top = bottom + height
        Ext3Axes1.set_xlabel("Total Flux: {:.4f}".format(np.sum(Average)))

        YMax = np.sum(Average)
        Ext3Axes1.set_ylim((0.0,YMax))
        Breaks = np.unique(np.array(SortedPathways.keys())[:,0],return_counts=True,return_index=True)[2]
        Breaks = np.cumsum(Breaks[::-1])[:-1]
        Average = np.array(Average)
        Marks = map(np.sum,np.split(Average, Breaks))
        Breaks = list(Breaks)
        Breaks.append(FluxIndex[-1])
        Breaks = np.array(Breaks)+1.0-0.475
        for i in range(len(Breaks)):
            Ext3Axes1.axvline(Breaks[i], 0.0,1.0,color="grey")
            Ext3Axes1.axvline(Breaks[i], 0.0,Marks[i]/YMax,color="red")

        # -------------------------------------------------------------
        ConcLog = np.log10(ConcRange)
        FluxIndex = range(len(FinalFlux))
        plt.figure("ExtConc")
        for i, T in enumerate(ConcLog):
            color = ["b"]*len(FluxIndex)
            color[np.argmax(Flux3D[i,:])] = "r"
            ExtAxes1.bar(
                FluxIndex,
                Flux3D[i,:],
                zs=T,
                zdir="y",
                color=color,
                alpha = 0.9)
            ConcIndex = range(Frac.shape[1])
            color = ["b"]*len(ConcIndex)
            color[np.argmax(Frac[A0Start+i,:])] = "r"
            ExtAxes2.bar(
                ConcIndex,
                Frac[A0Start+i,:],
                zs=T,
                color=color,
                zdir="y",
                alpha = 0.9
                )
        if not Redraw:
            plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
        #ExtAxes1.set_zbound(0.0,1.0)
        ExtAxes1.set_ybound(ConcLog[0],ConcLog[-1])
        ExtAxes2.set_ybound(ConcLog[0],ConcLog[-1])
        ExtAxes2.set_zbound(0.0,1.0)
        
        plt.draw()
        plt.figure("Pathway Average")
        plt.draw()
        plt.figure("Main")
        return ExtFig, ExtAxes1, ExtAxes2, ExtFig3, Ext3Axes1

    ExtConcFig, ExtConcAxes1, ExtConcAxes2, ExtConcFig3, Ext3ConcPlots1 = ExtConcPlots(
                     FullOutArray[KD1Range==KD1,KD2Range==KD2,:,TimeTarget,:][0],
                     Redraw = False
                     )

# --------------------------------------------------------------------------
if ShowSecondary3DScreen:
    from mpl_toolkits.mplot3d import Axes3D
    from time import time
    def Ext3DPlots(
                        Frac,
                        Time,
                        ExtFig = None,
                        ExtAxes1 = None,
                        ExtAxes2 = None,
                        Redraw = True,
                        A0 = A0Range[A0Target]
                     ):

        TimeSpace = np.diff(Time)
        FinalFlux = np.zeros(len(SortedPathways.values()))
        WFluxSum = np.zeros(len(SortedPathways.values()))
        Flux3D = np.zeros([len(Time)-1,len(SortedPathways.values())])

        for i in range(Frac.shape[0]-1):
            P.UpdateConcentrations(Frac[i,:], A0)
            #ODE.SetConc(A0Range[A0Target]*Frac[i,:]/float(P.ProteinsInSpecies[-1]))
            #ODE.GenCoef(
            #                Kd1   = KD1,
            #                Kd2   = KD2,
            #                Kp    = Kp,
            #                Delta = 0.0
            #                )

            Total = 0.0
            WTotal = 0.0
            for k, Path in enumerate(SortedPathways.values()):
                T = P.PathIndependentFlux(Path, Probability = False)
                WT = T*TimeSpace[i]
                FinalFlux[k] = T
                WFluxSum[k] += WT
                Total += T
                WTotal += WT
            Flux3D[i] = FinalFlux#*(ODE.Der(ODE.Conc_0, None)[-1])
        if not Redraw:
            ExtFig = plt.figure("Ext",facecolor="w", figsize=(10*1.8,10))
            ExtAxes1 = ExtFig.add_subplot(121, projection='3d')
            ExtAxes2 = ExtFig.add_subplot(122, projection='3d')
        else:
            ExtAxes1.cla()
            ExtAxes2.cla()
        plt.figure("Ext")
        ExtFig.suptitle("Time At Redraw: {}".format(datetime.fromtimestamp(time())))
        plt.draw()
        
        TimeLog = np.log10(Time[:-1])
        FluxIndex = range(1,len(FinalFlux)+1)
        for i, T in enumerate(TimeLog):
            color = ["b"]*len(FluxIndex)
            color[np.argmax(Flux3D[i,:])] = "r"
            ExtAxes1.bar(
                FluxIndex,
                Flux3D[i,:],
                zs=T,
                zdir="y",
                color=color,
                alpha = 0.9)
            ConcIndex = range(Frac.shape[1])
            color = ["b"]*len(ConcIndex )
            color[np.argmax(Frac[i,:])] = "r"
            ExtAxes2.bar(
                ConcIndex,
                Frac[i,:],
                zs=T,
                color=color,
                zdir="y",
                alpha = 0.9
                )
        if not Redraw:
            plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
        #ExtAxes1.set_zbound(0.0,1.0)
        ExtAxes1.set_ybound(np.log10(TimePoints)[0],np.log10(TimePoints)[-1])
        ExtAxes2.set_ybound(np.log10(TimePoints)[0],np.log10(TimePoints)[-1])
        ExtAxes2.set_zbound(0.0,1.0)
        plt.draw()
        plt.figure("Main")
        return ExtFig, ExtAxes1, ExtAxes2

    ExtFig, ExtAxes1, ExtAxes2 = Ext3DPlots(
                     FullOutArray[KD1Range==KD1,KD2Range==KD2,A0Target,:,:][0],
                     TimePoints,
                     Redraw = False,
                     A0 = A0Range[A0Target]
                     )
# ---------------------------------------------------------------------------------
import matplotlib.image as mpimg

def PlotPathway(Img, Fig, PathwayNumber):
    plt.figure("Pathway")
    img=mpimg.imread('./Paths/Pathway{}.png'.format(PathwayNumber))
    Img.clear()
    Img.axis('off')
    Img.imshow(img,
                interpolation='none',
                )
    
    Fig.suptitle("Pathway: {}".format(PathwayNumber))
    plt.draw()
    plt.figure("Main")

if ShowPaths:
    PathsFig = plt.figure("Pathway",facecolor="w", figsize=(10*1.8,10))
    Img = PathsFig.add_axes([0,0,1,1])
    PlotPathway(Img, PathsFig, 1)

# ---------------------------------------------------------------------------------
Ext = ((KD1Range[0],KD1Range[-1],KD2Range[0],KD2Range[-1]))    # KD1 vs KD2 Plotbounds

STitle = MainFigure.suptitle("DataSet: {}\nKD1: {:.2e} KD2: {:.2e} Time: {:.1e} Conc0: {:.2e} Species: {} Tau: {}".format(WorkingName, KD1, KD2, TimePoints[TimeTarget], A0Range[A0Target], Species, Tau))

ShowContourLevels = False
ShowMaxPoints = False
PlotType = "T99"

def UpdateP():
    global KD1, KD2
    Tau = Kp * A0Range[A0Target] * TimePoints[TimeTarget]

    ContourPlt.clear()
    ContourPlt.set_yscale('log')
    ContourPlt.set_xscale('log')
    ContourPlt.set_xlim(KD1Range[0],KD1Range[-1])
    ContourPlt.set_ylim(KD2Range[0],KD2Range[-1])
    
    xtics = matplotlib.ticker.LogLocator(base=10.0, subs=[1.0], numdecs=4, numticks=5)
    ContourPlt.xaxis.set_major_locator(xtics)
    ytics = matplotlib.ticker.LogLocator(base=10.0, subs=[1.0], numdecs=4, numticks=5)
    ContourPlt.yaxis.set_major_locator(ytics)

    
#    Min = 10.0**-6.0; Max = 10.0**-10.0;
#    ContourPlt.set_xlim(Min,Max)
#    ContourPlt.set_ylim(Min,Max)

    Z = np.zeros((len(KD2Range), len(KD1Range)), dtype='f8')
    ContourPlt.set_title(PlotType+"\n")
    ZMax = np.max(FullOutArray[:,:,A0Target,TimeTarget,Species])
    for x, TKD1 in enumerate(KD1Range):
        for y, TKD2 in enumerate(KD2Range):
            IP = InflectPoints(TimePoints, FullOutArray[x,y,A0Target,:,Species], np.sum(FullOutArray[x,y,A0Target,:,:],-1))
            try:
                if PlotType == "IP Count":
                    Z[y,x] = len(IP)
                elif PlotType == "LogWidth":
                    Z[y,x] = IP[-1] - IP[0]
                elif PlotType == "ScaledMax":
                    IPEnd = IP[-1]
                    IPStart = IP[0]
                    IPMid = float((IPEnd - IPStart)/2.0)
                    IPMidTP = find_nearest(TimePoints, 10.0**IPMid)
                    T = int(np.argmax(TimePoints == IPMidTP))
                    if len(IP) > 1:
                        Z[y,x] = (FullOutArray[x,y,A0Target,T ,Species])*np.max(FullOutArray[x,y,A0Target,:,Species])
                    else:
                        Z[y,x] = np.max(FullOutArray[x,y,A0Target,:,Species])
                elif PlotType == "FracMaxAnyTime":
                    Z[y,x] = np.max(FullOutArray[x,y,A0Target,:,Species])
                elif PlotType == "FracMaxAnyA0":
                    Z[y,x] = np.max(FullOutArray[x,y,:,TimeTarget,Species])
                elif PlotType == "FracMaxAny":
                    Z[y,x] = np.max(FullOutArray[x,y,:,:,Species])
                elif PlotType == "Plateau Frac":
                    IPEnd = IP[-1]
                    IPStart = IP[0]
                    IPMid = float((IPEnd - IPStart)/2.0)
                    IPMidTP = find_nearest(TimePoints, 10.0**IPMid)
                    T = int(np.argmax(TimePoints == IPMidTP))
                    if len(IP) > 1:
                        Z[y,x] = FullOutArray[x,y,A0Target,T ,Species]
                    else:
                        Z[y,x] = 0.0
                elif PlotType == "Flux":
                    Fraction = FullOutArray[x,y,A0Target,TimeTarget,:]
                    Fraction = A0Range[A0Target]*Fraction/float(P.ProteinsInSpecies[-1])
                    ODE.SetConc(Fraction)
                    ODE.GenCoef(
                                    Kd1   = TKD1,
                                    Kd2   = TKD2,
                                    Kp    = Kp,
                                    Delta = 0.0
                                    )
                    Z[y,x] = ODE.Der(ODE.Conc_0, None)[Species]
                elif PlotType == "Frac":
                    #Z[y,x] = ((np.max(FullOutArray[x,y,A0Target,TimeTarget,Species]))/ZMax)-1.0
                    Z[y,x] = np.max(FullOutArray[x,y,A0Target,TimeTarget,Species])
                elif PlotType == "T99":
                    if len(IP) >= 1:
                        Z[y,x] = IP[0]
                    else:
                        Z[y,x] = 0.0

                    # t99Index = np.argmax(FullOutArray[x,y,A0Target,:,Species] > ZMax*0.99)
                    # if t99Index == 0:
                    #     t99Index = -1
                    # Z[y,x] = np.log10(TimePoints[t99Index])


                    #Z[y,x] = ((np.max(FullOutArray[x,y,A0Target,TimeTarget,Species]))/ZMax)-1.0
                else:
                    Z[y,x] = 1.0
            except:
                Z[y,x] = 0.0

    # ----------------------
    # For data output from Spyder
    global SaveContourData
    SaveContourData = Z
    np.save("./Working/ContourData.npy",SaveContourData)
    np.save("./Working/KDRange.npy",KD1Range)

    # ----------------------

    #ContourImage = ContourPlt.pcolormesh(KD1Range, KD2Range, Z, cmap=cm.gray)
    #ColorMap = cm.bwr
    #ColorMap = cm.nipy_spectral
    #ColorMap = cm.Spectral
    #ColorMap = cm.seismic
    ColorMap = cm.gray
    # if PlotType == "Frac":
    #     ContourImage = ContourPlt.pcolormesh(KD1Range, KD2Range, Z, cmap=ColorMap, vmin=0.0, vmax=1.0)
    # else:
    ContourImage = ContourPlt.pcolormesh(KD1Range, KD2Range, Z, cmap=ColorMap)
    ContourImage.autoscale()
    #if np.min(Z) < 0.0 or np.max(Z) > 1.0:
        
    if ShowContourLevels:
        ContourLevels = np.hstack([np.arange(np.min(Z), np.max(Z), np.max(Z)/4.0),np.max(Z)*0.99])
        CS = ContourPlt.contour(
                                KD1Range,
                                KD2Range,
                                Z,
                                ContourLevels,
                                origin='lower',
                                linewidths=2,
                                extent=Ext
                                )
        ContourPlt.clabel(CS,
                       ContourLevels,
                       inline=1,
                       fmt='%.3f',
                       fontsize=14
                       )

    
    # Deal with Colorbar Initial Drawing Issues
    global CB
    try:
        CB.ax.get_position().bounds
    except:
        #CB = plt.colorbar(ContourImage, shrink=0.8, extend='both')
        CB = plt.colorbar(ContourImage, shrink=0.8)
        l,b,w,h = plt.gca().get_position().bounds
        ll,bb,ww,hh = CB.ax.get_position().bounds
        CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])

    CB.update_bruteforce(ContourImage)
    CB.draw_all()
    ContourPlt.set_xlabel('KD1')
    ContourPlt.set_ylabel('KD2')

    STitle.set_text("DataSet: {}\nKD1: {:.2e} KD2: {:.2e} Time: {:.1e} Conc0: {:.2e} Species: {} Tau: {}".format(WorkingName, KD1, KD2, TimePoints[TimeTarget], A0Range[A0Target], Species, Tau))

    ContourPlt.plot(KD1Range, KD2Range, color='r')
    
    #L = np.abs(-10.0**-8.0+KD1Range)

    #V = L > 0.0
    #L = np.log10(L[V])

    #ContourPlt.loglog(KD1Range,L)

    #ContourPlt.plot(np.logspace(-12,-3,len(KD1Range))**2.0,np.logspace(-12,-3,len(KD1Range)))

    if ShowMaxPoints:
        A=np.where(Z==np.max(Z))
        print "Dev: {:.2e} Max {} @ KD1:{:.2e} KD2:{:.2e}".format(1.0-Z[A][0],Z[A][0], KD1Range[A[0]][0],KD2Range[A[1]][0])
    
        A=np.where(Z>np.max(Z)*0.9999)
        try:
            ContourPlt.plot(KD1Range[A[:,1]],KD2Range[A[:,0]], marker="o", color="r",linestyle='None', ms=10.0)
        except:
            ContourPlt.plot(KD1Range[A[1]],KD2Range[A[0]], marker="o", color="r",linestyle='None', ms=10.0)
        AM=np.where(Z==np.max(Z))
        ContourPlt.plot(KD1Range[AM[1]],KD2Range[AM[0]], marker="o", color="b",linestyle='None', ms=10.0)

    ContourPlt.plot(KD1,KD2, marker="o", color="cyan",linestyle='None', ms=10.0)
    plt.draw()

if ShowRingCurve:
    import Ring as RingData
    Ring = RingData.Ring3()

PlotA0 = "MaxFrac"


def UpdateAxes():
    global KD1, KD2

    if ShowRingCurve:
        global Ring
        RingKD = np.sqrt(KD1*KD2)
        Ring.SetKD(RingKD)
        Ring.SetDeg(Deg)
        RingConc = A0Range[A0Target]
        Ring.SetConc0(RingConc)

    KD1Target = np.argmax(KD1Range==KD1)
    KD2Target = np.argmax(KD2Range==KD2)


    Tau = Kp * A0Range[A0Target] * TimePoints[TimeTarget]
    FvsTime.clear()
    FvsConc.clear()

    CurrentArray = FullOutArray[KD1Target,KD2Target,:,:,Species]
    AllSpeciesCurrentArray = FullOutArray[KD1Target,KD2Target,:,:,:]
    TotalConcArray = np.sum(FullOutArray[KD1Target,KD2Target,:,:,:],-1)
    FracAtA0 = CurrentArray[A0Target, :]
    if PlotA0 == "Delta":
        Delta = np.zeros(len(A0Range), dtype='f8')
        for i in range(len(Delta)):
            IP = InflectPoints(TimePoints, CurrentArray[i, :], TotalConcArray[i,:])
            try:
                Delta[i] = IP[-1] - IP[0]
            except:
                pass
        FvsConc.semilogx(A0Range, Delta)
        FvsConc.set_xlabel("Conc (M)")
        FvsConc.set_ylabel("Delta (s)")
        FvsConc.set_ylim(0.0,np.log10(TimePoints[-1])-np.log10(TimePoints[0]))
        FvsConc.axvline(A0Range[A0Target], 0.0,1.0,color="red")
    elif PlotA0 == "MaxFrac":
        MaxFrac = np.zeros((len(A0Range),SpeciesCount), dtype='f8')
        #Valid = ValidbyTotalConc(TotalConcArray[:,TimeTarget])
        for i in range(len(A0Range)):
            for j in range(SpeciesCount):
                MaxFrac[i,j] = AllSpeciesCurrentArray[i, TimeTarget, j]
        #for j in range(SpeciesCount):
        #    FvsConc.semilogx(A0Range, MaxFrac[:,j])

        FvsConc.semilogx(A0Range, MaxFrac[:,7]+MaxFrac[:,8], alpha = 0.8, color="green")
        FvsConc.semilogx(A0Range, MaxFrac[:,10], alpha = 0.8, color="red")
        FvsConc.semilogx(A0Range, MaxFrac[:,11], alpha = 0.8, color="blue")
        FvsConc.semilogx(A0Range,
                        MaxFrac[:,0] + \
                        MaxFrac[:,1] + \
                        MaxFrac[:,2] + \
                        MaxFrac[:,3] + \
                        MaxFrac[:,4] + \
                        MaxFrac[:,5] + \
                        MaxFrac[:,6] + \
                        MaxFrac[:,7] + \
                        MaxFrac[:,8] + \
                        MaxFrac[:,9] + \
                        MaxFrac[:,10] + \
                        MaxFrac[:,11],
                        alpha = 0.8, color="orange")
        if ShowRingCurve:
            FvsConc.semilogx(A0Range, Ring.GetA0Course(TimeTarget), color="purple", alpha = 0.8)
        FvsConc.set_xlabel("Conc (M)")
        FvsConc.set_ylabel("Fraction Assembled")
        FvsConc.set_ylim(0.0, 1.05)
        FvsConc.axvline(A0Range[A0Target], 0.0,1.0,color="red")

        # # --------------------------------
        # # For variable saving in Spyder
        # global SaveTimePoints
        # global SaveStackedTimeFraction
        # global SaveRingTimeFraction
        # SaveTimePoints = TimePoints[Valid]
        # SaveStackedTimeFraction = FracAtA0[Valid]
        # SaveRingTimeFraction = Ring.GetTimeCourse(RingConc)
        np.save("./Working/FracVsConc - ConcPoints.npy",A0Range)
        np.save("./Working/FracVsConc - Frac - Stacked.npy",MaxFrac)
        if ShowRingCurve:
            np.save("./Working/FracVsConc - Frac - Ring.npy",Ring.GetA0Course(TimeTarget))
        # np.save(".\\Working\\FracVsTime - Frac - Ring.npy",SaveRingTimeFraction)
        # # --------------------------------


    elif PlotA0 == "Pathways":
        Frac = FullOutArray[KD1Target,KD2Target,A0Target,TimeTarget,:]
        PathwayPlot(Frac, A0Range[A0Target], FvsConc)
    IP = InflectPoints(TimePoints, CurrentArray[A0Target,:], TotalConcArray[A0Target,:])
    try:
        FvsTime.axvspan(10.0**IP[-1], 10.0**IP[0], facecolor='g', alpha=0.5)
    except:
        pass
    FvsTime.axvline(TimePoints[TimeTarget], 0.0,1.0,color="red")

    Valid = ValidbyTotalConc(TotalConcArray[A0Target,:])
    InvalidPointCount = np.sum(Valid != True)
    print "Invalid Points: {}".format(InvalidPointCount)
    print TimePoints[Valid!=True]
    s5 = InterpolatedUnivariateSpline(np.log10(TimePoints[Valid]), CurrentArray[A0Target,Valid], k=5)
    print "# of Inflection Points: {}".format(len(IP))
    FvsTime.plot(TimePoints, s5(np.log10(TimePoints)), color='g', alpha=0.5)
    FvsTime.plot(10.0**IP, s5(IP), "^", markersize=10.0, color='r')
    
    # --------------------------------
    # For variable saving in Spyder
    global SaveTimePoints
    global SaveStackedTimeFraction
    if ShowRingCurve:
        global SaveRingTimeFraction
        SaveRingTimeFraction = Ring.GetTimeCourse(RingConc)
        np.save("./Working/FracVsTime - Frac - Ring.npy",SaveRingTimeFraction)

    SaveTimePoints = TimePoints[Valid]
    SaveStackedTimeFraction = FracAtA0[Valid]
    np.save("./Working/FracVsTime - TimePoints.npy",SaveTimePoints)
    np.save("./Working/FracVsTime - Frac - Stacked.npy",SaveStackedTimeFraction)

    # --------------------------------
    
    FvsTime.plot(TimePoints[Valid], FracAtA0[Valid], 'x', color="b")
    FvsTime.plot(TimePoints[Valid], FullOutArray[KD1Target,KD2Target,:,:,-1][A0Target, :][Valid], color="r")
    FvsTime.plot(TimePoints[Valid], 
                FullOutArray[KD1Target,KD2Target,:,:,0][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,1][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,2][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,3][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,4][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,5][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,6][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,7][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,8][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,9][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,10][A0Target, :][Valid] + \
                FullOutArray[KD1Target,KD2Target,:,:,11][A0Target, :][Valid],
                color="orange")
    FvsTime.plot(TimePoints[Valid!=True], s5(np.log10(TimePoints[Valid!=True])), 'x', color="r") # Show in red where the missing points would be according to the spline

    if ShowRingCurve:
        FvsTime.plot(Ring.Time, Ring.GetTimeCourse(RingConc), color="purple")
    
    #FvsTime.set_ylim(np.min(FracAtA0[Valid]), np.max(FracAtA0[Valid])*1.05)
    FvsTime.set_ylim(0.0, 1.05)
    FvsTime.set_xscale("log")
    FvsTime.set_xlabel("Time (s)")
    FvsTime.set_ylabel("Fraction")
    STitle.set_text("DataSet: {}\nKD1: {:.2e} KD2: {:.2e} Time: {:.1e} Conc0: {:.2e} Species: {} Tau: {}".format(WorkingName, KD1Target, KD1Target, TimePoints[TimeTarget], A0Range[A0Target], Species, Tau))
    plt.draw()


# Make room for Sliders
plt.subplots_adjust(left=0.15, bottom=0.25)

# Draw First Time
UpdateP()
UpdateAxes()

# KDSlider causes all other sliders to reset the selected KD1,KD2 point
KDSlider = False

# Add Sliders
axSpecies = plt.axes([0.15, 0.05, 0.65/2.0, 0.03])
axTime = plt.axes([0.15, 0.05+0.03, 0.65/2.0, 0.03])
axConc_0 = plt.axes([0.15, 0.05+0.03*2, 0.65/2.0, 0.03])
if KDSlider:
    axKD = plt.axes([0.15, 0.05+0.03*3, 0.65/2.0, 0.03])
    KD = Slider(axKD, 'KD', 0, len(KD1Range), valinit=np.argmax(KD1==KD1Range),dragging =False)

SpeciesSlider = Slider(axSpecies, 'Species', 0, SpeciesCount+1, valinit=Species,dragging =False)
TimeSlider = Slider(axTime, 'Time', 0, len(TimePoints), valinit=TimeTarget,dragging =False)
Conc_0 = Slider(axConc_0, 'Conc_0', 0, len(A0Range), valinit=A0Target,dragging =False)


# Show working dataset above the sliders
#plt.title(WorkingName)

OS = 0
OA = 0
OT = 0
OK = 0

def update(val):

    global OS, OA, OT, OK

    S = find_nearest(np.arange(SpeciesCount), SpeciesSlider.val)
    T = find_nearest(np.arange(len(TimePoints)), TimeSlider.val)
    A = find_nearest(np.arange(len(A0Range)), Conc_0.val)
    if KDSlider:
        K = find_nearest(np.arange(len(KD1Range)), KD.val)
    else:
        K = OK
    if not ((OS == S) and (OA == A) and (OT == T) and (OK == K)):
        print "--------------------------------------------"
        print "Species Entry: %i"%S
        print "Initial Concentration: %e"%A0Range[A]
        print "TimePoint: %e"%TimePoints[T]
        if KDSlider:
            print "TimePoint: %e"%KD1Range[K]
        print "Update Time: {}".format(datetime.fromtimestamp(time()))

    OS = S
    OA = A
    OT = T

    global Species, A0Target, TimeTarget, KD1, KD2

    Species = S
    A0Target = A
    TimeTarget = T
    if KDSlider:
        OK = K
        KD1 = KD1Range[K]
        KD2 = KD2Range[K]
    UpdateAxes()
    UpdateP()

def onpick(event):

    global KD1, KD2

    print "Update Time: {}".format(datetime.fromtimestamp(time()))
    if event.inaxes == ContourPlt:
        print "x: {} y: {}".format(event.xdata, event.ydata)
        KD1 = event.xdata
        KD2 = event.ydata
        KD1 = find_nearest(KD1Range, KD1)
        KD2 = find_nearest(KD2Range, KD2)
        print "KD1: {:2e} KD2: {:2e}".format(KD1, KD2)
        ContourPlt.plot(KD1,KD2, marker="o", color="cyan",linestyle='None', ms=10.0)
        UpdateAxes()
    elif PlotA0 == "Pathways" and event.inaxes == FvsConc and event.xdata < 46.5:
        print "Pathway Selected: {} : {}: Showing {}".format(event.xdata,event.ydata, int(round(event.xdata)))
        #FvsConc.set_xlabel("Pathway: {}".format(int(round(event.xdata))))
        PlotPathway(Img, PathsFig, int(round(event.xdata)))



# Button Functions
def FlipCL(event):
    global ShowContourLevels
    if ShowContourLevels:
        ShowContourLevels = False
    else:
        ShowContourLevels = True
    UpdateP()

def FlipCM(event):
    global ShowMaxPoints
    if ShowMaxPoints:
        ShowMaxPoints = False
    else:
        ShowMaxPoints = True
    UpdateP()

def SetCPlotType(event, NewPlotType):
    global PlotType
    Update = False
    if PlotType != NewPlotType:
        Update = True
    PlotType = NewPlotType
    if Update:
        UpdateP()

def SetA0Plot(event, NewPlotA0):
    global PlotA0
    Update = False
    if PlotA0 != NewPlotA0:
        Update = True
    PlotA0 = NewPlotA0
    if Update:
        UpdateAxes()
        
def SetDS(event, NewDS):
    global WorkingName
    Update = False
    if WorkingName != NewDS[0]:
        Update = True
    WorkingName = NewDS[0]
    global N
    if N != NewDS[1]:
        Update = True
    N = NewDS[1]
    if Update:
        ImportData()
        UpdateP()
        UpdateAxes()

def ReExt(event):
    global KD1, KD2
    KD1Target = np.argmax(KD1Range==KD1)
    KD2Target = np.argmax(KD2Range==KD2)
    if ShowSecondary3DScreen:
        Ext3DPlots(
            FullOutArray[KD1Target,KD2Target,A0Target,:,:],
            TimePoints,
            ExtFig,
            ExtAxes1,
            ExtAxes2,            
            Redraw=True,
            )
def ReExtConc(event):
    global KD1, KD2
    KD1Target = np.argmax(KD1Range==KD1)
    KD2Target = np.argmax(KD2Range==KD2)
    if ShowSecondaryPathwaysConcScreen:
        ExtConcPlots(
            FullOutArray[KD1Target,KD2Target,A0Target,:,:],
            ExtConcFig,
            ExtConcAxes1,
            ExtConcAxes2,
            ExtConcFig3,
            Ext3ConcPlots1,
            Redraw=True,
            )
        ExtConcPlots
def ReContour(event):
    UpdateP()


ButtonQueue = []
"""
Original:
ButtonQueue.append(("CL On/Off",FlipCL))
ButtonQueue.append(("CM On/Off",FlipCM))
ButtonQueue.append(("C:T99",lambda event: SetCPlotType(event, "T99")))
#ButtonQueue.append("") # Blank Button Placeholder
ButtonQueue.append(("C:Frac",lambda event: SetCPlotType(event, "Frac")))
#ButtonQueue.append("") # Blank Button Placeholder
ButtonQueue.append(("C:IP Count",lambda event: SetCPlotType(event, "IP Count")))
ButtonQueue.append(("C:LogWidth",lambda event: SetCPlotType(event, "LogWidth")))
ButtonQueue.append(("C:FracMaxAnyTime",lambda event: SetCPlotType(event, "FracMaxAnyTime")))
ButtonQueue.append(("C:FracMaxAnyA0",lambda event: SetCPlotType(event, "FracMaxAnyA0")))
ButtonQueue.append(("C:FracMaxAny",lambda event: SetCPlotType(event, "FracMaxAny")))
ButtonQueue.append(("C:ScaledMax",lambda event: SetCPlotType(event, "ScaledMax")))
ButtonQueue.append(("C:Plateau Frac",lambda event: SetCPlotType(event, "Plateau Frac")))
ButtonQueue.append(("C:Flux",lambda event: SetCPlotType(event, "Flux")))

#ButtonQueue.append("") # Blank Button Placeholder
ButtonQueue.append(("A0 Plot:Delta",lambda event: SetA0Plot(event, "Delta")))
ButtonQueue.append(("A0 Plot:Max Frac",lambda event: SetA0Plot(event, "MaxFrac")))
ButtonQueue.append(("A0 Plot:Pathways",lambda event: SetA0Plot(event, "Pathways")))

#ButtonQueue.append("") # Blank Button Placeholder
"""
ButtonQueue.append(("CL On/Off",FlipCL))
ButtonQueue.append(("CM On/Off",FlipCM))
ButtonQueue.append(("C:Frac",lambda event: SetCPlotType(event, "Frac")))
ButtonQueue.append(("C:IP Count",lambda event: SetCPlotType(event, "IP Count")))
ButtonQueue.append(("C:ScaledMax",lambda event: SetCPlotType(event, "ScaledMax")))
ButtonQueue.append(("Redraw:Contour",lambda event: ReContour(event)))
ButtonQueue.append(("Redraw:Ext",lambda event: ReExt(event)))
ButtonQueue.append(("Redraw:ExtConc",lambda event: ReExtConc(event)))
#ButtonQueue.append(("A0 Plot:Delta",lambda event: SetA0Plot(event, "Delta")))
#ButtonQueue.append(("A0 Plot:Max Frac",lambda event: SetA0Plot(event, "MaxFrac")))
ButtonQueue.append(("A0 Plot:Pathways",lambda event: SetA0Plot(event, "Pathways")))
#ButtonQueue.append(("DS:Static",lambda event: SetDS(event, "TriStackedRingODELargeStaticT")))
#ButtonQueue.append(("DS:SynDeg",lambda event: SetDS(event, "TriStackedRingODELargeSynDeg28e-04T")))
#ButtonQueue.append(("DS:Before",lambda event: SetDS(event, "TriStackedRingODELargeStaticTBeforeFix")))
#ButtonQueue.append(("DS:After",lambda event: SetDS(event, "TriStackedRingODELargeStaticT")))
#ButtonQueue.append(("DS:Before",lambda event: SetDS(event, "TriStackedRingODELargeSynDeg28e-04TBeforeFix")))
#ButtonQueue.append(("DS:After",lambda event: SetDS(event, "TriStackedRingODELargeSynDeg28e-04T")))
'''
ButtonQueue.append(("DS:Static",lambda event: SetDS(event, ("TriStackedRingODELargeStaticT",1))))
ButtonQueue.append(("DS:1.5e-8",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg15e-08T",1))))
ButtonQueue.append(("DS:9.0e-6",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg90e-06T",1))))
#ButtonQueue.append(("DS:5.9e-6",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg59e-06T",1))))
#ButtonQueue.append(("DS:3.0e-6",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg30e-06T",1))))
#ButtonQueue.append(("DS:1.5e-6",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg15e-06T",1))))
#ButtonQueue.append(("DS:2.8e-5",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-05T",1))))
# ButtonQueue.append(("DS:2.4e-5",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg24e-05T",1))))

#ButtonQueue.append(("DS:1.2e-5",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg12e-05T",1))))
#ButtonQueue.append(("DS:2.8e-4",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",1))))
ButtonQueue.append(("DS:1.7e-3",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg17e-03T",5))))
# ButtonQueue.append(("Static:5",lambda event: SetDS(event, ("TriStackedRingODELargeStaticT",5))))
# ButtonQueue.append(("Static:4",lambda event: SetDS(event, ("TriStackedRingODELargeStaticT",4))))
'''

#ButtonQueue.append(("Static:NP",lambda event: SetDS(event, ("NoPosStaticT",2))))
'''
ButtonQueue.append(("1.7e-3:V2",lambda event: SetDS(event, ("STV2SynDeg17e-03T",1))))
ButtonQueue.append(("2.8e-4:V2",lambda event: SetDS(event, ("STV2SynDeg28e-04T",1))))
ButtonQueue.append(("1.7e-5:V2",lambda event: SetDS(event, ("STV2SynDeg17e-05T",1))))
ButtonQueue.append(("3.5e-5:V2",lambda event: SetDS(event, ("STV2SynDeg35e-05T",1))))
ButtonQueue.append(("1.1e-6:V2",lambda event: SetDS(event, ("STV2SynDeg11e-06T",1))))
ButtonQueue.append(("2.2e-6:V2",lambda event: SetDS(event, ("STV2SynDeg22e-06T",1))))
ButtonQueue.append(("4.3e-6:V2",lambda event: SetDS(event, ("STV2SynDeg43e-06T",1))))
ButtonQueue.append(("8.7e-6:V2",lambda event: SetDS(event, ("STV2SynDeg87e-06T",1))))
ButtonQueue.append(("5.4e-7:V2",lambda event: SetDS(event, ("STV2SynDeg54e-07T",1))))
ButtonQueue.append(("2.7e-7:V2",lambda event: SetDS(event, ("STV2SynDeg27e-07T",1))))
ButtonQueue.append(("Static:V2",lambda event: SetDS(event, ("STV2StaticT",1))))
'''

#STV3SynDeg14e-04T
#STV3SynDeg28e-04T
#STV3SynDeg56e-04T
#STV3SynDeg17e-05T
#STV3SynDeg35e-05T
#STV3SynDeg69e-05T
#STV3SynDeg11e-06T
#STV3SynDeg22e-06T
#STV3SynDeg43e-06T
#STV3SynDeg87e-06T
#STV3SynDeg27e-07T
#STV3SynDeg54e-07T

#ButtonQueue.append(("1.4e-4:V3",lambda event: SetDS(event, ("STV3SynDeg14e-04T",1))))
#ButtonQueue.append(("2.8e-4:V3",lambda event: SetDS(event, ("STV3SynDeg28e-04T",1))))
#ButtonQueue.append(("5.6e-4:V3",lambda event: SetDS(event, ("STV3SynDeg56e-04T",1))))
#ButtonQueue.append(("1.7e-5:V3",lambda event: SetDS(event, ("STV3SynDeg17e-05T",1))))
#ButtonQueue.append(("3.5e-5:V3",lambda event: SetDS(event, ("STV3SynDeg35e-05T",1))))
ButtonQueue.append(("6.9e-5:V3",lambda event: SetDS(event, ("STV3SynDeg69e-05T",1))))
ButtonQueue.append(("1.1e-6:V3",lambda event: SetDS(event, ("STV3SynDeg11e-06T",1))))
ButtonQueue.append(("2.2e-6:V3",lambda event: SetDS(event, ("STV3SynDeg22e-06T",1))))
ButtonQueue.append(("4.3e-6:V3",lambda event: SetDS(event, ("STV3SynDeg43e-06T",1))))
ButtonQueue.append(("8.7e-6:V3",lambda event: SetDS(event, ("STV3SynDeg87e-06T",1))))
ButtonQueue.append(("2.7e-7:V3",lambda event: SetDS(event, ("STV3SynDeg27e-07T",1))))
ButtonQueue.append(("5.4e-7:V3",lambda event: SetDS(event, ("STV3SynDeg54e-07T",1))))
ButtonQueue.append(("3.2e-8:V3",lambda event: SetDS(event, ("STV3SynDeg32e-08T",1))))
ButtonQueue.append(("1.0e-9:V3",lambda event: SetDS(event, ("STV3SynDeg10e-09T",1))))
ButtonQueue.append(("1.0e-10:V3",lambda event: SetDS(event, ("STV3SynDeg10e-10T",1))))
ButtonQueue.append(("1.0e-11:V3",lambda event: SetDS(event, ("STV3SynDeg10e-11T",1))))
ButtonQueue.append(("1.0e-12:V3",lambda event: SetDS(event, ("STV3SynDeg10e-12T",1))))
ButtonQueue.append(("1.0e-13:V3",lambda event: SetDS(event, ("STV3SynDeg10e-13T",1))))
ButtonQueue.append(("Static:V3",lambda event: SetDS(event, ("STV3StaticT",1))))


#ButtonQueue.append(("Static:ND",lambda event: SetDS(event, ("NoDisStaticT",1))))
#ButtonQueue.append(("2.8e-4:ND",lambda event: SetDS(event, ("NoDisSynDeg28e-04T",1))))
#ButtonQueue.append(("Static:NBiP",lambda event: SetDS(event, ("NoBondsInPosStaticT",1))))
#ButtonQueue.append(("2.8e-4:NBiP",lambda event: SetDS(event, ("NoBondsInPosSynDeg28e-04T",1))))
#ButtonQueue.append(("Static:NBiP2",lambda event: SetDS(event, ("NoBondsInPosStaticT",2))))
#ButtonQueue.append(("2.8e-4:NBiP2",lambda event: SetDS(event, ("NoBondsInPosSynDeg28e-04T",2))))
#ButtonQueue.append(("Static:NBiP3",lambda event: SetDS(event, ("NoBondsInPosStaticT",3))))
#ButtonQueue.append(("2.8e-4:NBiP3",lambda event: SetDS(event, ("NoBondsInPosSynDeg28e-04T",3))))
#ButtonQueue.append(("2.8e-4:NP",lambda event: SetDS(event, ("NoPosSynDeg28e-04T",2))))
#ButtonQueue.append(("2.8e-4:8/7",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",11))))
#ButtonQueue.append(("2.8e-4:5",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",5))))
# ButtonQueue.append(("2.8e-4:4",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",4))))
#ButtonQueue.append(("2.8e-4:3",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",3))))
#ButtonQueue.append(("2.8e-4:2",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",2))))

# ButtonQueue.append(("Static:8/12",lambda event: SetDS(event, ("TriStackedRingODELargeStaticT",2))))
# ButtonQueue.append(("2.8e-4:8/",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",3))))
# ButtonQueue.append(("N2.8e-4:8/7",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",1))))
# ButtonQueue.append(("2.8e-4:8/12",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",2))))
# ButtonQueue.append(("2.8e-4:8/12:2",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",3))))
# ButtonQueue.append(("2.8e-4:8/13",lambda event: SetDS(event, ("TriStackedRingODELargeSynDeg28e-04T",4))))

def ButtonClick(event):
    print "-------------------------"
    print "Clicked!!"
    print event
    print "-------------------------"

# Add Button Grid to the Lower Right
ButtonGrid = []
C = 4
R = 6
T = 0.15
CW = 0.65/2.0+0.65/4.0
L = 0.15 + 0.65/2.0
H = 0.03

ButtonQueue.reverse()
for Row in range(R):
    TRow = []
    for Col in range(C):
        A = plt.axes([L+(CW/float(C))*Col, T-H*Row, CW/float(C), H])
        if len(ButtonQueue) > 0:
            B = ButtonQueue.pop()
            if B != "":
                TRow.append(Button(A, B[0]))
                TRow[-1].on_clicked(B[1])
                TRow[-1].drawon = False
            else:
                TRow.append(Button(A, ""))
                TRow[-1].on_clicked(ButtonClick)
                TRow[-1].drawon = False
        else:
            TRow.append(Button(A, ""))
            TRow[-1].on_clicked(ButtonClick)
            TRow[-1].drawon = False
    ButtonGrid.append(TRow)

#PlotType = "FracMaxAnyA0"
PlotType = "T99"

# Setup WX events
cid = MainFigure.canvas.mpl_connect('button_press_event', onpick)
SpeciesSlider.on_changed(update)
TimeSlider.on_changed(update)
Conc_0.on_changed(update)

if KDSlider:
    KD.on_changed(update)

plt.figure("Main")
plt.show()

#if ShowSecondary3DScreen:
#    plt.figure("Ext")
#    plt.show(ExtFig)
#    plt.figure("Main")

