# General Imports
import numpy as np
from os import curdir

# Initial MPL setup
import matplotlib
import matplotlib.font_manager as font_manager
matplotlib.use('PDF')
FontFile = './Fonts/tahomabd.ttf'
FontInstance = font_manager.FontProperties(fname=FontFile)
# FontFile = './Fonts/tahomait.ttf'
# FontInstance = font_manager.FontProperties(fname=FontFile, style="italic")
matplotlib.rcParams["pdf.fonttype"] = 42 # Not sure why to do this, but apparently it helps with font embeding in PDFs
matplotlib.rcParams['font.family'] = FontInstance.get_name()

#matplotlib.use("pgf")
#pgf_with_pdflatex = {
#    "pgf.texsystem": u'pdflatex',
#    "pgf.preamble": [
#        r"\usepackage[utf8x]{inputenc}",
#        r"\usepackage[T1]{fontenc}",
#        r"\usepackage{cmbright}",
#        r"\usepackage{fontspec}",
#        r'\setmathfont{'+FontInstance.get_name()+r'}',
#        ],
##    "font.family": "san-serif",
#    "font.sans-serif": "Helvetica",
#}
#matplotlib.rcParams.update(pgf_with_pdflatex)

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = FontInstance.get_name()
matplotlib.rcParams['mathtext.it'] = FontInstance.get_name()

FontInstanceMBF = font_manager.FontProperties(fname='./Fonts/achilles3expandital.ttf')
#FontInstanceMBF = font_manager.FontProperties(fname='./Fonts/achilles3boldital.ttf')
matplotlib.rcParams['mathtext.bf'] = FontInstanceMBF.get_name()



# matplotlib.rcParams['lines.antialiased'] = True
# matplotlib.rcParams['patch.antialiased'] = True
# matplotlib.rcParams['text.antialiased'] = True


# matplotlib.rcParams['mathtext.bf'] = FontInstance.get_name() # Causes a warning about font fallback


import matplotlib.pyplot as mpl

#mpl.switch_backend("pdf")


# Default DataFiles
DefaultDataDir = './FigureData'
DefaultBreakPoints = np.load("{}/{}".format(DefaultDataDir,"PathwayPlot-BreakPoints.npy"))

# Main Default Settings
# DPI = 1200
DPI = 300
matplotlib.rcParams['font.size'] = 10

# Import Images
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

from PIL import Image
ImageScalingFactor = 1.0

def ImportImages():
    ImageDirectory = "./Images/"
    ImageFileNames = [
        "32  1 0 0.0 0 0.png",
        "48  1 1 0.0 0 0.png",
        "36  1 0 0.1 0 0.png",
        "56  1 1 1.0 0 0.png",
        "52  1 1 0.1 0 0.png",
        "50  1 1 0.0 1 0.png",
        "60  1 1 1.1 0 0.png",
        "53  1 1 0.1 0 1.png",
        "54  1 1 0.1 1 0.png",
        "51  1 1 0.0 1 1.png",
        "62  1 1 1.1 1 0.png",
        "63  1 1 1.1 1 1.png",
        "Trimer.png",
        ]
    ImageObjects = []
    for ImageFileName in ImageFileNames:
        ImageFileName = "{}{}".format(ImageDirectory,ImageFileName)
        ImageObjects.append(Image.open(ImageFileName))
        OrigWidth, OrigHeight= ImageObjects[-1].size
        ScaledX = int(ImageScalingFactor*OrigWidth)
        ScaledY = int(ImageScalingFactor*OrigHeight)
        if ImageScalingFactor != 1:
            ImageObjects[-1] = ImageObjects[-1].resize((ScaledX, ScaledY), Image.ANTIALIAS)
    return ImageObjects

"""
Image.NEAREST      # use nearest neighbour
Image.BILINEAR     # linear interpolation in a 2x2 environment
Image.BICUBIC      # cubic spline interpolation in a 4x4 environment
Image.ANTIALIAS
"""

def OpenPDF(filename, directory = "Current"):
    if directory == "Current":
        directory = curdir
    import subprocess
    subprocess.Popen(["{}/{}".format(directory ,filename)],shell=True)

class MPLGPU(object):
    def __init__(self, Blank=False, hpad=0.3, wpad=0.4, FontSize = 10, Width = "double", Aspect = 4.0/3.0, axes_aspect=4.0/3.0, subplot=111, PrimaryFig = False, *args, **kwargs):
#        super(MPLGPU, self).__init__(*args, **kwargs)
        super(MPLGPU, self).__init__()

        SizeCorrectionFactor = 1.2

        # For Nature: 89 mm (single column) and 183 mm (double column) and the full depth of the page is 247 mm
        # For Nature Communications (per the checklist): 85 mm (single column) and 180 mm (double column) and the full depth of the page is 247 mm
        if Width == "double":
            Width = SizeCorrectionFactor * 180
        elif Width == "single":
            Width = 85
        elif Width == "wide":
            Width = inches2mm(6.0)
            Aspect = 16.0/9.0
            

        #hpad = 0.3
        #wpad = 0.4
        
        mpl.rcParams['text.usetex'] = False

        if PrimaryFig:
            self._MainFig = PrimaryFig
            #CSize = self._MainFig.get_size_inches()
#            self._MainFig.set_size_inches(
#                mm2inches(Width)*float(str(subplot)[1]) + wpad*(float(str(subplot)[1])-1.0),
#                mm2inches(Width/Aspect)*float(str(subplot)[0]) + hpad*(float(str(subplot)[0])-1.0)
#            )
            if axes_aspect == "equal":
                self._Plot = self._MainFig.add_subplot(subplot, aspect = 1.0, adjustable='box')
                #self._Plot.axis("equal")
            else:
                self._Plot = self._MainFig.add_subplot(subplot, adjustable='box')
            self._MainFig.subplots_adjust(
#                left=None,
#                bottom=None,
#                right=None,
#                top=None,
                wspace=wpad,
                hspace=hpad
            )
            self.SetAxesProperties(self._Plot)
        elif not Blank:
            matplotlib.rcParams['font.size'] = FontSize
            mpl.rcParams['figure.figsize'] = (
                mm2inches(Width)*float(str(subplot)[1]) + wpad*(float(str(subplot)[1])-1.0),
                mm2inches(Width/Aspect)*float(str(subplot)[0]) + hpad*(float(str(subplot)[0])-1.0)
            )
            self._MainFig = mpl.figure(
                facecolor = 'w',
                dpi = DPI,
            )
            if axes_aspect == "equal":
                self._Plot = self._MainFig.add_subplot(subplot, aspect = 1.0, adjustable='box')
                #self._Plot.axis('equal')
            else:
                self._Plot = self._MainFig.add_subplot(subplot, adjustable='box')
            
            self.SetAxesProperties(self._Plot)
        elif Blank:
            matplotlib.rcParams['font.size'] = FontSize
            self._MainFig = mpl.figure(
                facecolor = 'w',
                dpi = DPI,
            )
            self._MainFig.set_size_inches(
                mm2inches(Width)*float(str(subplot)[1]) + wpad*(float(str(subplot)[1])-1.0),
                mm2inches(Width/Aspect)*float(str(subplot)[0]) + hpad*(float(str(subplot)[0])-1.0)
            )
            
    def SetAxesProperties(self, axes, *args, **kwargs):
        axes.tick_params(
            which = 'both', # Major and Minor ticks are affected
            direction='out',
            bottom = True,
            top = False,
            left = True,
            right = False,
            *args,
            **kwargs
        )
        axes.tick_params(
            which = 'major',
            length = 5.0,
            width = 1.0,
            pad = 6.0, # How far are the labels from the mark
            *args,
            **kwargs
        )
        axes.tick_params(
            which = 'minor',
            length = 2.5,
            width = 0.75,
            pad = 5.0,
            *args,
            **kwargs
            )

    def FracToPercent(self, value, position=False):
        """
        Function for the beautification for Axes labels 
        """
        return "{:.0f}%".format(100*value)

    def ChangeYAxis(self, func):
        self._Plot.yaxis.set_major_formatter(mpl.FuncFormatter(func))

    def SaveFig(self, filename, dir = curdir, *args, **kargs):
        self._MainFig.savefig(
            "{}/{}".format(dir, filename),
            dpi=DPI,
            *args,
            **kargs
        )

    def SetAxisLabel(self, string, axis, *args, **kwargs):
        if axis == "x":
            self._Plot.set_xlabel(string, *args, **kwargs)
        if axis == "y":
            self._Plot.set_ylabel(string, *args, **kwargs)

    def Title(self, string, *args, **kargs):
        self._Plot.set_title(string, y = 1.04, *args, **kargs)

    def XLabel(self, string, *args, **kwargs):
        self.SetAxisLabel(string, "x", *args, **kwargs)

    def YLabel(self, string, *args, **kwargs):
        self.SetAxisLabel(string, "y", *args, **kwargs)

    def SetXRange(self, minimum, maximum, *args, **kwargs):
        self._Plot.set_xlim([minimum, maximum], *args, **kwargs)

    def SetYRange(self, minimum, maximum, *args, **kwargs):
        self._Plot.set_ylim([minimum, maximum], *args, **kwargs)

    def AddLine(self, x, y, xlog = True, *args, **kwargs):
        if xlog:
            self._Plot.semilogx(x, y, *args, **kwargs)
        else:
            self._Plot.plot(x, y, *args, **kwargs)

    def AddLegend(self, *args, **kargs):

        if not ("loc" in kargs.keys()):
            kargs["loc"] = "lower right"
        if not ("frameon" in kargs.keys()):
            kargs["frameon"] = False
        self.Legend = self._Plot.legend(*args, **kargs)
        
def mm2inches(size_in_mm):
    return float(size_in_mm)/25.4

def inches2mm(size_in_in):
    return float(size_in_in)*25.4

class LinePlot(MPLGPU):
    def __init__(self, YAxisInPercent = True, OF = True, *args, **kwargs):
        super(LinePlot, self).__init__(*args, **kwargs)
        self.SetAxesProperties(self._Plot)
        if YAxisInPercent:
            self.ChangeYAxis(self.FracToPercent)
        if OF:
            self.OpenFrame()

#    def AddLegend(self, *args, **kargs):
#
#        if not ("loc" in kargs.keys()):
#            kargs["loc"] = "lower right"
#        if not ("frameon" in kargs.keys()):
#            kargs["frameon"] = False
#        self.Legend = self._Plot.legend(*args, **kargs)

    def OpenFrame(self):
        self._Plot.spines['top'].set_visible(False)
        self._Plot.spines['right'].set_visible(False)
        self._Plot.spines['bottom'].set_visible(True)
        self._Plot.spines['left'].set_visible(True)

class RankWBinaryImages(MPLGPU):
    def __init__(self, Width = "double", PrimaryFig=False, subplot=111, *args, **kwargs):
        super(RankWBinaryImages, self).__init__(Width=Width, subplot=subplot, PrimaryFig=PrimaryFig, *args, **kwargs)
        #matplotlib.rcParams['font.size'] = FontSize
        if PrimaryFig:
            self._MainFig = PrimaryFig
        #self._Plot = self._MainFig.add_subplot(subplot)
        self.SetAxesProperties(self._Plot)
        self._Plot.yaxis.tick_right()
        self.PathwayCount = 46
        self.SetPathwayData(
            np.zeros(self.PathwayCount, dtype = "float64"),
            np.array(0.0,dtype='float64')
            )
        #self._Plot.set_yticklabels([])
        #self._Plot.set_xticklabels([])

#        self.TopSeqStartHeight = 1.2
        self.TopSeqStartHeight = 1.2
        self.BottomSeqStartHeight = self.TopSeqStartHeight-0.2
        self._Plot.set_ylim(0.0,self.TopSeqStartHeight+1.0)

        self._Plot.set_xlim(0.5,self.PathwayCount+0.5)
        self.AddLabels()
        self.AddTypeSeperators()
        self.ChangeYAxis(self.FracToPercent)

    def SetPathwayData(self, values, error = 0.0, Min = False, Max = False, NoError=False):
        self.PathwayData  = np.array(values, dtype='float64')
        try:
            Min[5]
            Singular = False
        except:
            Singular = True
        if not NoError:
            if not Singular:
                self.PathwayError = np.array([Min, Max])
            else:
                self.PathwayError = np.array(error, dtype='float64')
        self.SetBreakpoints()
        self.UpdatePathwayPlot(NoError)

    def SetBreakpoints(self, break_points = DefaultBreakPoints):
        self.BreakPoints = np.array(break_points,dtype='float64')
        self.PathwayTypeSums = np.array(map(np.sum,np.split(self.PathwayData, self.BreakPoints)),dtype='float64') # 0.99 is just to make sure the color dosn't leak from behind
        self.BreakPoints = np.hstack((self.BreakPoints ,self.PathwayCount))
        

    def UpdatePathwayPlot(self, NoError=False):
        for i, Tick in enumerate(self._Plot.get_yticks()):
            if Tick > 1.0:
                self._Plot.yaxis.get_major_ticks()[i].set_visible(False)

        LastBreakPoint = 0.0
        for Break, Height in np.column_stack((self.BreakPoints, self.PathwayTypeSums)):
            self._Plot.axvspan(
                LastBreakPoint+0.5,
                Break+0.5,
                0.0,
                Height/self._Plot.get_ylim()[1],
                color=(0.9,0.0,0.0),
                alpha = 0.8,
                capstyle = "projecting",
                linewidth = 0.0,
                fill = True,
                rasterized = False,
                )
            LastBreakPoint = Break
        if NoError:
            self._Plot.bar(
                left = np.arange(1,self.PathwayCount+1),
                height = self.PathwayData,
                width = 1.0,
                linewidth = 0.05,
                edgecolor = 'black',
                color = (25.0/255.0,80.0/255.0,210.0/255.0),
                ecolor = "black",
                align = "center",
                capsize = 2.0,
                alpha = 1.0
                )
        else:
            self._Plot.bar(
                left = np.arange(1,self.PathwayCount+1),
                height = self.PathwayData,
                width = 1.0,
                linewidth = 0.05,
                edgecolor = 'black',
                color = (25.0/255.0,80.0/255.0,210.0/255.0),
                yerr = self.PathwayError,
                ecolor = "black",
                align = "center",
                capsize = 2.0,
                alpha = 1.0
                )
            

    def AddLabels(self, *args, **kwargs):
        self._Plot.text(
            -1.0,
            self._Plot.get_ylim()[1]-0.5,
            'Final Reaction',
            ha='center',
            va='center',
            rotation='vertical',
            *args,
            **kwargs
            )

        self._Plot.text(
            -1.0,
            0.5,
            'Pathway Contribution',
            ha='center',
            va='center',
            rotation='vertical',
            *args,
            **kwargs
            )

        self._Plot.set_xlabel(
            "Assembly Pathway",
            *args,
            **kwargs
            )
        
    def AddTypeSeperators(self):

        LineWidth = 0.75

        self._Plot.plot(
            self._Plot.get_xlim(),
            [self.BottomSeqStartHeight,self.BottomSeqStartHeight],
            ":",
            color = "black",
            linewidth = LineWidth/2.0,
            alpha = 0.5,
            )
        self._Plot.plot(
            self._Plot.get_xlim(),
            [self.TopSeqStartHeight,self.TopSeqStartHeight],
            ":",
            color = "black",
            linewidth = LineWidth/2.0,
            alpha = 0.5,
            )

        Index = 1.0
        EqualSegWidth = self._Plot.get_xlim()[1]/float(len(self.BreakPoints))
        for Break, Height in np.column_stack((self.BreakPoints, self.PathwayTypeSums)):
            self._Plot.vlines(
                (Index*EqualSegWidth),
                ymin=self.TopSeqStartHeight,
                ymax=self._Plot.get_ylim()[1],
                color="black",
                linewidth = LineWidth
                )
            self._Plot.vlines(
                Break+0.5,
                ymin=0.0,
                ymax=self.BottomSeqStartHeight,
                color="black",
                linewidth = LineWidth
                )
            self._Plot.plot(
                [Break+0.5,(Index*EqualSegWidth)],
                [self.BottomSeqStartHeight,self.TopSeqStartHeight],
                linewidth = LineWidth,
                color="black",
                )
            Index += 1.0

class ContourPlot(MPLGPU):
    def __init__(self, Aspect=1.0, axes_aspect = "equal", *args, **kwargs):
        super(ContourPlot, self).__init__(Aspect = Aspect, axes_aspect = axes_aspect, *args, **kwargs)
        self.SetLogLog()

    def SetLogLog(self):
        self._Plot.loglog()

    def ContourData(self):
        pass

import matplotlib.pyplot as plt
import numpy as np

#Function to add ticks
def AddTicks(axis, newLocs, newLabels, pos='x'):
    # Draw to get ticks
    mpl.draw()

    # Get existing ticks
    if pos=='x':
        locs = axis.get_xticks().tolist()
        labels=[x.get_text() for x in axis.get_xticklabels()]
    elif pos =='y':
        locs = axis.get_yticks().tolist()
        labels=[x.get_text() for x in axis.get_yticklabels()]
    else:
        print("WRONG pos. Use 'x' or 'y'")
        return

    # Build dictionary of ticks
    Dticks=dict(zip(locs,labels))

    # Add/Replace new ticks
    for Loc,Lab in zip(newLocs,newLabels):
        Dticks[Loc]=Lab

    # Get back tick lists
    locs=list(Dticks.keys())
    labels=list(Dticks.values())

    # Generate new ticks
    if pos=='x':
        axis.set_xticks(locs)
        axis.set_xticklabels(labels)
    elif pos =='y':
        axis.set_yticks(locs)
        axis.set_yticklabels(labels)

if __name__ == "__main__":
    pass