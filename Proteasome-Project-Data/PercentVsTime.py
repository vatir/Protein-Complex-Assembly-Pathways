__author__ = 'Koan'

import numpy as np
import MPLGUI as MP
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

Images = MP.ImportImages()

if __name__ ==	"__main__":

	DataSet = "In Vivo"
	DataSet = "In Vitro"

	ImgFileName = "FracVsTime-{}.pdf".format(DataSet)

	TimePoints	= np.load("{}\\{}".format(MP.DefaultDataDir, "{} - FracVsTime - TimePoints.npy".format(DataSet)))
	FracRing	= np.load("{}\\{}".format(MP.DefaultDataDir, "{} - FracVsTime - Frac - Ring.npy".format(DataSet)))
	FracStacked	= np.load("{}\\{}".format(MP.DefaultDataDir, "{} - FracVsTime - Frac - Stacked.npy".format(DataSet)))

	XLim = [10.0**0.0, 10.0**9.0]

	Bounds = (XLim[0] <= TimePoints) & (TimePoints<= XLim[1])

	
	MainPlot = MP.LinePlot(Width="single")
	MainPlot.AddLine(TimePoints[Bounds], FracRing[Bounds], label="      ", color="red", clip_on=False, zorder=0)
	MainPlot.AddLine(TimePoints[Bounds], FracStacked[Bounds], label="      ", color="Blue", clip_on=False, zorder=0)

	MainPlot.AddLegend(
		loc = "lower right",
		#loc = "upper right",
		labelspacing = 2.5,
	)

	MainPlot.SetXRange(XLim[0],XLim[1])

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	
	#T1 =  MainPlot.Legend.get_texts()[0].get_transform()

	#ImageLoc = [9.0*10.0**5.0,0.123]
	LegendAdjust = (15.0, -10.0)
	ImageZoom = 0.07
	LegendEdges = MainPlot.Legend.get_window_extent()
	ImageLoc = [LegendEdges.x0 + LegendAdjust[0], LegendEdges.y1 - np.abs(LegendEdges.y0-LegendEdges.y1)/2.0 + LegendAdjust[1]+5.0]
	ab = AnnotationBbox(
		OffsetImage(
		Images[-1],
		zoom = ImageZoom,
		resample = True,
		dpi_cor = False,
		# Interpolation Options : 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'
		interpolation = "gaussian",
		),
		ImageLoc,
		xybox= ImageLoc,
		xycoords='axes points',
		#boxcoords="offset points",
		pad=0.0,
		frameon = False,
		)
	MainPlot._Plot.add_artist(ab)

	#ImageLoc = [9.0*10.0**5.0,0.23]
	
	
	ImageLoc = [LegendEdges.x0 + LegendAdjust[0], LegendEdges.y0 + LegendAdjust[1]]
	
	#ImageLoc = [LegendEdges.x0 + LegendAdjust[0], LegendEdges.y1 + LegendAdjust[1]]
	ab = AnnotationBbox(
		OffsetImage(
		Images[-2],
		zoom = ImageZoom,
		resample = True,
		dpi_cor = False,
		# Interpolation Options : 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'
		interpolation = "gaussian",
		),
		ImageLoc,
		xybox= ImageLoc,
		xycoords='axes points',
		
		#boxcoords="offset points",
		pad=0.0,
		frameon = False,
		)
	MainPlot._Plot.add_artist(ab)

	MainPlot._Plot.text(
		0.07,
		1.0,
		DataSet,
		transform=MainPlot._Plot.transAxes,
		fontsize=14,
		verticalalignment='top',
		)

	MainPlot.XLabel("Time (s)")
	MainPlot.YLabel("Percent in Complete Complex")

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	MP.OpenPDF(ImgFileName)
