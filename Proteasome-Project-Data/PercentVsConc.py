__author__ = 'Koan'

import numpy as np
import MPLGUI as MP
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

Images = MP.ImportImages()

if __name__ ==	"__main__":

	ImgFileName = "FracVsConc.pdf"

	ConcPoints	= np.load("{}\\{}".format(MP.DefaultDataDir, "FracVsConc - ConcPoints.npy"))
	FracRing	= np.load("{}\\{}".format(MP.DefaultDataDir, "FracVsConc - Frac - Ring.npy"))
	FracStacked	= np.load("{}\\{}".format(MP.DefaultDataDir, "FracVsConc - Frac - Stacked.npy"))

	"""
	Species order for FracStacked is: (Note image order is not the same as this)
	0  : 32;;1,0,0;0,0,0
	1  : 36;;1,0,0;1,0,0
	2  : 48;;1,1,0;0,0,0
	3  : 50;;1,1,0;0,1,0
	4  : 51;;1,1,0;0,1,1
	5  : 52;;1,1,0;1,0,0
	6  : 53;;1,1,0;1,0,1
	7  : 54;;1,1,0;1,1,0
	8  : 56;;1,1,1;0,0,0
	9  : 60;;1,1,1;1,0,0
	10 : 62;;1,1,1;1,1,0
	11 : 63;;1,1,1;1,1,1
	"""

	FracStable	= np.sum(FracStacked[:,(7,8,10)], -1)
	FracStacked	= FracStacked[:,-1]

	XLim = [10.0**-11.0, 10.0**-5.0]

	Bounds = (XLim[0] <= ConcPoints) & (ConcPoints<= XLim[1])

	MainPlot = MP.LinePlot(Width="single")
	MainPlot.AddLine(ConcPoints[Bounds], FracRing[Bounds], label="      ", color="red", clip_on=False, zorder=1)
	MainPlot.AddLine(ConcPoints[Bounds], FracStacked[Bounds], label="      ", color="Blue", clip_on=False, zorder=1)
	MainPlot.AddLine(ConcPoints[Bounds], FracStable[Bounds], label="      ", color="Orange", clip_on=False, zorder=0)

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
	ImageZoom = 0.05
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
	for ImageIndex in [3,8,10,11]:
		#ImageLoc = [LegendEdges.x0 + LegendAdjust[0], LegendEdges.y1 + LegendAdjust[1]]
		ab = AnnotationBbox(
			OffsetImage(
				Images[ImageIndex],
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

	# MainPlot._Plot.text(
	# 	0.07,
	# 	1.0,
	# 	DataSet,
	# 	transform=MainPlot._Plot.transAxes,
	# 	fontsize=14,
	# 	verticalalignment='top',
	# 	)

	MainPlot.XLabel("Time (s)")
	MainPlot.YLabel("Percent in Complete Complex")

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	MP.OpenPDF(ImgFileName)
