from MPLGUI import ContourPlot

__author__ = 'Koan'

import numpy as np
import MPLGUI as MP
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

Images = MP.ImportImages()

if __name__ ==	"__main__":


	ImgFileName = "IP w t99 - Contour.pdf"
	ContourData	= np.load("{}\\{}".format(MP.DefaultDataDir, "t99 - Static - Conc0 - 1e-7 - ContourData.npy"))
	ContourLineData	= np.load("{}\\{}".format(MP.DefaultDataDir, "Static-IP-ContourData.npy"))
	ContourData[ContourData == 0.0] = np.max(ContourData)
	KD1Range	= np.load("{}\\{}".format(MP.DefaultDataDir, "KDRange.npy"))
	KD2Range	= np.load("{}\\{}".format(MP.DefaultDataDir, "KDRange.npy"))

	XLim = [np.max(KD1Range),np.min(KD1Range)]
	YLim = [np.max(KD2Range),np.min(KD2Range)]

	XBounds = (XLim[0] <= KD1Range) & (KD1Range<= XLim[1])
	XBounds = (YLim[0] <= KD2Range) & (KD2Range<= YLim[1])


	MainPlot = MP.ContourPlot(Width="single")

	MainPlot.SetXRange(XLim[0],XLim[1])
	MainPlot.SetYRange(YLim[0],YLim[1])

	XTicks = MainPlot._Plot.xaxis.get_major_ticks()
	YTicks = MainPlot._Plot.yaxis.get_major_ticks()
	for Index in np.arange(0, len(XTicks),2):
		XTicks[Index].label1On = ""
		YTicks[Index].label1On = ""

	import matplotlib.cm as cm
	ColorMap = cm.gray_r

	ContourImage = MainPlot._Plot.pcolormesh(
		KD1Range,
		KD2Range,
		ContourData,
		cmap=ColorMap,
		# vmin=0.0,
		# vmax=1.0,
	)
	ContourImage.autoscale()

	ContourLevels = [
			1.0,
			2.0,
			3.0,
			4.0,
			]

	ContourLines =  MainPlot._Plot.contour(
		KD1Range,
		KD2Range,
		ContourLineData,
		ContourLevels,
		origin='lower',
		linewidths=1.0,
		colors=(
			'purple',
			'green',
			'blue',
			'r',
		)
	)

	# MainPlot._Plot.clabel(
	# 	ContourLines,
	# 	ContourLevels,
	# 	inline=1,
	# 	fmt='%d',
	# 	fontsize=7,
	# 	orientation='vertical',
	# 	inline_spacing = 30.0,
	# )

	MainPlot._Plot.tick_params(
			which = 'major',
			pad = 1.5, # How far are the labels from the mark
		)

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	MainPlot.XLabel(r"Intra-Ring K$_d$ (M)")
	MainPlot.YLabel(r"Inter-Ring K$_d$ (M)")

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	MP.OpenPDF(ImgFileName)
