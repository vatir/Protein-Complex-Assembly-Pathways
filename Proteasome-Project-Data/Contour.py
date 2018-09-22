from MPLGUI import ContourPlot

__author__ = 'Koan'

import numpy as np
import MPLGUI as MP
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

Images = MP.ImportImages()

if __name__ ==	"__main__":


	IP = False
	InVivo = False


	if InVivo:
		ImgFileName = "InVivo-Contour.pdf"
		ContourData	= np.load("{}\\{}".format(MP.DefaultDataDir, "In Vivo-Contour-Data.npy"))
	else:
		ImgFileName = "InVitro-Contour.pdf"
		ContourData	= np.load("{}\\{}".format(MP.DefaultDataDir, "Static-FracMaxAnyTime-ContourData.npy"))

	if IP:
		if InVivo:
			ImgFileName = "IP-InVivo-Contour.pdf"
			ContourData	= np.load("{}\\{}".format(MP.DefaultDataDir, "InVivo - IP - 1 - ContourData.npy"))
			ContourData	= np.load("{}\\{}".format(MP.DefaultDataDir, "InVivo - IP - 2 - ContourData.npy"))
		else:
			ImgFileName = "IP-InVivo-Contour.pdf"
			ContourData	= np.load("{}\\{}".format(MP.DefaultDataDir, "Static-IP-ContourData.npy"))

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
	if IP:
		ColorMap = cm.gray_r
	else:
		ColorMap = cm.gray

	ContourImage = MainPlot._Plot.pcolormesh(
		KD1Range,
		KD2Range,
		ContourData,
		cmap=ColorMap,
		# vmin=0.0,
		# vmax=1.0,
	)
	if not IP:
		ContourImage.autoscale()
		if InVivo:
			ContourLevels = [
					0.40,
					0.60,
					0.63,
					0.65,
					]
		else:
			ContourLevels = [
					0.50,
					0.59,
					0.75,
					0.99,
					]

		ContourLines =  MainPlot._Plot.contour(
			KD1Range,
			KD2Range,
			ContourData,
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

		MainPlot._Plot.clabel(
			ContourLines,
			ContourLevels,
			inline=1,
			fmt='%.2f',
			fontsize=7,
			orientation='vertical',
			inline_spacing = 30.0,
		)

	MainPlot._Plot.tick_params(
			which = 'major',
			pad = 1.5, # How far are the labels from the mark
		)

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	MainPlot.XLabel(r"Intra-Ring K$_d$ (M)")
	MainPlot.YLabel(r"Inter-Ring K$_d$ (M)")

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	MP.OpenPDF(ImgFileName)
