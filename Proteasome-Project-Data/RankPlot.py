__author__ = 'Koan'

import numpy as np
import MPLGUI as MP
from matplotlib.offsetbox import AnnotationBbox, OffsetImage

if __name__ ==	"__main__":

	ImgFileName = "RankPlot.pdf"

	SynDegRatePoints	= np.load("{}\\{}".format(MP.DefaultDataDir, "RankData - SynDegRatePoints.npy"))
	Results	= np.load("{}\\{}".format(MP.DefaultDataDir, "RankData.npy"))
	KD1Range	= np.load("{}\\{}".format(MP.DefaultDataDir, "KDRange.npy"))

	Data = []
	T = np.empty([len(KD1Range),len(KD1Range)],dtype=object)
	for i, p in enumerate(KD1Range):
		for j, k in enumerate(KD1Range):
			T[i,j] = (p,k)

	KDs = np.ravel(T)

	MainPlot = MP.LinePlot(
		Width="single",
		# Width="double",
		YAxisInPercent=False,
		# Aspect=(16.0/9.0),
	)

	for i in range(len(Results[0]))[::10]:
		MainPlot.AddLine(
			SynDegRatePoints,
			Results[:,i],
			color="grey",
			linewidth=0.05,
			clip_on=False,
			alpha=0.1,
			zorder=0
		)
	for i in range(len(Results[0])):
		LineWidth=2.0
		if KDs[i] == (1e-12,1e-12):
			MainPlot.AddLine(
				SynDegRatePoints,
				Results[:,i],
				color="green",
				linewidth=LineWidth,
				clip_on=True,
				zorder=1
			)
		elif KDs[i] == (1e-3,1e-3):
			MainPlot.AddLine(SynDegRatePoints, Results[:,i],color="purple",linewidth=LineWidth, clip_on=True, zorder=1)
		elif KDs[i] == (1e-3,1e-12):
			MainPlot.AddLine(SynDegRatePoints, Results[:,i],color="blue",linewidth=LineWidth, clip_on=False, zorder=1)
		elif KDs[i] == (1e-12,1e-3):
			MainPlot.AddLine(SynDegRatePoints, Results[:,i],color="red",linewidth=LineWidth, clip_on=True, zorder=1)
		else:
			pass


	XLim = [SynDegRatePoints[0],SynDegRatePoints[-1]]


	# Add Static Label
	MP.AddTicks(MainPlot._Plot, [XLim[1]],[unichr(0x221E)],'x')

	Bounds = (XLim[0] <= SynDegRatePoints) & (SynDegRatePoints<= XLim[1])

	MainPlot.SetXRange(XLim[0],XLim[1])

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	
	MainPlot.XLabel("Cell Division Rate (s)")
	MainPlot.YLabel("Assembly Efficiency Rank")

	MainPlot.SaveFig(ImgFileName, bbox_inches = 'tight')
	MP.OpenPDF(ImgFileName)
