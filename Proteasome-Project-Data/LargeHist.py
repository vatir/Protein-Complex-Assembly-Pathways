__author__ = 'Koan'

from MPLGUI import *

if __name__ ==	"__main__":
	import numpy as np

	DataSet = "WW"
	MainPlot = RankWBinaryImages(Aspect=(16.0/9.0))
	MainPlot.SetPathwayData(
		np.load("{}/{}".format(DefaultDataDir, "PathwayPlot-Average-{}.npy".format(DataSet))),
		np.load("{}/{}".format(DefaultDataDir, "PathwayPlot-StdDev-{}.npy".format(DataSet)))
		)

	#  0 :	"32  1 0 0.0 0 0.png"
	#  1 : 	"48  1 1 0.0 0 0.png"
	#  2 : 	"36  1 0 0.1 0 0.png"
	#  3 : 	"56  1 1 1.0 0 0.png"
	#  4 : 	"52  1 1 0.1 0 0.png"
	#  5 : 	"50  1 1 0.0 1 0.png"
	#  6 : 	"60  1 1 1.1 0 0.png"
	#  7 : 	"53  1 1 0.1 0 1.png"
	#  8 : 	"54  1 1 0.1 1 0.png"
	#  9 : 	"51  1 1 0.0 1 1.png"
	# 10 : 	"62  1 1 1.1 1 0.png"
	# 11 : 	"63  1 1 1.1 1 1.png"

	ImageOrderTop = [
		3,4,5,6,8,10,
		]
	ImageOrderBottom = [
		3,4,5,1,2,0,
		]
	Images = ImportImages()
	StepWidth = 46.0/6.0
	#ImageZoom = 0.09
	ImageZoom = 0.065
	XOffset = 0.15
	TopRowPos		=	1.95
	# BottomRowPos	=	1.45
	BottomRowPos	=	1.425
	for Index in np.arange(1.0,7.0):
		MainPlot._Plot.plot(
			StepWidth*Index-StepWidth/2.0+XOffset,
			(TopRowPos+BottomRowPos)/2.0,
			'+',
			color	=	"black",
			mew		=	1.5,
			ms		=	10.0,
		)
		xy = [StepWidth*Index-StepWidth/2.0+XOffset, TopRowPos]
		ab = AnnotationBbox(
			OffsetImage(
				Images[ImageOrderTop[int(Index-1)]],
				zoom = ImageZoom,
				resample = True,
				dpi_cor = False,
				# Interpolation Options : 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'
				interpolation = "gaussian",
				),
			xy,
			xybox=(1.0, 1.0),
			xycoords='data',
			boxcoords="offset points",
			pad=0.0,
			frameon = False,

			)
		MainPlot._Plot.add_artist(ab)

	for Index in np.arange(1.0,7.0):
		xy = [StepWidth*Index-StepWidth/2.0+XOffset, BottomRowPos]
		ab = AnnotationBbox(
			OffsetImage(
				Images[ImageOrderBottom[int(Index-1)]],
				zoom = ImageZoom,
				resample = True,
				dpi_cor = False,
				# Interpolation Options : 'none', 'nearest', 'bilinear', 'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos'
				interpolation = "gaussian",
				),
			xy,
			xybox=(1.0, 1.0),
			xycoords='data',
			boxcoords="offset points",
			pad=0.0,
			frameon = False,

			)
		MainPlot._Plot.add_artist(ab)


	#mpl.tight_layout(rect = [0, 0, 0.4, 1])
	ImgFileName = "LargeHist-{}.pdf".format(DataSet)

	MainPlot.SaveFig(
		ImgFileName,
		bbox_inches = 'tight',
	)
	OpenPDF(ImgFileName)
