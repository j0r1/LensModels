from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.plotutil as plotutil
from grale.cosmology import Cosmology
import grale.images as images
import numpy as np
import os

V = lambda x, y: np.array([x,y], dtype=np.double)

renderers.setDefaultLensPlaneRenderer("threads")
renderers.setDefaultMassRenderer("threads")
inversion.setDefaultInverter("threads")
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

z_lens = 0.4
iws = inversion.InversionWorkSpace(z_lens, 250*ANGLE_ARCSEC, cosmology=Cosmology(0.7, 0.3, 0, 0.7))

imgList = images.readInputImagesFile("images.txt", True) 
for i in imgList:
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
    iws.addImageDataToList(images.ImagesData.load("null_large.imgdata"), i["z"], "pointnullgrid")

iws.setDefaultInversionArguments(sheetSearch = "genome")

def setBasisFunctions(lens, minSub, maxSub):
    iws.clearBasisFunctions()
    iws.setUniformGrid(15) if not lens else iws.setSubdivisionGrid(lens, minSub, maxSub)
    iws.addBasisFunctionsBasedOnCurrentGrid()

bestStep, bestLens, bestFitness = None, None, None
prevLens = None
subDiv = 100
for i in range(1,6):
    if prevLens is None:
        setBasisFunctions(None, None, None)
    else:
        setBasisFunctions(prevLens, subDiv, subDiv+100)
    subDiv += 200

    lens, fitness, fitdesc = iws.invertBasisFunctions(512)
    lens.save(f"inv{i}.lensdata")
    prevLens = lens

    if bestFitness is None or fitness < bestFitness:
        bestFitness = fitness
        bestLens = lens
        bestStep = i

bestLens.save(f"best_step{bestStep}.lensdata")

