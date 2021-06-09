from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.plotutil as plotutil
from grale.cosmology import Cosmology
import grale.images as images
import numpy as np
import os

V = lambda x, y: np.array([x,y], dtype=np.double)

subRegSize = 2.5*ANGLE_ARCSEC
subNumUniform = 15

subTotalMass = 1e10*MASS_SUN

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

# The last images data set is the one that's split by an extra galaxy
# For this one, we'll try an extra null space constraint as well, over
# a smaller area
iws.addImageDataToList(images.ImagesData.load("null_splitimg_small.imgdata"), 1.5, "pointnullgrid")

iws.setDefaultInversionArguments(sheetSearch = "genome")

def setBasisFunctions(lens, minSub, maxSub):
    iws.clearBasisFunctions()
    
    iws.setUniformGrid(15) if not lens else iws.setSubdivisionGrid(lens, minSub, maxSub)
    iws.addBasisFunctionsBasedOnCurrentGrid()

    iws.setUniformGrid(subNumUniform, regionSize = subRegSize, regionCenter = V(-5,-18.0)*ANGLE_ARCSEC)
    iws.addBasisFunctionsBasedOnCurrentGrid(initialParameters = { "totalmass": subTotalMass })

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

