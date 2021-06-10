from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.plotutil as plotutil
import grale.images as images
import grale.lenses as lenses
import grale.cosmology as cosmology
import numpy as np
import os
import glob
import pprint

z_lens = 0.4
cosm = cosmology.Cosmology(0.7, 0.27, 0, 0.73)
cosmology.setDefaultCosmology(cosm)

iws = inversion.InversionWorkSpace(z_lens, 200*ANGLE_ARCSEC)

nullSubDiv = 48
nullSize = 400 * ANGLE_ARCSEC
nullHW = nullSize/2
nullImgData = images.createGridTriangles([-nullHW, -nullHW], [nullHW, nullHW], nullSubDiv, nullSubDiv)

V = lambda x, y: np.array([x,y], dtype=np.double)

renderers.setDefaultLensPlaneRenderer("threads")
renderers.setDefaultMassRenderer("threads")
inversion.setDefaultInverter("threads")
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

# Use the point image data to get the redshifts
for i in images.readInputImagesFile("images.txt", True):
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
    iws.addImageDataToList(nullImgData, i["z"], "pointnullgrid")

iws.setDefaultInversionArguments(sheetSearch = "genome")

def setBasisFunctions(lens, minSub = None, maxSub = None):
    iws.clearBasisFunctions()
    iws.setUniformGrid(15) if not lens else iws.setSubdivisionGrid(lens, minSub, maxSub)
    iws.addBasisFunctionsBasedOnCurrentGrid()

prevLens, subDivStart = None, 100
bestLens, bestFitness, bestIdx = None, None, 0
for idx in range(1,6):

    setBasisFunctions(prevLens, subDivStart, subDivStart+100) # subdiv is ignored if prevLens is None
    lens, fitness, fitdesc = iws.invertBasisFunctions(512)
    
    fileName = f"inv{idx}.lensdata" 
    lens.save(fileName)
    print(f"LENSFITNESS: {fileName} {fitness}")

    prevLens = lens
    subDivStart += 200

    if bestLens is None or fitness < bestFitness:
        bestLens = lens
        bestFitness = fitness
        bestIdx = idx

print(f"BESTFITNESS: {bestFitness}")
bestLens.save(f"best_is_{bestIdx}.lensdata")

