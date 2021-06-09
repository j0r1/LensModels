from grale.constants import *
import grale.inversion as inversion
import grale.renderers as renderers
import grale.plotutil as plotutil
from grale.cosmology import Cosmology
import grale.images as images
import grale.lenses as lenses
import numpy as np
import pprint
import os

#tdType = "Paper2009"
tdType = "NoSrc"

V = lambda x, y: np.array([x,y], dtype=np.double)

renderers.setDefaultLensPlaneRenderer("threads")
renderers.setDefaultMassRenderer("threads")
inversion.setDefaultInverter("threads")
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

# First two lines contain redshift and cosmology
with open("images_a1689.txt", "rt") as f:
    l1 = f.readline().strip()
    l2 = f.readline().strip()

if not l1.startswith("#") or not l2.startswith("#"):
    raise Exception("Expecting first lines to start with '#'")
l1, l2 = l1[1:].strip(), l2[1:].strip()

if not l1.startswith("zd="):
    raise Exception("Couldn't detect lens redshift")

z_lens = float(l1[3:])

cosm = { }
for p in l2.split():
    key, value = p.split("=")
    cosm[key] = float(value)

cosm = Cosmology(**cosm)
print(f"INFO: z_lens = {z_lens}")
print(f"INFO: cosmology = {cosm.getParameters()}")

def la(line):
    x, y, td, z = map(float, line.split())
    e = { "x": x*ANGLE_ARCSEC, "y": y*ANGLE_ARCSEC, "z": z }
    e["timedelay"] = td
    return e

imgList = images.readInputImagesFile("images_a1689.txt", True, lineAnalyzer=la)
allPts = np.array([ p["position"] for i in imgList for img in i["imgdata"].getAllImagePoints() for p in img ])
#pprint.pprint(allPts)

minX, maxX, minY, maxY = allPts[:,0].min(), allPts[:,0].max(), allPts[:,1].min(), allPts[:,1].max()
inversionSize = max(maxX-minX, maxY-minY)*1.3
nullSize = inversionSize * 3
center = V((minX+maxX)/2.0, (minY+maxY)/2.0)

print(f"INFO: inversion size: {inversionSize/ANGLE_ARCSEC}")
print(f"INFO: center: {center/ANGLE_ARCSEC}")

nullGrid = images.createGridTriangles(center-inversionSize/2, center+inversionSize/2, 48, 48)

iws = inversion.InversionWorkSpace(z_lens, inversionSize, regionCenter=center, cosmology=cosm)

for i in imgList:
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")
    iws.addImageDataToList(nullGrid, i["z"], "pointnullgrid")

iws.setDefaultInversionArguments(geneticAlgorithmParameters = { "selectionpressure": 2 }, sheetSearch = "genome", 
                                 fitnessObjectParameters = { "fitness_timedelay_type": tdType },
                                 
                                 maximumGenerations=2)

prevLens, baseStep, baseLens, baseFitness = None, None, None, [float("inf"), float("inf")]
subDiv = 100

for i in range(1,6):

    iws.setUniformGrid(15) if prevLens is None else iws.setSubdivisionGrid(prevLens, subDiv, subDiv+100)
    lens, fitness, fitdesc = iws.invert(512)
    print(f"Step {i} has fitness:", fitness)
    print(f"Fitness description:", fitdesc)
    lens.save(f"inv{i}.lensdata")
    
    prevLens = lens

    if fitness < baseFitness:
        baseFitness = fitness
        baseLens = lens
        baseStep = i

    subDiv += 200

lensInfo = plotutil.LensInfo(baseLens, size=60*ANGLE_ARCSEC)
corrMassScale = lensInfo.getIntegratedMass()/10.0
print("Mass scale for corrections: {:g} solar masses".format(corrMassScale/MASS_SUN))

iws.setUniformGrid(64)
corrections, corrFitness, fitdesc = iws.invert(512, baseLens=baseLens, allowNegativeValues=True, massScale=corrMassScale)
corrections.save("corrections.lensdata")

correctedLens = lenses.CompositeLens(baseLens.getLensDistance(), [
                                     { "factor": 1, "x": 0, "y": 0, "angle": 0, "lens": baseLens },
                                     { "factor": 1, "x": 0, "y": 0, "angle": 0, "lens": corrections } ])
correctedLens.save(f"correctedlens_{baseStep}.lensdata")

print("Base fitness:", baseFitness)
print("Corrected fitness:", corrFitness)

