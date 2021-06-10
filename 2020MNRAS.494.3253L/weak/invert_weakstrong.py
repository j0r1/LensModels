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
import fmtutil

z_lens = 0.4
cosm = cosmology.Cosmology(0.7, 0.27, 0, 0.73)
cosmology.setDefaultCosmology(cosm)

iws = inversion.InversionWorkSpace(z_lens, 200*ANGLE_ARCSEC)

V = lambda x, y: np.array([x,y], dtype=np.double)

renderers.setDefaultLensPlaneRenderer("threads")
renderers.setDefaultMassRenderer("threads")
inversion.setDefaultInverter("threads")
plotutil.setDefaultAngularUnit(ANGLE_ARCSEC)

# Use the point image data to get the redshifts
for i in images.readInputImagesFile("images.txt", True):
    iws.addImageDataToList(i["imgdata"], i["z"], "pointimages")

# The "Z" is used as a wildcard
weakImgData = "shear_mr_48x48_exact_zZ.imgdata"
#weakImgData = "shear_ell_48x48_pointavg_25_zZ.imgdata"
#weakImgData = "shear_zZ_gaussian_48x48_sigma60.000_nocentral120.imgdata"

fileNamesAndRedshifts = fmtutil.getFileNamesAndRedshifts(weakImgData)
if not fileNamesAndRedshifts:
    raise Exception("No weak lensing data found")

for fn, z_shear in fileNamesAndRedshifts:
    print(f"INFO: adding WL info for {fn} with z {z_shear}")
    iws.addImageDataToList(images.ImagesData.load(fn), z_shear, "sheardata", { "threshold": 0.1 })

#sheetType = "genome"
sheetType = "nosheet"
iws.setDefaultInversionArguments(sheetSearch = sheetType)

def setBasisFunctions(lens, minSub = None, maxSub = None):

    iws.clearBasisFunctions()
    # First the strong lensing grid
    iws.setUniformGrid(15) if not lens else iws.setSubdivisionGrid(lens, minSub, maxSub)
    iws.addBasisFunctionsBasedOnCurrentGrid()

    # Then, a uniform grid for weak lensing
    iws.setUniformGrid(48, regionSize = 30*ANGLE_ARCMIN)
    
    def myLensModelFunction(operation, operationInfo, parameters):
        r = inversion.defaultLensModelFunction(operation, operationInfo, parameters)
        if operation == "add":
            return (r[0], 0.1) # Set the mass that's counted in determining the scale factor in the GA to (almost) zero
        return r

    iws.addBasisFunctionsBasedOnCurrentGrid(myLensModelFunction, { "totalmass": 1e15*MASS_SUN })

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

