from grale.all import *

tdType = "Paper2009" # Set to "NoSrc" to use the newer time delay fitness from 2020MNRAS.494.3253L
maxSub = 5 # Produces inv1, inv2, ... inv5, like in the grale1 version (you can go higher as well)
corrMassScale = 0.5e13*MASS_SUN # Mass scale for the corrections, this is what was used in 2009
pressure = 2.0 # The default selection pressure parameter is 2.5, in the 2009 work this was set
               # slightly lower, to try to avoid getting stuck in a local optimum

# Write the RNG state, in case we want to reproduce the run exactly
# (note that the GRALE_DEBUG_SEED environment variable will need
# to be restored as well)
import random
print("RNG State:")
print(random.getstate())

z_lens = 0.68

z_Q = 1.734
z_A = 3.332
z_B = 2.74

invCenter = V(0, 5)*ANGLE_ARCSEC
slWidth = 35*ANGLE_ARCSEC

cosm = cosmology.Cosmology(0.71, 0.27, 0, 0.73)
cosmology.setDefaultCosmology(cosm)

# As in the grale1 version, the inversion consists of two steps: first a base
# lens is sought, then small scale corrections are applied. In the grale1
# version, manual inspection of the fitness values would reveal if we got
# stuck in a local optimum - here we can detect this in the script, and
# simply restart the whole inversion if this is the case.
# 
# Ideally of course, the GA would be better and avoid getting stuck in such
# a local optimum.
#
# This function does the refinement steps inv1, inv2, ... and is used below
def findBaseLens(attempt):

    iws = inversion.InversionWorkSpace(z_lens, slWidth, invCenter)
    iws.addImageDataToList(images.ImagesData.load("source1points.imgdata"), z_Q, "extendedimages", { "userectangles": False }) # Only use the actual points
    iws.addImageDataToList(images.ImagesData.load("source1null.imgdata"), z_Q, "extendednullgrid")
    iws.addImageDataToList(images.ImagesData.load("source1critgrid.imgdata"), z_Q, "causticgrid")

    iws.addImageDataToList(images.ImagesData.load("source2points.imgdata"), z_A, "extendedimages")
    iws.addImageDataToList(images.ImagesData.load("source2null.imgdata"), z_A, "extendednullgrid")
    iws.addImageDataToList(images.ImagesData.load("source2critgrid.imgdata"), z_A, "causticgrid")

    iws.addImageDataToList(images.ImagesData.load("source3points.imgdata"), z_B, "extendedimages")
    iws.addImageDataToList(images.ImagesData.load("source3null.imgdata"), z_B, "extendednullgrid")
    iws.addImageDataToList(images.ImagesData.load("source3critgrid.imgdata"), z_B, "causticgrid")


    iws.setDefaultInversionArguments(sheetSearch = "genome",
        geneticAlgorithmParameters = { "selectionpressure": pressure },
        fitnessObjectParameters = { "fitness_timedelay_type": tdType })

    prevLens, subDivStart = None, 100
    bestLens, bestFitness = None, None

    for idx in range(1, maxSub+1):

        iws.setUniformGrid(17) if prevLens is None else iws.setSubdivisionGrid(prevLens, subDivStart, subDivStart+100)

        lens, fitness, fitdesc = iws.invert(1024)
        lens.save(f"inv{idx}-attempt{attempt}.lensdata")

        if bestLens is None or fitness < bestFitness:
            bestLens = lens
            bestFitness = fitness

        prevLens = lens
        subDivStart += 200

    print(f"BESTFITNESS-attempt{attempt}: {bestFitness}")
    bestLens.save(f"best-{attempt}.lensdata")

    return bestLens, bestFitness

# Look for base lens

maxAttempts = 5
for i in range(maxAttempts):
    bestLens, bestFitness = findBaseLens(i)
    caustFitness, nullFitness, overlapFitness, tdFitness = bestFitness
    if caustFitness == 0 and nullFitness == 0:
        print("Found correct optimum, continuing")
        break
    else:
        print("Seems like local optimum, retrying")

else:
    raise Exception(f"Couldn't find a base lens in {maxAttempts} tries")

print(f"BESTFITNESS-baselens: {bestFitness}")
bestLens.save("baselens.lensdata")

# Add small scale corrections to improve this, with more detailed grids/images

iws = inversion.InversionWorkSpace(z_lens, slWidth, invCenter)
iws.addImageDataToList(images.ImagesData.load("source1points.imgdata"), z_Q, "extendedimages", { "userectangles": False }) # Only use the actual points
iws.addImageDataToList(images.ImagesData.load("source1detailednull.imgdata"), z_Q, "extendednullgrid")
iws.addImageDataToList(images.ImagesData.load("source1detailedcritgrid.imgdata"), z_Q, "causticgrid")

iws.addImageDataToList(images.ImagesData.load("source2detailedpoints.imgdata"), z_A, "extendedimages", { "userectangles": False })
iws.addImageDataToList(images.ImagesData.load("source2detailednull.imgdata"), z_A, "extendednullgrid")
iws.addImageDataToList(images.ImagesData.load("source2detailedcritgrid.imgdata"), z_A, "causticgrid")

iws.addImageDataToList(images.ImagesData.load("source3detailedpoints.imgdata"), z_B, "extendedimages")
iws.addImageDataToList(images.ImagesData.load("source3detailednull.imgdata"), z_B, "extendednullgrid")
iws.addImageDataToList(images.ImagesData.load("source3detailedcritgrid.imgdata"), z_B, "causticgrid")

iws.setDefaultInversionArguments(sheetSearch = "nosheet", # Don't look for a mass sheet when applying corrections! 
        baseLens = bestLens,
        allowNegativeValues = True,
        massScale = corrMassScale,
        fitnessObjectParameters = { "fitness_timedelay_type": tdType },
    )

iws.setUniformGrid(64)

corrections, fitness, fitdesc = iws.invert(1024)

correctedLens = lenses.CompositeLens(corrections.getLensDistance(), [
        { "x": 0, "y": 0, "factor": 1, "angle": 0, "lens": corrections },
        { "x": 0, "y": 0, "factor": 1, "angle": 0, "lens": bestLens },
    ])

print(f"CORRECTEDFITNESS: {fitness}")
correctedLens.save("correctedlens.lensdata")

