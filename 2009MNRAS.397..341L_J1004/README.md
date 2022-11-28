
[results.ipynb](results.ipynb) shows the results based on the models in the
article, and [resultstd.ipynb](resultstd.ipynb) shows the time delays
based on the second model.

Directory [inversion_grale1](inversion_grale1) shows the inversion script
and input that was used with grale1/graleshell. A recreation using grale2/pygrale
can be found in [inversion_grale2](inversion_grale2). In the inversion script
it is indicated which settings you can easily vary this way, but the current
settings are a translation of the grale1 ones. Results of this grale2 inversion,
both with the defaults from the script and a few other settings, can be
viewed in [results-with-grale2.ipynb](results-with-grale2.ipynb). The github
viewer won't show the 3D interactive maps, but when
[opened with nbviewer](https://nbviewer.org/github/j0r1/LensModels/blob/master/2009MNRAS.397..341L_J1004/results-with-grale2.ipynb)
they should be visible.

Note that in grale1, the inversion grids that are used always have the same
center. In the grale2 version a small random offset is applied each time, to
increase the variation in results (this can be turned off again, if needed for
some reason).

