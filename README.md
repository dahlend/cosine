# Analysis for comets observed by WISE


Comet orbits are obtained by querying JPL Horizons for a set of all known cometary orbits
JPL fits comet orbits typically on each apparation of the orbit, so comets such as 1P have
*many* orbit fits, in some cases going back millennia. Many of these orbits, especially the
older ones, have poor observational data, and as a result have extremely poor orbit fits.

Steps in this analysis:
1) Query horizons for the known set of cometary orbits, the result of this is saved in the
   "jpl_comets.csv" file, which defines JPL's knowledge of all known comets, and the epochs
   of the orbit fits. This is considered the master list of all cometary names.
2) Down select this list to only contain orbit epoch fits which obey the following criterion:
      a) Epoch of orbit fit is after 1980
      b) Name does not include "D/" for destroyed comets
      c) The standard deviation of position knowledge is greater than 0.01 AU
      d) during the period of WISE, from cryo through end of 2023, the comet at some point
         must be closer than 11.5 AU from the sun. (See figure below for a histogram of this)
      e) explicitly exclude 5 objects observed by WISE but not on this list. 
   This down selection results in the "keep.txt" file of comet names.
3) Using the list of "good" names, download the SPICE kernels from Horizons for these objects.
   This step is very time consuming, it is recommended to ask Dar for a copy of these files.
   There is a single object, 2013 A1 (siding spring) which is not available via the JPL website,
   and had to be manually requested via email. This object had a *very* close encounter with Mars.
4) Run spice based visibility checks for every field of view for WISE over the mission (2 hours).
5) Load all minor planet center observations for comets submitted by WISE, dropping everything
   submitted after the end of 2023 data release. (total of 7125)
6) Match all MPC detection to a neospy predicted observation, restricted by time to less
   than 5 seconds. (There are some slight irregularities with how time was calculated in the
   WISE data)


