# Jacobs2022-eLife
Code for Jacobs et al., 2022 Research Advance in eLife

Note the source data for figures can be found through the Dryad repository - doi:10.5061/dryad.9s4mw6mkn. This repo is still in the process of being set up by eLife as of 9/26/22.

## Required packages


For python (v3+) scripts - many of these will be included in your base environment if you use anaconda:
  - numpy, pandas, os, csv, glob, re, math, scipy, mapplotlib, pylab, datetime

For the one R code - 
  - dplyr, tidyr, reshape2


## Descriptions for each code:

**AUCautomation** - used to get area under the curve values around a particular window using the subject averages. Takes in a file where each column is a subject and each row is a sample for that subject over some epoch of interest. 

**CCF_processor**- used to get  peak (either min peak ,max peak, or absolute peak), n, and SD for an average cross correlation function. Note this is not used until data are already processed using the **Runcrosscorr_alltrials** and **CrossCorr_getavgs codes**.

**CrossCorr_getavgs** - the only R code which is a quick way to get the average of all cross correlation function over all trials. This could be done in Python obviously. But at the time of writing I was working more in R. Takes a csv where each column is a time-lag value in a cross correlation function and each row is a trial. 

**FR1_NPM_gob** and **JPTask_NPM_gob** - these codes do much of the heavy lifting in terms of analyzing fiber photometry data with some slight differences in things because the FR1 data has to be parsed differently than the PRT data (e.g. shock trials are only in the PRT). These take in both timestamps from bonsai and raw bonsai pixel value readings and parses the data based on each trial. Will also do all the rescaling, delta-F calculations and peri event z scores to get an average for each subjects for each session and block. This is able to parse all this by using the csv's returned from the **BehavioralProcessing** code that contains the metadata about each trial. 

**Runcrosscorr_alltrials** - does the heavy lifting regarding cross correlation analysis. Takes a file which contains peri-event z scores for ALL trials (this is returned from **JPTask_NPM_gob**. Takes each trial, gets full covariance values, normalizes those values in between -1 and 1, and returns them as a csv where each row is a trial and each column is a time point in the cross correlation function (excluding the initial columns which provide trial descriptors (block, session, subject, etc.)

**corrs_ephys** - short code that takes in the averaged single unit and averaged photometry data for a given epoch and normalizes them together by bringing them to the same sample rate (20 Hz) and min-max normalizes them. Will also get the correlations between unit data and photometry data as well as unit data and shuffled photometry data.

**permtesting** - Runs permutation tests and consecutive threshold analysis on subject averaged data. Generally, this needs to be tweaked depending on the comparison being run (e.g. session comparison or risk block comparison). Will return a csv of p values where each row is a p value at a given time in the epoch and each column is the group the control was compared to. Will also return a csv where each a row can be a Nan or a 1 and each column is the corresponding comparison. If a string of 1's is seen it means the consecutive threshold was reached at that time point (default threshold is 14 due to the low pass filter of 3Hz - see paper), this isn’t really useful except for if one wanted to use it for plotting a bar over 'significant' time periods.


## Questions: contact jacobsd@ohsu.edu

Elife paper: In press - https://elifesciences.org/articles/78912

Earlier BioArchive paper: https://www.biorxiv.org/content/10.1101/2022.03.29.486234v1
