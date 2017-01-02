# Kepler's Skye Planet Candidate Metric

***Identifies suspicious clusters of transits in time and space.***

## Introduction

Some of the CCD's on board the planet-hunting Kepler spacecraft
are known to be sensitive to thermal changes in the on-board electronics;
producing 'rolling band' artefacts which may introduce transit-like
signals in the lightcurves of stars.

An effective way to identify the transits that are likely due to this
type of artefact, is to look for clusters of transits at the same time
on the same CCD chip. This is what the so-called 'Skye metric' does
as part of the Kepler DR25 planet candidate vetting effort.

This repository contains the scripts used to identify the (time, ccd) pairs
during which the Kepler pipeline detected a suspicious number of transits,
where 'suspicious' is quantified by means of the inferred binomial distribution
of the observed frequency of transits produced by reliable planet candidates.


## Usage

The metric is implemented as a 4-step procedure, i.e. as 4 Python scripts,
located in the `scripts` folder of this repository.  The folder contains a
`Makefile` which captured how the scripts were ran.


## Output

The key output is the file called `output/ops-bin0.50-p1e-04-definition.txt`
which looks like this:

```
# This file specifies the times (floored bkjd) and skygroups during which
# an anomalous number of long-period (>50 days) transits were detected.
#
# bkjd skygroup
131.50 36
131.50 59
132.50 56
133.50 41
133.50 43
133.50 67
134.00 43
134.00 78
135.00 54
etc...
```

This file was used as input into the Kepler DR25 RoboVetter.
