This directory contains a copy of the output of the Kepler planet detection
pipeline, version DR25, copied on Oct 25, 2017.

The pipeline was run on different sets of lightcurves to assess completeness
and reliability. These runs include actual Kepler lightcurves ('OPS'),
inverted lightcurves ('INV'), and season-scrambled lightcurves ('SS1').

The TCEs detected by these different runs are located in this directory
under the following names:

* OPS-TCEs.txt (copied from /soc/nfs/so-nfs/DR25/OPS/DATA/TCEs.txt)
* INV-TCEs.txt (copied from /soc/nfs/so-nfs/DR25/INV/DATA/TCEs.txt)
* SS1-TCEs.txt (copied from /soc/nfs/so-nfs/DR25/SS1/DATA/TCEs.txt)

The table called `kepler-ccd-info-by-kic.csv` was created using the
following query at MAST/CasJobs:

    SELECT
        sci_kepler_id, sci_data_quarter, sci_season, sci_skygroup_id,
        sci_module, sci_channel, sci_output
    INTO mydb.kepler_ccd_info
    FROM kepler.kepler_science
    WHERE sci_archive_class = 'TPL';

Finally, the table `kepler-quarters.csv` details the quarter numbers
and time intervals for all Kepler Quarters. 
It was created by hand based on the data release notes.
