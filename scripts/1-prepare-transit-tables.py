"""Creates a table of long-period/low-mes TCE transits detected by the Kepler pipeline.

This script takes a 'TCEs.txt' table created by the Kepler pipeline and turns
it into a table of transits, detailing the time, quarter, season, ccd channel,
skygroup, etc.  This table is intended to allow the frequency of transits as a
function of time and CCD to be investigated, allowing time intervals and detectors
to be identified that produces a spurious number of transits at a given time,
in particular due to thermal changes in the spacecraft electronics.
"""
import numpy as np
import pandas as pd
from tqdm import tqdm

# A few useful constants
KEPLER_BEGIN_BK, KEPLER_END_BK = 130, 1582
KEPLER_QUARTERS = pd.read_csv('../data/kepler-quarters.csv')


def mjd2quarter(mjd):
    """Returns the Kepler quarter for a given Modified Julian Day."""
    mask = (KEPLER_QUARTERS.first_lc_mjd < mjd+0.01) & (KEPLER_QUARTERS.last_lc_mjd > mjd-0.01)
    if mask.any():
        return KEPLER_QUARTERS.loc[mask, 'quarter'].values[0]
    return None

def mjd2season(mjd):
    """Returns the Kepler season for a given Modified Julian Day."""
    mask = (KEPLER_QUARTERS.first_lc_mjd < mjd+0.01) & (KEPLER_QUARTERS.last_lc_mjd > mjd-0.01)
    if mask.any():
        return KEPLER_QUARTERS.loc[mask, 'season'].values[0]
    return None

def bkjd_to_mjd_approximate(bkjd):
    """Inexact conversion from Barycentric Kepler Julian Date (BKJD) to Modified Julian Date (MJD).

    Inexact because it ignores the TIMECORR and TIMSLICE corrections.
    """
    return bkjd + 2454833 - 2400000.5


def make_transit_table(tce_input_fn, robovetter_input_fn=None,
                       min_period=50, max_mes=9e99):
    # Read the TCE table
    columns = ['tce', 'kic', 'pn', 'period', 'epoch', 'mes', 'depth', 'duration',
               'rplanet', 'rstar', 'tstar', 'logg', 'a', 'radratio', 'arstar', 'snr', 'srad']
    tcedf_tmp = pd.read_fwf(tce_input_fn, comment='#', names=columns)

    if robovetter_input_fn is None:
        tcedf = tcedf_tmp
        tcedf['disposition'] = [None] * len(tcedf)
    else:
        # Add preliminary robovetter output to have best-effort dispositions
        robovetter_columns = ['tce', 'score', 'disposition', 'not_transit_like_flag',
                              'significant_secondary_flag', 'centroid_offset_flag',
                              'ephemeris_match_flag', 'minor_descriptive_flag']
        robovetterdf = pd.read_fwf(robovetter_input_fn, comment='#', names=robovetter_columns)

        # Both tcedf_tmp and robovetterdf should have the same number of TCEs
        tcedf = pd.merge(tcedf_tmp, robovetterdf, on='tce')
        # Sanity checks
        assert(len(tcedf_tmp) == len(robovetterdf))
        assert(len(tcedf) == len(tcedf_tmp))

    # Convert the catalog of TCEs into a catalog of LONG-PERIOD TRANSITS
    mask = (tcedf.period > min_period) & (tcedf.mes < max_mes)
    print('Selected {} out of {} TCEs using period and mes cut.'.format(mask.sum(), len(tcedf)))
    transitrows = []
    for mytce in tqdm(tcedf[mask].itertuples(), desc='Identifying transits from {}'.format(tce_input_fn)):
        mytime = mytce.epoch
        while mytime < KEPLER_END_BK:
            mjd = bkjd_to_mjd_approximate(mytime)
            newrow = {'transit_time': mytime,
                      'transit_time_mjd': mjd,
                      'tce': mytce.tce,
                      'kic': mytce.kic,
                      'period': mytce.period,
                      'mes': mytce.mes,
                      'disposition': mytce.disposition,
                      'not_transit_like_flag': mytce.not_transit_like_flag,
                      'significant_secondary_flag': mytce.significant_secondary_flag,
                      'centroid_offset_flag': mytce.centroid_offset_flag,
                      'ephemeris_match_flag': mytce.ephemeris_match_flag,
                      'minor_descriptive_flag': mytce.minor_descriptive_flag,
                      'quarter': mjd2quarter(mjd),
                      'season': mjd2season(mjd)}
            if newrow['quarter'] is not None:  # Ignore unobserved transits (falling between quarter boundaries)
                transitrows.append(newrow)
            mytime += mytce.period
    transits = pd.DataFrame(transitrows)
    print('Found {} long-period transits in {}'.format(len(transits), tce_input_fn))

    # Add the CCD channel number and mod/out for each transit
    ccdinfo = pd.read_csv('../data/kepler-ccd-info-by-kic.csv')
    # This join will remove unobserved transits, i.e. transits in a quarter on a dead module
    transits_with_ccdinfo = pd.merge(transits, ccdinfo,
                                     left_on=['kic', 'quarter'],
                                     right_on=['sci_kepler_id', 'sci_data_quarter'])
    assert((transits_with_ccdinfo.season == transits_with_ccdinfo.sci_season).all())  # Sanity check

    print('Found {} observed long-period transits'.format(len(transits_with_ccdinfo)))
    transits_with_ccdinfo.rename(columns={'sci_skygroup_id': 'skygroup',
                                          'sci_channel': 'channel',
                                          'sci_module': 'module',
                                          'sci_output': 'output'},
                                 inplace=True)
    export_columns = ['tce', 'kic', 'transit_time', 'period', 'mes', 'disposition',
                      'not_transit_like_flag', 'significant_secondary_flag', 'centroid_offset_flag',
                      'ephemeris_match_flag', 'minor_descriptive_flag',
                      'quarter', 'season', 'skygroup', 'channel', 'module', 'output']
    return transits_with_ccdinfo[export_columns]


if __name__ == '__main__':
    # Real 'OPS' TCE detection results
    ops = make_transit_table('../data/OPS-TCEs.txt', robovetter_input_fn='../data/RoboVetterOut-OPS.txt')
    ops.to_csv('ops-tces-transits.csv', index=False)

    # Inversion run
    #inv = make_transit_table('../data/INV-TCEs.txt', robovetter_input_fn=None)
    #inv.to_csv('inv-tces-transits.csv', index=False)

    # Season scrambling run
    #ss1 = make_transit_table('../data/SS1-TCEs.txt', robovetter_input_fn=None)
    #ss1.to_csv('ss1-tces-transits.csv', index=False)
