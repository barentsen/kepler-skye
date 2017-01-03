"""
Takes the sample of TCEs detected in the INV (inversion) and
SS1 (season scrambling) runs and add the sample of 'reliable' TCEs
from the OPS (i.e. real data) run to them.  The purpose of doing this is
to obtain a TCE sample to which the Skye false-positive identification
algorithm can meaningfully be applied.
"""

import pandas as pd

TCE_TABLE_COLUMNS = ['tce', 'kic', 'pn', 'period', 'epoch', 'mes', 'depth', 'duration',
                     'rplanet', 'rstar', 'tstar', 'logg', 'a', 'radratio', 'arstar', 'snr', 'srad']
ROBOVETTER_COLUMNS = ['tce', 'score', 'disposition', 'not_transit_like_flag',
                      'significant_secondary_flag', 'centroid_offset_flag',
                      'ephemeris_match_flag', 'minor_descriptive_flag']


def xmatch_robovetter(tce_input_fn, robovetter_input_fn):
    """Adds the columns from the robovetter output to the TCE list."""
    tce_df = pd.read_fwf(tce_input_fn, comment='#', names=TCE_TABLE_COLUMNS)
    robovetter_df = pd.read_fwf(robovetter_input_fn, comment='#', names=ROBOVETTER_COLUMNS)
    xmatch_df = pd.merge(tce_df, robovetter_df, on='tce', how='inner')
    # Both tcedf_tmp and robovetterdf should have the same number of TCEs
    assert(len(xmatch_df) == len(robovetter_df))  # sanity check
    return xmatch_df


def xmatch_clean(tce_df, clean_input_fn):
    """Only use the TCEs which are in the 'clean' table."""
    clean_df = pd.read_csv(clean_input_fn)
    xmatch_df = pd.merge(tce_df, clean_df, on='tce', how='inner')
    # Sanity check: are SS1-TCEs.txt and the clean table aligned?
    #assert(len(df_clean) == len(tce_df))
    #assert((df_clean.tce == tce_df.tce).sum() == len(df_clean))
    # Make a new dataframe containing only the TCEs to keep
    return xmatch_df[xmatch_df.keep]


def add_real_planet_candidates(tce_df):
    """Take a sample of TCEs and add reliable OPS TCEs to them."""
    ops_df = xmatch_robovetter('../data/OPS-TCEs.txt', '../data/RoboVetterOut-OPS.txt')
    planet_candidates_mask = (
                                (ops_df.not_transit_like_flag == 0) &
                                (ops_df.ephemeris_match_flag == 0)
                              )
    return pd.concat((tce_df, ops_df[planet_candidates_mask]))


if __name__ == '__main__':
    df = xmatch_robovetter('../data/SS1-TCEs.txt', '../data/RoboVetterOut-SS1.txt')
    df = xmatch_clean(df, '../data/ss1TCEClean-900day-7mes-Dec2016.csv')
    df = add_real_planet_candidates(df)
    df.to_csv('intermediate-output/ss1-tces-for-skye.csv')

    df = xmatch_robovetter('../data/INV-TCEs.txt', '../data/RoboVetterOut-INV.txt')
    df = xmatch_clean(df, '../data/invTCEClean-100day-9mes-Nov2016-no-header.csv')
    df = add_real_planet_candidates(df)
    df.to_csv('intermediate-output/inv-tces-for-skye.csv')
