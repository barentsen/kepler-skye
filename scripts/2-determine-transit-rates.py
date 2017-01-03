"""Produces a table detailing the frequency of transits observed for
long-period/low-mes TCEs as a function of skygroup and season.
This table will help us identify suspicious peaks in the number of transits
at a given epoch/skygroup in the next step."""
import pandas as pd
from tqdm import tqdm

# The list of quarters in each season
SEASON_QUARTERS = {
                    0: [2, 6, 10, 14],
                    1: [3, 7, 11, 15],
                    2: [4, 8, 12, 16],
                    3: [1, 5, 9, 13, 17]
                  }


def count_days_per_season():
    """Returns a pandas dataframe specifying the number of days per season.

    The returned dataframe should look like this:

                n_days
    season
    0       369.398637
    1       373.424100
    2       325.282400
    3       347.636795
    """
    KEPLER_QUARTERS = pd.read_csv('../data/kepler-quarters.csv')[1:]
    DAYS_PER_QUARTER = KEPLER_QUARTERS.last_lc_mjd - KEPLER_QUARTERS.first_lc_mjd
    seasons = [0, 1, 2, 3]
    n_days = []
    for season in seasons:
        quarters = SEASON_QUARTERS[season]
        n_days.append(sum([DAYS_PER_QUARTER.ix[q] for q in quarters]))
    df_n_days = pd.DataFrame({'n_days': n_days}, index=seasons)
    df_n_days.index.name = 'season'
    return df_n_days


def get_transit_rates(transit_table_fn, planet_candidates_only=True):
    """Returns a pandas dataframe detailing the observed frequency of transits
    grouped by skygroup and season.
    """
    # Read the list of transits
    transits = pd.read_csv(transit_table_fn)
    if planet_candidates_only:
        planet_candidates_mask = (
                                    (transits.not_transit_like_flag == 0) &
                                    (transits.ephemeris_match_flag == 0)
                                  )
        transits = transits[planet_candidates_mask]

    # Group the transits by skygroup and season and count their number;
    # store the counts in a dataframe
    groupby = transits.groupby(['skygroup', 'season'])
    df_n_transits = pd.DataFrame(groupby.size(), columns=['n_transits'])
    df_n_transits['n_tces'] = groupby.tce.nunique()  # number of unique TCEs

    # To determine the number of transits per day for a given season,
    # we first need to count the number of days contained in each season.
    df_n_days = count_days_per_season()
    # How many transits per day are in each skygroup/season?
    df_rates = df_n_transits.merge(df_n_days, left_index=True, right_index=True)
    df_rates['transits_per_day'] = df_rates.n_transits / df_rates.n_days

    # Add the channel number
    channel = []
    for skygroup, season in df_rates.index:
        channel.append(transits[(transits.skygroup == skygroup) & (transits.season == season)].channel.iloc[0])
    df_rates['channel'] = channel

    return df_rates


if __name__ == '__main__':
    ops = get_transit_rates('intermediate-output/ops-tces-transits.csv', planet_candidates_only=True)
    ops.to_csv('intermediate-output/ops-tces-transit-rates.csv')

    inv = get_transit_rates('intermediate-output/inv-tces-transits.csv', planet_candidates_only=True)
    inv.to_csv('intermediate-output/inv-tces-transit-rates.csv')

    ss1 = get_transit_rates('intermediate-output/ss1-tces-transits.csv', planet_candidates_only=True)
    ss1.to_csv('intermediate-output/ss1-tces-transit-rates.csv')
