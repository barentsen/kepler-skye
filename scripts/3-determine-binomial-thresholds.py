"""This final step identifies the (epoch, skygroup) pairs that show a
'suspiciously large' number of transits from long-period/low-mes TCEs;
such clusters of transits are indicative of noisy data.

The term 'suspiciously large' is quantified by means of a probability
to see the given number of transits, as inferred from the appropriate
binomial distribution.
"""
import pandas as pd
import numpy as np
from scipy.stats import binom, poisson
from tqdm import tqdm


KEPLER_BEGIN_BK, KEPLER_END_BK = 130, 1582

# Exclude the channels below from the statistics; these are known to show
# a large spurious population of TCEs due to rolling band artefacts.
OUTLIER_CHANNELS = [26, 44, 58]


class SkyeMetric(object):
    """This class reads in transit data and identifies (epoch, skygroup) pairs
    that show a suspiciously large number of transits."""

    def __init__(self, transit_table_fn, rates_table_fn, binsize=1.,
                 probability_threshold=1e-3, probability_threshold_combined=None):
        if probability_threshold_combined is None:
            probability_threshold_combined = probability_threshold
        print('Reading {} and {}'.format(transit_table_fn, rates_table_fn))
        self.transits = pd.read_csv(transit_table_fn)
        self.rates = pd.read_csv(rates_table_fn, index_col=['skygroup', 'season'])
        self.binsize = binsize
        self.probability_threshold = probability_threshold
        self.probability_threshold_combined = probability_threshold_combined

        self.transits['bin_id'] = binsize * np.floor(self.transits.transit_time / binsize)
        self._compute()

    def _compute(self):
        self.thresholds = self._compute_binomial_thresholds()
        self.bad_bins = self._compute_bad_bins()
        self.bad_bin_ids_all_channels = self._compute_bad_bins_all_channels()

    def _compute_binomial_thresholds(self):
        rates = self.rates
        tce_expect_col, rate_expected_col, rate_threshold_col = [], [], []
        for skygroup, season in rates.index:
            mask_reference = ((rates.index.get_level_values('skygroup') == skygroup) &
                              (rates.index.get_level_values('season') != season) &
                              ~rates.channel.isin(OUTLIER_CHANNELS))

            # Compute the probability for a TCE to produce a transit in a given bin
            n_transits_per_bin = self.binsize * rates[mask_reference].transits_per_day
            n_tces = rates[mask_reference].n_tces
            mean_transit_probability = (n_transits_per_bin / n_tces).mean()

            rate_expected_col.append(n_transits_per_bin.median())
            tce_expect_col.append(n_tces.median())
            rate_threshold_col.append(int(binom.ppf(1 - self.probability_threshold,
                                                    int(n_tces.median()),
                                                    mean_transit_probability)))

        rates['n_tces_expected'] = tce_expect_col
        rates['transit_rate_expected'] = rate_expected_col
        rates['transit_rate_threshold'] = rate_threshold_col
        return rates

    def _compute_bad_bins(self):
        transit_rate_groupby = self.transits.groupby(['bin_id', 'skygroup', 'season'])
        transit_rate = pd.DataFrame(transit_rate_groupby.size(), columns=['transit_rate_observed'])
        transit_rate['channel'] = transit_rate_groupby.first().channel
        transit_rate['quarter'] = transit_rate_groupby.first().quarter
        assert((transit_rate_groupby.first().channel == transit_rate_groupby.mean().channel).all())
        assert((transit_rate_groupby.first().quarter == transit_rate_groupby.mean().quarter).all())

        threshold_col = []
        for bin_id, skygroup, season in transit_rate.index:
            threshold_col.append(self.thresholds.ix[skygroup, season].transit_rate_threshold)
        transit_rate['transit_rate_threshold'] = threshold_col

        # Identify the bad bins
        mask_bad_bins = transit_rate.transit_rate_observed > transit_rate['transit_rate_threshold']
        print('Identified {} bad bins for single channels.'.format(mask_bad_bins.sum()))

        mask_transits_in_bad_bins = self.transits.channel > 999  # All False
        for bin_id, skygroup, season in transit_rate[mask_bad_bins].index:
            mask_transits_in_bad_bins |= (
                                    (self.transits.bin_id == bin_id) &
                                    (self.transits.skygroup == skygroup) &
                                    (self.transits.season == season)
                                 )
        self.mask_transits_in_bad_bins = mask_transits_in_bad_bins
        print('Flagged {} out of {} transits as suspicious.'.format(mask_transits_in_bad_bins.sum(),
                                                                    len(mask_transits_in_bad_bins)))

        return transit_rate[mask_bad_bins]

    def _compute_bad_bins_all_channels(self):
        tces_to_remove = self.transits[self.mask_transits_in_bad_bins].tce.unique()
        mask_transits_to_remove = self.transits.tce.isin(tces_to_remove)

        # Make a cut across all channels
        total_transit_count = self.transits[~mask_transits_to_remove].groupby('bin_id').size()
        n_tces = self.transits[~mask_transits_to_remove].tce.unique().size
        p_transit = total_transit_count.median() / n_tces
        total_count_threshold = binom.ppf(1 - self.probability_threshold_combined,
                                          int(n_tces),
                                          p_transit)
        # print('n={} p={}'.format(n_tces, p_transit))

        bins_to_remove = total_transit_count[total_transit_count > total_count_threshold].index
        print('Identified {} bad bins for all channels.'.format(len(bins_to_remove)))

        self.mask_transits_in_bad_bin_ids = self.transits['bin_id'].isin(bins_to_remove)

        return bins_to_remove

    def percent_removed(self):
        """Returns the fraction of the total Kepler data set removed by the metric."""
        n_bins_total = len(self.transits['bin_id'].unique()) * 80  # assume 80 channels
        n_bins_removed = len(self.bad_bins)  + 80 * len(self.bad_bin_ids_all_channels)
        return 100 * n_bins_removed / n_bins_total

    def print_summary(self):
        print('Number of bad bins per channel:')
        print(self.bad_bins.reset_index().groupby('channel').size().sort_values(ascending=False).head(10))
        print('Metric removes {:.1f}% of data.'.format(self.percent_removed()))


    def write_definition(self, output_fn):
        print('Writing skye metric definition to {}'.format(output_fn))
        with open(output_fn, 'w') as out:
            out.write('# This file specifies the times (floored bkjd) and skygroups during which\n')
            out.write('# an anomalous number of long-period (>50 days) transits were detected.\n')
            out.write('#\n')
            out.write('# bkjd skygroup\n')
            for row in self.bad_bins.reset_index().itertuples():
                out.write('{:.2f} {}\n'.format(row.bin_id, row.skygroup))
            for bin_id in self.bad_bin_ids_all_channels:
                for channel in range(1, 85):
                    out.write('{} {}\n'.format(bin_id, channel))

    def plot(self, output_fn):
        mask_transits_flagged = self.mask_transits_in_bad_bins | self.mask_transits_in_bad_bin_ids
        tces_to_remove = self.transits[mask_transits_flagged].tce.unique()
        mask_transits_to_remove = self.transits.tce.isin(tces_to_remove)

        pc_before = (self.transits.groupby('tce').first().disposition == 'PC').sum()
        pc_after = (self.transits[~mask_transits_to_remove].groupby('tce').first().disposition == 'PC').sum()
        print('Retain {} out of {} PCs.'.format(pc_after, pc_before))

        q10, q50, q90 = np.percentile(self.transits[~mask_transits_to_remove].groupby('tce').first().mes, [10, 50, 90])
        print('mes removed [{}, {}, {}]'.format(q10, q50, q90))

        import matplotlib.pyplot as pl
        from matplotlib import rcParams
        rcParams["savefig.dpi"] = 100
        import seaborn as sns
        print('Writing {}'.format(output_fn))
        pl.figure(figsize=(10, 7))
        pl.subplot(211)
        _ = pl.hist(self.transits.groupby('tce').first().period,
                    bins=100, label='TCEs considered', facecolor='red')
        _ = pl.hist(self.transits[~mask_transits_to_remove].groupby('tce').first().period,
                    bins=100, label='Not affected by cut', facecolor='blue')
        pl.xlabel('Period [days]', fontsize=14)
        pl.legend()
        pl.title('Binomial threshold (binsize={}, p={})'.format(self.binsize, self.probability_threshold), fontsize=18)
        #pl.text(700, 400, 'Only showing\nperiod > 100 days\nmax_mult_ev < 20', ha='left')
        pl.text(600, 400, 'Removes {:.1f}% of data\nRetains {} out of {} transit-like TCEs\nMedian MES removed: {:.1f}'.format(self.percent_removed(), pc_after, pc_before, q50), ha='left')

        pl.subplot(212)
        _ = pl.hist(self.transits.transit_time,
                    lw=0, facecolor='red', label='All long-period/low-mes transits',
                    bins=(KEPLER_END_BK - KEPLER_BEGIN_BK),
                    range=(KEPLER_BEGIN_BK, KEPLER_END_BK))
        _ = pl.hist(self.transits[~mask_transits_flagged].transit_time,
                    lw=0, facecolor='blue', label='Transits deemed ok',
                    bins=(KEPLER_END_BK - KEPLER_BEGIN_BK),
                    range=(KEPLER_BEGIN_BK, KEPLER_END_BK))
        pl.ylim([0, 200])
        pl.xlim([KEPLER_BEGIN_BK, KEPLER_END_BK])
        pl.legend()
        pl.xlabel('Transit time', fontsize=14)
        #pl.text(1390, 120, 'Only showing\nperiod > 100 days\nmax_mult_ev < 20', ha='left')
        pl.tight_layout()
        pl.savefig(output_fn)
        pl.close()


def run_skye(prefix, binsize=0.5, p_threshold=1e-6, p_threshold_global=1e-6):
    output_prefix = 'output/{}-bin{:.2f}-p{:.0e}'.format(prefix, binsize, p_threshold)
    transit_table_fn = 'intermediate-output/' + prefix + '-tces-transits.csv'
    transit_rates_fn = 'intermediate-output/' + prefix + '-tces-transit-rates.csv'
    skye = SkyeMetric(transit_table_fn,
                      transit_rates_fn,
                      binsize=binsize,
                      probability_threshold=p_threshold,
                      probability_threshold_combined=p_threshold_global)
    skye.print_summary()
    skye.write_definition(output_prefix + '-definition.txt')
    skye.plot(output_prefix + '-skye.png')
    skye.thresholds.to_excel(output_prefix + '-thresholds.xls')
    skye.bad_bins.to_excel(output_prefix + '-badbins.xls')
    return skye


if __name__ == '__main__':

    for prefix in ['ops', 'inv', 'ss1']:
        for binsize in [0.5]:
            for p_threshold in [1e-4, 5e-4]:
                skye = run_skye(prefix,
                                binsize=binsize,
                                p_threshold=p_threshold,
                                p_threshold_global=1e-12)
