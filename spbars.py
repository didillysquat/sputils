#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from datetime import datetime
import re
import matplotlib.gridspec as gridspec
import spcolors
import itertools

class SPBars:
    """
    A utility script for creating stacked bar plots from SymPortal output documents.
    It will produce sequence only, ITS2 type profile only, or a combination type plot.
    It will return a tuple of a matplotlib.pyplot.figure and matplotlib.axes.Axes objects.
    It will have the option to plot the output or not, i.e. to be used as a script.
    It should have a few basic options such as orientation.

    Parameters:

        seq_count_table_path (str): The full path to the 'absolute.abund_and_meta' post-MED sequence count table
        to use for plotting. N.B. Plotting of the pre-MED count table is currently not supported.

        profile_count_table_path (str): The full path to the absolute 'abund_and_meta' ITS2 type profile count table
        to use for plotting.

        plot_type (str): A string denoting what type of plot to produce. 'seq_only, 'profile_only', 'seq_and_profile'.
        ['seq_only']

        ax (matplotlib.axes.Axes): if passed, the plot will be done on the provided ax object.
        If None, a fig and ax will be generated and returned. [None]

        figsize (tuple): A tuple that will be passed to the plt.figure figure creation. Units should be in mm.
        If not provided this will automatically deduced. [None]

        sample_uids_included (list<int>): A list of the sample uids that should be plotted.
        The samples will be plotted in this order.[None]

        sample_uids_excluded (list<int>): A list of the sample uids that should be excluded from plotting. [None]

        sample_name_compiled_re_included (compiled regex): A regular expression, which if matched by a sample name,
        will denote that the sample should be included in plotting. [None]

        sample_name_compiled_re_included (compiled regex): A regular expression, which if matched by a sample name,
        will denote that the sample should be included in plotting. [None]

        orientation (str): 'vertical' or 'horizontal'. 'v' or 'h' may be passed. ['h']

        legend (bool): whether to plot a legend

        relative_abundnce (bool): whether to plot the abundances as relative or absolute [True]

        num_seq_leg_cols (int): The number of sequences to plot in the legend. Sequences will be plotted in order
        of abundance. [20]

        num_profile_leg_cols (int): The number of profiles to plot in the legend. Sequences will be plotted in order
        of abundance. [20]

    Returns:
        tuple(matplotlib.pyplot.figure, matplotlib.axes.Axes): The figure and axes object that contain the plot

    """
    def __init__(
            self, seq_count_table_path=None, profile_count_table_path=None, plot_type='seq', figsize=None,
            sample_uids_included=None, sample_uids_excluded=None,
            sample_name_compiled_re_included=None, sample_name_compiled_re_excluded=None,
            orientation='h', legend=True,
            relative_abundnce=True, num_seq_leg_cols=20, num_profile_leg_cols=20):

        # Check the inputs
        self._check_path_exists([seq_count_table_path, profile_count_table_path])
        if plot_type not in ['seq_only', 'profile_only', 'seq_and_profile']:
            raise RuntimeError("plot_type must be one of 'seq_only', 'profile_only', or 'seq_and_profile'")
        if orientation not in ['h', 'v', 'horizontal', 'vertical']:
            raise RuntimeError("orientation must be one of 'h', 'v', 'horizontal', or 'vertical' ")
        if not type(legend) == bool:
            raise RuntimeError("legend should be a bool")
        if (
                (sample_uids_included is not None or sample_uids_excluded is not None) and
                (sample_name_compiled_re_included is not None or sample_name_compiled_re_excluded is not None)
        ):
            raise RuntimeError("please provide either uids OR sample_names to include/exclude. Not both.")
        if sample_uids_included is not None and sample_uids_excluded is not None:
            raise RuntimeError(
                "Please provide only one of sample_uids_included or sample_uids_excluded"
            )
        if sample_name_compiled_re_included is not None and sample_name_compiled_re_excluded is not None:
            raise RuntimeError(
                "Please provide only one of sample_name_compiled_re_included or sample_name_compiled_re_excluded"
            )
        self._check_valid_params(plot_type, profile_count_table_path, seq_count_table_path)

        if seq_count_table_path:

            (
                self.sample_name_to_sample_uid_dict,
                self.sample_uid_to_sample_name_dict,
                self.seq_count_df
            ) = self._make_seq_count_df(
                seq_count_table_path, relative_abundnce,
                sample_uids_included, sample_uids_excluded,
                sample_name_compiled_re_included, sample_name_compiled_re_excluded
            )

        # Figure setup
        self._setup_fig_and_ax(figsize, legend, orientation, plot_type)

        # Plotting colors setup
        self.grey_iterator = itertools.cycle(spcolors.greys)
        if plot_type in ['seq_only', 'seq_and_profile']:
            self.colour_hash_iterator = iter(spcolors.colour_list)
            self.pre_def_seq_colour_dict = spcolors.pre_def_color_dict
            self.seq_color_dict = self._make_seq_colour_dict()

        # TODO setup colors for profiles
        foo = 'bar'

    def _make_seq_colour_dict(self):
        seq_color_dict = {}
        for seq_name in list(self.seq_count_df):
            if seq_name in self.pre_def_seq_colour_dict:
                seq_color_dict[seq_name] = self.pre_def_seq_colour_dict[seq_name]
            else:
                try:
                    seq_color_dict[seq_name] = next(self.colour_hash_iterator)
                except StopIteration:
                    seq_color_dict[seq_name] = next(self.grey_iterator)
        return seq_color_dict

    def _setup_fig_and_ax(self, figsize, legend, orientation, plot_type):
        """
        Setup up the ax objects that the plot and legends will be plotted to.

        :param figsize: user supplied figure size tuple(<int>, <int>) in mm
        :param legend: bool whether to plot a legend
        :param orientation: orientation of the plot either vertical or horizontal
        :param plot_type: 'seq_only', 'profile_only', 'seq_and_profile'

        :return: None. But, self.bar_ax, self.leg_ax_one and self.leg_ax_two (if plot_type is 'seq_and_profile')
        will be set.
        """
        if legend:
            # Then we need to produce axes for the legends
            if plot_type == 'seq_and_profile':
                # Then we need to have two legend axes
                self._setup_seq_and_profile(figsize, orientation)
            else:
                # Then we are working with a single legend
                self._setup_seq_or_profile_only_plot(figsize, orientation)

    def _setup_seq_or_profile_only_plot(self, figsize, orientation):
        if orientation in ['v', 'vertical']:
            # Then we want the bar plot to sit next to the legend plots
            # The bars will span the full height of the figure
            if figsize:
                # TODO enforce a minimum size
                fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be one column so 89mm
                # height will dependent on sample number
                # We will need a min height for the one legend of 30
                # and 5mm for each sample
                height = max(30, (5 * len(self.seq_count_df.index)))
                width = 89
                fig = plt.figure(figsize=self._mm2inch((height, width)))
            gs = gridspec.GridSpec(1, 2)
            self.bar_ax = plt.subplot(gs[:, :1])
            self.leg_ax_one = plt.subplot(gs[:, 1:2])
        else:
            # Then we want the bar plot to sit above the legend plot
            if figsize:
                # TODO enforce a minimum size
                fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be two columns so 183
                # height will be fixed, 30 for legend 50 for bars
                height = 80
                width = 183
                fig = plt.figure(figsize=self._mm2inch((height, width)))
                # TODO implement multiple rows of barplots
            gs = gridspec.GridSpec(2, 1)
            self.bar_ax = plt.subplot(gs[:1, :])
            self.leg_ax_one = plt.subplot(gs[1:2, :])

    def _setup_seq_and_profile(self, figsize, orientation):
        if orientation in ['v', 'vertical']:
            # Then we want the bar plot to sit next to the legend plots
            # The bars will span the full height of the figure
            # The seq legend will span the top half
            # The profile legend will span the bottom half
            if figsize:
                # TODO enforce a minimum size
                fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be one column so 89mm
                # height will dependent on sample number
                # We will need a min height for the legends of 30 * 2
                # and 5mm for each sample
                height = max(60, (5 * len(self.seq_count_df.index)))
                width = 89
                fig = plt.figure(figsize=self._mm2inch((height, width)))
            gs = gridspec.GridSpec(2, 2)
            self.bar_ax = plt.subplot(gs[:2, :1])
            self.leg_ax_one = plt.subplot(gs[:1, 1:2])
            self.leg_ax_two = plt.subplot(gs[1:2, 1:2])

        else:
            # Then we want the bar plot to sit above the legend plots
            # The bars will span the full width of the figure
            # The seq legend will span the left half
            # The profile legend will span the right half
            if figsize:
                # TODO enforce a minimum size
                fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be two columns so 183
                # height will be fixed, 30 for legend 50 for bars
                height = 80
                width = 183
                fig = plt.figure(figsize=self._mm2inch((height, width)))
                # TODO implement multiple rows of barplots
            gs = gridspec.GridSpec(2, 2)
            self.bar_ax = plt.subplot(gs[:1, :2])
            self.leg_ax_one = plt.subplot(gs[:2, :1])
            self.leg_ax_two = plt.subplot(gs[:2, 1:2])

    def _check_path_exists(self, paths):
        for path in paths:
            if path is not None:
                if not os.path.exists(path):
                    raise FileNotFoundError(f'Cannot find {path}')

    def _mm2inch(self, *tupl):
        inch = 25.4
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    def _make_seq_count_df(
            self, seq_count_table_path, relative_abundnce,
                sample_uids_included, sample_uids_excluded,
                sample_name_compiled_re_included, sample_name_compiled_re_excluded
    ):
        """
        make the seq_count_df, and produce dicts linking sample_uid to sample_name in both directions
        seq_count_df will be only the absolute abundance, with sequence_name as columns and sample_uid as index
        seq_count_df columns will be ordered by abundance across all samples
        samples will be excluded according to the user specified arguments
        """
        seq_count_df = pd.read_csv(seq_count_table_path, sep='\t', skipfooter=1, engine='python')
        sample_uid_to_sample_name_dict = {
            uid: s_name for uid, s_name in
            zip(seq_count_df['sample_uid'].values, seq_count_df['sample_name'].values)
        }
        sample_name_to_sample_uid_dict = {
            s_name: uid for uid, s_name in
            zip(seq_count_df['sample_uid'].values, seq_count_df['sample_name'].values)
        }
        seq_count_df.set_index(keys='sample_uid', drop=True, inplace=True)
        first_seq_index = self._find_first_seq_position(seq_count_df)
        seq_count_df = seq_count_df.iloc[:, first_seq_index:]

        seq_count_df = self._exclude_samples_from_seq_count_df(
            seq_count_df, sample_uids_included, sample_uids_excluded,
            sample_name_compiled_re_included, sample_name_compiled_re_excluded
        )

        # Convert to relative abundances
        if relative_abundnce:
            seq_count_df = seq_count_df.div(
                seq_count_df.sum(axis=1), axis=0
            )

        # sort by seq abundance, highest first
        sorted_seq_index = seq_count_df.sum(axis=0).sort_values(ascending=False).index
        seq_count_df = seq_count_df.reindex(sorted_seq_index, axis=1)

        return (sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict, seq_count_df)

    def _exclude_samples_from_seq_count_df(
            self, seq_count_df, sample_uids_included, sample_uids_excluded,
            sample_name_compiled_re_included, sample_name_compiled_re_excluded
    ):
        """
        Due to the checks we have already performed on the include/exclude arguments
        only one of them will be not None. This is the set we will use to exclude samples from the df.
        If all of them are none, then there are no samples to exclude.
        """
        if sample_uids_included is not None:
            # Return df containing only the included samples sorted according to that sampple order
            print('Excluding samples according to user supplied sample_uids_included list')
            # Check that all uids are found in the seq_count_df
            if not set(sample_uids_included).issubset(set(seq_count_df.index.values)):
                diff_set = set(sample_uids_included).difference(set(seq_count_df.index.values))
                raise RuntimeError('The following uids that were specified in the sample_uids_included list are not'
                                   f'found in the SymPortal count table: {diff_set}')
            else:
                # All samples were found in the df and those not listed can be exculded
                diff_set = set(sample_uids_included).difference(set(seq_count_df.index.values))
                seq_count_df = seq_count_df.drop(index=diff_set, inplace=False)
                print(f'Excluding {len(diff_set)} samples')
                return seq_count_df.reindex(sample_uids_included, axis=1)
        elif sample_uids_excluded is not None:
            print('Excluding samples according to user supplied sample_uids_excluded list')
            # Check that all the uids are found in the seq_count df
            if not set(sample_uids_excluded).issubset(set(seq_count_df.index.values)):
                diff_set = set(sample_uids_excluded).difference(set(seq_count_df.index.values))
                raise RuntimeError('The following uids that were specified in the sample_uids_excluded list are not'
                                   f'found in the SymPortal count table: {diff_set}')
            else:
                # All samples were found in the df and those not listed can be exculded
                print(f'Excluding {len(sample_uids_excluded)} samples')
                return seq_count_df.drop(index=sample_uids_excluded, inplace=False)
        elif sample_name_compiled_re_included is not None:
            print('Excluding samples according to user supplied sample_name_compiled_re_included regular expression')
            keep = []
            for sample_name in seq_count_df.index:
                if sample_name_compiled_re_included.match(sample_name):
                    keep.append(sample_name)
            if not keep:
                raise RuntimeError('No sample names matched the user supplied sample_name_compiled_re_included.\n')
            diff_set = set(keep).difference(seq_count_df.index.values)
            seq_count_df = seq_count_df.drop(index=diff_set)
            print(f'Excluding {len(diff_set)} samples')
            return seq_count_df.reindex(keep, axis=1)
        elif sample_name_compiled_re_excluded is not None:
            print('Excluding samples according to user supplied sample_name_compiled_re_excluded regular expression')
            exclude = []
            for sample_name in seq_count_df.index:
                if sample_name_compiled_re_excluded.match(sample_name):
                    exclude.append(sample_name)
            if not exclude:
                raise RuntimeError('No sample names matched the user supplied sample_name_compiled_re_excluded.\n')
            print(f'Excluding {len(exclude)} samples')
            return seq_count_df.drop(index=exclude)
        else:
            # Then they are all none and there are no samples to be excluded.
            return seq_count_df

    def _find_first_seq_position(self, seq_count_df):
        """
        Get the index of the first sequence column using a regex match rather
        than looking for a given metainfo column
        to ensure that changing the metadata down the line won't break this.

        :param seq_count_df (pandas.DataFrame): The dataframe we are finding the sequencing column for

        :returns index (int) of the column of the first sequence
        """
        seq_name_re = re.compile("^[A-I]{1}[0-9]{1}|^[0-9]+_[A-I]{1}")
        for i, col_str in enumerate(list(seq_count_df)):
            if seq_name_re.match(col_str):
                first_seq_index = i
                return first_seq_index
        raise RuntimeError('No valid sequence column found in df')

    def _check_valid_params(self, plot_type, profile_count_table_path, seq_count_table_path):
        if seq_count_table_path is None and profile_count_table_path is None:
            raise RuntimeError("Please provide either a seq_count_table_path or a profile_count_table_path")
        if plot_type is None:
            raise RuntimeError(
                "Please provide either a valid plot_type. Either seq_only, profile_only, or seq_and_profile"
            )

SPBars(
    seq_count_table_path='/Users/benjaminhume/Documents/projects/20210113_buitrago/sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt',
    plot_type='seq_only'
)