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

        seqs_right_bottom (bool): When True, for seq_and_profile plots the sequences will be plotted to the right of the
        profiles for vertical plots, and below the profiles for horizontal plots. [False]

        reverse_seq_abund (bool): When True, the order in which the sequences are plotted will be reversed. [False]

        reverse_profile_abund (bool): When True, the order in which the profiles are plotted will be reversed. [False]

        color_by_genus (bool): If True, sequences and profiles will be colored according to their genus. [False]

        sample_outline (bool): If True, each sample will be separated by a black line. [False]

    Returns:
        tuple(matplotlib.pyplot.figure, matplotlib.axes.Axes): The figure and axes object that contain the plot

    """
    def __init__(
            self, seq_count_table_path=None, profile_count_table_path=None, plot_type='seq', figsize=None,
            sample_uids_included=None, sample_uids_excluded=None,
            sample_name_compiled_re_included=None, sample_name_compiled_re_excluded=None,
            orientation='h', legend=True,
            relative_abundnce=True, num_seq_leg_cols=20, num_profile_leg_cols=20, seqs_right_bottom=False,
            reverse_seq_abund=False, reverse_profile_abund=False, color_by_genus=False, sample_outline=False
    ):

        # arguments that will be used throughout the class
        self.plot_type = plot_type
        if orientation in ['v', 'vertical']:
            self.orientation = 'v'
        elif orientation in ['h', 'horizontal']:
            self.orientation = 'h'
        self.legend = legend
        self.seqs_right_bottom = seqs_right_bottom
        self.reverse_seq_abund = reverse_seq_abund
        self.reverse_profile_abund = reverse_profile_abund
        self.color_by_genus = color_by_genus
        self.num_profile_leg_cols = num_profile_leg_cols
        self.num_seq_leg_cols = num_seq_leg_cols
        self.relative_abundance = relative_abundnce
        self.sample_outline = sample_outline

        # Check the inputs
        self._check_path_exists([seq_count_table_path, profile_count_table_path])
        if self.plot_type not in ['seq_only', 'profile_only', 'seq_and_profile']:
            raise RuntimeError("plot_type must be one of 'seq_only', 'profile_only', or 'seq_and_profile'")
        if self.orientation not in ['h', 'v', 'horizontal', 'vertical']:
            raise RuntimeError("orientation must be one of 'h', 'v', 'horizontal', or 'vertical' ")
        if not type(self.legend) == bool:
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
        self._check_valid_params(profile_count_table_path, seq_count_table_path)

        if seq_count_table_path:
            (
                self.sample_name_to_sample_uid_dict,
                self.sample_uid_to_sample_name_dict,
                self.seq_count_df
            ) = self._make_seq_count_df(
                seq_count_table_path,
                sample_uids_included, sample_uids_excluded,
                sample_name_compiled_re_included, sample_name_compiled_re_excluded
            )
        if profile_count_table_path:
            (
                self.sample_name_to_sample_uid_dict, self.sample_uid_to_sample_name_dict,
                self.profile_name_to_profile_uid_dict, self.profile_uid_to_profile_name_dict,
                self.profile_count_df
            ) = self._make_profile_count_df(
                    profile_count_table_path,
                    sample_uids_included, sample_uids_excluded,
                    sample_name_compiled_re_included, sample_name_compiled_re_excluded
            )

        # Figure setup
        self._setup_fig_and_ax(figsize)

        # Patch setup
        # List for holding the rectangle patches used in making the bars
        self.bar_patches = []

        # Plotting colors setup
        self.grey_iterator = itertools.cycle(spcolors.greys)
        self.genus_color_dict = spcolors.genus_color_dict

        if plot_type in ['seq_only', 'seq_and_profile']:
            # TODO we will only need one legend if we are coloring by genus
            if not self.color_by_genus:
                self.color_hash_iterator = iter(spcolors.color_list)
                self.pre_def_seq_color_dict = spcolors.pre_def_color_dict
            self.seq_color_dict = self._make_seq_color_dict()

        if plot_type in ['profile_only', 'seq_and_profile']:
            if not self.color_by_genus:
                self.profile_col_generator = (
                    '#%02x%02x%02x' % rgb_tup for rgb_tup in
                    spcolors.create_color_list(
                        mix_col=(255, 255, 255),
                        sq_dist_cutoff=5000,
                        num_cols=self.num_profile_leg_cols,
                        time_out_iterations=10000,
                        warnings_off=True)
                )
            self.profile_color_dict = self._make_profile_color_dict()

    def _make_profile_color_dict(self):
        if self.color_by_genus:
            return spcolors.genus_color_dict
        else:
            prof_color_dict = {}
            for prof_uid in list(self.profile_count_df):
                try:
                    prof_color_dict[prof_uid] = next(self.profile_col_generator)
                except StopIteration:
                    prof_color_dict[prof_uid] = next(self.grey_iterator)
            return prof_color_dict

    def plot(self):
        self._plot_bars()
        plt.show()
        plt.savefig('figure_1.svg')
        plt.savefig('figure_1.png', dpi=1200)
        foo = 'bar'

    def _plot_bars(self):
        """
            Plots every instance of a sequence or profile on a per sample basis.
            Each sample is either a vertical or horizontal stacked bar depending on the orientation
            with each rectangle of the stack representing
            a profile or sequences instance in that sample.

            :param color_dict: color dict used for either the sequence or profiles being plotted
            :param df: the abundance df used for the plotting
            :return: None
            """
        color_list = []
        if self.plot_type == 'seq_only':
            self._make_rect_and_cols(color_dict=self.seq_color_dict, color_list=color_list, df=self.seq_count_df)
            self._set_ax_lims_and_outline_only_plot(df=self.seq_count_df)
        elif self.plot_type == 'profile_only':
            self._make_rect_and_cols(
                color_dict=self.profile_color_dict, color_list=color_list, df=self.profile_count_df
            )
            self._set_ax_lims_and_outline_only_plot(df=self.profile_count_df)
        else:
            if self.seqs_right_bottom:
                self._make_rect_and_cols(color_dict=self.seq_color_dict, color_list=color_list, df=self.seq_count_df,
                                         pos_neg='positive')
                self._make_rect_and_cols(
                    color_dict=self.profile_color_dict, color_list=color_list, df=self.profile_count_df,
                    pos_neg='negative'
                )
                self._set_ax_lims_both_plot()
            else:
                self._make_rect_and_cols(color_dict=self.seq_color_dict, color_list=color_list, df=self.seq_count_df, pos_neg='negative')
                self._make_rect_and_cols(
                    color_dict=self.profile_color_dict, color_list=color_list, df=self.profile_count_df, pos_neg='positive'
                )
                self._set_ax_lims_both_plot()

        listed_color_map = ListedColormap(color_list)
        patches_collection = PatchCollection(self.bar_patches, cmap=listed_color_map)
        patches_collection.set_array(np.arange(len(self.bar_patches)))
        self.bar_ax.add_collection(patches_collection)
        self.bar_ax.autoscale_view()
        self.fig.canvas.draw()

    def _set_ax_lims_and_outline_only_plot(self, df):
        if self.orientation == 'v':
            self.bar_ax.set_ylim(-0.5, len(df) - 0.5)
            self.bar_ax.set_xlim(0, df.to_numpy().max())
            if self.sample_outline:
                for y in np.arange(0.5, len(df.index) - 1.5):
                    # TODO dynamically program the lw
                    self.bar_ax.axhline(y=y, lw=0.1, c='black')
        else:
            self.bar_ax.set_ylim(0, df.to_numpy().max())
            self.bar_ax.set_xlim(-0.5, len(df) - 0.5)
            if self.sample_outline:
                for x in np.arange(0.5, len(df.index) - 1.5):
                    #TODO dynamically program the lw
                    self.bar_ax.axvline(x=x, lw=0.1, c='black')

    def _set_ax_lims_both_plot(self):
        if self.orientation == 'v':
            self.bar_ax.set_ylim(-0.5, len(self.seq_count_df.index) - 0.5)
            if self.seqs_right_bottom:
                self.bar_ax.set_xlim(-1 * self.profile_count_df.to_numpy().max(), self.seq_count_df.to_numpy().max())
            else:
                self.bar_ax.set_xlim(-1 * self.profile_count_df.to_numpy().max(), self.seq_count_df.to_numpy().max())
        else:
            self.bar_ax.set_xlim(-0.5, len(self.seq_count_df.index) - 0.5)
            if self.seqs_right_bottom:
                self.bar_ax.set_ylim(-1 * self.profile_count_df.to_numpy().max(), self.seq_count_df.to_numpy().max())
            else:
                self.bar_ax.set_ylim(-1 * self.profile_count_df.to_numpy().max(), self.seq_count_df.to_numpy().max())

    def _make_rect_and_cols(self, color_dict, color_list, df, pos_neg='positive'):
        index_for_plot = 0
        for sample_uid in df.index:
            bottom = 0
            non_zero_seq_abundances = df.loc[sample_uid][
                df.loc[sample_uid] > 0]
            if pos_neg == 'negative':
                non_zero_seq_abundances = non_zero_seq_abundances * -1
            for obj_uid, abund in non_zero_seq_abundances.iteritems():
                if self.orientation in ['horizontal', 'h']:
                    self.bar_patches.append(Rectangle(
                        (index_for_plot - 0.5, bottom),
                        1,
                        abund, color=color_dict[obj_uid]))
                else:
                    coords = (bottom, index_for_plot - 0.5)
                    self.bar_patches.append(Rectangle(
                        (bottom, index_for_plot - 0.5),
                        abund,
                        1, color=color_dict[obj_uid]))
                bottom += abund
                color_list.append(color_dict[obj_uid])
            index_for_plot += 1

    def _make_seq_color_dict(self):
        if self.color_by_genus:
            return spcolors.genus_color_dict
        else:
            seq_color_dict = {}
            for seq_name in list(self.seq_count_df):
                if seq_name in self.pre_def_seq_color_dict:
                    seq_color_dict[seq_name] = self.pre_def_seq_color_dict[seq_name]
                else:
                    try:
                        seq_color_dict[seq_name] = next(self.color_hash_iterator)
                    except StopIteration:
                        seq_color_dict[seq_name] = next(self.grey_iterator)
            return seq_color_dict

    def _setup_fig_and_ax(self, figsize):
        """
        Setup up the ax objects that the plot and legends will be plotted to.

        :param figsize: user supplied figure size tuple(<int>, <int>) in mm

        :return: None. But, self.fig, self.bar_ax, self.leg_ax_one and self.leg_ax_two (if plot_type is 'seq_and_profile')
        will be set.
        """
        if self.plot_type == 'seq_and_profile':
            # Then we need to have two legend axes if legend plotted
            self._setup_seq_and_profile(figsize)
        else:
            # Then we are working with a single legend if legend plotted
            self._setup_seq_or_profile_only_plot(figsize)

    def _setup_seq_or_profile_only_plot(self, figsize):
        if self.orientation == 'v':
            # Then we want the bar plot to sit next to the legend plots
            # The bars will span the full height of the figure
            if figsize:
                # TODO enforce a minimum size
                self.fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be one column so 89mm
                # height will dependent on sample number
                # We will need a min height for the one legend of 30
                # and 5mm for each sample
                height = max(30, (5 * len(self.seq_count_df.index)))
                width = 89
                self.fig = plt.figure(figsize=self._mm2inch((width, height)))
            gs = gridspec.GridSpec(1, 2)
            if self.legend:
                self.bar_ax = plt.subplot(gs[:, :1])
                self.leg_ax_one = plt.subplot(gs[:, 1:2])
            else:
                self.bar_ax = plt.subplot(gs[:, :])
        else:
            # Then we want the bar plot to sit above the legend plot
            if figsize:
                # TODO enforce a minimum size
                self.fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be two columns so 183
                # height will be fixed, 30 for legend 50 for bars
                height = 80
                width = 183
                self.fig = plt.figure(figsize=self._mm2inch((width, height)))
                # TODO implement multiple rows of barplots
            gs = gridspec.GridSpec(2, 1)
            if self.legend:
                self.bar_ax = plt.subplot(gs[:1, :])
                self.leg_ax_one = plt.subplot(gs[1:2, :])
            else:
                self.bar_ax = plt.subplot(gs[:, :])

    def _setup_seq_and_profile(self, figsize):
        if self.orientation == 'v':
            # Then we want the bar plot to sit next to the legend plots
            # The bars will span the full height of the figure
            # The seq legend will span the top half
            # The profile legend will span the bottom half
            if figsize:
                # TODO enforce a minimum size
                self.fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be one column so 89mm
                # height will dependent on sample number
                # We will need a min height for the legends of 30 * 2
                # and 5mm for each sample
                height = max(60, (5 * len(self.seq_count_df.index)))
                width = 89
                self.fig = plt.figure(figsize=self._mm2inch((width, height)))
            gs = gridspec.GridSpec(2, 2)
            if self.legend:
                self.bar_ax = plt.subplot(gs[:2, :1])
                self.leg_ax_one = plt.subplot(gs[:1, 1:2])
                self.leg_ax_two = plt.subplot(gs[1:2, 1:2])
            else:
                self.bar_ax = plt.subplot(gs[:, :])
        else:
            # Then we want the bar plot to sit above the legend plots
            # The bars will span the full width of the figure
            # The seq legend will span the left half
            # The profile legend will span the right half
            if figsize:
                # TODO enforce a minimum size
                self.fig = plt.figure(figsize=self._mm2inch(figsize))
            else:
                # Deduce the fig size automatically
                # width will be two columns so 183
                # height will be fixed, 30 for legend 50 for bars
                height = 80
                width = 183
                self.fig = plt.figure(figsize=self._mm2inch((width, height)))
                # TODO implement multiple rows of barplots
            gs = gridspec.GridSpec(2, 2)
            if self.legend:
                self.bar_ax = plt.subplot(gs[:1, :2])
                self.leg_ax_one = plt.subplot(gs[1:2, :1])
                self.leg_ax_two = plt.subplot(gs[1:2, 1:2])
            else:
                self.bar_ax = plt.subplot(gs[:, :])

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

    def _make_profile_count_df(
            self, profile_count_table_path,
                sample_uids_included, sample_uids_excluded,
                sample_name_compiled_re_included, sample_name_compiled_re_excluded
    ):
        """
        make the seq_count_df, and produce dicts linking sample_uid to sample_name in both directions
        seq_count_df will be only the absolute abundance, with sequence_name as columns and sample_uid as index
        seq_count_df columns will be ordered by abundance across all samples
        samples will be excluded according to the user specified arguments
        """
        profile_count_df_abund = pd.read_csv(profile_count_table_path, sep='\t', skiprows=[1,2,3,4,5,6], skipfooter=2, engine='python')
        profile_count_df_meta = pd.read_csv(profile_count_table_path, sep='\t', index_col=0)
        profile_count_df_meta = profile_count_df_meta.drop(profile_count_df_meta.columns[0], axis=1)
        profile_uid_to_profile_name_dict = {
            p_uid: p_name for p_uid, p_name in profile_count_df_meta.loc['ITS2 type profile'].items()
        }
        profile_name_to_profile_uid_dict = {
            p_name: p_uid for p_uid, p_name in profile_count_df_meta.loc['ITS2 type profile'].items()
        }

        sample_uid_to_sample_name_dict = {
            uid: p_name for uid, p_name in
            zip(profile_count_df_abund[profile_count_df_abund.columns[0]].values, profile_count_df_abund[profile_count_df_abund.columns[1]].values)
        }
        sample_name_to_sample_uid_dict = {
            p_name: uid for uid, p_name in
            zip(profile_count_df_abund[profile_count_df_abund.columns[0]].values, profile_count_df_abund[profile_count_df_abund.columns[1]].values)
        }

        profile_count_df_abund.set_index(keys='ITS2 type profile UID', drop=True, inplace=True)
        profile_count_df_abund = profile_count_df_abund.drop(profile_count_df_abund.columns[0], axis=1)

        profile_count_df_abund = self._exclude_samples_from_count_df(
            profile_count_df_abund, sample_uid_to_sample_name_dict, sample_name_to_sample_uid_dict,
            sample_uids_included, sample_uids_excluded,
            sample_name_compiled_re_included, sample_name_compiled_re_excluded
        )

        # If color_by_genus convert to genera abunances
        if self.color_by_genus:
            present_genera = set([prof_name[0] for prof_name in profile_name_to_profile_uid_dict.keys()])
            for i, clade in enumerate(present_genera):
                profile_count_df_abund[clade] = profile_count_df_abund[
                    [
                        profile_uid for profile_uid in list(profile_count_df_abund)[:(-1 * i)] if
                        profile_uid_to_profile_name_dict[profile_uid].startswith(clade)
                    ]
                ].sum(axis=1)
            # drop non-clade cols
            profile_count_df_abund = profile_count_df_abund.drop(list(profile_count_df_abund)[:(-1 * len(present_genera))], axis=1)

        # Convert to relative abundances
        if self.relative_abundance:
            profile_count_df_abund = profile_count_df_abund.div(
                profile_count_df_abund.sum(axis=1), axis=0
            )

        # sort by column abundances, highest first
        if self.reverse_profile_abund:
            sorted_cols = profile_count_df_abund.sum(axis=0).sort_values(ascending=False).index
        else:
            sorted_cols = profile_count_df_abund.sum(axis=0).sort_values(ascending=False).index
        profile_count_df_abund = profile_count_df_abund.reindex(sorted_cols, axis=1)

        return (
            sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict,
            profile_name_to_profile_uid_dict, profile_uid_to_profile_name_dict,
            profile_count_df_abund
        )

    def _make_seq_count_df(
            self, seq_count_table_path,
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

        seq_count_df = self._exclude_samples_from_count_df(
            seq_count_df, sample_uid_to_sample_name_dict, sample_name_to_sample_uid_dict,
            sample_uids_included, sample_uids_excluded,
            sample_name_compiled_re_included, sample_name_compiled_re_excluded
        )

        # If color_by_genus convert to genera abunances
        if self.color_by_genus:
            present_genera = set([_[0] if _[0] in list('ABCDEFGHI') else _[-1] for _ in list(seq_count_df)])
            for clade in present_genera:
                seq_count_df[clade] = seq_count_df[
                    [
                        seq_name for seq_name in list(seq_count_df) if
                        (seq_name.startswith(clade) or seq_name.endswith(clade))
                    ]
                ].sum(axis=1)
            # drop non-clade cols
            seq_count_df = seq_count_df.drop(list(seq_count_df)[:(-1 * len(present_genera))], axis=1)

        # Convert to relative abundances
        if self.relative_abundance:
            seq_count_df = seq_count_df.div(
                seq_count_df.sum(axis=1), axis=0
            )

        # sort by seq abundance, highest first
        if self.reverse_seq_abund:
            sorted_seq_index = seq_count_df.sum(axis=0).sort_values(ascending=True).index
        else:
            sorted_seq_index = seq_count_df.sum(axis=0).sort_values(ascending=False).index
        seq_count_df = seq_count_df.reindex(sorted_seq_index, axis=1)

        return (sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict, seq_count_df)

    def _exclude_samples_from_count_df(
            self, count_df, sample_uid_to_sample_name_dict, sample_name_to_sample_uid_dict,
            sample_uids_included, sample_uids_excluded,
            sample_name_compiled_re_included, sample_name_compiled_re_excluded
    ):
        """
        Due to the checks we have already performed on the include/exclude arguments
        only one of them will be not None. This is the set we will use to exclude samples from the df.
        If all of them are none, then there are no samples to exclude.
        """
        if sample_uids_included is not None:
            return self._exclude_samples_sample_uids_included(sample_uids_included, count_df)
        elif sample_uids_excluded is not None:
            return self._exclude_samples_sample_uids_excluded(count_df, sample_uids_excluded)
        elif sample_name_compiled_re_included is not None:
            return self._exclude_samples_sample_name_compiled_re_included(
                count_df, sample_name_compiled_re_included,
                sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict
            )
        elif sample_name_compiled_re_excluded is not None:
            return self._exclude_samples_sample_name_compiled_re_excluded(
                count_df, sample_name_compiled_re_excluded,
                sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict
            )
        else:
            # Then they are all none and there are no samples to be excluded.
            return count_df

    def _exclude_samples_sample_name_compiled_re_excluded(self, count_df, sample_name_compiled_re_excluded,
                                                          sample_name_to_sample_uid_dict,
                                                          sample_uid_to_sample_name_dict):
        print('Excluding samples according to user supplied sample_name_compiled_re_excluded regular expression')
        exclude = []
        for sample_name in [sample_uid_to_sample_name_dict[sample_uid] for sample_uid in count_df.index]:
            if sample_name_compiled_re_excluded.match(sample_name):
                exclude.append(sample_name_to_sample_uid_dict[sample_name])
        if not exclude:
            raise RuntimeError('No sample names matched the user supplied sample_name_compiled_re_excluded.\n')
        print(f'Excluding {len(exclude)} samples')
        return count_df.drop(index=exclude)

    def _exclude_samples_sample_name_compiled_re_included(
            self, count_df, sample_name_compiled_re_included,
            sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict
    ):
        print('Excluding samples according to user supplied sample_name_compiled_re_included regular expression')
        keep = []
        for sample_name in [sample_uid_to_sample_name_dict[sample_uid] for sample_uid in count_df.index]:
            if sample_name_compiled_re_included.match(sample_name):
                keep.append(sample_name_to_sample_uid_dict[sample_name])
        if not keep:
            raise RuntimeError('No sample names matched the user supplied sample_name_compiled_re_included.\n')
        diff_set = set(keep).difference(count_df.index.values)
        count_df = count_df.drop(index=diff_set)
        print(f'Excluding {len(diff_set)} samples')
        return count_df.reindex(keep, axis=1)

    def _exclude_samples_sample_uids_excluded(self, count_df, sample_uids_excluded):
        print('Excluding samples according to user supplied sample_uids_excluded list')
        # Check that all the uids are found in the seq_count df
        if not set(sample_uids_excluded).issubset(set(count_df.index.values)):
            diff_set = set(sample_uids_excluded).difference(set(count_df.index.values))
            raise RuntimeError('The following uids that were specified in the sample_uids_excluded list are not'
                               f'found in the SymPortal count table: {diff_set}')
        else:
            # All samples were found in the df and those not listed can be exculded
            print(f'Excluding {len(sample_uids_excluded)} samples')
            return count_df.drop(index=sample_uids_excluded, inplace=False)

    def _exclude_samples_sample_uids_included(self, sample_uids_included, seq_count_df):
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

    def _check_valid_params(self, profile_count_table_path, seq_count_table_path):
        if seq_count_table_path is None and profile_count_table_path is None:
            raise RuntimeError("Please provide either a seq_count_table_path or a profile_count_table_path")
        if self.plot_type is None:
            raise RuntimeError(
                "Please provide either a valid plot_type. Either seq_only, profile_only, or seq_and_profile"
            )

SPBars(
    seq_count_table_path='/Users/benjaminhume/Documents/projects/20210113_buitrago/sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt',
    profile_count_table_path='/Users/benjaminhume/Documents/projects/20210113_buitrago/sp_output/its2_type_profiles/131_20201203_DBV_20201207T095144.profiles.absolute.abund_and_meta.txt',
    plot_type='seq_only', orientation='v', legend=False, relative_abundnce=True,
    color_by_genus=True, sample_outline=True
).plot()