#!/usr/bin/env python3
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
    ['seq_only]'

    figsize (tuple): A tuple that will be passed to the plt.figure figure creation. Units should be in mm.
    If not provided this will automatically deduced. [None]

    sample_uids_included (list<int>): A list of the sample uids that should be plotted. [None]

    sample_uids_excluded (list<int>): A list of the sample uids that should be excluded from plotting. [None]

    sample_name_compiled_re_included (compiled regex): A regular expression, which if matched by a sample name,
    will denote that the sample should be included in plotting. [None]

    sample_name_compiled_re_included (compiled regex): A regular expression, which if matched by a sample name,
    will denote that the sample should be included in plotting. [None]

    orientation (str): 'vertical' or 'horizontal'. 'v' or 'h' may be passed. ['h']

    legend (bool): whether to plot a legend

Returns:
    tuple(matplotlib.pyplot.figure, matplotlib.axes.Axes): The figure and axes object that contain the plot

"""

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

class SPBars:
    def __init__(
            self, seq_count_table_path=None, profile_count_table_path=None, plot_type='seq', figsize=None,
            sample_uids=None, sample_uids_included=None, sample_uids_excluded=None,
            sample_name_compiled_re_included=None, sample_name_compiled_re_excluded=None, orientation='h', legend=True):

        self._check_valid_params(plot_type, profile_count_table_path, seq_count_table_path)

        if seq_count_table_path:

            (
                self.sample_name_to_sample_uid_dict,
                self.sample_uid_to_sample_name_dict,
                self.seq_count_df
            ) = self._make_seq_count_df(seq_count_table_path)


        # Figure setup
        if figsize is not None:
            self.fig, self.ax = plt.figure(figsize=figsize)
        else:
            # if vertical, 5mm height per sample on top of
            self.fif, self.ax = plt.subplots(1, 1, figsize=self._mm2inch(89, 89))
        foo = 'bar'

    def _mm2inch(self, *tupl):
        inch = 25.4
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)

    def _make_seq_count_df(self, seq_count_table_path):
        """make the seq_count_df, and produce dicts linking sample_uid to sample_name in both directions"""
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
        return (sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict, seq_count_df)

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