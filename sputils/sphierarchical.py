import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy
import scipy.spatial.distance
from datetime import datetime
from sputils import sputils
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import re
import matplotlib.gridspec as gridspec
from sputils import spcolors
import itertools


class SPHierarchical(sputils.SPUtils):
    """
    A utility script for creating hierarchical clustering dendrograms from SymPortal outputs
    Will take in a distance matrix SymPortal output and produce a dendrogram either as a new figure
    or on a specified axis object.

    dist_output_path (str) = Full path to a SymPortal distance output.

    labels (bool) = Whether to plot labels for each of the leafs [False]

    color_dict (dict) = A dict that maps sample uid or sample name to a matplotlib color. [None]

    tickness_dict (dict) = A dict maping sample_uid or sample name to a thickness in XXX. [None]

    orientation (str) = Whether to plot the dendrogram vetically or horizontally.
    Value must be 'vertical' or 'horizontal'. [horizontal]

    ax (matplotlib.axes.Axes) = If passed, the dedrogram will be plotted to this axis.
    Else a figure and axis will be created. [None]

    no_plotting (bool): If True, then no figure or ax object will be created. This is useful when you want to
    run the __init__ of the class to get access to its objects without producing figures. [False]

    save_fig (bool): If passed the figure will be output. This will only have an effect if ax is None. [False]

    output_dir (str): Directory to save the figure to. This will only have an effect if save_fig is True. [cwd]

    sample_uids_included (list<int>): A list of the sample uids that should be plotted.
    For bar plots, the samples will be plotted in this order.[None]

    sample_uids_excluded (list<int>): A list of the sample uids that should be excluded from plotting. [None]

    sample_names_included (list<str>): A list of the sample names that should be plotted.
    For bar plots, the samples will be plotted in this order.[None]

    sample_names_excluded (list<str>): A list of the sample names that should be excluded from plotting. [None]

    sample_name_compiled_re_included (compiled regex): A regular expression, which if matched by a sample name,
    will denote that the sample should be included in plotting. [None]

    sample_name_compiled_re_included (compiled regex): A regular expression, which if matched by a sample name,
    will denote that the sample should be included in plotting. [None]
    """


    def __init__(
            self, dist_output_path=None, sample_profile=None, labels=False,
            color_dict=None, thickness_dict=None, orientation='h', ax=None, save_fig=False, output_dir=None,
            sample_uids_included=None, sample_uids_excluded=None,
            sample_names_included=None, sample_names_excluded=None,
            sample_name_compiled_re_included=None, sample_name_compiled_re_excluded=None, no_plotting=False
    ):

        self._check_exclude_arguments(
            sample_name_compiled_re_excluded, sample_name_compiled_re_included,
            sample_names_excluded, sample_names_included,
            sample_uids_excluded, sample_uids_included
        )

        self.labels = labels
        self.color_dict = color_dict
        self.thickness_dict = thickness_dict
        self.date_time = str(datetime.now()).split('.')[0].replace('-', '').replace(' ', 'T').replace(':', '')
        self.save_fig = save_fig
        self.output_dir = output_dir
        self.no_plotting = no_plotting

        if orientation in ['h', 'horizontal']:
            self.orientation = 'h'
        elif orientation in ['v', 'vertical']:
            self.orientation = 'v'
        else:
            raise RuntimeError("Orientation must be either 'horizontal' or 'vertical'.")

        if sample_profile is not None:
            if not sample_profile in ['sampe', 'profile']:
                raise RuntimeError("sample_profile must be either 'sample' or 'profile'.")
        else:
            if 'sample' in dist_output_path:
                print('sample distance data detected')
            elif 'profile' in dist_output_path:
                print('profile distance data deteceted')
            else:
                raise RuntimeError('Unable to automatically detect whether this is sample of profile data.\n'
                                   'Please explicitly set sample_profile to "sample" or "profile".\n')

        with open(dist_output_path, 'r') as f:
            f_list = [_.split('\t') for _ in f]
        # In the older .dist files there is a single line at the top that holds the number of samples
        # We need to delete this
        if len(f_list[0]) == 1:
            f_list = f_list[1:]
        self.obj_name_to_obj_uid_dict = {_[0]: int(_[1]) for _ in f_list}
        self.obj_uid_to_obj_name_dict = {v: k for k, v in self.obj_name_to_obj_uid_dict.items()}
        self.dist_df = pd.DataFrame(
            data=[_[2:] for _ in f_list],
            columns=list(self.obj_uid_to_obj_name_dict.keys()),
            index=list(self.obj_uid_to_obj_name_dict.keys())
        ).astype(float)

        self.dist_df = self._exclude_samples_from_count_df(
            df=self.dist_df,
            uid_to_name_dict=self.obj_name_to_obj_uid_dict, name_to_uid_dict=self.obj_name_to_obj_uid_dict,
            sample_uids_included=sample_uids_included, sample_uids_excluded=sample_uids_excluded,
            sample_names_included=sample_names_included, sample_names_excluded=sample_names_excluded,
            sample_name_compiled_re_included=sample_name_compiled_re_included,
            sample_name_compiled_re_excluded=sample_name_compiled_re_excluded, dists=True
        )


        # Figure setup
        if not self.no_plotting:
            self.fig, self.ax = self._setup_axis(ax)
        foo = 'bar'

    def plot(self):
        """Plot the dendrogram"""
        if self.no_plotting:
            raise RuntimeError("no_plotting is enabled. Please disable this option to enable plotting.")
        condensed_dist = scipy.spatial.distance.squareform(self.dist_df)
        # this creates the linkage df that will be passed into the dendogram_sp function
        linkage = scipy.cluster.hierarchy.linkage(y=condensed_dist, optimal_ordering=True)
        self.dendrogram = scipy.cluster.hierarchy.dendrogram(
            linkage, labels=(list(self.dist_df)), ax=self.ax, leaf_rotation=90, no_labels=not
            self.labels, link_color_func=lambda k: '#000000', orientation='top', get_leaves=True,
            show_leaf_counts=True
        )

        # Out put the sample order of the dendrogram for easy access
        self.dendrogram_sample_order_uid = self.dendrogram['ivl']
        self.dendrogram_sample_order_name = [self.obj_uid_to_obj_name_dict[uid] for uid in self.dendrogram['ivl']]

        if self.save_fig:
            if self.output_dir is None:
                plt.savefig(f'{self.date_time}_dendro.svg')
                plt.savefig(f'{self.date_time}_dendro.png', dpi=1200)
            else:
                plt.savefig(os.path.join(self.output_dir, f'{self.date_time}_dendro.svg'))
                plt.savefig(os.path.join(self.output_dir, f'{self.date_time}_dendro.png'), dpi=1200)

        foo = 'bar'

    def _setup_axis(self, ax):
        if ax is None:
            if self.orientation == 'v':
                # One col wide = 89mm
                # Height is 5mm per sample
                height = (1 * len(self.dist_df.index))
                width = 89
                return plt.subplots(
                    ncols=1, nrows=1, figsize=(self._mm2inch(width, height))
                )
            else:
                # Standard 80mm high
                # Width is max of 89mm or 5*num_samples
                width = max(89, (1 * len(self.dist_df.index)))
                height = 80
                return plt.subplots(
                    ncols=1, nrows=1, figsize=(self._mm2inch(width, height))
                )
        else:
            return ax.figure, ax

    def _mm2inch(self, *tupl):
        inch = 25.4
        if isinstance(tupl[0], tuple):
            return tuple(i / inch for i in tupl[0])
        else:
            return tuple(i / inch for i in tupl)
