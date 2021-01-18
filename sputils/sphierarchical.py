import os
import pandas as pd
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
import re
import matplotlib.gridspec as gridspec
from sputils import spcolors
import itertools
from datetime import datetime

class SPHierarchical:
    """
    A utility script for creating hierarchical clustering dendrograms from SymPortal outputs
    Will take in a distance matrix SymPortal output and produce a dendrogram either as a new figure
    or on a specified axis object.

    dist_output_path (str) = Full path to a SymPortal distance output.

    labels (bool) = Whether to plot labels for each of the leafs [False]

    color_dict (dict) = A dict that maps sample uid or sample name to a matplotlib color

    tickness_dict (dict) = A dict maping sample_uid or sample name to a thickness in XXX
    """
    def __init__(self, dist_output_path, labels=False, color_dict=None, thickness_dict=None):
        self.dist_path = dist_output_path
        self.labels = labels
        self.color_dict = color_dict
        self.thickness_dict = thickness_dict