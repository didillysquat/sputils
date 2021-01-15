#!/usr/bin/env python3

from sputils.spbars import SPBars
import argparse


class RunSPBars:
    def __init__(self):
        parser = argparse.ArgumentParser(
                    description='A script for plotting stacked bar plots from SymPortal outputs.',
                    epilog='For support email: didillysquat@gmail.com')
        parser.add_argument(
            '--seq-count-table-path',
            help="The full path to the 'absolute.abund_and_meta' post-MED sequence count table to use for plotting."
                 "N.B. Plotting of the pre-MED count table is currently not supported.",
            required=False, type=str, default=None
        )
        parser.add_argument(
            '--profile-count-table-path',
            help="The full path to the absolute 'absolute.abund_and_meta' "
                 "ITS2 type profile count table to use for plotting.",
            required=False, type=str, default=None
        )
        parser.add_argument(
            '--profile-count-table-path',
            help="The full path to the absolute 'absolute.abund_and_meta' "
                 "ITS2 type profile count table to use for plotting.",
            required=False, type=str, default=None
        )
        raise NotImplementedError()


        self.args = parser.parse_args()

    def run(self):
        SPBars()


RunSPBars()