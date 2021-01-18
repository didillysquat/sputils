class SPUtils:
    """Base class for the sputil modules"""
    def __init__(self):
        pass

    @staticmethod
    def _check_exclude_arguments(sample_name_compiled_re_excluded, sample_name_compiled_re_included,
                                 sample_names_excluded, sample_names_included, sample_uids_excluded,
                                 sample_uids_included):
        c_count = 0
        for constraint in [
            sample_uids_included, sample_uids_excluded,
            sample_names_included, sample_names_excluded,
            sample_name_compiled_re_included, sample_name_compiled_re_excluded
        ]:
            if constraint is not None:
                c_count += 1
                if c_count > 1:
                    raise RuntimeError('Please provide only one of:\n'
                                       '\tsample_uids_included\n'
                                       '\tsample_uids_excluded\n'
                                       '\tsample_names_included\n'
                                       '\tsample_names_excluded\n'
                                       '\tsample_name_compiled_re_included\n'
                                       '\tsample_name_compiled_re_excluded\n')

    def _exclude_samples_from_count_df(
            self, df, uid_to_name_dict, name_to_uid_dict,
            sample_uids_included, sample_uids_excluded, sample_names_included, sample_names_excluded,
            sample_name_compiled_re_included, sample_name_compiled_re_excluded, dists=False
    ):
        """
        Due to the checks we have already performed on the include/exclude arguments
        only one of them will be not None. This is the set we will use to exclude samples from the df.
        If all of them are none, then there are no samples to exclude.

        :param dists: if True then we are working with the distance df and we should remove
        columns and index.
        """
        if sample_uids_included is not None:
            return self._exclude_samples_sample_uids_included(sample_uids_included, df, dists)
        elif sample_uids_excluded is not None:
            return self._exclude_samples_sample_uids_excluded(df, sample_uids_excluded, dists)
        elif sample_names_included is not None:
            return self._exclude_samples_samples_names_included(df, name_to_uid_dict,
                                                                sample_names_included, dists)
        elif sample_names_excluded is not None:
            return self._exclude_samples_samples_names_excluded(df, name_to_uid_dict,
                                                                sample_names_excluded, sample_names_included, dists)
        elif sample_name_compiled_re_included is not None:
            return self._exclude_samples_sample_name_compiled_re_included(
                df, sample_name_compiled_re_included,
                name_to_uid_dict, uid_to_name_dict, dists
            )
        elif sample_name_compiled_re_excluded is not None:
            return self._exclude_samples_sample_name_compiled_re_excluded(
                df, sample_name_compiled_re_excluded,
                name_to_uid_dict, uid_to_name_dict, dists
            )
        else:
            # Then they are all none and there are no samples to be excluded.
            return df

    def _exclude_samples_samples_names_excluded(
            self, count_df, sample_name_to_sample_uid_dict, sample_names_excluded, sample_names_included, dist
    ):
        # Return df containing only the included samples sorted according to that sampple order
        print('Excluding samples according to user supplied sample_names_included list')
        # Check that all uids are found in the seq_count_df
        diff_set = set(sample_names_excluded).difference(set(sample_name_to_sample_uid_dict.keys()))
        if not set(sample_names_included).issubset(set(sample_name_to_sample_uid_dict.keys())):
            raise RuntimeError('The following uids that were specified in the sample_uids_excluded list are not'
                               f'found in the SymPortal count table: {diff_set}')
        else:
            # All samples were found in the df and those listed can be exculded
            print(f'Excluding {len(diff_set)} samples')
            return self._drop_from_df(
                df=count_df,
                labels=[sample_name_to_sample_uid_dict[sample_name] for sample_name in sample_names_excluded],
                dist=dist)

    def _exclude_samples_samples_names_included(
            self, count_df, sample_name_to_sample_uid_dict, sample_names_included, dist
    ):
        # Return df containing only the included samples sorted according to that sampple order
        print('Excluding samples according to user supplied sample_names_included list')
        # Check that all uids are found in the seq_count_df

        if not set(sample_names_included).issubset(set(sample_name_to_sample_uid_dict.keys())):
            diff_set = set(sample_names_included).difference(set(sample_name_to_sample_uid_dict.keys()))
            raise RuntimeError('The following uids that were specified in the sample_uids_included list are not'
                               f'found in the SymPortal count table: {diff_set}')
        else:
            diff_set = set(sample_names_included).symmetric_difference(set(sample_name_to_sample_uid_dict.keys()))
            # All samples were found in the df and those not listed can be exculded
            count_df = self._drop_from_df(
                df=count_df, labels=[sample_name_to_sample_uid_dict[sample_name] for sample_name in diff_set],
                dist=dist)
            print(f'Excluding {len(diff_set)} samples')
            return self._reindex_df(
                df=count_df,
                labels=[sample_name_to_sample_uid_dict[sample_name] for sample_name in sample_names_included],
                dist=dist)

    def _exclude_samples_sample_name_compiled_re_excluded(self, count_df, sample_name_compiled_re_excluded,
                                                          sample_name_to_sample_uid_dict,
                                                          sample_uid_to_sample_name_dict, dist):
        print('Excluding samples according to user supplied sample_name_compiled_re_excluded regular expression')
        exclude = []
        for sample_name in [sample_uid_to_sample_name_dict[sample_uid] for sample_uid in count_df.index]:
            if sample_name_compiled_re_excluded.match(sample_name):
                exclude.append(sample_name_to_sample_uid_dict[sample_name])
        if not exclude:
            raise RuntimeError('No sample names matched the user supplied sample_name_compiled_re_excluded.\n')
        print(f'Excluding {len(exclude)} samples')
        return self._drop_from_df(df=count_df, labels=exclude, dist=dist)

    def _exclude_samples_sample_name_compiled_re_included(self,
            count_df, sample_name_compiled_re_included,
            sample_name_to_sample_uid_dict, sample_uid_to_sample_name_dict, dist
    ):
        print('Excluding samples according to user supplied sample_name_compiled_re_included regular expression')
        keep = []
        for sample_name in [sample_uid_to_sample_name_dict[sample_uid] for sample_uid in count_df.index]:
            if sample_name_compiled_re_included.match(sample_name):
                keep.append(sample_name_to_sample_uid_dict[sample_name])
        if not keep:
            raise RuntimeError('No sample names matched the user supplied sample_name_compiled_re_included.\n')
        diff_set = set(keep).symmetric_difference(count_df.index.values)
        count_df = self._drop_from_df(df=count_df, labels=diff_set, dist=dist)
        print(f'Excluding {len(diff_set)} samples')
        return self._reindex_df(df=count_df, labels=keep, dist=dist)

    def _exclude_samples_sample_uids_excluded(self, count_df, sample_uids_excluded, dist):
        print('Excluding samples according to user supplied sample_uids_excluded list')
        # Check that all the uids are found in the seq_count df
        if not set(sample_uids_excluded).issubset(set(count_df.index.values)):
            diff_set = set(sample_uids_excluded).difference(set(count_df.index.values))
            raise RuntimeError('The following uids that were specified in the sample_uids_excluded list are not'
                               f'found in the SymPortal count table: {diff_set}')
        else:
            # All samples were found in the df and those not listed can be exculded
            print(f'Excluding {len(sample_uids_excluded)} samples')
            return self._drop_from_df(df=count_df, labels=sample_uids_excluded, dist=dist)

    def _exclude_samples_sample_uids_included(self, sample_uids_included, count_df, dist):
        # Return df containing only the included samples sorted according to that sampple order
        print('Excluding samples according to user supplied sample_uids_included list')
        # Check that all uids are found in the seq_count_df
        if not set(sample_uids_included).issubset(set(count_df.index.values)):
            diff_set = set(sample_uids_included).difference(set(count_df.index.values))
            raise RuntimeError('The following uids that were specified in the sample_uids_included list are not'
                               f'found in the SymPortal count table: {diff_set}')
        else:
            # All samples were found in the df and those not listed can be exculded
            diff_set = set(sample_uids_included).symmetric_difference(set(count_df.index.values))
            count_df = self._drop_from_df(df=count_df, labels=diff_set, dist=dist)
            print(f'Excluding {len(diff_set)} samples')
            return self._reindex_df(df=count_df, labels=sample_uids_included, dist=dist)

    @staticmethod
    def _drop_from_df(df, labels, dist):
        if dist:
            return df.drop(index=labels, columns=labels)
        else:
            return df.drop(index=labels)

    @staticmethod
    def _reindex_df(df, labels, dist):
        if dist:
            return df
        else:
            return df.reindex(labels, axis=0)
