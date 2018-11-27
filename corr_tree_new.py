class CorrTree(object):
    def __init__(self, feature_table, tree, metadata):
        self.feature_table = biom.load_table(feature_table)
        self.tree = tree
        self.metadata = metadata
    def compute_no(self, obj_col):
        """Compute correlation for every node.

        Args:
            obj_col: object column whose values are correlated to every node.(df.series)

        Return:

        """
        try:
            tmp = float(obj_col[0])
        except:
            print('it is not a value array')

    def compute_stats_test(self, obj_col):
        pass

        