# vim: fdm=indent
# author:     Fabio Zanini
# date:       17/06/19
# content:    Test the algorithm on same artificial data
import numpy as np
import pandas as pd
from northstar import Averages, AtlasFetcher


def test_constructor():
    sa = Averages(
            'Enge_2017',
            n_neighbors=2,
            threshold_neighborhood=0.8,
            n_pcs=2,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1
            )

    assert(sa is not None)


def test_arguments():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )

    sa = Averages(
            'Enge_2017',
            n_neighbors=2,
            threshold_neighborhood=0.8,
            n_pcs=2,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1
            )

    sa.new_data = matrix
    sa._check_init_arguments()
    assert(sa is not None)


def test_merge_small():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )
    n_pcs = 2
    k = 1
    threshold = 0.8

    sa = Averages(
            'Enge_2017',
            n_neighbors=k,
            threshold_neighborhood=threshold,
            n_pcs=n_pcs,
            distance_metric='correlation',
            n_neighbors_out_of_atlas=1
            )

    sa.new_data = matrix
    sa._check_init_arguments()
    sa.fetch_atlas_if_needed()
    sa.compute_feature_intersection()
    sa._check_feature_intersection()
    sa.prepare_feature_selection()
    sa.select_features()
    sa._check_feature_selection()
    sa.merge_atlas_newdata()

    assert(sa is not None)


def test_neighbors_small():
    matrix = pd.DataFrame(
        index=['INS', 'GCG', 'PPY'],
        columns=['cell1', 'cell2', 'cell3'],
        data=[
            [2302, 123, 0],
            [0, 5034, 6453],
            [0, 0, 1]],
        )
    n_pcs = 2
    k = 1
    threshold = 0.8

    sa = Averages(
            'Enge_2017',
            n_neighbors=k,
            threshold_neighborhood=threshold,
            n_pcs=n_pcs,
            distance_metric='correlation',
            n_cells_per_type=20,
            n_neighbors_out_of_atlas=1
            )

    sa.new_data = matrix
    sa._check_init_arguments()
    sa.fetch_atlas_if_needed()
    sa.compute_feature_intersection()
    sa._check_feature_intersection()
    sa.prepare_feature_selection()
    sa.select_features()
    sa._check_feature_selection()
    sa.merge_atlas_newdata()
    sa.compute_pca()
    sa.compute_similarity_graph()
    neis = sa.neighbors

    assert(isinstance(neis, list))
    assert(len(neis) == int(np.sum(sa.sizes)))
    for nei in neis:
        assert(isinstance(nei, list))


def test_run_mock_cells():
    aname = 'Darmanis_2015_nofetal'
    atlas = AtlasFetcher().fetch_atlas(
            aname, kind='average')

    ind = [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 4, 5]
    matrix = atlas[ind]
    cell_types = atlas.obs_names.values[ind]

    sa = Averages(
            aname,
            n_pcs=10,
            )
    sa.fit(matrix)

    assert((cell_types == sa.membership).mean() > 0.9)


def test_run_within_atlas():
    aname = 'Darmanis_2015_nofetal'
    atlas = AtlasFetcher().fetch_atlas(
            aname, kind='subsample')

    ind = [
        0, 2, 5, 8, 10, 15, 20, 25, 28, 30, 35, 38,
        40, 45, 50, 60, 70, 75, 80, 90]
    matrix = atlas[ind]
    cell_types = atlas.obs['CellType'].values[ind]

    sa = Averages(
            aname,
            n_features_per_cell_type=2,
            n_features_overdispersed=5,
            n_pcs=9,
            )
    sa.fit(matrix)

    for i in range(len(cell_types)):
        print('{:10s}    {:10s}'.format(cell_types[i], sa.membership[i]))
    print((cell_types == sa.membership).mean())
    assert((cell_types == sa.membership).mean() >= 0.5)


if __name__ == '__main__':

    test_run_within_atlas()
