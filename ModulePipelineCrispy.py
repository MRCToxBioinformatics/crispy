from cgatcore import pipeline as P
import cgatcore.iotools as iotools
import os
from cgatcore.pipeline import cluster_runnable

import numpy as np
import pandas as pd
from collections import Counter

@cluster_runnable
def resampleTallies(infile, sample_frac, outfile):
    '''
    infile = tab-seperated table of gRNA counts. Assumes columns 0 & 1 are index
    sample_frac = resampling fraction. Can be >1
    outfile = filename to save resampled counts
    '''

    counts = pd.read_csv(infile, sep='\t', index_col=[0, 1])

    # integer values for counts.index
    choice_index = range(0, len(counts))

    new_cols = []
    for column in range(0, len(counts.columns)):

        sample_counts = counts.iloc[:,column]

        choice_size = round((sample_counts.sum())*sample_frac)

        choice_p = sample_counts/sample_counts.sum()

        test_subsample = Counter(np.random.choice(a = choice_index,
                                                  size = choice_size,
                                                  p=choice_p))
        new_cols.append(dict(test_subsample))

    df = pd.DataFrame.from_dict(new_cols).transpose()
    df.columns = counts.columns

    # re-order dataframe to be the same as the input, allowing df.index to be a
    # subset of counts.index
    df = df.loc[[x for x in choice_index if x in df.index],]

    # rename df.index to be use the index columns in input
    df.index = counts.index[df.index]

    # NAs are zero counts from re-sampling
    df.fillna(0, inplace=True)

    df.to_csv(outfile, sep='\t')
