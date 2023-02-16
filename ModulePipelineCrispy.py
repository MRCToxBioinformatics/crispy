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


@cluster_runnable
def normaliseCounts(prenorm_table_filepath, norm_table_filepath):

    # mageckcount_getmediannormfactor is unmodified copy from
    # https://bitbucket.org/liulab/mageck/src/5a8503bbd864c9785f3d139e4d6317404b986cdd/mageck/mageckCountNorm.py#lines-15:28
    #
    # Copyright (c) 2013, Xiaole Shirley Liu lab at DFCI
    # All rights reserved.

    # Redistribution and use in source and binary forms, with or without
    # modification, are permitted provided that the following conditions are met:
    #     * Redistributions of source code must retain the above copyright
    #       notice, this list of conditions and the following disclaimer.
    #     * Redistributions in binary form must reproduce the above copyright
    #       notice, this list of conditions and the following disclaimer in the
    #       documentation and/or other materials provided with the distribution.
    #     * Neither the name of Xiaole Shirley Liu lab nor the
    #       names of its contributors may be used to endorse or promote products
    #       derived from this software without specific prior written permission.

    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    # "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    # LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    # A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT
    # HOLDERs AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
    # INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    # BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
    # OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
    # TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
    # USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
    # DAMAGE.

    def mageckcount_getmediannormfactor(ctable):
      """
      Get the factor by median normalization
      """
      n=len(ctable[list(ctable.keys())[0]]) # samples
      m=len(ctable) # sgRNAs
      meanval={k:math.exp( (sum( [ math.log(v2+1.0) for v2 in v])*1.0/n) ) for (k,v) in ctable.items() if sum(v)>0}  # geometric mean
      meanval={k:(lambda x: x if x>0 else 1)(v) for (k,v) in meanval.items()} # delete those with all 0 read counts
      #samplefactor=[0]*n
      medianfactor=[0.0]*n
      for ni in range(n):
        meanfactor=[ v[ni]/meanval[k] for (k,v) in ctable.items() if k in meanval]
        #print(str(sorted(meanfactor)))
        xfactor=sorted(meanfactor)[len(meanfactor)//2] # corrected
        if xfactor>0.0:
          medianfactor[ni]=1.0/xfactor
          #logging.debug('xfactor:'+str(xfactor))
      return medianfactor


    counts = pd.read_csv(prenorm_table_filepath, sep='\t')

    ctable = {}
    for colname in counts_table.columns[2:]:
        ctable[colname] = counts_table[colname]


    medianfactor = mageckcount_getmediannormfactor(ctable)

    samplefactor = medianfactor
    ntable={ k: [ samplefactor[i]*v[i] for i in range(n)] for (k,v) in ctable.items()}

    norm_counts_table = counts_table.copy()
    for k,v in ntable.items():
        norm_counts_table[k] = v

    norm_counts_table.to_csv(norm_table_filepath, sep='\t')