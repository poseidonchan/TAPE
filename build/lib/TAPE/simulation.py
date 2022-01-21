import anndata
import numpy as np
import pandas as pd
from tqdm import tqdm
from numpy.random import choice

def generate_simulated_data(sc_data, outname=None,
                            d_prior=None,
                            n=500, samplenum=5000):
    # sc_data should be a cell*gene matrix, no null value, txt file, sep='\t'
    # index should be cell names
    # columns should be gene labels
    print('Reading single-cell dataset, this may take 1 min')
    if type(sc_data) is str:
        sc_data = pd.read_csv(sc_data,index_col=0,sep='\t')
    elif type(sc_data) is pd.DataFrame:
        pass
    else:
        raise Exception("Please check the format of single-cell data!")
    print('Reading dataset is done')
    sc_data.dropna(inplace=True)
    sc_data['celltype'] = sc_data.index
    sc_data.index = range(len(sc_data))

    num_celltype = len(sc_data['celltype'].value_counts())
    genename = sc_data.columns[:-1]
    celltype_groups = sc_data.groupby('celltype').groups
    sc_data.drop(columns='celltype' ,inplace=True)

    # use ndarray to accelerate
    # change to C_CONTIGUOUS, 10x speed up
    sc_data = sc_data.values
    sc_data = np.ascontiguousarray(sc_data, dtype=np.float32)
    # make random cell proportions

    if d_prior is None:
        print('Generating cell fractions using Dirichlet distribution without prior info (actually random)')
        prop = np.random.dirichlet(np.ones(num_celltype),samplenum)
        print('RANDOM cell fractions is generated')
    elif d_prior is not None:
        print('Using prior info to generate cell fractions in Dirichlet distribution')
        assert len(d_prior) == num_celltype, 'dirichlet prior is a vector, its length should equals ' \
                                             'to the number of cell types'
        prop = np.random.dirichlet(d_prior, samplenum)
        print('Dirichlet cell fractions is generated')

    # make the dictionary
    for key, value in celltype_groups.items():
        celltype_groups[key] = np.array(value)
    # precise number for each celltype
    cell_num = np.floor(n*prop)

    # start sampling
    sample = np.zeros((prop.shape[0] ,sc_data.shape[1]))
    allcellname = celltype_groups.keys()
    print('Sampling cells to compose pseudo-bulk data')
    for i ,sample_prop in tqdm(enumerate(cell_num)):
        for j, cellname in enumerate(allcellname):
            select_index = choice(celltype_groups[cellname] ,size=int(sample_prop[j]) ,replace=True)
            sample[i] += sc_data[select_index].sum(axis=0)
            
    prop = pd.DataFrame(prop ,columns=celltype_groups.keys())
    simulated_dataset = anndata.AnnData(X=sample,
                                        obs=prop,
                                        var=pd.DataFrame(index=genename))
    print('Sampling is done')
    if outname is not None:
        simulated_dataset.write_h5ad(outname+'.h5ad')
    return simulated_dataset