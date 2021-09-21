import anndata
import numpy as np
import pandas as pd
from tqdm import tqdm

def generate_simulated_data(sc_data, outname=None,
                            prop = None,
                            n=500, samplenum=4000):
    # sc_data should be a cell*gene matrix, no null value, txt file, sep='\t'
    # index should be cell names
    # columns should be gene labels

    sc_data = pd.read_csv(sc_data,index_col=0,sep='\t')
    sc_data.dropna(inplace=True)
    sc_data['celltype'] = sc_data.index
    sc_data.index = range(len(sc_data))

    num_celltype = len(sc_data['celltype'].value_counts())
    genename = sc_data.columns[:-1]
    celltype_groups = sc_data.groupby('celltype').groups
    sc_data.drop(columns='celltype' ,inplace=True)

    # use ndarray to accelerate
    sc_data = sc_data.values

    # make random cell proportions
    if prop is None:
        prop = np.random.rand(samplenum ,num_celltype)
        prop = prop /np.sum(prop ,axis=1).reshape(-1 ,1)

        for i in range(prop.shape[0]):
            indices = np.random.choice(np.arange(prop.shape[1]), replace=False, size=int(prop.shape[1] * 0.3))
            prop[i ,indices] = 0
        prop = prop /np.sum(prop ,axis=1).reshape(-1 ,1)

    # make the dictionary
    for key, value in celltype_groups.items():
        celltype_groups[key] = np.array(value)
    # precise number for each celltype
    cell_num = np.floor(n*prop)

    # start sampling
    sample = np.zeros((prop.shape[0] ,sc_data.shape[1]))
    for i ,sample_prop in tqdm(enumerate(cell_num)):
        for j, cellname in enumerate(celltype_groups.keys()):
            select_index = np.random.choice(celltype_groups[cellname] ,size=int(sample_prop[j]) ,replace=True)
            sample[i] = sample[i] + sc_data[select_index].sum(axis=0)

    prop = pd.DataFrame(prop ,columns=celltype_groups.keys())
    simulated_dataset = anndata.AnnData(X=sample,
                                        obs=prop,
                                        var=pd.DataFrame(index=genename))
    if outname is not None:
        simulated_dataset.write_h5ad(outname+'.h5ad')
    return simulated_dataset