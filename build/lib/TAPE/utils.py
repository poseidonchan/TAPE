import os
import anndata
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
#### NEEDED FILES
# 1. GeneLength.txt


def counts2FPKM(counts,genelen):
    genelen = pd.read_csv(genelen, sep=',')
    genelen['TranscriptLength'] = genelen['Transcript end (bp)'] - genelen['Transcript start (bp)']
    genelen = genelen[['Gene name', 'TranscriptLength']]
    genelen = genelen.groupby('Gene name').max()
    # intersection
    inter = counts.columns.intersection(genelen.index)
    samplename = counts.index
    counts = counts[inter].values
    genelen = genelen.loc[inter].T.values
    # transformation
    totalreads = counts.sum(axis=1)
    counts = counts*1e9/(genelen*totalreads.reshape(-1,1))
    counts = pd.DataFrame(counts,columns=inter,index=samplename)
    return counts

def FPKM2TPM(fpkm):
    genename = fpkm.columns
    samplename = fpkm.index
    fpkm = fpkm.values
    total = fpkm.sum(axis=1).reshape(-1,1)
    fpkm = fpkm*1e6/total
    fpkm = pd.DataFrame(fpkm,columns=genename,index=samplename)
    return fpkm

def counts2TPM(counts,genelen):
    fpkm = counts2FPKM(counts,genelen)
    tpm = FPKM2TPM(fpkm)
    return tpm


def ProcessInputData(train_data, test_data, sep=None, datatype='TPM', variance_threshold=0.5,
                     genelenfile=None):

    ### read train data
    print('Reading training data')
    if type(train_data) is anndata.AnnData:
        pass
    elif type(train_data) is str:
        train_data = anndata.read_h5ad(train_data)
    train_data.var_names_make_unique()
    train_x = pd.DataFrame(train_data.X, columns=train_data.var.index)
    train_y = train_data.obs
    print('Reading is done')
    ### read test data
    print('Reading test data')
    if type(test_data) is str:
        test_x = pd.read_csv(test_data, index_col=0, sep=sep)
    elif type(test_data) is pd.DataFrame:
        test_x = test_data
    print('Reading test data is done')
    ### transform to datatype
    if datatype == 'FPKM':
        if genelenfile is None:
            raise Exception("Please add gene length file!")
        print('Transforming to FPKM')
        train_x = counts2FPKM(train_x,genelenfile)
    elif datatype == 'TPM':
        if genelenfile is None:
            raise Exception("Please add gene length file!")
        print('Transforming to TPM')
        train_x = counts2TPM(train_x,genelenfile)
    elif datatype == 'counts':
        print('Using counts data to train model')
    ### variance cutoff
    print('Variance Cutoff')
    test_x = test_x.loc[:, test_x.var(axis=0) > variance_threshold]

    ### find intersected genes
    print('Find intersected genes')
    inter = train_x.columns.intersection(test_x.columns)
    train_x = train_x[inter]
    test_x = test_x[inter]
    genename = list(inter)
    celltypes = train_y.columns
    samplename = test_x.index
    
    print('Intersected gene number is ',len(inter))
    ### MinMax process
    print('Log2 & MinMax scale')
    train_x = np.log2(train_x.values + 1)
    test_x = np.log2(test_x.values + 1)
    mms = MinMaxScaler()
    train_x = mms.fit_transform(train_x.T).T
    test_x = mms.fit_transform(test_x.T).T

    return train_x, train_y.values, test_x, genename, celltypes, samplename


def L1error(pred, true):
    return np.mean(np.abs(pred - true))

def CCCscore(y_pred, y_true, mode='all'):
    # pred: shape{n sample, m cell}
    if mode == 'all':
        y_pred = y_pred.reshape(-1,1)
        y_true = y_true.reshape(-1,1)
    elif mode == 'avg':
        pass
    ccc_value = 0
    for i in range(y_pred.shape[1]):
        r = np.corrcoef(y_pred[:, i], y_true[:, i])[0, 1]
        # Mean
        mean_true = np.mean(y_true[:, i])
        mean_pred = np.mean(y_pred[:, i])
        # Variance
        var_true = np.var(y_true[:, i])
        var_pred = np.var(y_pred[:, i])
        # Standard deviation
        sd_true = np.std(y_true[:, i])
        sd_pred = np.std(y_pred[:, i])
        # Calculate CCC
        numerator = 2 * r * sd_true * sd_pred
        denominator = var_true + var_pred + (mean_true - mean_pred) ** 2
        ccc = numerator / denominator
        ccc_value += ccc
    return ccc_value / y_pred.shape[1]

def score(pred, label):
    print('L1 error is', L1error(pred,label))
    print('CCC is ', CCCscore(pred, label))

def showloss(loss):
    sns.set()
    plt.plot(loss)
    plt.xlabel('iteration')
    plt.ylabel('loss')
    plt.show()



