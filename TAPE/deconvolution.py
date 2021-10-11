from Simulation import generate_simulated_data
from utils import ProcessInputData
from train import train_model, predict
from anndata import read_h5ad
from torch.utils.data import DataLoader
from model import reproducibility, simdatset, AutoEncoder, device

def OverallDeconvolution(sc_reference, real_bulk, bulkDataType,
                         sep = '\t',
                         mode='overall', adaptive=True,
                         save_model_name = None,):
    """
    :param sc_reference: a txt expression file path index is cell type name, columns is gene name
    :param real_bulk: an expression file path, index is sample, columns is gene name
    :param bulkDataType: FPKM or TPM, if type is RPKM, please just use FPKM
    :param sep: used to read bulk data, depends on the format
    :param mode: 'high-resolution' means this will apply adaptive stage to every single sample to generate signature matrix,
                 'overall' means that it will deconvolve all the samples at the same time
    :param adaptive: it has to be True, if model is 'high-resolution'
    :param save_model_name: the name used to save model, if it was not provided, it would not be saved
    :return: depends on the mode or adaptive
             there are three combinations:
             1. high-resolution and adaptive deconvolution
                this will return a dictionary and predicted cell fractions in pandas dataframe format
                the keys of the dict are the pre-defined cell type names in the single cell reference data
                the values of the dict are the dataframe of gene expression and samples
             2. overall and adaptive deconvolution
                this will return a signature matrix and a cell fraction
                the rows of the signature matrix is the gene expression in each cell types
                both of the variables are in dataframe format
             3. overall and non-adaptive deconvolution
                this will return a cell fraction directly
                the signature matrix in this mode is None
    """
    simudata = generate_simulated_data(sc_data=sc_reference, samplenum=5000)
    train_x, train_y, test_x, genename, celltypes, samplename = \
        ProcessInputData(simudata, real_bulk, sep=sep, datatype=bulkDataType)
    print('training data shape is ',train_x.shape,'\ntest data shape is ',test_x.shape)
    if save_model_name is not None:
        model = train_model(train_x,train_y,save_model_name,batch_size=128,iteration=5000)
    else:
        model = train_model(train_x,train_y,batch_size=128,iteration=5000)

    print('Notice that you are using parameters: mode='+str(mode)+' and adaptive='+str(adaptive))
    if adaptive is True:
        if mode == 'high-resolution':

            CellTypeSigm, Pred = \
                predict(test_x=test_x,genename=genename,celltypes=celltypes,samplename=samplename,
                        model=model, model_name=save_model_name,
                        adaptive=adaptive,mode=mode)
            return CellTypeSigm, Pred
        elif mode == 'overall':
            Sigm, Pred = \
                predict(test_x=test_x, genename=genename, celltypes=celltypes, samplename=samplename,
                        model=model, model_name=save_model_name,
                        adaptive=adaptive,mode=mode)
            return Sigm, Pred
    else:
        Pred = predict(test_x=test_x, genename=genename, celltypes=celltypes, samplename=samplename,
                       model=model, model_name=save_model_name,
                       adaptive=adaptive,mode=mode)
        Sigm = None
        return Sigm, Pred




