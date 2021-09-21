from Simulation import generate_simulated_data
from utils import ProcessInputData
from train import train_model, predict
from anndata import read_h5ad
from torch.utils.data import DataLoader
from model import reproducibility, simdatset, AutoEncoder, device

def OverallDeconvolution(sc_reference, real_bulk, bulkDataType,
                         sep = '\t',
                         save_simu = False, save_model = False,
                         mode='single', adaptive=True):
    """
    :param sc_reference: a txt expression file path index is cell type name, columns is gene name
    :param real_bulk: an expression file path, index is sample, columns is gene name
    :param bulkDataType: FPKM or TPM, if type is RPKM, please just use FPKM
    :param sep: used to read bulk data, depends on the format
    :param save_simu: save simulated training data or not
    :param save_model: save model or not
    :param mode: 'single' means this will apply adaptive stage to every single sample to generate signature matrix
    :param adaptive: it has to be True, if model is 'single'
    :return: depends on the mode or adaptive
    """
    simudata = generate_simulated_data(sc_data=sc_reference, samplenum=5000)
    train_x, train_y, test_x, genename, celltypes, samplename = \
        ProcessInputData(simudata, real_bulk, sep=sep, datatype=bulkDataType)
    print(train_x.shape,test_x.shape)
    if save_model is True:
        model = train_model(train_x,train_y,'model',batch_size=128,iteration=5000)
    else:
        model = train_model(train_x,train_y,batch_size=128,iteration=5000)
    if adaptive is True:
        if mode == 'single':
            CellTypeSigm, Pred = \
                predict(test_x=test_x,genename=genename,celltypes=celltypes,samplename=samplename,
                        model=model,
                        adaptive=adaptive,mode=mode)
            return CellTypeSigm, Pred
        elif mode == 'all':
            Sigm, Pred = \
                predict(test_x=test_x, genename=genename, celltypes=celltypes, samplename=samplename,
                        model=model,
                        adaptive=adaptive,mode=mode)
            return Sigm, Pred
    else:
        Pred = predict(test_x=test_x, genename=genename, celltypes=celltypes, samplename=samplename,
                       model=model,
                       adaptive=adaptive,mode=mode)
        return Pred




