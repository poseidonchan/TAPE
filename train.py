import numpy as np
import pandas as pd
from tqdm import tqdm
import torch
from torch.optim import Adam
from torch.utils.data import DataLoader
from model import reproducibility, simdatset, AutoEncoder, device
from utils import showloss
import torch.nn.functional as F

def training_stage(model, train_loader, optimizer, epochs=10):
    model.train()
    model.state = 'train'
    loss = []
    recon_loss = []
    for i in tqdm(range(epochs)):
        for k, (data, label) in enumerate(train_loader):
            optimizer.zero_grad()
            x_recon, cell_prop, sigm = model(data)
            batch_loss = F.l1_loss(cell_prop, label) + F.l1_loss(x_recon, data)
            batch_loss.backward()
            optimizer.step()
            loss.append(F.l1_loss(cell_prop, label).cpu().detach().numpy())
            recon_loss.append(F.l1_loss(x_recon, data).cpu().detach().numpy())

    return model, loss, recon_loss

def adaptive_stage(model, data, optimizerD, optimizerE, step=10, max_iter=5):
    data = torch.from_numpy(data).float().to(device)
    loss = []
    model.eval()
    model.state = 'train'
    ori_pred = model.encode(data).detach()
    ori_sigm = model.sigmatrix().detach()
    for k in range(max_iter):
        model.eval()
        for i in range(step):
            optimizerD.zero_grad()
            x_recon, _, sigm = model(data)
            batch_loss = F.l1_loss(x_recon, data)+F.l1_loss(ori_sigm,sigm)
            batch_loss.backward()
            optimizerD.step()
            loss.append(F.l1_loss(x_recon, data).cpu().detach().numpy())

        for i in range(step):
            optimizerE.zero_grad()
            x_recon, pred, _ = model(data)
            batch_loss = F.l1_loss(x_recon, data) + F.l1_loss(ori_pred, pred)
            batch_loss.backward()
            optimizerE.step()
            loss.append(F.l1_loss(x_recon, data).cpu().detach().numpy())

    model.eval()
    model.state = 'test'
    _, pred, sigm = model(data)
    print(model.state, pred.sum(dim=1))
    return sigm.cpu().detach().numpy(), loss, pred.cpu().detach().numpy()

def train_model(train_x, train_y,
                model_name=None,
                batch_size=128, iteration=10000):

    
    reproducibility(9)
    
    train_loader = DataLoader(simdatset(train_x, train_y), batch_size=batch_size, shuffle=True)
    model = AutoEncoder(train_x.shape[1], train_y.shape[1]).to(device)
    optimizer = Adam(model.parameters(), lr=1e-4)
    model, loss, reconloss = training_stage(model, train_loader, optimizer, epochs=int(iteration /(len(train_x)/128)))
    print('prediction loss is:')
    showloss(loss)
    print('reconstruction loss is:')
    showloss(reconloss)
    if model_name is not None:
        torch.save(model, model_name+".pth")
    return model



# a = np.random.randn(3,10)
# print(a)
# mms = MinMaxScaler()
# a = mms.fit_transform(a.T).T
# print(a)
# data_min = mms.data_max_.reshape(-1,1)
# data_max = mms.data_min_.reshape(-1,1)
# scale = np.concatenate((data_max,data_min),axis=1)
# for i in range(a.shape[0]):
#     a[i,:] = a[i,:]*(scale[i,1]-scale[i,0])+scale[i,0]
# print(a)

def predict(test_x, genename, celltypes, samplename,
            model_name=None, model=None,
            adaptive=True, mode='single'):
    reproducibility(9)
    
    if model is not None and model_name is None:
        torch.save(model, 'model.pth')
    if adaptive is True:
        if mode == 'single':
            TestSigmList = np.zeros((test_x.shape[0], len(celltypes), len(genename)))
            TestPred = np.zeros((test_x.shape[0], len(celltypes)))
            
            """
            scale the data
            """
#             mms = MinMaxScaler()
#             test_x = mms.fit_transform(test_x.T).T
#             data_min = mms.data_max_.reshape(-1,1)
#             data_max = mms.data_min_.reshape(-1,1)
#             scale = np.concatenate((data_max,data_min),axis=1)
            
            
            for i in tqdm(range(len(test_x))):
                x = test_x[i,:].reshape(1,-1)
                if model_name is not None and model is None:
                    model = torch.load(model_name + ".pth")
                elif model is not None and model_name is None:
                    model = torch.load("model.pth")
                decoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'decoder' in n]}]
                encoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'encoder' in n]}]
                optimizerD = torch.optim.Adam(decoder_parameters, lr=1e-4)
                optimizerE = torch.optim.Adam(encoder_parameters, lr=1e-4)
                test_sigm, loss, test_pred = adaptive_stage(model, x, optimizerD, optimizerE, step=250, max_iter=3)
                showloss(loss)
                TestSigmList[i, :, :] = test_sigm
                TestPred[i,:] = test_pred
            TestPred = pd.DataFrame(TestPred,columns=celltypes,index=samplename)
            CellTypeSigm = {}
            for i in range(len(celltypes)):
                cellname = celltypes[i]
                sigm = TestSigmList[:,i,:]
                
                """
                restore the data
                """
#                 for i in range(sigm.shape[0]):
#                     sigm[i,:] = sigm[i,:]*(scale[i,1]-scale[i,0])+scale[i,0]


                sigm = pd.DataFrame(sigm,columns=genename,index=samplename)
                CellTypeSigm[cellname] = sigm

            return CellTypeSigm, TestPred

        elif mode == 'all':
            if model_name is not None and model is None:
                model = torch.load(model_name + ".pth")
            elif model is not None and model_name is None:
                model = torch.load("model.pth")
            decoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'decoder' in n]}]
            encoder_parameters = [{'params': [p for n, p in model.named_parameters() if 'encoder' in n]}]
            optimizerD = torch.optim.Adam(decoder_parameters, lr=1e-4)
            optimizerE = torch.optim.Adam(encoder_parameters, lr=1e-4)
            test_sigm, loss, test_pred = adaptive_stage(model, test_x, optimizerD, optimizerE, step=250, max_iter=3)
            test_sigm = pd.DataFrame(test_sigm,columns=genename,index=celltypes)
            test_pred = pd.DataFrame(test_pred,columns=celltypes,index=samplename)
            return test_sigm, test_pred

    else:
        if model_name is not None and model is None:
            model = torch.load(model_name+".pth")
        elif model is not None and model_name is None:
            model = model
        model.eval()
        data = torch.from_numpy(test_x).float().to(device)
        _, pred, _ = model(data)
        pred = pred.cpu().detach().numpy()
        pred = pd.DataFrame(pred, columns=celltypes, index=samplename)
        return pred



