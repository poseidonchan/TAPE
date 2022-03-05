# TAPE: Tissue-AdaPtive autoEncoder for accurate deconvolution and gene expression analysis

![scTAPE](https://img.shields.io/badge/scTAPE-v0.2.6-blue)![GPL](https://img.shields.io/github/license/poseidonchan/TAPE)[![Downloads](https://static.pepy.tech/personalized-badge/sctape?period=total&units=international_system&left_color=grey&right_color=brightgreen&left_text=downloads)](https://pepy.tech/project/sctape)

**This model is able to accurately deconvolve bulk RNA-seq data into cell fractions and predict cell-type-specific gene expression at cell-type level based on scRNA-seq data**.

related article ***Deep autoencoder for interpretable tissue-adaptive deconvolution and cell-type-specific gene analysis*** is available on [*bioRxiv*](https://doi.org/10.1101/2021.10.26.465846)

## Setup
TAPE uses PyTorch as its Deep-learning framework, so the suitable version of PyTorch will accelerate the model training process. We recommend users to install PyTorch(>=1.8.0) with ***right*** compute platform (CUDA, CPU or ROCm) from its official [website](https://pytorch.org) in advance.

For example, we used NVIDIA GPU RTX 3090, so we choose the CUDA version 11.1 and the command is:

```bash
pip install torch==1.8.0+cu111 torchvision==0.9.0+cu111 torchaudio==0.8.0 -f https://download.pytorch.org/whl/torch_stable.html
```

If PyTorch is successfully installed, then TAPE could be installed from PyPI directly:

```bash
pip install scTAPE
```
## Usage
Required Files:
1. single-cell reference: txt format, indices are cell types, columns are gene names
2. bulk data: tabular format, needed to specify the seperation ('\t',','or others), indices are sample names, columns are gene names
3. gene length file: used to scale the expression value, columns should contain: [Gene name, Transcript start (bp), Transcript end (bp)]. This is provided in ./data/ directory.
```python
# basic example
from TAPE import Deconvolution
SignatureMatrix, CellFractionPrediction = \
    Deconvolution(sc_ref, bulkdata, sep='\t',
                  datatype='counts', genelenfile='./GeneLength.txt',
                  mode='overall', adaptive=True,
                  save_model_name=None,
                  batch_size=128, epochs=128)
```
parameters:

1. datatype: use '**TPM**', '**FPKM**' or '**counts**'. Users can choose different normalization method based on your single-cell seq technique, if single-cell data is from 10X Genomics, users should use '**counts**' to maintain a resonable procedure. The explanation could be found from the [webpage](https://kb.10xgenomics.com/hc/en-us/articles/115003684783-Should-I-calculate-TPM-RPKM-or-FPKM-instead-of-counts-for-10x-Genomics-data-).
2. mode: '**overall**' or '**high-resolution**'. If you need signature matrix for each sample, use 'high-resolution' mode.
3. adaptive: **True** or **False**. If this is False, then it would not predict signature matrix, the return will be ***None***
4. batch_size: **int**, related to training result. 32-128 are recommended. Smaller batch_size leads to more time consumption.
5. epochs: **int**, related to training result. Typically, *5000-10000* iterations are enough for TAPE, the relation is $epochs=\frac{iteration \times batch\_size}{sampleing\_num}$

Since the original implementation of Scaden [[repository](https://github.com/KevinMenden/scaden)] [[paper](https://www.science.org/doi/10.1126/sciadv.aba2619)] is not easy for us to test, we implemented the PyTorch version of Scaden. If you want to use Scaden to deconvolve bulk RNA-seq data, you can use the following code:

```python
from TAPE.deconvolution import ScadenDeconvolution
Pred = ScadenDeconvolution(sc_ref, bulkdata, sep='\t',
                           batch_size=128, epochs=128)
```

The parameters in this function are the same as the parameters in the *Deconvolution* function. 


## Example
An example is placed in the **Test** directory. Please run the example to get familiar with TAPE.

## Issues
If you find any bugs or have problems when you are using scTAPE, feel free to raise issues.

## Citation
If you feel TAPE is useful and want to cite TAPE, please use the following format:
[Bib TeX](https://cloud.tsinghua.edu.cn/f/6351d49841934fa79bc5/?dl=1)
[EndNote](https://cloud.tsinghua.edu.cn/f/88780d85281c40c0848f/?dl=1)

## Acknowledgement
Special thanks to [*Mengyue Sun*](https://github.com/sunmy2019), for his help to accelerate the sampling process (in the simulation.py).

Much thanks to [*Yibo Liu*](https://github.com/jedibobo), for his advice on building such a nice repository.