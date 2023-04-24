# TAPE: Tissue-AdaPtive autoEncoder for accurate deconvolution and gene expression analysis

![scTAPE](https://img.shields.io/badge/scTAPE-v1.1.2-blue)![GPL](https://img.shields.io/github/license/poseidonchan/TAPE)[![DOI](https://zenodo.org/badge/386163911.svg)](https://zenodo.org/badge/latestdoi/386163911)

**This model is able to accurately deconvolve bulk RNA-seq data into cell fractions and predict cell-type-specific gene expression at cell-type level based on scRNA-seq data**.

related article ***Deep autoencoder for interpretable tissue-adaptive deconvolution and cell-type-specific gene analysis*** is accepted by *Nature Communications*

## Setup

TAPE uses PyTorch as its Deep-learning framework, so the suitable version of PyTorch will accelerate the model training process. We recommend users to install PyTorch(>=1.8.0) with ***right*** compute platform (CUDA, CPU or ROCm) from its official [website](https://pytorch.org) in advance.

For example, we used NVIDIA GPU RTX 3090, so we choose the CUDA version 11.1 and the command is:

```bash
pip install torch==1.8.0+cu111 torchvision==0.9.0+cu111 torchaudio==0.8.0 -f https://download.pytorch.org/whl/torch_stable.html
```

If PyTorch is successfully installed, then TAPE could be installed from PyPI directly:

***Update: I relax the dependece requirements to enable the compatibility with current packages***.

```bash
pip install scTAPE==1.1.2
```
Usually, the installation time depends on your downloading speed.
## Usage

Required Files:
1. single-cell reference: txt format, indices are cell types, columns are gene names
2. bulk data: tabular format, needed to specify the seperation ('\t',','or others), indices are sample names, columns are gene names
3. gene length file: used to scale the expression value, columns should contain: [Gene name, Transcript start (bp), Transcript end (bp)]. This is provided in ./data/ directory.

***Warning: single-cell reference and bulk samples should contain the same cell types***

```python
# basic example
from TAPE import Deconvolution
SignatureMatrix, CellFractionPrediction = \
    Deconvolution(sc_ref, bulkdata, sep='\t', scaler='mms',
                  datatype='counts', genelenfile='./GeneLength.txt',
                  mode='overall', adaptive=True, variance_threshold=0.98,
                  save_model_name=None,
                  batch_size=128, epochs=128, seed=1)
```
parameters:

1. scaler: use '**mms**' or '**ss**' scaler to preprocess datasets, 'mms' stands for min-max scaler, 'ss' stands for standard scaler. In the paper, all datasets were tested using 'mms'.
2. datatype: use '**TPM**', '**FPKM**' or '**counts**'. Users can choose different normalization method based on your single-cell seq technique, if single-cell data is from 10X Genomics, users should use '**counts**' to maintain a resonable procedure. The explanation could be found from the [webpage](https://kb.10xgenomics.com/hc/en-us/articles/115003684783-Should-I-calculate-TPM-RPKM-or-FPKM-instead-of-counts-for-10x-Genomics-data-).
3. mode: '**overall**' or '**high-resolution**'. If you need signature matrix for each sample, use 'high-resolution' mode.
4. adaptive: **True** or **False**. If this is False, then it would not predict signature matrix, the return will be ***None***
5. variance_threshold: Float number from 0 to 1, it means how many genes you want to keep (in proportion) according to variance from high to low.
6. batch_size: **int**, related to training result. 32-128 are recommended. Smaller batch_size leads to more time consumption.
7. epochs: **int**, related to training result. Typically, *5000-10000* iterations are enough for TAPE, the relation is $epochs=\frac{iteration \times batch\_size}{sampleing\_num}$
8. seed: now, TAPE supports pinning the random seed to make results being reproducible.

Since the original implementation of Scaden [[repository](https://github.com/KevinMenden/scaden)] [[paper](https://www.science.org/doi/10.1126/sciadv.aba2619)] is not easy for us to test, we implemented the PyTorch version of Scaden. If you want to use Scaden to deconvolve bulk RNA-seq data, you can use the following code:

```python
from TAPE.deconvolution import ScadenDeconvolution
Pred = ScadenDeconvolution(sc_ref, bulkdata, sep='\t',
                           batch_size=128, epochs=128)
```


## Example
An example is placed in the **Test** directory. Please run the example to get familiar with TAPE.

Run the demo may takes 2 to 3 mins with GPU acceleration or 10 mins with CPU.


## Issues
If you find any bugs or have problems when you are using scTAPE, feel free to raise issues.

## Citation
```bibtex
@article{TAPE,
   author = {Chen, Yanshuo and Wang, Yixuan and Chen, Yuelong and Cheng, Yuqi and Wei, Yumeng and Li, Yunxiang and Wang, Jiuming and Wei, Yingying and Chan, Ting-Fung and Li, Yu},
   title = {Deep autoencoder for interpretable tissue-adaptive deconvolution and cell-type-specific gene analysis},
   journal = {Nature Communications},
   volume = {13},
   number = {1},
   pages = {6735},
   ISSN = {2041-1723},
   DOI = {10.1038/s41467-022-34550-9},
   url = {https://doi.org/10.1038/s41467-022-34550-9},
   year = {2022},
   type = {Journal Article}
}
```

## Acknowledgement
Special thanks to [*Mengyue Sun*](https://github.com/sunmy2019), for his help to accelerate the sampling process (in the simulation.py).

Much thanks to [*Yibo Liu*](https://github.com/jedibobo), for his advice on building such a nice repository.