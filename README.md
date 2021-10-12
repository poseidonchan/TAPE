# TAPE: Tissue-AdaPtive autoEncoder for accurate deconvolution and gene expression analysis

**This model is able to accurately deconvolve bulk RNA-seq data into cell fractions based on scRNA-seq data**.

More ***IMPORTANTLY***, TAPE enables fast **cell-type-specific gene expression prediction**, users could take this advantage to analyze gene expression at cell-type level. 

## Setup
TAPE uses PyTorch as its backbones, so the suitable version of PyTorch will accelerate the model. We recommend users to install PyTorch(>=1.8.0) with ***right*** compute platform (CUDA, CPU or ROCm) from its official [website](https://pytorch.org) in advance.

For example, we used the latest released Nvidia GPU RTX3090, so we choose the CUDA version 11.1 and the command is:

```bash
pip install torch==1.8.0+cu111 torchvision==0.9.0+cu111 torchaudio==0.8.0 -f https://download.pytorch.org/whl/torch_stable.html
```

If PyTorch is successfully installed, then the TAPE could be installed from PyPI directly:

```bash
pip install scTAPE
```
## Usage
Needed Files:
1. single-cell reference: txt format, indices are cell types, columns are gene names
2. bulk data: tabular format, needed to specify the seperation ('\t',','or others), indices are sample names, columns are gene names
3. gene length file: used to scale the expression value, columns should contain: [Gene name, Transcript start (bp), Transcript end (bp)]. This is provided in ./data/ directory.
```python
# basic example
from TAPE import Deconvolution
SignatureMatrix, CellFractionPrediction = \
    Deconvolution(sc_ref, bulkdata, sep='\t',
                  datatype='TPM', genelenfile='./GeneLength.txt',
                  mode='overall', adaptive=True,
                  save_model_name=None)
```
parameters:

1. datatype: use '**TPM**', '**FPKM**' or '**counts**', this should be the same with your bulk datatype.
2. mode: '**overall**' or '**high-resolution**'. If you need signature matrix for each sample, use 'high-resolution' mode.
3. adaptive: **True** or **False**. If this is False, then it would not predict signature matrix, the return will be ***None***

## Example
The example is listed in the Test directory. Please run the example to get familiar with TAPE.

## Acknowledgement
Special thanks to *Mengyue Sun*, for his help to accelerate the sampling process (in the simulation.py)