# TAPE: Tissue-AdaPtive autoEncoder for accurate deconvolution and gene expression analysis

**This model is able to deconvolve bulk RNA-seq data into cell fractions based on scRNA-seq data**

## Installation
TAPE uses PyTorch as its backbones, so the suitable version of PyTorch will accelerate the model. We recommend users to install PyTorch(>=1.8.0) with ***right*** compute platform (CUDA, CPU or ROCm) from its official [website](https://pytorch.org) in advance.

Then the TAPE could be installed from PyPI directly:

```bash
pip install scTAPE
```
## Usage
Needed Files:
1. single-cell reference: txt format, indices are cell types, columns are gene names
2. bulk data: tabular format, needed to specify the seperation ('\t',','or others), indices are sample names, columns are gene names
3. gene length file: used to scale the expression value, columns should contain: [Gene name, Transcript start (bp), Transcript end (bp)]. This is provided in ./data/ directory.
```python

```