from setuptools import setup, find_packages
setup(
    name = 'scTAPE',
    version = '0.1.3',
    description = 'deep learning tools for bulk RNA-seq deconvolution and gene expression analysis',
    author = 'Yanshuo CHEN',
    author_email = 'poseidonchan@icloud.com',
    url = 'https://github.com/poseidonchan/TAPE',
    license = 'GPL-3.0 License',
    packages = find_packages(),
    python_requires='>=3.7',
    platforms = 'any',
    install_requires = [
        'torch>=1.8.0',
        'numpy>=1.19',
        'pandas>=1.3.1',
        'matplotlib>=3.4',
        'anndata==0.7.6',
        'tqdm>=4.6',
        'scikit-learn>=0.23',
        'seaborn>=0.11'
    ],

)