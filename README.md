# CCWT_EXAFS
# Continuous Cauchy Wavelet Transform of EXAFS

This repository provides a Python implementation of the Continuous Cauchy Wavelet Transform (CCWT) applied to EXAFS χ(k). It enables CCWT analysis of EXAFS data, visualizing Fourier transforms and generating CCWT 2D maps (k vs R) for qualitative EXAFS data analysis.

## 📌 Features

Interpolation of EXAFS χ(k) data

Fourier Transform computation and visualization

Continuous Cauchy Wavelet Transform (CCWT) calculation

High-quality plots with Matplotlib

Option to save figures automatically in /figures folder

Simple, standalone Python script

📁 Repository Structure
exafs-cauchy-wavelet/
│
├── exafsCCWT.py            # Main Python script
├── znfoil_k.dat            # Example EXAFS data
├── requirements.txt        # Python dependencies
├── README.md               # This file
├── LICENSE                 # MIT License
└── figures/                # Generated plots

📖 Usage

Run the main script:

```bash
python exafsCCWT_02.py
```

The script will generate:

Interpolated χ(k) data plot

Fourier Transform magnitude plot

Continuous Cauchy Wavelet Transform (CCWT) 2D map

Save plots automatically:

The script can be modified to save plots using Matplotlib:

plt.savefig(f'figures/plot_name.png', dpi=300)

📦 Dependencies

numpy >= 1.21

scipy >= 1.7

matplotlib >= 3.5

Install all dependencies with:


```bash
pip install -r requirements.txt
```

🔬 References
1. Khan, L. U., Jabeen, N., & et al. (2021). Investigating local structure of ion-implanted ($\text{Ni}^{2+}$) and thermally annealed rocksalt CoO film by EXAFS simulation using evolutionary algorithm. ACS Applied Energy Materials, 4(3), 2049–2055. https://doi.org/10.1021/acsaem.0c02740
2. Khan, L. U., Khan, Z. U., Blois, L., Tabassam, L., Brito, H. F., & Figueroa, S. J. A. (2023). Strategy to probe the local atomic structure of luminescent rare earth complexes by X-ray absorption near-edge spectroscopy simulation using a machine learning-based PyFitIt approach. Inorganic Chemistry, 62(6), 2738–2750. https://doi.org/10.1021/acs.inorgchem.2c03823

🔬 Method

Munoz M., Argoul P., & Farges F. (2003).
Continuous Cauchy wavelet transform analyses of EXAFS spectra: a qualitative approach.
American Mineralogist, 88, 694–700.


✨ Author

Latif Ullah Khan
Beamline Scientist — BM08-XAFS/XRF Beamline, SESAME Light Source

⚡ License

MIT License

