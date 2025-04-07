# Spectral Line Fitting Tool

Interactive tool for measuring spectral lines using single or double Gaussian fitting. Designed for astronomical spectra saved as `.csv` files.

## ✨ Features

- Interactive selection of the spectral line region
- Fitting with 1 or 2 Gaussian curves
- Calculation of FWHM and Equivalent Width (EW)
- Saves fit plots and results to `.csv`

## 📂 Input Data Format

Input spectra must be **normalized** and saved in `.csv` files with **two required columns**:

- `wave`: wavelength (in Angstroms)
- `flux`: normalized flux

Example:

```csv
wave,flux
4860.1,0.95
4860.2,0.98
4860.3,1.00
...

🧪 Example Data
The example spectra used in this project were taken from the Sloan Digital Sky Survey (SDSS). These spectra were pre-processed and normalized prior to analysis.

## How to Use

```bash
pip install -r requirements.txt

python src/interactive_fitting.py

