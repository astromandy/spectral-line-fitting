#!/usr/bin/env python
# coding: utf-8

# Spectral Line Measurement with Interactive Double Gaussian Fitting

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simpson
from specutils import Spectrum1D
from astropy import units as u
from tkinter import Tk, filedialog, simpledialog, messagebox

# Gaussian functions
def gaussian(x, amp, mu, sigma):
    return amp * np.exp(-(x - mu)**2 / (2 * sigma**2))

def double_gaussian(x, amp1, mu1, sigma1, amp2, mu2, sigma2):
    return gaussian(x, amp1, mu1, sigma1) + gaussian(x, amp2, mu2, sigma2)

# Let user select folder
def select_folder():
    root = Tk()
    root.withdraw()
    folder_selected = filedialog.askdirectory(title="Select folder with CSV spectra")
    return folder_selected

# Load CSV as numpy arrays
def load_csv_spectrum(filepath):
    df = pd.read_csv(filepath)
    return df['wave'].values, df['flux'].values

# Semi-automatic peak fitting
def measure_line_interactively(wave, flux, filename, output_plot_dir):
    while True:
        # Ask user to select region
        plt.figure()
        plt.plot(wave, flux, label='Spectrum')
        plt.title(f'Select region of interest: {filename}\nClick 2 points (start & end)')
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux')
        points = plt.ginput(2, timeout=0)
        plt.close()

        if len(points) < 2:
            print("Insufficient points selected.")
            return None

        x1, x2 = sorted([p[0] for p in points])
        mask = (wave >= x1) & (wave <= x2)
        wave_sel = wave[mask]
        flux_sel = flux[mask]

        # Ask number of peaks
        root = Tk()
        root.withdraw()
        num_peaks = simpledialog.askinteger("Número de picos", "Quantos picos tem a linha? (1 ou 2)", minvalue=1, maxvalue=2)

        try:
            if num_peaks == 1:
                amp_guess = np.max(flux_sel) - np.min(flux_sel)
                mu_guess = wave_sel[np.argmax(flux_sel)]
                sigma_guess = (x2 - x1) / 6
                popt, pcov = curve_fit(gaussian, wave_sel, flux_sel, p0=[amp_guess, mu_guess, sigma_guess])
                fit_flux = gaussian(wave_sel, *popt)
                residual = flux_sel - fit_flux
                fwhm = 2.355 * popt[2]
                fwhm_err = 2.355 * np.sqrt(np.diag(pcov))[2]
                fit_label = "Gaussian Fit"

            else:
                amp1 = np.max(flux_sel) - np.min(flux_sel)
                mu1 = wave_sel[np.argmax(flux_sel)]
                sigma1 = (x2 - x1) / 10
                amp2 = amp1 / 2
                mu2 = mu1 + sigma1 * 2
                sigma2 = sigma1
                p0 = [amp1, mu1, sigma1, amp2, mu2, sigma2]
                popt, pcov = curve_fit(double_gaussian, wave_sel, flux_sel, p0=p0)
                g1 = gaussian(wave_sel, popt[0], popt[1], popt[2])
                g2 = gaussian(wave_sel, popt[3], popt[4], popt[5])
                fit_flux = g1 + g2
                residual = flux_sel - fit_flux
                fwhm = np.nan  # Não aplicável para curva composta
                fwhm_err = np.nan
                fit_label = "Soma de 2 Gaussianas"

                # Se quiser mostrar as gaussianas separadas:
                plt.plot(wave_sel, g1, '--', color='orange', label='Gaussiana 1')
                plt.plot(wave_sel, g2, '--', color='green', label='Gaussiana 2')

            # Plot ajuste
            plt.figure()
            plt.plot(wave_sel, flux_sel, label='Dados')
            plt.plot(wave_sel, fit_flux, label=fit_label, color='r')
            plt.fill_between(wave_sel, fit_flux, alpha=0.2, color='red', label='Área do Ajuste')
            plt.title(f'{filename} - Ajuste da Linha')
            plt.xlabel('Wavelength (Å)')
            plt.ylabel('Flux')
            plt.legend()
            plt.tight_layout()
            plt.show()

            # Confirmação do ajuste
            retry = messagebox.askyesno("Confirmação", "O ajuste está bom?")
            if retry:
                # Calcular área como medida_ew_flux
                ew = simpson(fit_flux, wave_sel)
                ew_err = np.trapz(residual**2, wave_sel)**0.5

                # Salvar plot
                plt.figure()
                plt.plot(wave_sel, flux_sel, label='Dados')
                plt.plot(wave_sel, fit_flux, label=fit_label, color='r')
                plt.fill_between(wave_sel, fit_flux, alpha=0.2, color='red', label='Área do Ajuste')
                plt.title(f'{filename} - Ajuste Final')
                plt.xlabel('Wavelength (Å)')
                plt.ylabel('Flux')
                plt.legend()
                plt.tight_layout()
                plot_path = os.path.join(output_plot_dir, f"{filename}_fit.png")
                plt.savefig(plot_path)
                plt.close()

                return {
                    "filename": filename,
                    "num_peaks": num_peaks,
                    "fwhm": fwhm,
                    "fwhm_err": fwhm_err,
                    "ew": ew,
                    "ew_err": ew_err
                }
        except RuntimeError:
            messagebox.showerror("Erro", "O ajuste falhou. Tente novamente.")

# Continuação: Main script

def main():
    folder = select_folder()
    if not folder:
        print("Nenhuma pasta selecionada.")
        return

    output_dir = os.path.join(folder, "results")
    os.makedirs(output_dir, exist_ok=True)

    output_plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(output_plot_dir, exist_ok=True)

    results = []
    for file in os.listdir(folder):
        if file.endswith(".csv"):
            path = os.path.join(folder, file)
            print(f"\n>>> Processando {file}...")
            wave, flux = load_csv_spectrum(path)
            result = measure_line_interactively(wave, flux, file, output_plot_dir)
            if result:
                results.append(result)

    if results:
        df = pd.DataFrame(results)
        df.to_csv(os.path.join(output_dir, "line_measurements.csv"), index=False)
        print(f"\nResultados salvos em: {os.path.join(output_dir, 'line_measurements.csv')}")
    else:
        print("Nenhuma medida válida realizada.")

if __name__ == "__main__":
    main()

