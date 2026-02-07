# NMR Diffusion Analysis

MATLAB tools for estimating diffusion coefficients from NMR (Nuclear Magnetic Resonance) data processed using Bruker TopSpin.

## Overview

**Problem:** Extracting accurate diffusion coefficients from NMR signal attenuation data requires fitting complex exponential models with robust initial guess strategies.

**Solution:** A suite of MATLAB scripts that perform mono-exponential and bi-exponential curve fitting with logarithmic initial guess iteration, signal normalization, and automated data export.

## Scripts

| Script | Purpose | Method |
|--------|---------|--------|
| `ADC_m_Fit_IteratingGuess.m` | Main diffusion coefficient estimator | Mono-exponential fit with iterating initial guesses |
| `Bi_exponential_ADC_plot.m` | Dual-model comparison | Mono + bi-exponential fitting with visualization |
| `Fitting_datapoints_through_equations.m` | General curve fitting | Least-squares fitting for custom equations |
| `Normalize_report.m` | Data normalization | Normalizes NMR peak intensities and exports to CSV |

## How It Works

### 1. Data Loading

Reads NMR data exported from TopSpin as `.xlsx` files containing b-values (gradient strengths) and signal intensities for multiple peaks.

```matlab
[file, path] = uigetfile('*.xlsx*', 'Select Data File');
filename = fullfile(path, file);
og = xlsread(filename);
```

### 2. Mono-Exponential Fitting

Fits signal attenuation to the Stejskal-Tanner equation:

**S(b) = S(0) * exp(-b * D)**

Where D is the apparent diffusion coefficient (ADC).

```matlab
funm = @(D, x) exp(-(B{k} - B{k}(1)) .* D(1));
[fit_m{k}, residual_m(k,1)] = lsqcurvefit(funm, D_guess_m(g_m), B{k}, S{k}, lb_m, ub_m, options);
```

### 3. Bi-Exponential Fitting

Models restricted and unrestricted diffusion components:

**S(b) = A * exp(-b * D1) + (1-A) * exp(-b * D2)**

```matlab
funb = @(D, x) D(1)*exp(-(B{k}-B{k}(1)).*D(2)) + (1-D(1))*(exp(-(B{k}-B{k}(1)).*D(3)));
[fit_b{k}, residual_b(k,1)] = lsqcurvefit(funb, D_guess_b{g_b}, B{k}, S{k}, lb, ub, options);
```

### 4. Initial Guess Strategy

Uses logarithmic spacing across the parameter space to avoid local minima:

```matlab
guessesPerDecade = 4;
logRange_m = log10(upperBound_m) - log10(lowerBound_m);
numIntervals_m = ceil(logRange_m) * guessesPerDecade;
```

The algorithm tests multiple initial guesses per decade of the parameter range and keeps the best fit (lowest residual).

### 5. Output

- Semi-log signal attenuation plots (`.fig`)
- ADC values table exported to CSV
- MATLAB workspace saved for reproducibility (`.mat`)

| Output | Description |
|--------|-------------|
| `Signal_Attenuation_*.fig` | Publication-ready signal decay plots |
| `ADC_Table_*.csv` | Diffusion coefficients per peak |
| `work_ADC_fit_*.mat` | Full workspace for reloading |

## Parameter Bounds

| Parameter | Lower Bound | Upper Bound | Unit |
|-----------|-------------|-------------|------|
| D (mono) | 10^-9 | 0.1 | mm^2/s |
| A (fraction) | 0.01 | 1.0 | dimensionless |
| D1, D2 (bi) | 10^-9 | 0.01 | mm^2/s |

## Requirements

- MATLAB R2019b or later
- Optimization Toolbox (`lsqcurvefit`)
- Input data: Excel files with b-values and normalized signal intensities

## Usage

1. Run `ADC_m_Fit_IteratingGuess.m` in MATLAB
2. Select your data file when prompted
3. Enter experiment parameters (number of peaks, points per scan, plot layout)
4. Results are automatically saved and displayed

```matlab
% Example: Run the main script
run('ADC_m_Fit_IteratingGuess.m')
% Follow the dialog prompts to select data and enter parameters
```

## Context

This project was developed for NMR diffusion-ordered spectroscopy (DOSY) research, where accurate diffusion coefficient estimation is critical for:

- Molecular size determination
- Mixture component identification
- Polymer characterization
- Drug delivery studies

## License

MIT License
