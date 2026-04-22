# iSTARMOD — Spectral Subtraction for Chromospheric Activity

**iSTARMOD** is a Python implementation of the spectral subtraction technique designed to quantify chromospheric activity in late-type stars using high-resolution spectra.

It is the modern evolution of the original **STARMOD** code, rewritten to improve usability, modularity, and integration with contemporary scientific workflows.

---

## 🔬 Scientific Background

The spectral subtraction technique isolates chromospheric emission by subtracting a synthetic reference spectrum (constructed from inactive stars) from an observed spectrum.

This allows measurement of **excess emission equivalent widths (EW)** in key activity indicators such as:

- Hα and Balmer series  
- Ca II H & K  
- Ca II infrared triplet (IRT)  
- He I D3, Na I D1/D2  
- Near-IR lines (Paβ, He I λ10830)  

These measurements can be converted into **chromospheric fluxes** using χ-factor calibrations included with the code.

---

## ⚙️ Features

- Iterative least-squares spectral fitting  
- Doppler shift + rotational broadening (Gray profile)  
- Support for **single stars and SB2 binaries**  
- Automatic **equivalent width (EW)** computation  
- Error estimation (Cayrel formula)  
- Batch processing of large spectral datasets  
- Compatible with FITS and simple ASCII spectra  

The code operates **order-by-order on echelle spectra**, improving robustness against instrumental effects and continuum uncertainties.

---
## 🚀 Installation
Option 1 — From GitHub 

Clone the repository and install locally:
```bash
git clone https://github.com/<USER>/iSTARMOD.git
cd iSTARMOD
git checkout v11

python -m venv .venv
source .venv/bin/activate   # On Windows: .venv\\Scripts\\activate

pip install .
```

Option 2 — From GitHub Release Archive (recommended)
1. Download release v11 from the GitHub Releases page
2. Extract the archive
3. Open a terminal in the extracted folder

Then run:
```bash
python -m venv .venv
source .venv/bin/activate
pip install .
```

Dependencies include:

- numpy  
- scipy  
- astropy
- matplotlib

> ⚠️ Note: This project uses Tkinter also for visualization. On some Linux systems,
> you may need to install it separately (e.g., `sudo apt-get install python3-tk`).

---

## ▶️ Usage
Run iSTARMOD from the repository root directory:

```bash
istarmod configsm/example.sm
```
Common options

Enable plotting:
```bash
istarmod configsm/example.sm --plot
```
Enable debugging output:
```bash
istarmod configsm/example.sm --debug
```
## 🖥️ GUI
Launch the graphical interface:
```bash
istarmod-gui
```

## 📁 Important Notes
The .sm configuration file defines the full workflow (input spectra, references, fitting parameters, etc.)
The code expects access to repository folders such as:
configsm/     → configuration files (.sm)
target/       → target spectra
referencesp/  → reference star spectra
data/         → auxiliary files (e.g., lambdas.dat)

➡️ For this reason, it is recommended to run iSTARMOD from the repository root.

## 🧪 Minimal Example

```bash
istarmod configsm/pwand_n3_cah_34_wvl.sm --plot
```

This will:

- build the synthetic spectrum
- subtract it from the observed spectrum
- compute equivalent widths
- optionally display plots

## ▶️ Quick Start

1. Prepare your observed spectrum (FITS or ASCII)
2. Select appropriate **reference star spectrum(s)**
3. Create a configuration file (`.sm`)
4. Run:

```bash
istarmod anyfolder/target.sm --plot
```

>⚠️ Note: Within iStarmod.py you must include the name of the .sm file as an input parameter in the call to the function starmod
>There are two additional default parameters: 'plot' to plot the figure (set to True by default), and debugging to print additional output messages (set to False by default)
>When running a batch of spectra, the plot parameter should be put to 'False'
As noted, in the command starmod("filename.sm") the parameter "filename.sm" can specify its absolute or relative path. If no path is specified the "filename.sm" input file must be found in the same directory where the executable is located. 

---

## 🧾 Configuration File (`.sm`)

The `.sm` file defines the full fitting setup using:

```
KEYWORD = value
```

### Key Sections

#### 1. General Information
```ini
IM_PATH = ./
RES_PATH = ./
OBJ_NAME = target_spectrum.fits
```
>⚠️ Note:
>The name of target spectrum can contain the wildcard "*" to process batch of spectra, with similar names.
>The RES_PATH parameter allows defining the path to the output files

#### 2. Output
```ini
SYN_SPEC = YES
SYN_NAME = synthetic_output
SUB_NAME = subtracted_output
```

#### 3. Fitting Parameters
```ini
N_ITER = 8
PIX_ZONE = 6520 6600
PIX_EXCL = 6520 6524
PIX_EXCL = 6558 6570
```

#### 4. Reference Star (Primary)
```ini
PRM_NAME = ref_star.fits
PRM_RAD = -20.0 var
PRM_ROT = 5.0 var
PRM_WGT = 1.0 fix
```

#### 5. Secondary Star (optional, SB2)
```ini
SEC_NAME = ref2.fits OR NONE
SEC_RAD = 20.0 var
SEC_ROT = 10.0 var
SEC_WGT = 0.3 var
```

#### 6. Spectral Region
```ini
APERTURE = 34
LINE = Halpha
```

#### 7. Additional Algorithm & Visualization Parameters
```ini
WVL_DISPLAY_RANGE = 6550 6580            
USE_RV_VALUES = NO  
```
>⚠️ Note:
>The WVL_DISPLAY_RANGE additional parameter allows defining the wavelength zone to be displayed in the output figure
>The USE_RV_VALUES additional parameter allows storing the RV values calculated in previous runs of the algorithm.

---

## 🔄 Workflow

1. Continuum normalization  
2. Doppler alignment  
3. Rotational broadening  
4. Synthetic spectrum construction  
5. Subtraction  
6. EW computation  

A successful subtraction yields near-zero residuals outside activity lines.

---

## Output

- Synthetic spectrum (`SYN_NAME`)
- Subtracted spectrum (`SUB_NAME`)
- Equivalent widths (EW)
- Diagnostic plots

---

## 🔁 Reproducibility & Paper Figures

This repository includes all components required to reproduce the main figures and results from the associated publication.

### Included

- Example spectra (or retrieval instructions)
- Reference star spectra
- Configuration files (`.sm`)
- Figure-generation scripts
- χ-factor based flux calculation routines

---

### ▶️ Reproducing Spectral Subtraction Figures

Example (Hα, single star):

```bash
python iStarmod.py 
```

Outputs:
- synthetic spectrum
- subtracted spectrum
- diagnostic plots

Expected:
- RMS ≲ 0.5
- clear Hα excess emission

---

### Reproducing Binary Star Cases (SB2)

```bash
python iStarmod.py 
```

Expected:
- correct two-component decomposition
- EW uncertainties ~1–2%

---

### 📁 Paper Configurations

Several configurations used in the publication are stored in:

```
configs/paper/
```

Naming convention:

```
<target>_<line>_<instrument>.sm
```


---

### Reproducibility Checklist

- Install dependencies  
- Download or link datasets  
- Run `.sm` configurations  
- Execute plotting scripts  
- Compare outputs with paper figures  

---

### ⚠️ Notes

Results depend critically on:
- reference star selection
- spectral type matching  
- S/N and normalization quality  

Small variations may affect RMS and EW values.

---

## Applications

- Stellar activity characterization  
- Time-series analysis  
- Spectroscopic surveys  
- Radial velocity activity mitigation  

---

## Citation

If you use this code, please cite:

```
Labarga, F. & Montes, D. (2026)
iSTARMOD: A Python Code to Quantify Chromospheric Activity
The Astronomical Journal, 171, 15
```

---

## 🤝 Contributing

Pull requests are welcome. Please include:
- clear description
- test cases
- documentation updates

---

## 📄 License

Open-source (see LICENSE file).

---

## 🔭 Future Work

- Extended χ calibrations and analysis
- Enhanced visualization tools
- Time series analysis (e.g. Lyapunov exponents, activity period, etc.)
- RV verification
