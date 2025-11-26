# 2026 ISRO–ESA Heliophysics Workshop — Hands-On Materials

This repository contains all notebooks, data-access examples, and analysis workflows for the hands-on component of the **Joint ISRO–ESA Heliophysics Workshop on Aditya-L1, Solar Orbiter, and Proba-3** (January 19–23, 2026, Thiruvananthapuram, India).

The goal is to provide a unified, Python-based environment where participants can explore and analyse heliophysics data from multiple missions using SunPy, mission-specific tools, and coordinated multi-instrument workflows.

---

## Contents

### **Instrument tutorials (Mon–Wed)**  
Introductory notebooks contributed by each instrument team, including:

- Data access examples  
- Basic workflows and visualisations  
- Instrument-specific context  

**Instruments covered include:**

- **Solar Orbiter:** EUI, STIX, SPICE, PHI, Metis, SoloHI, MAG, SWA, EPD, RPW  
- **Aditya-L1:** SUIT, VELC, HEL1OS, ASPEX, SoLEXS, MAG, PAPA  
- **Proba-3:** ASPIICS  

---

### **Cross-instrument science workflows (Thu–Fri)**  
Focused, multi-mission analysis examples exploring:

- Solar flares (multi-wavelength & multi-mission)  
- Large & small-scale coronal structures  
- In-situ heliospheric variability (solar wind & energetic particles)  

---

### **Common resources**

- Shared Python environment (`env/`)  
- Utility functions and documentation links  

---

# Contributing Instrument Notebooks

Each instrument team is invited to contribute one or more Jupyter notebooks for the hands-on sessions. This section explains how to add your notebook to the repository.

---

## What should be in the notebooks?

A typical notebook should include:

- A short instrument introduction  
- How to access data (examples + links)  
- A basic plotting / quicklook example  
- Suggested exercises or exploration ideas  
- Links to documentation  

Feel free to adapt or extend your tutorial as needed.

---

# 0. Quick Start: How to Add Your Notebook (Recommended)

Most contributors will add notebooks **locally** and then open a Pull Request on GitHub.  
Here is the simplest workflow:

---

## **Step 1 — Fork the repository**

Open the main repo:

**https://github.com/ISRO-ESA-Heliophysics-Workshop/hands-on-materials**

Click **Fork** (top right).  
This creates a personal copy of the repository under your GitHub account.

---

## **Step 2 — Clone your fork**

On your own computer:

```bash
git clone https://github.com/<your-username>/hands-on-materials.git
cd hands-on-materials
```

## Step 3 — Create a new branch

Use a new branch for your contribution:

`git checkout -b add-<instrument>-notebook`

Example:

`git checkout -b add-stix-notebook`


## Step 4 — Add your notebook to the correct folder

Place your notebook inside the appropriate directory:

- `solar-orbiter/<instrument>/`
- `aditya-l1/<instrument>/`
- `proba-3/aspiics/`

Example:

`solar-orbiter/stix/stix_intro.ipynb`

If your instrument folder is missing, feel free to create it.


## Step 5 — Commit and push your changes

Run:
```
git add .
git commit -m "Add <instrument> tutorial notebook"
git push origin add-<instrument>-notebook
```



## Step 6 — Open a Pull Request (PR)

On GitHub, open your fork.  
A banner should appear offering **“Open a Pull Request”** — click it.

Your PR description should include:

- A brief overview of what the notebook covers  
- Any environment notes or extra dependencies  
- Links to relevant documentation  
- Any known issues or potential follow-ups  

We review PRs regularly and merge them as they come in.

---

# Additional Notes for Contributors

## Do not commit large data files

To keep the repository lightweight:

- **Do NOT upload large datasets.**
- Instead:
  - provide download links, or  
  - include code to retrieve data from archives (SOAR, SUIT, ISRO, Proba-3), or  
  - include only very small sample files (<1–2 MB) if necessary.

Large data should be downloaded at runtime.


## Test your notebook (easy version)

Before submitting, please:

- Run your notebook from start to finish  
- Ensure all cells execute without errors  
- Confirm that plots and outputs appear correctly  
- Avoid obscure or heavy Python packages if possible  

In your Pull Request, **please list the Python packages you used**, e.g.:

`Packages used: numpy, matplotlib, astropy, sunpy, ndcube`

This helps the organisers keep the shared environment up to date.

> You do *not* need to understand Conda or environments — just list the packages, and we will handle the rest.


## Optional: Add helper scripts

If your notebook uses repeated helper functions, you may add:

`solar-orbiter/<instrument>/helpers.py`

(or the equivalent folder for your mission).  
This keeps notebooks clean and easier for participants to follow.


## Need help?

If you run into any issues:
- Email me: laura.hayes@dias.ie
- Open an Issue in this repository, or  
- Contact the workshop organisers  

We are happy to help ensure your material integrates smoothly into the hands-on programme.


---


