# HEL1OS Instrument Data Analysis

Welcome to the hands-on material for the HEL1OS instrument onboard Aditya-L1. Open the **[Instrument Guide](HEL1OS.md)** for a brief overview of the HEL1OS instrument.  

<!--Run the demo notebook [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/msudh974/hands-on-materials/main?labpath=aditya-l1%2Fhel1os%2FTiming-tools-demo.ipynb) to learn how to plot spectrally resolved light curves from the HEL1OS Level-1 Event data. Loading Binder may take a few minutes. 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/msudh974/hands-on-materials/blob/main/aditya-l1/hel1os/Timing-tools-demo.ipynb)-->  

### Running Locally using Jupyter  

1. **Install and Launch Jupyter**
   ```bash
    jupyter lab
    ```
   
2.  **Install requirements in the first cell:**
    ```bash
    !pip install -r requirements.txt
    ```

3.  **Clone the repository:**
    ```bash
    !git clone [https://github.com/msudh974/hands-on-materials.git](https://github.com/msudh974/hands-on-materials.git)
    !cd hands-on-materials/aditya-l1/hel1os
    ```

4.  **Create a new folder to store the event data:**
    ```bash
    !mkdir data
    !cd data
    ```

<!--Run the **[Demo Notebook](demo1.ipynb)** to see how to load data.
Try the **[Analysis Tools](analysis_tools/)** folder for advanced processing.-->  

**Scripts and Notebooks Prepared and Tested by:**
**Manju Sudhakar**<sup>1</sup>, **Kiran Kumar Kalisetti**<sup>1,2</sup>

<sup>1</sup> Space Astrnomy Group, U.R. Rao Satellite Centre (URSC), ISRO, Bangalore, India <br>
<sup>2</sup> Department of Physics, Indian Institute of Science (IISc), Bangalore, India

**Contact:** msudh2021@gmail.com, manjus@ursc.gov.in



---
*Note: Make sure to install the libraries in `requirements.txt` before running the notebooks.*


