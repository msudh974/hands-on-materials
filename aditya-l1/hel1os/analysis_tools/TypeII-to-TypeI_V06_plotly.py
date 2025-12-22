"""
Version 6.0 (December 09, 2025) Manju Sudhakar
Code to generate multiple Type I files from Type II spectra. Plotly-based visualization for Colab/Jupyter notebooks, with:
    (1) Science spectra with common CdTe/CZT time interval for spectral fitting
    (2) Detector-specific background spectra with separate CdTe, CZT time intervals if user wishes to select a diferent background period. 
HISTORY: 

Modified from Version 3.0 (Dec 03, 2025) Kiran Kumar Kalisetti 
    - GUI-based tool to be executed from a terminal. 
      Modifications include 
       (1) enabling user to change the end time after right-clicking 
       (2) Removing background file name hard code 

Modified from Version 1.0 (Jan 04, 2025) Manju Sudhakar 
    - GUI-based tool to be executed from a terminal. Purpose
"""


import numpy as np
import pandas as pd
from astropy.io import fits
from itertools import cycle

import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ---------------------------------------------------------------------
# Global paths (set in main) so filter_and_process_files can see them
# ---------------------------------------------------------------------
fits1_path = None  # CdTe Type II PHA
fits2_path = None  # CZT Type II PHA

color_cycle = cycle(["red", "green", "orange", "purple", "black", "cyan", "magenta"])


# ---------------------------------------------------------------------
# Engine Function 1: generate a Type I PHA for a given interval and one detector
# ---------------------------------------------------------------------
def generate_typeI_for_interval(fits_file, output_name, start_time, end_time,
                                is_background=False):
    """
    Generate a Type I PHA file for a single detector over a single time interval.

    Parameters
    ----------
    fits_file : str
        Path to the Type II PHA FITS file.
    output_name : str
        Base name for the detector (e.g. 'cdte', 'czt').
    start_time, end_time : str or pandas.Timestamp or datetime
        Absolute times (same system as DATE_OBS/TIME_OBS in the FITS header).
    is_background : bool
        If True, include 'back' in the output filename.
    """
    start_time = pd.to_datetime(start_time)
    end_time = pd.to_datetime(end_time)

    with fits.open(fits_file) as hdul:
        spectrum_data = hdul["SPECTRUM"].data
        header = hdul["SPECTRUM"].header

        tstart = spectrum_data["TSTART"]
        tstop = spectrum_data["TSTOP"]

        # Convert start/end times to seconds relative to DATE_OBS/TIME_OBS
        base_time = pd.to_datetime(header["DATE_OBS"] + " " + header["TIME_OBS"])
        start_time_seconds = (start_time - base_time).total_seconds()
        end_time_seconds = (end_time - base_time).total_seconds()

        # Find the rows that cover the requested times
        start_indices = np.where(
            (tstart <= start_time_seconds) & (tstop >= start_time_seconds)
        )[0]
        end_indices = np.where(
            (tstart <= end_time_seconds) & (tstop >= end_time_seconds)
        )[0]

        if len(start_indices) == 0 or len(end_indices) == 0:
            print(f"WARNING: No rows found in {fits_file} for the selected interval.")
            return

        start_row = start_indices[0]
        end_row = end_indices[0]

        # Sum counts and errors across the selected rows
        summed_counts = np.sum(
            spectrum_data["COUNTS"][start_row : end_row + 1], axis=0
        )
        summed_stat_err = np.sqrt(
            np.sum(
                spectrum_data["STAT_ERR"][start_row : end_row + 1] ** 2,
                axis=0,
            )
        )
        summed_exposure = np.sum(
            spectrum_data["EXPOSURE"][start_row : end_row + 1]
        )

        # Set ANCRFILE and RESPFILE based on the detector type
        if output_name.lower() == "cdte":
            ancrfile = "hel1os_cdte_arf_v03.fits"
            respfile = "hel1os_cdte_srf_v04.fits"
        elif output_name.lower() == "czt":
            ancrfile = "hel1os_czt_arf_v03.fits"
            respfile = "hel1os_czt_srf_v04.fits"
        else:
            ancrfile = "none"
            respfile = "none"

        # Time strings for filenames
        start_time_str = start_time.strftime("%H%M%S")
        end_time_str = end_time.strftime("%H%M%S")

        if is_background:
            file_name = f"{output_name.lower()}_back_{start_time_str}_{end_time_str}.fits"
        else:
            file_name = f"{output_name.lower()}_{start_time_str}_{end_time_str}.fits"

        # Build the Type I PHA table
        cols = [
            fits.Column(
                name="CHANNEL",
                format="J",
                array=np.arange(len(summed_counts)) + 1,
            ),
            fits.Column(name="COUNTS", format="J", array=summed_counts),
            fits.Column(name="STAT_ERR", format="E", array=summed_stat_err),
        ]
        pha_hdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))

        # OGIP-compliant header keywords (same as science; only filename differs)
        headers = {
            "EXTNAME": "SPECTRUM",
            "TELESCOP": "ADITYA-L1",
            "INSTRUME": "HEL1OS",
            "RESPFILE": respfile,
            "ANCRFILE": ancrfile,
            "BACKFILE": "none",
            "CORRFILE": "none",
            "CORRSCAL": 0.0,
            "AREASCAL": 1.0,
            "BACKSCAL": 1.0,
            "HDUCLASS": "OGIP",
            "HDUCLAS1": "SPECTRUM",
            "HDUCLAS2": "TYPE:I",
            "HDUVERS": "1.2.1",
            "POISSERR": "false",
            "CHANTYPE": "PHA",
            "DETCHANS": len(summed_counts),
            "TLMIN": 1,
            "TLMAX": len(summed_counts),
            "EXPOSURE": float(summed_exposure),
            "DATE_OBS": start_time.strftime("%Y-%m-%d"),
            "TIME_OBS": start_time.strftime("%H:%M:%S.%f"),
            "DATE_END": end_time.strftime("%Y-%m-%d"),
            "TIME_END": end_time.strftime("%H:%M:%S.%f"),
            "TSTART": float(tstart[start_row]),
            "TSTOP": float(tstop[end_row]),
        }

        for key, value in headers.items():
            pha_hdu.header[key] = value

        pha_hdu.writeto(file_name, overwrite=False)
        tag = "background" if is_background else "science"
        print(f"{tag.capitalize()} Type I file generated: {file_name}")


# ---------------------------------------------------------------------
# Engine Function 2 Process both detectors for a common time interval and generate Science Spectrum
# ---------------------------------------------------------------------
def filter_and_process_files(start_time, end_time):
    """
    Filter both Type II FITS files between start_time and end_time
    and generate science Type I PHA files for CdTe and CZT. This is a wrapper. It calls "generate_typeI_for_interval TWICE. Once for CdTe and once for CZT. 
    """
    global fits1_path, fits2_path

    if fits1_path is None or fits2_path is None:
        raise RuntimeError(
            "fits1_path and fits2_path must be set before calling filter_and_process_files()."
        )

    start_time = pd.to_datetime(start_time)
    end_time = pd.to_datetime(end_time)

    # Science spectra: same interval for both
    generate_typeI_for_interval(fits1_path, "cdte", start_time, end_time,
                                is_background=False)
    generate_typeI_for_interval(fits2_path, "czt", start_time, end_time,
                                is_background=False)


# ---------------------------------------------------------------------
# Load Type II PHA as a time profile (for plotting)
# ---------------------------------------------------------------------
def load_time_profile(fits_file):
    """
    Load a Type II PHA FITS file and return a DataFrame
    with mid-bin times and total counts per time bin. This function takes a Type II FITS file 
    and produces a Pandas dataframe with with two columns (1) Time - is the midpoint of each time bin 
    in actual date/time UTC (2) Counts - the total counts in that time bin summed over all energy channels.
    """
    with fits.open(fits_file) as hdul:
        spectrum_data = hdul["SPECTRUM"].data
        tstart = spectrum_data["TSTART"]
        tstop = spectrum_data["TSTOP"]
        counts = np.sum(spectrum_data["COUNTS"], axis=1)

        time_midpoints = tstart + (tstop - tstart) / 2.0

        header = hdul["SPECTRUM"].header
        base_time = pd.to_datetime(header["DATE_OBS"] + " " + header["TIME_OBS"])
        times = [base_time + pd.to_timedelta(tm, unit="s") for tm in time_midpoints]

    return pd.DataFrame({"Time": times, "Counts": counts})


# ---------------------------------------------------------------------
# Plotly lightcurve - we make a figure with two rows.  
# ---------------------------------------------------------------------
def make_lightcurve_figure(cdte_df, czt_df):
    """
    Create a Plotly Figure with CdTe and CZT lightcurves. Top row is CdTe light curve and bottom 
    figure is CZT light curve. Both y-axes are on log scale (counts). This figure helps select the flare
    and background time periods. 
    """
    base_fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        subplot_titles=("CdTe Lightcurve", "CZT Lightcurve"),
    )

    fig = go.Figure(base_fig)

    # CdTe
    fig.add_scatter(
        row=1,
        col=1,
        x=cdte_df["Time"],
        y=cdte_df["Counts"],
        mode="lines",
        name="CdTe",
    )

    # CZT
    fig.add_scatter(
        row=2,
        col=1,
        x=czt_df["Time"],
        y=czt_df["Counts"],
        mode="lines",
        name="CZT",
    )

    fig.update_yaxes(type="log", title_text="Counts", row=1, col=1)
    fig.update_yaxes(type="log", title_text="Counts", row=2, col=1)
    fig.update_xaxes(title_text="Time [UTC]", row=2, col=1)

    date_str = cdte_df["Time"].iloc[0].strftime("%Y-%m-%d")
    fig.update_layout(
        title_text=f"HEL1OS Type II Lightcurves ({date_str})",
        hovermode="x unified",
    )

    return fig


# ---------------------------------------------------------------------
# Detector-specific background helper function
# ---------------------------------------------------------------------
def generate_detector_specific_background():
    """
    Ask user for CdTe and CZT background intervals (separately),
    and generate background Type I PHA files with 'back' in the names. 
    """
    global fits1_path, fits2_path

    # CdTe background
    while True:
        try:
            s_cd = input(
                "Enter CdTe BACKGROUND START time (YYYY-mm-dd HH:MM:SS) "
                "or press ENTER to skip CdTe background: "
            ).strip()
        except EOFError:
            print("\nEOF received, skipping CdTe background.")
            s_cd = ""
        if not s_cd:
            print("Skipping CdTe background.")
            break

        try:
            e_cd = input(
                "Enter CdTe BACKGROUND END time   (YYYY-mm-dd HH:MM:SS): "
            ).strip()
        except EOFError:
            print("\nEOF received, skipping CdTe background.")
            break

        if not e_cd:
            print("Empty CdTe end time; skipping CdTe background.")
            break

        try:
            start_cd = pd.to_datetime(s_cd)
            end_cd = pd.to_datetime(e_cd)
        except Exception as exc:
            print(f"Could not parse CdTe times: {exc}")
            continue

        if end_cd <= start_cd:
            print("CdTe END time must be after START time. Try again.")
            continue

        print(
            f"\nGenerating CdTe BACKGROUND for interval:\n  {start_cd}  -->  {end_cd}"
        )
        generate_typeI_for_interval(
            fits1_path, "cdte", start_cd, end_cd, is_background=True
        )
        break

    # CZT background
    while True:
        try:
            s_czt = input(
                "Enter CZT BACKGROUND START time  (YYYY-mm-dd HH:MM:SS) "
                "or press ENTER to skip CZT background: "
            ).strip()
        except EOFError:
            print("\nEOF received, skipping CZT background.")
            s_czt = ""
        if not s_czt:
            print("Skipping CZT background.")
            break

        try:
            e_czt = input(
                "Enter CZT BACKGROUND END time    (YYYY-mm-dd HH:MM:SS): "
            ).strip()
        except EOFError:
            print("\nEOF received, skipping CZT background.")
            break

        if not e_czt:
            print("Empty CZT end time; skipping CZT background.")
            break

        try:
            start_czt = pd.to_datetime(s_czt)
            end_czt = pd.to_datetime(e_czt)
        except Exception as exc:
            print(f"Could not parse CZT times: {exc}")
            continue

        if end_czt <= start_czt:
            print("CZT END time must be after START time. Try again.")
            continue

        print(
            f"\nGenerating CZT BACKGROUND for interval:\n  {start_czt}  -->  {end_czt}"
        )
        generate_typeI_for_interval(
            fits2_path, "czt", start_czt, end_czt, is_background=True
        )
        break


# ---------------------------------------------------------------------
# The Main Function: User entry point (works with %run and plain python)
# ---------------------------------------------------------------------
def main(cdte_fits, czt_fits):
    """
    (1) Stores the two input filenames (CdTe and CZT) in fits1_path and fits2_path. 
    (2) Use "load_time_profile" to get lightcurve dataframes. 
    (3) Call "make_lightcurve_figure" and display the figure. 
    (4) Enter a "while True:" loop 
        (a) Asks for science start time
            (i) If you press "Enter" - exit the program
        (b) Ask for science end time.
        (c) Convert to datetime and check that end > start. 
        (d) Call "filter_process_files(start,end) which produces science Type I files for 
        both CdTe and CZT
        (e) Print a separator line and show a small menu 
        Next action? [N]ew science spectrum, [B]ackground (detector specific), [Q]uit:
            (i)  N/n - loops back and asks for a new science interval.
            (ii) B/b - call "generate_detector_specific_background()" and then return to menu
            (iii) Q/q - exit the function
    """
    global fits1_path, fits2_path, color_cycle

    fits1_path = cdte_fits
    fits2_path = czt_fits

    color_cycle = cycle(["red", "green", "orange", "purple", "black", "cyan", "magenta"])

    cdte_df = load_time_profile(fits1_path)
    czt_df = load_time_profile(fits2_path)

    fig = make_lightcurve_figure(cdte_df, czt_df)

    # Display
    try:
        from IPython.display import display

        display(fig)
    except ImportError:
        fig.show()

    print("\nLightcurves displayed.")
    print("You can now generate SCIENCE Type I PHA files by typing time ranges.")
    print("Times must be in the same absolute system as the lightcurve (DATE_OBS / TIME_OBS).")
    print("Example format:  2023-12-15 07:10:00\n")

    while True:
        # --- SCIENCE interval ---
        try:
            s = input(
                "Enter START time for SCIENCE spectrum (YYYY-mm-dd HH:MM:SS) "
                "or press ENTER to quit: "
            ).strip()
        except EOFError:
            print("\nEOF received, exiting.")
            break

        if not s:
            print("No start time entered. Exiting.")
            break

        try:
            e = input(
                "Enter END time   for SCIENCE spectrum (YYYY-mm-dd HH:MM:SS): "
            ).strip()
        except EOFError:
            print("\nEOF received, aborting this interval.")
            break

        if not e:
            print("Empty end time; skipping this interval.")
            continue

        try:
            start = pd.to_datetime(s)
            end = pd.to_datetime(e)
        except Exception as exc:
            print(f"Could not parse times: {exc}")
            continue

        if end <= start:
            print("END time must be after START time. Try again.")
            continue

        print(f"\nGenerating SCIENCE Type I PHA files for interval:\n  {start}  -->  {end}")
        filter_and_process_files(start, end)
        print("-" * 60)

        # --- Next action menu ---
        while True:
            try:
                choice = input(
                    "\nNext action? [N]ew science spectrum, "
                    "[B]ackground (detector specific), [Q]uit: "
                ).strip().lower()
            except EOFError:
                print("\nEOF received, exiting.")
                return fig

            if not choice:
                print("No choice entered. Exiting.")
                return fig

            if choice.startswith("n"):
                # Back to outer loop: ask for next SCIENCE interval
                break
            elif choice.startswith("b"):
                generate_detector_specific_background()
                # After background, ask again what to do next
                continue
            elif choice.startswith("q"):
                print("Exiting.")
                return fig
            else:
                print("Please enter N, B, or Q.")

    return fig


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage:")
        print("  python TypeII-to-TypeI_V06_plotly.py <cdte_TypeII.fits> <czt_TypeII.fits>")
        sys.exit(1)

    cdte_path = sys.argv[1]
    czt_path = sys.argv[2]
    main(cdte_path, czt_path)

