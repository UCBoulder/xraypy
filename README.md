# XRDpy

Tools to assist in the operation and analysis of the WAXS/SAXS system in CU Physics. This software, installed on the computer running the system, can be used to assist in GIWAXS experiments, and can be installed on any computer to assist in data stitching, transformations, and reductions.

The system utilizes a Dectris Eiger R 1M CCD array to produce images of the X-ray scattering for WAXS/SAXS/GIWAXS/GISAXS experiments. The detector has 1065 rows and 1030 columns of pixels. However, there is a band of 37 pixel rows that are not active in the middle of the array. The eiger_stitch command for the saxs program produces raw TIFF files with the file format R-C-S where R represents a row, C a column, and S is a number either 1 or 2 for the two images used to stitch together data in the band of inactive pixels. A motor moves the detector to two positions per stitch so that the two images overlap in a way to fill in the band. Any additional rows or columns of stitching will have a 10 pixel overlap. This causes 3 distinct regions in the stitch:

- Regions that get half the counts (54 rows of pixels in the middle, and 37 rows at the top and bottom)
- Regions that get the normal amount of counts (parts where the two stitched images overlap)
- Regions that get double counts (10 columns of pixels where separate stitch images overlap)

eiger_stitch will produce a stitched image which accounts for these regions by multiplying the half-count regions by 2 and the double-count regions by 0.5; however, there are some issues with this if you care about details. If you save the raw R-C-S files, XRDpy is able to stitch them together in a slightly different way. XRDpy will:

- locate dead, hot, and problematic pixels and mask them
- stitch the images together
- instead of adjusting the counts, it will track how many images contributed to a single pixel in the overall stitch.
- produce a stitch without exposure time adjustment
- produce a "flat field" which can tell your analysis software how sensitive each pixel is
    - Pixels which only had 1 image contribute to it are considered "half-senstive"
    - Pixels which have 2 images contributing to it are considered "full-sensitive"
    - Pixels which have 3 images contributing to it are considered "three-halves-sensitive"
    - 

# Installation

It is recommended to create a unique Anaconda environment for the installation because the install will generate globally callable commands from the terminal, and it is nice to isolate these within an Anaconda environment. This guide assumes you will call the environment "XRD", but you can choose any name you wish.

<!--There are two options to install: from PyPI (with pip) or from source.-->

## Dependencies

XRDpy depends on several other packages. All methods of installation below will automatically install the required dependencies in your conda environment. These are:

- numpy
- matplotlib
- pyFAI
- pyyaml
- fabio
- pathlib
- pyside6
- pyopencl

An optional dependency is [gixpy](https://github.com/etortorici/gixpy) which is a C-compiled package for the transform. In order to install gixpy, you must be running on Windows and have MSVC installed to compile the package. It can be installed with pip if so. Having gixpy installed will allow the GIX transform to run very fast; if it isn't installed, it will run a pure-Python implementation of the transform which is a few orders of magnitude slower.

## Setting up Conda Environment

In Windows, open an Anaconda prompt and run

`conda create --new xrd`

`conda activate xrd`

## Install from source

To install from source code, clone the repository and `cd` into the repository and run (make sure you have the conda environment activated)

`pip install -U .`

# Included CLI and GUI programs

## pyFAI-calib2

This is included from pyFAI. Run `pyFAI-calib2` as a first step to generate a PONI file and call it "cal.poni" to place in the directory with the raw data. Using the auto-stitched TIFF file from a AgBh exposure is preferred.

PONI refers to *point of normal incidence* and is the point on the detector that intersects the shortest line between the sample and the detector. In the orientation where the incident beam is normal to the detector and the sample is brought into the beam path, the PONI is the same as the beam center on the detector.

This GUI can also be used to create and save masks for images.

### Calibration procedure

- In the toolbar, select "Y" icon and select "Orient Y-axis downward"
- In the right panel:
    - Set wavelength to 1.5418 Angstrom
    - Select AgBh from calibrants
    - Select the ... button under "Detector" and then load from file
        - the install should save .h5 files for each eiger_stitch size in "Documents/XRDpy/detectors/". Select the file matching the calibration's stitch size.
    - Designate the calibration TIFF
    - Optionally set a mask file if you have one saved (otherwise, you can make one in the next step)
- Press "Next" to get to the mask making step. Once a mask is made, Press "Next to get to the ring selection tool
- Select the first few rings and then "extract more rings"
- Select "Next" to begin fitting. Select the "SAXS constraints" and then fit.
- Save the PONI file as "cal.poni" with your data.

The PONI file contains all the relevant parameters (except the exposure time, incident angle, and the detector tilt angle).

## `XRD-stitch`

This will stitch all the raw files in the current working directory. It won't adjust the intensity of pixels with different exposure times; instead it will output two TIFF files:

- *raw-stitched-data.tif*: The stitched data without adjustment of different pixel's exposure times.
- *stitched-exposure-time.tif*: A TIFF where each pixel's intensity represents the exposure time (in seconds) that that pixel in the data set had due to the stitching. This can be used to adjust the intensity of the data in *raw-stitched-data.tif* while keeping track of weights.

It is not necessary to run `XRD-stitch` before transforming GIWAXS data with `XRD-film`, because, that command will also stitch. This command is meant for powder patterns using WAXS, and will copy a reduction Jupyter Notebook template to the data directory.

### Arguments

- Exposure time (in seconds)
    - This is exposure time given to `eiger-stitch` so it represents the total exposure time of two overlapping images (the exposure time of a single raw image is half of `<time in seconds>`).
    - `<time in seconds>` must be an integer.
    - optional if `"exposure"` is defined in *params.yaml*.
    - eg: `XRD-stitch 1200`.

### Flags
- `--incident INCIDENT_ANGLE` or `-I INCIDENT_ANGLE`
    - This is the angle of incidence of the beam with the sample; i.e. how many degrees the sample is tipped wrt to the beam.
    - Optional to give to XRD-stitch, and can be defined later. This will also get saved in *params.yaml*
    - eg: `XRD-stitch 1200 --incident 0.2`
- `--tilt TILT_ANGLE` or `-T TILT_ANGLE`
    - This is fully optional.
    - This will rotate the detector clockwise (positive direction) relative to the sample so that you can align the sample normal to the detector's y-direction if the sample is tilted.
- `--dir "PATH TO DATA"` or `-D "PATH TO DATA"`
    - By default, it will look for the data in the current working directory, but you can specify another location. All the saved files will be saved here as well.
- `--plot` or `-P`
    - Produce plots of the stitched data using matplotlib.
- `--override` or `-O`
    - By default, if `XRD-stitch` has already been ran and has saved data, running it again will just load the data found in the working directory. This tag will delete these first and re-stitch the data.

## `XRD-film`

Running this will look for data produced by `XRD-stitch`, and if that isn't found will first stitch the raw data. It will also copy a GIWAXS reduction template Jupyter Notebook in the data directory.

This will transform data taken on aligned thin films (e.g. GIWAXS). An incident angle must be defined for the transform to work. This works for grazing incidence as well as non-grazing incidence.

If `"incident"` isn't already defined in *params.yaml*, it must be provided with a flag. By default, both a cake ($|\mathbf{q}|$ vs $\Omega$) and a reciprocal lattice ($q_z$ vs $q_{xy}$) plot will be made and saved as TIFFs. The plots will have intensity adjusted due to the relative exposure at that location, but the TIFF will save the data unadjusted.

For both the cake and reciprocal lattice plots, two individual TIFF files will be produced (one for the intensity per bin and one for exposure time per bin) as well as yaml file containing the first and last value of each axis as well as the number of bins on that axis.

### Flags

- `--exposure EXPOSURE_TIME` or `-E EXPOSURE_TIME`
    - This is exposure time given to `eiger-stitch` so it represents the total exposure time of two overlapping images (the exposure time of a single raw image is half of `EXPOSURE_TIME`).
    - `EXPOSURE_TIME` must be an integer (and is in seconds).
    - optional if `"exposure"` is defined in *params.yaml*.
    - eg: `XRD-stitch --exposure 1200`.

- `--incident INCIDENT_ANGLE` or `-I INCIDENT_ANGLE`
    - This is the angle of incidence of the beam with the sample; i.e. how many degrees the sample is tipped wrt to the beam.
    - Optional to give to XRD-stitch, and can be defined later. This will also get saved in *params.yaml*
    - eg: `XRD-stitch 1200 --incident 0.2`
- `--tilt TILT_ANGLE` or `-T TILT_ANGLE`
    - This is fully optional.
    - This will rotate the detector clockwise (positive direction) relative to the sample so that you can align the sample normal to the detector's y-direction if the sample is tilted.
- `--dir "PATH TO DATA"` or `-D "PATH TO DATA"`
    - By default, it will look for the data in the current working directory, but you can specify another location. All the saved files will be saved here as well.
- `--plot` or `-P`
    - Plot everything via matplotlib
- `--override` or `-O`
    - By default, it will first look for previously transformed data and load that, instead of doing the transform from scratch. This tag will delete previously saved data and redo the transform.
- `--cake` or `-C`
    - This will produce a cake plot (2D integration) via pyFAI.
- `--reduce` or `-R`
    - This will produce a reduction (1D integration) via pyFAI

## `XRD-imacro`

Running this will produce a macro in Documents/XRDpy/Macros to scan om in order to tune the angle of incidence for thin film studies. The scan takes approximately 1 minute for every 10 angles (6 seconds / angle).

The program will print the command you need to feed to the SAXS program to run the macro.

Note: Macro names include date and hour, so if you save two macros in the same hour, they will overwrite.

### Arguments

- Starting angle (in degrees)
    - The first angle in the scan (smallest or most negative)
- End angle (in degrees)
    - The last angle in the scan (largest or least negative)
- Angle step size (in degrees)
    - how large of a step to take between each angle in the scan

e.g. `XRD-imacro 0 1 .01`

### Flags
- `--name STRING_TO_ADD_TO_NAME` or `-N STRING_TO_ADD_TO_NAME`
    - Optionally add a short tag to the TIFF file names generated
- `--clear` or `-C`
    - Clear all previously saved macros to free up space

## `XRD-imove`

All the TIFFs produced by `XRD-imacro` will be saved in the DATA directory, to avoid clutter, run this command to move them to a personal directory: Documents/XRDpy/\<user-name\>. It is recommended to use your identikey (i.e. abcd1234) as your user-name.

### Arguments
- User-name
    - This will be where your data gets saved. Use your identikey

e.g. `XRD-imove abcd1234`

### Flags
- `--append` or `-A`
    - By default, a new directory will be made in your personal directory every time you use `XRD-imove`, but using this flag will move the data to the last directory created in the same day. This is for if you need to append an om scan with more data.

## `XRD-iplot`

This will find the latest om scan you have performed in your personal directory and plot the data.

### Arguments
- User-name
    - This is where the program will look for your latest scan

e.g. `XRD-iplot abcd1234`

### Flags
- `--animate` or `-A`
    - Produce an animation (movie) playing through each image in the scan.
- `--title TITLE` or `-T TITLE`
    - Give the main plot this TITLE.
- `--date DATE` or `-D DATE`
    - Look for a directory from a previous date from today. Must give the following format YYYYMMDD.
- `--number NUM` or `-N NUM`
    - Look for a directory from a previous scan that wasn't the most recent in today's (or a previous date if `--date` was used). Indexed from 1. If this isn't specified, it will automatically find the most recent.
- `--crop_width` or `-CW`
    - The TIFFs are all cropped around the beam center. This specifies how many pixels horizontally to keep. 13 is the default value
- `--crop_above` or `-CA`
    - How many pixels to keep above the beam center. 100 is the default.
- `--crop_below` or `-CB`
    - How many pixels to keep below the beam center. 20 is the default.
- `--crop_offset` or `-CO`
    - Move the center of the crop in the horizontal direction. By default, this value is 0.
- `--color` or `-C`
    - An integer value you can give to change the color scheme of the plots. If you give too large a value, it will just cycle through the available color schemes.
- `--dir PATH`
    - Override directory search entirely and look in `PATH` for the files. Can give `-D CWD` to get current working directory. Note: `-D` is for `--date`, not `--dir`.
