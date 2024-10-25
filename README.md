# XRDpy

Tools to assist in the operation and analysis of the WAXS/SAXS system in CU Physics. This software, installed on the computer running the system, can be used to assist in GIWAXS experiments, and can be installed on any computer to assist in data stitching, transformations, and reductions.

# Installation

It is recommended to create a unique Anaconda environment for the installation because the install will generate globally callable commands from the terminal, and it is nice to isolate these within an Anaconda environment. [Detailed installation instructions can be found on the wiki](https://github.com/UCBoulder/XRDpy/wiki/Detailed-Installation-Instructions-(Windows)).

You can simply download the repository, and in an Anaconda terminal, `cd` into the directory and call `python setup.py install` or `pip install -U .` 

# Included CLI and GUI programs

## pyFAI-calib2

`pyFAI-calib2` is a simple GUI included from pyFAI. This can generate PONI files from fititng AgBh calibration data, and it can be used to generate masks.

PONI refers to *point of normal incidence* and is the point on the detector that intersects the shortest line between the sample and the detector. In the orientation where the incident beam is normal to the detector and the sample is brought into the beam path, the PONI is the same as the beam center on the detector.

### Calibration procedure

- In the toolbar, select "Y" icon and select "Orient Y-axis downward"
- In the right panel:
    - Set wavelength to 1.54185 Angstrom
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

## `STITCH`

This will stitch the raw image files from the /DATA/ directory on the Kubuntu control computer and generate .edf files of the data and the flat field so as to preserve counting statistics.

### Arguments

- Experiment type
    - Description of experiment. E.g. WAXS, SAXS, GIWAXS, or GISAXS.
- Stitch Rows
    - Number of rows in stitch (1 or 2)
- Stitch Columns
    - Number of columns in stitch (1, 2, or 3)
- Exposure time
    - Exposure time of each image in seconds (this time gets split between each 'half' image that fill in the gap between pixel banks).

### Flags

These optional flags allow you to add more information the .edf file header or allow you to optionally save the data as .tif files as well.
- `--incident INCIDENT_ANGLE` or `-I INCIDENT_ANGLE`
    - This is the angle of incidence of the beam with the sample; i.e. how many degrees the sample is tipped wrt to the beam.
- `--om OMEGA_ANGLE` or `-O OMEGA_ANGLE`
    - This is the motor position of the omega motor (aka om).
- `--dist DET_DIST` or `-D DET_DIST`
    - The sample detector distance in meters.
- `--tif` or `-T`
    - If you use this flag, .tif files will also be generated.

## `SCAN`

Running this will produce a macro in Documents/XRDpy/Macros to scan om or z in order to tune the angle of incidence for thin film studies. The scan takes approximately 1 minute for every 10 angles (6 seconds / angle).

The program will print the command you need to feed to the SAXS program to run the macro.

Note: Macro names include date and hour, so if you save two macros in the same hour, they will overwrite.

### Arguments

- Motor to scan.
    - Either 'om' for angular scan or 'z' for linear scan.
- Starting angle (in degrees) or position (in mm).
    - The first position in the scan (smallest or most negative).
    - How far below the current position (in degrees or mm) for relative scans.
- End angle (in degrees) or position (in mm).
    - The last angle in the scan (largest or least negative).
    - How far above the current position (in degrees or mm) for relative scan.
- Angle step size (in degrees).
    - how large of a step to take between each position in the scan.

### Flags
- `--relative` or `-R`
    - Makes scan relative to current position.
- `--name TAG_FOR_FILE` or `-N TAG_FOR_FILE`
    - Add an optional note to the file names.
- `--clear` or `-C`
    - Clear out macro files built up in macro folder.

e.g. `SCAN om 0 1 .01`

## `MOVE`

All the TIFFs produced by `SCAN` will be saved in the DATA directory, to avoid clutter, run this command to move them to a personal directory: Documents/XRDpy/\<user-name\>. It is recommended to use your identikey (i.e. abcd1234) as your user-name.

### Arguments
- Motor that was scanned
    - Either "om" or "z".
- User-name
    - This will be where your data gets saved. Use your identikey.
- Position of other motor
    - Use `wm om` or `wm z` in SPEC to find the position.
    - If doing an om-scan, specify `wm z`
    - If doing a z-scan, specify `wm om`

e.g. `MOVE om abcd1234 [OTHER_POSITION]`

### Flags
- `--append` or `-A`
    - By default, a new directory will be made in your personal directory every time you use `XRD-imove`, but using this flag will move the data to the last directory created in the same day. This is for if you need to append an om scan with more data.

## `PLOT`

This will find the latest om scan you have performed in your personal directory and plot the data.

### Arguments
- User-name
    - This is where the program will look for your latest scan

e.g. `PLOT PATH_TO_DATA -C .25,.19,.19`

### Flags
- `--animate` or `-A`
    - Produce an animation (movie) playing through each image in the scan.
- `--title TITLE` or `-T TITLE`
    - Give the main plot this TITLE.
- `--mod -1` or `-M -1`
    - move $z_0$ in om-scan or change cut off in z-scan (in standard deviations of the beam width; default is 4).
- `--crit COMMA,SEPARATED,CRITICAL,ANGLES` or `-C COMMA,SEPARATED,CRITICAL,ANGLES`
    - comma separated critical angles. Double up on numbers for waveguiding in that layer.
- `--z` or `-Z`
    - Default is an om-specular-scan. This will change to z-specular-scan.
- `--range 1.0` or `-R 1.0`
    - The angular range for an om-specular-scan. This will change where the plot is cut off at the top.
- `--beamwidth HORIZONTAL_BEAMSIZE` or `-B HORIZONTAL_BEAMSIZE`
    - Sets the horizontal crop of the beam. Should be roughly the same size as the horizontal shutter.
- `--save IMAGE_DPI` or `-S IMAGE_DPI`
    - Save the plot.