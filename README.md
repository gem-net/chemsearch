
# chemsearch

This python package is used to create a local Molecules repository, 
with all the information necessary for performing similarity and substructure 
matching across the collected molecules.

The package assumes that all molecules are stored in a Shared Drive, whose top 
level directory contains category folders (one folder per category, e.g. `aramids`),
and each category folder contains molecule folders (one folder per molecule).


## Setup instructions

First, download the package: https://github.com/gem-net/chemsearch/archive/master.zip. 
Next, unzip it and move the unzipped folder somewhere you're happy storing it, 
e.g. your Documents folder. 

### Python package dependencies

The code in this repository assumes a Python 3 environment. The easiest way to 
install all required dependencies is via Conda:
1. download and install Miniconda (if you don't have conda already): 
  https://docs.conda.io/en/latest/miniconda.html
2. In your shell (Terminal app on Mac), navigate to the unzipped `chemsearch` package.
3. Create a conda virtual environment with all dependencies:
   ```bash
   conda env create -n chemsearch -f environment.yaml
   ```
    - or to create an environment precisely mimicking the tested version, specify 
    `environment.lock.yaml` instead of `environment.yaml` above. 
4. Activate the new environment with `conda activate chemsearch`.
5. Install the chemsearch package by running `python setup.py develop` (to 
   install in development mode, allowing the code to be modified.)


### Configure using env file

An env file is required to provide configuration detail that is not suitable for hard-coding. 
A demo env file, `demo.env` has been provided, which you should update and rename to `.env`. 

The env file is used to specify: 
- the location of a service account credentials JSON file. For more on this see 
  the [wiki](https://github.com/gem-net/chemsearch/wiki/Get-API-Credentials-File)
- the Shared Drive alphanumeric ID 
   - you can get this from the end of the Shared Drive URL, e.g. 
   https://drive.google.com/drive/u/1/folders/ABCDEFG123456789
- a username in your domain that has Drive access.

`.env_demo`:
```bash
SERVICE_ACCOUNT_FILE=/path/to/credentials/file.json
SHARED_DRIVE_ID=ABCDEFG123456789
CREDENTIALS_AS_USER=someuser@gem-net.net
```

### (Optional) Run demo using notebooks

- Make sure the chemsearch conda environment is active (`conda activate chemsearch`).
- In the chemsearch directory, run `jupyter lab`, which should open a Jupyter Lab 
  window in your browser.
- Double click on `notebooks` in the navigation pane, then open one of the 
  notebooks, e.g. `step1_scan_gdrive.ipynb`.
- Run the notebook one cell at a time using the '▶' button or with Shift-Return.


## Build your local repository

When you ran `python seutp.py develop` above, it made an executable command 
called `build_molecules` available in your shell. This executable creates a local
folder of molecule files in whatever folder you tell it to use. To save this 
local 'repository' to ~/Desktop/local_db, for example, you can run:
```bash
conda activate chemsearch  # to make sure you're in the correct coda environment
build_molecules --verbose ~/Desktop/local_db
```
