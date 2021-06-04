# chemsearch

Chemsearch is web application used to create, browse, and search a chemical 
compounds library. Powered by the cheminformatics toolkit `rdkit`, Chemsearch 
extracts molecular structure from a collection of MOL files, gathers metadata, 
and provides an interface for exploring the library of compounds. To search 
within the library, Chemsearch offers 'search by similarity' and 'search by 
substructure', taking queries in the form of molecular structure drawings 
(using a molecule editor from Kekule.js) or as SMILES or SMARTS strings. Each
compound has a dedicated page which shows molecular structure, summary data, a 
file listing for the corresponding data folder, and custom information (provided
as a Google Doc, markdown file or other text file).

Chemsearch can work with a Google Workspace / G Suite account for extended functionality.
The app can work with libraries stored in a Google Shared Drive 
(in DRIVE mode) or with a local folder (default mode). Optional access control 
is provided in AUTH mode, with a Google OAuth login page that restricts access 
to members of a specified Google Group.

![FIGURE](src/chemsearch/app/static/schematic.png)

## Contents

<!-- MarkdownTOC autolink="true" -->

- [Getting started](#getting-started)
  - [Quick start, using Docker image](#quick-start-using-docker-image)
  - [Manual installation](#manual-installation)
    - [Installing the dependencies and `chemsearch` executable](#installing-the-dependencies-and-chemsearch-executable)
  - [Starting the server](#starting-the-server)
- [Customization](#customization)
  - [The env file](#the-env-file)
  - [Custom queries: substructure search shortcuts](#custom-queries-substructure-search-shortcuts)
  - [Providing custom data for a compound](#providing-custom-data-for-a-compound)
    - [Example markdown for rich formatting](#example-markdown-for-rich-formatting)
  - [Google Workspace setup](#google-workspace-setup)
- [Building your library](#building-your-library)

<!-- /MarkdownTOC -->


## Getting started

### Quick start, using Docker image

The quickest way to try out Chemsearch is to use the Docker container image on 
[Docker Hub](https://hub.docker.com/r/cgemcci/chemsearch), which has all code and
dependencies bundled in. With Docker [installed and running](https://www.docker.com/get-started) 
on your machine, run the following command in the terminal to download the image 
and browse the demo library with default settings:
```shell script
docker run --rm -p 5000:5000 cgemcci/chemsearch
```
The app will be accessible in your browser at http://localhost:5000/.

Customizing the app requires additional command arguments that bind two local 
folders:
1. a folder of configuration files (see the [Customization](#customization) 
   section below).
2. a folder to contain the local data archive.
```shell script
docker run -p 5000:5000 --rm \
 -v /path/to/config/dir:/app/config \
 -v /path/to/archive/dir:/app/demo_db cgemcci/chemsearch
``` 


### Manual installation

#### Installing the dependencies and `chemsearch` executable

Chemsearch requires Python 3.6 and above, plus a number of third-party Python packages, 
and rdkit. The easiest way to install everything is with Conda. Conda can set up a
suitable Python virtual environment for use with Chemsearch and install rdkit, all
with a single command. 
1. If you don't have `conda` on your system, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).  
2. Clone or [download](https://github.com/gem-net/chemsearch/archive/master.zip) 
  the Chemsearch repository, and unzip it somewhere convenient for long-term storage
  —your Python environment will store a link to the folder as part of installation.
3. In your terminal, navigate to the folder and run the following to create a 
  conda virtual environment with all dependencies:
   ```bash
   conda env create -n chemsearch -f environment.yaml
   ```
    - or to create an environment precisely mimicking a tested version, specify 
    `environment.lock.yaml` instead of `environment.yaml` above. 
4. Activate the new environment with `conda activate chemsearch`. (This environment
  will need to be activated each time you want to run Chemsearch.)
5. Install the chemsearch package by running `python setup.py develop`. This 
    will add a `chemsearch` CLI executable to your PATH that you can use to 
    configure and run your Chemsearch server.
6. (Optional) To bypass app configuration and import a demo library for use with the app, run 
  `chemsearch build`. This demo library will be available for browsing when you 
  start the server.


### Starting the server

The WSGI server executables `gunicorn` and `flask` are available after installing 
the package and dependencies as above and activating the chemsearch virtual environment. 
Either executable can be used to serve the application. `gunicorn` is the recommended, 
more production-ready option, while `flask` is a 'development' server, best suited 
to testing code changes.

 **GUNICORN**:
To serve using `gunicorn`, run the following :
```shell script
gunicorn -w 1 -k gevent chemsearch.chemsearch:app
```
- Note that you can add [command arguments](https://docs.gunicorn.org/en/stable/run.html#commonly-used-arguments) 
to further customize the application server, to modify e.g. listening port and 
number of worker processes. The command above uses a single 'gevent' worker.
 
**FLASK**:

To serve using `flask`, run the following:
```shell script
chemsearch flask run
```
- This can be simplified to `flask run` if you set `FLASK_APP` in the configuration 
file `config/.env`.


## Customization

Chemsearch has a command-line interface (CLI) that will help you configure the 
app. Get started by running `chemsearch setup prompt` in your terminal. This will 
walk you through the creation of the 'env' file that stores the configuration settings. 

Run `chemsearch setup` to see a complete list of setup subcommands, including: 
```
  creds   Copy Google JSON credentials to config folder.
  edit    Open configuration .env file in an editor.
  import  Load variables from specified .env path.
  prompt  Create configuration .env file via command prompt.
  revert  Revert to previous .env file.
  show    Print configuration path and contents.
```

### The env file

The env file is used to specify non-default options via environment variables:
- root directory for local compounds library (LOCAL_DB_PATH)
- app title, for banner (APP_TITLE, default 'Chemsearch')
- DRIVE mode status, for use with Google Shared Drive (USE_DRIVE, default 'off')
- AUTH mode status (on/off), for authentication with Google OAuth (USE_AUTH, default 'off')
- Flask environment mode (development vs production), if serving using 
  `flask` application server (FLASK_ENV, default 'production').
- Fingerprint for similarity matching (SIM_FINGERPRINT, default 'Morgan').
- Coefficient for similarity matching (SIM_COEFFICIENT, default 'Dice').
  
To use Shared Drive and/or Authentication modes, you will need to [set up a 
Google Workspace](#google-workspace-setup) and provide some additional env 
variables. In the Google-specific variable list below, a tick below means the 
variable is required:

variable | DRIVE mode | AUTH mode | example | notes
-------- | ----------------------- | ---------------------- | ------- | -----
CREDENTIALS_AS_USER | ✅ | ✅ | admin@example.org | (1)
SHARED_DRIVE_ID | ✅ |  | 0AbCdEfG_abc123 | (2)
GOOGLE_CLIENT_ID |  | ✅ | 123456789.apps.googleusercontent.com | (3)
GOOGLE_SECRET |  | ✅ | 98765_abcde | (3)
GROUP_KEY |  | ✅ | 012345abcde | (4)

Notes
1. User email in your domain, whose privileges will be mirrored by service 
   account. For DRIVE mode, this user must have access to the Shared Drive. For 
   AUTH mode, the user must have access to the Google Group. This user will be 
   granted admin status.
2. You can get the Drive ID from the end of the Shared Drive URL, e.g. 
   https://drive.google.com/drive/u/1/folders/ABCDEFG123456789
3. Client and secret correspond to OAuth credentials generated at 
   https://console.developers.google.com/
4. The Group key/ID is shown at the end of the group's URL in the G Suite admin 
   panel, e.g. https://admin.google.com/ac/groups/012345abcde


### Custom queries: substructure search shortcuts

The search page can be customized with shortcut buttons that run a specific
substructure search. These buttons are created automatically (at app launch time) 
when a `custom_queries.yaml` file is added to the `config` folder. The `chemsearch` 
CLI can walk you through the creation of this file with the command:
```shell script
chemsearch shortcuts prompt
```

An example specification:
```yaml
# FORMAT EACH SHORTCUT AS <Name for button>: '<SMARTS STRING>'
Aliphatic amines: '[$([NH2][CX4]),$([NH]([CX4])[CX4]),$([NX3]([CX4])([CX4])[CX4])]'
Bicyclic: '[$([*R2]([*R])([*R])([*R]))].[$([*R2]([*R])([*R])([*R]))]'
```

### Providing custom data for a compound

Chemsearch can display 'custom data' on a compound page, extracting text
from (in decreasing order of priority) a Google Doc (in DRIVE mode), markdown, 
or other text file in the compound folder. Rich formatting in Google Docs and 
Markdown files are preserved by conversion to HTML.

Markdown filenames should end in `.md` (e.g. `custom.md`), while plain text files
should end in `.txt`.


#### Example markdown for rich formatting

An example Markdown file, showing syntax for lists, tables, and blockquotes, 
is provided in `demo_db/base/CHEMBL1195529/custom.md`, and replicated below:

```
<h4>Hello!</h4>

This is some custom information for the molecule page.

1. first item in enumerated list
2. second item
3. final list item.

- bullet point.
- another bullet point.

First Header  | Second Header
------------- | -------------
Content Cell with a long long LONG value  | Content Cell
Content Cell  | Content Cell


> Blockquote
>
> Here's a second paragraph.

See our team website at http://gem-net.net.
```



### Google Workspace setup

To use Shared Drive and/or Authentication modes, you will need a 
[Google Workspace](https://workspace.google.com/) / G Suite account. 
The 'G Suite for Education' edition is free for academic institutions, otherwise 
a Business Standard plan (or higher) will be required to create Shared Drives 
(for DRIVE mode). The premium plans (which offer a free trial period) charge 
based on the number of users, but for Chemsearch only one official user is required: 
a user to create a Google Group (for AUTH mode) and a Shared Drive (for DRIVE mode).
The members of the Google Group and Shared Drive can all use free Google accounts.   

To sign up for a Workspace account, you'll need to register a domain 
(e.g. example.com) and verify your ownership of the domain with Google. The 
remaining steps are as follows:
1. Create a [Google Group](https://admin.google.com/ac/groups), with suitable 
  address e.g. 'everyone@example.com', adding members by email address (either free 
  Google-associated accounts or addresses in your domain).
2. (For DRIVE mode) Create a [Shared Drive](https://support.google.com/a/users/answer/9310249) 
  at [https://drive.google.com/], with a suitable name (e.g. 'Compounds Library'). 
  Add the Google Group as a 'member' of the Drive with read/write access. This 
  means that everyone in the Group will have immediate access, and permissions can be 
  granted and revoked by updating the Google Group members list.
3. At http://console.cloud.google.com/, create a new project, called e.g. Chemsearch.
  Here you'll enable Google APIs and create the credentials used to access them. 
  Under APIs & Services > Library, find and enable the following APIs: 
  Google Drive API (for DRIVE), Google People API (for AUTH), Admin SDK API (for AUTH).
4. Follow [instructions to create a service account](https://developers.google.com/identity/protocols/oauth2/service-account)
   with 'domain-wide authority', with a JSON credentials file. In the process, 
   you'll need to specify 'scopes' to use with the account, which are as follows:
    ```
    https://www.googleapis.com/auth/drive.readonly
    https://www.googleapis.com/auth/drive.metadata.readonly
    https://www.googleapis.com/auth/admin.directory.user.readonly
    https://www.googleapis.com/auth/admin.directory.group.member.readonly
    ```
   - only the first two are required for DRIVE mode, while only the last two are 
     required for AUTH mode.
5. Next [follow instructions](https://developers.google.com/identity/protocols/oauth2/web-server#creatingcred) 
  to create OAuth client ID credentials. Under 'Authorized redirect URI', use:
  - `http://localhost:5000/callback/google` (for testing from `localhost`)
  - `https://<yourdomainname>/callback/google` (if you're serving the app at a 
    custom domain)
 
You should now be able to use the Workspace account with Chemsearch. Just place 
the JSON file in the app `config` directory (renaming it creds.json), and update
the env file as described [above](#the-env-file).

## Building your library

Chemsearch is minimally prescriptive in how you structure your library, with 
only the following requirements (which apply to both the 'local' and 'DRIVE' 
versions:
- the top level folder contains 'category' folders — a way to offer basic 
  categorization and filtering in the Chemsearch browser. 
- each compound gets its own folder within a category folder, and the name of 
  this folder determines the compound name.
- Each compound folder must have a MOL file, to store its 2-D structure.

Additional files and folders associated with compounds can be nested within 
compound folders as desired.
