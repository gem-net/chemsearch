
# chemsearch

This Python web application is used to create, browse, and search a compounds 
library.


## Manual setup instructions

After you clone or [download](https://github.com/gem-net/chemsearch/archive/master.zip) 
and unzip the repository, you must install the dependencies. The easiest way to 
install the required Python 3 environment and all dependencies is via Conda.
1. If you don't have `conda` on your system, install [Miniconda](https://docs.conda.io/en/latest/miniconda.html).  
2. In the terminal, navigate to the cloned/downloaded folder.
3. Create a conda virtual environment with all dependencies:
   ```bash
   conda env create -n chemsearch -f environment.yaml
   ```
    - or to create an environment precisely mimicking a tested version, specify 
    `environment.lock.yaml` instead of `environment.yaml` above. 
4. Activate the new environment with `conda activate chemsearch`.
5. [Optional] Install the chemsearch package by running `python setup.py develop` (to 
   install in development mode, allowing the code to be modified.)


### Configuration

With default settings, Chemsearch will serve a demo library, but the app can be 
customized by modifying configuration files in the `config` subdirectory.  

An 'env' file provides the primary configuration detail. A demo env file, 
`demo.env` has been provided, which you should update, rename to `.env`, and place 
in the `config` folder.  

The env file is used to specify non-default options for:
- root directory for local compounds library
- app title, for banner
- Drive mode status (on/off). If 'on':
  - Shared Drive ID. You can get this from the end of the Shared Drive URL, e.g. 
   https://drive.google.com/drive/u/1/folders/ABCDEFG123456789
  - User email to use with service account. (A username in your domain that has Drive access.)
- Authentication mode status (on/off). If 'on':
  - Oauth client ID and secret
  - Google Group ID for membership list
  - User email to use with service account (same as Drive mode)
- Flask settings: development vs production mode, port (if serving with `flask`), 
 WSGI script path)

Drive and Authentication modes require a service account credentials JSON file 
named `creds.json` in the config directory. For information on how to generate 
the credentials file, see the [wiki](https://github.com/gem-net/chemsearch/wiki/Get-API-Credentials-File).


## Run the server

The WSGI server executables `flask` and `gunicorn` are available after installing 
the dependencies as above, and either can be used to serve the application. 
`flask` is a 'development' server, best suited to testing code changes, while    
`gunicorn` is more production-ready.
 
**FLASK**:
To serve using `flask`, run the following from the top directory of the repository:
```shell script
FLASK_APP=src/chemsearch/chemsearch.py flask run
```
- This can be simplified to `flask run` if you set `FLASK_APP` in the configuration 
file `config/.env`.

**GUNICORN**:
To serve using `gunicorn`, run the following from the top directory of the repository:
```shell script
gunicorn src.chemsearch.chemsearch:app
```
- Note that you can add command arguments to further customize the application 
server, to modify e.g. listening port and number of worker processes.
