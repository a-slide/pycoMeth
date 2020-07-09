# Installation

## Create a clean virtual environment (optional but recommended)

Ideally, before installation, create a clean **python3.6+** virtual environment to deploy the package.
Earlier version of Python3 should also work but **Python 2 is not supported**.
For example one can use conda or virtualenvwrapper.

With [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/install.html):

```bash
mkvirtualenv pycoMeth -p python3.7
workon pycoMeth
```

With [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -n pycoMeth python=3.7
conda activate pycoMeth
```

You might also want to install [Nanopolish](https://github.com/jts/nanopolish) in the same virtual environment so you can pipe nanopolish output directly into `pycoMeth`

## Dependencies

[Nanopolish 0.10+](https://github.com/jts/nanopolish) is not a direct dependency but is required to generate the files used by several commands from this package

`pycoMeth` relies on a the following robustly maintained third party python libraries:

* `numpy>=1.14.0`
* `scipy>=1.4.1`
* `statsmodels>=0.11.1`
* `pandas>=1.0.3`
* `Jinja2>=2.11.1`
* `plotly>=4.6.0`
* `pyfaidx>=0.5.8`
* `tqdm>=4.45.0`
* `colorlog>=4.1.0`

The correct versions of packages are installed together with the software when using pip

in addition if you want to enable static images export you will also need the following dependencies:

* `plotly-orca==1.2.1`
* `psutil>=5.7.0`
* `requests>=2.24.0`

Those dependencies will only be installed with conda as `orca` is not a python package.

## Option 1: Installation with pip from pypi

Install or upgrade the package with pip from pypi

```bash
pip install pycoMeth
```

You can also update to the **unstable** development version from test.pypi repository

```bash
pip install --index-url https://test.pypi.org/simple/ pycoMeth -U
```

## Option 2: Installation with conda from Anaconda cloud

This option will also install `orca` which is required for static images export

```bash
# First installation
conda install -c aleg -c plotly pycometh
```

You can also get the **unstable** development version from the dev channel

```bash
conda update -c aleg_dev -c plotly pycometh
```

## Option 3: Installation with pip from Github

Or from github to get the last version

```bash
# First installation
pip install git+https://github.com/a-slide/pycoMeth.git

# First installation bleeding edge
pip install git+https://github.com/a-slide/pycoMeth.git@dev

# Update to last version
pip install git+https://github.com/a-slide/pycoMeth.git --upgrade
```

## Option 4: Clone the repository and install locally in develop mode

With this option, the package will be locally installed in *editable* or *develop mode*. This allows the package to be both installed and editable in project form. This is the recommended option if you wish to modify the code and/or participate to the development of the package (see [contribution guidelines](contributing.md)).

```bash
# Clone repo localy
git clone https://github.com/a-slide/pycoMeth.git

# Enter in repo directory
cd pycoMeth

# Make setup.py executable
chmod u+x setup.py

# Install with pip3
pip3 install -e ./
```
