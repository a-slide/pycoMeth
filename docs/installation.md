# Installation

## Create a clean virtual environment

Ideally, before installation, create a clean **python3.6+** virtual environment to deploy the package. **Python 2 is not supported**.
For example one can use conda or virtualenvwrapper.

With [virtualenvwrapper](https://virtualenvwrapper.readthedocs.io/en/latest/install.html):

```bash
mkvirtualenv pycoMeth -p python3.6
```

With [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```bash
conda create -n pycoMeth python=3.6
```

You might also want to install [Nanopolish](https://github.com/jts/nanopolish) in the same virtual environment so you can pipe nanopolish output directly into pycoMeth

## Dependencies

[Nanopolish 0.10+](https://github.com/jts/nanopolish) is not a direct dependency but is required to generate the files used by several commands from this package

Nanocompore relies on a the following robustly maintained third party python libraries:

* numpy>=1.14.0
* tqdm>=4.23.4

The correct versions of packages are installed together with the software when using pip.

## Option 1: Installation with pip from pypi

Install or upgrade the package with pip from pypi

```bash
# First installation
pip install pycoMeth

# Update to last version
pip install pycoMeth --upgrade
```

## Option 2: Installation with pip from Github

Or from github to get the last version

```bash
# First installation
pip install git+https://github.com/a-slide/pycoMeth.git

# First installation bleeding edge
pip install git+https://github.com/a-slide/pycoMeth.git@dev

# Update to last version
pip install git+https://github.com/a-slide/pycoMeth.git --upgrade
```

## Option 3: Clone the repository and install locally in develop mode

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
