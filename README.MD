OCEAN - Optimized Cross rEActivity estimatioN
--------

So what is OCEAN?
--------
OCEAN is a web-tool for target-prediction of chemical structures which uses ChEMBL as datasource.
For more details please feel free to read the [paper](http://www.paper.foo).

Software dependencies
--------

We recommend to use [myChembl20-VM](http://chembl.blogspot.de/2015/07/mychembl-20-has-landed.html) which has almost everything you need to run OCEAN.

Nevertheless here are the core-depencies:

* `RDKit (2015.03.1)` http://www.rdkit.org/
* `Django (1.7.1)`
* `numpy (1.9.0)`
* `sciPy (0.15.1)`
* `your preferred database-server (Oracle or PostgreSQL, we didn't test MySQL or sqlite but they are supported by Django as well)`
* `python module for your database-server (cx_oracle for Oracle, psycopg2 for PostgreSQL)`

Optional dependencies
--------

* `matplotlib (1.3.0)`
* `pillow (2.6.1)`
* `images2gif (1.0.1)`

Installation using myCHEMBL-VM
--------
Prepare and activate a virtual environment for OCEAN:
```
# create a virtualEnv for OCEAN
virtualenv --system-site-packages ocean_env

# activate ocean_env
source ocean_env/bin/activate

# install Django 1.7.1
pip install Django==1.7.1
```
Get a clone of OCEAN from the Repo:
```
git clone https://www.github.com/rdkit/OCEAN.git
cd OCEAN
```
Create a PostgreSQL-User, the database for OCEAN and fill the OCEAN-DB:
```
# create PostgreSQL-User for OCEAN
# password 'ocean_pw' should be used when asked,
# according to database-entry in settings.py
createuser -P -s -d -r -e ocean_user

# create PostgreSQL-Database for OCEAN
createdb --owner=ocean_user ocean
```
If you want to use ChEMBL17-data for OCEAN you should:


>change line `CHEMBL_VERSION = "chembl_%d" % (20)` in file `ocean/settings.py` to `CHEMBL_VERSION = "chembl_%d" % (17)`

```
createdb --owner=ocean_user chembl_17

# create and populate OCEAN-DB (duration <2 min)
psql chembl_17 < postgres_createOceanDB.sql
./importOcean17DB.sh

# Chembl17 is not a part of myChembl20-VM, so we need to create and populate
# a few tables (duration <10 sec)
psql chembl_17 < postgres_createTargetDict.sql
./importChembl17Tables.sh
```
For ChEMBL20-data you should:
```
psql chembl_20 < postgres_createOceanDB.sql
./importOcean20DB.sh
```
It's time to start OCEAN:
```
# On first start, you have to sync the django-server to create all missing tables.
# When asked, you should create a superuser and set a password.
# This is necessary to have access to the django-admin portal of OCEAN.
python manage.py syncdb

# Start OCEAN
python manage.py runserver 0.0.0.0:8080
```
The next steps are necessary to create the fingerprints for all structures:
1. Login to ocean-admin portal http://0.0.0.0:8080/admin

2. Select *CHEMBL* and *6* (section a), press *createFPs* to create Morgan-Fingerprints for all entries. This may take a few minutes. You should see a message at the bottom of the site when the process is done.

The last step is calculating the statistical parameters for OCEAN's scoring-function:
1. In the Django-Admin portal select *CHEMBL*, *6* and activate Checkbox *rebuild Pairs* (section b), press *recalc Statistics for Fingerprint and Datasource* to calculate statistic-parameters for the scoring-function. This will also take a few minutes. Wait until you see a new message at the bottom of the site.

If you have a MarvinJS licence, you can copy the MarvinJS-folder to `ocean/media/MarvinJS` and change the line `url(r'^$', main_smiles, name='home')` in `ocean/urls.py` to `url(r'^$', main_marvin, name='home')` to use MarvinJS as sketcher.

Now lets try it:

1. Open ocean-website http://0.0.0.0:8080

2. Write a smiles-string into the *smiles*-field or sketch a molecule (if you have MarvinJS) and press Button *Start OCEAN-Search*. You should see the results after a few seconds.

Questions or Problems?
--------
If you have any questions or problems with OCEAN please do not hesitate to contact [Paul Czodrowski](mailto:Paul.Czodrowski@merckgroup.com) or [Wolf-Guido Bolick](mailto:wolf-guido.bolick@gmx.de).
