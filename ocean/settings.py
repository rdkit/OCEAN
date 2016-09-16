# Django settings for ocean project.

import os,sys
from platform import node
from os.path import dirname
from enum import Enum
from rdkit import Chem
import SimilarityMaps as SM
from rdkit.Chem import MACCSkeys


PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
PROJECT_NAME = os.path.split(PROJECT_ROOT)[-1]

SECRET_KEY = 'wq6sbq%h0quizegiy@-s^g24_a+atmlfl%+mca)+n==^o58v93'

DEBUG = True
TEMPLATE_DEBUG = DEBUG

ADMINS = (
)

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',  # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'ocean',  # Or path to database file if using sqlite3.
#         # The following settings are not used with sqlite3:
        'USER': 'ocean_user',
        'PASSWORD': 'ocean_pw',
#         'HOST': '',                      # Empty for localhost through domain sockets or '127.0.0.1' for localhost through TCP.
#         # 'PORT': '',                      # Set to empty string for default.
    },
    'chembl_17':{
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
    # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'chembl_17',  # Or path to database file if using sqlite3.
        #         # The following settings are not used with sqlite3:
        'USER': 'ocean_user',
        'PASSWORD': 'ocean_pw',
    },
    'chembl_20':{
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
    # Add 'postgresql_psycopg2', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': 'chembl_20',  # Or path to database file if using sqlite3.
        #         # The following settings are not used with sqlite3:
        'USER': 'ocean_user',
        'PASSWORD': 'ocean_pw',
    }
}

FP_METHODS = [Chem.RDKFingerprint,
              Chem.LayeredFingerprint,
              Chem.PatternFingerprint,
              MACCSkeys.GenMACCSKeys,
              SM.GetAPFingerprint,
              SM.GetTTFingerprint,
              SM.GetMorganFingerprint,
              ]

SCORING_PARAMS = {
    'CHEMBL':
        {'FP': 6,
         'THRESHOLD': 0.30}}

DATASOURCES = Enum('DataSource', 'CHEMBL')

DATASOURCE_LINK_TARGET = {'CHEMBL': "https://www.ebi.ac.uk/chembl/target/inspect/{0}"}
DATASOURCE_LINK_COMPOUND = {'CHEMBL': "https://www.ebi.ac.uk/chembl/compound/inspect/{0}"}

CALC_OCEAN_PARAMETER_REPEATS = 10
CALC_OCEAN_PARAMETER_START = 10
CALC_OCEAN_PARAMETER_END = 316
CALC_OCEAN_PARAMETER_STEPS = 1
CALC_OCEAN_PARAMETER_THRESH_START = 0.01
CALC_OCEAN_PARAMETER_THRESH_END = 1.00
CALC_OCEAN_PARAMETER_THRESH_STEPS = 0.01

OCEAN_DB_TABLE = 'ocean_db'

# change this to False to accept TC 1.0-Hits for Scoring, for validation it is useful to ignore exact matches,
# in most cases these are the same structures (duplicates) with other id
VALIDATING_PROCESS = False

ADMIN_HELPTEXT = """
On a fresh server you should:<br><br>
\ta)\tcreate FPs (especially for the default FP's of Datasets, which is defined in SCORING_PARAMS (in custom_settings.py or settings.py)<br><br>
\tb)\trecalc Statistics for your Datasources and Fingerprints (used threshold is defined by SCORING_PARAMS-Variable (in custom_settings.py or settings.py))
"""

# this should be a path/file to a line-seperated list of molregnos/compound_ids
# this variable is replaced by a set of its entries later
# or set it to None
drop_compounds = None

CSV_HEADER_ITEMS = ["ID",
                    "DRUG_NAME",
                    "MOLECULE_ID",
                    "CANONICAL_SMILES",
                    "MECHANISM_OF_ACTION",
                    "TARGET_NAME",
                    "TARGET_ID",
                    "CLASSIFICATION",
                    "RESULT_VALUE"]
CSV_SEP = ";"

# Masquerading of CSV-Entries, None should be possible but not tested
CSV_MASQ = '"'

csv_header_tmp = []
for header_item in CSV_HEADER_ITEMS:
    csv_header_tmp.append('{0}{1}{0}'.format(CSV_MASQ,header_item) if CSV_MASQ else header_item)
CSV_HEADER = CSV_SEP.join(csv_header_tmp)
del(csv_header_tmp)

CMPD_COUNT_CUTOFF = 10
CMPD_NM_CUTOFF = 10000
PARALLEL_PROCESSES =  4 # 45 # for server, 4 for workstation
ORGANISM = "Homo sapiens"

CHEMBL_VERSION = "chembl_%d" % (20)

SUBDOMAIN = 'ocean'
HTTP_PORT = '8081'
URL_BASE = '%s.my_subdomain.com:%s' % (SUBDOMAIN, HTTP_PORT)
# Hosts/domain names that are valid for this site, required if DEBUG is False
# See https://docs.djangoproject.com/en/1.5/ref/settings/#allowed-hosts
ALLOWED_HOSTS = ['localhost', '127.0.0.1', URL_BASE.split(":")[0]]
HTTP_HOSTNAME = node()

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# In a Windows environment this must be set to your system time zone.
TIME_ZONE = 'Europe/Berlin'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-us'
LANGUAGES = (('en', 'English'),
             ('de', 'German'))

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = True

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale.
USE_L10N = True

# If you set this to False, Django will not use timezone-aware datetimes.
USE_TZ = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/var/www/example.com/media/"
# MEDIA_ROOT = ''
MEDIA_ROOT = os.path.join(PROJECT_ROOT, '../media')
# MEDIA_ROOT = '/home/m207409/python/djangoprojects/mmp/media'
# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash.
# Examples: "http://example.com/media/", "http://media.example.com/"
MEDIA_URL = '/ocean_media/'

# Absolute path to the directory static files should be collected to.
# Don't put anything in this directory yourself, store your static files
# in apps' "static/" subdirectories and in STATICFILES_DIRS.
# Example: "/var/www/example.com/static/"
STATIC_ROOT = ''

# URL prefix for static files.
# Example: "http://example.com/static/", "http://static.example.com/"
STATIC_URL = '/static/'

# Additional locations of static files
STATICFILES_DIRS = (
    # Put strings here, like "/home/html/static" or "C:/www/django/static".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

# List of finder classes that know how to find static files in
# various locations.
STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
    #    'django.contrib.staticfiles.finders.DefaultStorageFinder',
)



# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
    #     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    # Uncomment the next line for simple clickjacking protection:
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'ocean.urls'

# Python dotted path to the WSGI application used by Django's runserver.
WSGI_APPLICATION = 'ocean.wsgi.application'

TEMPLATE_DIRS = (
    os.path.join(PROJECT_ROOT, '../templates'),
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.sites',
    'django.contrib.admin',
    'django.contrib.admindocs',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'ocean',
)

# A sample logging configuration. The only tangible logging
# performed by this configuration is to send an email to
# the site admins on every HTTP 500 error when DEBUG=False.
# See http://docs.djangoproject.com/en/dev/topics/logging for
# more details on how to customize your logging configuration.
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse'
        }
    },
    'handlers': {
        'mail_admins': {
            'level': 'ERROR',
            'filters': ['require_debug_false'],
            'class': 'django.utils.log.AdminEmailHandler'
        }
    },
    'loggers': {
        'django.request': {
            'handlers': ['mail_admins'],
            'level': 'ERROR',
            'propagate': True,
        },
    }
}

TEST_RUNNER = 'django.test.runner.DiscoverRunner'

if drop_compounds is not None:
    if os.path.isfile(drop_compounds):
        with open(drop_compounds, 'r') as fh:
            drop_compounds = map(lambda x: int(x.rstrip()) if not x.startswith('#') else -1, fh.readlines())
            drop_compounds = set(drop_compounds)
else:
    drop_compounds = set()

if os.path.exists('ocean/custom_settings.py'):
    import custom_settings
    msee = {k:v for k,v in custom_settings.__dict__.items() if not k.startswith('__')}
    current_locals = locals()
    for entry,value in msee.items():
        current_locals.update({entry:value})
        print >> sys.stderr, "monkey patch custom setting", entry, value
