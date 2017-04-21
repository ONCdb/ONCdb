"""
Code to generate the ONCdb from Robberto et al. (2013) data
"""
from astrodbkit import astrodb
from astropy.io import ascii

# Load the empty database
db = astrodb.Database('orion.sql')

# Read ALL the data into one astropy table
table = ascii.read()
