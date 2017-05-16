"""
Code to generate the ONCdb from Robberto et al. (2013) data
"""
# from astrodbkit import astrodb
from astropy.io import ascii

# Load the empty database
# db = astrodb.Database('orion.db')

# Read in NICMOS data and fix columns
# nicmos = ascii.read('raw_data/onc_nicmos.CDS.dat.txt', delimiter='\t')
# for row in nicmos:
# 	if row['m_110']>99 and row['dm_110']<99:
# 		row['m_110'] = row['dm_110']
# 		row['dm_110'] = 99.9999
# 	if row['m_160']>99 and row['dm_160']<99:
# 		row['m_160'] = row['dm_160']
# 		row['dm_160'] = 99.9999

def f():
	acs = ascii.read('raw_data/onc_acs.CDS.dat.txt', delimiter='\t')
	print(acs)
