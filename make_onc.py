"""
Code to generate the ONCdb from Robberto et al. (2013) data
"""
# from astrodbkit import astrodb
from astropy.io import ascii
import astropy.table as at
import numpy as np
from astrodbkit import astrodb

# Load the empty database
db = astrodb.Database('/Users/jfilippazzo/Desktop/orion.db')

# Clear the photometry table
db.modify("DELETE FROM photometry")

def acs_data(file='raw_data/acs.txt'):
    acs = ascii.read(file, data_start=3)
    
    acs = fix_cols(acs)
    
    # Add the ACS objects to the database with their ONCacs number as their source_id
    db.modify("DELETE FROM sources")
    sources = [['ra','dec','publication_shortname']]
    ids = []
    for row in acs:
        id = row['ONCacs']
        if id not in ids:
            sources.append([row['ra'],row['dec'],'Robb13'])
            ids.append(id)
    db.add_data(sources, table='sources')
    
    # Rename some columns
    acs.rename_column('ONCacs', 'source_id')
    acs.rename_column('Obs', 'epoch')
    bands = [c for c in acs.colnames if c.startswith('F')]
    for b in bands:
        acs.rename_column('e_'+b, b+'_unc')
        
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    acs['publication_shortname'] = ['Robb13']*len(acs)
    acs['telescope_id'] = [1]*len(acs)
    acs['instrument_id'] = [1]*len(acs)
    acs['system_id'] = [1]*len(acs)
        
    # Add the photometry to the database
    for b in bands:
        try:
            db.add_data(acs, table='photometry', bands=[b])
        except:
            pass
    
    return acs

def wpc_data(file='raw_data/wpc.txt'):
    wpc = ascii.read(file, data_start=3)
    
    wpc = fix_cols(wpc)
    
    wpc.rename_column('Obs', 'epoch')
    bands = [c for c in wpc.colnames if c.startswith('F')]
    for b in bands:
        wpc.rename_column('e_'+b, b+'_unc')
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    wpc['publication_shortname'] = ['Robb13']*len(wpc)
    wpc['telescope_id'] = [1]*len(wpc)
    wpc['instrument_id'] = [1]*len(wpc)
    wpc['system_id'] = [1]*len(wpc)
    
    # Get the source_id
    source_ids = []
    delta = 0.0001
    for row in wpc:
        matches = db.query("SELECT * FROM sources WHERE (ra BETWEEN ? AND ?) AND (dec BETWEEN ? AND ?)", \
                    (row['ra']-delta,row['ra']+delta,row['dec']-delta,row['dec']+delta), format='table')
        source_id = matches[0]['id']
        source_ids.append(source_id)
    wpc['source_id'] = source_ids
    
    # Add the photometry to the database
    for b in bands:
        try:
            db.add_data(wpc, table='photometry', bands=[b])
        except:
            pass
    
    return wpc

# Read in NICMOS data and fix columns
def nic_data(file='raw_data/nic.txt'):
    nic = ascii.read(file, data_start=3)
    
    nic = fix_cols(nic)
    
    nic.rename_column('Obs', 'epoch')
    bands = [c for c in nic.colnames if c.startswith('F')]
    for b in bands:
        nic.rename_column('e_'+b, b+'_unc')
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    nic['publication_shortname'] = ['Robb13']*len(nic)
    nic['telescope_id'] = [1]*len(nic)
    nic['instrument_id'] = [1]*len(nic)
    nic['system_id'] = [1]*len(nic)
    
    # Get the source_id
    source_ids = []
    delta = 0.0001
    for row in nic:
        matches = db.query("SELECT * FROM sources WHERE (ra BETWEEN ? AND ?) AND (dec BETWEEN ? AND ?)", \
                    (row['ra']-delta,row['ra']+delta,row['dec']-delta,row['dec']+delta), format='table')
        source_id = matches[0]['id']
        source_ids.append(source_id)
    nic['source_id'] = source_ids
    
    # Add the photometry to the database
    for b in bands:
        try:
            db.add_data(nic, table='photometry', bands=[b])
        except:
            pass
        
    return nic
    
def fix_cols(table):
    """
    The error and magnitude columns are reversed in rows where there is no error. This fixes them.
    """
    # Get magnitude columns ad make them floats
    cols = [c for c in table.colnames if c.startswith('F')]
    
    # For each row
    for row in table:
        
        # For each column to correct
        for col in cols:
            
            # Fix misaligned NULL rows
            if isinstance(row[col], np.ma.core.MaskedConstant):
                val = row[col]
                row[col] = row['e_'+col]
                row['e_'+col] = val
                
    # Rename some columns
    table.rename_column('_RAJ2000','ra')
    table.rename_column('_DEJ2000','dec')
    
    return table