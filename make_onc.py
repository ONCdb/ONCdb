"""
Code to generate the ONCdb from Robberto et al. (2013) data
"""
# from astrodbkit import astrodb
from astropy.io import ascii
import astropy.table as at
import numpy as np
from astrodbkit import astrodb

def generate_ONCdb():
    """
    Generate the database from a list of unique sources and the Robberto+2013 data
    """
    # Make an empty database
    astrodb.create_database('orion.db', 'orion.sql', overwrite=True)

    # Load the empty database
    db = astrodb.Database('orion.db')
    
    # Load the source list
    source_list = ascii.read('raw_data/generated_source_list.txt', delimiter='\t')[:20]
    
    # Populate the SOURCES table (must have 'ra' and 'dec' columns)
    db.add_data(source_list, 'sources')
    
    # Populate the SYSTEMS, INSTRUMENTS, TELESCOPES and PUBLICATIONS tables
    db.add_data([['name'],['Vega']], 'systems')
    db.add_data([['name','publication_shortname'],['HST','']], 'telescopes')
    db.add_data([['name','publication_shortname'],['ACS',''],['NICMOS',''],['WFPC2',''],['WFC3','']], 'instruments')
    db.add_data([['bibcode','shortname','DOI','description'],['2013yCat..22070010R','Robb13','','VizieR Online Data Catalog: HST Treasury Program on the ONC']], 'publications')
    
    db.save()
    
    # Add the ACS photometry
    # add_acs_data()
    #
    # # Add the NICMOS photometry
    # add_nicmos_data()
    #
    # # Add the WPC photometry
    # add_wpc2_data()
    
    return

def add_acs_data(db='orion.db', file='raw_data/acs.txt'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    db = astrodb.Database(db)
    
    # Read in the data
    acs = ascii.read(file, data_start=3)
    
    # Rename some columns
    acs.rename_column('Obs', 'epoch')
    acs.rename_column('_RAJ2000', 'ra')
    acs.rename_column('_DEJ2000', 'dec')
    
    # Find the source_id for each row
    IDS = []
    for ra,dec in np.array(acs[['ra','dec']]):
        try:
            id = int(db.query("select id from sources where ra={} and dec={}".format(ra,dec), fmt='table')['id'][0])
        except:
            id = 99999999
        IDS.append(id)
            
    # Add the source_id column
    acs.add_column(at.Column(IDS, name='source_id'), index=0)
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    acs['publication_shortname'] = ['Robb13']*len(acs)
    acs['telescope_id'] = [1]*len(acs)
    acs['instrument_id'] = [1]*len(acs)
    acs['system_id'] = [1]*len(acs)

    # Add the photometry to the database one band at a time
    # Rename the columns to match svo_filters
    bands = [c for c in acs.colnames if c.startswith('F')]
    for b in bands:
        try:
            # Change the column names to add the band
            acs.rename_column(b, 'magnitude')
            acs.rename_column('e_'+b, 'magnitude_unc')
            acs.rename_column('f_'+b, 'flag')
            acs.add_column(at.Column(['ACS_HRC.'+b]*len(acs), 'band'))
            
            # Add the data
            db.query("pragma foreign_keys=OFF")
            db.add_data(acs, table='photometry')
            db.query("pragma foreign_keys=ON")
            
            # Change the column name back
            acs.rename_column('magnitude', b)
            acs.rename_column('magnitude_unc', 'e_'+b)
            acs.rename_column('flag', 'f_'+b)
            acs.remove_column('band')
        except IOError:
            pass
            
    db.save()
    db.close()
    
    return

def add_nicmos_data(db='orion.db', file='raw_data/nic.txt'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    db = astrodb.Database(db)
    
    # Read in the data
    nic = ascii.read(file, data_start=3)
    
    # Rename some columns
    # nic.rename_column('Obs', 'epoch')
    nic.rename_column('_RAJ2000', 'ra')
    nic.rename_column('_DEJ2000', 'dec')
    
    # Find the source_id for each row
    IDS = []
    for ra,dec in np.array(nic[['ra','dec']]):
        try:
            id = int(db.query("select id from sources where ra={} and dec={}".format(ra,dec), fmt='table')['id'][0])
        except:
            id = 99999999
        IDS.append(id)
            
    # Add the source_id column
    nic.add_column(at.Column(IDS, name='source_id'), index=0)
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    nic['publication_shortname'] = ['Robb13']*len(nic)
    nic['telescope_id'] = [1]*len(nic)
    nic['instrument_id'] = [2]*len(nic)
    nic['system_id'] = [1]*len(nic)

    # Add the photometry to the database one band at a time
    # Rename the columns to match svo_filters
    bands = [c for c in nic.colnames if c.startswith('F')]
    for b in bands:
        try:
            # Change the column names to add the band
            nic.rename_column(b, 'magnitude')
            nic.rename_column('e_'+b, 'magnitude_unc')
            nic.rename_column('f_'+b, 'flag')
            nic.add_column(at.Column(['NICMOS3.'+b]*len(nic), 'band'))
            
            # Add the data
            db.query("pragma foreign_keys=OFF")
            db.add_data(nic, table='photometry')
            db.query("pragma foreign_keys=ON")
            
            # Change the column name back
            nic.rename_column('magnitude', b)
            nic.rename_column('magnitude_unc', 'e_'+b)
            nic.rename_column('flag', 'f_'+b)
            nic.remove_column('band')
        except IOError:
            pass
            
    db.save()
    db.close()
    
    return

def add_wpc2_data(db='orion.db', file='raw_data/wpc.txt'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    db = astrodb.Database(db)
    
    # Read in the data
    wpc = ascii.read(file, data_start=3)
    
    # Rename some columns
    # wpc.rename_column('Obs', 'epoch')
    wpc.rename_column('_RAJ2000', 'ra')
    wpc.rename_column('_DEJ2000', 'dec')
    
    # Find the source_id for each row
    IDS = []
    for ra,dec in np.array(wpc[['ra','dec']]):
        try:
            id = int(db.query("select id from sources where ra={} and dec={}".format(ra,dec), fmt='table')['id'][0])
        except:
            id = 99999999
        IDS.append(id)
            
    # Add the source_id column
    wpc.add_column(at.Column(IDS, name='source_id'), index=0)
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    wpc['publication_shortname'] = ['Robb13']*len(wpc)
    wpc['telescope_id'] = [1]*len(wpc)
    wpc['instrument_id'] = [2]*len(wpc)
    wpc['system_id'] = [1]*len(wpc)

    # Add the photometry to the database one band at a time
    # Rename the columns to match svo_filters
    bands = [c for c in wpc.colnames if c.startswith('F')]
    for b in bands:
        try:
            # Change the column names to add the band
            wpc.rename_column(b, 'magnitude')
            wpc.rename_column('e_'+b, 'magnitude_unc')
            wpc.rename_column('f_'+b, 'flag')
            wpc.add_column(at.Column(['WFPC2.'+b]*len(wpc), 'band'))
            
            # Add the data
            db.query("pragma foreign_keys=OFF")
            db.add_data(wpc, table='photometry')
            db.query("pragma foreign_keys=ON")
            
            # Change the column name back
            wpc.rename_column('magnitude', b)
            wpc.rename_column('magnitude_unc', 'e_'+b)
            wpc.rename_column('flag', 'f_'+b)
            wpc.remove_column('band')
        except IOError:
            pass
            
    db.save()
    db.close()
    
    return