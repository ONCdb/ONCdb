"""
Code to generate the ONCdb from Robberto et al. (2013) data
"""
# from astrodbkit import astrodb
from astropy.io import ascii
import astropy.table as at
import numpy as np
from astrodbkit import astrodb
import pandas as pd

# Photometry flags
# --------------------------------------
# 1:good                            => A
# 2:saturated but recoverd          => B
# 3:saturatated                     => C
# 4:from drizzle/coadded image      => D
# 0:undetected                      => E
# --------------------------------------
flags = {0:'E', 1:'A', 2:'B', 3:'C', 4:'D'}

def generate_ONCdb(cat):
    """
    Generate the database from a list of unique sources and the Robberto+2013 data
    
    Parameters
    ----------
    cat: pandas.DataFrame
         The assembled catalog
    """
    # Make an empty database
    astrodb.create_database('orion.db', 'orion.sql', overwrite=True)
    
    # Load the empty database
    db = astrodb.Database('orion.db')
    
    # Load the source list
    source_list = at.Table(cat.catalog.values, names=cat.catalog.columns)
    
    # Rename some columns
    source_list.rename_column('oncID', 'id')
    source_list.rename_column('oncflag', 'comments')
    source_list.rename_column('ra_corr', 'ra')
    source_list.rename_column('dec_corr', 'dec')
    
    # Populate the SOURCES table (must have 'ra' and 'dec' columns)
    db.add_data(source_list, 'sources')
    
    # Populate the SYSTEMS, INSTRUMENTS, TELESCOPES and PUBLICATIONS tables
    db.add_data([['name'],['Vega']], 'systems', clean_up=False)
    db.add_data([['name','publication_shortname'],['HST','']], 'telescopes', clean_up=False)
    db.add_data([['name','publication_shortname'],['ACS',''],['NICMOS',''],['WFPC2',''],['WFC3','']], 'instruments', clean_up=False)
    db.add_data([['bibcode','shortname','DOI','description'],['2013yCat..22070010R','Robb13','','VizieR Online Data Catalog: HST Treasury Program on the ONC']], 'publications', clean_up=False)
    
    # Populate the other tables with dummy data
    db.query("pragma foreign_keys=OFF")
    db.modify("INSERT INTO spectra (source_id, spectrum) VALUES(0,'/foobar/test.fits')")
    db.modify("INSERT INTO spectral_types (source_id, spectral_type) VALUES(0,0)")
    db.modify("INSERT INTO parallaxes (source_id, parallax) VALUES(0,0)")
    db.query("pragma foreign_keys=ON")
    
    # Add the ACS photometry
    try:
        add_acs_data(db)
    except IOError:
        pass
    
    # # Add the NICMOS photometry
    # try:
    #     add_nicmos_data()
    # except:
    #     pass
    #
    # # Add the WPC photometry
    # try:
    #     add_wpc2_data()
    # except:
    #     pass
        
    return db

def add_acs_data(db, file='raw_data/viz_acs_with_IDs.tsv'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    # Read in the data
    acs = ascii.read(file)
    
    # Rename some columns
    acs.rename_column('Obs', 'epoch')
    acs.rename_column('_RAJ2000', 'ra')
    acs.rename_column('_DEJ2000', 'dec')
    
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
            acs.rename_column('f_'+b, 'flags')
            acs.add_column(at.Column(['ACS_HRC.'+b]*len(acs), 'band'))
            
            # Convert flag integer to string
            acs['flags'] = at.Column([flags[i] for i in acs['flags']], 'flags')
            
            # Move the magnitudes into the correct column
            for row in acs:
                if row['flags'] in ['D','E']:
                    row['magnitude'] = row['magnitude_unc']
                    row['magnitude_unc'] = np.nan
                if row['flags']=='C':
                    row['magnitude_unc'] = np.nan
                if not str(row['magnitude']).strip():
                    row['magnitude'] = np.nan
                if not str(row['magnitude_unc']).strip():
                    row['magnitude_unc'] = np.nan
                    
            # Make sure the magntiudes are floats
            acs['magnitude'] = at.Column(acs['magnitude'], 'magnitude', dtype=float)
            acs['magnitude_unc'] = at.Column(acs['magnitude_unc'], 'magnitude_unc', dtype=float)
            
            # Add the data
            db.query("pragma foreign_keys=OFF")
            db.add_data(acs, table='photometry', clean_up=False)
            db.query("pragma foreign_keys=ON")
            
            # Change the column name back
            acs.rename_column('magnitude', b)
            acs.rename_column('magnitude_unc', 'e_'+b)
            acs.rename_column('flags', 'f_'+b)
            acs.remove_column('band')
        except IOError:
            pass
            
    db.save()

def add_nicmos_data(db='orion.sql', file='raw_data/viz_nicmos_with_IDs.tsv'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    db = astrodb.Database(db)
    
    # Read in the data
    nic = ascii.read(file)
    
    # Rename some columns
    nic.rename_column('Obs', 'epoch')
    nic.rename_column('_RAJ2000', 'ra')
    nic.rename_column('_DEJ2000', 'dec')
    nic.rename_column('oncID', 'id')
    
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
            nic.rename_column('f_'+b, 'flags')
            nic.add_column(at.Column(['NICMOS3.'+b]*len(nic), 'band'))
            
            # Convert flag integer to string
            nic['flags'] = at.Column([flags[i] for i in nic['flags']], 'flags')
            
            # Move the magnitudes into the correct column
            for row in nic:
                if row['flags'] in ['D','E']:
                    row['magnitude'] = row['magnitude_unc']
                    row['magnitude_unc'] = np.nan
                if row['flags']=='C':
                    row['magnitude_unc'] = np.nan
                if not row['magnitude'].strip():
                    row['magnitude'] = np.nan
                if not row['magnitude_unc'].strip():
                    row['magnitude_unc'] = np.nan
                    
            # Make sure the magntiudes are floats
            nic['magnitude'] = at.Column(nic['magnitude'], 'magnitude', dtype=float)
            nic['magnitude_unc'] = at.Column(nic['magnitude_unc'], 'magnitude_unc', dtype=float)
            
            # Add the data
            db.query("pragma foreign_keys=OFF")
            db.add_data(nic, table='photometry', clean_up=False)
            db.query("pragma foreign_keys=ON")
            
            # Change the column name back
            nic.rename_column('magnitude', b)
            nic.rename_column('magnitude_unc', 'e_'+b)
            nic.rename_column('flags', 'f_'+b)
            nic.remove_column('band')
        except IOError:
            pass
            
    db.save()
    db.close()

def add_wpc2_data(db='orion.sql', file='raw_data/viz_wfpc2_with_IDs.tsv'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    db = astrodb.Database(db)
    
    # Read in the data
    wpc = ascii.read(file)
    
    # Rename some columns
    wpc.rename_column('Obs', 'epoch')
    wpc.rename_column('_RAJ2000', 'ra')
    wpc.rename_column('_DEJ2000', 'dec')
    wpc.rename_column('oncID', 'id')
    
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
            wpc.rename_column('f_'+b, 'flags')
            wpc.add_column(at.Column(['WFPC2.'+b.lower()]*len(wpc), 'band'))
            
            # Convert flag integer to string
            wpc['flags'] = at.Column([flags[i] for i in wpc['flags']], 'flags')
            
            # Move the magnitudes into the correct column
            for row in wpc:
                if row['flags'] in ['D','E']:
                    row['magnitude'] = row['magnitude_unc']
                    row['magnitude_unc'] = np.nan
                if row['flags']=='C':
                    row['magnitude_unc'] = np.nan
                if not row['magnitude'].strip():
                    row['magnitude'] = np.nan
                if not row['magnitude_unc'].strip():
                    row['magnitude_unc'] = np.nan
                    
            # Make sure the magntiudes are floats
            wpc['magnitude'] = at.Column(wpc['magnitude'], 'magnitude', dtype=float)
            wpc['magnitude_unc'] = at.Column(wpc['magnitude_unc'], 'magnitude_unc', dtype=float)
            
            # Add the data
            db.query("pragma foreign_keys=OFF")
            db.add_data(wpc, table='photometry', clean_up=False)
            db.query("pragma foreign_keys=ON")
            
            # Change the column name back
            wpc.rename_column('magnitude', b)
            wpc.rename_column('magnitude_unc', 'e_'+b)
            wpc.rename_column('flags', 'f_'+b)
            wpc.remove_column('band')
        except IOError:
            pass
            
    db.save()
    db.close()

