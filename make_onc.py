"""
Code to generate the ONCdb from Robberto et al. (2013) data
"""
# from astrodbkit import astrodb
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.table as at
import astropy.units as q
import numpy as np
from astrodbkit import astrodb, astrocat
import pandas as pd
path = '/Users/jfilippazzo/Documents/Modules/ONCdb/'
# Photometry flags
# --------------------------------------
# 1:good                            => A
# 2:saturated but recoverd          => B
# 3:saturatated                     => C
# 4:from drizzle/coadded image      => D
# 0:undetected                      => E
# --------------------------------------
flags = {0:'E', 1:'A', 2:'B', 3:'C', 4:'D'}

def ONC_catalogs_to_database(ra=83.81775*q.deg, dec=-5.38788889*q.deg, radius=0.0001, count=-1):
    """
    Generated the SQL database from the input catalogs
    """
    # Empty instance
    onc, db = astrocat.Catalog(), None
    
    # ACS catalog from Robberto+2013
    onc.Vizier_query('J/ApJS/207/10/table5', 'ACS', ra, dec, 10*q.deg, group=False)
    
    # ACS catalog from Robberto+2013
    onc.Vizier_query('J/ApJS/207/10/table6', 'WFPC2', ra, dec, 10*q.deg, group=False)

    # ACS catalog from Robberto+2013
    onc.Vizier_query('J/ApJS/207/10/table7', 'NICMOS', ra, dec, 10*q.deg, group=False)
    
    # Group sources
    onc.group_sources(radius)
    
    # Get the radius from the ONC center which includes all Robberto+2013 sources
    center = coord.SkyCoord(ra=ra, dec=dec, frame='icrs')
    radec = coord.SkyCoord(ra=onc.catalog['ra'], dec=onc.catalog['dec'], unit=(q.deg, q.deg), frame='icrs')
    onc_radius = np.max(radec.separation(center)).value*q.deg
    
    # Get 2MASS
    onc.Vizier_query('II/246/out', 'TMASS', ra, dec, onc_radius, group=False)

    # Get GAIA (Only 19 sources! Placeholder for DR2!)
    onc.Vizier_query('I/337/tgas', 'TGAS', ra, dec, onc_radius, ra_col='_RA.icrs', dec_col='_DE.icrs', group=False)

    # Get ALLWISE
    onc.Vizier_query('II/328/allwise', 'ALLWISE', ra, dec, onc_radius, group=False)
   
    # Get spectral types from Hillenbrand+2013
    onc.Vizier_query('J/AJ/146/85/table2', 'Hill13', ra, dec, onc_radius, group=False)
    
    # Get the LAMOST spectra
    onc.Vizier_query('V/149/dr2', 'LAMOST', ra, dec, onc_radius, group=False, column_filters={"Class":"=STAR", "objType":"=Star"})
    
    # Get SDSS spectra
    onc.Vizier_query('V/147/sdss12', 'SDSS', ra, dec, onc_radius, ra_col='RA_ICRS', dec_col='DE_ICRS', column_filters={"class":"=6"}, group=False)
    
    onc.group_sources(radius)

    # Generate SQL database
    db = generate_ONCdb(onc)
        
    return onc, db

def generate_ONCdb(cat):
    """
    Generate the database from a list of unique sources and the Robberto+2013 data
    
    Parameters
    ----------
    cat: astrodbkit.astrocat.Catalog
         The assembled catalog
    """
    # Make an empty database
    astrodb.create_database(path+'orion.db', path+'orion.sql', overwrite=True)
    
    # Load the empty database
    db = astrodb.Database(path+'orion.db')
    
    # Load the source list
    source_list = at.Table(cat.catalog.values, names=cat.catalog.columns)

    # Rename some columns
    source_list.rename_column('flag', 'comments')

    # Populate the SOURCES table (must have 'ra' and 'dec' columns)
    db.add_data(source_list, 'sources')

    # Populate the SYSTEMS, INSTRUMENTS, TELESCOPES and PUBLICATIONS tables
    db.add_data([['name'],['Vega'],['AB']], 'systems', clean_up=False)
    db.add_data([['name','publication_shortname'],['HST',''],['2MASS',''],['WISE',''],['GAIA',''],['LAMOST',''],['SDSS','']], 'telescopes', clean_up=False)
    db.add_data([['name','publication_shortname'],['ACS',''],['NICMOS',''],['WFPC2',''],['WFC3',''],['2MASS',''],['WISE',''],['GAIA',''],['LAMOST',''],['SDSS','']], 'instruments', clean_up=False)
    db.add_data([['bibcode','shortname','DOI','description'],\
                 ['2013yCat..22070010R','Robb13','','VizieR Online Data Catalog: HST Treasury Program on the ONC'],\
                 ['2003yCat.2246....0C','Cutr03','','VizieR Online Data Catalog: 2MASS All-Sky Catalog of Point Sources'],\
                 ['2014yCat.2328....0C','Cutr13','','VizieR Online Data Catalog: AllWISE Data Release '],\
                 ['2016A&A...595A...1G','Gaia16','','Gaia DR1'],\
                 ['2016yCat.5149....0L','Luo_16','','The second data release (DR2) of the LAMOST regular survey'],\
                 ['2015ApJS..219...12A','Alam15','','The Eleventh and Twelfth Data Releases of the Sloan Digital Sky Survey: Final Data from SDSS-III']\
                ], 'publications', clean_up=False)

    # Populate the other tables with dummy data
    db.query("pragma foreign_keys=OFF")
    db.modify("INSERT INTO spectra (source_id, spectrum) VALUES(0,'/foobar/test.fits')")
    db.modify("INSERT INTO spectral_types (source_id, spectral_type) VALUES(0,0)")
    db.modify("INSERT INTO parallaxes (source_id, parallax) VALUES(0,0)")
    db.query("pragma foreign_keys=ON")

    # Add the ACS photometry
    try:
        add_acs_data(db, cat.ACS)
    except:
        print('No ACS data added')

    # Add the NICMOS photometry
    try:
        add_nicmos_data(db, cat.NICMOS)
    except:
        print('No NICMOS data added')

    # Add the WPC photometry
    try:
        add_wpc2_data(db, cat.WFPC2)
    except:
        print('No WFPC2 data added')

    # Add the 2MASS photometry
    try:
        add_2MASS_data(db, cat.TMASS)
    except:
        print('No 2MASS data added')

    # Add the ALLWISE photometry
    try:
        add_WISE_data(db, cat.ALLWISE)
    except:
        print('No WISE data added')

    # Add the GAIA photometry and parallaxes
    try:
        add_TGAS_data(db, cat.TGAS)
    except:
        print('No TGAS data added')

    # Add the spectral types
    try:
        add_Hill13_data(db, cat.Hill13)
    except:
        print('No Hill13 data added')
        
    # Add LAMOST spectra
    try:
        add_LAMOST_data(db, cat.LAMOST)
    except:
        print('No LAMOST data added')
        
    # Add SDSS photometry
    try:
        add_SDSS_data(db, cat.SDSS)
    except:
        print('No SDSS data added')
        
    return db

def add_LAMOST_data(db, cat):
    lamo = at.Table.from_pandas(cat)
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    lamo['publication_shortname'] = ['Luo_15']*len(lamo)
    lamo['telescope_id'] = [5]*len(lamo)
    lamo['instrument_id'] = [8]*len(lamo)
    lamo['flux_units'] = ['Wm-2um-1']*len(lamo)
    lamo['wavelength_units'] = ['A']*len(lamo)
    lamo.rename_column('Obs.Date', 'obs_date')
    lamo['spectrum'] = ['http://cdsarc.u-strasbg.fr/ftp/cats/V/146/LAMOST/fits/{0}/spec-{1}-{0}%5Fsp{2}-{3}.fits'.format(lamo[n]['PlanId'],lamo[n]['LMJD'],lamo[n]['spId'],lamo[n]['FiberId']) for n in range(len(lamo))]
    db.query("pragma foreign_keys=OFF")
    db.add_data(lamo, table='spectra', clean_up=False)
    db.query("pragma foreign_keys=ON")
            
    # Collect the spectral types
    spts = lamo['SubClass']
    
    # Convert the spectral types to integers
    spts = [specType(s) if isinstance(s, str) else [np.nan,''] for s in spts]
    typ, lc = np.array(spts).T
    
    # Add to the table
    lamo['spectral_type'] = typ
    lamo['luminosity_class'] = lc
    
    # Add the data
    db.query("pragma foreign_keys=OFF")
    db.add_data(lamo, table='spectral_types', clean_up=False)
    db.query("pragma foreign_keys=ON")
            
    db.save()
    
def add_TGAS_data(db, cat):
    """
    Read in the GAIA photometry and astrometry
    """
    # Read in the data
    gaia = at.Table.from_pandas(cat)
    
    # Rename some columns
    gaia.rename_column('Epoch', 'epoch')
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    gaia['publication_shortname'] = ['Gaia16']*len(gaia)
    gaia['telescope_id'] = [4]*len(gaia)
    gaia['instrument_id'] = [7]*len(gaia)
    gaia['system_id'] = [1]*len(gaia)
    
    # Rename and make sure the magntiudes are floats
    # gaia.rename_column('__Gmag_', 'magnitude')
    # gaia.rename_column('e__Gmag_', 'magnitude_unc')
    # gaia['band'] = ['G']*len(gaia)
    # gaia['magnitude'] = at.Column(gaia['magnitude'], 'magnitude', dtype=float)
    # gaia['magnitude_unc'] = at.Column(gaia['magnitude_unc'], 'magnitude_unc', dtype=float)
    # db.query("pragma foreign_keys=OFF")
    # db.add_data(gaia, table='photometry', clean_up=False)
    # db.query("pragma foreign_keys=ON")
    
    # Add the parallaxes
    gaia.rename_column('Plx', 'parallax')
    gaia.rename_column('e_Plx', 'parallax_unc')
    db.query("pragma foreign_keys=OFF")
    db.add_data(gaia, table='parallaxes', clean_up=False)
    db.query("pragma foreign_keys=ON")
            
    db.save()
    
def add_acs_data(db, cat):#, file=path+'raw_data/viz_acs_with_IDs.tsv'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    # Read in the data
    # acs = ascii.read(file)
    acs = at.Table.from_pandas(cat)
    
    # Rename some columns
    acs.rename_column('Obs', 'epoch')
    
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

def add_nicmos_data(db, cat):#file=path+'raw_data/viz_nicmos_with_IDs.tsv'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    # Read in the data
    # nic = ascii.read(file)
    nic = at.Table.from_pandas(cat)
    
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
                if not str(row['magnitude']).strip():
                    row['magnitude'] = np.nan
                if not str(row['magnitude_unc']).strip():
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

def add_wpc2_data(db, cat):#file=path+'raw_data/viz_wfpc2_with_IDs.tsv'):
    """
    Read in the Robberto+2013 ACS data and match objects by RA and Dec
    """
    # Read in the data
    # wpc = ascii.read(file)
    wpc = at.Table.from_pandas(cat)
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    wpc['publication_shortname'] = ['Robb13']*len(wpc)
    wpc['telescope_id'] = [1]*len(wpc)
    wpc['instrument_id'] = [3]*len(wpc)
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
            wpc.add_column(at.Column(['WFPC2.'+b.upper()]*len(wpc), 'band'))
            
            # Convert flag integer to string
            wpc['flags'] = at.Column([flags[i] for i in wpc['flags']], 'flags')
            
            # Move the magnitudes into the correct column
            for row in wpc:
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

def add_SDSS_data(db, cat):
    """
    Read in the cross matched SDSS data
    """
    # Read in the data
    sdss = at.Table.from_pandas(cat)
    
    # Rename some columns
    sdss.rename_column('ObsDate', 'epoch')
    sdss.rename_column('flags', 'oflags')
    sdss.rename_column('Q', 'flags')
    sdss.rename_column('e_umag','umag_unc')
    sdss.rename_column('e_gmag','gmag_unc')
    sdss.rename_column('e_rmag','rmag_unc')
    sdss.rename_column('e_imag','imag_unc')
    sdss.rename_column('e_zmag','zmag_unc')
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    sdss['publication_shortname'] = ['Alam15']*len(sdss)
    sdss['telescope_id'] = [6]*len(sdss)
    sdss['instrument_id'] = [9]*len(sdss)
    sdss['system_id'] = [2]*len(sdss)
    
    # Add the data
    db.add_data(sdss, table='photometry', bands=['umag','gmag','rmag','imag','zmag'], clean_up=False)
    
    # Change band names
    for i,j in zip(['umag','gmag','rmag','imag','zmag'],['SDSS.u','SDSS.g','SDSS.r','SDSS.i','SDSS.z']):
        db.modify("UPDATE photometry SET band='{}' WHERE band='{}'".format(j,i))
            
    db.save()

def add_2MASS_data(db, cat):
    """
    Read in the cross matched 2MASS data
    """
    # Read in the data
    tmass = at.Table.from_pandas(cat)
    
    # Rename some columns
    tmass.rename_column('JD', 'epoch')
    tmass.rename_column('Qflg', 'flags')
    tmass.rename_column('e_Jmag','Jmag_unc')
    tmass.rename_column('e_Hmag','Hmag_unc')
    tmass.rename_column('e_Kmag','Kmag_unc')
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    tmass['publication_shortname'] = ['Cutr03']*len(tmass)
    tmass['telescope_id'] = [2]*len(tmass)
    tmass['instrument_id'] = [5]*len(tmass)
    tmass['system_id'] = [2]*len(tmass)
    
    # Add the data
    db.add_data(tmass, table='photometry', bands=['Jmag','Hmag','Kmag'], clean_up=False)
    
    # Change band names
    for i,j in zip(['Jmag','Hmag','Kmag'],['2MASS.J','2MASS.H','2MASS.Ks']):
        db.modify("UPDATE photometry SET band='{}' WHERE band='{}'".format(j,i))
        
    db.save()

def add_WISE_data(db, cat):
    """
    Read in the cross matched WISE data
    """
    # Read in the data
    wise = at.Table.from_pandas(cat)
    
    # Rename some columns
    wise.rename_column('qph', 'flags')
    wise.rename_column('e_W1mag','W1mag_unc')
    wise.rename_column('e_W2mag','W2mag_unc')
    wise.rename_column('e_W3mag','W3mag_unc')
    wise.rename_column('e_W4mag','W4mag_unc')
    
    # Add columns for telescope_id, instrument_id, system_id, and publication_shortname
    wise['publication_shortname'] = ['Cutr13']*len(wise)
    wise['telescope_id'] = [3]*len(wise)
    wise['instrument_id'] = [6]*len(wise)
    wise['system_id'] = [2]*len(wise)
    
    # Add the data
    db.add_data(wise, table='photometry', bands=['W1mag','W2mag','W3mag','W4mag'], clean_up=False)
    
    # Change band names
    for i,j in zip(['W1mag','W2mag','W3mag','W4mag'],['WISE.W1','WISE.W2','WISE.W3','WISE.W4']):
        db.modify("UPDATE photometry SET band='{}' WHERE band='{}'".format(j,i))
        
    db.save()

def add_Hill13_data(db, cat):
    """
    Read in the cross matched Hillenbrand+13 data
    """
    # Read in the data
    hill13 = at.Table.from_pandas(cat)
    
    # Collect the spectral types
    spts = hill13['SpT2']
    
    # Convert the spectral types to integers
    spts = [specType(s) if isinstance(s, str) else [np.nan,''] for s in spts]
    typ, lc = np.array(spts).T
    
    # Add to the table
    hill13['spectral_type'] = typ
    hill13['luminosity_class'] = lc
    
    # Add column publication_shortname
    hill13['publication_shortname'] = ['Hill13']*len(hill13)
    
    # Add the data
    db.query("pragma foreign_keys=OFF")
    db.add_data(hill13, table='spectral_types', clean_up=False)
    db.query("pragma foreign_keys=ON")
    
    db.save()

def specType(SpT, types=[i for i in 'OBAFGKMLTY'], verbose=False):
    """
    Converts between float and letter/number spectral types (e.g. 14.5 => 'B4.5' and 'A3' => 23).
    
    Parameters
    ----------
    SpT: float, str
        Float spectral type or letter/number spectral type between O0.0 and Y9.9
    types: list
        The MK spectral type letters to include, e.g. ['M','L','T','Y']
      
    Returns
    -------
    list, str
        The converted spectral type string or (spectral type, luminosity class) numbers
    """
    result = [np.nan, '']
    try:
        # String input
        if isinstance(SpT, (str,bytes)):
            SpT = SpT.replace("'",'').replace('b','')
            val, LC = np.nan, ''
            if SpT[0] in types and SpT!='':
                MK, LC = SpT[0], 'V'
                suf = SpT[1:].replace('n','').replace('e','').replace('w','')\
                             .replace('m','').replace('a','').replace('Fe','')\
                             .replace('-1','').replace(':','').replace('?','')\
                             .replace('-V','').replace('p','').replace('<','')\
                             .replace('>','')
                
                if suf.replace('.','').isdigit():
                    val = float(suf)
                    
                else:
                
                    for cl in ['III','V','IV']:
                        try:
                            idx = suf.find(cl)
                            val = float(suf[:idx].split('/')[0])
                            LC = suf[idx:].split('/')[0].split(',')[0]
                            break
                        except:
                            try:
                                val = float(suf)
                            except:
                                continue
                                
                # return [types.index(MK)*10+val-(4. if MK in ['M','L','T','Y'] else 0), LC]
                return [types.index(MK)*10+val, LC]

        # Numerical input
        elif isinstance(SpT, float) or isinstance(SpT, int) and 0.0 <= SpT < len(types)*10:
            letter = ''.join(types)[int(SpT // 10)]
            number = int(SpT % 10) if SpT%10==int(SpT%10) else SpT%10
            result = '{}{}'.format(letter, number)

        # Bogus input
        else:
            if verbose:
                print('Sir, Spectral type',SpT,'must be a float between 0 and',len(types)*10,'or a string of class',types)
        
    except:
        pass
        
    return result