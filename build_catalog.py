"""
Module to ingest arbitrary VizieR catalogs into a custom cross-matched database

Authors: Andrea Lin, Joe Filippazzo
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import astropy.units as u
import datetime
from scipy.stats import norm
from sklearn.externals import joblib
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude

class Database(object):
    
    def __init__(self, name='Test'):
        """
        Initialize a database object
        
        Parameters
        ----------
        name: str
            The name of the database
        """
        self.name = name
        self.catalog = pd.DataFrame(columns=('oncID','oncflag','catname','catID','ra_corr','dec_corr'))
        self.n_sources = len(self.catalog)
        self.history = "{}: Database created".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.grouped = False
        
    @property
    def info(self):
        """
        Print the history
        """
        print(self.history)
        
    def ingest_VizieR(self, viz_output, cat_name, id_col, radius=0.25):
        """
        Ingest a data file from Vizier and regroup sources
        
        Parameters
        ----------
        viz_output: str
            The path to the exported VizieR data
        cat_name: str
            The name of the added catalog
        id_col: str
            The name of the column containing the unique ids
        radius: float
            The distance in arcsec to be used for cross matching
        """
        # Check if the catalog is already ingested
        if cat_name in self.catalog['catname'].tolist():
            
            print('Catalog {} already ingested.'.format(cat_name))
            
        else:
            
            # Read in the TSV file and rename some columns
            raw_data = pd.read_csv(viz_output, sep='\t', comment='#', engine='python')[:10]
            data = raw_data[[id_col,'_RAJ2000','_DEJ2000']].groupby(id_col).agg(lambda x: np.mean(x))
            data.insert(0,'dec_corr', data['_DEJ2000'])
            data.insert(0,'ra_corr', data['_RAJ2000'])
            data.insert(0,'catID', data.index)
            data.insert(0,'catname', cat_name)
            data.insert(0,'oncflag', '')
            data.insert(0,'oncID', np.nan)
            data = data.reset_index(drop=True)
            
            print('{} has {} sources in {} rows of data.'.format(cat_name,len(data),len(raw_data)))
            
            # Get sky coordinates of new sources from RA and Dec
            new_coords = SkyCoord(data['_RAJ2000'], data['_DEJ2000'], unit='degree')
            
            # If it is the initial table...
            if self.catalog.empty:
                
                # Generate a pandas data frame
                build_dist = pd.DataFrame()
                print("Measuring pairwise distances...")
                for i,c in enumerate(new_coords):
                    sep = new_coords.separation(c).arcsecond
                    build_dist.loc[:,i] = sep
                    progress_meter((i+1)*100./len(new_coords))
                    
                print('\n')
                
                # Add the distances to each row
                new_dist = pd.concat([data, build_dist], axis=1)
                new_dist.columns = new_dist.columns.astype(str)
                new_dist.index = new_dist.index.astype(str)
                
                # Update the attribute
                self.catalog = new_dist
                self.n_sources = len(self.catalog)
                
            # If not first ingested table...
            else:
                
                # Get sky coordinaes of existing sources
                coords = SkyCoord(self.catalog['_RAJ2000'], self.catalog['_DEJ2000'], unit='degree')
                
                # Calculate separations
                cross_dist = pd.DataFrame()
                self_dist = pd.DataFrame()
                print("Measuring pairwise distances...")
                for i,c in enumerate(new_coords):
                    
                    # sep between new_cat and existing oncdb objects
                    sep_cross = coords.separation(c).arcsecond
                    
                    # internal sep between new_cat objects
                    sep_self = new_coords.separation(c).arcsecond
                    
                    cross_dist.loc[:,i] = sep_cross
                    self_dist.loc[:,i] = sep_self
                    
                    progress_meter((i+1)*100./len(new_coords))
                    
                print('\n')
                
                # offsetting indices to make it join properly
                nc_join = data.rename(index = lambda x: (int(x) + self.n_sources), inplace=False)
                cross_dist.rename(columns = lambda x: (int(x) + self.n_sources), inplace=True)
                self_dist.rename(columns = lambda x: (int(x) + self.n_sources), inplace=True)
                self_dist.rename(index = lambda x: (int(x) + self.n_sources), inplace=True)
                
                # join
                pw1 = self.catalog.join(cross_dist)
                pw2 = nc_join.join(cross_dist.transpose().join(self_dist))
                
                pw1.columns = pw1.columns.astype(str)
                pw2.columns = pw2.columns.astype(str)
                
                # Update the catalog
                self.catalog = pd.concat([pw1,pw2], ignore_index=True)
                self.n_sources = len(self.catalog)
                
                # Symmetrize the matrix
                mtx = self.catalog[[str(i) for i in range(self.n_sources)]].as_matrix()
                self.catalog[[str(i) for i in range(self.n_sources)]] = np.tril(mtx)+np.tril(mtx,-1).T 
                
            # Update the history
            self.history += "\n{}: Catalog {} ingested.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),cat_name)

    def group_sources(self, dist_crit=0.25):
        """
        Function to group sources by some critical distance
        
        Parameters
        ----------
        dist_crit: float
            The cross-matching radius (observations within this radius will be grouped together as a source)
        """
        # Clear the source grouping
        self.catalog['oncID'] = np.nan
        self.catalog['oncflag'] = ''
        
        # Make a copy of the pairwise distances
        onc_df = self.catalog.copy()
        
        # 'new source' numbering starts at highest ACS number + 1
        new_source = max(onc_df.loc[onc_df['catname'] == 'ACS', 'catID'].values) + 1
        
        exclude = set()
        
        print("Grouping sources with {} arcsec radius...".format(dist_crit))
        for k in range(len(onc_df)):
            
            if k not in exclude:
                
                # find where dist < dist_crit
                m = onc_df.loc[onc_df[str(k)] < dist_crit]
                
                mindex = set(m[str(k)].index.tolist())
                mindex_updated = mindex.copy()
                
                # initially set False to ensure it goes through the loop at least once
                mindex_same = False
                
                # keep adding match values until no new values are added
                while mindex_same == False:
                    for x in mindex:
                        y = onc_df.loc[onc_df[str(x)] < dist_crit]
                        
                        yindex = set(y[str(x)].index.tolist())
                        
                        mindex_updated.update(yindex)
                        
                    # drops out of loop if no new things are added
                    mindex_same = (mindex == mindex_updated)
                    
                    mindex.update(mindex_updated)
                    
                # if already grouped, don't need to do it again
                exclude.update(mindex)
                
                num_group = len(mindex)
                
                match = onc_df.loc[mindex,['catname','catID']]
                
                # check for multiple objects in same catalog (any duplicates will flag as True)
                if match.duplicated(subset='catname', keep=False).any() == True:
                    
                    onc_df.loc[mindex,'oncflag'] = 'd'
                    
                # check for one-to-one matches between ACS sources and new_cat sources (when new_cat is not ACS)
                elif 'ACS' in match['catname'].values:
                    
                    onc_df.loc[mindex,'oncflag'] = 'o'
                    
                onc_df.loc[mindex,'oncflag'] += str(num_group)
                
                # use ACS number if it exists -- if multiple, use lowest
                if ('ACS' in match['catname'].values) == True:
                    onc_df.loc[mindex,'oncID'] = min(match.loc[match['catname'] == 'ACS','catID'].values)
                    
                # otherwise give it a new number
                else:
                    onc_df.loc[mindex,'oncID'] = new_source
                    new_source += 1
                    
                progress_meter((k+1)*100./len(onc_df))
                
        print('\n')
        
        # change id columns to ints (defaults to floats...)
        onc_df.loc[:,'catID'] = onc_df.loc[:,'catID'].astype(int)
        onc_df.loc[:,'oncID'] = onc_df.loc[:,'oncID'].astype(int)
        
        # Update the catalog
        self.catalog = onc_df
        
        # Update the history
        self.history += "\n{}: Catalog grouped with radius {} arcsec.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),dist_crit)
        
        self.grouped = True
        
    def load(self, path):
        """
        Load the catalog from file
        
        Parameters
        ----------
        path: str
            The path to the file
        """
        # Get the object
        DB = joblib.load(path)
        
        # Load the attributes
        self.catalog   = DB.catalog
        self.n_sources = DB.n_sources
        self.name      = DB.name
        self.history   = DB.history
        
        del DB
        
    def save(self, path):
        """
        Save the catalog to file for faster loading next time
        
        Parameters
        ----------
        path: str
            The path to the file
        """
        joblib.dump(self, path)
        
    def correct_offsets(self, catname, truth='ACS'):
        """
        Function to determine systematic, linear offsets between catalogs
        
        FUTURE -- do this with TweakReg, which also accounts for rotation/scaling
        See thread at https://github.com/spacetelescope/drizzlepac/issues/77
        
        Parameters
        ----------
        catname: str
            Name of catalog to correct
        truth: str
            The catalog to measure against
        """
        # Must be grouped!
        if not self.grouped:
            
            print("Please run group_sources() before running correct_offsets().")
            
        else:
            
            onc_gr = self.catalog.copy()
            
            # restrict to one-to-one matches, sort by oncID so that matches are paired
            o2o_new = onc_gr.loc[(onc_gr['oncflag'].str.contains('o')) & (onc_gr['catname'] == catname) ,:].sort_values('oncID')
            o2o_old = onc_gr.loc[(onc_gr['oncID'].isin(o2o_new['oncID']) & (onc_gr['catname'] == truth)), :].sort_values('oncID')
            
            # get coords
            c_o2o_new = SkyCoord(o2o_new.loc[o2o_new['catname'] == catname, 'ra_corr'],\
                                 o2o_new.loc[o2o_new['catname'] == catname, 'dec_corr'], unit='degree')
            c_o2o_old = SkyCoord(o2o_old.loc[o2o_old['catname'] == truth, 'ra_corr'],\
                                 o2o_old.loc[o2o_old['catname'] == truth, 'dec_corr'], unit='degree')
                             
            print(len(c_o2o_old), 'one-to-one matches found!')
            
            if len(c_o2o_old)>0:
                
                delta_ra = []
                delta_dec = []
                
                for i in range(len(c_o2o_old)):
                    # offsets FROM ACS TO new catalog
                    ri, di = c_o2o_old[i].spherical_offsets_to(c_o2o_new[i])
                    
                    delta_ra.append(ri.arcsecond)
                    delta_dec.append(di.arcsecond)
                    
                    progress_meter((i+1)*100./len(c_o2o_old))
                    
                delta_ra = np.array(delta_ra)
                delta_dec = np.array(delta_dec)
                
                print('\n')
                
                # fit a gaussian
                mu_ra, std_ra = norm.fit(delta_ra)
                mu_dec, std_dec = norm.fit(delta_dec)
                
                # Update the coordinates of the appropriate sources
                self.catalog.loc[self.catalog['catname']==catname]['ra_corr'] += mu_ra
                self.catalog.loc[self.catalog['catname']==catname]['dec_corr'] += mu_dec
                
                print('Delta RA (arcsec):', mu_ra)
                print('Delta DEC (arcsec):', mu_dec)
                
                # Update history
                now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.history += "\n{}: {} sources shifted by {} deg in RA and {} deg in Declination.".format(now, catname, mu_ra, mu_dec)
                
                # return (delta_ra, delta_dec, mu_ra, mu_dec, std_ra, std_dec)
                
            else:
                
                print('Cannot correct offsets in {} sources.'.format(catname))

def progress_meter(progress):
    """
    Print nice progress update
    
    Parameters
    ----------
    progress: float
        Some fraction of the completed job
    """
    sys.stdout.write("\rloading... %.1f%%" % progress)
    sys.stdout.flush()
