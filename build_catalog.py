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
from sklearn.cluster import DBSCAN
from collections import Counter
from scipy.stats import norm
from sklearn.externals import joblib
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude

class Dataset(object):
    
    def __init__(self, name='Test'):
        """
        Initialize a database object
        
        Parameters
        ----------
        name: str
            The name of the database
        """
        self.name = name
        self.catalog = pd.DataFrame(columns=('source_id','ra','dec','flag','cat_name','catID'))
        self.n_sources = len(self.catalog)
        self.history = "{}: Database created".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.catalogs = {}
        self.xmatch_radius = 0
        
    @property
    def info(self):
        """
        Print the history
        """
        print(self.history)
        
    def ingest_VizieR(self, path, cat_name, id_col, count=-1):
        """
        Ingest a data file from Vizier and regroup sources
        
        Parameters
        ----------
        path: str
            The path to the exported VizieR data
        cat_name: str
            The name of the added catalog
        id_col: str
            The name of the column containing the unique ids
        count: int
            The number of table rows to add
            (This is mainly for testing purposes)
        """
        # Check if the catalog is already ingested
        if cat_name in self.catalogs:
            
            print('Catalog {} already ingested.'.format(cat_name))
            
        else:
            
            # Read in the TSV file and rename some columns
            data = pd.read_csv(path, sep='\t', comment='#', engine='python')[:count]
            # data = data.groupby(id_col).agg(lambda x: np.mean(x))
            data.insert(0,'catID', ['{}_{}'.format(cat_name,n) for n in range(len(data))])
            data.insert(0,'dec_corr', data['_DEJ2000'])
            data.insert(0,'ra_corr', data['_RAJ2000'])
            data.insert(0,'source_id', np.nan)
            data = data.reset_index(drop=True)
            
            print('Ingesting {} rows from {} catalog...'.format(len(data),cat_name))
            
            # Save the raw data as an attribute
            setattr(self, cat_name, data)
                
            # Update the history
            self.history += "\n{}: Catalog {} ingested.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),cat_name)
            self.catalogs.update({cat_name:(path,id_col)})
            
    def group_sources(self, radius=0.001, plot=False):
        """
        Calculate the centers of the point clusters given the
        radius and minimum number of points

        Parameters
        ----------
        coords: array-like
            The list of (x,y) coordinates of all clicks
        radius: int
            The distance threshold in degrees for cluster membership

        Returns
        -------
        np.ndarray
            An array of the cluster centers
        """
        # Gather the catalogs
        cats = pd.concat([getattr(self, cat_name) for cat_name in self.catalogs])
        
        # Clear the source grouping
        cats['oncID'] = np.nan
        cats['oncflag'] = ''
        self.xmatch_radius = radius
        
        # Make a list of the coordinates of each catalog row
        coords = cats[['ra_corr','dec_corr']].values
        
        # Perform DBSCAN to find clusters
        db = DBSCAN(eps=radius, min_samples=1).fit(coords)
        
        # Group the sources
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        source_ids = db.labels_
        unique_source_ids = list(set(source_ids))
        n_sources = len(unique_source_ids)
        
        # Get the average coordinates of all clusters
        unique_coords = np.asarray([np.mean(coords[source_ids==id], axis=0) for id in list(set(source_ids))])
        
        # Generate a source catalog
        self.catalog = pd.DataFrame(columns=('source_id','ra','dec','flag'))
        self.catalog['id'] = unique_source_ids
        self.catalog[['ra','dec']] = unique_coords
        self.catalog['flag'] = ['d{}'.format(i) for i in Counter(source_ids).values()]
        
        # Update history
        self.history += "\n{}: Catalog grouped with radius {} arcsec.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.xmatch_radius)
        
        # Update the source_ids in each catalog
        cats['source_id'] = source_ids
        for cat_name in self.catalogs:
            
            # Get the source_ids for the catalog
            cat_source_ids = cats.loc[cats['catID'].str.startswith(cat_name)]['source_id']
            
            # Get the catalog
            cat = getattr(self, cat_name)
            
            # Update the source_ids and put it back
            cat['source_id'] = cat_source_ids
            setattr(self, cat_name, cat)
            
            del cat, cat_source_ids
            
        del cats
        
        # Plot it
        if plot:
            plt.figure()
            plt.title('{} clusters for {} sources'.format(n_sources,len(coords)))
            
            colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, n_sources)]
            for k, col in zip(unique_source_ids, colors):
                
                class_member_mask = (source_ids == k)
                xy = coords[class_member_mask & core_samples_mask]
                
                marker = 'o'
                if len(xy)==1:
                    col = [0,0,0,1]
                    marker = '+'
                    
                plt.plot(xy[:, 0], xy[:, 1], color=tuple(col), marker=marker, markerfacecolor=tuple(col))

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
        
    def correct_offsets(self, cat_name, truth='ACS'):
        """
        Function to determine systematic, linear offsets between catalogs
        
        FUTURE -- do this with TweakReg, which also accounts for rotation/scaling
        See thread at https://github.com/spacetelescope/drizzlepac/issues/77
        
        Parameters
        ----------
        cat_name: str
            Name of catalog to correct
        truth: str
            The catalog to measure against
        """
        # Must be grouped!
        if not self.xmatch_radius:
            
            print("Please run group_sources() before running correct_offsets().")
            
        else:
            
            # First, remove any previous catalog correction
            self.catalog.loc[self.catalog['cat_name']==cat_name, 'ra_corr'] = self.catalog.loc[self.catalog['cat_name']==cat_name, '_RAJ2000']
            self.catalog.loc[self.catalog['cat_name']==cat_name, 'dec_corr'] = self.catalog.loc[self.catalog['cat_name']==cat_name, '_DEJ2000']
            
            # Copy the catalog
            onc_gr = self.catalog.copy()
            
            # restrict to one-to-one matches, sort by oncID so that matches are paired
            o2o_new = onc_gr.loc[(onc_gr['oncflag'].str.contains('o')) & (onc_gr['cat_name'] == cat_name) ,:].sort_values('oncID')
            o2o_old = onc_gr.loc[(onc_gr['oncID'].isin(o2o_new['oncID']) & (onc_gr['cat_name'] == truth)), :].sort_values('oncID')
            
            # get coords
            c_o2o_new = SkyCoord(o2o_new.loc[o2o_new['cat_name'] == cat_name, 'ra_corr'],\
                                 o2o_new.loc[o2o_new['cat_name'] == cat_name, 'dec_corr'], unit='degree')
            c_o2o_old = SkyCoord(o2o_old.loc[o2o_old['cat_name'] == truth, 'ra_corr'],\
                                 o2o_old.loc[o2o_old['cat_name'] == truth, 'dec_corr'], unit='degree')
                             
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
                
                # Fix precision
                mu_ra = round(mu_ra, 6)
                mu_dec = round(mu_dec, 6)
                
                # Update the coordinates of the appropriate sources
                print('Shifting {} sources by {}" in RA and {}" in Dec...'.format(cat_name,mu_ra,mu_dec))
                self.catalog.loc[self.catalog['cat_name']==cat_name, 'ra_corr'] += mu_ra
                self.catalog.loc[self.catalog['cat_name']==cat_name, 'dec_corr'] += mu_dec
                
                # Update history
                now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.history += "\n{}: {} sources shifted by {} deg in RA and {} deg in Declination.".format(now, cat_name, mu_ra, mu_dec)
                
                # Regroup the sources since many have moved
                self.group_sources(self.xmatch_radius)
                
            else:
                
                print('Cannot correct offsets in {} sources.'.format(cat_name))


class Dataset_old(object):
    
    def __init__(self, name='Test'):
        """
        Initialize a database object
        
        Parameters
        ----------
        name: str
            The name of the database
        """
        self.name = name
        self.catalog = pd.DataFrame(columns=('oncID','oncflag','cat_name','catID','ra_corr','dec_corr'))
        self.n_sources = len(self.catalog)
        self.history = "{}: Database created".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        self.catalogs = {}
        self.xmatch_radius = 0
        
    @property
    def info(self):
        """
        Print the history
        """
        print(self.history)
        
    def ingest_VizieR(self, path, cat_name, id_col, count=-1):
        """
        Ingest a data file from Vizier and regroup sources
        
        Parameters
        ----------
        path: str
            The path to the exported VizieR data
        cat_name: str
            The name of the added catalog
        id_col: str
            The name of the column containing the unique ids
        count: int
            The number of table rows to add
            (This is mainly for testing purposes)
        """
        # Check if the catalog is already ingested
        if cat_name in self.catalog['cat_name'].tolist():
            
            print('Catalog {} already ingested.'.format(cat_name))
            
        else:
            
            # Read in the TSV file and rename some columns
            raw_data = pd.read_csv(path, sep='\t', comment='#', engine='python')[:count]
            data = raw_data[[id_col,'_RAJ2000','_DEJ2000']].groupby(id_col).agg(lambda x: np.mean(x))
            data.insert(0,'dec_corr', data['_DEJ2000'])
            data.insert(0,'ra_corr', data['_RAJ2000'])
            data.insert(0,'catID', data.index)
            data.insert(0,'cat_name', cat_name)
            data.insert(0,'oncflag', '')
            data.insert(0,'oncID', np.nan)
            data = data.reset_index(drop=True)
            
            print('{} has {} unique sources in {} rows of data.'.format(cat_name,len(data),len(raw_data)))
            
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
                
                # # Symmetrize the matrix (not necessary!)
                # mtx = self.catalog[[str(i) for i in range(self.n_sources)]].as_matrix()
                # self.catalog[[str(i) for i in range(self.n_sources)]] = np.tril(mtx)+np.tril(mtx,-1).T
                
            # Update the history
            self.history += "\n{}: Catalog {} ingested.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),cat_name)
            
            self.catalogs.update({cat_name:(path,id_col)})

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
        self.xmatch_radius = dist_crit
        
        # Make a copy of the pairwise distances
        onc_df = self.catalog.copy()
        
        # 'new source' numbering starts at highest ACS number + 1
        exclude = set()
        
        print("Grouping sources with {} arcsec radius...".format(self.xmatch_radius))
        for k in range(len(onc_df)):
            
            if k not in exclude:
                
                # find ROWS where dist < dist_crit
                m = onc_df.loc[onc_df[str(k)] < self.xmatch_radius]
                
                mindex = set(m[str(k)].index.tolist())
                mindex_updated = mindex.copy()
                
                # initially set False to ensure it goes through the loop at least once
                mindex_same = False
                
                # keep adding match values until no new values are added
                while mindex_same == False:
                    for x in mindex:
                        y = onc_df.loc[onc_df[str(x)] < self.xmatch_radius]
                        
                        yindex = set(y[str(x)].index.tolist())
                        
                        mindex_updated.update(yindex)
                        
                    # drops out of loop if no new things are added
                    mindex_same = (mindex == mindex_updated)
                    
                    mindex.update(mindex_updated)
                    
                # if already grouped, don't need to do it again
                exclude.update(mindex)
                
                num_group = len(mindex)
                
                match = onc_df.loc[mindex,onc_df.columns[:len(onc_df.columns)-len(onc_df)]]
                
                # Get average RA and Dec for source if multiple matches
                if num_group>1:
                    
                    onc_df.loc[mindex,'ra_corr'] = np.mean(match['_RAJ2000'])
                    onc_df.loc[mindex,'dec_corr'] = np.mean(match['_DEJ2000'])
                    match = onc_df.loc[mindex,onc_df.columns[:len(onc_df.columns)-len(onc_df)]]
                    
                # check for multiple objects in same catalog (any duplicates will flag as True)
                if match.duplicated(subset='cat_name', keep=False).any() == True:
                    
                    onc_df.loc[mindex,'oncflag'] = 'd'
                    
                # check for one-to-one matches between ACS sources and new_cat sources (when new_cat is not ACS)
                elif 'ACS' in match['cat_name'].values:
                    
                    onc_df.loc[mindex,'oncflag'] = 'o'
                    
                onc_df.loc[mindex,'oncflag'] += str(num_group)
                
                # use ACS number if it exists -- if multiple, use lowest
                if ('ACS' in match['cat_name'].values) == True:
                    onc_df.loc[mindex,'oncID'] = min(match.loc[match['cat_name'] == 'ACS','catID'].values)
                    
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
        self.history += "\n{}: Catalog grouped with radius {} arcsec.".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), self.xmatch_radius)
                
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
        
    def export_IDs(self):
        """
        Put the new source IDs of the cross-matched catalog back into the raw data
        
        Parameters
        """
        for cat_name,(path,id_col) in self.catalogs.items():
            
            # Change the filename
            filename = path.replace('.tsv','_with_IDs.tsv')
            
            # Get object metadata
            cat = self.catalog.loc[:,['oncID','oncflag','cat_name','catID','ra_corr','dec_corr','_RAJ2000','_DEJ2000']]
            
            # Rename columns to fit the ONCdbWeb schema
            cat.rename(columns={'oncID':'source_id', 'oncflag':'comments', '_RAJ2000':'ra', '_DEJ2000':'dec'}, inplace=True)
            
            # Get the relevant rows
            cat = cat.loc[self.catalog['cat_name'] == cat_name, ['source_id','catID']]
            
            # Rename for easy merging
            cat.rename(columns={'catID':id_col}, inplace=True)
            
            # Get the raw data
            raw = pd.read_csv(path, sep='\t', comment='#', engine='python')
            
            # Merge the lists
            final = cat.merge(raw)
            
            # Write it to file
            final.to_csv(filename, sep='\t', index=False)
            
            print('{} IDs exported to {}.'.format(len(cat),filename))
        
    def correct_offsets(self, cat_name, truth='ACS'):
        """
        Function to determine systematic, linear offsets between catalogs
        
        FUTURE -- do this with TweakReg, which also accounts for rotation/scaling
        See thread at https://github.com/spacetelescope/drizzlepac/issues/77
        
        Parameters
        ----------
        cat_name: str
            Name of catalog to correct
        truth: str
            The catalog to measure against
        """
        # Must be grouped!
        if not self.xmatch_radius:
            
            print("Please run group_sources() before running correct_offsets().")
            
        else:
            
            # First, remove any previous catalog correction
            self.catalog.loc[self.catalog['cat_name']==cat_name, 'ra_corr'] = self.catalog.loc[self.catalog['cat_name']==cat_name, '_RAJ2000']
            self.catalog.loc[self.catalog['cat_name']==cat_name, 'dec_corr'] = self.catalog.loc[self.catalog['cat_name']==cat_name, '_DEJ2000']
            
            # Copy the catalog
            onc_gr = self.catalog.copy()
            
            # restrict to one-to-one matches, sort by oncID so that matches are paired
            o2o_new = onc_gr.loc[(onc_gr['oncflag'].str.contains('o')) & (onc_gr['cat_name'] == cat_name) ,:].sort_values('oncID')
            o2o_old = onc_gr.loc[(onc_gr['oncID'].isin(o2o_new['oncID']) & (onc_gr['cat_name'] == truth)), :].sort_values('oncID')
            
            # get coords
            c_o2o_new = SkyCoord(o2o_new.loc[o2o_new['cat_name'] == cat_name, 'ra_corr'],\
                                 o2o_new.loc[o2o_new['cat_name'] == cat_name, 'dec_corr'], unit='degree')
            c_o2o_old = SkyCoord(o2o_old.loc[o2o_old['cat_name'] == truth, 'ra_corr'],\
                                 o2o_old.loc[o2o_old['cat_name'] == truth, 'dec_corr'], unit='degree')
                             
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
                
                # Fix precision
                mu_ra = round(mu_ra, 6)
                mu_dec = round(mu_dec, 6)
                
                # Update the coordinates of the appropriate sources
                print('Shifting {} sources by {}" in RA and {}" in Dec...'.format(cat_name,mu_ra,mu_dec))
                self.catalog.loc[self.catalog['cat_name']==cat_name, 'ra_corr'] += mu_ra
                self.catalog.loc[self.catalog['cat_name']==cat_name, 'dec_corr'] += mu_dec
                
                # Update history
                now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                self.history += "\n{}: {} sources shifted by {} deg in RA and {} deg in Declination.".format(now, cat_name, mu_ra, mu_dec)
                
                # Regroup the sources since many have moved
                self.group_sources(self.xmatch_radius)
                
            else:
                
                print('Cannot correct offsets in {} sources.'.format(cat_name))

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
