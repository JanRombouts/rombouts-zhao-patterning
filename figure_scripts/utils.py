### utility functions for making the figures

import matplotlib.pyplot as plt
import h5py
import skimage.io
import dill
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec
import glob

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors

from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

######## some global parameters
fig_width = 8.6 / 2.54 # in inch. this is single col width, we can use double width or 1.5 width if necessary

## folder to save
savefolder = '../figures/'
datafolder = '../data_for_figures/'

# Set some other things
pxpermicron = 1.084

sigma0=35 ### space scale
tau0=sigma0**2 / 3 ##### time scale

for font_attribute in ['font.size','axes.labelsize', 'axes.titlesize', 'legend.fontsize', 'xtick.labelsize','ytick.labelsize']:
    plt.rcParams[font_attribute] = 8

plt.rcParams['font.sans-serif']='Arial'
plt.rcParams['lines.linewidth']=1
plt.rcParams['legend.frameon']=False

##### 
data_marker = dict(marker='*', markeredgecolor='none') #### for indicating things on diagrams

#### color maps
cmap_bra = plt.colormaps['plasma']
cmap_rho = mcolors.LinearSegmentedColormap.from_list('density colormap', ['#ffebd6', '#f29d21'])
cmap_h2b = plt.colormaps['Greys_r'] # mcolors.LinearSegmentedColormap.from_list('H2B color map', [[0,0,0], [1,0,0]])
cmap_eigenfunctions = plt.colormaps['tab20b']

### intensity values to use for color maps

df_intensityvalues = pd.read_csv(datafolder + 'experimental_data/intensityvalues_perexperiment_fromsegmentation.csv')
df_intensityvalues['Experiment'] = df_intensityvalues['Experiment'].astype(str)
df_intensityvalues.set_index(['Experiment', 'delta_frame'], inplace=True)

##### for the bifurcation diagrams
color_unstable='gray'
color_stable='k'
linestyle_stable='-'
linestyle_unstable=':'
color_highlight='k' ##### for the branch corresponding to k=1.5

### for plotting cell positions
color_bra_high_cells = 'black'
color_bra_low_cells = '#ebebeb'


def roman_numeral(i):
    if i<4:
        return i*'i'
    elif i==4:
        return 'iv'
    elif i<9:
        return 'v'+ (i-5)*'i'
    elif i==9:
        return 'ix'
    elif i==10:
        return 'x'
    else:
        return '/'

# function for the panel indications
def lf(l):
    # returns a string with the panel letter and a dictionary with options to pass to the text function 
    return '('+l.lower()+')', {'fontweight': 'bold'}

#### colormap and normalization to use for the Lengths
## don't use all colors per se because some are too faint
def get_color_length(L):
    maxL=350
    length_normalized = 0.25 + (L/maxL)/2
    return plt.cm.Blues_r(length_normalized)

#### Color scheme for the profile classification
def get_color_pattern(pattern_group):
    """ Pattern group is:
    0 for no pattern, 1 for polarized, 2 for middle peak, 3 for two peaks, 4 for three or more peaks

    use gray for the no pattern case
    """
    if pattern_group==0:
        return '#ebebeb'
    else:
        return plt.cm.viridis((pattern_group)/(max_pattern_n+1))

def classify_pattern_from_mode(k):
    """ k is wavenumber corresponding to dominant cosine mode (corrected for sign)"""
    if k==0:
        return 0 # no pattern
    if k==0.5:
        return 1 # polarized
    elif k==1: 
        return 2 # middle bump
    elif k==1.25 or k==1.5 or k==2:
        return 3 # two bumps
    else: # three or more
        return 4


def classify_pattern_from_peaks(r, peaktype='smooth'):
    ### r should be a series, with entry npeaks_plateau and entry xpeak_frac
    if r[f'npeaks_{peaktype}'] == 0:
        return 0
    elif r[f'npeaks_{peaktype}'] ==1:
        if r[f'xpeak_frac_{peaktype}'] < .4 or r[f'xpeak_frac_{peaktype}']>.6:
            return 1
        else:
            return 2
    elif r[f'npeaks_{peaktype}'] == 2:
        return 3
    else:
        return 4

max_pattern_n = 4 #### for using directly

#### load the pattern classification
### function returns two dataframes, one with the dominant mode and pscore per frame (delta_frame), 
### one 'reduced', where an overall patterning score and classification has been determined (so after selecting frames)
def load_pattern_classification(manual = True, pattern_threshold = 0.5, std_threshold=0.15 , density_threshold=0.15, delta_final=160, exclude_verylongbars = True):
    if manual:
        dfc_reduced = pd.read_csv(datafolder + 'experimental_data/experiment_pattern_classification_semimanual.csv')
        dfc_reduced['Experiment'] = dfc_reduced['Experiment'].astype(str)
        dfc_reduced['Bar'] = dfc_reduced['Bar'].astype(str)
        dfc_reduced = dfc_reduced.rename(columns={'pgroup': 'pattern_group'}).set_index(['Experiment', 'Bar'])
        return None, dfc_reduced
    else:
        dfc = pd.read_csv(datafolder + 'experimental_data/experiment_pattern_classification.csv')
        dfc['Experiment'] = dfc['Experiment'].astype(str)
        dfc['Bar'] = dfc['Bar'].astype(str)

        #### apply patterning criteria
        #dfc.loc[dfc.index[dfc['pscore']<pattern_threshold], ['mode_k', 'npeaks_smooth']] = 0
        ## set mode to zero for bars with low density amplitude
        #inds_nopattern = dfc.index[(dfc['std_rho']<std_threshold) | (dfc['avg_rho']<density_threshold)]
        inds_nopattern = dfc.index[(dfc['pscore']<1)]

        dfc.loc[inds_nopattern, ['mode_k', 'pgroup_sp']] = 0
        

        #### keep only the pattern group per bar, based on a selected number of frames
        def reduce(df):
            #### df will be the dataframe for one experiment, bar
            df_selectedframes = df.query("delta_frame>{} & delta_frame <= {}".format(delta_final-10, delta_final))
            pgroup_mode = pd.Series.mode(df_selectedframes['pgroup_sp'])[0]

            # keep L too and avg rho
            avg_rho = df_selectedframes['avg_rho'].mean()
            return pd.Series([df['L'].values[0], pgroup_mode, avg_rho], index=['L', 'pattern_group', 'avg_rho'])

        dfc_reduced = dfc.groupby(['Experiment','Bar']).apply(reduce)

        if exclude_verylongbars:
            dfc.drop(dfc.query("L>400").index, inplace=True)
            dfc_reduced.drop(dfc_reduced.query("L>400").index, inplace=True)

        return dfc, dfc_reduced



##### Loading general data
def load_bar_data(exclude_verylongbars = True):
    allbardata = pd.read_csv("../data/bardata_allexperiments_extended.csv")
    allbardata['Experiment'] = allbardata['Experiment'].astype(str)
    allbardata['Bar'] = allbardata['Bar'].astype(str)
    allbardata.set_index(['Experiment', 'Bar'], inplace=True)
    allbardata['Length_micron'] = np.round((allbardata['Length']/1.084/20))*20

    fit_data = pd.read_csv("../data/data_bar_fitalphabar.csv")
    fit_data['Experiment'] = fit_data['Experiment'].astype(str)
    fit_data['Bar'] = fit_data['Bar'].astype(str)
    fit_data.set_index(['Experiment', 'Bar'], inplace=True)
    #### add length to fit data
    fit_data = fit_data.join(allbardata['Length_micron'])

    if exclude_verylongbars:
        allbardata.drop(allbardata.query("Length_micron>400").index, inplace=True)
        fit_data.drop(fit_data.query("Length_micron>400").index, inplace=True)
    return allbardata, fit_data


############ For reformatting the branch dataset (separating stable/unstable)

#### Split branches containing both stable and unstable parts
### solution from chatgpt

def reformat_branches(df):
    """ Assuming that df is pandas dataframe with columns branch_index, stable and pattern_group."""
    df["group"] = (df["stable"] != df["stable"].shift()).cumsum()
    df["group2"] = (df["pattern_group"] != df["pattern_group"].shift()).cumsum()
    df["new_branch_index"] = df.groupby(["branch_index", "group", "group2"]).ngroup() + 1
    df["new_branch_index"] = df["branch_index"].astype(str) + "." + df["new_branch_index"].astype(str)
    df.drop(columns=['group', 'group2'], inplace=True)
    return df

##### order parameter
def get_mu(x, u):
    
    dx = np.diff(x)[0]
    L = x[-1]+dx/2
    xr = x/L #### because bifurcation diagrams are all on [0,1]
    dxr = np.diff(xr)[0] 

    #print(x, u, L, dx)
    L2norm = np.sum((u-np.mean(u))**2*dxr)**0.5
    prefactor = np.sum((u-np.mean(u))*((xr-1/2)*0.5 + (xr-1/2)**2)*dxr)
    return prefactor*L2norm
