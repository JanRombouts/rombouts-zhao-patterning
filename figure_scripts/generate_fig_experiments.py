
from utils import *

###################### Load experimental data we will show
allbardata, fit_data = load_bar_data()


######## Linear stability fits
files_init = glob.glob(datafolder + "/experimental_data/init/linear_profile_data*")
bars_example_init = []


#### keep only three (more in the supp figure)
L_to_use = [320., 200., 120.]

for f in files_init:
    bar = f.split('.')[-2].split('bar')[-1]
    exp = f.split('_')[-2]
    if allbardata.loc[(exp, bar), 'Length_micron'] in L_to_use:
        bars_example_init.append((exp, bar))
print(len(bars_example_init))
### rearrange on length
bars_example_init = [bars_example_init[i] for i in np.argsort([allbardata.loc[eb,'Length_micron'] for eb in bars_example_init])[-1::-1]]
bar_lengths_init = [allbardata.loc[ee, 'Length_micron'] for ee in bars_example_init]


############ Final bars to show
files_final = glob.glob(datafolder + "/experimental_data/final/densityprofiles*")
bars_example_final = []
for f in files_final:
    bar = f.split('.')[-2].split('bar')[-1]
    exp = f.split('_')[-2]
    bars_example_final.append((exp, bar))

#### Load the data from the parameter scan
df_paramscan = pd.read_csv(datafolder + "data_alphaLscan.csv")

dfpa_counts = df_paramscan.groupby(['L','pattern_group'])[['alpha']].count()

dfpa_counts['fraction'] = dfpa_counts['alpha'] / dfpa_counts.groupby('L')['alpha'].transform('sum')

#### add zeros
newidx = pd.MultiIndex.from_product([dfpa_counts.index.get_level_values(0).unique(), range(max_pattern_n+1)], names = dfpa_counts.index.names)
dfpa_counts = dfpa_counts.reindex(newidx, fill_value=0)


############## Final profile classification experiments
### Load data
dfc, dfc_r = load_pattern_classification()
# drop the 40 microns
dfc_r.drop(dfc_r[dfc_r['L']==40].index, inplace=True)
print(dfc_r.head())
## keep only those with avg_rho in range??
dfc_r.drop(dfc_r.index[ (dfc_r['avg_rho']<0.25) | (dfc_r['avg_rho']>0.55)], inplace=True)

#### get the counts per pattern type
df_exp_counts = dfc_r.groupby(['L','pattern_group'])[['avg_rho']].count()

#### add zeros if there are pattern groups that don't exist for given length
newidx = pd.MultiIndex.from_product([df_exp_counts.index.get_level_values(0).unique(), range(0, max_pattern_n+1)], names = df_exp_counts.index.names)
df_exp_counts = df_exp_counts.reindex(newidx, fill_value=0)
df_exp_counts['fraction'] = df_exp_counts['avg_rho'] / df_exp_counts.groupby('L')['avg_rho'].transform('sum')


################### Theoretical (simulated) examples of pattern groups
profiles_pgroups_th = {}
for pgroup in [1, 2, 3, 4]:
    xrho = np.loadtxt(f"../data_for_figures/simulations/density_evolution_pgroup{pgroup}.txt")
    profiles_pgroups_th[pgroup] = xrho
   

###############

######### Figure setup

fig = plt.figure(figsize=(2*fig_width, .9*fig_width))

leftc, rightc, topc, bottomc = 0.08, 0.97, 0.93, 0.1

gs_init = gridspec.GridSpec(nrows = 3, ncols = 1, left=leftc, right=0.3, top=topc, bottom=0.58, hspace=0.3)

gs_diagrams = gridspec.GridSpec(nrows = 1, ncols =3, width_ratios = [1, 1,0.6],  left=gs_init.right+0.07, right=rightc,\
 top=topc, bottom=gs_init.bottom, wspace=0.5)
gs_final = gridspec.GridSpec(nrows = 3, ncols = 4,  left=leftc, right=rightc, 
top=gs_diagrams.bottom-0.13, bottom=bottomc, hspace=0.1, height_ratios=[1, 1, 0.7])


##########################
#### Initial regime fit

#### new only the second profile
for i, (exp, bar) in enumerate(bars_example_init):
    L = allbardata.loc[(exp, bar), 'Length_micron']
    ### profiles
    Xdata = np.loadtxt(datafolder + f'experimental_data/init/linear_profile_data_{exp}_bar{bar}.txt')
    Xmodel = np.loadtxt(datafolder + f'experimental_data/init/linear_profile_model_{exp}_bar{bar}.txt')



    #ax_start = fig.add_subplot(gs_init[0,i])
    ax_end = fig.add_subplot(gs_init[i])

    #### start and end
    fstart, fend = fit_data.loc[(exp, bar), ['fstart', 'fend']]
    fstart = int(fstart)
    fend = int(fend)

    for j, profile_index, frame, ax in zip([1], [-1], [fend-1], [ax_end]):
       
        ax_rho = ax  ### remnant of old code where we split this axis
        
        #### cell density
        y_data = Xdata[profile_index]
        y_model = Xmodel[profile_index]

        ax_rho.plot(Xdata[0], y_data - np.mean(y_data), '.', color='k', markersize=3, label='Data')
        ax_rho.plot(Xmodel[0], y_model - np.mean(y_model), '-', color='gray', label='Model fit')

        ax_rho.set_ylim(-0.55, 0.55)
        ax_rho.spines['left'].set_bounds(-0.5, 0.5)
        ax_rho.spines['bottom'].set_position(('outward', 1.5))
        for s in ['top','right']:
            ax_rho.spines[s].set_visible(0)

        
        ax_rho.set_xlim(0, np.max(bar_lengths_init))
        ax_rho.set_xticks([0, 80, 160, 240, 320])
        if i==2:
            ax_rho.set_xticklabels(ax_rho.get_xticks())
            ax_rho.set_xlabel("Position (µm)", labelpad=0)
            ax_rho.set_ylabel(r"$\rho-\bar\rho$")
            ax_rho.legend(handlelength=1, labelspacing=0.2)
        else:
            ax_rho.set_xticklabels([])
            #ax_rho.set_yticklabels([])
        
        
##### add a little arrow in the middle
#arrow = plt.annotate('', (0.12, .95), (0.12,.8), xycoords='figure fraction', textcoords='figure fraction', arrowprops=dict(arrowstyle='<-'))
#fig.text(arrow.xy[0]-0.02,arrow.xy[1]-0.075,'250 min', va='center', ha='right')



##################
######### Diagrams, theoretical and experimental, with the fractions of observed patterns

#################### For getting the colors right etc.
values = [0,1,2,3,4]
cmap = mcolors.ListedColormap([get_color_pattern(v) for v in values])
bounds = [v-0.5 for v in values] + [values[-1]+0.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

c_offset = 1
npatterngroups = len(dfpa_counts.index.get_level_values(1).unique())
colorlist = cmap(norm(values))
pattern_labels = {0: "No pattern", 1: 'Polarized', 2: 'Middle peak', 3: 'Two peaks',4: "More peaks"}




ax_observedpatterns_th = fig.add_subplot(gs_diagrams[:,0])
ax_observedpatterns_exp = fig.add_subplot(gs_diagrams[:,1])

##### Simulations
Lv = dfpa_counts.index.get_level_values(0).unique()

bottoms = np.zeros(len(Lv))
dL = Lv[1]-Lv[0]
for i, (pat, df) in enumerate(dfpa_counts.groupby('pattern_group')):
    ax_observedpatterns_th.bar(df.index.get_level_values(0), df['fraction'], bottom=bottoms, width=dL, label=pattern_labels[pat], color=colorlist[i])
    bottoms+=df['fraction'].values

ax_observedpatterns_th.set_xlabel("System size $\ell$", labelpad=1)
ax_observedpatterns_th.set_xticks(np.arange(0, np.max(Lv), 1))
ax_observedpatterns_th.set_ylabel("Fraction of simulations", labelpad=0)

#### bg color: no pattern
ax_observedpatterns_th.set_facecolor(get_color_pattern(0))
## lim x
ax_observedpatterns_th.set_xlim(0,10)
ax_observedpatterns_th.set_xticks(np.arange(0,12,2))

##### add little arrows corresponding to the ell values in the saddles figure (CHECK if alpha is the same)
ell_values = [1.3, 2.3, 3.5] # manually determined

for ell in ell_values:
    ax_observedpatterns_th.axvline(ell, color='k')
    
ax_observedpatterns_th.set_xticks([0] + ell_values + [5, 10])
ax_observedpatterns_th.set_xticklabels([0] +  [r'$\ell_\mathrm{PF}$', r'$\ell_1^*$', r'$\ell_2^*$'] + [5, 10], fontsize=6)

##### Experiments

Lv = df_exp_counts.index.get_level_values(0).unique()
bottoms = np.zeros(len(Lv))
for i, (pat, df) in enumerate(df_exp_counts.groupby('pattern_group')):
    ax_observedpatterns_exp.bar(df.index.get_level_values(0), df['fraction'], bottom=bottoms, width=16, label=pattern_labels[pat], color=colorlist[i])
    bottoms+=df['fraction'].values


### add number of repeats for each length
reps_per_length = dfc_r.groupby("L").count().iloc[:,0].to_dict()
for L,n in reps_per_length.items():
    if L<45:
        continue
    else:
        ax_observedpatterns_exp.text(L, 1.03, f'{n}', clip_on=False, ha='center', fontsize=6)
# add a small 'n'
ax_observedpatterns_exp.text(80-20, 1.03, 'N:', clip_on=False, ha='center', fontsize=6)


ax_observedpatterns_exp.set_ylabel("Fraction of experiments", labelpad=0)
ax_observedpatterns_exp.set_xlabel("Bar length (µm)", labelpad=1)
ax_observedpatterns_exp.set_xticks(df_exp_counts.index.get_level_values(0).unique())
plt.setp(ax_observedpatterns_exp.get_xticklabels()[1], visible=False)
plt.setp(ax_observedpatterns_exp.get_xticklabels(), fontsize=6)

#### y ticks
for ax in [ax_observedpatterns_th, ax_observedpatterns_exp]:
    ax.set_yticks([0, 1/2, 1])
    
### small axes with the pattern groups (simulation examples)

### new: add a colored rectangle
ax_states=  fig.add_subplot(gs_diagrams[2])
x0, y0, w, h = ax_states.get_position().bounds
ax_states.set_axis_off()
gs_thpatterns = gridspec.GridSpec(nrows=5, ncols=1, left=x0+0.048, right=x0+w, bottom=y0, top=y0+h, hspace=0.3, height_ratios=[1, 1, 1, 1, 1])

for i in range(5):
    ax = fig.add_subplot(gs_thpatterns[i])
    
    tt= ax.text(-1.65, 0.5, pattern_labels[i], ha='left', transform=ax.transAxes, fontsize=8, va='center')
    wr =0.3
    hr=0.2
    rect =mpatches.Rectangle((tt._x-wr-0.05, 0.5-hr/2), wr, hr, fc=get_color_pattern(i), ec='none', clip_on=False, transform=ax.transAxes)
    ax.add_artist(rect)
    ax.set_axis_off()

    if i>0:
        xrho = profiles_pgroups_th[i]

        ax.plot(xrho[0], xrho[-1],color = get_color_pattern(i))
        ax.set_ylim(-0.1, 1.1)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_yticklabels([])
    #ax.set_title(pattern_labels[i+1], color=get_color_pattern(i+1), fontsize=6, y=1., pad=1)


################ Bar examples final
delta_frame_to_use=160
for i, (exp, bar) in enumerate(bars_example_final):
    ax_image = fig.add_subplot(gs_final[0,i])
    ax_cellpos = fig.add_subplot(gs_final[1,i])
    ax_rho = fig.add_subplot(gs_final[2,i])

    ##### Load data from file
    finalframe = delta_frame_to_use + int(allbardata.loc[(exp, bar), 'frame0'])
    L = allbardata.loc[(exp, bar), 'Length_micron']

    vmin_bra, vmax_bra = df_intensityvalues.loc[(exp, delta_frame_to_use), ['Bra_p0.5', 'Bra_p99.5']]

    ## Bra image
    im = skimage.io.imread(datafolder + f'experimental_data/final/{exp}_bar{bar}_brachyury_maxz2-4.tif')[finalframe]
    ## cell positions dataframe
    cellposdf = pd.read_csv(datafolder + f'experimental_data/final/cellpositions_{exp}_bar{bar}.csv').query("frame==@finalframe")
    ## density profile
    xrho = np.loadtxt(datafolder + f'experimental_data/final/densityprofiles_{exp}_bar{bar}.txt')

    df_high = cellposdf.query("type=='high' &  Class==2")
    df_low = cellposdf.query("type=='low' &  Class==2")


    x = np.arange(im.shape[1])/pxpermicron
    y = np.arange(im.shape[0])/pxpermicron

    for ax in [ax_image, ax_cellpos]:
        ax.set_xlim(0,L)
        ax.set_ylim(0,40)
        ax.set_xticks([0, L])
        ax.set_yticks([0, 40])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if i==0:
            ax.set_yticklabels([0, 40])
    if i==0:
        ax_cellpos.set_ylabel("Position\n(µm)", y=1.5)
    
    ## Cell positions high Bra cells
    x_cells_high = df_high['x']-df_high['min_col']
    y_cells_high = df_high['y']-df_high['min_row']
    # low brachyury cells
    x_cells_low = df_low['x']-df_low['min_col']
    y_cells_low = df_low['y']-df_low['min_row']
    
    if allbardata.loc[(exp, bar), 'Orientation'] == 'horizontal':
        ax_image.pcolormesh(x, y, im, cmap='plasma', shading='nearest', rasterized=True, vmax=vmax_bra, vmin=vmin_bra)
        ax_cellpos.plot(x_cells_low/pxpermicron, y_cells_low/pxpermicron, '.', color=color_bra_low_cells, markersize=3)
        ax_cellpos.plot(x_cells_high/pxpermicron, y_cells_high/pxpermicron, '.', color=color_bra_high_cells, markersize=3)
    else: # switch x and y
        ax_image.pcolormesh(y, x, im.T, cmap='plasma', shading='nearest', rasterized=True, vmax=vmax_bra, vmin=vmin_bra)
        ax_cellpos.plot(y_cells_low/pxpermicron, x_cells_low/pxpermicron, '.', color=color_bra_low_cells, markersize=3)
        ax_cellpos.plot(y_cells_high/pxpermicron, x_cells_high/pxpermicron, '.', color=color_bra_high_cells, markersize=3)

    ax_rho.plot(xrho[0], xrho[1+finalframe], color=get_color_pattern(i+1)) #### shortcut since patterns are ordered left-right
    ax_rho.set_ylim(-0.1, 1.3)
    ax_rho.spines['left'].set_bounds(0,1)
    ax_rho.set_yticks([0,1])
    ax_rho.set_xticks([0, L])
    if i==0:
        ax_rho.set_ylabel("Density\n(cells/µm)")
    else:
        ax_rho.set_yticklabels([])
    ax_rho.set_xlim([0,L])
    ax_rho.set_xlabel("Position (µm)", labelpad=-5)
    ax_rho.spines.top.set_visible(0)
    ax_rho.spines.right.set_visible(0)
    ax_rho.spines.bottom.set_position(('outward', 2))



    ax_image.set_aspect(1, anchor='W')
    ax_cellpos.set_aspect(1, anchor='W')

    ### make rho axis same width as the others
    new_w = ax_image.get_position().bounds[2]
    xr, yr, wr, hr = ax_rho.get_position().bounds
    ax_rho.set_position((xr, yr, new_w, hr))

# panel letters
letter_positions = {'a': (0.03, 0.96), 'b': (0.32, 'a'), 'c': ('a', 0.46)}

for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')


fig.savefig(savefolder + 'main_figures/fig_experiments.svg')
fig.savefig(savefolder + 'main_figures/fig_experiments.pdf')