from utils import *


##########
######## Supplemental Figure
### related to experimental fig in main. Shows extra details of the linear stability fit.
## This figure shows: time evolutoin of the profiles, where in main figure only final frame is shown
# - The dots for individual fitted alphabar (compared to mean +- std dev in main figure)
################


######################### Load data

allbardata, fit_data = load_bar_data()

######## data for the example
files_examplebars = glob.glob(datafolder + "/experimental_data/init/linear_profile_data*")
bars_example = []
for f in files_examplebars:
    bar = f.split('.')[-2].split('bar')[-1]
    exp = f.split('_')[-2]
    bars_example.append((exp, bar))
    
profiles = {}

for exp, bar in bars_example:
    Xdata = np.loadtxt(datafolder + f'/experimental_data/init/linear_profile_data_{exp}_bar{bar}.txt')
    Xmodel = np.loadtxt(datafolder + f'/experimental_data/init/linear_profile_model_{exp}_bar{bar}.txt')

    profiles[(exp, bar)]= (Xdata, Xmodel)

cellpos_bar = []

### sort on length (long to short)
bars_example = [bars_example[i] for i in np.argsort([allbardata.loc[eb,'Length_micron'] for eb in bars_example])[-1::-1]]


###################
############### Make the figure

fig = plt.figure(figsize=(2*fig_width, 1.7*fig_width))

leftc, rightc, bottomc, topc  = 0.1, 0.95, 0.07, 0.95

gs_regimediag = gridspec.GridSpec(nrows = 1, ncols = 1, bottom=0.57, top=topc, left=leftc, right=.55)
gs_histograms = gridspec.GridSpec(nrows = 2, ncols = 1, bottom=gs_regimediag.bottom, top=topc, left=gs_regimediag.right+0.1, right=rightc, hspace=0.4)

lengths = np.array([allbardata.loc[(eb,'Length_micron')] for eb in bars_example])

gs_profiles = gridspec.GridSpec(nrows = 1, ncols = len(bars_example), width_ratios=lengths/ np.sum(lengths), left=leftc+0.04, right=rightc, bottom=bottomc, top=gs_regimediag.bottom-0.1)

ax_regimediag = fig.add_subplot(gs_regimediag[0])
ax_hist_alphabar = fig.add_subplot(gs_histograms[0])
ax_hist_alpha = fig.add_subplot(gs_histograms[1])

################ Regime diagram
#### Lines from theory, only show exact
alphav = np.linspace(0,20,201)
Lv = np.linspace(0,350/sigma0,201) # multiples of sigma in scaled units

aa,LL = np.meshgrid(alphav,Lv)
# draw them using parametrized curves
cmap = plt.colormaps['gray']
ncurves = 20
for n in range(1, ncurves):
    thetav = np.linspace(-np.pi, 0, 200)
    xv = -1j *(1+np.exp(1j*thetav))/(1-np.exp(1j*thetav))
    Lv = (thetav - n*np.pi)/xv
    av = 0.5*(xv**2+1)
    l, = ax_regimediag.plot(Lv*sigma0,av, label='$n={}$'.format(n), color = cmap((n+1)/ncurves), alpha=0.5)
 
ax_regimediag.set_ylabel(r"$\bar\alpha$", rotation=0)
ax_regimediag.set_xlabel(r'System size $L$ (µm)')

###### Add the fits
remove_40 = True
if remove_40:
    df_toplot = fit_data.query("Length_micron > 40")
else:
    df_toplot = fit_data

if 1: #### add jitter to better distinguish points
    L_offset = np.random.normal(size=fit_data.shape[0])*4
else:
    L_offset=0

ax_regimediag.scatter(df_toplot['Length_micron']+L_offset, df_toplot['alphabar'],\
                      c = get_color_length(df_toplot['Length_micron']), marker='.', s=10, alpha=0.7)

ax_regimediag.set_ylim((0, 20))
ax_regimediag.set_xlim((0, 10*sigma0))

# legend
ld = mlines.Line2D([],[], marker='.', linestyle='none', markersize=5, color='gray', alpha=0.3)
#lme = ax_regimediag.errorbar([-1],[0], yerr=[0], capsize=2, marker='.', markersize=5, linestyle='none', color='k') ### put outside of bounds (cheat)

ax_regimediag.legend([ld], [r'Fitted $\bar\alpha$'],\
                     bbox_to_anchor=(0.5, 0.97), loc='lower center', fontsize=8)

##### add markers for the example bars
for exp, bar in bars_example:
    alphabar_example, L_example = fit_data.loc[(exp, bar), ['alphabar','Length_micron']]
    mark_offset=2
    lx, = ax_regimediag.plot(L_example-mark_offset,alphabar_example, **data_marker, color='k', markersize=6)
    
################### Histograms

#### axes with histograms for the inferred parameters
rho_max = 1.2 #### to update based on the data -- doesn't actually matter that much (it's 1-1.4)

alphabarv = fit_data['alphabar']
rhov = fit_data['rho_avg']
alphav = alphabarv  / rhov/(1-rhov/rho_max)

N = len(alphabarv)

for dat, ax in zip([alphabarv, alphav], [ax_hist_alphabar, ax_hist_alpha]):
    n, binedges= np.histogram(dat, bins=20)
    binmids = 0.5*(binedges[1:]+binedges[:-1])
    ax.bar(x=binmids,height=n/N, width=binedges[1]-binedges[0], color='gray')
    
    ax.axvline(np.median(dat), color='tab:red')
    ax.set_xticks([0, np.median(dat), int(np.max(binedges))+1 ])
    ax.set_xticklabels([0, '{:.1f}'.format(np.median(dat)), int(np.max(binedges))+1 ])

    ax.set_ylim(0, 0.3)
    ax.set_yticks([0, 0.1, 0.2, 0.3])
ax_hist_alphabar.set_ylabel("Fraction")
ax_hist_alphabar.set_xlabel(r'$\bar\alpha$', labelpad=0)
ax_hist_alpha.set_xlabel(r'$\alpha$', labelpad=0)
ax_hist_alpha.set_ylabel("Fraction")


######################### Profiles

axes_profiles = []
for i in range(len(bars_example)):
    ax = fig.add_subplot(gs_profiles[i])
    axes_profiles.append(ax)
    
    exp, bar = bars_example[i]

    L = allbardata.loc[(exp, bar), 'Length_micron']

    
    Xdata, Xmodel = profiles[(exp, bar)] # includes both data and model
    ntimepoints = (Xdata.shape[0]-1)

    xdata = Xdata[0]
    xmodel = Xmodel[0]

    nsnaps = 10
    snaptimes = np.int16(np.linspace(0, ntimepoints-1, nsnaps))

    ### create a new axis for every profile
    x0,y0,w,h = ax.get_position().bounds
    ih = .01*h
    space_bottom=0.01
    hh = (h-(nsnaps-1)*ih-space_bottom)/nsnaps


    for j, st in enumerate(snaptimes):
        y_data = Xdata[st+1]
        y_model = Xmodel[st+1]
        current_ax = fig.add_axes((x0, y0+h-(j+1)*(hh+ih), w, hh))
        
        current_ax.plot(xdata, y_data-np.mean(y_data), linestyle = 'none', marker= '.', \
                        markersize=3,color='k', alpha=0.3+j/nsnaps*0.7)
        current_ax.plot(xmodel, y_model-np.mean(y_model), color='gray', alpha=0.1+j/nsnaps*0.9)

        current_ax.set_yticks([-0.3 ,0.3])
        current_ax.set_ylim(-0.4, 0.4)
        current_ax.set_xticks([])

        for s in ['top','right','bottom']:
            current_ax.spines[s].set_visible(0)
        if j < nsnaps-1 or i>0:
            current_ax.set_yticklabels([])
    
    for s in ['top','left','right']:
        ax.spines[s].set_visible(0)
    ax.set_xticks([0,L])
    ax.set_yticks([])

### add the x label directly to the figure
fig.text((gs_profiles.right+gs_profiles.left)/2, 0.01, "Position (µm)", ha='center', va='bottom')
axes_profiles[0].set_ylabel(r"Density deviation from mean $(\rho - \bar\rho)/\rho_\text{max}$", labelpad=10, y=0.55)
### time arrow
time_arrow = axes_profiles[0].annotate('', (-0.3, 1), (-0.3,0), xycoords='axes fraction', textcoords='axes fraction', arrowprops=dict(arrowstyle='<-'))
axes_profiles[0].text(time_arrow.xy[0]-0.1,0.5,'Time (min)', transform=axes_profiles[0].transAxes, rotation=90)
### ticks on the time arrow
axes_profiles[0].text(time_arrow.xy[0]-0.1,1.,'0', transform=axes_profiles[0].transAxes, va='top', ha='right')
axes_profiles[0].text(time_arrow.xy[0]-0.1,0.,'250', transform=axes_profiles[0].transAxes, va='bottom', ha='right')


# panel letters
letter_positions = {'a': (0.03, 0.97), 'b': (0.58, 'a'), 'c': ('a', 0.5)}
for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')


fig.savefig(savefolder + '/supplemental_figures/sfig_experimental_linearfit.pdf')
fig.savefig(savefolder + '/supplemental_figures/sfig_experimental_linearfit.svg')
