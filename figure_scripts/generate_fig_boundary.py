from utils import *


######## Figure with boundary effects. Mainly bifurcation diagrams. 
## Sketches added later
### Load data

df_bd_beta = pd.read_csv(datafolder + 'bifurcation_data/bd_data_beta.csv')


#### stability and pattern classification
stab_thr = 1e-4
df_bd_beta['stable'] = df_bd_beta['l0']<stab_thr
df_bd_beta['pattern_group'] = df_bd_beta['mode_k'].apply(classify_pattern_from_mode)


#### Reformat bifurcation data
df_bd_beta = reformat_branches(df_bd_beta)

### 2 par BD
df_2par = pd.read_csv(datafolder + 'bifurcation_data/2par_bifurcationdiagram_alphabeta.csv')

#### profiles to show
df_pointstoshow = pd.read_csv(datafolder + 'bifurcation_data/data_bifurcation_pointstoshow.csv')
df_pointstoshow['pattern_group'] = df_pointstoshow['mode_k'].apply(classify_pattern_from_mode)

#### drop the last two --- tmp
df_pointstoshow.drop(df_pointstoshow.index[-2:], inplace=True)



### Setup figure
fig = plt.figure(figsize=(fig_width, 1*fig_width))


topc, bottomc, leftc, rightc = 0.9, 0.12, 0.18, 0.94

gs_bd = gridspec.GridSpec(nrows = 2, ncols=1, hspace=0.35, wspace=0.3, \
        right=0.6, left=leftc, top=topc, bottom=bottomc)
gs_profiles = gridspec.GridSpec(nrows = 3, ncols=1, hspace=0.3, wspace=0.3, \
        right=rightc, left=gs_bd.right+0.1, top=topc, bottom=bottomc)

#### setup axes
## two par bifurcation diagram
ax_diag_2p = fig.add_subplot(gs_bd[0])

### bifurcation diag beta
ax_diag_beta = fig.add_subplot(gs_bd[1], sharex=ax_diag_2p)

#### axes for profiles
axes_profiles = [fig.add_subplot(gs_profiles[i]) for i in range(3)]

#### The two parameter bifurcation diagram beta alpha

### fill
polarized_region_color = plt.cm.Grays(0.3)
ax_diag_2p.fill_between(df_2par['beta'], df_2par['alpha'], 30*np.ones_like(df_2par['alpha']), color=polarized_region_color)
ax_diag_2p.fill_between([15,30], [15,30], [30,30] , color=polarized_region_color)

l,=ax_diag_2p.plot(df_2par['beta'], df_2par['alpha'], 'k')

## add asymptotes at beta=0 and beta=alpha
#### dashed extending to origin
ax_diag_2p.plot([0,0], [0,30], color=l.get_color(), linestyle='--')
ax_diag_2p.plot([0,30], [0,30], color=l.get_color(), linestyle='--')
###solid higher up
ax_diag_2p.plot([0,0], [15,30], color=l.get_color(), linestyle='-')
ax_diag_2p.plot([15,30], [15,30], color=l.get_color(), linestyle='-')

# add lines to indicate where the other diagrams are located
ax_diag_2p.axhline(20, color='tab:red', linestyle='-') 

# add some text
ax_diag_2p.text(0.55, 0.8, 'Polarization', transform=ax_diag_2p.transAxes, ha='center', bbox=dict(facecolor=polarized_region_color, edgecolor='none', pad=0))

ax_diag_2p.set_xlim(-10, 30)
ax_diag_2p.set_ylim(0, 30)
ax_diag_2p.set_xlabel(r"Boundary interaction $\beta$", labelpad=1)

ax_diag_2p.set_ylabel(r"Interparticle interaction $\alpha$", labelpad=5)

ax_diag_2p.set_aspect(1)


######## Diagram beta

#### boundaries of region of polarization... approximate, just for colors
beta1 = 0
beta2 = 19



# draw the bifurcation diagram
for bi, df in df_bd_beta.groupby('new_branch_index'): # iterate over the branches
    pv = df['par']
    muv = df['mu']
    stab = df['stable'].iloc[0]
    pg = df['pattern_group'].iloc[0]
    modek = df['mode_k'].iloc[0]

    ###### a bit ad hoc to get the colors right
    #### the 1.5 branch has label 5,6
    if bi[0]=='5' or bi[0]=='6':
        current_color = color_highlight
    elif bi[0]=='0' and stab:
        if np.max(pv)<beta1:
            current_color = get_color_pattern(2)
        if np.min(pv)>beta2:
            current_color = get_color_pattern(3)

    elif (bi[0] in '23') or (bi[0] in '14'):
        current_color = get_color_pattern(1)

    elif stab:
        current_color = color_stable
    else:
        current_color = color_unstable

    if stab:
        ax_diag_beta.plot(pv,muv,color=current_color, ls=linestyle_stable)
    else:
        ax_diag_beta.plot(pv,muv,color=current_color, ls=linestyle_unstable)

ax_diag_beta.set_xlabel("Boundary " + r"interaction $\beta$", labelpad=1)
ax_diag_beta.set_ylabel(r"Order parameter $\mu(\rho)$", labelpad=-1)

##### Profiles

for i, (_, r) in enumerate(df_pointstoshow.iterrows()):
    ax = axes_profiles[i]
    xu = np.loadtxt(datafolder + 'bifurcation_data/bifurcationdiagram_point_{}_branch{}_point{}.txt'.format(r['diagram'], r['branch_index'], r['point_index']))

    stab = r['l0']<stab_thr
    if stab:
        ax.plot(xu[0], xu[1], color=get_color_pattern(r['pattern_group']), linestyle=linestyle_stable)
    else:
        ax.plot(xu[0], xu[1], color='gray', linestyle=linestyle_unstable)

    ##### Add markers
    msize=7
    ax.plot(-0.1, .5, **data_marker, color=get_color_pattern(r['pattern_group']), markersize=msize, clip_on=False)
    if r['show_on_beta']:
        ax_diag_beta.plot(r['beta'], r['mu'], **data_marker, color=get_color_pattern(r['pattern_group']), markersize=msize)
    if r['show_on_2']:
        ax_diag_2p.plot(r['beta'], r['alpha'], **data_marker, color=get_color_pattern(r['pattern_group']), markersize=msize)

### labels etc
for i, ax in enumerate(axes_profiles):
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xlim(0,1)
    ax.set_ylim(-0.04,1.04)
    ax.spines.top.set_visible(False)
    ax.spines.right.set_visible(False)
    ax.spines.bottom.set_position(('outward',1))
    ax.spines.left.set_position(('outward',1))
    ax.spines.left.set_bounds((0,1))
    
    if i<2:
        ax.set_xticklabels([])
    else:
        ax.set_xticklabels([0, '$\ell$'])
        ax.set_xlabel("Position", labelpad=0)
    if i==0:
        ax.text(0., 1.05, r'Density $\rho/\rho_0$', transform=  ax.transAxes)

### align the y labels
for ax in [ax_diag_2p, ax_diag_beta]:
    ax.yaxis.set_label_coords(-0.25, 0.5)


# panel letters
letter_positions = {'a': (0.05, 0.95), 'b': (0.65, 'a'), 'c': ('a', 0.47)}

for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')


fig.savefig(savefolder + 'main_figures/fig_boundary_withoutsketches.svg')
fig.savefig(savefolder + 'main_figures/fig_boundary_withoutsketches.pdf')