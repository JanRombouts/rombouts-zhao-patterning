from utils import *

####

df_bd_L = pd.read_csv(datafolder + 'bifurcation_data/bd_data_L.csv')
stab_thr = 1e-4
df_bd_L['stable'] = df_bd_L['l0']<stab_thr
df_bd_L['pattern_group'] = df_bd_L['mode_k'].apply(classify_pattern_from_mode)

df_bd_L = reformat_branches(df_bd_L)

#### profiles to show
df_pointstoshow = pd.read_csv(datafolder + 'bifurcation_data/data_bifurcation_L_pointstoshow_supplemental.csv')
df_pointstoshow['pattern_group'] = df_pointstoshow['mode_k'].apply(classify_pattern_from_mode)


df_pointstoshow = df_pointstoshow.sort_values("L").reset_index()

#######

fig = plt.figure(figsize=(1.5*fig_width, 1.5*fig_width))

### in total 9 profiles

#### gs is just for diagram, we add the profile axes manually (because sizes are varying)
gs = gridspec.GridSpec(nrows = 1, ncols = 1, left=0.11, top=0.96, right=0.96, bottom=0.5)


ax_diag_L = fig.add_subplot(gs[0,:])
ax_diag_L.set_xlim(0, 7)
for bi, df in df_bd_L.groupby('new_branch_index'): # iterate over the branches
    pv = df['par']
    muv = df['mu']

    ##### after reformatting the branches, stability and patterngroup should be the same for the whole branch
    stab = df['stable'].iloc[0]
    pg = df['pattern_group'].iloc[0]
    modek = df['mode_k'].iloc[0]
    
    if stab:
        current_color = color_stable
    else:
        current_color = color_unstable

    if stab:
        ax_diag_L.plot(pv,muv,color=current_color, ls=linestyle_stable)
    else:
        ax_diag_L.plot(pv,muv,color=current_color, ls=linestyle_unstable)

ax_diag_L.set_xlabel(r'System size $\ell$', labelpad=0)
ax_diag_L.set_ylabel(r'Order parameter $\mu(\rho)$', labelpad=0)

######### profiles

### positioning is a bit fiddly
Lscale = 18
h = 0.1
hi=0.03
yv_rows = [0.32-i*(h+hi) for i in range(3)]
xv_col = [gs.left, 0.3, gs.right - df_pointstoshow['L'].max()/Lscale]
axes_profiles = []
for i, r in df_pointstoshow.iterrows():


    # adapt the length to L
    L = r['L']
    ax = fig.add_axes((xv_col[i//3], yv_rows[i%3], L/Lscale, h))

    x0, y0, w, h = ax.get_position().bounds
    #ax.set_position((x0, y0, w*L/Lref, h))

    axes_profiles.append(ax)

    xu = np.loadtxt(datafolder + 'bifurcation_data/bifurcationdiagram_point_{}_branch{}_point{}.txt'.format(r['diagram'], r['branch_index'], r['point_index']))

    stab = r['l0']<stab_thr
    if stab:
        ax.plot(xu[0], xu[1], color=color_stable, linestyle=linestyle_stable)
    else:
        ax.plot(xu[0], xu[1], color=color_unstable, linestyle=linestyle_unstable)

    ##### Add roman numeral
    ax_diag_L.text(r['L'], r['mu'], roman_numeral(i+1), ha='center', va='center', bbox=dict(facecolor='w', edgecolor='none', pad=1), fontstyle='italic')
    ax.text(0.03, 0.95, roman_numeral(i+1), ha='left', va='top', bbox=dict(facecolor='w', edgecolor='none', pad=1), fontstyle='italic', transform=ax.transAxes)
    #msize=10
    #ax_diag_alpha.plot(r['alpha'], r['mu'], **data_marker, color=get_color_pattern(r['pattern_group']), markersize=msize)
   
### labels etc
for i, ax in enumerate(axes_profiles):
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    if i//3!=0:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(r"Density $\rho$", labelpad=1)

    if i%3==2:
        ax.set_xticklabels([0, '$\ell$'])
        ax.set_xlabel("Position", labelpad=0)

    else:
        ax.set_xticklabels([])
   



fig.savefig(savefolder + '/supplemental_figures/sfig_bifurcationL.svg')
fig.savefig(savefolder + '/supplemental_figures/sfig_bifurcationL.pdf')

