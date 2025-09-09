from utils import *


####
##### Supp figure, to go with the bifurcation diagram main figure
#### We show the one parameter diagram with alpha and add many different profiles
#### later addition: put the eigenvalue diagram here too


df_bd_alpha = pd.read_csv(datafolder + 'bifurcation_data/bd_data_alpha.csv')
stab_thr = 1e-4
df_bd_alpha['stable'] = df_bd_alpha['l0']<stab_thr
df_bd_alpha['pattern_group'] = df_bd_alpha['mode_k'].apply(classify_pattern_from_mode)

df_bd_alpha = reformat_branches(df_bd_alpha)

#### profiles to show
df_pointstoshow = pd.read_csv(datafolder + 'bifurcation_data/data_bifurcation_pointstoshow_supplemental.csv')
df_pointstoshow['pattern_group'] = df_pointstoshow['mode_k'].apply(classify_pattern_from_mode)


### manually give some new branch labels that will be colored
def new_branch_index(s):
    # s is series, this functino should be applied to rows of dataframe
    bi = s['branch_index']
    mu = s['mu']
    if bi == 0:
        if s['stable']:
            return '0a'
        else:
            return '0b'
    if bi == 1 or bi==2:
        if mu>0:
            return '1.a'
        else:
            return '1.b'
    if bi==3 or bi==4:
        if mu>0:
            return '3.a'
        else:
            return '3.b'
    if bi==5 or bi==6:
        if mu>0:
            return '5.a'
        else:
            return '5.b'
    if bi==7 or bi==8:
        if mu>0:
            return '7.a'
        else:
            return '7.b'
    if bi==9 or bi==10:
        if mu>0:
            return '9.a'
        else:
            return '9.b'

df_bd_alpha['new_branch_index2'] = df_bd_alpha.apply(new_branch_index, axis=1)

#######

fig = plt.figure(figsize=(2*fig_width, 1.2*fig_width))

### in total 9 profiles
gs_diag = gridspec.GridSpec(nrows = 2, ncols = 1, height_ratios= [3, 2], top=0.96, right=0.5)
gs_profiles = gridspec.GridSpec(nrows = 3, ncols = 3, left= gs_diag.right+0.1, top=gs_diag.top, right=0.97)

ax_diag_alpha = fig.add_subplot(gs_diag[0])
ax_diag_eig = fig.add_subplot(gs_diag[1], sharex=ax_diag_alpha)
for bi, df in df_bd_alpha.groupby('new_branch_index2'): # iterate over the branches
    df=df.sort_values('par')
    pv = df['par']
    muv = df['mu']

    ##### after reformatting the branches, stability and patterngroup should be the same for the whole branch
    stab = df['stable'].iloc[0]
    pg = df['pattern_group'].iloc[0]
    modek = df['mode_k'].iloc[0]

    if stab:
        l, = ax_diag_alpha.plot(pv,muv,ls=linestyle_stable, lw=1)
    else:
        l, = ax_diag_alpha.plot(pv,muv,ls=linestyle_unstable, lw=1)

    # plot the eigenvalue
    l0v = df['l0']
    ax_diag_eig.plot(pv, l0v, color=l.get_color())

ax_diag_eig.set_xlabel(r'Cell-cell interaction $\alpha$')
ax_diag_alpha.set_ylabel(r'Order parameter $\mu(\rho)$')
ax_diag_eig.set_ylabel("Largest eigenvalue")

######### profiles
axes_profiles = []
for i, r in df_pointstoshow.iterrows():
    ax = fig.add_subplot(gs_profiles[i//3, i%3])
    axes_profiles.append(ax)

    xu = np.loadtxt(datafolder + 'bifurcation_data/bifurcationdiagram_point_{}_branch{}_point{}.txt'.format(r['diagram'], r['branch_index'], r['point_index']))

    stab = r['l0']<stab_thr
    if stab:
        ax.plot(xu[0], xu[1], color=color_stable, linestyle=linestyle_stable)
    else:
        ax.plot(xu[0], xu[1], color=color_unstable, linestyle=linestyle_unstable)

    ##### Add roman numeral
    ax_diag_alpha.text(r['alpha'], r['mu'], roman_numeral(i+1), ha='center', va='center', bbox=dict(facecolor='w', edgecolor='none', pad=1), fontstyle='italic')
    ax.text(0.03, 0.95, roman_numeral(i+1), ha='left', va='top', bbox=dict(facecolor='w', edgecolor='none', pad=1), fontstyle='italic', transform=ax.transAxes)
    #msize=10
    #ax_diag_alpha.plot(r['alpha'], r['mu'], **data_marker, color=get_color_pattern(r['pattern_group']), markersize=msize)
   
### labels etc
for i, ax in enumerate(axes_profiles):
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    if i%3!=0:
        ax.set_yticklabels([])
    if i>=6:
        ax.set_xticklabels([0, '$\ell$'])
    else:
        ax.set_xticklabels([])
    if i==3:
        ax.set_ylabel(r"Density $\rho$")
    if i==7:
        ax.set_xlabel("Position")

letter_positions = {'a': (0.03, 0.96), 'b':('a', 0.42), 'c':(0.55, 'a') }
for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')


fig.savefig(savefolder + '/supplemental_figures/sfig_bifurcationalpha_all.pdf')
fig.savefig(savefolder + '/supplemental_figures/sfig_bifurcationalpha_all.svg')

