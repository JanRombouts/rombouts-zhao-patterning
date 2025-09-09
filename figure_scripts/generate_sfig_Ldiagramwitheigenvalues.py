from utils import *


####
##### Supp figure, to go with the bifurcation diagram main figure
#### We show the one parameter diagram with L
### now with extra panel showing eigenvalues



df_bd_L = pd.read_csv(datafolder + 'bifurcation_data/bd_data_L.csv')
stab_thr = 1e-4
df_bd_L['stable'] = df_bd_L['l0']<stab_thr
df_bd_L['pattern_group'] = df_bd_L['mode_k'].apply(classify_pattern_from_mode)

df_bd_L = reformat_branches(df_bd_L)

### manually give some new branch labels that will be colored
### (manual, based on a plot of all branches)
def new_branch_index(s):
    # s is series, this functino should be applied to rows of dataframe
    bi = s['branch_index']
    mu = s['mu']
    if bi == 0:
        if s['stable']:
            return '0a'
        else:
            return '0b'
    if bi == 1 or bi==2 or bi==17:
        if mu>0:
            return '1.a'
        else:
            return '1.b'
    if bi==3 or bi==4 or bi==18:
        if mu>0:
            return '3.a'
        else:
            return '3.b'
    if bi==5 or bi==6 or bi==19:
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
    if bi==11 or bi==12:
        if mu>0:
            return '11.a'
        else:
            return '11.b'
    else:
        return 'other'
df_bd_L['new_branch_index2'] = df_bd_L.apply(new_branch_index, axis=1)
#######

fig = plt.figure(figsize=(2*fig_width, 0.8*fig_width))

gs = gridspec.GridSpec(nrows = 1, ncols = 2, top=0.9, right=0.95, bottom=0.15, wspace=0.4)

### left: the diagram
ax_diag_L = fig.add_subplot(gs[0])
ax_diag_eig = fig.add_subplot(gs[1])
for bi, df in df_bd_L.sort_values('new_branch_index2').groupby('new_branch_index2'): # iterate over the branches
    df=df.sort_values('par')
    pv = df['par']
    muv = df['mu']

    ##### after reformatting the branches, stability and patterngroup should be the same for the whole branch
    stab = df['stable'].iloc[0]
    pg = df['pattern_group'].iloc[0]
    modek = df['mode_k'].iloc[0]

    if stab:
        l, = ax_diag_L.plot(pv,muv,ls=linestyle_stable, lw=1)
    else:
        l, = ax_diag_L.plot(pv,muv,ls=linestyle_unstable, lw=1)

    # on the right plot the eigenvalue
    l0v = df['l0']
    ax_diag_eig.plot(pv, l0v, color=l.get_color())

ax_diag_L.set_ylabel(r'Order parameter $\mu(\rho)$')
ax_diag_eig.set_ylabel(r'Largest eigenvalue')

for ax in [ax_diag_L, ax_diag_eig]:
    ax.set_xlim(0,7)
    ax.set_xlabel(r'System size $\ell$')

ax_diag_eig.set_ylim(-20, 300)

##### panel letters
y_toprow = 0.94
letter_positions = {'a': (0.05, y_toprow), 'b': (0.55, y_toprow)}
for letter, (xl, yl) in letter_positions.items():
    fig.text(xl, yl, f'({letter})', ha='center')


fig.savefig(savefolder + '/supplemental_figures/sfig_bifurcationL_eigenvalues.pdf')
fig.savefig(savefolder + '/supplemental_figures/sfig_bifurcationL_eigenvalues.svg')

