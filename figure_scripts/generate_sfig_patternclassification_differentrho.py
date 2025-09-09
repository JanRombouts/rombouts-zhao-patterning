##### Figure that accompanies Fig. 4 e

#### we show the pattern groups per length for the experiments, for different selected density ranges

from utils import *


allbardata, fit_data = load_bar_data()


############## Final profile classification experiments
### Load data
dfc, dfc_r = load_pattern_classification()
# drop the 40 microns
dfc_r.drop(dfc_r[dfc_r['L']==40].index, inplace=True)

##in the main figure we show density from 0.25 to 0.55. Here we show all, and lower/higher
rho_threshold_low = 0.25
rho_threshold_high = 0.55
dfc_r_low = dfc_r.drop(dfc_r.index[dfc_r['avg_rho']>=rho_threshold_low])
dfc_r_high = dfc_r.drop(dfc_r.index[dfc_r['avg_rho']<=rho_threshold_high])



##### For the colors
### the colormap

values = [0,1,2,3,4]
cmap = mcolors.ListedColormap([get_color_pattern(v) for v in values])
bounds = [v-0.5 for v in values] + [values[-1]+0.5]
norm = mcolors.BoundaryNorm(bounds, cmap.N)

colorlist = cmap(norm(values))

pattern_labels = {0: "No pattern", 1: 'Polarized', 2: 'Middle peak', 3: 'Two peaks',4: "More peaks"}

#### make figure

fig = plt.figure(figsize=(2*fig_width, 0.8*fig_width))

gs = gridspec.GridSpec(nrows = 1, ncols=3, wspace=0.4,right=0.97, left=0.09, top=0.8, bottom=0.15)

for i, df, ax_title in zip(range(3), [dfc_r, dfc_r_low, dfc_r_high], ['All', r'$\rho_0 < {}$ cells/µm'.format(rho_threshold_low),
                                                                       r'$\rho_0 > {}$ cells/µm'.format(rho_threshold_high)]):
    #### get the counts per pattern type
    df_exp_counts = df.groupby(['L','pattern_group'])[['avg_rho']].count()

    #### add zeros if there are pattern groups that don't exist for given length
    newidx = pd.MultiIndex.from_product([df_exp_counts.index.get_level_values(0).unique(), range(0, max_pattern_n+1)], names = df_exp_counts.index.names)
    df_exp_counts = df_exp_counts.reindex(newidx, fill_value=0)
    df_exp_counts['fraction'] = df_exp_counts['avg_rho'] / df_exp_counts.groupby('L')['avg_rho'].transform('sum')

    ax = fig.add_subplot(gs[i])

    Lv = df_exp_counts.index.get_level_values(0).unique()
    bottoms = np.zeros(len(Lv))
    for i, (pat, dff) in enumerate(df_exp_counts.groupby('pattern_group')):
        ax.bar(dff.index.get_level_values(0), dff['fraction'], bottom=bottoms, width=16, label=pattern_labels[pat], color=colorlist[i])
        bottoms+=dff['fraction'].values


    ### add number of repeats for each length
    reps_per_length = df.groupby("L").count().iloc[:,0].to_dict()
    for L,n in reps_per_length.items():
        if L<45:
            continue
        else:
            ax.text(L, 1.03, f'{n}', clip_on=False, ha='center')
    # add a small 'n'
    ax.text(80-20, 1.03, 'n:', clip_on=False, ha='center')

    ax.set_title(ax_title, y=1.07)

    ax.set_ylabel("Fraction of experiments")
    ax.set_xlabel("Bar length (µm)")
    ax.set_xlim(60, 340)
    ax.set_xticks(df_exp_counts.index.get_level_values(0).unique())
    plt.setp(ax.get_xticklabels()[1], visible=False)

    ### add legend

## legend
pgroups_legend = [0, 1,2,3,4]
patches_legend = [mpatches.Rectangle([0,0],0,0, fc=get_color_pattern(pp)) for pp in pgroups_legend]
labels_legend = [pattern_labels[pp] for pp in pgroups_legend]

fig.axes[1].legend(patches_legend, labels_legend, 
fontsize=8, ncols=len(pgroups_legend), bbox_to_anchor=(0.5, 1.15), columnspacing=.7, handlelength=1, handletextpad=0.25, loc='lower center')

# panel letters

def gtl(ax):
    #### returns (x, y) of top left point of the axis
    axpos = ax.get_position()
    return axpos.bounds[0], axpos.bounds[1] + axpos.bounds[3]

letter_positions = {'a': (gtl(fig.axes[0])[0]-0.05, 0.9), 'b': (gtl(fig.axes[1])[0]-0.05, 'a'), 'c': (gtl(fig.axes[2])[0]-0.05, 'a') }
for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')

fig.savefig(savefolder + '/supplemental_figures/sfig_exp_classification_diffrho.pdf')
fig.savefig(savefolder + '/supplemental_figures/sfig_exp_classification_diffrho.svg')
