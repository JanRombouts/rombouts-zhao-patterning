from utils import *

### data polarization time
df = pd.read_csv(datafolder + "poltimes_forhistograms.csv")

#### one extra with large value of L only betaf = 1
df2 = pd.read_csv(datafolder + "poltimes_forhistograms_largeLbetaf1.csv")


#########################################################
### Setup figure
fig = plt.figure(figsize=(2*fig_width, 1*fig_width))


topc, bottomc, leftc, rightc = 0.9, 0.15, 0.2, 0.95

Lv = np.sort(df['L'].unique())
betafv = np.sort(df['beta_factor'].unique())
nL = len(Lv)
nbeta = len(betafv)

gs = gridspec.GridSpec(nrows = nbeta, ncols=nL+1, hspace=0.35, wspace=0.5,
        right=rightc, left=leftc, top=topc, bottom=bottomc)

for i, betaf in enumerate(betafv):
    for j, L in enumerate(Lv):
        ax = fig.add_subplot(gs[i,j])
        if i==0:
            bins = np.arange(0,500,20)
        elif i==1:
            bins = np.arange(0,300,20)
        else:
            bins=np.arange(0,1500,50)
        
        ax.hist(df.query("L==@L & beta_factor==@betaf & alpha==20.")['poltime_mode'], bins=bins)

        if i==0:
            ax.set_title('$\ell={:.2f}$'.format(L))
        if j==0:
            ax.text(-0.7, 0.5, r'$\beta = {}\alpha\rho_0$'.format(betaf) + '\n{}'.format(['repelling', 'neutral', 'attracting'][i]), transform=ax.transAxes, va='center', ha='right')
            ax.set_ylabel("Count")

        if i==nbeta-1:
            ax.set_xlabel('Polarization time\n'+r'$[1/\tau]$')

#### one with large L for betaf = 1
axn = fig.add_subplot(gs[1,-1])
axn.hist(df2.query("alpha==20.")['poltime_mode'], bins=np.arange(0, 5000, 250))
axn.set_title('$\ell={:.2f}$'.format(df2['L'].iloc[0]))

axn.set_xlabel('Polarization time\n'+r'$[1/\tau]$')
###########

fig.savefig(savefolder +  'supplemental_figures/sfig_distr_poltime.pdf')
fig.savefig(savefolder +  'supplemental_figures/sfig_distr_poltime.svg')