from utils import *


##########
######## Supplemental Figure
#### showing how we calculate density profiles

### use one of the bars of the intro figure
allbardata, fit_data = load_bar_data()

bar_to_use  = ('20240320','23')
exp, bar = bar_to_use

L = allbardata.loc[(exp, bar), 'Length_micron']

#### load all the info
delta_frame_to_use = 160
frame_to_show = int(allbardata.loc[(exp, bar), 'frame0'])+delta_frame_to_use
print(frame_to_show)

#### Load the Brachyury and H2B channels
im_bra = skimage.io.imread(datafolder + f"experimental_data/{exp}_bar{bar}_brachyury_maxz2-4.tif")[frame_to_show,:,:]
im_h2b = skimage.io.imread(datafolder + f"experimental_data/{exp}_bar{bar}_h2b_maxz2-4.tif")[frame_to_show,:,:]

### load cell positions
cell_positions = pd.read_csv(datafolder + f"experimental_data/cellpositions_{exp}_bar{bar}.csv").query('frame==@frame_to_show')
## load density profile
density_profile = np.loadtxt(datafolder + f"experimental_data/densityprofile_{exp}_bar{bar}.txt")

###########
#### Build the figure
fig = plt.figure(figsize=(1.5*fig_width, 1.7*fig_width))

leftc, rightc, bottomc, topc  = 0.1, 0.97, 0.08, 0.96

gs_all = gridspec.GridSpec(nrows=5, ncols=2, left=leftc, right=rightc, top=topc, bottom=bottomc,hspace=.5)

ax_im_bra = fig.add_subplot(gs_all[0,0])
ax_im_h2b = fig.add_subplot(gs_all[0,1])

ax_deadalive = fig.add_subplot(gs_all[1,1])
ax_pos_bra = fig.add_subplot(gs_all[2,0])
ax_pos_highlow = fig.add_subplot(gs_all[3,0])

for ax in [ax_im_bra, ax_im_h2b, ax_deadalive, ax_pos_bra, ax_pos_highlow]:
    ax.set_aspect(1, anchor='W')
    ax.set_ylim(0,40)
    ax.set_yticks([0,40])
    ax.set_xlim(0, L)
    ax.set_xticks([0,L])

ax_densityprofile = fig.add_subplot(gs_all[4,0])

#### images
horizontal = allbardata.loc[(exp, bar), 'Orientation']=='horizontal'

vmin_bra, vmax_bra = df_intensityvalues.loc[(exp, delta_frame_to_use), ['Bra_p0.5', 'Bra_p99.5']]
vmin_h2b, vmax_h2b = df_intensityvalues.loc[(exp, delta_frame_to_use), ['H2B_p0.5', 'H2B_p99.5']]

x = np.arange(im_bra.shape[1])/pxpermicron
y = np.arange(im_bra.shape[0])/pxpermicron 

if horizontal:
    ppbra=ax_im_bra.pcolormesh(x, y, im_bra, cmap=cmap_bra, shading='nearest', rasterized=True, vmax=vmax_bra, vmin=vmin_bra)
    pph2b=ax_im_h2b.pcolormesh(x, y, im_h2b, cmap=cmap_h2b, shading='nearest', rasterized=True, vmax=vmax_h2b, vmin=vmin_h2b)

else: # switch x and y
    ppbra=ax_im_bra.pcolormesh(y, x, im_bra.T, cmap=cmap_bra, shading='nearest', rasterized=True, vmax=vmax_bra, vmin=vmin_bra)
    pph2b=ax_im_h2b.pcolormesh(y, x, im_h2b.T, cmap=cmap_h2b, shading='nearest', rasterized=True, vmax=vmax_h2b, vmin=vmin_h2b)

#### dead/alive
df_alive = cell_positions.query("Class==2")
df_dead = cell_positions.query("Class==1")

row_offset = allbardata.loc[(exp, bar), 'min_row']
col_offset = allbardata.loc[(exp, bar), 'min_col']

xalive = (df_alive['x']-col_offset)/pxpermicron
yalive = (df_alive['y']-row_offset)/pxpermicron
xdead = (df_dead['x']-col_offset)/pxpermicron
ydead = (df_dead['y']-row_offset)/pxpermicron

if horizontal:
    ax_deadalive.plot(xalive, yalive, '.k', markersize=3)
    ax_deadalive.plot(xdead, ydead, 'x', color='tab:red', markersize=3)
else:
    ax_deadalive.plot(yalive, xalive, '.k', markersize=3)
    ax_deadalive.plot(ydead, xdead, 'x', color='tab:red', markersize=3)

### cell positions with Brachyury level (scatter plot)
if horizontal:
    ppposbra = ax_pos_bra.scatter(xalive, yalive, c = df_alive['Bra_ratio'], s=3, cmap=cmap_bra)
else:
    ppposbra = ax_pos_bra.scatter(yalive, xalive, c = df_alive['Bra_ratio'], s=3, cmap=cmap_bra)

## cell positions high low
df_high = df_alive.query('type=="high"')
df_low = df_alive.query('type=="low"')
xhigh = (df_high['x']-col_offset)/pxpermicron
yhigh = (df_high['y']-row_offset)/pxpermicron
xlow = (df_low['x']-col_offset)/pxpermicron
ylow = (df_low['y']-row_offset)/pxpermicron

if horizontal:
    ax_pos_highlow.plot(xlow, ylow, '.', markersize=3, color="#ebebeb")
    ax_pos_highlow.plot(xhigh, yhigh, 'k.', markersize=3)

else:
    ax_pos_highlow.plot(ylow, xlow, '.', markersize=3, color="#ebebeb")
    ax_pos_highlow.plot(yhigh, xhigh, 'k.', markersize=3)


### density profile
x = density_profile[0,:]
rho = density_profile[1+frame_to_show, :]
ax_densityprofile.plot(x, rho, '-ok', markersize=3)
ax_densityprofile.set_ylim(-0.1, 1.1)

#### do a check..
#rho_check,_ = np.histogram(yhigh, bins=np.linspace(0, L, int(L/10)+1))

#ax_densityprofile.plot(x, rho_check/10, 'tab:red')


### labels etc.
for ax in [ax_im_bra, ax_pos_bra, ax_pos_highlow]:
    ax.set_ylabel("Position (µm)")
for ax in fig.axes[:-1]:
    ax.set_xlabel("Position (µm)", labelpad=-5)

ax_densityprofile.set_xlabel("Position (µm)")
ax_densityprofile.set_ylabel("Density (cells/µm)")

### color bars
for ax, pp, label in zip((ax_im_bra, ax_im_h2b,ax_pos_bra),\
                                      (ppbra, pph2b, ppposbra), (r'Bra/Bra$_0$', 'H2B/H2B$_0$', 'Bra/H2B (norm.)')):
    x0,y0,w,h = ax.get_position().bounds
    vmin, vmax =pp.get_clim()
    cbarax = fig.add_axes([x0+w/2, y0+h+0.01, w/2, 0.01])
    fig.colorbar(pp, cax= cbarax, label=label, orientation='horizontal', ticks=[vmin, vmax])
    cbarax.set_xticklabels([0,1])
    cbarax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, pad=0)
    cbarax.xaxis.set_label_coords(0.5, 3.5)


### text
ax_im_bra.text(0.05, 0.95, 'Brachyury', fontweight='bold',color= plt.cm.plasma(1.), \
                        ha='left', va='top', transform=ax_im_bra.transAxes)
ax_im_h2b.text(0.05, 0.95, 'H2B', fontweight='bold',color= 'w', \
                       ha='left', va='top', transform=ax_im_h2b.transAxes)

### legend for the dead/alive
ax_deadalive.legend(['Alive', 'Dead'], bbox_to_anchor=(0.5, 1), loc='lower center', ncols=2, handletextpad=0.1)
### legend for the high/low
ax_pos_highlow.legend(['Low Bra', 'High Bra'], bbox_to_anchor=(0.5, 1), loc='lower center', ncols=2, handletextpad=0.1)



fig.savefig(savefolder + '/supplemental_figures/sfig_pipeline_withoutannotations.pdf')
fig.savefig(savefolder + '/supplemental_figures/sfig_pipeline_withoutannotations.svg')
