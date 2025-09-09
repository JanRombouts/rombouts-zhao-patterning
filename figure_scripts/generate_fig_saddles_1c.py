from utils import *


######## Figure with the appearance of all the saddle points, also with polarization timing

###### L bifurcation diagram
df_bd_L = pd.read_csv(datafolder + 'bifurcation_data/bd_data_L.csv')

#### stability and pattern classification
stab_thr = 1e-4
df_bd_L['stable'] = df_bd_L['l0']<stab_thr
df_bd_L['pattern_group'] = df_bd_L['mode_k'].apply(classify_pattern_from_mode)

### reformat branches
df_bd_L = reformat_branches(df_bd_L)

##### projection on 2D
vectorfield_data = np.loadtxt(datafolder + 'griddata_2dvectorfield.txt') # x y vx vy
traj_files = glob.glob(datafolder + 'traj_2d_coefficients*')
trajectories = {}
for tf in traj_files:
    traj_ind = int(tf.split('.')[-2].split('_')[-1])
    trajectories[traj_ind] = np.loadtxt(tf)

#### Density trajectories
traj_files_profiles = glob.glob(datafolder + 'simulations/traj_2d_fulldata*')
alldata_profiles = {}
for tf in traj_files_profiles:
    traj_ind = int(tf.split('.')[-2].split('_')[-1])
    if traj_ind!=1 and traj_ind!=2: # keep 2 only to show
        continue
    with h5py.File(tf, 'r') as f:
        alldata_profiles[traj_ind] = (f['sim_bounded']['x'][:], f['sim_bounded']['t'][:], f['sim_bounded']['u'][:]) 



### data polarization time
df_poltime = pd.read_csv(datafolder + "polarization_time.csv")


#########################################################
### Setup figure
fig = plt.figure(figsize=(fig_width, 1.15*fig_width))


topc, bottomc, leftc, rightc = 0.94, 0.08, 0.16, 0.95

gs_length = gridspec.GridSpec(nrows = 2, ncols=1, hspace=0.45, 
        right=0.5, left=leftc, top=topc, bottom=0.3)

gs_pp = gridspec.GridSpec(nrows = 2, ncols=1, hspace=gs_length.hspace, 
        right=rightc, left=gs_length.right+0.1, top=topc, bottom=gs_length.bottom)
ax_diag_L = fig.add_subplot(gs_length[0,0])
ax_poltime = fig.add_subplot(gs_length[1,0]) ### dummy axis

ax_pp = fig.add_subplot(gs_pp[0,0])
ax_equilibria = fig.add_subplot(gs_pp[1,0])

gs_profiles = gridspec.GridSpec(nrows=1, ncols=2, left=leftc, right=rightc-0.12, bottom=bottomc, top = gs_pp.bottom-0.12, wspace=0.3)

axes_traj = [fig.add_subplot(gs_profiles[0,i]) for i in [0,1]]


#############
#### Bifurcation diagram
for bi, df in df_bd_L.groupby('new_branch_index'): # iterate over the branches
    pv = df['par']
    muv = df['mu']
    stab = df['stable'].iloc[0]
    pg = df['pattern_group'].iloc[0]
    modek = df['mode_k'].iloc[0]

    if modek==1.5:
        current_color = color_highlight
    elif stab:
        current_color = color_stable
    else:
        current_color = color_unstable
    if stab:
        ax_diag_L.plot(pv,muv,color=current_color, ls=linestyle_stable)
    else:
        ax_diag_L.plot(pv,muv,color=current_color, ls=linestyle_unstable)

ax_diag_L.set_ylabel(r"Order parameter $\mu(\rho)$", labelpad=0)
#ax_diag_L.set_xlabel("Normalized domain\n"+r"length $\ell$")
## trim the x axis. 
ax_diag_L.set_xlim((0, 7))
ax_diag_L.tick_params(axis='y', pad=0)

################### Add the approximations for transitions points from LSA
ell_transition = np.loadtxt(datafolder + "approx_L_transitions_alpha20_rho00.4.txt")
for ell in ell_transition:
    #ax_diag_L.plot(ell, 0, 'x', color='tab:red', markersize=2)
    ax_diag_L.annotate('', (ell, 0.01), (ell, 0), arrowprops=dict(arrowstyle='<|-', mutation_aspect=1., color='tab:red'))



#######
#### Polarization time, mean and std

x0, y0, w, h = ax_poltime.get_position().bounds
ax_poltime.set_axis_off()

ax_poltime_mean = fig.add_axes((x0, y0+h/2, w, h/2), sharex=ax_diag_L)
ax_poltime_std = fig.add_axes((x0, y0, w, h/2),  sharex=ax_diag_L, sharey=ax_poltime_mean)


alpha=20
for betafactor,la in zip([0.5, 1, 1.5],['Repelling','Neutral', 'Attracting']) :
    fdat = df_poltime.query("alpha=={} & beta_factor=={}".format(alpha,betafactor)).reset_index()
    
    l, = ax_poltime_mean.semilogy(fdat['L'], fdat['poltime_mode_mean'], label=la, color=cmap_eigenfunctions(betafactor/2))
    ax_poltime_std.semilogy(fdat['L'], fdat['poltime_mode_std'], color=l.get_color(), ls = ':', label=la)


ax_poltime_mean.legend(handlelength=0, fontsize=6, bbox_to_anchor=(1, -0.05), borderpad=0, labelcolor='linecolor',
labelspacing=0.2, loc='lower right')

for ax in [ax_poltime_mean, ax_poltime_std]:
    ax.minorticks_off()
    ax.set_xticks([2,3,4,5,6])
    ax.set_xticklabels(ax.get_xticks())
    ax.set_yticks([1, 1000])

ax_poltime_std.set_xlabel("System size $\ell$", labelpad=1)

#ax_scalingL_mean.text(-0.43, -0.7, 'Polarization time', rotation=90, clip_on=False, transform=ax_scalingL_mean.transAxes)
ax_poltime_mean.set_ylabel(r"Mean", labelpad=0)
ax_poltime_std.set_ylabel(r"St.dev.", labelpad=0)
ax_poltime_std.tick_params(axis='y', pad=0)
ax_poltime_mean.tick_params(axis='y', pad=0)

### extra text 
ax_poltime.text(-0.37, 0.5, r'Polarization time $t/\tau$', transform = ax_poltime.transAxes, ha='center', va='center', rotation=90)

##### add vertical lines crossing the left panels (critical lengths)
ell_values = [1.3, 2.3, 3.5] # manually determined

for ell in ell_values:
    #### bottom of line:
    bottom_co = fig.transFigure.inverted().transform((ax_poltime_std.transData.transform((ell, ax_poltime_std.get_ylim()[0]))))
    top_co = fig.transFigure.inverted().transform((ax_diag_L.transData.transform((ell,ax_diag_L.get_ylim()[1]))))
    fig.add_artist(mlines.Line2D(*zip(bottom_co, top_co), color='gray', alpha=0.25))


ax_poltime_std.set_xticks([0] + ell_values + [5])
ax_poltime_std.set_xticklabels([0, r'$\ell_\mathrm{PF}$', r'$\ell_1^*$', r'$\ell_2^*$', 5 ])
plt.setp(ax_diag_L.get_xticklabels() + ax_poltime_mean.get_xticklabels(), visible=False)


###############################

#### Phase plane
######### 2D representation
x, y, vx, vy = vectorfield_data.T
speed = (vx**2 + vy**2)**0.5*0+1
ax_pp.quiver(x, y, vx/speed, vy/speed, color='gray', scale=3., width=.005, angles='xy')
#ax_pp.set_aspect(1)

### Add trajectories

#### get speed to normalize with
maxvs=[]
for tr in trajectories.values():
    tvx = tr[:,2]
    tvy = tr[:,3]
    tspeed = np.sqrt(tvx**2 + tvy**2)
    maxvs.append(np.max(tspeed))
v_norm = np.max(maxvs)

## now plot
for tr in trajectories.values():
    tx = tr[:,0]
    ty = tr[:,1]
    tvx = tr[:,2]
    tvy = tr[:,3]
    tspeed = np.sqrt(tvx**2 + tvy**2)
    tspeed_norm = mcolors.Normalize(vmin=0, vmax=v_norm)

    ## line collection
    points = np.array([tx, ty]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    # Create a LineCollection
    lc = LineCollection(segments, cmap='Reds', norm=tspeed_norm)
    lc.set_array(tspeed)
    ax_pp.add_collection(lc)


### small axis for colorbar
data_to_fig = ax_pp.transData + fig.transFigure.inverted()
axc_x, axc_y = data_to_fig.transform((0.37, 0.6))
ax_for_colorbar = fig.add_axes([axc_x, axc_y, 0.1, 0.01])
cbspeed = colorbar_trajspeed = fig.colorbar(lc, cax = ax_for_colorbar, orientation='horizontal')
cbspeed.set_label(r"speed $[1/\tau]$", fontsize=6)
ax_for_colorbar.set_xticks([0, np.max(tspeed)])
ax_for_colorbar.set_xticklabels([0, 'max'], fontsize=6)
ax_for_colorbar.tick_params(axis='x', pad=0)
ax_for_colorbar.xaxis.set_label_coords(0.5, 3.5)

### add markers and small axes and profiles

### create four small axes
x0, y0, w, h = ax_equilibria.get_position().bounds
ax_equilibria.set_axis_off()
gs_equi = gridspec.GridSpec(nrows = 2, ncols = 2, left=x0, right=x0+w, bottom=y0, top=y0+h, hspace=0.3)

axes_eq = [fig.add_subplot(gs_equi[i,j]) for (i, j) in [(0,0), (0,1), (1,0), (1,1)]]

### ordering: left polar, right polar, top side aggregates, bottom middle aggregate
for co_point, axn, marker_fillstyle, pg, axtitle in zip(\
    [ (0.59, 0.168), (-0.59, 0.168), (0,0.54), (0, -0.54), (0,0)],\
     axes_eq, ['full','full', 'left', 'left'], [1, 1, 3, 2], ['left polar','right polar', 'saddle', 'saddle']):
    ax_pp.plot(*co_point, 'o',color=get_color_pattern(pg), markersize=4, \
               markeredgewidth=.3,fillstyle=marker_fillstyle)
    ## add a small axis with the corresponding profile
    xv = np.linspace(0, 1, 50)
    phi1 = np.cos(np.pi*xv)
    phi2 = np.cos(2*np.pi*xv)
    a1, a2 = co_point
    uv = 0.5 + a1*phi1 + a2*phi2
    
    axn.plot(xv,uv, color=get_color_pattern(pg))
    axn.set_xticks([])
    axn.set_yticks([])
    axn.set_xlim(-0.1,1.1)
    axn.set_ylim(uv.min()-0.2, uv.max()+0.2)
    for spine in axn.spines.values():
        spine.set_linewidth(0.5)  
        spine.set_color('gray')

    ### title text + marker
    xpmarker=0
    tt = axn.text(xpmarker+0.15, 1.08, axtitle, fontsize=6, transform=axn.transAxes, ha='left', va='center')
    
    axn.plot(xpmarker,1.08, 'o',color=get_color_pattern(pg), markersize=4, \
               markeredgewidth=.3,fillstyle=marker_fillstyle, clip_on=False, transform=axn.transAxes)
    
## add marker at zero
ax_pp.plot(0,0, 'o',color='k', markersize=5, \
               markeredgewidth=.3,fillstyle='full', markerfacecolor='w')

#ax_pp.set_axis_off()
ax_pp.spines.top.set_visible(0)
ax_pp.spines.right.set_visible(0)

ax_pp.set_xticks([-0.5, 0, 0.5])
ax_pp.set_yticks([-0.5, 0, 0.5])
ax_pp.set_xlabel("Cosine mode $a_1$", labelpad=0)
ax_pp.set_ylabel("Cosine mode $a_2$", labelpad=-70, rotation=0, loc='top')
ax_pp.tick_params(axis='y', pad=0)

ax_pp.set_xlim(-0.7, 0.7)
ax_pp.set_ylim(-0.7, 0.7)


#### Trajectory profile plots on the bottom
snapshots = 4
for i, (ax, data) in enumerate(zip(axes_traj, alldata_profiles.values())):
    x, t, u = data

    for j, f in enumerate(np.linspace(0, len(t)//6*5-1, snapshots).astype(int)):
        ax.plot(x, u[f], color=plt.cm.Greys((j+1)/(snapshots+1)), label='{:.1f}'.format(t[f]))


    L = x[-1] + 0.5*np.diff(x)[0]
    ax.set_xlim([0, L])
    ax.set_xticks([0, L])

    ax.spines.top.set_visible(0)
    ax.spines.right.set_visible(0)
    ax.spines.bottom.set_position(('outward', 2))
    ax.spines.left.set_bounds((0, 1))
    ax.set_xlabel('Position $x/\sigma$', labelpad=-4)
axes_traj[0].set_ylabel("Density\n"+r"$\rho/\rho_\mathrm{max}$")
axes_traj[1].legend(title=r'Time $t/\tau$', labelspacing=0.3, borderpad=0, \
bbox_to_anchor=(1, 0), loc='lower left', handlelength=1, fontsize=6, title_fontsize=6)

#### add little arrow to the bifurcation diagram with the L for this panel
ax_diag_L.annotate('', (L, -0.04), (L, ax_diag_L.get_ylim()[0]), arrowprops=dict(arrowstyle='<-'), annotation_clip=False)

### add markers for trajectories
for pos_pp, ax, mmcolor in zip([(0.32, 0.85), (0.73, 0.2)], axes_traj, ['k','tab:red']):
    
    posx, posy = pos_pp
    ax_pp.plot(posx, posy, **data_marker, color=mmcolor, transform=ax_pp.transAxes)
    ax.plot(0., 1.05, **data_marker, color=mmcolor, transform=ax.transAxes, clip_on=False)

# panel letters
letter_positions = {'a': (0.04, 0.97), 'b': (0.53, 'a'), 'c': ('a', 0.21)}
for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')


###########

fig.savefig(savefolder +  'main_figures/fig_saddles.pdf')
fig.savefig(savefolder +  'main_figures/fig_saddles.svg')