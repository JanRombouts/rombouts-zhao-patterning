from utils import *


##### Introduction figure
### contains diagrams of interactions, plot of kernel, and two example simulations

### load data for kymo simulation
simulation_datafiles = glob.glob(datafolder + "simulations/data_simulation_intro_*")
print(simulation_datafiles)
simulation_data = []
simulation_domainlength = []
for filename in simulation_datafiles:
    with h5py.File(filename, 'r') as f:
        simulation_data.append((f['sim_bounded']['x'][:] , f['sim_bounded']['t'][:], f['sim_bounded']['u'][:]))
        simulation_domainlength.append(f['sim_bounded'].attrs['L'])

## sort from shortest to longest
ind_sorted = np.argsort(simulation_domainlength)
simulation_data = [simulation_data[i] for i in ind_sorted]
simulation_domainlength = [simulation_domainlength[i] for i in ind_sorted]


#######################

fig = plt.figure(figsize=(fig_width, .6*fig_width))

topc, bottomc, leftc, rightc = 0.93, 0.13, 0.15, 0.97

gs_sketches = gridspec.GridSpec(nrows = 1, ncols=2, right=rightc,\
                       top=0.95, bottom=0.57, left=leftc, wspace=0.3, width_ratios=[1, 0.65])

gs_simulation = gridspec.GridSpec(nrows = 1, ncols=2,right=0.8,\
                       top=gs_sketches.bottom-0.22, bottom=bottomc, \
                        left=leftc, wspace=.2, width_ratios = simulation_domainlength)


### axis for sketch of model and interaction function
ax_model = fig.add_subplot(gs_sketches[0])
ax_kernel_cc = fig.add_subplot(gs_sketches[1])


## simulations
axes_simulations = [fig.add_subplot(gs_simulation[i]) for i in range(len(simulation_data))]
ax_model.set_axis_off()

##### drawing of interaction kernel
xv = np.linspace(0, 4, 50)
yv = np.exp(-xv)
ax_kernel_cc.plot(xv, yv, color='k')

ax_kernel_cc.set_ylabel("Interaction\nstrength $g(s/\sigma)$")
ax_kernel_cc.set_xlabel(r"Distance $s/\sigma$", labelpad=0)

ax_kernel_cc.set_yticks([0,1])
ax_kernel_cc.set_xticks(np.arange(5))

for s in ['top', 'right']:
    ax_kernel_cc.spines[s].set_visible(0)

######### show simulation data: snapshots on top of each other

for i, d in enumerate(simulation_data):
    x, t, u = d
    ax = axes_simulations[i]


    snapshots = 5
    dt = t[1]-t[0]
    tmax = np.max(t)
    timepoints = np.array([0, tmax/3, tmax/3*2, tmax])
    frames = (timepoints / dt).astype(int)# np.linspace(0, u.shape[0]-1, snapshots).astype(int)
    for j, frame in enumerate(frames):
        ax.plot(x, u[frame], color='k', alpha= (j+1)/(snapshots+1), label='{:.2f}'.format(timepoints[j]))
    
    ax.set_xlabel("Position $x/\sigma$", labelpad=-3)
    for s in ['top','right']:
        ax.spines[s].set_visible(0)
    if i==1:
        ax.spines['left'].set_visible(0)
        ax.set_yticks([])
    else:
        ax.set_yticks([0,1])
        ax.set_ylim([0,1])
        ax.spines['left'].set_position(('outward', 5))

    L = x[-1]+x[0] #### because of discretization
    ax.set_xticks([0,L])
    ax.set_xlim([0,L])

    

    if i==1:
        ax.legend(title=r'Time $t/\tau$', bbox_to_anchor=(1, 0.5), loc='center left', fontsize=6, title_fontsize=6)

axes_simulations[0].set_ylabel(r"Density $\rho/\rho_\text{max}$", labelpad=0)


# panel letters
letter_positions = {'a': (0.04, 0.94), 'b': (0.5, 'a'), \
'c': ('a',0.45)}

for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')

fig.savefig(savefolder + "/main_figures/fig_intro_withoutsketches.svg")
fig.savefig(savefolder + "/main_figures/fig_intro_withoutsketches.pdf")
