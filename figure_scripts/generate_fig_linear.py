from utils import *



########### Linear stability with eigenvalues


########################## Setup for the eigenfunction plots 
def zeta(z, L, n):
    return np.exp(L*z/2)/(z-1) - (-1)**n *np.exp(-L*z/2)/(z+1)
def get_eigenfunction(w1,w2,l,alpha,L,n):
    ## returns a function (callable)
    z1 = np.sqrt(w1+0*1j)
    z2 = np.sqrt(w2 + 0*1j)
    if n==0:
        z = np.sqrt(1-2*alpha +0*1j)
        return lambda x: zeta(z,L,0) + np.exp(z*(x-L/2)) + np.exp(-z*(x-L/2))
    def rf(z,n):
        if n % 2==0:
            return 1
        elif abs(z.imag)>1e-10:
            return 1j
        else:
            return 1

    A = 1
    B = -A *rf(z1,n) * zeta(z1,L,n)/rf(z2,n)/zeta(z2,L,n)
    return lambda x: A*rf(z1, n)*(np.exp(z1*(x-L/2)) + (-1)**n*np.exp(-z1*(x-L/2))) + B*rf(z2,n)*(np.exp(z2*(x-L/2)) + (-1)**n*np.exp(-z2*(x-L/2)))

## Load in the branch data for the exact eigenvalues
with open(datafolder + 'linear_stability/branches_exact_eigenvalues_exp.pkl', 'rb') as f:
    b1 = dill.load(f)
#### lists to 2D arrays
for n,(X,v) in b1.items():
    X2 = np.array(X)
    v2 = np.array(v)
    b1[n] = (X2, v2)

##### manually add the zero branch (dummy values for w1, w1,...)
XXX0 = b1[1][0].copy()
vvv0 = b1[1][1].copy()
XXX0[:,-2]=0
b1[0] = (XXX0, vvv0)



########################

fig = plt.figure(figsize=(fig_width, .9*fig_width))

topc, bottomc, leftc, rightc = 0.93, 0.03, 0.12, 0.97


n_ef = 4 # how many eigenfunctions on one branch do we show
n_b = 4 # of how many branches

### eigenvalue loci
gs_omega = gridspec.GridSpec(nrows = 1, ncols=1, hspace=0.3, wspace=0.3, \
        right=0.48, left=leftc, top=topc, bottom=0.53)
ax_ev = fig.add_subplot(gs_omega[0])
ax_zoom = ax_ev.inset_axes([0.04, 0.6, 0.35, 0.35])

##### regime diagram
gs_regimediag = gridspec.GridSpec(nrows = 1, ncols=1, \
        right=rightc, left=gs_omega.right+0.14, top=topc, bottom=gs_omega.bottom)

ax_regimediag = fig.add_subplot(gs_regimediag[0])

gs_ef = gridspec.GridSpec(nrows = n_ef, ncols=n_b, figure=fig, hspace=0.3, wspace=0.2, \
                          top=gs_omega.bottom-0.18, left=leftc+0.05, bottom=bottomc, right=rightc)

## axes for eigenfunctions
axes_ef = [[] for i in range(n_ef)]
for i in range(n_ef):
    for c in range(n_b):
        axes_ef[i].append(fig.add_subplot(gs_ef[i,c]))
    
axes_ef = np.array(axes_ef)

## alphas at which we draw an eigenfunction
alphabar_toshow = [0., 4, 8,12]

L=5

maxmode_toshow = 7

maxn_cmap = 10
for n in range(maxmode_toshow+1): ### iterate over branches
    skip=10 # careful with this!
    ##### load data from the branch structures
    XXX = np.array(b1[n][0])[::skip,:]
    vvv = np.array(b1[n][1])[::skip,:]
    alphabarv = XXX[:,-1]
    omegav = XXX[:,-2]



    l,=ax_ev.plot(alphabarv,omegav,color=cmap_eigenfunctions(n/maxn_cmap), label=r'$n={}$'.format(n))

    # add results from perturbation
    l0 = -(n*np.pi/L)**2
    l1 = -2*l0/(1-l0)**2*(1-l0 - 2/L*(1 - (-1)**n*np.exp(-L)))
    ax_ev.plot(alphabarv, l0+alphabarv*l1, ':', color=l.get_color())

    ### inset: only selected modes
    if n==1 or n==3:
        l,=ax_zoom.plot(alphabarv,omegav,color=cmap_eigenfunctions(n/maxn_cmap), label=r'$n={}$'.format(n))
    
        # add results from perturbation
        l0 = -(n*np.pi/L)**2
        l1 = -2*l0/(1-l0)**2*(1-l0 - 2/L*(1 - (-1)**n*np.exp(-L)))
        ax_zoom.plot(alphabarv, l0+alphabarv*l1, ':', color=l.get_color())


    # add some eigenfunctions for some branches
    if n<n_b:
        xx = np.linspace(0, L, 100)
        for j in range(n_ef):
            ## get index
            index_alpha = np.nonzero(alphabarv>alphabar_toshow[j])[0][0]
            ef = get_eigenfunction(*XXX[index_alpha], L, n)

            yv = ef(xx)
            #normalize and align directions
            NN = np.max(abs(yv))
            dx = xx[1]-xx[0]
            ## a little bit ad hoc... center of mass for odd and peak in middle for even
            if n % 2 == 0:
                aa = np.sign(np.sum(yv*np.cos(2*np.pi/L*xx)*dx))
            else:
                aa = np.sign(np.sum(yv*xx*dx)-L/2)
            axes_ef[j][n].plot(xx, yv/NN*aa, color=l.get_color())

        axes_ef[0,n].set_title("$n={}$".format(n), color=l.get_color(), y=0.8)
    
#### add zero eigenfunction manually

ax_ev.axhline(0, color=cmap_eigenfunctions(0))
ax_ev.set_xlabel(r"Effective interaction strength $\bar\alpha$", labelpad=1)

ax_ev.set_ylim(-5, 20)
ax_ev.set_ylabel(r"Eigenvalue $\omega_n$", labelpad=1)
ax_ev.set_xticks(range(0, 20, 5))

ax_zoom.set_xlim(1,5)
ax_zoom.set_ylim(-1,3)
ax_zoom.set_xticks([])
ax_zoom.set_yticks([])

ax_ev.indicate_inset_zoom(ax_zoom)


###### edit the eigenfunction axes
for ax in axes_ef.flatten():
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylim(-1.5, 1.5)
    ax.axhline(0, color='gray', lw=.5)

## add alpha values on the left and arrows to the top diagram
for i in range(n_ef):
    axes_ef[i,0].text(-0.45, 0.5, r'$\bar\alpha={:.0f}$'.format(alphabar_toshow[i]), ha='center', va='center', transform=axes_ef[i,0].transAxes)
    yarrow = 0
    ax_ev.annotate('', [alphabar_toshow[i], yarrow], [alphabar_toshow[i], yarrow-3], arrowprops=dict(arrowstyle='->'), clip_on=False)



######### regime diagram

#### Lines from theory
alphav = np.linspace(0,15,201)
ellv = np.linspace(0,7,201) # multiples of sigma in scaled units

aa,LL = np.meshgrid(alphav,ellv)
# draw them using parametrized curves
ncurves = 20
for n in range(1, ncurves):
    thetav = np.linspace(-np.pi, 0, 200)
    xv = -1j *(1+np.exp(1j*thetav))/(1-np.exp(1j*thetav))
    Lv = (thetav - n*np.pi)/xv
    av = 0.5*(xv**2+1)
    l, = ax_regimediag.plot(Lv,av, label='$n={}$'.format(n), color = cmap_eigenfunctions(n/maxn_cmap))
    if 1:
        # add perturbation
        l0v = -(n*np.pi/LL)**2
        l1v = -2*l0v/(1-l0v)**2*(1-l0v - 2/LL*(1 - (-1)**n*np.exp(-LL)))
        
        lamdav = l0v+aa*l1v
        ax_regimediag.contour(LL,aa,lamdav,levels=[0], linestyles=[':'], colors=[l.get_color()], alpha=.5)

    #### shading below first line
    if n==1:
        ind_sh = (Lv>0.001) & (Lv<7)
        ax_regimediag.fill_between(Lv[ind_sh], av[ind_sh], color=get_color_pattern(0), alpha=0.5)


## no patterns text
ax_regimediag.text(0.05, 0.05, 'No\npatterns', transform=ax_regimediag.transAxes, fontsize=6)

ax_regimediag.set_ylabel(r"Eff. interaction strength $\bar\alpha$", labelpad=0)
ax_regimediag.set_xlabel(r'System size $\ell$', labelpad=1)

ax_regimediag.set_xlim(ellv.min(), ellv.max())
ax_regimediag.set_ylim(alphav.min(), alphav.max())
ax_regimediag.set_yticks([0, 5, 10, 15])
ax_regimediag.set_xticks([0, 2, 4, 6])



#### Add two legends
legend_y = 1.08

### add the word 'modes'
mode_t = ax_ev.text(0, legend_y, 'Modes:', transform=ax_ev.transAxes, fontsize=6, ha='left', va='center')

handles, labels = ax_ev.get_legend_handles_labels()
dummy = mlines.Line2D([], [], color='w', alpha=0)
legend_modes = ax_ev.legend([handles[0], handles[1]], [labels[0] + ',', labels[1] + ', ...'],ncols=2, columnspacing=0.5, \
                              bbox_to_anchor=(mode_t.get_position()[0]+0.18, legend_y), loc='center left', handlelength=1, fontsize=6, handletextpad=0.4)
ax_ev.add_artist(legend_modes)

#### legend with black lines making distinction exact/approximate
ls = mlines.Line2D([], [], color='k')
ld = mlines.Line2D([], [], color='k', linestyle='dotted')

ax_regimediag.legend([ls, ld], ['Exact', 'Perturbation'], ncols=2, columnspacing=0.5,\
                              bbox_to_anchor=(0, legend_y), loc='center left', handlelength=1.5, fontsize=6)


### panel letters

# panel letters
letter_positions = {'a': (0.03, 0.96), 'b':('a', 0.38) }
for letter, (xl, yl) in letter_positions.items():
    while isinstance(xl, str):
        xl = letter_positions[xl][0]
    while isinstance(yl, str):
        yl = letter_positions[yl][1]
    fig.text(xl, yl, f'({letter})', ha='center')


###### Add the fits
fig.savefig(savefolder + 'main_figures/fig_linear.svg')
fig.savefig(savefolder + 'main_figures/fig_linear.pdf')