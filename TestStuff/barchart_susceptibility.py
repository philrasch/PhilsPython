
import numpy as np
import matplotlib.pyplot as plt


from matplotlib.backends.backend_pdf import PdfPages

# bar graphs
x = np.arange(5)
#y1, y2 = np.random.randint(1, 25, size=(2, 5))
y1 = [0.05,0.46,0.45,0.10,0.03]
width = 0.35
fig, ax = plt.subplots()
ax.bar(x, y1, width)
#ax.bar(x + width, y2, width, color=plt.rcParams['axes.color_cycle'][2])
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)
ax.set_xticks(x + width/2)
ax.set_xticklabels(['OBS', 'Mod$\mathrm{_{clim}}$', 'Mod$\mathrm{_{orb}}$', 'Mod$\mathrm{_{orb,cld}}$', 'Mod$\mathrm{_{orb,cld,aer}}$'],rotation='45')
#plt.setp(ax.yaxis.get_ticklines(), 'markersize', 25)
#plt.setp(ax.yaxis.get_ticklines(), 'markeredgewidth', 3)
ax.set_title('Cloud Susceptibility to Aerosols',fontsize=20)
#plt.xlabel('xlabel',fontsize=15)
plt.ylabel('$\mathrm{\partial R_{e} / \partial AOD}$',fontsize=20)
plt.subplots_adjust(bottom=0.25)
#plt.axes().set_aspect('equal')
#plt.tight_layout()
plt.savefig('test.pdf')
plt.show()
