import matplotlib.pyplot as plt
import pickle
import numpy as np

a = open('dg.pickle', 'rb')
dGs = pickle.load(a)
a.close()

a = open('dg_error.pickle', 'rb')
dG_error = pickle.load(a)

#specify how many ps to discard for equilibration
discard = 1000
start = discard+1000

#Plot the cumulative dG estimate as a function of simulation time in ns
x_e = np.arange(start/1000,len(dGs)+start/1000,1)
plt.errorbar(x_e, dGs, yerr=dG_error, alpha=0.5, fmt = '-o', mfc='w', mew=2.5, lw=5, ms=10)
plt.xlabel('simulation time in ns',fontsize=20)
plt.ylabel('ΔG in kcal/mol',fontsize=20)
plt.grid(True)    
plt.savefig('dG.png', bbox_inches = "tight")
plt.show()

