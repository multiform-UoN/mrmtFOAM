from plt_settings import *
import scipy.io as spio
plt.style.use('seaborn')


# Import matlab stuff
mat = spio.loadmat('resultsChebfun.mat', squeeze_me=True)
Cheb = mat['M']

# Import data from breakthrough
codeRes = np.loadtxt(open("breakthrough.dat", "rb"), delimiter=" ", skiprows=0)

plt.figure(figsize=(250 /25.4, 200 / 25.4))
plt.plot(Cheb[:,0],Cheb[:,1],'-k',label='Chebfun')
plt.plot(codeRes[:,0],codeRes[:,1],'bo',label='mrmtScalarTransportFoam')
plt.xlabel('$t (s)$',fontsize=45)
plt.ylabel('$BT$',fontsize=45)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.ylim(0, 0.25)
plt.xlim(0, 10)
#plt.show()
plt.legend(loc='upper right', frameon=False,fontsize=30)

#    plt.rc('grid', linestyle="-", color='black')
plt.grid(True)
#    plt.tight_layout()
plt.savefig('verification.pdf')
