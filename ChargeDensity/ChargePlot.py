import pandas as pd, matplotlib.pyplot as plt, glob

f=glob.glob("*.dat")

df=pd.read_csv(f[1], sep="\s+")
dfr=pd.read_csv(f[0], sep="\s+")

x,p=df.iloc[:,0], df.iloc[:,1]
i,r=dfr.iloc[:,0], dfr.iloc[:,1]

Inp=input("Enter n for density plot \n or \n q for charge plot \n")

if(Inp=='n'):
	plt.semilogy(x,p, label=r"$n_{inv}$")
	plt.xlabel("distance (cm)", fontsize=18)
	plt.ylabel(r"$n_{inv} (cm^{-3})$", fontsize=18)
	plt.legend()
	plt.show()

if(Inp=='q'):
	plt.semilogy(i[:-2],r[:-2], label=r" $Q_{inv} (cm^{-2})$")
	plt.xlabel("Vgs (V)", fontsize=18)
	plt.ylabel("Residual", fontsize=18)
	plt.legend()
	plt.show()
