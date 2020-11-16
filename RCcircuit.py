#! python3.7
"""
RCcircuit.py
version 1
Charge d'un circuit RC, tension aux bornes du condensateur (Uc(t)) et aux bornes de la resistance (Ur(t))
Représentation graphique pour un echelon de tension

Auteur : Kenneth Maussang, kenneth.maussang@umontpellier.fr
Licence : GNU General Public License
"""
# importation des modules utiles
import numpy as np
import scipy.integrate as itg
import matplotlib.pyplot as plt
from matplotlib import rc # permet d'ecrire des formules en LaTeX dans les labels

# parametres du circuit RC etudie
R=1000 # R value in Ohm
C=1 # 1 value in microFarad

# Calcul de la constante de temps pertinente tau=RC
tau=R*C*1e-6
# Definition d'un vecteur temporel t pour l'abscisse
Npts = 250 # nb de points pour les calculs numeriques
t=np.linspace(-tau,10*tau,Npts)
# Conditions initiales sur la tension Uc(t) aux bornes du condensateur
U0=0
# Amplitude maximale de la tension excitrice
Emax=10

# Definition d'une tension d'excitation
def E(t):
    if t<0 :
        return 0.0
    else :
        return Emax

# definition de du/dt a partir de l'equation differentielle
def dUc(Uc, t):
    return (E(t)-Uc)/tau

# integration de l'equation differentielle sur Uc(t), par rapport a la variable t
Uc = itg.odeint(dUc, [U0], t)

# representation graphique
plt.close('all') # fermeture de toutes les fenetres de trace de courbes matplotlib
fig=plt.figure(1, figsize=(7, 5)) #c reation de la figure 1 de taille definie
plt.plot([ -1000*tau, 10000*tau], [ Emax, Emax ], "y--", linewidth=2, label=None) # construction des asymptotes (non labelises dans legend)
plt.plot([ 1000*tau, 1000*tau ], [-Emax/10, 11.0*Emax ], "y--", linewidth=2, label=None)
plt.plot([ 3*1000*tau, 3*1000*tau ], [-Emax/10, 11.0*Emax ], "y--", linewidth=2, label=None)
plt.plot(1000*t,[E(x) for x in t],'-',color='green', linewidth=4, label=r'$E(t)$ (V)') # trace de E(t)
plt.plot(1000*t,Uc,'-',color='blue', linewidth=4, label=r'$U_c(t)$ (V)') # trace de Uc(t)
plt.plot([ 0, 1000*tau ], [ 0.0, Emax], "r-.", linewidth=2, label=None) # tangente a l'origine
plt.plot([0],[0],'o',color='red',markersize=6, label=None)
plt.plot([1000*tau],[Emax],'o',color='red',markersize=6, label=None)
plt.plot([1000*tau],[0.63*Emax],'o',color='black',markersize=6, label=None) # t=tau
plt.plot([3000*tau],[0.95*Emax],'o',color='black',markersize=6, label=None) # t=3*tau
plt.xlabel(r'$t$ (ms)')
plt.ylabel(r'$E(t)$ (V) et $U_c(t)$ (V) ')
plt.title(r'Circuit RC, R='+str(R/1000)+r'k$\Omega$, C='+str(C)+r'$\mu$F, $E_{max}$='+str(Emax)+r'V')
plt.legend(loc='best') # legende, localisation 'best'
plt.grid() # ajout d'une grille
plt.annotate(r"$t=\tau$", (1000*tau, 0), (1.05*1000*tau, -0.05*Emax))
plt.annotate(r"$t=3\tau$", (3000*tau, 0), (3.05*1000*tau, -0.05*Emax))
plt.annotate("63%", (1000*tau, 0.63*Emax), (1500*tau, 0.52*Emax), arrowprops={"arrowstyle":"->"})
plt.annotate("95%", (3000*tau, 0.95*Emax), (3500*tau, 0.84*Emax), arrowprops={"arrowstyle":"->"})
plt.xlim(-1000*tau,5000*tau) # axe des abscisses sur [-tau,5*tau]
plt.ylim(-0.1*Emax,1.1*Emax) # axe des ordonnees sur [-0.1*Emax, 1.1*Emax]
plt.tight_layout()
plt.savefig('circuit_RC.png',dpi=200) # export de la figure en .png, resolution = 200dpi
plt.show()