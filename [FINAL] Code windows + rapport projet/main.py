# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:00:33 2019

@author: Adrien
"""

from prog_lin import *
from prog_dyn import *
from instance_montagne import *
import time

# On génère une liste d'instances, et on sauvegarde les positions des puits et du versant de chaque instance
# (pour pouvoir les tracer plus tard si on le souhaite)
posPuits = []
limVersant = []
for i in range(10): # On prend 10 instances ici, mais l'utilisateur peut choisir la valeur qu'il souhaite
    # Format: voir instance_montagne.py, on utilise random pour plus d'aléatoire dans la génération d'instance
    posPuits_i, limVersant_i = creeInstance(15,16,i,0.1,0.3,30)
    posPuits.append(posPuits_i)
    limVersant.append(limVersant_i)

# On détermine les politiques optimales pour chacune des instances créées, et on mesure le temps que ça prend
L = []
optiPoli = []
optiPoli2 = []
optiPoli3 = []
epsilon = 0.1
for k in range(10):
    actions,redirections,N_noeud,N_actions,N_but=lire_graphe('montagne' + str(k) + '.txt')
    liste_actions, cout = modelMDP(actions,redirections,N_noeud,N_actions)
    m1,m2 = modelSTAUFFER(actions,redirections,N_noeud,N_actions,N_but)
    start = time.time()
    optiPoli.append(policyIteration(m1,m2,cout))
    middle = time.time()
    L.append(middle-start)
    optiPoli2.append(sSSPdual(m1,m2,cout))
    end = time.time()
    L.append(end-middle)
#    optiPoli3.append(valueIteration(m1,m2,cout,epsilon))
#    L.append(time.time()-end)
print(L)

# On dessine plusieurs instanciations d'une politique (ici, la politique obtenue sur la 8ème instance générée
# précédemment, par méthode de policyIteration, en partant du noeud 1 à chaque fois)

k = 5
dessinInstanceM(N_noeud,posPuits[k-1],limVersant[k-1],1,optiPoli[k-1],redirections)
turtle.reset()