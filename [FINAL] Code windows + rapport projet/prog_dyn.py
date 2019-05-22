# -*- coding: utf-8 -*-
"""
Created on Tue May 21 21:14:29 2019

@author: Adrien
"""

from lecture_fichier import *

def norme(u,v):
    res = 0
    for i in range(len(u)):
        res += (u[i] - v[i])**2
    return np.sqrt(res)

def valueIteration(J,P,cout,epsilon):
    nLignes = len(J)
    nCol = len(J[0])
    M = max(cout)
    V = np.ndarray.tolist(np.zeros(nCol))
    V1 = []
    for i in range(nCol):
        V1.append(M)
    while norme(V1,V)>=epsilon:
        V=list.copy(V1)         
        for i in range(nCol):
            mini = M
            for j in range(nLignes):
                if J[j][i] == 1:
                    temp = cout[j]
                    for k in range(nCol):
                        temp += P[j][k]*V[k]
                    if temp<mini:
                        mini = temp
            V1[i] = mini
        
    optiPoli = []
    for i in range(nCol):
        mini = max(V1)*nCol
        for j in range(nLignes):
            if J[j][i] == 1:
                temp = cout[j]
                for k in range(nCol):
                    temp += P[j][k]*V1[k]
                if temp<=mini:
                    actionOpti = j
                    mini = temp
        optiPoli.append(actionOpti+1)
    
    return optiPoli

def policyIteration(J,P,cout):
    nLignes = len(J)
    nCol = len(J[0])
    M = max(cout)
    # On initialise la politique par une politique "simple": on va toujours de l'avant quand on le peut
    # (on rappelle que par construction le noeud objectif est le noeud avec l'indice le plus élevé)
    poli = []
    for i in range(nCol):
        suivant = -1
        for j in range(nLignes):
            if J[j][i] == 1:
                if suivant == -1:
                    action = j
                for k in range(nCol):
                    if P[j][k]>0:
                        if k>suivant:
                            suivant = k
                            action = j  
        poli.append(action)
    poli0 = list.copy(poli)
    poli0[0] -= 1
    while poli != poli0:
        poli0 = list.copy(poli)
        matriceCoeff = []
        for i in range(nCol):
            matriceCoeff.append(np.ndarray.tolist(np.zeros(nCol)))
        vecteurConst = []
        for i in range(nCol):
            vecteurConst.append(cout[poli0[i]])
            matriceCoeff[i][i] += 1
            for j in range(nCol):
                matriceCoeff[i][j] -= P[poli0[i]][j]
        vecteurConst = np.array(vecteurConst)
        matriceCoeff = np.array(matriceCoeff)
        V = np.ndarray.tolist(np.linalg.solve(matriceCoeff,vecteurConst))
        
        for i in range(nCol):
            mini = max(V)*nCol
            actionOpti = []
            for j in range(nLignes):
                if J[j][i] == 1:
                    temp = cout[j]
                    for k in range(nCol):
                        temp += P[j][k]*V[k]
                    if temp == mini:
                        actionOpti.append(j)
                    if temp < mini:
                        actionOpti = [j]
                        mini = temp
            if poli0[i] in actionOpti:
                poli[i] = poli0[i]
            else:
                poli[i] = actionOpti[0]
        
    for i in range(nCol):
        poli[i] += 1
        
    return poli