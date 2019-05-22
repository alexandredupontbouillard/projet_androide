# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from lecture_fichier import *
import gurobipy as grb


def sSSPdual(J,P,cout):
    nLignes = len(J)
    nCol = len(J[0])
    y = []
    try:
        
        # On crée le modèle
        m = grb.Model("PL dual")
        
        # On ajoute les variables
        for i in range(nCol):
            y.append(m.addVar(vtype=grb.GRB.CONTINUOUS,name=str(["y",i+1])))
        
        # On écrit la fonction objectif
        obj = grb.LinExpr()
        for i in range(nCol):
            obj += y[i]*1/nCol
        m.setObjective(obj, grb.GRB.MAXIMIZE)
        
        # On ajoute les contraintes
        for i in range(nLignes):
            contr = grb.LinExpr()
            for j in range(nCol):
                contr += (J[i][j]-P[i][j])*y[j]
            m.addConstr(contr <= cout[i])
        
        # On résout
        m.optimize()
        
        # On déduit la stratégie optimale
        strat = []
        M = max(cout)
        for i in range(nCol):
            opti = M
            action = -1
            for j in range(nLignes):
                minPot = 0
                if J[j][i] == 1:
                    if action<0:
                        action = j
                    minPot += cout[j]
                    for k in range(nCol):
                        minPot += P[j][k]*m.getVars()[k].x
                    if minPot < opti:
                        opti = minPot
                        action = j
            strat.append(action+1)
        
        return strat
    
    except grb.GurobiError:
        print('Error reported')