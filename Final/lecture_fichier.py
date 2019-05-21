# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 13:32:00 2019
@author: alex0
"""
import numpy as np
import graphviz as pgv


#Pour la lecture de fichier on utilisera le format suivant :
#1) c chaine_de_caractère  ->  sera ignoré par la lecture de fichier, interprété comme un commentaire
#2) p N_noeud N_stoc N_actions N_redirection -> permettra de donner le nombre de noeud,d'actions, et de noeuds chance dans le graphe
#3) a n1 n2 poid  -> représentera une action pouvant être faite en n1 amenant au noeud stochastique n2 avec un certain poids
#4) e n1 n2 proba -> représentera un tirage chance fait en n1 amenant au choix n2 avec une certaine probabilité
# on fera attention à suivre l'ordre c puis p puis a puis e
    

#prend une chaine de caractères et un séparateur et coupe la chaine de caractère à chaque séparateur pour en former une liste
def ma_partition(sep,chaine):
    result=[]
    ch=chaine
    for i in range(4):
        c=ch.partition(sep)
        result.append(c[0])
        ch=c[2]
    return result

#renvoit une listes d'adjacences représentant les actions et une autre représentant les redirections possibles
def lire_graphe(filename):
    data = np.loadtxt ( filename, delimiter='\n', dtype=np.str )
    compteur=0
    ligne = data[0]
    while(ligne[0]=='c'):
        compteur= compteur+1
        ligne = data[compteur]
    actions=[]
    redirections=[]
    l=ma_partition(' ',ligne)
    N_noeud = int(l[1])
    N_actions = int(l[2])
    N_stoc = int(l[3])
    compteur=compteur+1
    for i in range(0,N_actions):
        ligne=data[compteur]
        l=ma_partition(' ',ligne)
        actions.append([int(l[1]),int(l[2]),int(l[3])])
        compteur = compteur+1
    for i in range(0,N_stoc):
        ligne=data[compteur]
        l=ma_partition(' ',ligne)
        redirections.append([int(l[1]),int(l[2]),float(l[3])])
        compteur = compteur+1
    N_util = []
    N_but = 0
    for i in actions:
        N_util.append(i[0])
    N_util = list(set(N_util))
    for i in range(N_util[-1],N_noeud):
        N_but+=1
    return actions,redirections,N_noeud,N_actions,N_but

def afficheGraphe(l_actions,l_redirections, N_noeud, N_actions):
    graph = pgv.Digraph()
    for i in range(N_noeud):
        graph.node('n'+str(i+1))
    for i in range(N_actions):
        graph.node('a'+str(i+1),shape='Msquare')
    for i in l_actions:
        graph.edge('n'+str(i[0]),'a'+str(i[1]),label = str(i[2]))
    for i in l_redirections:
        graph.edge('a'+str(i[0]),'n'+str(i[1]),label = str(i[2]))
    
    graph.view()

#renvoit une liste de matrices, une matrice pour chaque action
#et un vecteur de coût pour chaques actions
def modelMDP(actions,redirections,N_noeud,N_actions):
    liste_actions = np.zeros((N_actions,N_noeud,N_noeud))
    cout = np.zeros(N_actions)
    for i in range(N_actions):
        cout[actions[i][1]-1] = actions[i][2]
        for j in range(len(redirections)):
            if(redirections[j][0]-1==i):
            
                liste_actions[i][actions[i][0]-1][redirections[j][1]-1]=redirections[j][2]
    return liste_actions,cout

#retourne une matrice a_ij = 1 si l'action i est accessible depuis le sommet j
#et une seconde matrice m_sa = p(s|a)     
def modelSTAUFFER(actions,redirections,N_noeud,N_actions,N_but):
    matrice_action = np.zeros((N_actions,N_noeud))
    matrice_p_s_a = np.zeros((N_actions,N_noeud))
    for i in range(N_actions):
        matrice_action[actions[i][1]-1][actions[i][0]-1]=1
    for j in range(len(redirections)):
        matrice_p_s_a[redirections[j][0]-1][redirections[j][1]-1] = redirections[j][2]
    matrice_action = np.ndarray.tolist(matrice_action)
    matrice_p_s_a = np.ndarray.tolist(matrice_p_s_a)
    matrice_arcs = np.zeros((N_noeud,N_noeud))
    for i in range(N_actions):
        for j in range(N_noeud):
            if matrice_action[i][j]>0:
                for k in range(N_noeud):
                    if matrice_p_s_a[i][k]>0:
                        matrice_arcs[j][k] = 1
    for i in range(N_but):
        for j in range(N_actions):
            del(matrice_action[j][-1])
            del(matrice_p_s_a[j][-1])
    return matrice_action,matrice_p_s_a,matrice_arcs

#x= sparseMatrix(4)
#l_x = sparseMatrixVector(4)
#l_x[0]= x
#l_x[1] = x
#l_x[2] = x
#l_x[3] = x
#y = sparseMatrix(4)
#m1 = marmoteInterval(0,4-1)
#m2 = marmoteInterval(0,4-1)
#mdp1 = totalRewardMDP("min", m1, m2, l_x, y)
