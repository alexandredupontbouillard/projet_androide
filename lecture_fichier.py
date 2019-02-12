# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 13:32:00 2019

@author: alex0
"""
import numpy as np
import graphviz as pgv


#Pour la lecture de fichier on utilisera le format suivant :
#1) c chaine_de_caractère  ->  sera ignoré par la lecture de fichier, interprété comme un commentaire
#2) p N_noeud N_stoc N_actions N_redirection -> permettra de donner le nombre de noeud,d'actions, et de transitions stochastiques dans le graphe
#3) a n1 n2 poid  -> représentera une action pouvant être faites en n1 amenant au noeud stochastique n2 avec un certain poid
#4) e n1 n2 proba -> représentera une redirection stochastique faites en n1 amenant au choix n2 avec une certaine probabilité
# on fera attention à suivre l'orde c puis p puis a puis e
    

#prend une chaine de caractère et un séparateur et coupe la chaine de caractère à chaque séparateurs pour en former une liste
def ma_partition(sep,chaine):
    result=[]
    ch=chaine
    for i in range(4):
        c=ch.partition(sep)
        result.append(c[0])
        ch=c[2]
    return result

#renvoi une listes d'adjacences représentant les actions et une autre représentant les redirections possibles
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
        
    return actions,redirections,N_noeud,N_actions

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

actions,redirections,N_noeud,N_actions=lire_graphe('graphe_test.txt')

afficheGraphe(actions,redirections,N_noeud,N_actions)
    
    
    
    
    
    