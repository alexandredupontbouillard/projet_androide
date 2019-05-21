# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 16:18:18 2019

@author: Adrien
"""

import numpy as np
import random
import turtle

# Paramètres: n² = taille de la grille
#             npuits = nombre de noeuds à éviter, au maximum 3n
#             i = numéro de l'instance (pour pouvoir en générer plusieurs avec une boucle)
#             p1 = proba de ne pas "tomber" (quand on choisit une action avec un risque de descendre, p1 petit)
#             p2 = proba de "grimper" (quand on choisit une action de type "monter" ou gauche/droite, p2~0.5)
#             relief = mesure la linéarité de la frontière entre versant nord et sud.
#                      proche de 0 si peu de relief, proche de 100 si beaucoup de relief
def creeInstance(n,npuits,i,p1,p2,relief):
    instance = open("montagne" + str(i) + ".txt","w+")
    instance.write("c "+"Instance numero "+str(i)+" du type montagne")
    # On "dessine" le versant de la montagne de manière aléatoire
    limVersant = [int(n/2)]
    for i in range(n):
        courant = limVersant[-1]
        limVersant.append(round(courant+relief/100*1.5*np.random.normal(n/2-courant)))
    # On rajoute des points à éviter de manière aléatoire
    if npuits <= n:
        posPuits = random.sample(range(n),npuits)
        for i in range(npuits):
            courant = posPuits[i]
            posPuits[i] = n*random.randint(limVersant[courant],limVersant[courant]+1)+courant
    elif npuits <= 2*n:
        posPuits1 = random.sample(range(n),npuits//2)
        posPuits2 = random.sample(range(n),npuits//2+npuits%2)
        for i in range(npuits//2):
            courant = posPuits1[i]
            posPuits1[i] = n*random.randint(limVersant[courant]+2,limVersant[courant]+3)+courant
        for i in range(npuits//2+npuits%2):
            courant = posPuits2[i]
            posPuits2[i] = n*random.randint(limVersant[courant]-2,limVersant[courant]-1)+courant
        posPuits = posPuits1+posPuits2
    elif npuits <= 3*n:
        posPuits1 = random.sample(range(n),npuits//3)
        posPuits2 = random.sample(range(n),npuits//3)
        posPuits3 = random.sample(range(n),npuits//3+npuits%3)
        for i in range(npuits//3):
            courant = posPuits1[i]
            posPuits1[i] = n*random.randint(limVersant[courant]+4,limVersant[courant]+5)+courant
        for i in range(npuits//3):
            courant = posPuits2[i]
            posPuits2[i] = n*random.randint(limVersant[courant]-1,limVersant[courant])+courant
        for i in range(npuits//3+npuits%3):
            courant = posPuits3[i]
            posPuits3[i] = n*random.randint(limVersant[courant]-6,limVersant[courant]-5)+courant
        posPuits = posPuits1+posPuits2+posPuits3
        
    listeA = []
    listeE = []
    nAction = 1
    nNoeud = 1
    # On crée les listes associant noeud -> actions possibles (listeA)
    #                 et action -> noeuds d'arrivée possibles (listeE)
    for ligne in range(n):
        for colonne in range(n):
            # On rajoute les conditions au bord (si on essaye de sortir, on reste à sa place)
            # La condition poids est rajoutée pour les noeuds à éviter
            bordGauche,bordDroit,bordBas,bordHaut,poids = 1,1,n,n,1
            if n*ligne+colonne in posPuits:
                poids = n**10
            if colonne == 0:
                bordGauche = 0
            if colonne == n-1:
                if ligne == n-1: # case objectif
                    break
                else:
                    bordDroit = 0
            if ligne == 0:
                bordBas = 0
            if ligne == n-1:
                bordHaut = 0
            # On remplit les listes
            listeA.append([nNoeud,nAction,poids])
            listeA.append([nNoeud,nAction+1,poids])
            listeA.append([nNoeud,nAction+2,poids])
            listeA.append([nNoeud,nAction+3,poids])
            # Condition pour savoir sur quel versant on est (plus de chance de descendre que de monter)
            if limVersant[colonne]<=ligne:
                if bordDroit == bordHaut == 0:
                    listeE.append([nAction,nNoeud,1])
                else:
                    listeE.append([nAction,nNoeud+bordDroit,p2])
                    listeE.append([nAction,nNoeud+bordHaut,1-p2])
                if bordDroit == bordBas == 0:
                    listeE.append([nAction+1,nNoeud,1])
                else:
                    listeE.append([nAction+1,nNoeud+bordDroit,p1])
                    listeE.append([nAction+1,nNoeud-bordBas,1-p1])
                if bordGauche == bordHaut == 0:
                    listeE.append([nAction+2,nNoeud,1])
                else:
                    listeE.append([nAction+2,nNoeud-bordGauche,p2])
                    listeE.append([nAction+2,nNoeud+bordHaut,1-p2])
                if bordGauche == bordBas == 0:
                    listeE.append([nAction+3,nNoeud,1])
                else:
                    listeE.append([nAction+3,nNoeud-bordGauche,p1])
                    listeE.append([nAction+3,nNoeud-bordBas,1-p1])
            else:
                if bordDroit == bordHaut == 0:
                    listeE.append([nAction,nNoeud,1])
                else:
                    listeE.append([nAction,nNoeud+bordDroit,p1])
                    listeE.append([nAction,nNoeud+bordHaut,1-p1])
                if bordDroit == bordBas == 0:
                    listeE.append([nAction+1,nNoeud,1])
                else:
                    listeE.append([nAction+1,nNoeud+bordDroit,p2])
                    listeE.append([nAction+1,nNoeud-bordBas,1-p2])
                if bordGauche == bordHaut == 0:
                    listeE.append([nAction+2,nNoeud,1])
                else:
                    listeE.append([nAction+2,nNoeud-bordGauche,p1])
                    listeE.append([nAction+2,nNoeud+bordHaut,1-p1])
                if bordGauche == bordBas == 0:
                    listeE.append([nAction+3,nNoeud,1])
                else:
                    listeE.append([nAction+3,nNoeud-bordGauche,p2])
                    listeE.append([nAction+3,nNoeud-bordBas,1-p2])
            nAction += 4
            nNoeud += 1
                        
    # On n'a plus qu'à remplir le fichier .txt de manière à respecter le format d'instance choisi
    instance.write("\n"+"p "+str(nNoeud)+" "+str(len(listeA))+" "+str(len(listeE)))
    for i in range(len(listeA)):
        instance.write("\n"+"a "+str(listeA[i][0])+" "+str(listeA[i][1])+" "+str(listeA[i][2]))
    for i in range(len(listeE)):
        instance.write("\n"+"e "+str(listeE[i][0])+" "+str(listeE[i][1])+" "+str(listeE[i][2]))
    instance.close()
    
    return posPuits, limVersant
    
def dessinInstanceM(size,N_noeud,posPuits,limVersant,initial,policy,redirections):
    longueur = int(np.sqrt(N_noeud))
    width = 5
    screen = turtle.Screen()
    screen.bgcolor("black")
    turtle.pensize(width)
    turtle.speed(0)
    
    turtle.penup()
    turtle.left(180)
    turtle.forward(0.5*size*longueur)
    turtle.left(90)
    turtle.forward(0.5*size*longueur)
    turtle.left(90)
    turtle.pendown()
    
    for colonne in range(longueur):
        for ligne in range(longueur):
            if ligne<limVersant[colonne]:
                turtle.color("white")
            else:
                turtle.color("#8dade0")
            for i in range(4):
                turtle.forward(size)
                turtle.left(90)
            turtle.penup()
            turtle.left(90)
            turtle.forward(size)
            turtle.right(90)
            turtle.pendown()
        turtle.penup()
        turtle.forward(size)
        turtle.right(90)
        turtle.forward(size*longueur)
        turtle.left(90)
        turtle.pendown()
        
    turtle.color("#bc104d")
    for elem in posPuits:
        turtle.penup()
        turtle.home()
        turtle.left(180)
        turtle.forward(0.5*size*longueur)
        turtle.left(90)
        turtle.forward(0.5*size*longueur)
        turtle.left(90)
        quotient, reste = elem//longueur, elem%longueur
        turtle.forward(reste*size)
        turtle.left(90)
        turtle.forward(quotient*size)
        turtle.right(90)
        turtle.pendown()
        for i in range(4):
            turtle.forward(size)
            turtle.left(90)
        for i in range(int(size/width)):
            turtle.forward(size-width)
            turtle.left(90)
            turtle.forward(width/2)
            turtle.left(90)
            turtle.forward(size-width)
            turtle.right(90)
            turtle.forward(width/2)
            turtle.right(90)
    
    for j in range(5):    
        turtle.penup()
        turtle.home()
        turtle.left(180)
        turtle.forward(0.5*size*longueur-size/2)
        turtle.left(90)
        turtle.forward(0.5*size*longueur-size/2)
        turtle.left(90)
        courant = initial-1
        quotient,reste = courant//longueur, courant%longueur
        turtle.forward(reste*size)
        turtle.left(90)
        turtle.forward(quotient*size)
        turtle.right(90)
        turtle.pendown()
        turtle.pensize(10)
        if j == 0:
            turtle.color("#69e07f")
        elif j == 1:
            turtle.color("#3773d3")
        elif j == 2:
            turtle.color("#e51d49")
        elif j == 3:
            turtle.color("#efb837")
        elif j == 4:
            turtle.color("#e527d8")
                         
        turtle.speed(5)
        recap=[]
        while 1:
            old = courant+1
            courant = [elem for elem in redirections if elem[0] == policy[courant]]
            print(courant)
            ber = np.random.binomial(1,courant[0][2])
            courant = ber*courant[0][1] + (1-ber)*courant[1][1]
            quotient,reste = (courant-old)//longueur, (courant-old)%longueur
            recap.append([courant,quotient,reste])
            courant -= 1
            if courant == N_noeud-1:
                if old == N_noeud-longueur:
                    turtle.left(90)
                    turtle.forward(size)
                    break
                else:
                    turtle.forward(size)
                    break
            else:
                if reste == 1:
                    turtle.forward(size)
                elif reste == longueur-1:
                    turtle.left(180)
                    turtle.forward(size)
                    turtle.left(180)
                elif quotient == 1:
                    turtle.left(90)
                    turtle.forward(size)
                    turtle.right(90)
                elif quotient == -1:
                    turtle.right(90)
                    turtle.forward(size)
                    turtle.left(90)

    screen.exitonclick()
    
    return recap