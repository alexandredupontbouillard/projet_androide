Le programme peut prendre en entrée des fichiers du type :


c mon commentaire1
c mon commentaire2
p 7 9 13
a 1 2 2
a 1 3 1
a 1 4 45
a 2 1 2
a 3 6 35
a 4 7 20
a 5 8 30
a 6 9 70
a 2 5 3
e 1 1 1
e 2 2 0.1
e 2 3 0.9
e 3 4 0.2
e 3 5 0.7
e 3 6 0.1
e 4 7 1
e 5 2 0.1
e 5 3 0.9
e 6 7 1
e 7 7 1
e 8 7 1
e 9 7 1

la lettre c en début de ligne indique un commentaire
la lettre p indique la déclaration du problème, suivit du nombre de sommet du graphe, puis du nombre d'action, 
puis du nombre de redirections stochastiques
les lignes a n1 n2 p représentent les actions accessibles en n1 amenant à la redirection stochastique n2 de coût p
les lignes e n1 n2 p représentent les redirections stochastiques de n1 vers n2 de probabilités p.