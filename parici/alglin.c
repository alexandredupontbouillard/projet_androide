/* Marmote and MarmoteMDP are free softwares: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Marmote is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Marmote. If not, see <http://www.gnu.org/licenses/>.

Copyright 2018 EMmanuel Hyon, Alain Jean-Marie*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include "alglin.h"

/* Suite de fonctions de traitement de matrice */

/* affichages */
void  verifVd(double *Vecteur,int dim){
  int i;
  for(i=0;i<dim;i++){
    printf("%d  %6.10f\n",i,Vecteur[i]);
  }
  return;
}

void  verifVi(int *Vecteur,int dim){
  int i;
  for(i=0;i<dim;i++){
    printf("%d  %d\n",i,Vecteur[i]);
  }
  return;
}

void  verifMat(double **M,int dim){
int j,k;
  for(j=0;j<dim;j++){
    for(k=0;k<dim;k++){
      printf("%1.5f   ",M[j][k]);
    }
    printf("\n\n");
  }
  return;
}

void  verifMat2D(double **M,int dim1,int dim2){
int j,k;
  for(j=0;j<dim1;j++){
    for(k=0;k<dim2;k++){
      printf("%1.5f   ",M[j][k]);
    }
    printf("\n");
  }
  return;
}

void  verifMatP(double **M,int dim){
int j,k;
  for(j=0;j<dim;j++){
    for(k=0;k<dim;k++){
      printf("%1.10f   ",M[j][k]);
    }
    printf("\n\n");
  }
  return;
}

/* transposition d'une matrice carree */
void Transpose(double **Matrice,int dimension){
int i,j;
double temp;
for (i=0;i<dimension;i++){
  for(j=0;j<i;j++){
    temp=Matrice[i][j];
    Matrice[i][j]=Matrice[j][i];
    Matrice[j][i]=temp;
  }
}
return ;
}


/* Produit  matrice vecteur */
double *produitMatVect(double **A,double *B,int dim1){
  int i,l;
  double temp;
  double *Res;

  /*creation matrice */
  Res = (double *) calloc(dim1, sizeof(double));

  for(i=0;i<dim1;i++){
    temp=0;
    for(l=0;l<dim1;l++){
      temp=temp+(A[i][l]*B[l]);
    }
    Res[i]=temp;
  }
  return Res;
}

/* Produit  matrice matrice  carree*/
double **produitMatMat(double **A,double **B,int dim1){
  int i,j,l;
  double temp;
  double **Res;

  /*creation matrice */
  Res = (double **) calloc(dim1, sizeof(double*));

  for(i=0;i<dim1;i++)
    Res[i] = (double *) calloc(dim1, sizeof(double));

  for(i=0;i<dim1;i++){
    for(j=0;j<dim1;j++){
      temp=0;
      for(l=0;l<dim1;l++){
	temp=temp+(A[i][l]*B[l][j]);
      }
    Res[i][j]=temp;
    }
  }
  return Res;
}

double prodScalV(double *U,int Size){
  int  i;
  double res=0.0;
  for ( i=0; i<Size; i++ )
    res+=U[i];

  return res;
}


double Norm( double *U, double *V, int Size )
{
  int		i;
  double	Res;

  Res = 0.0;

  for ( i=0; i<Size; i++ )
    if ( fabs( U[i] - V[i] ) > Res ) {
      Res = fabs( U[i] - V[i] );
    }

  return Res;
}

/* calcul du span */
double Span( double *U, double *V, int Size )
{
  int       i;
  double    Min;
  double    Max;
  double    ecart;

  Min = U[0] - V[0];
  Max = U[0] - V[0];


  for ( i=1; i<Size; i++ ){
    ecart=U[i] - V[i];
    if ( ecart  > Max ){
      Max= ecart ;
    }
    if (  ecart  < Min ){
      Min= ecart;
    }
  }
  return (Max-Min);
}

/* calcul du span et recup des valeurs mex et min*/
double SpanRecup( double *U, double *V, int Size, double *maxi, double *mini)
{
  int       i;
  double    Min;
  double    Max;
  double    ecart;

  Min = U[0] - V[0];
  Max = U[0] - V[0];


  for ( i=1; i<Size; i++ ){
    ecart=U[i] - V[i];
    if ( ecart  > Max ){
      Max= ecart ;
    }
    if (  ecart  < Min ){
      Min= ecart;
    }
  }
  *maxi=Max;
  *mini=Min;
  return (Max-Min);
}




/******************************************************************************* */
/* resolution d'un systeme lineaire */
/* Pivot de gauss avec pivot le plus grand et pas de permutation effective (les permutations sont gardees dans un  vecteur */
double *ResolutionSysLin(double **Mat,double *Ve, int dim){
  int i,j,k,l;
  int inditemp; /* indice temporaire pour stocker valeur pivot */
  int inditab; /* indice temporaire de la ligne selectionnee pour pivot*/
  long double apiv; /* valeur pivot */
  int inpiv; /* indice pivot */
  long double temp; /* valeur temporaire pour les calculs des lignes*/

  int *permu; /* vecteur des permutations */
  long double *vtemp; /* vecteur pour la remontee */
  double *Vsol; /* solution */
  long double **M; /* matrice sur laquelle on travaille */
  long double *V;  /* vecteur sur lequel on travaille*/

  /* remplissage du vecteur de permutation*/
  permu=(int *) malloc( dim*sizeof(int) );
  for(i=0;i<dim;i++){
    permu[i]=i;
  }

   /*creation M */
  M= (long double **) calloc(dim, sizeof(long double*));
  for(i=0;i<dim;i++) {
    M[i] = (long double *) calloc(dim, sizeof(long double));
  }
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++)
      M[i][j]=(long double) Mat[i][j];
  }

  /* vecteur sol */
  V=(long double*) malloc(dim*sizeof(long double));
  for(i=0;i<dim;i++){
    V[i]=Ve[i];
  }


  /**** transformation de la matrice en matrice "diagonale" sans changement de ligne ****/
  for(i=0;i<dim-1;i++){
    inditab=i; /*ne sert pas */
    apiv=M[permu[i]][i];
    inpiv=permu[i];
    /*selection du pivot partiel sur les lignes*/
    for(j=i+1;j<dim;j++){
      if ((fabsl(M[permu[j]][i]))> fabsl(apiv)){
          inpiv=permu[j];
          inditab=j;
          apiv=M[permu[j]][i];
      }
    }
    // permutation si necessaire
    if (inpiv!= permu[i]){
      inditemp=permu[i];
      permu[i]=inpiv;
      permu[inditab]=inditemp;
    }
    // changement des termes (pivotage)
    for(j=i+1;j<dim;j++){
      temp=M[permu[j]][i]/apiv;
      /* Modif vecteur */
      // formule   V[tab[j,1],1]:= V[tab[j,1],1]-V[tab[i,1],1]*temp ;
      V[permu[j]]=(double) V[permu[j]]-(V[permu[i]]*temp);
      /* modif matrice */
      for(k=i;k<dim;k++){
	M[permu[j]][k]=M[permu[j]][k]-(M[permu[i]][k]*temp);
      }
    }
  }

  /******** forme triangulaire sup�rieure ******/

  /*verif matrice et vecteur permu dans resl sys lin*/
  /*printf("\n\n");
  printf("Permu\n");
  verifVi(permu,dim);
  printf("\n");
  printf("Vecteur\n");
  verifVD(V,dim);
  printf("\n");
  printf("Matrice\n");
  for(j=0;j<dim;j++){
    for(k=0;k<dim;k++){
      printf("%1.5Lf  ", M[j][k]);
    }
    printf("\n");
  }
  */


  /***** remontee *****/
  /* Creation du vecteur temporaire*/
  vtemp=(long double *) malloc(dim*sizeof(long double));
  for(i=0;i<dim;i++){
    vtemp[i]=0;
  }

  i= dim-1;
  vtemp[permu[i]]=V[permu[i]]/M[permu[i]][i];
  for(i=2;i<=dim;i++){
    temp=0;
    k=dim-i;
    for(j=1;j<=i;j++){
       l= dim-j;
       temp= temp+M[permu[k]][l]*vtemp[permu[l]];
    }
    temp=V[permu[k]]-temp;
    vtemp[permu[k]]=temp/M[permu[k]][k];
  }



  /*remise des valeurs dans l'ordre et creation du vecteur que l'on retournera */
  Vsol=(double *) malloc(dim*sizeof(double));
  for(i=0;i<dim;i++){
    Vsol[i]=(double) vtemp[permu[i]];
  }


  /** Nettoyage des malloc */
  for(i=0;i<dim;i++)
    free(M[i]);
  free(M);
  free(V);
  free(permu);
  free(vtemp);
 return Vsol;

}
/* Fin **************************************************************************************** */








/* ******************************************************************************************** */
/* **************** Inversion matrice *****************************/

/* Inversion par pivot  de gauss avec pivot le plus grand et pas de permutation effective (les permutations sont gardees dans un  vecteur */
void Inversion(double **Mat,double **Inverse, int dim){
  int i,j,k,t;
  int inditemp; /* indice temporaire pour stocker valeur pivot */
  int inditab; /* indice temporaire de la ligne selectionnee pour pivot*/
  long double apiv; /* valeur pivot */
  int inpiv; /* indice pivot */
  long double temp; /* valeur temporaire pour les calculs des lignes*/

  int *permu; /* vecteur des permutations */
  long double **M; /* matrice sur laquelle on travaille */
  long double **V;  /* inverse de la matrice sur laquel on travaille*/

  /* remplissage du vecteur de permutation*/
  permu=(int *) malloc( dim*sizeof(int) );
  for(i=0;i<dim;i++){
    permu[i]=i;
  }

   /*creation M et V */
  M= (long double **) calloc(dim, sizeof(long double*));
  V = (long double **) calloc(dim, sizeof(long double*));
  for(i=0;i<dim;i++) {
    M[i] = (long double *) calloc(dim, sizeof(long double));
    V[i] = (long double *) calloc(dim, sizeof(long double));
  }
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      M[i][j]=(long double) Mat[i][j];
      V[i][j]= 0.0;
    }
    V[i][i]=1.0;
  }



  /**** transformation de la matrice en matrice "diagonale" sans changement de ligne ****/
  for(i=0;i<dim-1;i++){
    inditab=i;
    apiv=M[permu[i]][i];
    inpiv=permu[i];
    //selection du pivot partiel sur les lignes
    for(j=i+1;j<dim;j++){
      if ((fabsl(M[permu[j]][i]))> fabsl(apiv)){
          inpiv=permu[j];
          inditab=j;
          apiv=M[permu[j]][i];
      }
    }
    // permutation si necessaire
    if (inpiv!= permu[i]){
      inditemp=permu[i];
      permu[i]=inpiv;
      permu[inditab]=inditemp;
    }
    // changement des termes (pivotage)
    for(j=i+1;j<dim;j++){
      temp=M[permu[j]][i];   /* ou encore  M[permu[j]][i]/apiv mais on divise plus parapiv ensuite*/
      /* Modif Matrice Identite */
      // formule   V[tab[j,1],k]:= V[tab[j,1],k]-V[tab[i,1],k]*temp ;
      for(k=0;k<dim;k++){
	// k est la colonne
	V[permu[j]][k]=V[permu[j]][k]-(V[permu[i]][k]*temp)/apiv;
      }
      for(k=i;k<dim;k++){
	/* modif matrice */
	M[permu[j]][k]=M[permu[j]][k]-(M[permu[i]][k]*temp)/apiv;
      }
    }
  }

  /*on a forme triangulaire sup�rieure */
  /*verif matrice triangulaire sup et vecteur permu
  printf("\n\n");
  printf("Permu\n");
  verifVi(permu,dim);
  printf("\n");
  printf("Matrice Inverse \n");
  verifMat(V,dim);
  printf("\n");
  printf("Matrice\n");
  for(j=0;j<dim;j++){
    for(k=0;k<dim;k++){
      printf("%1.4Lf  ",M[j][k]);
    }
    printf("\n");
  }
  */

 /***** remontee *****/
  for(i=1;i<=dim;i++){
    k=dim-i;
    // divise la ligne par pivot
     for(t=0;t<dim;t++){
       V[permu[k]][t]=V[permu[k]][t]/M[permu[k]][k];
     }
     // je met le pivot � 1
     M[permu[k]][k]=1.0;
     // je remplace chaque ligne sur la colonne
      for(j=0;j<k;j++){
         for(t=0;t<dim;t++){
	   V[permu[j]][t]= V[permu[j]][t]-V[permu[k]][t]*M[permu[j]][k];
         }
         // je remplace M par 0 c est juste pour debug on en a pas besoin apres
         M[permu[j]][k]=M[permu[j]][k]-M[permu[j]][k]*M[permu[k]][k];
      }
  }




  /*remise des valeurs dans l'ordre et creation du vecteur que l'on retournera */
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++)
      Inverse[i][j]=(double) V[permu[i]][j];
  }


  /** Nettoyage des malloc */
  for(i=0;i<dim;i++){
    free(M[i]);
    free(V[i]);
  }
  free(M);
  free(V);
  free(permu);
 return;
  }




/* ********************************************************************************************************************** */
/* ***** calcul d'une probabilit� invariante en fonction generateur pour processus continu ***** */
/* les equations redondantes vont  etre supprimees dedans*/
double *CalculPI(double **Generator,int dim){
  int i,j;
  double *pi; /* proba invariante */
  double *p;  /* proba initiale */
  double **M; /* matrice qui stocke le generateur */

   /*creation M */
  M= (double **) calloc(dim, sizeof(double*));
  for(i=0;i<dim;i++) {
    M[i] = (double *) calloc(dim, sizeof(double));
  }
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++)
      M[i][j]= Generator[i][j];
  }

  // Suppression des redondances et normalisation
  for(i=0;i<dim;i++){
     M[i][dim-1]= 1;
  }
  /* vecteur sol */
  p=(double *) malloc(dim*sizeof(double));
  for(i=0;i<dim;i++){
    p[i]=0.0;
  }
  p[dim-1]=1;

  /* transposition */
  Transpose(M,dim);
  /*printf("Matrice Transposee\n");*/
  // resolution
  pi=ResolutionSysLin(M,p, dim);

  for(i=0;i<dim;i++)
    free(M[i]);
  free(M);
  free(p);
  return pi;
}



/* ********************************************************************************************************************** */
/* ***** calcul d'une probabilit� invariante en fonction generateur pour processus continu ***** */
/* les equations redondantes ont a etre supprimees */
double *CalculPIChaine(double **Noyau,int dim){
  int i,j;
  double *pi; /* proba invariante */
  double *p;  /* proba initiale */
  double **M; /* matrice qui stocke le generateur */

   /*creation M */
  M= (double **) calloc(dim, sizeof(double*));
  for(i=0;i<dim;i++) {
    M[i] = (double *) calloc(dim, sizeof(double));
  }
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++)
      M[i][j]= Noyau[i][j];
  }


  // integration de soustraction identit�
  for(i=0;i<dim;i++){
     M[i][i]= M[i][i]-1.0;
  }

  // Suppression des redondances
  for(i=0;i<dim;i++){
     M[i][dim-1]= 1.0;
  }


  /* vecteur sol */
  p=(double*) malloc(dim*sizeof(double));
  for(i=0;i<dim;i++){
    p[i]=0.0;
  }
  p[dim-1]=1.0;

  /* pour debug
  printf("DANS alglin Verif valeur Mat\n");
  verifMat(M,dim);
  printf("\n \n Verif vecteur \n");
  verifVd(p,dim);
  fin pour debug */

  /* transposition */
  Transpose(M,dim);
  /*printf("Matrice Transposee\n");*/

  // resolution
  pi=ResolutionSysLin(M,p, dim);

  for(i=0;i<dim;i++)
    free(M[i]);
  free(M);
  free(p);
  return pi;
}


/* Pour debug
void main(){
  int i,di;
  di=3;
  double **Ma;
  double **Inv;

  Ma = (double **) calloc(di, sizeof(double*));
  Inv = (double **) calloc(di, sizeof(double*));
  for(i=0;i<di;i++) {
    Ma[i] = (double *) calloc(di, sizeof(double));
    Inv[i] = (double *) calloc(di, sizeof(double));
  }



   Ma[0][0]=1;  Ma[0][1]=3;  Ma[0][2]=3;
  Ma[1][0]=2;  Ma[1][1]=2;  Ma[1][2]=0;
  Ma[2][0]=3; Ma[2][1]=3; Ma[2][2]=6;

  printf("Dans main\n");
  verifMat(Ma,3);
  printf("\n");
  Inversion(Ma,Inv,3);
   printf("\ninverse\n");
   verifMat(Inv,3);


}
*/
