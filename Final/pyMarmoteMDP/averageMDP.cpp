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

/**
 * @brief Class averageMDP source: class for implementing an averageMDP class
 * @author Hyon, lip6.
 * @version 0.2
 * @date jan 2019
 *
 */

#include <values.h>
#include <iostream>
#include "header/averageMDP.h"
#include "alglin.h"


using namespace std;

// constructeurs

averageMDP::averageMDP(string r, marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, sparseMatrix* rews)
: genericMDP("discrete",r,states,actions,trans,rews), rho(0.0), type_c("average"), indexEtatSpecifique(0)
{
}


averageMDP::averageMDP(string r, marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, vector<sparseMatrix*> rews)
: genericMDP(string("discrete"),r,states,actions,trans,rews), rho(0), type_c("average"), indexEtatSpecifique(0)
{
}

// other utilities function
void averageMDP::writeMDP()
{
    genericMDP::writeMDP();
    printf("#############################################\n");
    printf("write Average MDP\n");
    printf("MDP Criteria : %s \n" , type_c.c_str());
}

//A function to solve average infinite  MDP using value iteration algorithm.

solutionMDP* averageMDP::valueIteration(double epsilon, int maxIter){
  //Create necessary variables.
  int i; // boucle etat
  int a; // boucle action
  int nbIter; // count the iteration number

  double optimal; //keep optimal value

  int indiceEtat;
  int indiceAction;
  int argopt;

  double distance;
  double rho_n;
  double res;
  double mini,maxi;

  double *newVal; //the iterates of the value function
  double *oldVal;
  double *tmpVal; // for exchanging the value

  feedbackSolutionMDP* optimum;

  //Allocate buffers (to be used to iterate state space and action space).
  int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
  int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

  //Allocate newVal, oldVal, res and optimum action arrays.
  newVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );
  oldVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );


  optimum=new feedbackSolutionMDP();
  optimum->setSize(stateSpace->cardinal());

  // tableau qui pour chaque etat enregistre l'action optimale
  int *action = (int*) calloc( stateSpace->cardinal(), sizeof(int) );
  optimum->setAction(action);

  //Initialisation.
  nbIter = 0;

  // fill in the value vector
  // we should iterate the state space, initialize statebuffer by the first state in the state space.
  stateSpace->firstState(statebuffer);
  //Initialize V0.
  for(i = 0 ; i < stateSpace->cardinal() ; i++)
  {
    oldVal[stateSpace->index(statebuffer)] = 0;
    //Move to next state.
    stateSpace->nextState(statebuffer);
  }

  //Start value iterations.
  do
  {
    //increment the number of iterations.
    nbIter++;

    // iteration avec la valeur
    stateSpace->firstState(statebuffer);
    for(i = 0 ; i < stateSpace->cardinal() ; i++)
    {
      if(type_r == "min")    //If it is a minimization problem, initialize optimal variable with a very big value.
      {
        optimal = MAXDOUBLE;
      }
      else            //If it is a maximization problem, initialize optimal variable with a very small value.
      {
        optimal = MINDOUBLE;
      }
      indiceEtat=stateSpace->index(statebuffer);
      //To iterate the action space, initialize actionbuffer by the first action in the action space.
      actionSpace->firstState(actionbuffer);
      for(a = 0 ; a < actionSpace->cardinal() ; a++)
      {
        indiceAction= actionSpace->index(actionbuffer);
        //Calculate the gain (evaluation function) for the current state at statebuffer.
        //gain = cost per stage + transitions.evaluateValueState(values vector, results vector, state)
        // transitions.at recupere la matrice de transition dépendant de l'action
        // evalueState calcule la ligne corresondant à l'etat dans cette fonction de valeur
        //the value of the multiplication evaluateValueState is stored at res.
        //oldval i used everywhere
        res = transitions.at(indiceAction)->evaluateValueState(oldVal,indiceEtat);
        res += rewards->getEntry(indiceEtat,indiceAction) ;
        //Check if we have a new optimal (min or max depending on the problem type).
        if(((res < optimal)&&(type_r == "min"))||((res > optimal)&&(type_r == "max")))
        {
          //New optimal found, update.
          optimal = res;
          argopt = indiceAction;
        }
        //Move to next action.
        actionSpace->nextState(actionbuffer);
      }
      // fin du parcours des actions
      newVal[indiceEtat] = optimal;
      optimum->setActionIndex(indiceEtat,argopt);
      //Move to next state.
      stateSpace->nextState(statebuffer);
    }
    distance = SpanRecup(newVal, oldVal, stateSpace->cardinal(),&maxi,&mini);
    tmpVal=newVal;
    newVal=oldVal;
    oldVal=tmpVal;
    // fin while
  }while( (distance > epsilon) && (nbIter < maxIter) );
  printf("# Value Iteration Done %d iterations. Final distance = %e\n", nbIter, distance);

  // on connait la politique optimale mais pas le gain moyen
  rho_n=(mini+maxi)/2;
  // remplissage de la solution
  stateSpace->firstState(statebuffer);
  for(i = 0; i < stateSpace->cardinal() ; i++)
  {
    indiceEtat=stateSpace->index(statebuffer);
    oldVal[indiceEtat]=rho_n;
    //Move to next state.
    stateSpace->nextState(statebuffer);
  }

  // affectation de la valeur de la politique
  optimum->setValue(oldVal);

  //Free unnecessary memory blocks.
  free(newVal);
  free(actionbuffer);
  free(statebuffer);

  return optimum;
}

//A function to solve average infinite  MDP using value iteration and Gauss Seidel Improvement
solutionMDP* averageMDP::valueIterationGS(double epsilon, int maxIter){
     std::cout << "Warning not yet implemented" << std::endl;
     return NULL;
}



//A function to solve average infinite  MDP using relative value iteration algorithm.
solutionMDP* averageMDP::relativeValueIteration(double epsilon, int maxIter){
  //Create necessary variables.
  int i; // boucle etat
  int a; // boucle action
  int nbIter; // count the iteration number

  double optimal; //keep optimal value
  int argopt; // keep arg opt

  int indiceEtat;
  int indiceAction;

  double distance;
  double gain;
  double rho_n;
  double res;

  double *newVal; //the iterates of the value function
  double *oldVal;
  double *tmpVal;

  feedbackSolutionMDP* optimum;

  //Allocate buffers (to be used to iterate state space and action space).
  int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
  int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

  //Allocate newVal, oldVal, res and optimum action arrays.
  newVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );
  oldVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );

  optimum=new feedbackSolutionMDP();
  optimum->setSize(stateSpace->cardinal());

  // tableau qui pour chaque etat enregistre l'action optimale
  int *action = (int*) calloc( stateSpace->cardinal(), sizeof(int) );
  optimum->setAction(action);

  //Initialisation.
  nbIter = 0;

  // fill in the value vector
  // we should iterate the state space, initialize statebuffer by the first state in the state space.
  stateSpace->firstState(statebuffer);
  //Initialize V0.
  for(i = 0 ; i < stateSpace->cardinal() ; i++)
  {
    oldVal[stateSpace->index(statebuffer)] = 0;
    //Move to next state.
    stateSpace->nextState(statebuffer);
  }

  //Start value iterations.
  do
  {
    //increment the number of iterations.
    nbIter++;
    // first compute the expected gain
    // calcul du rho optimal sur l'etat specifique
    if(type_r == "min")    //If it is a minimization problem, initialize optimal variable with a very big value.
    {
      optimal = MAXDOUBLE;
    }
    else            //If it is a maximization problem, initialize optimal variable with a very small value.
    {
      optimal = MINDOUBLE;
    }
    //To iterate the action space, initialize actionbuffer by the first action in the action space.
    actionSpace->firstState(actionbuffer);
    for(a = 0 ; a < actionSpace->cardinal() ; a++)
    {
      indiceAction= actionSpace->index(actionbuffer);
      //Calculate the gain (evaluation function) for the current state at statebuffer.
      //gain = cost per stage + transitions.evaluateValueState(values vector, results vector, state);
      //the value of the multiplication evaluateValueState is stored at res.
      //oldval i used everywhere
      res = transitions.at(indiceAction)->evaluateValueState(oldVal,indexEtatSpecifique );
      gain = rewards->getEntry(indexEtatSpecifique,indiceAction) + res;
      //Check if we have a new optimal (min or max depending on the problem type).
      if(((gain < optimal)&&(type_r == "min"))||((gain > optimal)&&(type_r == "max")))
      {
        //New optimal found, update.
        optimal = gain;
      }
      //Move to next action.
      actionSpace->nextState(actionbuffer);
    }
    // gain moyen temporaire
    rho_n=optimal;

    // iteration avec la valeur relative
    // avec une autre fonction de cout qui inclut la valeur moyenne
    stateSpace->firstState(statebuffer);
    for(i = 0 ; i < stateSpace->cardinal() ; i++)
    {
      if(type_r == "min")    //If it is a minimization problem, initialize optimal variable with a very big value.
      {
        optimal = MAXDOUBLE;
      }
      else            //If it is a maximization problem, initialize optimal variable with a very small value.
      {
        optimal = MINDOUBLE;
      }
      indiceEtat=stateSpace->index(statebuffer);
      //To iterate the action space, initialize actionbuffer by the first action in the action space.
      actionSpace->firstState(actionbuffer);
      for(a = 0 ; a < actionSpace->cardinal() ; a++)
      {
        indiceAction= actionSpace->index(actionbuffer);
        //Calculate the gain (evaluation function) for the current state at statebuffer.
        //gain = cost per stage + transitions.evaluateValueState(values vector, results vector, state)
        // transitions.at recupere la matrice de transition dépendant de l'action
        // evalueState calcule la ligne corresondant à l'etat dans cette fonction de valeur
        //the value of the multiplication evaluateValueState is stored at res.
        //oldval i used everywhere
        res = transitions.at(indiceAction)->evaluateValueState(oldVal,indiceEtat);
        gain = rewards->getEntry(indiceEtat,indiceAction) + res ;
        gain =gain - rho_n;
        //Check if we have a new optimal (min or max depending on the problem type).
        if(((gain < optimal)&&(type_r == "min"))||((gain > optimal)&&(type_r == "max")))
        {
          //New optimal found, update.
          optimal = gain;
          argopt = indiceAction;
        }
        //Move to next action.
        actionSpace->nextState(actionbuffer);
      }
      // fin du parcours des actions
      newVal[indiceEtat] = optimal;
      // move to next space
      stateSpace->nextState(statebuffer);
    }
    distance = Span(newVal, oldVal, stateSpace->cardinal());
    // echange des pointeurs des vecteurs pour recommencer en 1.
    tmpVal=newVal;
    newVal=oldVal;
    oldVal=tmpVal;
  }while( (distance > epsilon) && (nbIter < maxIter) );
  printf("# Relative Value Iteration Done %d iterations. Final distance = %e\n", nbIter, distance);
  // Fin de la boucle de value iteration

  // Calcul de arg min ou argmax pour determiner la politique
  // Se fait sur la politique non modifiee par le biais
  // parcours des etats
  stateSpace->firstState(statebuffer);
  for(i = 0 ; i < stateSpace->cardinal() ; i++)
  {
    if(type_r == "min")    //If it is a minimization problem, initialize optimal variable with a very big value.
    {
      optimal = MAXDOUBLE;
    }
    else            //If it is a maximization problem, initialize optimal variable with a very small value.
    {
      optimal = MINDOUBLE;
    }
    indiceEtat=stateSpace->index(statebuffer);
    //To iterate the action space, initialize actionbuffer by the first action in the action space.
    actionSpace->firstState(actionbuffer);
    for(a = 0 ; a < actionSpace->cardinal() ; a++)
    {
      indiceAction= actionSpace->index(actionbuffer);
      res = transitions.at(indiceAction)->evaluateValueState(oldVal,indiceEtat);
      gain = rewards->getEntry(indiceEtat,indiceAction) + res ;
      if(((gain < optimal)&&(type_r == "min"))||((gain > optimal)&&(type_r == "max")))
      {
        //New optimal found, update.
        optimal=gain;
        argopt = indiceAction;
      }
      //Move to next action.
      actionSpace->nextState(actionbuffer);
    }
    // affectation des valeurs optimales de l etat
    optimum->setActionIndex(indiceEtat,argopt);
    newVal[indiceEtat]=rho_n;
    stateSpace->nextState(statebuffer);
  }

  // affectation de la valeur de la politique
  optimum->setValue(newVal);

  //Free unnecessary memory blocks.
  free(oldVal);
  free(actionbuffer);
  free(statebuffer);

  return optimum;
}


// policyIteration
solutionMDP* averageMDP::policyIteration(int maxIter) {
    std::cout << "Warning NOT Implemented" << std::endl;
    return NULL;
}

// variante de policy iteration
solutionMDP* averageMDP::policyIterationModified(double epsilon, int maxIter, double delta, int maxInIter){
  int i, a;
  int indiceEtat;
  int indiceAction;
  int indiceArgOpt;
  int nbIter =0; //nb iterations totales
  int nbInIter=0; // nb iterations internes

  bool flag=true;
  double distance;
  double gain; // pour calculer la gain d'une action
  double optimal; // garde l'optimal
  double res; // res phase

  double mini, maxi;
  double rho_n;

  double *tmpVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur temporaire
  double *initVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur initiale
  double *newVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) );  // fonction de valeur
  double *echVal;

  //Allocate buffers (to be used to iterate state space and action space).
  int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
  int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

  feedbackSolutionMDP *optimum;
  optimum=new feedbackSolutionMDP();
  optimum->setSize(stateSpace->cardinal());

  // tableau qui pour chaque etat enregistre l'action optimale
  int *tmpAction = (int*) calloc( stateSpace->cardinal(), sizeof(int) );

  /* boucle generique */
  do {
      //increment the number of iterations.
      nbIter++;
      /*Policy iteration phasis. */
      //To iterate the state space, initialize statebuffer by the first state in the state space.
      stateSpace->firstState(statebuffer);
      for(i = 0 ; i < stateSpace->cardinal() ; i++)
      {
          if(type_r == "min")             //If it is a minimization problem, initialize optimal variable with a very big value.
          {
              optimal = MAXDOUBLE;
          }
          else                            //If it is a maximization problem, initialize optimal variable with a very small value.
          {
              optimal = MINDOUBLE;
          }
          indiceEtat=stateSpace->index(statebuffer);
          //To iterate the action space, initialize actionbuffer by the first action in the action space.
          actionSpace->firstState(actionbuffer);
          for(a = 0 ; a < actionSpace->cardinal() ; a++)
          {
              //Calculate the gain (evaluation function) for the current state at statebuffer.
              //gain = cost per stage + transitions.evaluateValueState(values vector, results vector, state);
              indiceAction =actionSpace->index(actionbuffer);
              res = transitions.at(indiceAction)->evaluateValueState(initVal,indiceEtat);
              gain = rewards->getEntry(indiceEtat,indiceAction) + res;
              //Check if we have a new optimal (min or max depending on the problem type).
              if(((gain < optimal)&&(type_r == "min"))||((gain > optimal)&&(type_r == "max")))
              {
                  //New optimal found, update.
                  optimal = gain;
                  indiceArgOpt = actionSpace->index(actionbuffer);
              }
              //Move to next action.
              actionSpace->nextState(actionbuffer);
          }
          //Save optimal action and its associated gain value for the current state (statebuffer).
          tmpVal[indiceEtat] = optimal;
          tmpAction[indiceEtat]=indiceArgOpt;
          //Move to next state.
          stateSpace->nextState(statebuffer);
      }
      //Calculate the new distance between Vn+1 and Vn.
      distance = SpanRecup(tmpVal,initVal, stateSpace->cardinal(),&maxi,&mini );
      //Optimality check.
      if(distance < epsilon)
      {
          //Optimal found.
          flag = false;
      }
      else
      //Optimal not found, run value iteration to improve.
      {
          nbInIter=0;
          //Value iteration phase.
          do
          {
              nbInIter++;
              //To iterate the state space, initialize statebuffer by the first state in the state space.
              stateSpace->firstState(statebuffer);
              for(i = 0; i < stateSpace->cardinal() ; i++)
              {
                  indiceEtat=stateSpace->index(statebuffer);
                  indiceAction=tmpAction[indiceEtat];
                  res = transitions.at(indiceAction)->evaluateValueState(tmpVal,indiceEtat);
                  newVal[indiceEtat] = rewards->getEntry(indiceEtat,indiceAction) + res;
                  //Move to next state.
                  stateSpace->nextState(statebuffer);
              }
              distance = Span(newVal, tmpVal, stateSpace->cardinal());
              // echange des pointeurs des vecteurs pour recommencer en 1.
              echVal=newVal;
              newVal=tmpVal;
              tmpVal=echVal;
          }while ( (distance > delta ) && ( nbInIter < maxInIter) );
          // fin du calcul de la valeur
          // echange des espaces pointes
          echVal=initVal;
          initVal=tmpVal;
          tmpVal=echVal;
      } // fin de la phase iteration valeur
  // fin du step
  }while( flag && (nbIter < maxIter));
  printf("# value Iteration Modified Done %d iterations. Final distance = %e\n", nbIter, distance);

  rho_n=(mini+maxi)/2;
  // remplissage de la solution
  stateSpace->firstState(statebuffer);
  for(i = 0; i < stateSpace->cardinal() ; i++)
  {
    indiceEtat=stateSpace->index(statebuffer);
    tmpVal[indiceEtat]=rho_n;
    //Move to next state.
    stateSpace->nextState(statebuffer);
  }
  //
  optimum->setValue(tmpVal);
  optimum->setAction(tmpAction);

  free(newVal);
  if (! flag) {free(initVal);}
  free(actionbuffer);
  free(statebuffer);

  return optimum;
}

//calcul du cout moyen sachant une politique donnee
double* averageMDP::policyCost(solutionMDP *policy,double epsilon,int maxIter){
  int i;
  int nbIter=0; // count the iteration number

  int indiceEtat;
  int indiceAction;

  double distance;
  double rho_n;
  double res;
  double mini, maxi;

  int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );
  double *oldVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur temporaire
  double *newVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ) ;
  double *tmpVal;

  feedbackSolutionMDP* optimum;

  // stockage des indices de la politique optimaledynamic_cast<B*>(aPtr) != nullptr
  if (  dynamic_cast <feedbackSolutionMDP *> (policy) != NULL ){
    optimum = dynamic_cast <feedbackSolutionMDP *> (policy);
  }
  else {
    std::cout << "Wrong type of policy" << std::endl;
    return NULL;
  }

  // on part d une fonction V nulle
  //Initialize V0.
  stateSpace->firstState(statebuffer);
  for(i = 0 ; i < stateSpace->cardinal() ; i++)
  {
    oldVal[stateSpace->index(statebuffer)] = 0;
    //Move to next state.
    stateSpace->nextState(statebuffer);
  }

  // debut iteration
  do
  {
      nbIter++;
      //To iterate the state space, initialize statebuffer by the first state in the state space.
      stateSpace->firstState(statebuffer);

      for(i = 0; i < stateSpace->cardinal() ; i++)
      {
        indiceEtat= stateSpace->index(statebuffer);
        indiceAction=optimum->getActionIndex(indiceEtat);
        res = transitions.at(indiceAction)->evaluateValueState(oldVal,indiceEtat);
        newVal[indiceEtat] = rewards->getEntry(indiceEtat,indiceAction) +res;

        //Move to next state.
        stateSpace->nextState(statebuffer);
      }

      distance = SpanRecup(newVal, oldVal, stateSpace->cardinal(),&maxi,&mini);
      // echange des pointeurs des vecteurs pour recommencer en 1.
      tmpVal=newVal;
      newVal=oldVal;
      oldVal=tmpVal;
  }while((distance > epsilon) && (nbIter < maxIter) );
  std::cout << "# Power method Done. With "<< nbIter << " iterations and  final distance: " << distance << std::endl;

  rho_n=(mini+maxi)/2;
  // remplissage de la solution
  stateSpace->firstState(statebuffer);
  for(i = 0; i < stateSpace->cardinal() ; i++)
  {
    indiceEtat=stateSpace->index(statebuffer);
    oldVal[indiceEtat]=rho_n;
    //Move to next state.
    stateSpace->nextState(statebuffer);
  }


  free (newVal);
  free(statebuffer);

  return oldVal;
}
