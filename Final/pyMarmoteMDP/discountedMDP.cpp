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
 * @brief Class dicountedMDP source
 * @author Hyon, lip6
 * @date 18 Jan 2018
 * @version 3
 *
 * This .cpp file is for implementing functions of the discountedMDP class
 */

#include <values.h>
#include <iostream>
#include "header/discountedMDP.h"
#include "alglin.h"

using namespace std;


// constructeurs

discountedMDP::discountedMDP(std::string r, marmoteSet* states, marmoteSet* actions, std::vector<sparseMatrix*> trans, sparseMatrix* rews, double b)
: genericMDP("discrete",r,states,actions,trans,rews), beta(b), type_c("infinite discounted")
{

}


discountedMDP::discountedMDP(string r, marmoteSet* states, marmoteSet* actions, vector<sparseMatrix*> trans, vector<sparseMatrix*> rews, double b)
: genericMDP(string("discrete"),r,states,actions,trans,rews), beta(b), type_c("infinite discounted")
{


}

// destructeur
discountedMDP::~discountedMDP(){
}

// other utilities function
void discountedMDP::writeMDP()
{
    genericMDP::writeMDP();
    printf("#############################################\n");
    printf("write Infinite Discounted MDP\n");
    printf("MDP Criteria : %s \n" , type_c.c_str());
    printf("Discount factor: %f\n" , beta );
}


// solving discounted MDP

//A function to solve discounted infinite  MDP using value iteration algorithm.
solutionMDP* discountedMDP::valueIteration(double epsilon, int maxIter)
{
  //Create necessary variables.
    int i;
    int a;
    int nbIter; // count the iteration number
    double optimal; //keep optimal value
    int argopt;
    int indiceEtat;
    int indiceAction;

    double distance;
    double gain;
    double res;
    double borne;

    double *newVal;
    double *oldVal;
    double *tempVal;
    
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



    //Initialization.
    nbIter = 0;
    // calcul de la borne
    // borne donnee par powell avec arrangement pour les cas grands
    if ( (1-beta)/(2*beta) > 2) {
        borne = 2* epsilon;
    }
    else {
      borne= epsilon* ((1-beta)/(2*beta));
    }
    //To iterate the state space, initialize statebuffer by the first state in the state space.
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
        //To iterate the state space, initialize statebuffer by the first state in the state space.
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
                //gain = cost per stage + transitions.evaluateValueState(values vector, results vector, state);
                //the value of the multiplication evaluateValueState is stored at res.
                //one alternates between oldval and newval to avoid the copy
                res = transitions.at(indiceAction)->evaluateValueState(oldVal, indiceEtat);
                gain = rewards->getEntry(indiceEtat,indiceAction) + (beta*res);
                
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
            //Save optimal action and its associated gain value for the current state (statebuffer).
            newVal[indiceEtat] = optimal;
            optimum->setActionIndex(indiceEtat,argopt);
            //Move to next state.
            stateSpace->nextState(statebuffer);
        }
        //Calculate the new distance between Vn+1 and Vn.
        distance = Norm(newVal, oldVal, stateSpace->cardinal());
        // exchange 
        tempVal=oldVal;
        oldVal=newVal;
        newVal=tempVal;
    }while((distance > borne) && (nbIter < maxIter));
    printf("# Value Iteration Done %d iterations. Final distance = %e\n", nbIter, distance);


    optimum->setValue(oldVal);
    //Free unnecessary memory blocks.
    free(newVal);
    

    free(actionbuffer);
    free(statebuffer);

    //return the resulting solution.
    return optimum;
}

solutionMDP* discountedMDP::valueIterationGS(double epsilon, int maxIter){
     std::cout << "Warning not yet implemented" << std::endl;
     return NULL;
}

// A function to solve (discrete time) MDP using policy iteration algorithm.
solutionMDP* discountedMDP::policyIteration(int maxIter) {
    int i, a;
    int indiceEtat;
    int indiceAction;
    int indiceArgOpt; // to keep the optimal action
    int nbIter =0; //nb iterations totales
    double res =0; // pour stocker des resultats temporaires
    double gain; // pour calculer la gain d'une action
    double optimal; // garde l'optimal
    bool differencePol; // indicates that there is a difference between the other term

    double *tmpVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur temporaire

    //Allocate buffers (to be used to iterate state space and action space).
    int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
    int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

    // solutions
    feedbackSolutionMDP *optimum;
    optimum=new feedbackSolutionMDP();
    optimum->setSize(stateSpace->cardinal());
    // tableau qui pour chaque etat enregistre l'action optimale
    // c'est lui qui sera utilisé  pour pi
    int *tmpAction = (int*) calloc( stateSpace->cardinal(), sizeof(int) );
    // ancien tableau qui stocke les valeurs optimale
    int *tmpActionPrec = (int*) calloc( stateSpace->cardinal(), sizeof(int) );

    /*Policy selection phase. */
    //One iterates the state space, one initializes statebuffer by the first state in the state space.
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
        res = transitions.at(indiceAction)->evaluateValueState(tmpVal,indiceEtat);
        gain = rewards->getEntry(indiceEtat,indiceAction) + (beta*res);
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
    // On commence la boucle
    do{
      nbIter++;
      // on teste si y a des differences entre politiques
      differencePol=false;
      stateSpace->firstState(statebuffer);
      for(i=0;i < stateSpace->cardinal() ; i++)
      {
        indiceEtat=stateSpace->index(statebuffer);
        if (tmpAction[indiceEtat]!=tmpActionPrec[indiceEtat]){
          // Il y a une difference entre politiques
          differencePol=true;
        }
        // changer la matrice de gain
      }


    } while((differencePol) && (nbIter<maxIter));
    //

    return NULL;
}

/* *********************************************************************** */
// A function to evaluate the cost of a policy with iteration of the power.
solutionMDP* discountedMDP::policyIterationModified(double epsilon, int maxIter, double delta, int maxInIter){
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
    double borne; // sur l'ecart;
    double res; // res phase

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

    // calcul de la borne
    // borne donnee par powell avec arrangement pour les cas grands
    if ( (1-beta)/(2*beta) > 2) {
        borne = 2* epsilon;
    }
    else {
      borne= epsilon* ((1-beta)/(2*beta));
    }

    // initialisation
    // on doit avoir LV0 > V0 pour max et LV0 < V0 pour min
    // choix (1-gamma)^{-1} min_s min_a pour un  pb de maximisation
    // choix (1-gamma)^{-1} max_s max_a pour un  pb de minimisation
    // attention optimal fonctionne differement que dans value iteration

    stateSpace->firstState(statebuffer);
    if(type_r == "min")    //If it is a minimization problem, initialize optimal variable with a very big value.
    {
        optimal = MINDOUBLE;
    }
    else            //If it is a maximization problem, initialize optimal variable with a very small value.
    {
        optimal = MAXDOUBLE;
    }
    for(i = 0 ; i < stateSpace->cardinal() ; i++)
    {
        //To iterate the action space, initialize actionbuffer by the first action in the action space.
        actionSpace->firstState(actionbuffer);
        for(a = 0 ; a < actionSpace->cardinal() ; a++)
        {
            gain = rewards->getEntry(stateSpace->index(statebuffer),actionSpace->index(actionbuffer));
            //Check if we have a new optimal (min or max depending on the problem type).
            if(((gain > optimal)&&(type_r == "min"))||((gain < optimal)&&(type_r == "max")))
            {
                //New optimal found, update.
                optimal = gain;
            }
        }
    }
    for(i = 0 ; i < stateSpace->cardinal() ; i++)
    {
        initVal[stateSpace->index(statebuffer)] = optimal/(1-beta);
    }

    // fin  init
    /* boucle generique */
    do {
        //increment the number of iterations.
        nbIter++;
        /*Policy iteration phase. */
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
                gain = rewards->getEntry(indiceEtat,indiceAction) + (beta*res);
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
        distance = Norm(tmpVal,initVal, stateSpace->cardinal());
        //Optimality check.
        if(distance < borne)
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
                nbInIter+=2;
                //To iterate the state space, initialize statebuffer by the first state in the state space.
                stateSpace->firstState(statebuffer);
                for(i = 0; i < stateSpace->cardinal() ; i++)
                {
                    indiceEtat=stateSpace->index(statebuffer);
                    indiceAction=tmpAction[indiceEtat];
                    res = transitions.at(indiceAction)->evaluateValueState(tmpVal,indiceEtat);
                    newVal[indiceEtat] = rewards->getEntry(indiceEtat,indiceAction) + (beta*res);
                    //Move to next state.
                    stateSpace->nextState(statebuffer);
                }
                // converse step to avoid a inefficient copy
                stateSpace->firstState(statebuffer);
                for(i = 0; i < stateSpace->cardinal() ; i++)
                {
                    indiceEtat=stateSpace->index(statebuffer);
                    indiceAction=tmpAction[indiceEtat];
                    res = transitions.at(indiceAction)->evaluateValueState(newVal,indiceEtat);
                    tmpVal[indiceEtat] = rewards->getEntry(indiceEtat,indiceAction) + (beta*res);
                    //Move to next state.
                    stateSpace->nextState(statebuffer);
                }
                //Calculate the new distance between Vmn+1 and Vmn.
                distance = Norm(newVal, tmpVal, stateSpace->cardinal());
            }while ( (distance > delta ) && ( nbInIter < maxInIter) );
            // fin du calcul de la valeur
            // echange des espaces pointes
            echVal=initVal;
            initVal=tmpVal;
            tmpVal=echVal;
        } // fin de la phase iteration valeur
    // fin du step
    }while((flag) && (nbIter < maxIter));
    printf("# Done %d iterations. Final distance = %e\n", nbIter, distance);

    //
    optimum->setValue(tmpVal);
    optimum->setAction(tmpAction);

    free(newVal);
    free(initVal);
    free(actionbuffer);
    free(statebuffer);

    return optimum;
}


//A function to evaluate the cost of a policy with iteration of the power
double *discountedMDP::policyCost(solutionMDP *policy,double epsilon,int maxIter){
    int i;
    int indice;
    int nbIter;

    double res=0;
    double borne;

    double norme1;
    double temponorme;
    double vieilleval;

    double *tmpVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur
    double *initVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur initiale

    
    int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

    feedbackSolutionMDP* optimum;


    // calcul de la borne
    // borne donnee par powell avec arrangement pour les cas grands
    if ( (1-beta)/(2*beta) > 2) {
        borne = 2* epsilon;
    }
    else {
      borne= epsilon* ((1-beta)/(2*beta));
    }


    // stockage des indices de la politique optimaledynamic_cast<B*>(aPtr) != nullptr
    if (  dynamic_cast <feedbackSolutionMDP *> (policy) != NULL ){
        optimum = dynamic_cast <feedbackSolutionMDP *> (policy);
    }
    else {
        std::cout << "Wrong type of policy" << std::endl;
        return NULL;
    }

    // calcul de la valeur optimale iteration 1
    //To iterate the state space, initialize statebuffer by the first state in the state space.
    stateSpace->firstState(statebuffer);
    for(i = 0; i < stateSpace->cardinal() ; i++)
    {
        indice=stateSpace->index(statebuffer);
        res = transitions.at(optimum->getActionIndex(indice))->evaluateValueState(initVal,indice);
        tmpVal[indice] = rewards->getEntry(indice,optimum->getActionIndex(indice)) + (beta*res);
        //Move to next state.
        stateSpace->nextState(statebuffer);
    }

    // calcul avec la methode de gauss seidel
    nbIter=1;
    do
    {  nbIter++;
       //To iterate the state space, initialize statebuffer by the first state in the state space.
       stateSpace->firstState(statebuffer);
       norme1=0;
       for(i = 0; i < stateSpace->cardinal() ; i++)
       {
          indice= stateSpace->index(statebuffer);
          vieilleval=tmpVal[indice];
          res = transitions.at(optimum->getActionIndex(indice))->evaluateValueState(tmpVal,indice);
          tmpVal[indice] = rewards->getEntry(indice,optimum->getActionIndex(indice)) + (beta*res);
          //Calculate the new distance between Vn+1 and Vn.
          temponorme= fabs(tmpVal[indice]-vieilleval);
          if (temponorme > norme1) { norme1=temponorme;}
          //Move to next state.
          stateSpace->nextState(statebuffer);
        }
    }while((norme1 > borne) && (nbIter < maxIter) );
    printf("# Cost policy Done %d iterations. Final distance = %e\n", nbIter, norme1);

    free(statebuffer);
    free(initVal);

    return tmpVal;

}

//A function to evaluate the cost of a policy with iteration of the power and a scan set with indexes
double *discountedMDP::policyCostbyIndex(solutionMDP *policy,double epsilon,int maxIter){
    int i;
    double res=0;
    double borne;
    int nbIter;

    double norme1;
    double temponorme;
    double vieilleval;

    double *tmpVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur
    double *initVal = (double*) calloc( stateSpace->cardinal(), sizeof(double) ); // fonction de valeur initiale

    int* actionbuffer = (int*) calloc( actionSpace->totNbDims(), sizeof(int) );
    int* statebuffer = (int*) calloc( stateSpace->totNbDims(), sizeof(int) );

    feedbackSolutionMDP* optimum;


    // calcul de la borne
    // borne donnee par powell avec arrangement pour les cas grands
    if ( (1-beta)/(2*beta) > 2) {
        borne = 2* epsilon;
    }
    else {
      borne= epsilon* ((1-beta)/(2*beta));
    }


    // stockage des indices de la politique optimaledynamic_cast<B*>(aPtr) != nullptr
    if (  dynamic_cast <feedbackSolutionMDP *> (policy) != NULL ){
        optimum = dynamic_cast <feedbackSolutionMDP *> (policy);
    }
    else {
        std::cout << "Wrong type of policy" << std::endl;
        return NULL;
    }

    // calcul de la valeur optimale iteration 1
    for(i = 0; i < stateSpace->cardinal() ; i++)
    {
        stateSpace->decodeState(i,statebuffer);
        if (i != stateSpace->index(statebuffer)){
            std::cout << "Wrong index translation" << std::endl;
        }
        res = transitions.at(optimum->getActionIndex(i))->evaluateValueState(initVal,i);
        tmpVal[i] = rewards->getEntry(i,optimum->getActionIndex(i)) + (beta*res);

    }

    // calcul avec la methode de gauss seidel
    nbIter=1;
    do
    {
       nbIter++;
       //To iterate the state space, initialize statebuffer by the first state in the state space.
       norme1=0;
       for(i = 0; i < stateSpace->cardinal() ; i++)
       {
          stateSpace->decodeState(i,statebuffer);
          vieilleval=tmpVal[i];
          res = transitions.at(optimum->getActionIndex(i))->evaluateValueState(tmpVal,i);
          tmpVal[i] = rewards->getEntry(i,optimum->getActionIndex(i)) + (beta*res);
          //Calculate the new distance between Vn+1 and Vn.
          temponorme= fabs(tmpVal[i]-vieilleval);
          if (temponorme > norme1) { norme1=temponorme;}
          //Move to next state.
        }
    }while((norme1 > borne) && (nbIter < maxIter) );
    printf("# Cost policy Done %d iterations. Final distance = %e\n", nbIter, norme1);

    free(actionbuffer);
    free(statebuffer);
    free(initVal);

    return tmpVal;

}
