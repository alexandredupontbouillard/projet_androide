class feedbackSolutionMDP : public solutionMDP
{
public:
   /**
   * @brief Constructor to create feedbackSolution object.
   * @author EH
   *
   * Constructeur qui initialise deux vecteurs null pour les deux attributs
   */
    feedbackSolutionMDP();

   /**
   * @brief destructor to delete feedbackSolution object.
   * @author AM EH
   */
    ~feedbackSolutionMDP();

   /**
   * @brief setter of the action vector
   * @author AM EH
   * @date feb 2018
   * @param a : int * . Pointeur sur un vecteur d'entiers de dimension le cardinal de l'espace d etat.
   *  Chaque entree représente represente l'index de l'action à
   *  effectuer quand on est dans l'etat dont l'index est la coordonnée.
   * @return none.
   */
    void setAction(int* a);

   /**
   * @brief setter of the value function vector
   * @author AM EH
   * @date feb 2018
   * @param t : double * Pointeur sur un vecteur de double de dimension le cardinal de l'espace d etat. Chaque entree représente  la valeur
   *  est la valeur de la fonction de valeur dans l etat dont l'index est la coordonnée du vecteur
   * @return none.
   */
    void setValue(double* t);

   /**
   * @brief setter of only one action at index
   * @author AM EH
   * @date feb 2018
   * @param indice : int index de la valeur de l'action à modifier
   * @param value : int nouvelle valeur
   * @return none.
   */
    void setActionIndex(int indice,int value);

    /**
   * @brief getter of only one action at index
   * @author AM EH
   * @date feb 2018
   * @param indice : int index de la valeur de l'action à retourner
   * @return int the index of the action
   */
    int getActionIndex(int indice);

  /**
   * @brief getter of only one value at index
   * @author EH
   * @date feb 2019
   * @param indice : int index de la valeur à retourner
   * @return double the value at indice
   */
    double getValueIndex(int indice);


    /**
    * @brief A function to print the solution object.
    * @author Abood Mourad.
    * @return none.
    */
     void writeSolution();

//variables specifiques de la classe
protected:
    int *action;      /**< A vector of actions representing the optimal policy. */
    double *value;    /**< A vector of associated values. */
};

