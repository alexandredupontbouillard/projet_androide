
class solutionMDP
{
public:

  /**
  * @brief Constructor to create solutionMDP object.
  * @author Hyon
  * @version 1
  * @date feb 2018
  */
  solutionMDP();

  /** 
  * @brief the destructor of the object solutionMDP
  * @author EH
  * @version 1
  * @date feb 2018
  */
  virtual ~solutionMDP();


  /**
  * @brief A function to print the solution object.
  * @author Abood Mourad.
  * @date feb 2018
  * @return none.
  */
   virtual void writeSolution();
   
   
  /**
  * @brief setter of the size.
  * @author EH
  * @date mar 2018
  * @return none.
  */ 
  void setSize(int s);
     

protected:
  int _size;                         /**< _size of the solution */
};






