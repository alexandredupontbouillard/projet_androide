class marmoteInterval : public marmoteSet {

 protected:
  int _min; /**< the lower end of the interval */
  int _max; /**< the higher end of the interval */

 public:
  /**
   * Constructor for an interval. By convention, if max < min, then the
   * interval is empty. Otherwise, both min and max are inside the interval.
   * @author Alain Jean-Marie
   * @param min the low end of the interval
   * @param max the high end of the interval
   */
  marmoteInterval( int min, int max );
  /**
   * Standard destructor for an interval.
   * @author Alain Jean-Marie
   */
  ~marmoteInterval();

 public:
  /**
   * @brief Test if the set is finite. These sets always are.
   *
   * @return true
   */
  bool isFinite() { return true; };
  /**
   * @copydoc marmoteSet::isZero(int*)
   */
  bool isZero(int* buffer);
  /**
   * @brief Initializes some state buffer with the first state of the set.
   * @param buffer the buffer to be set.
   */
  void firstState(int* buffer);
  /**
   * @copydoc marmoteSet::nextState()
   */
  void nextState(int *buffer);
  /**
   * @copydoc marmoteSet::decodeState(int,int*)
   */
  void decodeState(int index, int* buffer);
  /**
   * @copydoc marmoteSet::index(int*)
   */
  int index(int* buffer);
  /**
   * @copydoc marmoteSet::printState()
   */
  void printState(FILE* out, int *buffer);

 public:
  /**
   * @brief Enumeration procedure.
   *
   */
  void enumerate();

};


