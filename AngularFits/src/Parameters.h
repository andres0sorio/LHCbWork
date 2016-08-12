// $Id: $
#ifndef PARAMETERS_H 
#define PARAMETERS_H 1

// Include files

/** @class Parameters Parameters.h
 *  
 *
 *  @author Andres Osorio Oliveros
 *  @date   2006-10-31
 */
class Parameters {
public: 
  /// Standard constructor
  Parameters( ) {};
  
  Parameters( const char * pname , double pvalue ) 
  {
    m_name  = std::string( pname );
    m_value = pvalue;
  }
  
  virtual ~Parameters( ); ///< Destructor
  
  std::string name() 
  { 
    return m_name;
  }
  
  double value() 
  {
    return m_value;
  }
  
  void updateValue( double val ) 
  {
    m_value = val;
  }
  
  
protected:
  
private:
  
  std::string m_name;
  
  double m_value;
  
  
};
#endif // PARAMETERS_H
