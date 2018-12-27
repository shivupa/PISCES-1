#ifndef PISCES_TIMER_HPP
#define PISCES_TIMER_HPP

#include <string>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#else
#include <ctime>
#endif

/// \class timer \brief Provide simple timing facility (uses OpenMP_Wtime if available)

#if defined(_OPENMP)
class timer
{
public:
   timer() { restart(); }
   void   restart() { _start = omp_get_wtime(); }
   double elapsed() const { return  omp_get_wtime() - _start; }
private:
   double _start;
};


#else 


class timer
{
 public:
    timer() {restart();}
    void   restart() {_start = std::clock();}
    double elapsed() const { return double(std::clock()-_start)/CLOCKS_PER_SEC; }

 private:
    std::clock_t _start;
}; // timer

#endif // _OPENMP


/// Writes to std::cout the total time of scoping block
class progress_timer
{
public:
  progress_timer(const std::string& name="_", int verbose = 2)
      : m_name (name) 
      , m_verb (verbose)
      {m_name.resize(20,' ');}
   ~progress_timer() {
      if ( m_verb >= 0 )
         std::cout << "  time: " << m_name << ": " << m_t.elapsed()  << std::endl;
   }

private:
   timer m_t;
   std::string m_name;
  int m_verb;
};

#endif  // BOOST_TIMER_HPP
