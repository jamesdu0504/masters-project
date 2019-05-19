/**
 * @file debug.cpp
 * @author John Lees-Miller
 * @brief implementation of debug class
 */

#include <ANTL/debug.hpp>

using namespace std;

namespace ANTL
{
  /**
   * The debug class is a singleton. Only one instance exists, and it is
   * accessible through this private method.
   */
  debug& debug::singleton()
  {
    static debug debug_singleton;
    return debug_singleton;
  }

  /**
   * Mechanism for a streambuf that creates a stream that does nothing.
   */
  class nullbuf : public std::streambuf
  {
  public:
    int_type overflow(int_type c)
    {
      return traits_type::not_eof(c);
    }
  } debug_nullbuf; // A module-scope variable.

  /**
   * Default constructor for debug class (sets defaults according to NDEBUG).
   */ 
  debug::debug()
    : null_stream(&debug_nullbuf)
  {
    // Set stream pointers to defaults.
    this->perror = &cerr;
    this->pinfo = &cout;
    this->ptrace = &cout;
	
    // Enable streams depending on NDEBUG.
#ifdef NDEBUG
    this->error_enabled = true;
    this->info_enabled = false;
    this->trace_enabled = false;
#else
    this->error_enabled = true;
    this->info_enabled = true;
    this->trace_enabled = true;
#endif
  }

  //
  // Implementation of debug members.
  //
  void debug::enable_info() 	{ debug::singleton().info_enabled = true; }
  void debug::enable_trace() 	{ debug::singleton().trace_enabled = true; }

  void debug::disable_info() 	{ debug::singleton().info_enabled = false; }
  void debug::disable_trace() { debug::singleton().trace_enabled = false; }

  void debug::set_error(ostream& os)	{ debug::singleton().perror = &os; }
  void debug::set_info(ostream& os) 	{ debug::singleton().pinfo = &os; }
  void debug::set_trace(ostream& os)	{ debug::singleton().ptrace = &os; }

  ostream & debug::error()
  {
    assert(debug::singleton().perror);
    return *(debug::singleton().perror);
  }

  ostream & debug::info()
  { 
    if (debug::singleton().info_enabled)
      {
	assert(debug::singleton().pinfo);
	return *(debug::singleton().pinfo);
      }
    return debug::singleton().null_stream;
  }

  ostream & debug::trace()
  {
    if (debug::singleton().trace_enabled)
      {
	assert(debug::singleton().ptrace);
	return *(debug::singleton().ptrace);
      }
    return debug::singleton().null_stream;
  }

  //
  // Helpful Output Methods
  //
  void print_binary_expansion(ostream & os, const ZZ & z)
  {
    long n = NumBits(z);
    while(n > 0)
      os << NTL::bit(z, --n);
  }

  void print_binary_expansion(ostream & os, const RR & r)
  {
    ZZ m = r.mantissa();
    long n = NumBits(m);
    long bin_pt_index = -r.exponent() - 1;
    for (int i = n - 1; i > bin_pt_index; --i)
      os << NTL::bit(m, i);
    os << ".";
    for (int i = bin_pt_index; i >= 0; --i)
      os << NTL::bit(m, i);
  }

} // ANTL

