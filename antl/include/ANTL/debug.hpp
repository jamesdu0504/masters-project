/**
 * @file debug.hpp
 * @author John Lees-Miller
 * @brief Macros and functions useful for debugging.
 *
 * You can turn off errors and "info" messages by defining \c NDEBUG. This is
 * usually done in the makefile, and also disables C assertions. For gcc, you
 * would do this by passing the \c -DNDEBUG option on the command line.
 */

#ifndef ANTL_DEBUG_H
#define ANTL_DEBUG_H

#include <cassert>
#include <iostream>
#include <vector>
#include <NTL/ZZ.h>
#include <NTL/RR.h>


// We use the NTL namespace everywhere. Rather than have a using directive in
// every file, we just put it here, for convenience and clarity.
NTL_CLIENT

namespace ANTL {

  /**
   * @brief Controls where debug output goes, and allows global enable/disable
   * on debug output streams.
   *
   * You should use the macros defined in this file to generate debugging output.
   * This way, it's easy to remove all the debug messages at compile time, by
   * defining NDEBUG. You should only use this class directly when you want to
   * configure debug output from a driver application.
   *
   * By default, \c error is \c std::cerr, and \c pinfo and \c ptrace
   * are \c std::cout.
   */
  class debug
  {
  public:

    /// Info messages will be printed to the \c info() stream.
    static void enable_info();

    /// Trace messages will be printed to the \c trace() stream.
    static void enable_trace();

    /// \c info() will return a "null" stream that directs output to nowhere.
    static void disable_info();
	
    /// \c trace() will return a "null" stream that directs output to nowhere.
    static void disable_trace();

    static void set_error(std::ostream & os);
    static void set_info(std::ostream & os);
    static void set_trace(std::ostream & os);
	
    /// Gets the error stream.
    static std::ostream & error();
	
    /// Gets the info stream, or a "null" stream if info is disabled.
    static std::ostream & info();
	
    /// Gets the trace stream, or a "null" stream if trace is disabled.
    static std::ostream & trace();
	
  private:
    static debug& singleton();
	
    debug();

    bool error_enabled, info_enabled, trace_enabled;
    std::ostream* perror;
    std::ostream* pinfo;
    std::ostream* ptrace;
    std::ostream null_stream;
  };
}

#ifdef NDEBUG

//
// Disable the debugging macros. Note that the error macros still work.
//
#define ANTL_TRACE()							// nothing
#define ANTL_TRACE_ARG(arg)						// nothing
#define ANTL_TRACE_ARG_DETAIL(arg, detail)				// nothing
#define ANTL_INFO(message)						// nothing
#define ANTL_INFO_VAR(arg)						// nothing
#define ANTL_INFO_VAR_DETAIL(arg, detail)				// nothing
#define ANTL_INFO_IF(condition, message)				// nothing
#define ANTL_INFO_VAR_IF(condition, arg)				// nothing
#define ANTL_INFO_VAR_DETAIL_IF(condition, arg, detail)	// nothing

#else // NDEBUG is _not_ defined

// __PRETTY_FUNCTION__ is a gcc extension to the standard. It gives the full
// signature, which is better than the more widely supported __func__ symbol,
// which just gives the name.
#if defined(__GNUC__)
#	define ANTL_CURRENT_METHOD_NAME __PRETTY_FUNCTION__
#else
#	define ANTL_CURRENT_METHOD_NAME __func__
#endif

/**
 * @brief Prints the current method name to the trace stream, if \c NDEBUG is
 * not defined.
 */
#define ANTL_TRACE()							\
  ANTL_BEGIN_MACRO_FUNCTION						\
  ANTL::debug::trace()							\
    << "@@:" << __FILE__ << ":" << __LINE__ << ":"			\
					       << ANTL_CURRENT_METHOD_NAME << std::endl; \
  ANTL_END_MACRO_FUNCTION

/**
 * @brief Prints message of the form "arg = value_of_arg" to the trace stream,
 * if \c NDEBUG is not defined.
 */
#define ANTL_TRACE_ARG(arg)			\
  ANTL_BEGIN_MACRO_FUNCTION			\
  ANTL::debug::trace()				\
    << "@@:\t" #arg " = " << arg << std::endl;	\
  ANTL_END_MACRO_FUNCTION

/**
 * @brief Prints message of the form "detail (arg) = value_of_arg" to the
 * trace stream, if \c NDEBUG is not defined.
 */
#define ANTL_TRACE_ARG_DETAIL(arg, detail)				\
  ANTL_BEGIN_MACRO_FUNCTION						\
  ANTL::debug::trace()							\
    << "@@:\t" << detail << " (" #arg ") = " << arg << std::endl; 	\
  ANTL_END_MACRO_FUNCTION
	
/**
 * @brief Prints \p message to the info stream, if \c NDEBUG is not defined.
 */
#define ANTL_INFO(message)						\
  ANTL_BEGIN_MACRO_FUNCTION						\
  ANTL::debug::info()							\
    << "ii:" << __FILE__ << ":" << __LINE__ << ":"			\
					       << __func__ << ": " << message << std::endl; \
  ANTL_END_MACRO_FUNCTION

/**
 * @brief Prints message of the form "var = value_of_var" to the info stream,
 * if \c NDEBUG is not defined.
 */
#define ANTL_INFO_VAR(var)			\
  ANTL_INFO( #var " = " << var )

/**
 * @brief Prints message of the form "detail (var) = value_of_var" to the
 * info stream, if \c NDEBUG is not defined.
 */
#define ANTL_INFO_VAR_DETAIL(var, detail)		\
  ANTL_INFO( detail << " (" #var ") = " << var )

/**
 * @brief Prints \p message to the info stream, if \p condition evaluates to
 * true and \c NDEBUG is not defined.
 */
#define ANTL_INFO_IF(condition, message)	\
  ANTL_BEGIN_MACRO_FUNCTION			\
  if (condition)				\
    ANTL_INFO(message);				\
  ANTL_END_MACRO_FUNCTION

/**
 * @brief Prints message of the form "var = value_of_var" to the info stream,
 * if \c NDEBUG is not defined and \p condition evaluates to true.
 */
#define ANTL_INFO_VAR_IF(condition, var)	\
  ANTL_BEGIN_MACRO_FUNCTION			\
  if (condition)				\
    ANTL_INFO_VAR(var);				\
  ANTL_END_MACRO_FUNCTION

/**
 * @brief Prints message of the form "detail (var) = value_of_var" to the
 * info stream, if \c NDEBUG is not defined and \p condition evaluates to true.
 */
#define ANTL_INFO_VAR_DETAIL_IF(condition, var, detail)	\
  ANTL_BEGIN_MACRO_FUNCTION				\
  if (condition)					\
    ANTL_INFO_VAR_DETAIL(var, detail);			\
  ANTL_END_MACRO_FUNCTION

#endif // ifdef NDEBUG

/**
 * @brief Prints \p message to the error stream.
 *
 * Errors are always printed, regardless of whether \c NDEBUG is set.
 */
#define ANTL_ERROR(message)						\
  ANTL_BEGIN_MACRO_FUNCTION						\
  ANTL::debug::error()							\
    << "**:" << __FILE__ << ":" << __LINE__ << ":"			\
					       << __func__ << ": " << message << std::endl; \
  ANTL_END_MACRO_FUNCTION
	
/**
 * @brief Prints \p message to the error stream and exits the application with
 * the given \p returnCode. \p returnCode must be an integer.
 */
#define ANTL_FATAL_RC(returnCode, message)	\
  ANTL_BEGIN_MACRO_FUNCTION			\
  ANTL_ERROR(message);				\
  std::exit(returnCode);			\
  ANTL_END_MACRO_FUNCTION

/**
 * @brief Prints \p message to the error stream and exits the application with
 * the given \p returnCode. \p returnCode must be an integer.
 */
#define ANTL_FATAL(message)			\
  ANTL_FATAL_RC(1, message)

//
// Helpful Output Methods
//

namespace ANTL
{

  /// Print bits of \p z to \p os, from highest order to lowest order.
  void print_binary_expansion(std::ostream & os, const ZZ & z);

  /// Prints bits of \p r to \p os, from highest order to lowest order.
  void print_binary_expansion(std::ostream & os, const RR & r);

  /// Outputs an STL vector in the form [ v_0 v_1 ... v_n ].
  template <class T>
  std::ostream & operator << (std::ostream & os,
			      const std::vector<T> & vector)
  {
    os << "[ ";
    for (size_t i = 0; i < vector.size(); ++i) os << vector[i] << " ";
    os << "]";
    return os;
  }

} // ANTL

#endif // guard

