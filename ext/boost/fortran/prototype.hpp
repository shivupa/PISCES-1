#ifndef BOOST_FORTRAN_PROTOTYPE_H_
#define BOOST_FORTRAN_PROTOTYPE_H_


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include <boost/preprocessor/detail/is_unary.hpp>
#include <boost/preprocessor/if.hpp>
#include <boost/preprocessor/control/expr_iif.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/seq/cat.hpp>
#include <boost/preprocessor/logical/bitand.hpp>
#include <boost/preprocessor/facilities/identity.hpp>
#include <cstring>
#include <string>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// compiler specific configurations 
// TODO Come up with cannonical configs
#include <boost/config.hpp>
#ifdef BOOST_MSVC
#  pragma warning( push )
#  pragma warning( disable : 4003 )
#endif

//==============================================================================
// CONFIGURATION MACROS 
//    Canonical definition of portability macros based on platform.
//==============================================================================


// Casing
// Append underscore
// Extra underscore
// Split char
// Calling convention

#ifdef unix
#   define BOOST_FORTRAN_DECL_SPLITCHAR   1
#   define BOOST_FORTRAN_DECL_UNDERSCORE  1
#   define BOOST_FORTRAN_DECL_CAPS        0
#   define BOOST_FORTRAN_DECL_XUNDERSCORE 1
#   define BOOST_FORTRAN_DECL_CALLCONV 
#else 
#   define BOOST_FORTRAN_DECL_SPLITCHAR   0
#   define BOOST_FORTRAN_DECL_UNDERSCORE  0
#   define BOOST_FORTRAN_DECL_CAPS        1
#   define BOOST_FORTRAN_DECL_XUNDERSCORE 0
#   define BOOST_FORTRAN_DECL_CALLCONV 
#endif




//==============================================================================
// IMPLEMENTATION cruft
//==============================================================================
#define BOOST_FORTRAN_decorate(name,NAME,hasus) \
   BOOST_PP_SEQ_CAT( \
      (BOOST_PP_IIF( BOOST_FORTRAN_DECL_CAPS, NAME, name ))    \
      (BOOST_PP_EXPR_IIF(BOOST_FORTRAN_DECL_UNDERSCORE,_))     \
      (BOOST_PP_EXPR_IIF(BOOST_PP_BITAND( hasus, BOOST_FORTRAN_DECL_XUNDERSCORE ), _ ))\
      )

//  namespace for the extern prototypes
#define BOOST_FORTRAN_namespace boost_fortran_impl

//==============================================================================
// CHAR_PROXY   Proxy struct to enable splitting character content and length
//==============================================================================
namespace boost   {
namespace fortran {
namespace detail  {

   
   typedef char char_t;  // Character type in Fortran (usually 1 byte char unsigned)
   typedef unsigned int charlen_t;  // Character length type in Fortran (usually 4 byte unsigned)
   

   
   struct char_proxy
   {
      char_proxy( char_t* c, size_t sz=size_t(-1) ) : p(c) , n(static_cast<charlen_t>(sz==size_t(-1)? std::strlen(c):sz)) {}
      char_proxy( ::std::string& s, size_t sz=size_t(-1) ) : p(&s[0]) , n(static_cast<charlen_t>(s.size())) {}
      char_t* p;
      charlen_t n;
   };
   struct char_const_proxy
   {
      char_const_proxy( const char_t* c, size_t sz=size_t(-1) ) : p(c) , n(static_cast<charlen_t>(sz==size_t(-1)? std::strlen(c):sz)) {}
      char_const_proxy( const ::std::string& s, size_t sz=size_t(-1) ) : p(s.c_str()) , n(static_cast<charlen_t>(s.size())) {}
      const char_t* p;
      charlen_t n;
   };

}  // namespace detail
   
typedef detail::char_proxy char_ref;
typedef detail::char_const_proxy char_const_ref;
} // namespace fortran
} // namespace mytools




#if BOOST_FORTRAN_DECL_SPLITCHAR 
#  define BOOST_FORTRAN_extern_char_t ::boost::fortran::detail::char_t
#else 
#  define BOOST_FORTRAN_extern_char_t ::boost::fortran::detail::char_proxy
#endif


//==============================================================================
//  CHARACTER  Part of the interface, to be used to prototype character 
//        arguments.
//
//  The character inside the parenthesis indicates constness. See discussion here 
//  http://tinyurl.com/2u8dwt .
//
//==============================================================================
#define BOOST_FORTRAN_CHARACTER (0)
#define BOOST_FORTRAN_CHARACTER_CONST (1)

//==============================================================================
//  IS_CHAR(x)  will return 1 if x is FORTRAN_CHAR, 0 otherwise. 
// 
//  This works since x is supposed to be BOOST_FORTRAN_CHAR or some type 
//  declaration which never comes in the form of a tuple.
// 
//  If CHARACTER arguments need not be split, this will always return 0, which 
//  is OK.     
//==============================================================================
#define BOOST_FORTRAN_is_char BOOST_PP_IS_UNARY
#define BOOST_FORTRAN_char_is_const(x) BOOST_PP_IDENTITY x ()

//==============================================================================
// extern_arg_type  Convert `i'th element of the formal argument type array 
//                  `args' to argument type in extern declaration.      
//==============================================================================
#if BOOST_FORTRAN_DECL_SPLITCHAR

#  define BOOST_FORTRAN_extern_arg(r, data, i, elem)                    \
   BOOST_PP_COMMA_IF(i)                                                 \
   BOOST_PP_IIF( BOOST_FORTRAN_is_char(elem),                           \
                 BOOST_PP_IIF( BOOST_FORTRAN_char_is_const(elem),       \
                               const ::boost::fortran::detail::char_t*, \
                               ::boost::fortran::detail::char_t*),      \
                 elem ) 
                      
#  define BOOST_FORTRAN_actual_arg(r,data,i,elem) \
   BOOST_PP_COMMA_IF(i)                           \
   BOOST_PP_IIF( BOOST_FORTRAN_is_char(elem),     \
                 v##i.p ,                         \
                 v##i   )

#  define BOOST_FORTRAN_extern_arg_extra(r, data, i, elem)              \
   BOOST_PP_COMMA_IF( BOOST_FORTRAN_is_char(elem) )                     \
   BOOST_PP_EXPR_IIF( BOOST_FORTRAN_is_char(elem),                      \
                      ::boost::fortran::detail::charlen_t )

#  define BOOST_FORTRAN_actual_arg_extra(r,data,i,elem)        \
   BOOST_PP_COMMA_IF( BOOST_FORTRAN_is_char(elem) )            \
   BOOST_PP_EXPR_IIF( BOOST_FORTRAN_is_char(elem), v##i.n )

#else // no SPLITCHAR

#  define BOOST_FORTRAN_extern_arg(r, data, i, elem)                    \
   BOOST_PP_COMMA_IF( i )                                               \
   BOOST_PP_IIF( BOOST_FORTRAN_is_char(elem),                           \
                 BOOST_PP_IIF( BOOST_FORTRAN_char_is_const(elem),       \
                               ::boost::fortran::detail::char_const_proxy, \
                               ::boost::fortran::detail::char_proxy),   \
                 elem ) 

#  define BOOST_FORTRAN_actual_arg(r,data,i,elem) \
   BOOST_PP_COMMA_IF(i) v##i


#  define BOOST_FORTRAN_extern_arg_extra(r, data, i, elem) /*empty*/              
#  define BOOST_FORTRAN_actual_arg_extra(r, data, i, elem) /*empty*/

#endif


#define BOOST_FORTRAN_wrapper_arg(r,data,i,elem)                        \
   BOOST_PP_COMMA_IF(i)                                                 \
   BOOST_PP_IIF( BOOST_FORTRAN_is_char(elem),                           \
                 BOOST_PP_IIF( BOOST_FORTRAN_char_is_const(elem),       \
                               ::boost::fortran::char_const_ref,        \
                               ::boost::fortran::char_ref),             \
                 elem )  v##i





//==============================================================================
// Generate declaration of the extern (actual) routine. 
//==============================================================================
#define BOOST_FORTRAN_extern_decl(args) \
   BOOST_PP_SEQ_FOR_EACH_I( BOOST_FORTRAN_extern_arg, ~, args )  \
   BOOST_PP_SEQ_FOR_EACH_I( BOOST_FORTRAN_extern_arg_extra, ~, args )

//==============================================================================
// Generate declaration of the wrapper function.
//==============================================================================
#define BOOST_FORTRAN_wrapper_decl(args) \
   BOOST_PP_SEQ_FOR_EACH_I(BOOST_FORTRAN_wrapper_arg, ~, args)

//==============================================================================
// Generate a call to the extern routine. 
//==============================================================================
#define BOOST_FORTRAN_extern_call(args)                          \
   BOOST_PP_SEQ_FOR_EACH_I(BOOST_FORTRAN_actual_arg, ~, args)    \
   BOOST_PP_SEQ_FOR_EACH_I(BOOST_FORTRAN_actual_arg_extra, ~, args)








//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//  MACROS provided as a part of the library interface:  
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------




//==============================================================================
// PROTOTYPE  Generic prototyping macro
//==============================================================================
#define BOOST_FORTRAN_PROTOTYPE(cxxname, fname, FNAME, has_us, ret, not_void, args) \
   namespace BOOST_FORTRAN_namespace { \
   extern "C" BOOST_PP_IIF(not_void, ret, void)                         \
      BOOST_FORTRAN_DECL_CALLCONV                                       \
      BOOST_FORTRAN_decorate(fname,FNAME,has_us)                        \
      (BOOST_FORTRAN_extern_decl(args))                                 \
      ;                                                                 \
      }\
   inline BOOST_PP_IIF(not_void, ret, void)                             \
      cxxname (BOOST_FORTRAN_wrapper_decl(args))                        \
   {                                                                    \
      BOOST_PP_EXPR_IIF(not_void, return)                               \
      BOOST_FORTRAN_namespace :: \
	 BOOST_FORTRAN_decorate(fname,FNAME,has_us)                     \
	 (BOOST_FORTRAN_extern_call(args))                              \
	 ;                                                              \
   }

//==============================================================================
// Convenience macros (special case PROTOTYPEs)
//==============================================================================
#define BOOST_FORTRAN_FUNCTION(ret, cxxname, fname, FNAME, args )       \
   BOOST_FORTRAN_PROTOTYPE(cxxname, fname, FNAME, 0, ret, 1, args )

#define BOOST_FORTRAN_FUNCTION_WITH_UNDERSCORE(ret, cxxname, fname, FNAME, args ) \
   BOOST_FORTRAN_PROTOTYPE(cxxname, fname, FNAME, 1, ret, 1, args )

#define BOOST_FORTRAN_SUBROUTINE(cxxname, fname, FNAME, args )          \
   BOOST_FORTRAN_PROTOTYPE(cxxname, fname, FNAME, 0, void, 0, args )

#define BOOST_FORTRAN_SUBROUTINE_WITH_UNDERSCORE(cxxname, fname, FNAME, args ) \
   BOOST_FORTRAN_PROTOTYPE(cxxname, fname, FNAME, 1, void, 0, args )



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#endif // BOOST_FORTRAN_PROTOTYPE_H_
