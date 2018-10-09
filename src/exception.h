/******************************************************************************
 *                       SMART EXCEPTIONS                                     *
 ******************************************************************************/

#ifndef SMARTASTRO_EXCEPTIONS_H
#define SMARTASTRO_EXCEPTIONS_H

#include <exception>
#include <cassert>
#include <iostream>
#include <string>

#define _SMARTASTRO_EXCEPTION_QUOTEME(x) #x
#define SMARTASTRO_EXCEPTION_QUOTEME(x) _SMARTASTRO_EXCEPTION_QUOTEME(x)
#define SMARTASTRO_EXCEPTION_EXCTOR(s) ((std::string(__FILE__ "," SMARTASTRO_EXCEPTION_QUOTEME(__LINE__) ": ") + s) + ".")
#define SMARTASTRO_EX_THROW(s) (throw smartastro_exception(SMARTASTRO_EXCEPTION_EXCTOR(s)))

#define smartastro_throw(s) SMARTASTRO_EX_THROW(s)

namespace smartastro{
class smartastro_exception: public std::exception {
	public:
		smartastro_exception(const std::string &s):m_what(s) {}
			virtual const char *what() const throw() {	
				return m_what.c_str();
		}
	virtual ~smartastro_exception() throw() {}
	protected:
		std::string m_what;
};
}

#endif
