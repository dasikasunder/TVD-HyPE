/*
 * headers.hh
 *      Author: sunder
 */

#ifndef HEADERS_HH_
#define HEADERS_HH_

//----------------------------------------------------------------------------
// Commonly used C/C++ header files
//----------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <ctime>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <array>
#include <map>
#include <utility>

//----------------------------------------------------------------------------
// Preprocessors
//----------------------------------------------------------------------------

// Run in debug mode. The code runs slower but errors are easy to catch

#define DEBUG true

//----------------------------------------------------------------------------
// Error handling mechanism
//----------------------------------------------------------------------------

/* Enumerator defining whether to exit the program when an exception is caught
 * or to throw the error object using the C++ throw syntax
 */

enum ErrorHandling {
	Exit_No_Throw,
	Throw_No_Exit
};

//----------------------------------------------------------------------------

/* Base class for error handling */

class ErrorBase {
protected:
    const char* file;
    int line;
    const char* function;
    const char* cond;
    const char* exc;
public:
    ErrorBase();
    ErrorBase (const ErrorBase&);
    virtual ~ErrorBase() noexcept;

    void set_fields (const char*, const int, const char*, const char*);

    void print_log_data(std::ostream&) const;
    virtual void print_debug() const {std::cerr << "Generic Exception" << std::endl;}
};


//----------------------------------------------------------------------------

/* Class to check whether the subscript index is within proper range for
 * arrays
 */

class ErrIndexRange: public ErrorBase {
	int arg1, arg2, arg3;
public:
	ErrIndexRange(int, int, int);
	virtual ~ErrIndexRange() noexcept;
	virtual void print_debug() const {
		std::cerr << "The index " << arg1
				  << " is out of the range [" << arg2
				  << ", " << arg3
				  << ")" << std::endl;
	}
};

//----------------------------------------------------------------------------

/* Class to check if a number is zero */

class ErrZero: public ErrorBase {
public:
	ErrZero();
	virtual ~ErrZero() noexcept;
	virtual void print_debug() const {
		std::cerr << "A number is zero where it should not be" << std::endl;
	}
};



//----------------------------------------------------------------------------

/* Raise an error if an empty object is found, where it should not be */

class ErrEmptyObj: public ErrorBase {
public:
	ErrEmptyObj();
	virtual ~ErrEmptyObj() noexcept;
	virtual void print_debug() const {
		std::cerr << "You are trying to perform an operation on an empty object.\n"
				  << "It does not make any sense to perform such operation"
				  << std::endl;
	}
};

//----------------------------------------------------------------------------

/* Class to check whether two dimensions match */

class ErrDimMisMatch: public ErrorBase {
	int dim1, dim2;
public:
	ErrDimMisMatch(int, int);
	virtual ~ErrDimMisMatch() noexcept;
	virtual void print_debug() const {
		std::cerr << "The dimensions " << dim1 << " and " << dim2
				  << " do not match" << std::endl;
	}
};

//----------------------------------------------------------------------------

/* Raise an error if the pressure or density become zero in a riemann solver */

class ErrNegativePressureDensity: public ErrorBase {
	double d, p, x_loc, y_loc;
public:
	ErrNegativePressureDensity(double, double, double, double);
	virtual ~ErrNegativePressureDensity() noexcept;
	virtual void print_debug() const {
        std::cerr << "Negative Pressure/Density" << std::endl;
        std::cerr << "Density = " << d << ", Pressure = " << p << std::endl;
        std::cerr << "at (" << x_loc << ", " << y_loc << ")" << std::endl;
       }
};



//----------------------------------------------------------------------------

/* Raise an error if NaN or infinite number is encountered */

class ErrNotFinite: public ErrorBase {
public:
	ErrNotFinite();
	virtual ~ErrNotFinite() noexcept;
	virtual void print_debug() const {
		std::cerr << "A number that is not finite has been encountered in the program"
				  << std::endl;
	}
};



//----------------------------------------------------------------------------

/* Raise an error if the asked feature is not implemented */

class ErrNotImplemented: public ErrorBase {
public:
	ErrNotImplemented();
	virtual ~ErrNotImplemented() noexcept;
	virtual void print_debug() const {
		std::cerr << "The requested feature has not been implemented yet"
				  << std::endl;
	}
};

//----------------------------------------------------------------------------

/* Raise an error and print some message */

class ErrMessage: public ErrorBase {
	std::string msg;
public:
	ErrMessage(std::string);
	virtual ~ErrMessage() noexcept;
	virtual void print_debug() const {
		std::cerr << msg << std::endl;
	}
};



//----------------------------------------------------------------------------

/* Unable to open a file for reading or writing data */

class ErrIO: public ErrorBase {
	std::string filename;
public:
	ErrIO(std::string);
	virtual ~ErrIO() noexcept;
	virtual void print_debug() const {
		std::cerr <<  "Unable to open the file: " << filename << " for reading or writing data." << std::endl;
	}
};



//----------------------------------------------------------------------------

/* Singular matrix detected */

class ErrSingular: public ErrorBase {
public:
	ErrSingular();
	virtual ~ErrSingular() noexcept;
	virtual void print_debug() const {
		std::cerr <<  "Singular matrix encountered" << std::endl;
	}
};

//----------------------------------------------------------------------------

/* Templated function for setting various fields of the error object. Depending on
 * the ErrorHandling enumerator, either the program would be aborted or an exception
 * will be thrown
 */

namespace Internal {
	template <class exc>
	void issue_error (ErrorHandling handler, const char *file, int line, const char *function, const char *cond, exc e)
	{
		 // Fill the fields of the exception object
		e.set_fields (file, line, function, cond);

		switch (handler) {
		case Exit_No_Throw:
			e.print_log_data(std::cerr);
			e.print_debug();
			std::exit(EXIT_FAILURE);
			break;
		case Throw_No_Exit:
			throw(e);
			break;
		}
	}
}

/* C-style assert macro. If condition is false, useful information
 * about the condition that failed along with the file, line and
 * additional information is printed. The program will then exit immediately.
 */

#if DEBUG == true
#define Assert(cond, err) if (!(cond)) {\
	Internal::issue_error(Exit_No_Throw, __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, err);\
}
#else
#define Assert(cond, err)
#endif

/* C++-style assert macro. If condition is false, an exception of any of the
 * derived class of ErrorBase will be thrown. No information will be printed.
 * Should the user require, the necessary information can be printed using the
 * catch statement.
 */

#define AssertThrow(cond, err) if (!(cond)) {\
	Internal::issue_error(Throw_No_Exit, __FILE__, __LINE__, __PRETTY_FUNCTION__, #cond, err);\
}

//----------------------------------------------------------------------------
// Utilities namespace. Simple random functions used frequently
//----------------------------------------------------------------------------

namespace Utilities {

	// SIGN function. returns a with the sign of b

	template<class T>
	inline T SIGN(const T &a, const T &b) {
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}

	// SGN function. returns a with the sign of b

	template<class T>
	inline T sgn(const T &a) {
		return a >= 0 ? 1.0 : -1.0;
	}

	// Check if two floating point numbers are close

	template<class T>
	bool isclose(const T& a, const T& b, const double rtol=1.0e-05,
			const double atol=1.0e-08) {

		if (std::isnan(a) || std::isnan(b) )
			return false;

		else
			return std::abs(a-b) <= (atol + rtol * std::abs(b)) ? true : false;

	}
	
	// Allocate and free memory for dynamic arrays

	template<typename T>
	void allocate_mem_1d(T*&, size_t);

	template<typename T>
	void allocate_mem_2d(T**&, size_t, size_t);

	template<typename T>
	void allocate_mem_3d(T***&, size_t, size_t, size_t);

	template<typename T>
	void allocate_mem_4d(T****&, size_t, size_t, size_t, size_t);

	template<typename T>
	void free_mem_1d(T*&);

	template<typename T>
	void free_mem_2d(T**&, size_t);

	template<typename T>
	void free_mem_3d(T**&, size_t, size_t);

	template<typename T>
	void free_mem_4d(T****&, size_t, size_t, size_t);

	// Allocate memory for 1D array

	template<typename T>
	void allocate_mem_1d(T*& U, size_t m) {
		U = new T[m];
	}


	// Allocate memory for 2D array

	template<typename T>
	void allocate_mem_2d(T**& U, size_t m, size_t n) {
		U = new T*[m];

		for (unsigned int i = 0; i < m; ++i) {
			U[i] = new T[n];
		}
	}

	// Allocate memory for 3D array

	template<typename T>
	void allocate_mem_3d(T***& U, size_t m, size_t n, size_t o) {
		U = new T**[m];

		for(unsigned int i = 0 ; i < m; ++i) {
			U[i] = new T*[n];
		    for(unsigned int j = 0; j < n; ++j){
		    	U[i][j] = new T[o];
		    }
		}
	}

	// Allocate memory for 4D array

	template<typename T>
	void allocate_mem_4d(T****& U, size_t size1, size_t size2, size_t size3, size_t size4) {

		U = new T***[size1];

		for(size_t i = 0 ; i < size1; ++i) {

			U[i] = new T**[size2];

			for(size_t j = 0; j < size2; ++j){

				U[i][j] = new T*[size3];

				for (size_t k = 0; k < size3; ++k) {

					U[i][j][k] = new T[size4];

				}

			}
		}
	}

	// Release memory for 1D array

	template<typename T>
	void free_mem_1d(T*& U) {
		delete[] U;
	}

	// Release memory for 2D array

	template<typename T>
	void free_mem_2d(T**& U, size_t m) {

		for (unsigned int i = 0; i < m; ++i) {
			delete[] U[i];
		}

		delete[] U;
	}

	// Release memory for 3D array

	template<typename T>
	void free_mem_3d(T***& U, size_t m, size_t n) {
	    for (unsigned int i = 0; i < m; ++i) {
	        for (unsigned int j = 0; j < n; ++j) {
	            delete [] U[i][j];
	        }
	        delete [] U[i];
	    }
	    delete [] U;
	}

	// Release memory for 4D array

	template<typename T>
	void free_mem_4d(T****& U, size_t size1, size_t size2, size_t size3) {

		for (size_t i = 0; i < size1; ++i) {
			for (size_t j = 0; j < size2; ++j) {
				for (unsigned int k = 0; k < size3; ++k) {
					delete[] U[i][j][k];
				}
				delete[] U[i][j];
			}
			delete[] U[i];
		}

		delete[] U;
	}

	std::string int_to_string (unsigned int value, const unsigned int digits);
}

#endif /* HEADERS_HH_ */

