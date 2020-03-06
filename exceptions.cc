/*
 * exceptions.cc
 *      Author: sunder
 */

// Contains exception handling classes

#include "headers.hh"

//----------------------------------------------------------------------------
// Base error class
//----------------------------------------------------------------------------

// Default constructor

ErrorBase::ErrorBase() :
    file(""),
    line(0),
    function(""),
    cond(""),
    exc("")
{}

// Copy constructor

ErrorBase::ErrorBase(const ErrorBase &err) :
    file(err.file),
    line(err.line),
    function(err.function),
    cond(err.cond),
    exc(err.exc)
{}

// Destructor

ErrorBase::~ErrorBase() noexcept
{ // Does nothing
}

// Set various entries of the ErrorBase class

void ErrorBase::set_fields (const char *f,
                                const int  l,
                                const char *func,
                                const char *c)
{
	file = f;
    line = l;
    function = func;
    cond = c;
}

// Print general information about the error, this includes name of the file,
// line number, violated condition and the function


void ErrorBase::print_log_data (std::ostream& out) const
{
	// print a header for the error
    out << "An error occurred in line <" << line
        << "> of file <" << file
        << "> in function" << std::endl
        << "    " << function << std::endl
        << "The violated condition was: "<< std::endl
        << "    " << cond << std::endl;
    // finally print the additional information the error provides:
    out << "Additional information: " << std::endl;
}


//----------------------------------------------------------------------------
// Check array indices
//----------------------------------------------------------------------------

// Constructor

ErrIndexRange::ErrIndexRange(int a1, int a2, int a3) :
		arg1(a1),
		arg2(a2),
		arg3(a3)
{}

// Destructor

ErrIndexRange::~ErrIndexRange() noexcept
{ // Does nothing
}

//----------------------------------------------------------------------------
// Check if a number is zero
//----------------------------------------------------------------------------

// Constructor

ErrZero::ErrZero()
{}

// Destructor

ErrZero::~ErrZero() noexcept
{ // Does nothing
}

//----------------------------------------------------------------------------
// Check if an empty object is found
//----------------------------------------------------------------------------

// Constructor

ErrEmptyObj::ErrEmptyObj()
{}

// Destructor

ErrEmptyObj::~ErrEmptyObj() noexcept
{ // Does nothing
}


//----------------------------------------------------------------------------
// Class to check whether two dimensions match
//----------------------------------------------------------------------------

// Constructor

ErrDimMisMatch::ErrDimMisMatch(int d1, int d2) :
		dim1(d1),
		dim2(d2)
{}

// Destructor

ErrDimMisMatch::~ErrDimMisMatch() noexcept
{ // Does nothing
}

//----------------------------------------------------------------------------
// Class if pressure or density become negative
//----------------------------------------------------------------------------

// Constructor

ErrNegativePressureDensity::ErrNegativePressureDensity(double d_, double p_,
		double x_, double y_):

		d(d_),
		p(p_),
		x_loc(x_),
		y_loc(y_)
{}

// Destructor

ErrNegativePressureDensity::~ErrNegativePressureDensity() noexcept
{ // Does nothing
}


//----------------------------------------------------------------------------
// NaN is encountered
//---------------------------------------------------------------------------

// Constructor

ErrNotFinite::ErrNotFinite()
{}

// Destructor

ErrNotFinite::~ErrNotFinite() noexcept
{ // Does nothing
}

//----------------------------------------------------------------------------
// Requested feature is not implemented
//----------------------------------------------------------------------------

// Constructor

ErrNotImplemented::ErrNotImplemented()
{}

// Destructor

ErrNotImplemented::~ErrNotImplemented() noexcept
{ // Does nothing
}

//----------------------------------------------------------------------------
// Print some error message
//----------------------------------------------------------------------------

// Constructor

ErrMessage::ErrMessage(std::string message) :
		msg(message)
{}

// Destructor

ErrMessage::~ErrMessage() noexcept
{ // Does nothing
}

//----------------------------------------------------------------------------
// Unable to open file for reading/writing data
//----------------------------------------------------------------------------

// Constructor

ErrIO::ErrIO(std::string f) :
		filename(f)
{}


// Destructor

ErrIO::~ErrIO() noexcept
{ // Does nothing
}

//----------------------------------------------------------------------------
// Singular matrix encountered
//----------------------------------------------------------------------------

// Constructor

ErrSingular::ErrSingular() {
}

// Destructor

ErrSingular::~ErrSingular() noexcept
{ // Does nothing
}
