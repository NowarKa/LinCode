#pragma once

typedef signed long long int Fint;

typedef enum
{
   OK,                      // No error
   ErrorIncompatibleFields, // The two number or element are from two different fields.
   ErrorNoInverse,          // The number of element does not have a multiplicative inverse.
   ErrorDivideByZero,       // Divide by zero happened.
   ErrorInvalidArraySize    // One or more of arrays in function argument have invalid size.
} ErrorCode;
