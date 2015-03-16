// UTypes.h
#pragma once

#include "UConstants.h"

// integer types

// 8-bit
typedef unsigned char UInt8;
typedef signed char SInt8;

// 16-bit
typedef unsigned short UInt16;
typedef short SInt16;

// 32-bit
typedef unsigned int UInt32;
typedef int SInt32;

// 64-bit
#if defined(_MSC_VER)
typedef unsigned __int64 UInt64;
typedef __int64 SInt64;
#else
typedef unsigned long long UInt64;
typedef long long SInt64;
#endif

