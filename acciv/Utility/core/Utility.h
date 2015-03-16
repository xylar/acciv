// Utility.h
#pragma once

#include "UTypes.h"
#include <new>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _DEBUG
#define U_ASSERT(condition) logMsg(condition, #condition, __FILE__, __LINE__);
#else
#define U_ASSERT(condition)
#endif

#define u_inline __inline

u_inline bool logMsg(bool x, const char *msg, char* file, unsigned int line)
{
  if( false == x)
  {
    //do that extra logging information to a file etc.
    printf("On line %i in file %s : Error !! Assert %s  failed\n", line, file, msg);
    abort();
    return (true);
  }
  else
  {
    return(true);
  }
}

template<class T>
u_inline void UConstructElements(T * elements, UInt32 count)
{
	// first do bit-wise zero initialization
	memset((void*)elements, 0, count * sizeof(T));

	// then call the constructor(s)
	for (; count--; elements++)
		::new((void *)elements) T;
}

template<class T>
u_inline void UDestructElements(T * elements, UInt32 count)
{
	// call the destructor(s)
	for (; count--; elements++)
		elements->~T();
}

template<class T>
u_inline void UCopyElements(T * destination, const T * source, UInt32 count)
{
	while(count--)
		*destination++ = *source++;
}

template <class T>
u_inline void USwapElements( T & object1, T & object2 )
{
	struct Data
	{
		UInt8 buffer[ sizeof(T) ];
	};

	Data temp;

	temp = *((Data*)&object1);
	*((Data*)&object1) = *((Data*)&object2);
	*((Data*)&object2) = *((Data*)&temp);
}

template<class T>
u_inline void UStaticCopy( T& destination, const T& source )
{
	struct Data
	{
		UInt8 buffer[ sizeof(T) ];
	};

	*((Data*)&destination) = *((Data*)&source);
}
