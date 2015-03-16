// UArray.h
#pragma once

#include "../core/Utility.h"
#include <memory.h>
#include <string.h>

template<class T>
class UArray
{
public:
	typedef T Element;

// Construction
	UArray();
	~UArray();

	UArray(const UArray<T> & other);
	UArray( UInt32 initSize );

	UArray<T> operator =(const UArray<T> & other);

// Attributes
	UInt32 getSize() const;
	SInt32 getUpperBound() const;
	void setSize(UInt32 newSize, SInt32 growBy = -1);

// Operations
	// Clean up
	void freeExtra();
	void removeAll();

	// Accessing elements
	T getAt(UInt32 index) const;
	void setAt(UInt32 index, const T& newElement);
	T & elementAt(UInt32 index);

	// Direct Access to the element data (may return NULL)
	const T* getData() const;
	T* getData();

	// Potentially growing the array
	void setAtAndGrow(UInt32 index, const T& newElement);
	UInt32 add( const T& newElement);
	UInt32 appendArray(const UArray& source);
	void copy(const UArray& source);

	// overloaded operator helpers
	const T& operator[](UInt32 index) const;
	T& operator[](UInt32 index);

	// Operations that move elements around
	void insertAt(UInt32 index, const T& newElement, UInt32 count = 1);
	void removeAt(UInt32 index, UInt32 count = 1);
	void insertAt(UInt32 startIndex, UArray * newArray);

	// Sorting
	void quickSort();
	void quickSort( UInt32 lowIndex, UInt32 highIndex );
	void insertionSort();
	void insertionSort( UInt32 lowIndex, UInt32 highIndex );

	// removes duplicates
	void uniq();

// Implementation
protected:
	T * data;   // the actual array of data
	UInt32 size;     // # of elements (upperBound - 1)
	UInt32 maxSize;  // max allocated
	UInt32 growBy;   // grow amount
};

/////////////////////////////////////////////////////////////////////////////
// UArray<T> inline functions

template<class T>
u_inline UInt32 UArray<T>::getSize() const
{
	return size;
}

template<class T>
u_inline SInt32 UArray<T>::getUpperBound() const
{
	return size-1;
}

template<class T>
u_inline void UArray<T>::removeAll()
{
	setSize(0, -1);
}

template<class T>
u_inline T UArray<T>::getAt(UInt32 index) const
{
	U_ASSERT(index < size);
	return data[index];
}

template<class T>
u_inline void UArray<T>::setAt(UInt32 index, const T& newElement)
{
	U_ASSERT( index < size);
	data[index] = newElement;
}

template<class T>
u_inline T& UArray<T>::elementAt(UInt32 index)
{
	U_ASSERT(index < size);
	return data[index];
}

template<class T>
u_inline const T* UArray<T>::getData() const
{
	return (const T*)data;
}

template<class T>
u_inline T* UArray<T>::getData()
{
	return (T*)data;
}

template<class T>
u_inline UInt32 UArray<T>::add( const T& newElement)
{
	UInt32 index = size;
	setAtAndGrow(index, newElement);
	return index;
}

template<class T>
u_inline const T& UArray<T>::operator[](UInt32 index) const
{
	U_ASSERT(index < size);
	return data[index];
}

template<class T>
u_inline T& UArray<T>::operator[](UInt32 index)
{
	U_ASSERT(index < size);
	return data[index];
}

/////////////////////////////////////////////////////////////////////////////
// UArray<T> out-of-line functions

template<class T>
UArray<T>::UArray()
{
	data = NULL;
	size = maxSize = growBy = 0;
}

template<class T>
UArray<T>::~UArray()
{
	if (data != NULL)
	{
		UDestructElements<T>(data, size);
		delete[] (UInt8*)data;
	}
}

template<class T>
UArray<T>::UArray(const UArray<T> & other)
{
	data = NULL;
	size = maxSize = growBy = 0;
	copy(other);
}

template<class T>
UArray<T>::UArray( UInt32 initSize )
{
	data = NULL;
	size = maxSize = growBy = 0;
	setSize( initSize );
}

template<class T>
UArray<T> UArray<T>::operator =(const UArray<T> & other)
{
	removeAll();
	copy(other);
	return *this;
}

template<class T>
void UArray<T>::setSize(UInt32 newSize, SInt32 newGrowBy)
{
	if (newGrowBy != -1)
		growBy = newGrowBy;

	if (newSize == 0)
	{
		if (data != NULL)
		{
			UDestructElements<T>(data, size);
			delete[] (UInt8*)data;
			data = NULL;
		}
		size = maxSize = 0;
	}
	else if (data == NULL)
	{
		data = (T*) new UInt8[newSize * sizeof(T)];
		UConstructElements<T>(data, newSize);
		size = maxSize = newSize;
	}
	else if (newSize <= maxSize)
	{
		if (newSize > size)
		{
			UConstructElements<T>(&data[size], newSize-size);
		}
		else if (size > newSize)
		{
			UDestructElements<T>(&data[newSize], size-newSize);
		}
		size = newSize;
	}
	else
	{
		UInt32 newGrowBy = growBy;
		if (newGrowBy == 0)
		{
			newGrowBy = size / 8;
			newGrowBy = (newGrowBy < 4) ? 4 : ((newGrowBy > 1024) ? 1024 : newGrowBy);
		}
		UInt32 newMax;
		if (newSize < maxSize + newGrowBy)
			newMax = maxSize + newGrowBy;
		else
			newMax = newSize;

		U_ASSERT(newMax >= maxSize);
		T* newData = (T*) new UInt8[newMax * sizeof(T)];

		memcpy(newData, data, size * sizeof(T));

		U_ASSERT(newSize > size);
		UConstructElements<T>(&newData[size], newSize-size);

		delete[] (UInt8*)data;
		data = newData;
		size = newSize;
		maxSize = newMax;
	}
}

template<class T>
UInt32 UArray<T>::appendArray(const UArray& source)
{
	U_ASSERT(this != &source);

	UInt32 oldSize = size;
	setSize(size + source.size);
	UCopyElements<T>(data + oldSize, source.data, source.size);
	return oldSize;
}

template<class T>
void UArray<T>::copy(const UArray& source)
{
	U_ASSERT(this != &source);

	setSize(source.size);
	UCopyElements<T>(data, source.data, source.size);
}

template<class T>
void UArray<T>::freeExtra()
{
	if (size != maxSize)
	{
		T* newData = NULL;
		if (size != 0)
		{
			newData = (T*) new UInt8[size * sizeof(T)];
			memcpy(newData, data, size * sizeof(T));
		}

		delete[] (UInt8*)data;
		data = newData;
		maxSize = size;
	}
}

template<class T>
void UArray<T>::setAtAndGrow(UInt32 index, const T& newElement)
{
	if (index >= size)
		setSize(index+1, -1);
	data[index] = newElement;
}

template<class T>
void UArray<T>::insertAt(UInt32 index, const T& newElement, UInt32 count)
{
	if(index >= size)
	{
		setSize(index + count, -1);
	}
	else
	{
		UInt32 oldSize = size;
		setSize(size + count, -1);
		UDestructElements<T>(&data[oldSize], count);
		memmove(&data[index+count], &data[index],
			(oldSize-index) * sizeof(T));

		UConstructElements<T>(&data[index], count);
	}

	U_ASSERT(index + count <= size);
	while (count--)
		data[index++] = newElement;
}

template<class T>
void UArray<T>::removeAt(UInt32 index, UInt32 count)
{
	U_ASSERT(index + count <= size);

	SInt32 moveCount = size - (index + count);
	UDestructElements<T>(&data[index], count);
	if (moveCount)
		memmove(&data[index], &data[index + count], moveCount * sizeof(T));
	size -= count;
}

template<class T>
void UArray<T>::insertAt(UInt32 startIndex, UArray* newArray)
{
	U_ASSERT((newArray != NULL) && (newArray != this));

	if (newArray->getSize() > 0)
	{
		insertAt(startIndex, newArray->getAt(0), newArray->getSize());
		for (UInt32 i = 0; i < newArray->getSize(); i++)
			setAt(startIndex + i, newArray->getAt(i));
	}
}


template <class T>
u_inline void UArray<T>::quickSort()
{
	if( getSize() != 0 )
		quickSort(0, getUpperBound());
}

// Internal quicksort method that makes recursive calls.
// Uses median-of-three partitioning and a cutoff of 10.
// lowIndex is the left-most index of the subarray.
// highIndex is the right-most index of the subarray.
template <class T>
void UArray<T>::quickSort( UInt32 lowIndex, UInt32 highIndex )
{
	T* pivot;
    if( lowIndex + 10 > highIndex )
        insertionSort( lowIndex, highIndex );
    else
    {
          // Sort lowIndex, middleIndex, highIndex
        UInt32 middleIndex = ( lowIndex + highIndex ) / 2;
        if( data[middleIndex] < data[lowIndex] )
            USwapElements( data[lowIndex], data[middleIndex] );
        if( data[highIndex] < data[lowIndex] )
            USwapElements( data[lowIndex], data[highIndex] );
        if( data[highIndex] < data[middleIndex] )
            USwapElements( data[middleIndex], data[highIndex] );

          // Place pivot at position highIndex - 1
 //       T pivot = data[middleIndex];
        USwapElements( data[middleIndex], data[highIndex - 1] );
		pivot = &data[highIndex-1];

          // Begin partitioning
        SInt32 i = lowIndex;
		SInt32 j = highIndex - 1;
        while(true)
        {
            while( data[++i] < *pivot )
				;
            while( *pivot < data[--j] )
				;
            if( i < j )
                USwapElements( data[i], data[j] );
            else
                break;
        }
        USwapElements( data[i], data[highIndex - 1] );  // Restore pivot

        quickSort( lowIndex, i - 1 );     // Sort small elements
        quickSort( i + 1, highIndex );    // Sort large elements
    }
}

template <class T>
u_inline void UArray<T>::insertionSort()
{
	if( getSize() != 0 )
		insertionSort(0, getUpperBound());
}

// Internal insertion sort routine for subarrays
// that is used by quicksort.
// lowIndex is the left-most index of the subarray.
// highIndex is the right-most index of the subarray.
template <class T>
void UArray<T>::insertionSort( UInt32 lowIndex, UInt32 highIndex )
{
	struct Data
	{
		UInt8 buffer[ sizeof(T) ];
	};

	if( highIndex == lowIndex )
		return;

	Data tempData;
	T* temp;
	for( UInt32 i = lowIndex; i <= highIndex; i++ )
	{
//		UInt8 buffer[ sizeof(T) ];
//		memcpy( buffer, &data[i], sizeof(T) );
//		T& temp = *((T*)buffer);
//      T temp = data[i];

		tempData = *((Data*) &(data[i]));
		temp = ((T*)&tempData);

		UInt32 j;
		for( j = i; (j > lowIndex) && (*temp < data[j - 1]); j-- )
		{
//			memcpy( &data[j], &data[j-1], sizeof(T) );
//            data[j] = data[j - 1];
			*((Data*)&(data[j])) = *((Data*)&(data[j-1]));
		}
//		memcpy( &data[j], &buffer, sizeof(T) );
//		data[j] = temp;
		*((Data*)&(data[j])) = tempData;
	}
}

template <class T>
void UArray<T>::uniq()
{
	if( !getSize() ) return;
	UInt32 i, j;
	UArray<T> tmp = *this;
	tmp.quickSort();
	for( i = 0; i < tmp.getSize() - 1; i++ )
	{
		if( tmp[i] == tmp[i+1] )
		{
			tmp.removeAt(i+1);
			for( j = this->getUpperBound(); j >= 0; j-- )
			{
				if( (*this)[j] == tmp[i] )
				{
					this->removeAt( j );
					break;
				}
			}
			U_ASSERT( j != 0 ); // this means the function is incorrect
			i--;
		}
	}
}

