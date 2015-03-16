// UString.h
#pragma once

#include "Utility.h"
#include <stdio.h>
#include <stdarg.h>

class UString  
{
public:
	UString();
//	UString( const std::string inString );
	UString( bool inBoolean );
	UString( const UString& inString );
	UString( const char* inString );
	UString( char inCharacter );
	UString( int inInt );
	UString( float inFloat );

	~UString();
	
	UInt32 getLength() const;
	operator const char* () const;

	char getCharacterAt( UInt32 inIndex ) const;
	void setCharacterAt( UInt32 inIndex, char inCharacter );

	UString operator=( const UString& inString );
	UString operator=( const char* inString );

	UString operator+=( const UString& inString );
	UString operator+=( const char* inString ); 

	SInt32 compare( const char* inString ) const;
	SInt32 compareNoCase( const char* inString ) const;

	UString mid( UInt32 inFirstIndex ) const;
	UString mid( UInt32 inFirstIndex, UInt32 inCount ) const;

	UString left( UInt32 inCount ) const;
	UString right( UInt32 inCount ) const;

	UString spanIncluding( const char* inCharacterSet ) const;
	UString spanExcluding( const char* inCharacterSet ) const;

	void makeUpper();
	void makeLower();
	void makeReverse();

	UInt32 replace( char inCharacterFrom, char inCharacterTo );
	UInt32 replace( const char* inStringFrom, const char* inStringTo );

	UInt32 remove( char inCharacter );
	UInt32 remove( const char* inString );

	void remove( UInt32 inFirstIndex, UInt32 inCount = 1 );

	void insert( UInt32 inIndex, char inCharacter );
	void insert( UInt32 inIndex, const char* inString );

	static UString makeFormat( const char* inFormatString, ... );
	static UString makeFormatV( const char* inFormatString, va_list argList);
	void format( const char* inFormatString, ... );
	void formatV( const char* inFormatString, va_list argList);

	void trimLeft();
	void trimLeft( char inCharacter );
	void trimLeft( const char* inCharacterSet );

	void trimRight();
	void trimRight( char inCharacter );
	void trimRight( const char* inCharacterSet );

	SInt32 find( char inCharacter, SInt32 inFirstIndex = 0 ) const;
	SInt32 find( const char* inString, SInt32 inFirstIndex = 0 ) const;

	SInt32 reverseFind( char inCharacter ) const;
	SInt32 findOneOf( const char* inCharacterSet, SInt32 inFirstIndex = 0 ) const;

	void split2( UString& dst1, UString& dst2, char separator ) const;
	void split2( UString& dst1, UString& dst2 ) const;

	char* getBuffer( UInt32 inMinimumLength );
	char* getBufferOfLength( UInt32 inLength );
	char* getBuffer() { return getBufferOfLength( getLength() ); }

	void releaseBuffer( SInt32 inNewLength = -1 );

	void freeExtra();

private:
	struct Data
	{
		UInt32 stringLength;
		UInt32 allocatedLength;
		SInt32 referenceCount;
	};

	void setData( Data* inData );

	UInt32 getDataLength() const;
	UInt32 getDataAllocatedLength() const;

	u_inline char* getDataBuffer();
	u_inline const char * getDataBuffer() const;
	u_inline UString::Data* getData() const;

	void setDataLength(UInt32 inDataLength);

	void assureDataModifiable();
	void assureDataAvailable( UInt32 inMinimumLength );

	void reallocateData( UInt32 inStringLength );

	static void acquireData( Data* inData );
	static void releaseData( Data* inData );

	char* data;

	static Data * const kEmptyStringData;
};

u_inline UInt32 UString::getLength() const
{
	return getDataLength();
}

u_inline UString::operator const char* () const
{
	return getDataBuffer();
}

u_inline UInt32 UString::getDataLength() const
{
	return getData()->stringLength;
}

u_inline UInt32 UString::getDataAllocatedLength() const
{
	return getData()->allocatedLength;
}

u_inline char* UString::getDataBuffer()
{
	return data;
}

u_inline const char * UString::getDataBuffer() const
{
	return data;
}

u_inline UString::Data* UString::getData() const
{
	return ((Data*)data) - 1;
}

u_inline UString operator+( const UString& inLeft, const UString& inRight )
{
	UString result = inLeft;
	result += inRight;
	return result;
}

u_inline UString operator+( const UString& inLeft, const char* inRight )
{
	UString result = inLeft;
	result += inRight;
	return result;
}

u_inline UString operator+( const char* inLeft, const UString& inRight )
{
	UString result = inLeft;
	result += inRight;
	return result;
}

u_inline bool operator <( const UString& inLeft, const UString& inRight )
{
	return (inLeft.compare(inRight) < 0);
}

u_inline bool operator <( const UString& inLeft, const char* inRight )
{
	return (inLeft.compare(inRight) < 0);
}

u_inline bool operator <( const char* inLeft, const UString& inRight )
{
	return (inRight.compare(inLeft) >= 0);
}

u_inline bool operator >( const UString& inLeft, const UString& inRight )
{
	return (inLeft.compare(inRight) > 0);
}

u_inline bool operator >( const UString& inLeft, const char* inRight )
{
	return (inLeft.compare(inRight) > 0);
}

u_inline bool operator >( const char* inLeft, const UString& inRight )
{
	return (inRight.compare(inLeft) <= 0);
}

u_inline bool operator <=( const UString& inLeft, const UString& inRight )
{
	return (inLeft.compare(inRight) <= 0);
}

u_inline bool operator <=( const UString& inLeft, const char* inRight )
{
	return (inLeft.compare(inRight) <= 0);
}

u_inline bool operator <=( const char* inLeft, const UString& inRight )
{
	return (inRight.compare(inLeft) > 0);
}

u_inline bool operator >=( const UString& inLeft, const UString& inRight )
{
	return (inLeft.compare(inRight) >= 0);
}

u_inline bool operator >=( const UString& inLeft, const char* inRight )
{
	return (inLeft.compare(inRight) >= 0);
}

u_inline bool operator >=( const char* inLeft, const UString& inRight )
{
	return (inRight.compare(inLeft) < 0);
}

u_inline bool operator ==( const UString& inLeft, const UString& inRight )
{
	return (inLeft.compare(inRight) == 0);
}

u_inline bool operator ==( const UString& inLeft, const char* inRight )
{
	return (inLeft.compare(inRight) == 0);
}

u_inline bool operator ==( const char* inLeft, const UString& inRight )
{
	return (inRight.compare(inLeft) == 0);
}

u_inline bool operator !=( const UString& inLeft, const UString& inRight )
{
	return (inLeft.compare(inRight) != 0);
}

u_inline bool operator !=( const UString& inLeft, const char* inRight )
{
	return (inLeft.compare(inRight) != 0);
}

u_inline bool operator !=( const char* inLeft, const UString& inRight )
{
	return (inRight.compare(inLeft) != 0);
}

