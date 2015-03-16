// UString.cpp
#include "UString.h"

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

UInt32 gUString_EmptyStringDataObject[4] = { 0, 4, 1, 0 };
UString::Data * const UString::kEmptyStringData = (UString::Data*)(&gUString_EmptyStringDataObject);


UString::UString()
{
	acquireData( kEmptyStringData );
	data = (char*)(kEmptyStringData+1);
}

UString::UString( int inInt )
{
	acquireData( kEmptyStringData );
	data = (char*)(kEmptyStringData+1);

	format( "%d", inInt );
}

UString::UString( float inFloat )
{
	acquireData( kEmptyStringData );
	data = (char*)(kEmptyStringData+1);

	format( "%f", inFloat );
}

UString::~UString()
{
	releaseData( getData() );
}

UString::UString( bool inBoolean )
{
	acquireData( kEmptyStringData );
	data = (char*)(kEmptyStringData+1);

	if( inBoolean )
		*this = "true";
	else
		*this = "false";
}

//UString::UString( const std::string inString )
//{
//	U_ASSERT( inString != NULL );
//	acquireData( kEmptyStringData );
//	data = (char*)(kEmptyStringData+1);
//	*this = inString.c_str(); // Pany: Why is it telling me that this wants to be recursive or some crap?
//}

UString::UString( const UString& inString )
{
	acquireData( inString.getData() );
	data = inString.data;
}

UString::UString( const char* inString )
{
	U_ASSERT( inString != NULL );
	acquireData( kEmptyStringData );
	data = (char*)(kEmptyStringData+1);
	*this = inString;
}

UString::UString( char inCharacter )
{
	acquireData( kEmptyStringData );
	data = (char*)(kEmptyStringData+1);
	assureDataAvailable(1);
	assureDataModifiable();
	getDataBuffer()[0] = inCharacter;
	setDataLength( 1 );
}

char UString::getCharacterAt( UInt32 inIndex ) const
{
	U_ASSERT(inIndex < getLength());
	return getDataBuffer()[inIndex];
}

void UString::setCharacterAt( UInt32 inIndex, char inCharacter )
{
	U_ASSERT(inIndex < getLength());
	assureDataModifiable();
	getDataBuffer()[inIndex] = inCharacter;
	if( inCharacter == '\0' )
		setDataLength( inIndex );
}

UString UString::operator=( const UString& inString )
{
	setData( inString.getData() );
	return *this;
}

UString UString::operator=( const char* inString )
{
	UInt32 length = strlen( inString );
	
	assureDataAvailable( length );
	assureDataModifiable();

	memcpy( getDataBuffer(), inString, length );
	setDataLength( length );

	return *this;
}

UString UString::operator+=( const UString& inString )
{
	UInt32 length = getLength();
	UInt32 otherLength = inString.getLength();
	assureDataAvailable( length + otherLength );
	assureDataModifiable();

	memcpy( getDataBuffer() + length, inString.getDataBuffer(), otherLength );
	setDataLength( length + otherLength );
	return *this;
}

UString UString::operator+=( const char* inString )
{
	UInt32 length = getLength();
	UInt32 otherLength = strlen(inString);
	assureDataAvailable( length + otherLength );
	assureDataModifiable();
	
	memcpy( getDataBuffer() + length, inString, otherLength );
	setDataLength( length + otherLength );
	return *this;
}

SInt32 UString::compare( const char* inString ) const
{
	return strcmp(getDataBuffer(), inString);
}

SInt32 UString::compareNoCase( const char* inString ) const
{
	const char * thisString = getDataBuffer();
	const char * otherString = inString;
	SInt32 result = 0;
	do
   	{
		result = toupper(*thisString) - toupper(*otherString);
	} while (*(thisString++) && *(otherString++) && !result);
	return result;
}

UString UString::mid( UInt32 inFirstIndex ) const
{
	U_ASSERT(inFirstIndex <= getLength());
	return UString(getDataBuffer() + inFirstIndex);
}

UString UString::mid( UInt32 inFirstIndex, UInt32 inCount ) const
{
	U_ASSERT(inFirstIndex + inCount <= getLength());
	if((inFirstIndex == 0) && (inCount == getLength()))
		return *this;
	UString result;
	result.assureDataAvailable(inCount);
	result.assureDataModifiable();
	memcpy(result.getDataBuffer(), getDataBuffer() + inFirstIndex, inCount);
	result.setDataLength(inCount);
	return result;
}

UString UString::left( UInt32 inCount ) const
{
	return mid(0, inCount);
}

UString UString::right( UInt32 inCount ) const
{
	return mid(getLength() - inCount, inCount);
}

UString UString::spanIncluding( const char* inCharacterSet ) const
{
	return left(strspn(getDataBuffer(), inCharacterSet));
}

UString UString::spanExcluding( const char* inCharacterSet ) const
{
	return left(strcspn(getDataBuffer(), inCharacterSet));
}

void UString::makeUpper()
{
	assureDataModifiable();
	char * thisString = getDataBuffer();
	while(*thisString)
	{
		*thisString = toupper(*thisString);
		thisString++;
	}
}

void UString::makeLower()
{
	assureDataModifiable();
	char * thisString = getDataBuffer();
	while(*thisString)
	{
		*thisString = tolower(*thisString);
		thisString++;
	}
}

void UString::makeReverse()
{
	assureDataModifiable();
	char * string1 = getDataBuffer();
	char * string2 = string1+getLength()-1;
	while(string1 < string2)
	{
		char temp = *string1;
		*string1 = *string2;
		*string2 = temp;
		string1++;
		string2--;
	}
}

UInt32 UString::replace( char inCharacterFrom, char inCharacterTo )
{
	UInt32 count = 0;
	if(inCharacterFrom != inCharacterTo)
	{
		assureDataModifiable();
		char * buffer = getDataBuffer();
		char * end = buffer + getLength();
		while(buffer < end)
		{
			if(*buffer == inCharacterFrom)
			{
				*buffer = inCharacterTo;
				count++;
			}
			buffer++;
		}
	}
	return count;
}


UInt32 UString::replace( const char* inStringFrom, const char* inStringTo )
{
	U_ASSERT((inStringFrom != NULL) && (inStringTo != NULL));

	UInt32 fromLength = strlen(inStringFrom);
	if(fromLength == 0)
		return 0;

	UInt32 toLength = strlen(inStringTo);

	UInt32 count = 0;

	char * start = getDataBuffer();
	char * end = start + getLength();
	char * target;
	while(start < end)
	{
		while((target = strstr(start, inStringFrom)) != NULL)
		{
			count++;
			start = target + fromLength;
		}
		start += strlen(start) + 1;
	}

	if(count > 0)
	{
		UInt32 oldLength = getLength();
		UInt32 newLength = oldLength + (fromLength-toLength)*count;
		assureDataAvailable(newLength);
		assureDataModifiable();

		char * buffer = start = getDataBuffer();
		end = start + getLength();

		while(start < end)
		{
			while((target = strstr(start, inStringFrom)) != NULL)
			{
				UInt32 balance = oldLength - (target - buffer + toLength);
				memmove(target + toLength, target + oldLength, balance*sizeof(char));
				memcpy(target, inStringTo, toLength*sizeof(char));
				start = target + toLength;
				start[balance] = '\0';
				oldLength += (toLength - fromLength);
			}
			start += strlen(start) + 1;
		}
		U_ASSERT(getDataBuffer()[newLength] == '\0');
		setDataLength(newLength);
	}
	return count;
}


UInt32 UString::remove( char inCharacter )
{
	assureDataModifiable();

	char * source = getDataBuffer();
	char * dest = source;
	char * end = source + getLength();

	while(source < end)
	{
		if(*source != inCharacter)
		{
			*dest = *source;
			dest++;
		}
		source++;
	}
	*dest = '\0';
	UInt32 count = source - dest;
	setDataLength(getDataLength() - count);
	return count;
}


UInt32 UString::remove( const char* inString )
{
	return replace(inString, "");
}


void UString::remove( UInt32 inFirstIndex, UInt32 inCount)
{
	UInt32 length = getLength();
	U_ASSERT(inFirstIndex + inCount <= length);

	if(inCount == 0)
		return;

	assureDataModifiable();
	UInt32 numBytesToCopy = length - (inFirstIndex + inCount) + 1;
	char * buffer = getDataBuffer();
	memcpy(buffer + inFirstIndex, buffer + inFirstIndex + inCount, numBytesToCopy * sizeof(char));

	setDataLength(length - inCount);
}


void UString::insert( UInt32 inIndex, char inCharacter )
{
	UInt32 length = getLength();
	U_ASSERT(inIndex <= length);
	
	length++;
	assureDataAvailable(length);
	assureDataModifiable();

	char * buffer = getDataBuffer();
	memcpy(buffer + inIndex + 1, buffer + inIndex, (length-inIndex)*sizeof(char));
	buffer[inIndex] = inCharacter;
	setDataLength(length);
}

void UString::insert( UInt32 inIndex, const char* inString )
{
	UInt32 length = getLength();
	U_ASSERT(inIndex <= length);

	UInt32 insertLength = strlen(inString);
	if(insertLength == 0)
		return;

	length += insertLength;

	assureDataAvailable(length);
	assureDataModifiable();

	char * buffer = getDataBuffer();
	memcpy(buffer + inIndex + insertLength, buffer + inIndex, (length - inIndex - insertLength + 1)*sizeof(char));
	memcpy(buffer + inIndex, inString, insertLength*sizeof(char));

	setDataLength(length);
}

UString UString::makeFormat( const char* inFormatString, ... )
{
	va_list argList;
	va_start(argList, inFormatString);
	UString result;
	result.formatV(inFormatString, argList);
	va_end(argList);
	return result;
}

UString UString::makeFormatV( const char* inFormatString, va_list argList)
{
	UString result;
	result.formatV(inFormatString, argList);
	return result;
}

void UString::format( const char* inFormatString, ... )
{
	va_list argList;
	va_start(argList, inFormatString);
	formatV(inFormatString, argList);
	va_end(argList);
}

void UString::formatV( const char* inFormatString, va_list argList)
{
	const UInt32 FORCE_INT64 = 0x40000;

	va_list argListSave;
#ifdef _WIN32
	argListSave = argList;
#else
	va_copy(argListSave, argList);
#endif

	UInt32 maxLength = 0;
	for(const char * buffer = inFormatString; *buffer != '\0'; buffer++)
	{
		if ((*buffer != '%') || (*(++buffer) == '%'))
		{
			maxLength++;
			continue;
		}

		UInt32 itemLength = 0;

		UInt32 width = 0;
		while(*buffer != '\0')
		{
			if(*buffer == '#')
				maxLength += 2;
			else if(*buffer == '*')
				width = va_arg(argList, int);
			else if(*buffer == '-' || *buffer == '+' || *buffer == '0' ||
				*buffer == ' ')
				;
			else
				break;
			buffer++;
		}

		if(width == 0)
		{
			width = atoi(buffer);
			while((*buffer != '\0') && isdigit(*buffer))
				buffer++;
		}

		UInt32 precision = 0;
		if(*buffer == '.')
		{
			buffer++;

			if (*buffer == '*')
			{
				precision = va_arg(argList, int);
				buffer++;
			}
			else
			{
				precision = atoi(buffer);
				while((*buffer != '\0') && isdigit(*buffer))
					buffer++;
			}
		}

		UInt32 modifier = 0;
		if(strncmp(buffer, "I64", 3) == 0)
		{
			buffer += 3;
			modifier = FORCE_INT64;
		}
		else
		{
			switch (*buffer)
			{
			case 'h':
			case 'l':
			case 'F':
			case 'N':
			case 'L':
				buffer++;
				break;
			}
		}

		switch (*buffer | modifier)
		{
		case 'c':
		case 'C':
			itemLength = 2;
			va_arg(argList, int);
			break;
		case 's':
		case 'S':
			{
				char * nextArgString = va_arg(argList, char *);
				if (nextArgString == NULL)
				   itemLength = 6; // "(null)"
				else
				{
				   itemLength = strlen(nextArgString);
				   itemLength = std::max((UInt32)1, itemLength);
				}
			}
			break;
		}

		// adjust itemLength for strings
		if (itemLength != 0)
		{
			if (precision != 0)
				itemLength = std::min(itemLength, precision);
			itemLength = std::max(itemLength, width);
		}
		else
		{
			switch (*buffer)
			{
			// integers
			case 'd':
			case 'i':
			case 'u':
			case 'x':
			case 'X':
			case 'o':
				if (modifier & FORCE_INT64)
					va_arg(argList, SInt64);
				else
					va_arg(argList, int);
				itemLength = 32;
				itemLength = std::max(itemLength, width+precision);
				break;

			case 'e':
			case 'g':
			case 'G':
				va_arg(argList, double);
				itemLength = 128;
				itemLength = std::max(itemLength, width+precision);
				break;

			case 'f':
				va_arg(argList, double);
				itemLength = 128; // width isn't truncated
				// 312 == strlen("-1+(309 zeroes).")
				// 309 zeroes == max precision of a double
				itemLength = std::max(itemLength, 312+precision);
				break;

			case 'p':
				va_arg(argList, void*);
				itemLength = 32;
				itemLength = std::max(itemLength, width+precision);
				break;

			// no output
			case 'n':
				va_arg(argList, int*);
				break;

			default:
				U_ASSERT(false);  // unknown formatting option
			}
		}

		// adjust maxLength for output itemLength
		maxLength += itemLength;
	}

	char * dataBuffer = getBuffer(maxLength);
	SInt32 result = vsnprintf(dataBuffer, maxLength+1, inFormatString, argListSave);
	U_ASSERT(result <= (SInt32)getDataAllocatedLength());
	releaseBuffer();

	va_end(argListSave);
}

void UString::trimLeft()
{
	assureDataModifiable();

	char * buffer = getDataBuffer();

//printf("isspace(13): %i\n", isspace(13));
	while(isspace(*buffer))
		buffer++;

	if(buffer != getDataBuffer())
	{
		UInt32 length = getLength() - (buffer - getDataBuffer());
		memmove(getDataBuffer(), buffer, length*sizeof(char));
		setDataLength(length);
	}
}

void UString::trimLeft( char inCharacter )
{
	assureDataModifiable();

	char * buffer = getDataBuffer();

	while(inCharacter == *buffer)
		buffer++;

	if(buffer != getDataBuffer())
	{
		UInt32 length = getLength() - (buffer - getDataBuffer());
		memmove(getDataBuffer(), buffer, length*sizeof(char));
		setDataLength(length);
	}
}

void UString::trimLeft( const char* inCharacterSet )
{
	if(strlen(inCharacterSet) == 0)
		return;

	assureDataModifiable();

	char * buffer = getDataBuffer();

	while(buffer != '\0')
	{
		if(strchr(inCharacterSet, *buffer) == NULL)
			break;
		buffer++;
	}

	if(buffer != getDataBuffer())
	{
		UInt32 length = getLength() - (buffer - getDataBuffer());
		memmove(getDataBuffer(), buffer, length*sizeof(char));
		setDataLength(length);
	}
}

void UString::trimRight()
{
	assureDataModifiable();

	char * buffer = getDataBuffer();
	char * last = NULL;

	while(*buffer != '\0')
	{
		if(isspace(*buffer))
		{
			if(last == NULL)
				last = buffer;
		}
		else
		{
			last = NULL;
		}
		buffer++;
	}

	if(last != NULL)
		setDataLength(last-getDataBuffer());
}

void UString::trimRight( char inCharacter )
{
	assureDataModifiable();

	char * buffer = getDataBuffer();
	char * last = NULL;

	while(*buffer != '\0')
	{
		if(*buffer == inCharacter)
		{
			if(last == NULL)
				last = buffer;
		}
		else
		{
			last = NULL;
		}
		buffer++;
	}

	if(last != NULL)
		setDataLength(last-getDataBuffer());
}

void UString::trimRight( const char* inCharacterSet )
{
	assureDataModifiable();

	char * buffer = getDataBuffer();
	char * last = NULL;

	while(*buffer != '\0')
	{
		if(strchr(inCharacterSet, *buffer) != NULL)
		{
			if(last == NULL)
				last = buffer;
		}
		else
		{
			last = NULL;
		}
		buffer++;
	}

	if(last != NULL)
		setDataLength(last-getDataBuffer());
}

SInt32 UString::find( char inCharacter, SInt32 inFirstIndex) const
{
	if( inFirstIndex >= (SInt32)getDataLength() )
		return -1;

	const char * buffer = getDataBuffer();
	const char * firstCharacter = strchr(buffer + inFirstIndex, inCharacter);
	return (firstCharacter == NULL) ? -1 : (SInt32)(firstCharacter - buffer);
}

SInt32 UString::find( const char* inString, SInt32 inFirstIndex) const
{
	if( inFirstIndex >= (SInt32)getDataLength() )
		return -1;

	const char * buffer = getDataBuffer();
	const char * firstCharacter = strstr(buffer + inFirstIndex, inString);
	return (firstCharacter == NULL) ? -1 : (SInt32)(firstCharacter - buffer);
}

SInt32 UString::reverseFind( char inCharacter ) const
{
	const char * buffer = getDataBuffer();
	const char * lastCharacter = strrchr(buffer, inCharacter);
	return (lastCharacter == NULL) ? -1 : (SInt32)(lastCharacter - buffer);
}

SInt32 UString::findOneOf( const char* inCharacterSet, SInt32 inFirstIndex ) const
{
	if( inFirstIndex >= (SInt32)getDataLength() )
		return -1;

	const char * buffer = getDataBuffer();
	const char * character = strpbrk(buffer + inFirstIndex, inCharacterSet);
	return (character == NULL) ? -1 : (SInt32)(character - buffer);
}

void UString::split2( UString& dst1, UString& dst2, char separator ) const
{
	UString tmp(*this);
	tmp.trimLeft();
	tmp.trimRight();
	int splitIndex = -1;

	if( separator == ' ' || separator == '\t' )
	{
		int tabIndex = tmp.find('\t');
		int spaceIndex = tmp.find(' ');
		if( tabIndex == -1 )
			splitIndex = spaceIndex;
		else if( spaceIndex == -1 ) 
			splitIndex = tabIndex;
		else
			splitIndex = spaceIndex < tabIndex ? spaceIndex : tabIndex;
	} else
		splitIndex = tmp.find( separator );

	if( splitIndex == -1 )
	{
		dst1 = tmp;
		dst2 = "";
		return;
	}

	dst1 = tmp.left( splitIndex );
	dst1.trimRight();
	dst2 = tmp.mid( splitIndex + 1 );
	dst2.trimLeft();
}

void UString::split2( UString& dst1, UString& dst2 ) const
{
	split2( dst1, dst2, ' ' );
}

/*
void UString::split2( UString& dst1, UString& dst2 ) const
{
	UString tmp(*this);
	tmp.trimLeft();
	tmp.trimRight();
	int tabIndex = tmp.find('\t');
	int spaceIndex = tmp.find(' ');
	int splitIndex;
	if( tabIndex == -1 && spaceIndex == -1 )
	{
		dst1 = tmp;
		dst2 = "";
		return;
	}
	if( tabIndex == -1 )
		splitIndex = spaceIndex;
	else if( spaceIndex == -1 )
		splitIndex = tabIndex;
	else splitIndex = spaceIndex < tabIndex ? spaceIndex : tabIndex;

	dst1 = tmp.left( splitIndex );
	dst1.trimRight();
	dst2 = tmp.mid( splitIndex + 1 );
	dst2.trimLeft();
}
*/
char* UString::getBuffer( UInt32 inMinimumLength )
{
	assureDataAvailable(inMinimumLength);
	assureDataModifiable();
	return getDataBuffer();
}

char* UString::getBufferOfLength( UInt32 inLength )
{
	char * buffer = getBuffer(inLength);
	setDataLength(inLength);
	return buffer;
}

void UString::releaseBuffer( SInt32 inNewLength)
{
	assureDataModifiable();
	if(inNewLength == -1)
		inNewLength = strlen(getDataBuffer());

	U_ASSERT(inNewLength <= (SInt32)getDataAllocatedLength());
	setDataLength(inNewLength);
}

void UString::freeExtra()
{
	if( (getDataLength() + 1) < getDataAllocatedLength() )
		reallocateData( getDataLength() );
}

void UString::setData( UString::Data* inData )
{
	U_ASSERT( inData != NULL );

	acquireData( inData );

	releaseData( getData() );

	data = (char*)(inData + 1);
}

void UString::setDataLength(UInt32 inDataLength)
{
	U_ASSERT( inDataLength+1 <= getDataAllocatedLength() );
	getData()->stringLength = inDataLength;
	getDataBuffer()[inDataLength] = '\0';
}

void UString::assureDataModifiable()
{
	if( getData()->referenceCount > 1 )
		reallocateData( getDataAllocatedLength() - 1 );
	U_ASSERT(getData() != kEmptyStringData);
}

void UString::assureDataAvailable( UInt32 inMinimumLength )
{
	if( (inMinimumLength + 1) > getDataAllocatedLength() )
		reallocateData( inMinimumLength );
}

void UString::reallocateData( UInt32 inStringLength )
{
	Data* newData = (Data*)( new char[ sizeof(Data) + inStringLength + 1 ] );
	newData->allocatedLength = inStringLength + 1;
	newData->referenceCount = 0;
	newData->stringLength = std::min( inStringLength, getDataLength() );

	memcpy( (char*)(newData+1), getDataBuffer(), newData->stringLength );
	
	setData( newData );
	setDataLength( newData->stringLength );
}

void UString::acquireData( UString::Data* inData )
{
	inData->referenceCount++;
}

void UString::releaseData( UString::Data* inData )
{
	inData->referenceCount--;
	if(inData->referenceCount <= 0 )
		delete[] (char *)inData;
}
