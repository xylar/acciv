#include "ParameterFileReader.h"
#include <platform/UStreamFile.h>
#include <stdlib.h>

ParameterFileReader::ParameterFileReader(const UString & fileName)
{
	addFile(fileName);
}

ParameterFileReader::~ParameterFileReader(void)
{
}

void ParameterFileReader::addFile(const UString & fileName)
{
	try
	{
		UStreamFile file(fileName, UStreamFile::modeRead|UStreamFile::shareDenyWrite|UStreamFile::typeText);
		UString line, temp;
		while(file.readLine(line))
		{
			if((line.getLength() > 0) && !(line.getCharacterAt(0) == '#'))
			{
				SInt32 breakIndex = line.find('=');
				if(breakIndex != -1)
				{
					UString tag = line.left(breakIndex);
					while((tag.getLength() > 0) && (tag.getCharacterAt(tag.getLength()-1) == ' '))
						tag = tag.left(tag.getLength()-1);
					UString value = line.mid(breakIndex+1);
					while((value.getLength() > 0) && (value.getCharacterAt(0) == ' '))
						value = value.mid(1);

					// stupid DOS file problem: have to remove ^M characters
					if((value.getLength() > 0) && (value.getCharacterAt(value.getLength()-1) == 13))
						value = value.left(tag.getLength()-1);

					while((value.getLength() > 0) && (value.getCharacterAt(value.getLength()-1) == ' '))
						value = value.left(value.getLength()-1);

					// only insert the tag/value if there is not already a value associated with the tag
					if((tag.getLength() > 0) && (value.getLength() > 0) && !tagValueMap.lookup(tag,temp))
							tagValueMap.insert(tag, value);
				}
			}
		}
	}
	catch (UFileException & e)
	{
		fprintf(stderr, "ParameterFileReader::ParameterFileReader: Error while reading parameter file.");
	}
}


bool ParameterFileReader::getInteger(const UString & tag, SInt32 & outValue) const
{
	UString valueString;
	if(!tagValueMap.lookup(tag, valueString))
		return false;
	outValue = atoi(valueString);
	return true;
}

bool ParameterFileReader::getDouble(const UString & tag, double & outValue) const
{
	UString valueString;
	if(!tagValueMap.lookup(tag, valueString))
		return false;
	outValue = atof(valueString);
	return true;
}

bool ParameterFileReader::getString(const UString & tag, UString & outValue) const
{
	return tagValueMap.lookup(tag, outValue);
}

bool ParameterFileReader::getBool(const UString & tag, bool & outValue) const
{
	UString valueString;
	if(!tagValueMap.lookup(tag, valueString))
		return false;
	valueString.makeLower();
	outValue = (valueString == "true") || (valueString == "t") || (valueString == "1");
	return true;
}

bool ParameterFileReader::getIntegerList(const UString & tag, UArray<SInt32> & outValues) const
{
	UString valueString;
	if(!tagValueMap.lookup(tag, valueString))
		return false;
	UArray<UString> listStrings;
	if(!getListStrings(valueString, listStrings))
		return false;
	for(SInt32 index = 0; index < listStrings.getSize(); index++)
		outValues.add(atoi(listStrings[index]));
	return true;
}

bool ParameterFileReader::getDoubleList(const UString & tag, UArray<double> & outValues) const
{
	UString valueString;
	if(!tagValueMap.lookup(tag, valueString))
		return false;
	UArray<UString> listStrings;
	if(!getListStrings(valueString, listStrings))
		return false;
	for(SInt32 index = 0; index < listStrings.getSize(); index++)
		outValues.add(atof(listStrings[index]));
	return true;
}

bool ParameterFileReader::getBoolList(const UString & tag, UArray<bool> & outValues) const
{
	UString valueString;
	if(!tagValueMap.lookup(tag, valueString))
		return false;
	UArray<UString> listStrings;
	if(!getListStrings(valueString, listStrings))
		return false;
	for(SInt32 index = 0; index < listStrings.getSize(); index++)
		outValues.add((listStrings[index] == "true") || (listStrings[index] == "t") || (listStrings[index] == "1"));
	return true;
}

bool ParameterFileReader::getListStrings(const UString & valueString, UArray<UString> & outListStrings) const
{
	SInt32 openIndex = valueString.find('[');
	SInt32 closeIndex = valueString.find(']');
	if((openIndex == -1) || (closeIndex == -1))
		return false;
	UString subString = valueString.mid(openIndex+1,closeIndex-openIndex-1);
	UString left, right;
	char separator = ',';
	subString.split2(left, right, separator);
	if(right.getLength() == 0)
	{
		separator = ' ';
		subString.split2(left, right, separator);
                // if right is still empty then there is just one string, which is fine
	}
	outListStrings.add(left);
	while(right.getLength() != 0)
	{
		subString = right;
		subString.split2(left, right, separator);
		outListStrings.add(left);
	}
	return true;
}
