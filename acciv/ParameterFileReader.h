#pragma once

#include <core/UTypes.h>
#include <containers/UMap.h>
#include <core/UString.h>
#include <containers/UArray.h>

class ParameterFileReader
{
public:
	ParameterFileReader(const UString & fileName);
	~ParameterFileReader();

	void addFile(const UString & fileName);

	bool getInteger(const UString & tag, SInt32 & outValue) const;
	bool getDouble(const UString & tag, double & outValue) const;
	bool getString(const UString & tag, UString & outValue) const;
	bool getBool(const UString & tag, bool & outValue) const;

	bool getIntegerList(const UString & tag, UArray<SInt32> & outValues) const;
	bool getDoubleList(const UString & tag, UArray<double> & outValues) const;
	bool getBoolList(const UString & tag, UArray<bool> & outValues) const;

private:
	bool getListStrings(const UString & valueString, UArray<UString> & outListStrings) const;

	UMap<UString, UString> tagValueMap;
};
