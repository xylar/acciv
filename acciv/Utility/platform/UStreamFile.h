// UStreamFile.h

#pragma once

#include <stdio.h>
#include "../core/Utility.h"
#include "../core/UString.h"

class UFileException;

class UStreamFile  
{
public:
	enum OpenFlags {
		modeRead =          0x0000,
		modeWrite =         0x0001,
		modeReadWrite =     0x0002,
		shareCompat =       0x0000,
		shareExclusive =    0x0010,
		shareDenyWrite =    0x0020,
		shareDenyRead =     0x0030,
		shareDenyNone =     0x0040,
		modeNoInherit =     0x0080,
		modeCreate =        0x1000,
		modeNoTruncate =    0x2000,
		typeText =          0x4000,
		typeBinary =   (int)0x8000
		};

	enum SeekPosition { begin = 0x0, current = 0x1, end = 0x2 };

	UStreamFile();
	virtual ~UStreamFile();
	UStreamFile(FILE * inOpenStream);
	UStreamFile(const char * inFileName, UInt32 inOpenFlags);

	operator FILE *() const;
	FILE * getFileStream() const;
	FILE * releaseFileStream();

	// not necessarily the number of bytes in the file
	// because of the carriage return / line feed weirdness
	// but can be used to allocate a buffer to take the
	// characters from the file
	UInt32 getLength() const;

	void writeString(const char * inString) const;
	const char * readLine(char * inString, UInt32 maxSize) const;
	bool readLine(UString & inString) const;

	void open(const char * inFileName, UInt32 inOpenFlags);
	UInt32 read(void * inBuffer, UInt32 inCount) const;
	void write(const void * inBuffer, UInt32 inCount) const;
	SInt32 seek(SInt32 offset, UInt32 seekPosition) const;
	SInt32 getPosition() const;

	void abort();
	void flush() const;
	void close();
	bool endOfFile() const;

	static void rename(const char * oldFileName, const char * newFileName);
	static void remove(const char * fileName);

private:

	FILE * stream;
	UString fileName;

};

class UFileException
{
public:
	enum {
		none = 0,
		generic,
		fileNotFound,
		badPath,
		tooManyOpenFiles,
		accessDenied,
		invalidFile,
		removeCurrentDir,
		directoryFull,
		badSeek,
		hardIO,
		sharingViolation,
		lockViolation,
		diskFull,
		endOfFile,
		numCauses
	};

	UFileException(UInt32 cause = UFileException::none, const char * fileName = NULL);
	UFileException(const UFileException & inException);
	~UFileException();
	static UFileException createFromError(SInt32 errorNumber, const char * fileName = NULL);
	UFileException operator = (const UFileException & inException);

	UInt32 getCause();
	UString getFileName();
	UString getErrorMessage();
private:
	UInt32 cause;
	UString fileName;
};
