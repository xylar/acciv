// UStreamFile.cpp

#include "UStreamFile.h"
#include <stdlib.h>
#include <errno.h>

UStreamFile::UStreamFile()
	:stream(NULL)
{
}

UStreamFile::~UStreamFile()
{
	close();
}

UStreamFile::UStreamFile(FILE * inOpenStream)
{
	stream = inOpenStream;
}

UStreamFile::UStreamFile(const char * inFileName, UInt32 inOpenFlags)
	:stream(NULL)
{
	open(inFileName, inOpenFlags);
}

UStreamFile::operator FILE *() const
{
	return stream;
}

FILE * UStreamFile::getFileStream() const
{
	return stream;
}

FILE * UStreamFile::releaseFileStream()
{
	FILE * result = stream;
	stream = NULL;
	return result;
}

void UStreamFile::writeString(const char * inString) const
{
	U_ASSERT(inString != NULL);
	U_ASSERT(stream != NULL);

	if (fputs(inString, stream) == EOF)
		throw UFileException(UFileException::diskFull, fileName);
}

const char * UStreamFile::readLine(char * inString, UInt32 maxSize) const
{
	U_ASSERT(inString != NULL);
	U_ASSERT(stream != NULL);

	char * result = fgets(inString, maxSize, stream);
	if(ferror(stream))
	{
		clearerr(stream);
		throw UFileException::createFromError(errno, fileName);
	}
	return result;
}

bool UStreamFile::readLine(UString & inString) const
{
	inString = "";
	const UInt32 maxSize = 1024;
	char * buffer = inString.getBuffer(maxSize);
	char * result;
	UInt32 length = 0;
	while(true)
	{
		result = fgets(buffer, maxSize, stream);
		inString.releaseBuffer();

		if((result == NULL) && !endOfFile())
		{
			clearerr(stream);
			throw UFileException::createFromError(errno, fileName);
		}

		length = inString.getLength();
		if((result == NULL) || (length < maxSize)
			|| (inString.getCharacterAt(length-1) == '\n'))
			break;

		buffer = inString.getBuffer(maxSize + length) + length;
	}

	length = inString.getLength();
	if ((length != 0) && (inString[length-1] == '\n'))
		inString.setCharacterAt(length-1, '\0');

	return (result != NULL);
}

void UStreamFile::open(const char * inFileName, UInt32 inOpenFlags)
{
	U_ASSERT(inFileName != NULL);

	close();

	char mode[4];
	int modeIndex = 0;

	if (inOpenFlags & modeCreate)
	{
		if (inOpenFlags & modeNoTruncate)
			mode[modeIndex++] = 'a';
		else
			mode[modeIndex++] = 'w';
	}
	else if (inOpenFlags & modeWrite)
		mode[modeIndex++] = 'a';
	else
		mode[modeIndex++] = 'r';

	if (mode[0] == 'r' && (inOpenFlags & modeReadWrite) ||
		mode[0] != 'r' && !(inOpenFlags & modeWrite))
	{
		mode[modeIndex++] = '+';
	}

	if (inOpenFlags & typeBinary)
		mode[modeIndex++] = 'b';
	else
		mode[modeIndex++] = 't';
	mode[modeIndex++] = '\0';
//printf("file name: %s\n", inFileName);
//printf("mode: %s\n", mode);

	stream = fopen(inFileName, mode);
//printf("stream %i\n", stream);
	if(stream == NULL)
	{
		throw UFileException::createFromError(errno, inFileName);
	}
	else
	{
		fileName = inFileName; 
	}
}

UInt32 UStreamFile::read(void * inBuffer, UInt32 inCount) const
{
	U_ASSERT(stream != NULL);

	if(inCount == 0)
		return 0;

	UInt32 readCount = fread(inBuffer, sizeof(UInt8), inCount, stream);
	if((readCount == 0) && !endOfFile())
		throw UFileException::createFromError(errno, fileName);
	if(ferror(stream))
	{
		clearerr(stream);
		throw UFileException::createFromError(errno, fileName);
	}
	return readCount;
}

void UStreamFile::write(const void * inBuffer, UInt32 inCount) const
{
	U_ASSERT(stream != NULL);

	UInt32 writeCount = fwrite(inBuffer, sizeof(UInt8), inCount, stream);
	if(writeCount != inCount)
		throw UFileException::createFromError(errno, fileName);
}

SInt32 UStreamFile::seek(SInt32 offset, UInt32 seekPosition) const
{
	U_ASSERT(seekPosition == begin || seekPosition == end || seekPosition == current);
	U_ASSERT(stream != NULL);

	UInt32 result = fseek(stream, offset, seekPosition);
	if(result != 0)
		throw UFileException(UFileException::badSeek, fileName);

	return ftell(stream);
}

SInt32 UStreamFile::getPosition() const
{
	U_ASSERT(stream != NULL );
	return ftell(stream);
}

void UStreamFile::abort()
{
	if(stream == NULL)
		return;
	fclose(stream);
	stream = NULL;
}

void UStreamFile::flush() const
{
	if(stream == NULL)
		return;
	SInt32 result = fflush(stream);
	if(result != 0)
		throw UFileException(UFileException::diskFull, fileName);
}

void UStreamFile::close()
{
	if(stream == NULL)
		return;
	SInt32 error = fclose(stream);
	stream = NULL;
	if(error != 0)
		throw UFileException(UFileException::diskFull, fileName);
}

UInt32 UStreamFile::getLength() const
{
	UInt32 currentPosition = seek(0, current);
	UInt32 length = seek(0, end);
	UInt32 position = seek(currentPosition, begin);
	if(position != currentPosition)
		throw UFileException(UFileException::badSeek, fileName);
	return length;
}

void UStreamFile::rename(const char * oldFileName, const char * newFileName)
{
	SInt32 result = ::rename(oldFileName, newFileName);
	if(result != 0)
		throw UFileException::createFromError(errno, oldFileName);
}

void UStreamFile::remove(const char * fileName)
{
	SInt32 result = ::remove(fileName);
	if(result != 0)
		throw UFileException::createFromError(errno, fileName);
}

bool UStreamFile::endOfFile() const
{
	if(stream == NULL)
		return true;
	return (feof(stream) != 0);
}


static const char * gFileExceptionCauses[UFileException::numCauses] =
{
	"none",
	"generic",
	"fileNotFound",
	"badPath",
	"tooManyOpenFiles",
	"accessDenied",
	"invalidFile",
	"removeCurrentDir",
	"directoryFull",
	"badSeek",
	"hardIO",
	"sharingViolation",
	"lockViolation",
	"diskFull",
	"endOfFile"
};

UFileException::UFileException(UInt32 cause, const char * inFileName)
	:cause(cause)
{
	if(inFileName != NULL)
		fileName = inFileName;
}

UFileException::UFileException(const UFileException & inException)
	:cause(inException.cause), fileName(inException.fileName)
{
}

UFileException::~UFileException()
{
}

UFileException UFileException::createFromError(SInt32 errorNumber, const char * inFileName)
{
	UFileException result;
	switch ((UInt32)errorNumber)
	{
	case 0:
		result.cause = UFileException::none;
		break;
	case EPERM:
		result.cause = UFileException::accessDenied;
		break;
	case ENOENT:
		result.cause = UFileException::fileNotFound;
		break;
	case EIO:
		result.cause = UFileException::generic;
		break;
	case ENXIO:
		result.cause = UFileException::badPath;
		break;
	case EBADF:
		result.cause = UFileException::sharingViolation;
		break;
	case EACCES:
		result.cause = UFileException::accessDenied;
		break;
	case EISDIR:
		result.cause = UFileException::badPath;
		break;
	case ENOTDIR:
		result.cause = UFileException::badPath;
		break;
	case EEXIST:
		result.cause = UFileException::lockViolation;
		break;
	case ENFILE:
		result.cause = UFileException::tooManyOpenFiles;
		break;
	case EMFILE:
		result.cause = UFileException::tooManyOpenFiles;
		break;
	case EFBIG:
		result.cause = UFileException::invalidFile;
		break;
	case ENOSPC:
		result.cause = UFileException::diskFull;
		break;
	case ENAMETOOLONG:
		result.cause = UFileException::badPath;
		break;
	case ENOLCK:
		result.cause = UFileException::lockViolation;
		break;
	default:
		result.cause = UFileException::generic;
		break;
	}
	if(inFileName != NULL)
		result.fileName = inFileName;
	return result;
}

UFileException UFileException::operator = (const UFileException & inException)
{
	cause = inException.cause;
	fileName = inException.fileName;
	return *this;
}

UInt32 UFileException::getCause()
{
	return cause;
}

UString UFileException::getFileName()
{
	return fileName;
}

UString UFileException::getErrorMessage()
{
	if(cause < numCauses)
		return gFileExceptionCauses[cause];
	else
		return "";
}

