#pragma once

#ifdef _WIN32

#include <direct.h>
#include <errno.h>


inline bool makeDirectory(const char * directory)
{
	// make directory if it doesn't exist
	int status = _mkdir(directory);
	if((status != 0) && (errno != EEXIST))
	{
		return false;
	}
	return true;
}

#else

#include <sys/stat.h>
#include <errno.h>


inline bool makeDirectory(const char * directory)
{
	// make directory if it doesn't exist
	int status = mkdir(directory, S_IRWXU);
	if((status != 0) && (errno != EEXIST))
	{
		return false;
	}
	return true;
}

#endif

