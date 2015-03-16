//UMap.h
#pragma once

#include "../core/Utility.h"
#include "../core/UString.h"
#include "UAllocationPool.h"
#include <memory.h>


template< class Key, class Value >
class UMap
{
public:
	UMap( UInt32 inTableSize = 17, UInt32 inGrowBy = 10 );
	~UMap();

	UInt32 getSize() const;

	bool lookup( const Key& inKey, Value& outValue ) const;

	void insert( const Key& inKey, const Value& inValue );

	void remove( const Key& inKey );
	void removeAll();

	typedef void* Position;
	Position getFirstPosition() const;
	void getNextAssociation( Position& ioPosition, Key& outKey, Value& outValue ) const;

	UInt32 getTableSize() const;
	void setTableSize( UInt32 inSize );

protected:
	struct Association
	{
		Association* next;
		UInt32 hash;
		Key key;
		Value value;
	};

	Association* allocateAssociation( const Key& inKey, UInt32 inHash, const Value& inValue );
	void releaseAssociation( Association* inAssociation );

	Association* findAssociation( const Key& inKey, UInt32 inHash) const;
	Association* getAssociation( const Key& inKey, UInt32 inHash );

	UInt32 getHashKey( const Key& inKey) const;

	void allocateTable();

	Association** table;
	UInt32 tableSize;
	UInt32 size;

	Association* freeList;

	UAllocationPool<Association> allocationPool;
	UInt32 growBy;
};

#pragma warning(push)
#pragma warning(disable: 4311)
#pragma warning(disable: 4312)

template< class Key >
inline UInt32 UHashKey( const Key& key )
{
	// default identity hash - works for most primitive values
	return ((UInt32)(void*)(UInt32)key) >> 4;
}

#pragma warning(pop)

template<>
inline UInt32 UHashKey<const char*>( const char* const & key )
{
	const char* array = key;
	UInt32 hash = 0;
	while( *array )
		hash = (hash<<5) + hash + *array++;
	return hash;
}

template<>
inline UInt32 UHashKey<UString>( const UString& key )
{
	return UHashKey( (const char*) key );
}

template< class Key, class Value >
UMap<Key,Value>::UMap( UInt32 inTableSize, UInt32 inGrowBy )
{
	table = NULL;
	tableSize = inTableSize;
	size = 0;
	freeList = NULL;
	growBy = inGrowBy;
}

template< class Key, class Value >
UMap<Key,Value>::~UMap()
{
	removeAll();
	delete table;
}

template< class Key, class Value >
UInt32 UMap<Key,Value>::getSize() const
{
	return size;
}

template< class Key, class Value >
bool UMap<Key,Value>::lookup( const Key& inKey, Value& outValue ) const
{
	Association* association = findAssociation( inKey, getHashKey(inKey) );
	if( association )
	{
		outValue = association->value;
		return true;
	}
	return false;
}

template< class Key, class Value >
void UMap<Key,Value>::insert( const Key& inKey, const Value& inValue )
{
	if( !table )
		allocateTable();

	UInt32 hash = getHashKey(inKey);

	Association* association = table[hash];
	while( association )
	{
		if( association->key == inKey )
		{
			association->value = inValue;
			return;
		}
		association = association->next;
	}

	// allocate a new association and return it
	Association* newAssociation = allocateAssociation( inKey, hash, inValue );
	newAssociation->next = table[hash];
	table[hash] = newAssociation;
}

template< class Key, class Value >
void UMap<Key,Value>::remove( const Key& inKey )
{
	if( !table )
		return;

	UInt32 hash = getHashKey(inKey);
	Association* previous = NULL;
	Association* association = table[hash];
	while( association )
	{
		if( association->key == inKey )
		{
			if( previous )
				previous->next = association->next;
			else
				table[hash] = association->next;

			releaseAssociation( association );

			return;
		}
		previous = association;
		association = association->next;
	}
}

template< class Key, class Value >
void UMap<Key,Value>::removeAll()
{
	if( !table || (size == 0) )
		return;

	for( UInt32 i = 0; i < tableSize; i++ )
	{
		Association* association = table[i];
		while( association )
		{
			Association* nextAssociation = association->next;
			releaseAssociation( association );
			association = nextAssociation;
		}
	}

	allocationPool.releaseAll();
	size = 0;
	freeList = NULL;
	
	if( table )
	{
		for( UInt32 i = 0; i < tableSize; i++ )
			table[i] = NULL;
	}
}

template< class Key, class Value >
typename UMap<Key,Value>::Position UMap<Key,Value>::getFirstPosition() const
{
	if( !table || (size == 0) )
		return (Position)( NULL );

	for( UInt32 i = 0; i < tableSize; i++ )
	{
		Association* association = table[i];
		if( association )
			return (Position)( association );
	}
	return (Position)( NULL );
}

template< class Key, class Value >
void UMap<Key,Value>::getNextAssociation( typename UMap<Key,Value>::Position& ioPosition, Key& outKey, Value& outValue ) const
{
	U_ASSERT( table != NULL );

	Association* association = (Association*)( ioPosition );
	U_ASSERT( association != NULL );

	outKey = association->key;
	outValue = association->value;

	Association* nextAssociation = association->next;
	UInt32 i = association->hash + 1;
	while( (nextAssociation == NULL ) && (i < tableSize) )
	{
		if( table[i] )
			nextAssociation = table[i];
		i++;
	}
	ioPosition = (Position)( nextAssociation );
}

template< class Key, class Value >
UInt32 UMap<Key,Value>::getTableSize() const
{
	return tableSize;
}

template< class Key, class Value >
void UMap<Key,Value>::setTableSize( UInt32 inSize )
{
	Association* associationsToRehash = NULL;
	if( table )
	{
		for( UInt32 i = 0; i < tableSize; i++ )
		{
			Association* association = table[i];
			while( association )
			{
				Association* temp = association;
				association = association->next;

				temp->next = associationsToRehash;
				associationsToRehash = temp;
			}
		}
		delete table;
	}

	tableSize = inSize;
	allocateTable();

	while( associationsToRehash )
	{
		Association* next = associationsToRehash->next;

		UInt32 hash = getHashKey( associationsToRehash->key );
		associationsToRehash->hash = hash;

		associationsToRehash->next = table[hash];
		table[hash] = associationsToRehash;

		associationsToRehash = next;
	}
}

template< class Key, class Value >
typename UMap<Key,Value>::Association* UMap<Key,Value>::allocateAssociation( const Key& inKey, UInt32 inHash, const Value& inValue )
{
	if( !freeList )
	{
		Association* newAssociations = allocationPool.allocate( growBy );

		for( UInt32 i = 0; i < growBy; i++ )
		{
			newAssociations->next = freeList;
			freeList = newAssociations;
			newAssociations++;
		}
	}

	Association* association = freeList;
	freeList = freeList->next;

	::new((void*)&(association->key)) Key( inKey );
	::new((void*)&(association->value)) Value( inValue );
	association->hash = inHash;

	size++;

	return association;
}

template< class Key, class Value >
void UMap<Key,Value>::releaseAssociation( typename UMap<Key,Value>::Association* inAssociation )
{
	U_ASSERT( inAssociation != NULL );

	(inAssociation->key).~Key();
	(inAssociation->value).~Value();

	inAssociation->next = freeList;
	freeList = inAssociation;

	size--;
}

template< class Key, class Value >
typename UMap<Key,Value>::Association* UMap<Key,Value>::findAssociation( const Key& inKey, UInt32 inHash) const
{
	if( !table )
		return NULL;

	Association* association = table[inHash];
	while( association )
	{
		if( association->key == inKey )
			return association;
		association = association->next;
	}
	
	return NULL;
}

template< class Key, class Value >
typename UMap<Key,Value>::Association* UMap<Key,Value>::getAssociation( const Key& inKey, UInt32 inHash )
{
	if( !table )
		allocateTable();

	Association* association = table[inHash];
	while( association )
	{
		if( association->key == inKey )
			return association;
		association = association->next;
	}

	// allocate a new association and return it
	Association* newAssociation = allocateAssociation( inKey, inHash );
	newAssociation->next = table[inHash];
	table[inHash] = newAssociation;

	return newAssociation;
}

template< class Key, class Value >
void UMap<Key,Value>::allocateTable()
{
	table = (Association**) new UInt8[ tableSize*sizeof(Association*) ];
	memset( table, 0, tableSize*sizeof(Association*) );
}

template< class Key, class Value >
UInt32 UMap<Key,Value>::getHashKey( const Key& inKey) const
{
	return UHashKey(inKey) % tableSize;
}

