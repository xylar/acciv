// UAllocationPool.h
#pragma once

#include "../core/Utility.h"

template< class T >
class UAllocationPool
{
public:
	UAllocationPool();
	~UAllocationPool();

	T* allocate( UInt32 inCount );
	void releaseAll();

private:
	struct Block
	{
		Block* next;
	};

	Block* firstBlock;
};

template< class T >
UAllocationPool<T>::UAllocationPool()
{
	firstBlock = NULL;
}

template< class T >
UAllocationPool<T>::~UAllocationPool()
{
	releaseAll();
}

template< class T >
T* UAllocationPool<T>::allocate( UInt32 inCount )
{
	UInt32 size = sizeof(Block) + inCount*sizeof(T);
	Block* block = (Block *)(new UInt8[ size ]);
	block->next = firstBlock;
	firstBlock = block;

	T* data = (T*)( (UInt8*)block + sizeof(Block) );
	return data;
}

template< class T >
void UAllocationPool<T>::releaseAll()
{
	Block* block = firstBlock;
	while( block )
	{
		Block* next = block->next;
		delete block;
		block = next;
	}
	firstBlock = NULL;
}

