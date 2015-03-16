// UList.h
#pragma once

#include "../core/Utility.h"
#include "UAllocationPool.h"

template< class T >
class UList
{
public:
	typedef void* Position;

	UList( UInt32 inGrowBy = 10 );
	~UList();

// Attributes
	UInt32 getSize() const;

// Operations
	// Clean up
	void removeAll();

	// element access
	T getFirst() const;
	T removeFirst();
	Position addToFront( const T& inElement );
	void appendToFront( const UList<T>& inList );

	T getLast() const;
	T removeLast();
	Position addToBack( const T& inElement );
	void appendToBack( const UList<T>& inList );

	Position getFirstPosition() const;
	Position getLastPosition() const;

	T getNextElement( Position& ioPosition ) const;
	T getPreviousElement( Position& ioPosition ) const;
	
	T getAtPosition( Position inPosition ) const;
	void setAtPosition( Position inPosition, const T& inElement );
	T removeAtPosition( Position inPosition );

	Position insertBeforePosition( Position inPosition, const T& inElement );
	Position insertAfterPosition( Position inPosition, const T& inElement );

	Position find( const T& inElement, Position inStartAfter = NULL ) const;
	Position findIndex( UInt32 inIndex ) const;

private:
	class Node
	{
	public:
		Node* next;
		Node* prev;
		T element;
	};

	Node* allocateNode( const T& inElement );
	void releaseNode( Node* inNode );

	UInt32 size;
	Node* first;
	Node* last;
	Node* freeList;
	UAllocationPool<Node> allocationPool;
	UInt32 growBy;
};

template< class T >
UList<T>::UList<T>( UInt32 inGrowBy )
{
	size = 0;
	growBy = inGrowBy;
	first = last = freeList = NULL;
}

template< class T >
UList<T>::~UList<T>()
{
	removeAll();
}

template< class T >
UInt32 UList<T>::getSize() const
{
	return size;
}

template< class T >
void UList<T>::removeAll()
{
	while(first)
	{
		Node* temp = first;
		first = first->next;
		releaseNode(temp);
	}
	allocationPool.releaseAll();
	size = 0;
	first = last = freeList = NULL;
}

template< class T >
T UList<T>::getFirst() const
{
	U_ASSERT( first != NULL );
	return first->element;
}

template< class T >
T UList<T>::removeFirst()
{
	U_ASSERT( first != NULL );

	Node* node = first;
	first = first->next;
	if(first)
		first->prev = NULL;
	else
		last = NULL;

	T result = node->element;
	releaseNode( node );

	size--;

	return result;
}

template< class T >
typename UList<T>::Position UList<T>::addToFront( const T& inElement )
{
	Node* newNode = allocateNode( inElement );
	U_ASSERT( newNode != NULL );

	newNode->next = first;
	newNode->prev = NULL;
	if(first == NULL)
	{
		last = newNode;
	}
	else
	{
		first->prev = newNode;
	}

	first = newNode;

	size++;

	return (Position)newNode;
}

template< class T >
void UList<T>::appendToFront( const UList<T>& inList )
{
	Node* node = inList.last;

	while( node )
	{
		addToFront( node->element );
		node = node->prev;
	}
}

template< class T >
T UList<T>::getLast() const
{
	U_ASSERT( last != NULL );
	return last->element;
}

template< class T >
T UList<T>::removeLast()
{
	U_ASSERT( last != NULL );

	Node* node = last;
	last = last->prev;
	if(last)
		last->next = NULL;
	else
		first = NULL;

	T result = node->element;

	releaseNode( node );

	size--;

	return result;
}

template< class T >
typename UList<T>::Position UList<T>::addToBack( const T& inElement )
{
	Node* newNode = allocateNode( inElement );
	U_ASSERT( newNode != NULL );

	newNode->next = NULL;
	newNode->prev = last;
	if(last == NULL)
	{
		first = newNode;
	}
	else
	{
		last->next = newNode;
	}
	last = newNode;
	size++;
	return (Position)newNode;
}

template< class T >
void UList<T>::appendToBack( const UList<T>& inList )
{
	Node* node = inList.first;

	while( node )
	{
		addToBack( node->element );
		node = node->next;
	}
}

template< class T >
typename UList<T>::Position UList<T>::getFirstPosition() const
{
	return (Position)( first );
}

template< class T >
typename UList<T>::Position UList<T>::getLastPosition() const
{
	return (Position)( last );
}

template< class T >
T UList<T>::getNextElement( Position& ioPosition ) const
{
	Node* node = (Node*)(ioPosition);
	U_ASSERT( node != NULL );

	ioPosition = (Position)( node->next );
	return node->element;
}

template< class T >
T UList<T>::getPreviousElement( Position& ioPosition ) const
{
	Node* node = (Node*)(ioPosition);
	U_ASSERT( node != NULL );

	ioPosition = (Position)( node->prev );
	return node->element;
}

template< class T >
T UList<T>::getAtPosition( Position inPosition ) const
{
	Node* node = (Node*)(inPosition);
	U_ASSERT( node != NULL );

	return node->element;
}

template< class T >
void UList<T>::setAtPosition( Position inPosition, const T& inElement )
{
	Node* node = (Node*)(inPosition);
	U_ASSERT( node != NULL );

	node->element = inElement;
}

template< class T >
T UList<T>::removeAtPosition( Position inPosition )
{
	Node* node = (Node*)(inPosition);
	U_ASSERT( node != NULL );

	if( node->next )
		node->next->prev = node->prev;
	else
		last = node->prev;

	if( node->prev )
		node->prev->next = node->next;
	else
		first = node->next;

	T result = node->element;
	releaseNode( node );

	size--;

	return result;
}

template< class T >
typename UList<T>::Position UList<T>::insertBeforePosition( Position inPosition, const T& inElement )
{
	Node* node = (Node*)(inPosition);
	//U_ASSERT( node != NULL );
	if(node == NULL)
		return addToFront(inElement);

	Node* newNode = allocateNode( inElement );
	U_ASSERT( newNode != NULL );

	newNode->next = node;
	newNode->prev = node->prev;

	if( node->prev )
		node->prev->next = newNode;
	else
		first = newNode;

	node->prev = newNode;

	size++;

	return (Position)newNode;
}

template< class T >
typename UList<T>::Position UList<T>::insertAfterPosition( Position inPosition, const T& inElement )
{
	Node* node = (Node*)(inPosition);
//	U_ASSERT( node != NULL );
	if(node == NULL)
		return addToBack(inElement);

	Node* newNode = allocateNode( inElement );
	U_ASSERT( newNode != NULL );

	newNode->next = node->next;
	newNode->prev = node;

	if( node->next )
		node->next->prev = newNode;
	else
		last = newNode;

	node->next = newNode;

	size++;

	return (Position)newNode;
}

template< class T >
typename UList<T>::Position UList<T>::find( const T& inElement, Position inStartAfter ) const
{
	Node* node = (Node*)(inStartAfter);
	if( node == NULL )
		node = first;

	while( node )
	{
		if( node->element == inElement )
			return (Position)( node );
		node = node->next;
	}

	return (Position)( NULL );
}

template< class T >
typename UList<T>::Position UList<T>::findIndex( UInt32 inIndex ) const
{
	Node* node = first;

	while( inIndex )
	{
		U_ASSERT( node );
		node = node->next;
		inIndex--;
	}
	
	return (Position)( node );
}

template< class T >
typename UList<T>::Node* UList<T>::allocateNode( const T& inElement )
{
	if( !freeList )
	{
		Node* newNodes = allocationPool.allocate( growBy );
		U_ASSERT( newNodes != NULL );

		for( UInt32 i = 0; i < growBy; i++ )
		{
			newNodes->next = freeList;
			freeList = newNodes;
			newNodes++;
		}
	}

	Node* newNode = freeList;
	freeList = freeList->next;

	::new((void*)&(newNode->element)) T( inElement );

	return newNode;
}

template< class T >
void UList<T>::releaseNode( typename UList<T>::Node* inNode )
{
	U_ASSERT( inNode != NULL );
	
	(inNode->element).~T();
	
	inNode->next = freeList;
	freeList = inNode;
}