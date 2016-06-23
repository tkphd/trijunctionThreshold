// File: sparse_priority_queue.h
// Declaration and Implementation of SparseDistanceVoxel and SparsePriorityQueue classes

#include <iostream>
#include <climits>
#include <vector>
#include <map>
#include <cassert>
#include <cmath>
#include <omp.h>
#include "MMSP.sparse.hpp"

#ifndef _DISTANCE_VOXEL_H_
#define _DISTANCE_VOXEL_H_

class SparseDistanceVoxel {
public:
  SparseDistanceVoxel() {
  	x[0]=x[1]=x[2]=0;
  }
  // accessor
  int getPosition(const int d) const {return x[d];}
  double getValue(const int id) const {
  	for (int i=0; i<distance.length(); ++i)	{
  		int index=MMSP::index(distance,i);
  		if (index==id) return distance.value(i);
  	}
  	return std::sqrt(std::numeric_limits<double>::max());
  }
  int length() const {return distance.length();}
  int index(int j) const {return distance.index(j);}
  // modifier
  void setPosition(const int d, int i){x[d]=i;}
  void setValue(const int index, double v) {distance.set(index) = v;}
  omp_lock_t* ompLock() {return &omp_lock;}
private:
  // REPRESENTATION
  int x[3];      									// position in the grid
  MMSP::sparse<double> distance;	// how far away each vertex is
  omp_lock_t omp_lock;
};
#endif // _DISTANCE_VOXEL_H_

#ifndef _SPARSE_PRIORITY_QUEUE_H_
#define _SPARSE_PRIORITY_QUEUE_H_

// The SparsePriorityQueue is a customized, non-templated
// priority queue that stores SparseDistanceVoxel pointers in a heap.  The
// elements in the heap can be looked up in a map, to quickly find out
// the current index of the element within the heap.

// =========================================================================

class SparsePriorityQueue {

public:
	// --------------------------
	// CONSTRUCTORS
	SparsePriorityQueue() {} // Default Constructor
	SparsePriorityQueue(const int index, const std::vector<SparseDistanceVoxel*> &values ) {
		// Class Constructor
		// Build heap from input vector using PUSH member function
		for ( std::vector<SparseDistanceVoxel*>::const_iterator itr = values.begin(); itr != values.end(); ++itr)
			push(index, *itr );
	}
	// NOTE: No Copy Constructor or Assignment Operator are defined.

	// ------------------------
	// ACCESSORS
	int size() {return m_heap.size();}
	bool empty() {return m_heap.empty();}
	int last_non_leaf() {return ( size() - 1 ) / 2;}
	int get_parent( int i ) {
		assert ( i > 0 );
		assert( i <= size() );
		return ( i - 1 ) / 2;
	}
	int has_left_child( int i ) {return ( 2 * i ) + 1 < size();}
	int has_right_child( int i ) {return ( 2 * i ) + 2 < size();}
	int get_left_child( int i ) {
		assert ( i >= 0 && has_left_child( i ) );
		return 2 * i + 1;
	}
	int get_right_child( int i ) {
		assert ( i >= 0 && has_right_child( i ) );
		return 2*i + 2;
	}

	// read the top element
	const SparseDistanceVoxel* top() const {
		assert( !m_heap.empty() );
		return m_heap[0];
	}

	// is this element in the heap?
	bool in_heap( SparseDistanceVoxel* element ) const {
		std::map<SparseDistanceVoxel*, int>::const_iterator itr = backpointers.find( element );
		return ( itr != backpointers.end() );
	}

	// add an element to the heap
	void push(const int index, SparseDistanceVoxel* element ) {
		std::map<SparseDistanceVoxel*, int>::iterator itr = backpointers.find( element );
		assert ( itr == backpointers.end() );
		m_heap.push_back( element );
		backpointers[element] = m_heap.size() - 1;
		this->percolate_up( index, int( m_heap.size() - 1 ) );
	}

	// the value of this element has been edited, move the element up or down
	void update_position(const int index, SparseDistanceVoxel* element ) {
		std::map<SparseDistanceVoxel*, int>::iterator itr = backpointers.find( element );
		assert ( itr != backpointers.end() );
		this->percolate_up( index, itr->second );
		this->percolate_down( index, itr->second );
	}

	// remove the top (minimum) element
	void pop(const int index) {
		assert( !m_heap.empty() );
		int success = backpointers.erase( m_heap[0] );
		assert ( success == 1 );
		m_heap[0] = m_heap.back();
		m_heap.pop_back();
		this->percolate_down( index, 0 );
	}

private:
	// REPRESENTATION
	//  the heap is stored in a vector representation (the binary tree
	//  structure "unrolled" one row at a time)
	std::vector<SparseDistanceVoxel*> m_heap;

	// the map stores a correpondence between elements & indices in the heap
	std::map<SparseDistanceVoxel*,int> backpointers;

	// Helper function for percolators: swap contents of heap
	void swap( int i, int j) {
		// Change positions in heap
		SparseDistanceVoxel* temp = m_heap[j];
		m_heap[j] = m_heap[i];
		m_heap[i] = temp;
		// Update positions in map
		backpointers[m_heap[i]] = i;
		backpointers[m_heap[j]] = j;
	}

	// private helper functions
	void percolate_up(const int j, int i) {
		while (i > 0) {
			// Loop until you become the root node (i==0)
			if ( m_heap[i]->getValue(j) < m_heap[get_parent(i)]->getValue(j) ) {
				swap(i, get_parent(i));
				i = get_parent(i);
			} else {
				break;
			}
		}
	}

	void percolate_down(const int j, int i) {
		while ( has_left_child(i) ) {
			// Loop as long as there's at least one level below you
			int child;
			// Choose the child to compare against
			if ( has_right_child(i) && m_heap[get_right_child(i)]->getValue(j) < m_heap[get_left_child(i)]->getValue(j)) {
				child = get_right_child(i);
			} else child = get_left_child(i);
			if ( m_heap[child]->getValue(j) < m_heap[i]->getValue(j) ) {
				swap(child, i);
				i = child;
			} else {
				break;
			}
		}
	}

};

#endif
