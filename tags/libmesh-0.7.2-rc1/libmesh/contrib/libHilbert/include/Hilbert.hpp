/*
 * Copyright (C) 2006-2007 Chris Hamilton <chamilton@cs.dal.ca>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef _HILBERT_HPP_
#define _HILBERT_HPP_


#include "Hilbert/FixBitVec.hpp"
#include "Hilbert/BigBitVec.hpp"


// Description of parameters:
//
// FOR REGULAR HILBERT INDICES
//
// CFixBitVec/CBigBitVec *p
// Pointer to array of non-negative coordinate values.
//
// int m
// Precision of all coordinate values (number of bits required to
// represent the largest possible coordinate value).
//
// int n
// Number of dimensions (size of the array *p).
//
// CFixBitVec/CBigBitVec &h
// Hilbert index of maximum precision m*n.
//
// int *ms
// Array of precision values, one per dimension.
//
// FOR COMPACT HILBERT INDICES
//
// CFixBitVec/CBigBitVec &hc
// Compact Hilbert index of maximum precision M.
//
// int M
// Net precision value, corresponding to the size of the compact
// Hilbert code.  If not provided, defaults to zero and will be calculated
// by the function (sum_i { ms[i] }).
//
// int m
// Largest precision value (max_i { ms[i] }).  If not provided, defaults
// to zero and will be calculated by the function,


namespace Hilbert
{
	// fix -> fix
	void coordsToIndex( const CFixBitVec *p, int m, int n, CFixBitVec &h );
	void indexToCoords( CFixBitVec *p, int m, int n, const CFixBitVec &h );
	void coordsToCompactIndex( const CFixBitVec *p, const int *ms, int n,
		CFixBitVec &hc, int M = 0, int m = 0 );
	void compactIndexToCoords( CFixBitVec *p, const int *ms, int n,
		const CFixBitVec &hc, int M = 0, int m = 0 );

	// fix -> big
	void coordsToIndex( const CFixBitVec *p, int m, int n, CBigBitVec &h );
	void indexToCoords( CFixBitVec *p, int m, int n, const CBigBitVec &h );
	void coordsToCompactIndex( const CFixBitVec *p, const int *ms, int n,
		CBigBitVec &hc, int M = 0, int m = 0 );
	void compactIndexToCoords( CFixBitVec *p, const int *ms, int n,
		const CBigBitVec &hc, int M = 0, int m = 0 );
	
	// big -> big
	void coordsToIndex( const CBigBitVec *p, int m, int n, CBigBitVec &h );
	void indexToCoords( CBigBitVec *p, int m, int n, const CBigBitVec &h );
	void coordsToCompactIndex( const CBigBitVec *p, const int *ms, int n,
		CBigBitVec &hc, int M = 0, int m = 0 );
	void compactIndexToCoords( CBigBitVec *p, const int *ms, int n,
		const CBigBitVec &hc, int M = 0, int m = 0 );
}


#endif
