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

#ifndef _FIXBITVEC_HPP_
#define _FIXBITVEC_HPP_


#include <inttypes.h>
#include "Hilbert/Common.hpp"


// This must be an unsigned integer that is either
// 32 or 64 bits.  Otherwise, there are places in the
// code that simply will not work.
// For speed, this should be the native word size.
//typedef uint64_t FBV_UINT;
//#define FBV_BITS		64
typedef uint32_t FBV_UINT;
#define FBV_BITS		32

#define FBV0		((FBV_UINT)0)
#define FBV1		((FBV_UINT)1)
#define FBV1S		(~FBV0)
#define FBVN1S(n)	(n==FBV_BITS?FBV1S:(FBV1<<n)-1)
#define FBVMOD(i,m)	if((i)>=(m))(i)-=(m)*((i)/(m));


typedef enum 
{
	eFix,
	eBig
} EBitVecType;


class CFixBitVec
{
  public:

	static
	EBitVecType
	type();

	// Default constructor.  The bits parameter
	// is completely ignored, but accepted in order
	// to look and feel the same as a BigBitVec.
	CFixBitVec(
		int iBits = FBV_BITS
		);
	
	// Copy constructor.
	CFixBitVec(
		const CFixBitVec &cFBV
		);
	
	// Returns the current size in bits.
	int
	getSize();
	
	// Sets the size.  This is a dummy
	// function just for BigBitVec compatibility.
	CFixBitVec &
	setSize(
		int iBits
		);

	// Zeros the bit-vector.
	CFixBitVec &
	zero();

	// Truncates the bit-vector to a given precision in
	// bits (zeroes MSBs without shrinking the vector)
	CFixBitVec &
	truncate(
		int iBits
		);

	// Assignment operator.
	CFixBitVec &
	operator=(
		const CFixBitVec &cFBV
		);

	// Assignment operator.
	CFixBitVec &
	operator=(
		FBV_UINT i
		);

	// Returns the value of the nth bit.
	bool
	getBit(
		int iIndex
		) const;

	// Sets the value of the nth bit.
	CFixBitVec &
	setBit(
		int iIndex,
		bool bBit
		);

	// Toggles the value of the nth bit.
	CFixBitVec &
	toggleBit(
		int iIndex
		);

	// AND operation in place.
	CFixBitVec &
	operator&=(
		const CFixBitVec &cFBV
		);
	
	CFixBitVec &
	operator&=(
		FBV_UINT i
		);

	// AND operation.
	CFixBitVec
	operator&(
		const CFixBitVec &cFBV
		) const;
	CFixBitVec
	operator&(
		FBV_UINT i
		);
	
	// OR operation in place.
	CFixBitVec &
	operator|=(
		const CFixBitVec &cFBV
		);
	CFixBitVec &
	operator|=(
		FBV_UINT i
		);

	// OR operation.
	CFixBitVec
	operator|(
		const CFixBitVec &cFBV
		) const;
	CFixBitVec
	operator|(
		FBV_UINT i
		);

	// XOR operation in place.
	CFixBitVec &
	operator^=(
		const CFixBitVec &cFBV
		);
	CFixBitVec &
	operator^=(
		FBV_UINT i
		);

	// XOR operation.
	CFixBitVec
	operator^(
		const CFixBitVec &cFBV
		) const;
	CFixBitVec
	operator^(
		FBV_UINT i
		);

	// Shift left operation, in place.
	CFixBitVec &
	operator<<=(
		int iBits
		);

	// Shift left operation.
	CFixBitVec
	operator<<(
		int iBits
		) const;
	
	// Shift right operation, in place.
	CFixBitVec &
	operator>>=(
		int iBits
		);

	// Shift right operation.
	CFixBitVec
	operator>>(
		int iBits
		) const;

	// Right rotation, in place.
	CFixBitVec &
	rotr(
		int iBits,
		int iWidth = FBV_BITS
		);

	// Right rotation.
	CFixBitVec
	rotrCopy(
		int iBits,
		int iWidth = FBV_BITS
		) const;

	// Left rotation, in place.
	CFixBitVec &
	rotl(
		int iBits,
		int iWidth = FBV_BITS
		);

	// Left rotation.
	CFixBitVec
	rotlCopy(
		int iBits,
		int iWidth = FBV_BITS
		) const;

	// Is the bit rack zero valued?
	bool
	isZero() const;

	// Returns the number of trailing set bits.
	int
	tsb() const;

	// Returns the index of the most significant bit
	int
	msb() const;
                
	// Returns the index of the first set bit, numbered from
	// 1 to n.  0 means there were no set bits.
	int
	fsb() const;

	// Prefix decrement.  Returns true if a carry
	// was generated.
	bool
	operator--();

	// Calculates the Gray Code.
	CFixBitVec &
	grayCode();

	// Calculates the Gray Code Inverse
	CFixBitVec &
	grayCodeInv();

	// Ones-complements the rack
	CFixBitVec &
	complement();

	// Returns the first rack.
	FBV_UINT &
	rack();
	FBV_UINT
	rack() const;

	// Return a pointer to the racks
	FBV_UINT *
	racks();
	const FBV_UINT *
	racks() const;

	// Returns the number of racks.
	int
	rackCount();

        // Comparison operator
        bool
	operator<(const CFixBitVec& other) const
        { return this->m_uiRack < other.m_uiRack; }

        // Comparison operator
        bool
	operator==(const CFixBitVec& other) const
        { return this->m_uiRack == other.m_uiRack; }
  
  private:
	
	static
	void
	compileTimeAssertions();

	FBV_UINT	m_uiRack;
};


#endif
