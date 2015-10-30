/*
 * Copyright (C) 2006-2014 Chris Hamilton <chamilton@cs.dal.ca>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#ifndef _BIGBITVEC_HPP_
#define _BIGBITVEC_HPP_


#include "Hilbert/Common.hpp"
#include "Hilbert/FixBitVec.hpp"

#define BBV_MIN(a,b)		((a)<(b)?(a):(b))
#define BBV_MAX(a,b)		((a)>(b)?(a):(b))
#define FBVS_NEEDED(b)		((BBV_MAX(b,1)+FBV_BITS-1)/FBV_BITS)
#define BBV_MODSPLIT(r,b,k) { b=(k); r=b/FBV_BITS; b-=r*FBV_BITS; }

class CBigBitVec
{
  public:

	static
	EBitVecType
	type();

	// Constructor, with optional number of bits.
	CBigBitVec(
		int iBits = FBV_BITS
		);

	// Copy construct.  Creates duplicate.
	CBigBitVec(
		const CBigBitVec &cBBV
		);

	// Copy constructor.
	CBigBitVec(
		const CFixBitVec &cFBV
		);

	// Destructor
	~CBigBitVec();

	// Returns the current size in bits.
	int
	getSize() const;

	// Resize function.  Returns the number of bits
	// we can accomodate after resizing.
	CBigBitVec &
	setSize(
		int iBits
		);

	// Zeros the bit-vector.
	CBigBitVec &
	zero();

	// Truncates the bit-vector to a given precision in
	// bits (zeroes MSBs without shrinking the vector)
	CBigBitVec &
	truncate(
		int iBits
		);

	// Assignment operator.  No resizing.
	CBigBitVec &
	operator=(
		const CBigBitVec &cBBV
		);
	CBigBitVec &
	operator=(
		const CFixBitVec &cFBV
		);
	CBigBitVec &
	operator=(
		 FBV_UINT j
		);

	// Returns the value of the nth bit.
	bool
	getBit(
		int iIndex
		) const;

	// Sets the value of the nth bit.
	CBigBitVec &
	setBit(
		int iIndex,
		bool bBit
		);

	// Toggles the value of the nth bit.
	CBigBitVec &
	toggleBit(
		int iIndex
		);

	// In place AND.
	CBigBitVec &
	operator&=(
		const CBigBitVec &cBBV
		);
	CBigBitVec &
	operator&=(
		const CFixBitVec &r
		);
	CBigBitVec &
	operator&=(
		 FBV_UINT i
		);

	// AND operator.
	CBigBitVec
	operator&(
		const CBigBitVec &cBBV
		) const;
	CBigBitVec
	operator&(
		const CFixBitVec &r
		);
	CBigBitVec
	operator&(
		 FBV_UINT i
		);

	// In place OR.
	CBigBitVec &
	operator|=(
		const CBigBitVec &cBBV
		);
	CBigBitVec &
	operator|=(
		const CFixBitVec &r
		);
	CBigBitVec &
	operator|=(
		 FBV_UINT i
		);

	// OR operator.
	CBigBitVec
	operator|(
		const CBigBitVec &cBBV
		) const;
	CBigBitVec
	operator|(
		const CFixBitVec &r
		);
	CBigBitVec
	operator|(
		 FBV_UINT i
		);

	// In place XOR.
	CBigBitVec &
	operator^=(
		const CBigBitVec &cBBV
		);
	CBigBitVec &
	operator^=(
		const CFixBitVec &r
		);
	CBigBitVec &
	operator^=(
		FBV_UINT i
		);

	// XOR operator.
	CBigBitVec
	operator^(
		const CBigBitVec &cBBV
		) const;
	CBigBitVec
	operator^(
		const CFixBitVec &r
		);
	CBigBitVec
	operator^(
		 FBV_UINT i
		);

	// Shift left operation, in place.
	CBigBitVec &
	operator<<=(
		int iBits
		);

	// Shift left operation.
	CBigBitVec
	operator<<(
		int iBits
		) const;

	// Shift right operation, in place.
	CBigBitVec &
	operator>>=(
		int iBits
		);

	// Shift right operation.
	CBigBitVec
	operator>>(
		int iBits
		) const;

	// Right rotation, in place.
	CBigBitVec &
	rotr(
		int iBits,
		int iWidth = 0
		);

	// Right rotation.
	CBigBitVec
	rotrCopy(
		int iBits,
		int iWidth = 0
		) const;

	// Left rotation, in place.
	CBigBitVec &
	rotl(
		int iBits,
		int iWidth = 0
		);

	// Left rotation.
	CBigBitVec
	rotlCopy(
		int iBits,
		int iWidth = 0
		) const;

	// Returns true if the rack is zero valued.
	bool
	isZero() const;

	// Returns the number of trailing set bits.
	int
	tsb() const;

	// OB:
	// Returns the index of the most significant bit (numbered
	// 1 to n)
	int
	msb() const;

	// Returns the index of the first set bit.
	// (numbered 1 to n, with 0 meaning no bits were set)
	int
	fsb() const;

	// Prefix decrement.  Returns true if there
	// was a carry, false otherwise.
	bool
	operator--();

	// Gray Code
	CBigBitVec &
	grayCode();

	// Gray Code Inverse
	CBigBitVec &
	grayCodeInv();

	// Complement
	CBigBitVec &
	complement();

	// Returns the first rack.
	FBV_UINT &
	rack();
	FBV_UINT
	rack() const;

	// Returns the racks.
	FBV_UINT *
	racks();
	const FBV_UINT *
	racks() const;

	// Returns the number of racks.
	int
	rackCount() const;

        // Comparison operator
        bool
	operator<(const CBigBitVec& other) const;

        // Comparison operators
        bool
	operator==(const CBigBitVec& other) const;
        bool
        operator!=(const CBigBitVec& other) const;

  private:

	// Right rotates entire racks (in place).
	void
	rackRotr(
		int k
		);

	CFixBitVec *m_pcRacks;
	int m_iRacks;
};


#endif
