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


#include <Hilbert/FixBitVec.hpp>


// This must be an unsigned integer that is either
// 32 or 64 bits.  Otherwise, there are places in the
// code that simply will not work.
// For speed, this should be the native word size.
//typedef uint64_t FBV_UINT;
//#define FBV_BITS		64
typedef uint32_t FBV_UINT;
#define FBV_BITS		32


#define COMPILE_TIME_ASSERT(pred) switch(0){case 0:case pred:;}


EBitVecType
CFixBitVec::type()
{
	return eFix;
}

// Default constructor.  The bits parameter
// is completely ignored, but accepted in order
// to look and feel the same as a BigBitVec.
CFixBitVec::CFixBitVec(
		       int /* iBits */
		       )
{
	return;
}

// Copy constructor.
CFixBitVec::CFixBitVec(
	const CFixBitVec &cFBV
	)
{
	m_uiRack = cFBV.m_uiRack;
}

// Returns the current size in bits.
int
CFixBitVec::getSize()
{
	return FBV_BITS;
}

// Sets the size.  This is a dummy
// function just for BigBitVec compatibility.
CFixBitVec &
CFixBitVec::setSize(
		    int /* iBits */
		    )
{
	return (*this);
}

// Zeros the bit-vector.
CFixBitVec &
CFixBitVec::zero()
{
	m_uiRack = 0;
	return (*this);
}

// Truncates the bit-vector to a given precision in
// bits (zeroes MSBs without shrinking the vector)
CFixBitVec &
CFixBitVec::truncate(
	int iBits
	)
{
	assert( 0 <= iBits && iBits <= FBV_BITS );
	m_uiRack &= FBVN1S(iBits);
	return (*this);
}

// Assignment operator.
CFixBitVec &
CFixBitVec::operator=(
	const CFixBitVec &cFBV
	)
{
	m_uiRack = cFBV.m_uiRack;
	return (*this);
}

// Assignment operator.
CFixBitVec &
CFixBitVec::operator=(
	FBV_UINT i
	)
{
	m_uiRack = i;
	return (*this);
}

// Returns the value of the nth bit.
bool
CFixBitVec::getBit(
	int iIndex
	) const
{
	assert( 0 <= iIndex && iIndex < FBV_BITS );
	return ( (m_uiRack & (FBV1<<iIndex)) > 0 );
}

// Sets the value of the nth bit.
CFixBitVec &
CFixBitVec::setBit(
	int iIndex,
	bool bBit
	)
{
	assert( 0 <= iIndex && iIndex < FBV_BITS );
	FBV_UINT m = (FBV1<<iIndex);
	m_uiRack |= m;
	if ( !bBit ) m_uiRack ^= m;
	return (*this);
}

// Toggles the value of the nth bit.
CFixBitVec &
CFixBitVec::toggleBit(
	int iIndex
	)
{
	assert( 0 <= iIndex && iIndex < FBV_BITS );
	m_uiRack ^= (FBV1<<iIndex);
	return (*this);
}

// AND operation in place.
CFixBitVec &
CFixBitVec::operator&=(
	const CFixBitVec &cFBV
	)
{
	m_uiRack &= cFBV.m_uiRack;
	return (*this);
}
CFixBitVec &
CFixBitVec::operator&=(
	FBV_UINT i
	)
{
	m_uiRack &= i;
	return (*this);
}

// AND operation.
CFixBitVec
CFixBitVec::operator&(
	const CFixBitVec &cFBV
	) const
{
	CFixBitVec t(*this);
	t &= cFBV;
	return t;
}
CFixBitVec
CFixBitVec::operator&(
	FBV_UINT i
	)
{
	CFixBitVec t(*this);
	t &= i;
	return t;
}

// OR operation in place.
CFixBitVec &
CFixBitVec::operator|=(
	const CFixBitVec &cFBV
	)
{
	m_uiRack |= cFBV.m_uiRack;
	return (*this);
}
CFixBitVec &
CFixBitVec::operator|=(
	FBV_UINT i
	)
{
	m_uiRack |= i;
	return (*this);
}

// OR operation.
CFixBitVec
CFixBitVec::operator|(
	const CFixBitVec &cFBV
	) const
{
	CFixBitVec t(*this);
	t |= cFBV;
	return t;
}
CFixBitVec
CFixBitVec::operator|(
	FBV_UINT i
	)
{
	CFixBitVec t(*this);
	t |= i;
	return t;
}


// XOR operation in place.
CFixBitVec &
CFixBitVec::operator^=(
	const CFixBitVec &cFBV
	)
{
	m_uiRack ^= cFBV.m_uiRack;
	return (*this);
}
CFixBitVec &
CFixBitVec::operator^=(
	FBV_UINT i
	)
{
	m_uiRack ^= i;
	return (*this);
}

// XOR operation.
CFixBitVec
CFixBitVec::operator^(
	const CFixBitVec &cFBV
	) const
{
	CFixBitVec t(*this);
	t ^= cFBV;
	return t;
}
CFixBitVec
CFixBitVec::operator^(
	FBV_UINT i
	)
{
	CFixBitVec t(*this);
	t ^= i;
	return t;
}

// Shift left operation, in place.
CFixBitVec &
CFixBitVec::operator<<=(
	int iBits
	)
{
	m_uiRack <<= iBits;
	return (*this);
}

// Shift left operation.
CFixBitVec
CFixBitVec::operator<<(
	int iBits
	) const
{
	CFixBitVec t(*this);
	t <<= iBits;
	return t;
}

// Shift right operation, in place.
CFixBitVec &
CFixBitVec::operator>>=(
	int iBits
	)
{
	m_uiRack >>= iBits;
	return (*this);
}

// Shift right operation.
CFixBitVec
CFixBitVec::operator>>(
	int iBits
	) const
{
	CFixBitVec t(*this);
	t >>= iBits;
	return t;
}

// Right rotation, in place.
CFixBitVec &
CFixBitVec::rotr(
	int iBits,
	int iWidth
	)
{
	assert( iBits >= 0 );
	assert( iWidth >  0 );
	assert( iBits < iWidth );//#D FBVMOD(iBits,iWidth);
	m_uiRack &= FBVN1S(iWidth);
	m_uiRack = (m_uiRack>>iBits) | (m_uiRack<<(iWidth-iBits));
	m_uiRack &= FBVN1S(iWidth);
	return (*this);
}

// Right rotation.
CFixBitVec
CFixBitVec::rotrCopy(
	int iBits,
	int iWidth
	) const
{
	CFixBitVec t(*this);
	t.rotr(iBits,iWidth);
	return t;
}

// Left rotation, in place.
CFixBitVec &
CFixBitVec::rotl(
	int iBits,
	int iWidth
	)
{
	assert( iBits >= 0 );
	assert( iWidth > 0 );
	assert( iBits < iWidth );//#D FBVMOD(iBits,iWidth);
	m_uiRack &= FBVN1S(iWidth);
	m_uiRack = (m_uiRack<<iBits) | (m_uiRack>>(iWidth-iBits));
	m_uiRack &= FBVN1S(iWidth);
	return (*this);
}

// Left rotation.
CFixBitVec
CFixBitVec::rotlCopy(
	int iBits,
	int iWidth
	) const
{
	CFixBitVec t(*this);
	t.rotl(iBits,iWidth);
	return t;
}

// Is the bit rack zero valued?
bool
CFixBitVec::isZero() const
{
	return m_uiRack == 0;
}

// Returns the number of trailing set bits.
int
CFixBitVec::tsb() const
{
	FBV_UINT i = m_uiRack;
	int c = 0;

	#if FBV_BITS == 64
	if ( i == FBV1S ) return 64;
	if ( (i&FBVN1S(32)) == FBVN1S(32) ) { i>>=32; c^=32; }
	#elif FBV_BITS == 32
	if ( i == FBV1S ) return 32;
	#endif
	if ( (i&FBVN1S(16)) == FBVN1S(16) ) { i>>=16; c^=16; }
	if ( (i&FBVN1S( 8)) == FBVN1S( 8) ) { i>>= 8; c^= 8; }
	if ( (i&FBVN1S( 4)) == FBVN1S( 4) ) { i>>= 4; c^= 4; }
	if ( (i&FBVN1S( 2)) == FBVN1S( 2) ) { i>>= 2; c^= 2; }
	if ( (i&FBVN1S( 1)) == FBVN1S( 1) ) { i>>= 1; c^= 1; }

	return c;	
}

// Returns the index of the most significant bit
int
CFixBitVec::msb() const
{
	FBV_UINT i = m_uiRack;
	int c = 0;
	#if FBV_BITS == 64
	if ( i == FBV0 ) return 0;
	if ( i & (FBVN1S(32) << 32) ) { i >>= 32; c ^= 32; }
	#elif FBV_BITS == 32
	if ( i == FBV0 ) return 0;
	#endif
	if ( i & (FBVN1S(16) << 16) ) { i >>= 16; c ^= 16; }
	if ( i & (FBVN1S( 8) <<  8) ) { i >>=  8; c ^=  8; }
	if ( i & (FBVN1S( 4) <<  4) ) { i >>=  4; c ^=  4; }
	if ( i & (FBVN1S( 2) <<  2) ) { i >>=  2; c ^=  2; }
	if ( i & (FBVN1S( 1) <<  1) ) { i >>=  1; c ^=  1; }
	return ++c;
}
                
// Returns the index of the first set bit, numbered from
// 1 to n.  0 means there were no set bits.
int
CFixBitVec::fsb() const
{
	FBV_UINT i = m_uiRack;
	int c = 0;

	#if FBV_BITS == 64
	if ( i == FBV0 ) return 0;
	if ( (i&FBVN1S(32)) == FBV0 ) { i>>=32; c^=32; }
	#elif FBV_BITS == 32
	if ( i == FBV0 ) return 0;
	#endif
	if ( (i&FBVN1S(16)) == FBV0 ) { i>>=16; c^=16; }
	if ( (i&FBVN1S( 8)) == FBV0 ) { i>>= 8; c^= 8; }
	if ( (i&FBVN1S( 4)) == FBV0 ) { i>>= 4; c^= 4; }
	if ( (i&FBVN1S( 2)) == FBV0 ) { i>>= 2; c^= 2; }
	if ( (i&FBVN1S( 1)) == FBV0 ) { i>>= 1; c^= 1; }

	return ++c;
}

// Prefix decrement.  Returns true if a carry
// was generated.
bool
CFixBitVec::operator--()
{
	return ( m_uiRack-- == 0 );
}

// Calculates the Gray Code.
CFixBitVec &
CFixBitVec::grayCode()
{
	m_uiRack ^= (m_uiRack>>1);
	return (*this);
}

// Calculates the Gray Code Inverse
CFixBitVec &
CFixBitVec::grayCodeInv()
{
	m_uiRack ^= (m_uiRack>>1);
	m_uiRack ^= (m_uiRack>>2);
	m_uiRack ^= (m_uiRack>>4);
	m_uiRack ^= (m_uiRack>> 8);
	m_uiRack ^= (m_uiRack>>16);
	#if FBV_BITS == 64
	m_uiRack ^= (m_uiRack>>32);
	#endif
	return (*this);
}

// Ones-complements the rack
CFixBitVec &
CFixBitVec::complement()
{
	m_uiRack = ~m_uiRack;
	return (*this);
}

// Returns the first rack.
FBV_UINT &
CFixBitVec::rack()
{
	return m_uiRack;
}
FBV_UINT
CFixBitVec::rack() const
{
	return m_uiRack;
}

// Return a pointer to the racks
FBV_UINT *
CFixBitVec::racks()
{
	return &m_uiRack;
}
const FBV_UINT *
CFixBitVec::racks() const
{
	return &m_uiRack;
}

// Returns the number of racks.
int
CFixBitVec::rackCount()
{
	return 1;
}

void
CFixBitVec::compileTimeAssertions()
{
	COMPILE_TIME_ASSERT( 8*sizeof(FBV_UINT) == FBV_BITS );
	COMPILE_TIME_ASSERT( (sizeof(FBV_UINT) == 4) || (sizeof(FBV_UINT) == 8) );
}
