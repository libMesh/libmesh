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

#include "Hilbert/BigBitVec.hpp"
#include <string.h>


EBitVecType
CBigBitVec::type()
{
	return eBig;
}

// Constructor, with optional number of bits.
CBigBitVec::CBigBitVec(
	int iBits
	)
{
	// Determine number of racks required.
	m_iRacks = FBVS_NEEDED(iBits);

	// Allocate the memory.
	m_pcRacks = new CFixBitVec[m_iRacks];

	return;
}


// Copy construct.	Creates duplicate.
CBigBitVec::CBigBitVec(
	const CBigBitVec &cBBV
	)
{
	m_iRacks = cBBV.m_iRacks;
	m_pcRacks = new CFixBitVec[m_iRacks];

	// Copy the rack values.
	/*int i;
	for ( i = 0; i < m_iRacks; i++ )
		m_pcRacks[i] = cBBV.m_pcRacks[i];*/
	memcpy( static_cast<void*>(m_pcRacks),
		static_cast<const void*>(cBBV.m_pcRacks),
		sizeof(CFixBitVec)*m_iRacks );

	return;
}


// Copy constructor.
CBigBitVec::CBigBitVec(
	const CFixBitVec &cFBV
	)
{
	m_iRacks = 1;
	m_pcRacks = new CFixBitVec[m_iRacks];

	m_pcRacks[0] = cFBV;

	return;
}


// Destructor
CBigBitVec::~CBigBitVec()
{
	delete [] m_pcRacks;

	return;
}


// Returns the current size in bits.
int
CBigBitVec::getSize() const
{
	return m_iRacks*FBV_BITS;
}


// Resize function.	Returns the number of bits
// we can accomodate after resizing.
CBigBitVec &
CBigBitVec::setSize(
	int iBits
	)
{
	int i;

	// How many racks do we need?
	int iRacks = FBVS_NEEDED(iBits);

	// Same size?
	if ( iRacks == m_iRacks )
		return (*this);

	// Allocate new racks.
	CFixBitVec *pcRacks = new CFixBitVec[iRacks];

	// Copy over the old values.
	/*for ( i = 0; i < BBV_MIN(iRacks,m_iRacks); i++ )
		pcRacks[i] = m_pcRacks[i];*/
	i = BBV_MIN(iRacks,m_iRacks);
	memcpy( static_cast<void*>(pcRacks),
		static_cast<const void*>(m_pcRacks),
		sizeof(CFixBitVec)*i );

	// zero the new values
	/*for ( ; i < iRacks; i++ )
		pcRacks[i].zero();*/
	if ( iRacks > i )
		memset( static_cast<void*>(pcRacks + i), 0,
			sizeof(CFixBitVec)*(iRacks-i) );

	// Release the old stuff.
	delete [] m_pcRacks;

	// Replace old with new.
	m_iRacks = iRacks;
	m_pcRacks = pcRacks;

	return (*this);
}


// zeros the bit-vector.
CBigBitVec &
CBigBitVec::zero()
{

	/*int i;
	for ( i = 0; i < m_iRacks; i++ )
		m_pcRacks[i].zero();*/
	memset( static_cast<void*>(m_pcRacks), 0,
		sizeof(CFixBitVec)*m_iRacks );

	return (*this);
}


// truncates the bit-vector to a given precision in
// bits (zeroes MSBs without shrinking the vector)
CBigBitVec &
CBigBitVec::truncate(
	int iBits
	)
{
	assert( iBits >= 0 && iBits <= getSize() );
	int r, b, i;

	BBV_MODSPLIT(r,b,iBits);

	if ( r >= m_iRacks )
		return (*this);

	m_pcRacks[r].truncate(b);

	for ( i = r+1; i < m_iRacks; i++ )
		m_pcRacks[i].zero();

	return (*this);
}


// Assignment operator.	No resizing.
CBigBitVec &
CBigBitVec::operator=(
	const CBigBitVec &cBBV
	)
{
	if ( m_iRacks < cBBV.m_iRacks )
	{
		/*int i;
		for ( i = 0; i < m_iRacks; i++ )
			m_pcRacks[i] = cBBV.m_pcRacks[i];*/
		memcpy( static_cast<void*>(m_pcRacks),
			static_cast<const void*>(cBBV.m_pcRacks),
			sizeof(CFixBitVec)*m_iRacks );
	}
	else
	{
		/*int i;
		for ( i = 0; i < cBBV.m_iRacks; i++ )
			m_pcRacks[i] = cBBV.m_pcRacks[i];
		for ( ; i < m_iRacks; i++ )
			m_pcRacks[i].zero();*/
		memcpy( static_cast<void*>(m_pcRacks),
			static_cast<const void*>(cBBV.m_pcRacks),
			sizeof(CFixBitVec)*cBBV.m_iRacks );
		memset( static_cast<void*>(m_pcRacks+cBBV.m_iRacks),
			0, sizeof(CFixBitVec)*(m_iRacks-cBBV.m_iRacks) );
	}

	return (*this);
}
CBigBitVec &
CBigBitVec::operator=(
	const CFixBitVec &cFBV
	)
{
	m_pcRacks[0] = cFBV;
	/*int i;
	for ( i = 1; i < m_iRacks; i++ )
		m_pcRacks[i].zero();*/
	memset( static_cast<void*>(m_pcRacks+1),
		0, sizeof(CFixBitVec)*(m_iRacks-1) );
	return (*this);
}
CBigBitVec &
CBigBitVec::operator=(
	 FBV_UINT j
	)
{
	m_pcRacks[0] = j;
	/*int i;
	for ( i = 1; i < m_iRacks; i++ )
		m_pcRacks[i].zero();*/
	memset( static_cast<void*>(m_pcRacks+1),
		0, sizeof(CFixBitVec)*(m_iRacks-1) );
	return (*this);
}

// Returns the value of the nth bit.
bool
CBigBitVec::getBit(
	int iIndex
	) const
{
	assert( iIndex >= 0 && iIndex < getSize() );
	int r, b;
	BBV_MODSPLIT(r,b,iIndex);
	return m_pcRacks[r].getBit(b);
}


// Sets the value of the nth bit.
CBigBitVec &
CBigBitVec::setBit(
	int iIndex,
	bool bBit
	)
{
	assert( iIndex >= 0 && iIndex < getSize() );
	int r, b;
	BBV_MODSPLIT(r,b,iIndex);
	m_pcRacks[r].setBit(b,bBit);
	return (*this);
}


// Toggles the value of the nth bit.
CBigBitVec &
CBigBitVec::toggleBit(
	int iIndex
	)
{
	assert( iIndex >= 0 && iIndex < getSize() );
	int r, b;
	BBV_MODSPLIT(r,b,iIndex);
	m_pcRacks[r].toggleBit(b);
	return (*this);
}


// In place AND.
CBigBitVec &
CBigBitVec::operator&=(
	const CBigBitVec &cBBV
	)
{
	int i;

	for ( i = 0; i < BBV_MIN(m_iRacks,cBBV.m_iRacks); i++ )
		m_pcRacks[i] &= cBBV.m_pcRacks[i];

	return (*this);
}
CBigBitVec &
CBigBitVec::operator&=(
	const CFixBitVec &r
	)
{
	m_pcRacks[0] &= r;
	return (*this);
}
CBigBitVec &
CBigBitVec::operator&=(
	 FBV_UINT i
	)
{
	m_pcRacks[0] &= i;
	return (*this);
}


// AND operator.
CBigBitVec
CBigBitVec::operator&(
	const CBigBitVec &cBBV
	) const
{
	CBigBitVec t( *this );
	t &= cBBV;

	return t;
}
CBigBitVec
CBigBitVec::operator&(
	const CFixBitVec &r
	)
{
	CBigBitVec t( *this );
	t &= r;

	return t;
}
CBigBitVec
CBigBitVec::operator&(
	 FBV_UINT i
	)
{
	CBigBitVec t( *this );
	t &= i;

	return t;
}


// In place OR.
CBigBitVec &
CBigBitVec::operator|=(
	const CBigBitVec &cBBV
	)
{
	int i;

	for ( i = 0; i < BBV_MIN(m_iRacks,cBBV.m_iRacks); i++ )
		m_pcRacks[i] |= cBBV.m_pcRacks[i];

	return (*this);
}
CBigBitVec &
CBigBitVec::operator|=(
	const CFixBitVec &r
	)
{
	m_pcRacks[0] |= r;
	return (*this);
}
CBigBitVec &
CBigBitVec::operator|=(
	 FBV_UINT i
	)
{
	m_pcRacks[0] |= i;
	return (*this);
}


// OR operator.
CBigBitVec
CBigBitVec::operator|(
	const CBigBitVec &cBBV
	) const
{
	CBigBitVec t( *this );
	t |= cBBV;

	return t;
}
CBigBitVec
CBigBitVec::operator|(
	const CFixBitVec &r
	)
{
	CBigBitVec t( *this );
	t |= r;

	return t;
}
CBigBitVec
CBigBitVec::operator|(
	 FBV_UINT i
	)
{
	CBigBitVec t( *this );
	t |= i;

	return t;
}


// In place XOR.
CBigBitVec &
CBigBitVec::operator^=(
	const CBigBitVec &cBBV
	)
{
	int i;

	for ( i = 0; i < BBV_MIN(m_iRacks,cBBV.m_iRacks); i++ )
		m_pcRacks[i] ^= cBBV.m_pcRacks[i];

	return (*this);
}
CBigBitVec &
CBigBitVec::operator^=(
	const CFixBitVec &r
	)
{
	m_pcRacks[0] ^= r;
	return (*this);
}
CBigBitVec &
CBigBitVec::operator^=(
	 FBV_UINT i
	)
{
	m_pcRacks[0] ^= i;
	return (*this);
}


// XOR operator.
CBigBitVec
CBigBitVec::operator^(
	const CBigBitVec &cBBV
	) const
{
	CBigBitVec t( *this );
	t ^= cBBV;

	return t;
}
CBigBitVec
CBigBitVec::operator^(
	const CFixBitVec &r
	)
{
	CBigBitVec t( *this );
	t ^= r;

	return t;
}
CBigBitVec
CBigBitVec::operator^(
	 FBV_UINT i
	)
{
	CBigBitVec t( *this );
	t ^= i;

	return t;
}


// Shift left operation, in place.
CBigBitVec &
CBigBitVec::operator<<=(
	int iBits
	)
{
	assert( iBits >= 0 );
	int r, b, i;

	// No shift?
	if ( iBits == 0 )
		return (*this);

	BBV_MODSPLIT(r,b,iBits);

	// All racks?
	if ( r >= m_iRacks )
	{
		zero();
		return (*this);
	}

	// Do rack shifts.
	if ( r > 0 )
	{
		for ( i = m_iRacks-1; i >= r; i-- )
			m_pcRacks[i] = m_pcRacks[i-r];
		for ( ; i >= 0; i-- )
			m_pcRacks[i].zero();
	}

	// Do bit shifts.
	if ( b > 0 )
	{
		int bi = FBV_BITS - b;
		for ( i = m_iRacks-1; i >= r+1; i-- )
		{
			m_pcRacks[i] <<= b;
			m_pcRacks[i] |= m_pcRacks[i-1] >> bi;
		}
		m_pcRacks[i] <<= b;
	}

	return (*this);
}


// Shift left operation.
CBigBitVec
CBigBitVec::operator<<(
	int iBits
	) const
{
	CBigBitVec t( *this );
	t <<= iBits;
	return t;
}


// Shift right operation, in place.
CBigBitVec &
CBigBitVec::operator>>=(
	int iBits
	)
{
	assert( iBits >= 0 );
	int r, b, i;

	// No shift?
	if ( iBits == 0 )
		return (*this);

	BBV_MODSPLIT(r,b,iBits);

	// All racks?
	if ( r >= m_iRacks )
	{
		zero();
		return (*this);
	}

	// Do rack shifts.
	if ( r > 0 )
	{
		for ( i = 0; i < m_iRacks-r; i++ )
			m_pcRacks[i] = m_pcRacks[i+r];
		for ( ; i < m_iRacks; i++ )
			m_pcRacks[i].zero();
	}

	// Do bit shifts.
	if ( b > 0 )
	{
		int bi = FBV_BITS - b;
		for ( i = 0; i < m_iRacks-r-1; i++ )
		{
			m_pcRacks[i] >>= b;
			m_pcRacks[i] |= m_pcRacks[i+1] << bi;
		}
		m_pcRacks[i] >>= b;
	}

	return (*this);
}


// Shift right operation.
CBigBitVec
CBigBitVec::operator>>(
	int iBits
	) const
{
	CBigBitVec t( *this );
	t >>= iBits;
	return t;
}


// Right rotation, in place.
CBigBitVec &
CBigBitVec::rotr(
	int iBits,
	int iWidth
	)
{
	assert( iBits >= 0 );

	// Fill in the width, if necessary.
	if ( iWidth <= 0 )
		iWidth = getSize();

	// Modulo the number of bits.
	//FBVMOD(iBits,iWidth);
	assert( iBits < iWidth );
	if ( iBits == 0 ) return (*this);

	// Ensure we are truncated appropriately.
	truncate(iWidth);

	CBigBitVec t1( *this );
	(*this) >>= iBits;
	t1 <<= (iWidth-iBits);
	(*this) |= t1;

	truncate(iWidth);

	return (*this);
}


// Right rotation.
CBigBitVec
CBigBitVec::rotrCopy(
	int iBits,
	int iWidth
	) const
{
	CBigBitVec t( *this );
	t.rotr(iBits,iWidth);
	return t;
}


// Left rotation, in place.
CBigBitVec &
CBigBitVec::rotl(
	int iBits,
	int iWidth
	)
{
	assert( iBits >= 0 );

	// Fill in the width, if necessary.
	if ( iWidth <= 0 )
		iWidth = getSize();

	// Modulo the number of bits.
	//FBVMOD(iBits,iWidth);
	assert( iBits < iWidth );
	if ( iBits == 0 ) return (*this);

	// Ensure we are truncated appropriately.
	truncate(iWidth);

	CBigBitVec t1( *this );
	(*this) <<= iBits;
	t1 >>= (iWidth-iBits);
	(*this) |= t1;

	truncate(iWidth);

	return (*this);
}


// Left rotation.
CBigBitVec
CBigBitVec::rotlCopy(
	int iBits,
	int iWidth
	) const
{
	CBigBitVec t( *this );
	t.rotl(iBits,iWidth);
	return t;
}


// Returns true if the rack is zero valued.
bool
CBigBitVec::isZero() const
{
	int i;
	for ( i = 0; i < m_iRacks; i++ )
		if ( !m_pcRacks[i].isZero() ) return false;
	return true;
}


// Returns the number of trailing set bits.
int
CBigBitVec::tsb() const
{
	int c, i, j;
	c = 0;
	for ( i = 0; i < m_iRacks; i++ )
	{
		j = m_pcRacks[i].tsb();
		c += j;
		if ( j < FBV_BITS )
			break;
	}
	return c;
}

// OB:
// Returns the index of the most significant bit (numbered
// 1 to n)
int
CBigBitVec::msb() const
{
	int c, i, j = 0;
	c = FBV_BITS * m_iRacks;
	for ( i = m_iRacks - 1; i >= 0 && !j; i-- )
	{
		j = m_pcRacks[i].msb();
		if (j)
			return c - (FBV_BITS - j);
		else
			c -= FBV_BITS;
	}
	return 0;
}

// Returns the index of the first set bit.
// (numbered 1 to n, with 0 meaning no bits were set)
int
CBigBitVec::fsb() const
{
	int c, i, j;
	c = 0;
	for ( i = 0; i < m_iRacks; i++ )
	{
		j = m_pcRacks[i].fsb();
		if ( j < FBV_BITS )
			return c + j;
		else
			c += FBV_BITS;
	}
	return 0;
}


// Prefix decrement.	Returns true if there
// was a carry, false otherwise.
bool
CBigBitVec::operator--()
{
	int i = 0;
	bool b=false;
	while ( i < m_iRacks && (b = --m_pcRacks[i]) ) i++;

	return b;
}


// Gray Code
CBigBitVec &
CBigBitVec::grayCode()
{
	int i;
	FBV_UINT s = 0;

	for ( i = m_iRacks-1; i >= 0; i-- )
	{
		FBV_UINT t = m_pcRacks[i].rack() & 1;
		m_pcRacks[i].grayCode();
		m_pcRacks[i] ^= (s<<(FBV_BITS-1));
		s = t;
	}

	return (*this);
}


// Gray Code Inverse
CBigBitVec &
CBigBitVec::grayCodeInv()
{
	int i;
	FBV_UINT s = 0;

	for ( i = m_iRacks-1; i >= 0; i-- )
	{
		m_pcRacks[i].grayCodeInv();
		if ( s ) m_pcRacks[i].complement();
		s = m_pcRacks[i].rack() & 1;
	}
	return (*this);
}


// Complement
CBigBitVec &
CBigBitVec::complement()
{
	int i;
	for ( i = 0; i < m_iRacks; i++ )
		m_pcRacks[i].complement();
	return (*this);
}


// Returns the first rack.
FBV_UINT &
CBigBitVec::rack()
{
	return m_pcRacks[0].rack();
}
FBV_UINT
CBigBitVec::rack() const
{
	return m_pcRacks[0].rack();
}


// Returns the racks.
FBV_UINT *
CBigBitVec::racks()
{
	return reinterpret_cast<FBV_UINT*>(m_pcRacks);
}
const FBV_UINT *
CBigBitVec::racks() const
{
	return reinterpret_cast<FBV_UINT*>(m_pcRacks);
}


// Returns the number of racks.
int
CBigBitVec::rackCount() const
{
	return m_iRacks;
}


// Right rotates entire racks (in place).
void
CBigBitVec::rackRotr(
	int k
	)
{
	assert( 0 <= k && k < m_iRacks );

	int c, v;
	CFixBitVec tmp;

	if (k == 0) return;

	c = 0;
	for (v = 0; c < m_iRacks; v++)
	{
		int t = v, tp = v + k;
		tmp = m_pcRacks[v];
		c++;
		while (tp != v)
		{
			m_pcRacks[t] = m_pcRacks[tp];
			t = tp;
			tp += k;
			if (tp >= m_iRacks) tp -= m_iRacks;
			c++;
		}
		m_pcRacks[t] = tmp;
	}

	return;
}


// Comparison operator
bool
CBigBitVec::operator<(const CBigBitVec& b) const
{
  const CBigBitVec &a = *this;

  int rack;

  if (a.rackCount() > b.rackCount())
    {
      for (rack = (a.rackCount()-1); rack >= b.rackCount(); rack--)
	if (a.racks()[rack] > 0)
	  return false;
    }

  else if (a.rackCount() < b.rackCount())
    {
      for (rack = (b.rackCount()-1); rack >= a.rackCount(); rack--)
	if (0 < b.racks()[rack])
	  return true;
    }
  else
    rack = a.rackCount()-1;

  while (rack > 0 && a.racks()[rack] == b.racks()[rack]) rack--;

  return a.racks()[rack] < b.racks()[rack];
}


// Comparison operators
bool
CBigBitVec::operator==(const CBigBitVec& b) const
{
  const CBigBitVec &a = *this;

  if (a.rackCount() != b.rackCount())
    return false;

  int rack = a.rackCount()-1;

  while (rack > 0 && a.racks()[rack] == b.racks()[rack]) rack--;

  return a.racks()[rack] == b.racks()[rack];
}


bool
CBigBitVec::operator!=(const CBigBitVec& b) const
{
  return !(*this == b);
}
