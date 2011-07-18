#include "Hilbert.hpp"
#include "Hilbert/Algorithm.hpp"


namespace Hilbert
{

// fix -> fix

void coordsToIndex( const CFixBitVec *p, int m, int n, CFixBitVec &h )
{
	coordsToIndex<CFixBitVec,CFixBitVec>(p,m,n,h);
}

void indexToCoords( CFixBitVec *p, int m, int n, const CFixBitVec &h )
{
	indexToCoords<CFixBitVec,CFixBitVec>(p,m,n,h);
}

void coordsToCompactIndex( const CFixBitVec *p, const int *ms, int n,
	CFixBitVec &hc, int M, int m )
{
	coordsToCompactIndex<CFixBitVec,CFixBitVec>(p,ms,n,hc,M,m);
}

void compactIndexToCoords( CFixBitVec *p, const int *ms, int n,
	const CFixBitVec &hc, int M, int m )
{
	compactIndexToCoords<CFixBitVec,CFixBitVec>(p,ms,n,hc,M,m);
}

// fix -> big

void coordsToIndex( const CFixBitVec *p, int m, int n, CBigBitVec &h )
{
	coordsToIndex<CFixBitVec,CBigBitVec>(p,m,n,h);
}

void indexToCoords( CFixBitVec *p, int m, int n, const CBigBitVec &h )
{
	indexToCoords<CFixBitVec,CBigBitVec>(p,m,n,h);
}

void coordsToCompactIndex( const CFixBitVec *p, const int *ms, int n,
	CBigBitVec &hc, int M, int m )
{
	coordsToCompactIndex<CFixBitVec,CBigBitVec>(p,ms,n,hc,M,m);
}

void compactIndexToCoords( CFixBitVec *p, const int *ms, int n,
	const CBigBitVec &hc, int M, int m )
{
	compactIndexToCoords<CFixBitVec,CBigBitVec>(p,ms,n,hc,M,m);
}

// big -> big

void coordsToIndex( const CBigBitVec *p, int m, int n, CBigBitVec &h )
{
	coordsToIndex<CBigBitVec,CBigBitVec>(p,m,n,h);
}

void indexToCoords( CBigBitVec *p, int m, int n, const CBigBitVec &h )
{
	indexToCoords<CBigBitVec,CBigBitVec>(p,m,n,h);
}

void coordsToCompactIndex( const CBigBitVec *p, const int *ms, int n,
	CBigBitVec &hc, int M, int m )
{
	coordsToCompactIndex<CBigBitVec,CBigBitVec>(p,ms,n,hc,M,m);
}

void compactIndexToCoords( CBigBitVec *p, const int *ms, int n,
	const CBigBitVec &hc, int M, int m )
{
	compactIndexToCoords<CBigBitVec,CBigBitVec>(p,ms,n,hc,M,m);
}

}
