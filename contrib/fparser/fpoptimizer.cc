/***************************************************************************\
|* Function Parser for C++ v4.5.1                                          *|
|*-------------------------------------------------------------------------*|
|* Function optimizer                                                      *|
|*-------------------------------------------------------------------------*|
|* Copyright: Joel Yliluoma                                                *|
|*                                                                         *|
|* This library is distributed under the terms of the                      *|
|* GNU Lesser General Public License version 3.                            *|
|* (See lgpl.txt and gpl.txt for the license text.)                        *|
\***************************************************************************/

/* NOTE:
 This file contains generated code (from the optimizer sources) and is
 not intended to be modified by hand. If you want to modify the optimizer,
 download the development version of the library.
*/

#include "fpconfig.hh"
#ifdef FP_SUPPORT_OPTIMIZER
#include "fparser.hh"
#include "extrasrc/fptypes.hh"
#include "extrasrc/fpaux.hh"
#define lD4 )return
#define lC4 c0 tmp2
#define lB4 (lS i71
#define lA4 xZ a)
#define l94 eN cAbs);
#define l84 eN cMul
#define l74 if(n6 iJ3
#define l64 );if(
#define l54 tree.xJ1
#define l44 :if(&*e41){
#define l34 :{n21 r
#define l24 :if lP
#define l14 "Found "
#define l04 stackpos
#define iZ3 ,1,538 x0
#define iY3 "dup(%u) "
#define iX3 eX{assert
#define iW3 "%d, cost "
#define iV3 "PUSH " c72
#define iU3 "immed "<<
#define iT3 mFuncParsers
#define iS3 stderr
#define iR3 sep2=" "
#define iQ3 FPHASH_CONST
#define iP3 cache_needed[
#define iO3 fprintf
#define iN3 ::cout<<"Applying "
#define iM3 ||tree.GetOpcode
#define iL3 FUNCTIONPARSER_INSTANTIATE_OPTIMIZE
#define iK3 FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE
#define iJ3 HANDLE_UNARY_CONST_FUNC
#define iI3 ;tG1 eN
#define iH3 ;l71 yV1
#define iG3 tmp.nJ 0)
#define iF3 ;tmp eN
#define iE3 ,tC info
#define iD3 ;typedef
#define iC3 l2 2,2,
#define iB3 if(iX1==
#define iA3 const eM
#define i93 .n8 synth.
#define i83 within,
#define i73 nC==cLog2&&
#define i63 nC==cS3
#define i53 c_count
#define i43 s_count
#define i33 MaxOp
#define i23 2)lT 2*
#define i13 lQ3 val
#define i03 max.val
#define tZ3 iD);}if
#define tY3 for(xJ3
#define tX3 for xW1
#define tW3 sim.x11
#define tV3 b.Value)
#define tU3 b.Opcode
#define tT3 ].swap(
#define tS3 =synth.
#define tR3 codes[b
#define tQ3 whydump
#define tP3 c2 lU 2,
#define tO3 ;++a){
#define tN3 if(a>0){
#define tM3 );tmp tB1);
#define tL3 e41=r.specs;if(r.found){
#define tK3 info.
#define tJ3 ;GetParam
#define tI3 ;sim.x3 2,
#define tH3 .empty()
#define tG3 e41,info
#define tF3 =(*xB)[a]
#define tE3 for(;a<
#define tD3 nparams
#define tC3 l3 16,1,
#define tB3 l3 0,1,
#define tA3 0x12 nM
#define t93 nW 0,
#define t83 cAbs,nW
#define t73 iI;case
#define t63 ;tmp2.nJ
#define t53 fp_pow(
#define t43 false;}
#define t33 cAbsIf)
#define t23 ParamHolder
#define t13 FP_GetOpcodeName(
#define t03 .second
#define eZ3 ]t03
#define eY3 size_t>
#define eX3 l21++b)
#define eW3 =true;yZ1
#define eV3 =false;
#define eU3 ].first
#define eT3 Ne_Mask
#define eS3 Gt_Mask
#define eR3 Lt_Mask
#define eQ3 opcode,
#define eP3 resize(
#define eO3 public:
#define eN3 {data->
#define eM3 Immeds
#define eL3 switch(
#define eK3 pclone
#define eJ3 l8 a))
#define eI3 ,cEqual
#define eH3 cOr,l6
#define eG3 param.
#define eF3 AddFrom(
#define eE3 (p1.xJ1
#define eD3 middle2
#define eC3 ::string
#define eB3 );cZ1 eQ
#define eA3 };enum
#define e93 )const
#define e83 (e93{
#define e73 :start_at()
#define e63 NumConstant:
#define e53 eL3 nJ3.first cZ3
#define e43 (eG3
#define e33 ;}void
#define e23 xQ)e33
#define e13 std xV3<bool>
#define e03 ;}static t31
#define cZ3 ){case
#define cY3 newpow
#define cX3 (count
#define cW3 133,2,
#define cV3 Needs
#define cU3 byteCode
#define cT3 lW1 nC==
#define cS3 cPow&&tH
#define cR3 iC2 n01
#define cQ3 long iA1
#define cP3 factor_t
#define cO3 value1
#define cN3 Finite
#define cM3 a)l64!
#define cL3 fp_max(
#define cK3 else{if(
#define cJ3 xH());nE
#define cI3 max.known
#define cH3 known&&
#define cG3 eL3 lX1
#define cF3 xZ 0)nC
#define cE3 p2 nQ3
#define cD3 IsLogicalValue(xZ
#define cC3 eO p2;p2
#define cB3 cAbsNot
#define cA3 }switch
#define c93 stackptr
#define c83 cLog);xJ
#define c73 opcodes
#define c63 did_muli
#define c53 &Value){
#define c43 yL const
#define c33 used[b]
#define c23 sizeof(
#define c13 cAbsIf,
#define c03 cNotNot,
#define yZ3 cTan,y82
#define yY3 450998,
#define yX3 cExp2,nW
#define yW3 7168,
#define yV3 &param=*
#define yU3 xC3 eE2
#define yT3 default:
#define yS3 ){pow.cN2
#define yR3 range lO3
#define yQ3 range<yL2
#define yP3 range xH
#define yO3 Ge0Lt1
#define yN3 Gt0Le1
#define yM3 cAdd lT2
#define yL3 if(op1==
#define yK3 iterator
#define yJ3 begin();
#define yI3 TreeSet
#define yH3 parent
#define yG3 insert(i
#define yF3 comp.tP
#define yE3 change=
#define yD3 newrel
#define yC3 =*(cR*lH2
#define yB3 y73 eL3
#define yA3 ;if(half
#define y93 iZ2 size()
#define y83 ::ByteCodeSynth xH
#define y73 break;}
#define y63 void set
#define y53 b_needed
#define y43 cachepos
#define y33 grammar
#define y23 half=
#define y13 ,i7,1,l61+1);
#define y03 131,4,1,
#define xZ3 131,8,1,
#define xY3 4,1,2,1,
#define xX3 ){lM1 xQ=
#define xW3 ),child);
#define xV3 ::vector
#define xU3 FindPos(
#define xT3 src_pos
#define xS3 )i6|x91)
#define xR3 reserve(
#define xQ3 treeptr
#define xP3 yR lM1(
#define xO3 tS1 void
#define xN3 ImmedTag
#define xM3 a,const
#define xL3 RefCount
#define xK3 Birth();
#define xJ3 typename
#define xI3 mulgroup
#define xH3 unsigned
#define xG3 exponent
#define xF3 cost_t
#define xE3 fpdata
#define xD3 middle
#define xC3 result
#define xB3 cLog2by);
#define xA3 for e0 nU
#define x93 sqrt_cost
#define x83 const int
#define x73 mul_count
#define x63 maxValue1
#define x53 minValue1
#define x43 lQ3 known
#define x33 maxValue0
#define x23 minValue0
#define x13 ValueType
#define x03 xK n4 0),
#define nZ3 cB known
#define nY3 cB n4 0),
#define nX3 PlusInf
#define nW3 Value_t
#define nV3 (m xK val
#define nU3 >nW3(
#define nT3 abs_mul
#define nS3 pos_set
#define nR3 Rehash nV
#define nQ3 eQ ifp2);
#define nP3 sim.x3 1,
#define nO3 default_function_handling
#define nN3 yO 2,cAdd
#define nM3 subtree
#define nL3 invtree
#define nK3 MakeHash(
#define nJ3 parampair
#define nI3 rulenumit
#define nH3 l7 0,2,
#define nG3 l7 0,1,
#define nF3 (cond yW
#define nE3 goto eT
#define nD3 ){sim.Eat(
#define nC3 cJ1.SubTrees
#define nB3 cJ1.Others
#define nA3 )iI m;}case
#define n93 cAnd,l6
#define n83 cEqual,lB
#define n73 cNeg,lU 1
#define n63 MakeEqual
#define n53 n91,l4::
#define n43 n91,{l4::
#define n33 newbase
#define n23 fp_equal(
#define n13 branch1op
#define n03 branch2op
#define lZ3 overlap
#define lY3 truth_b
#define lX3 truth_a
#define lW3 found_dup
#define lV3 .tX1 n]
#define lU3 1)iM2(1));
#define lT3 nW3(2)));
#define lS3 ;TriTruthValue
#define lR3 lP3 eV3 if(
#define lQ3 xC3 xK
#define lP3 xC3 e61
#define lO3 xH xC3
#define lN3 tmp c0 tree
#define lM3 (xY l8 a)xL
#define lL3 if(n23
#define lK3 (iJ2 l8 1)x51
#define lJ3 (xI3
#define lI3 continue;
#define lH3 rangeutil
#define lG3 )val=iL1
#define lF3 );eB cU1);
#define lE3 synth lO2
#define lD3 Plan_Has(
#define lC3 StackMax)
#define lB3 iI true;}
#define lA3 xC3 xS
#define l93 namespace
#define l83 DelParam(
#define l73 ByteCode[
#define l63 inverted
#define l53 IsNever:
#define l43 depcodes
#define l33 explicit
#define l23 cCosh,nW
#define l13 VarBegin
#define l03 Params[a].
#define iZ2 Params.
#define iY2 ].data);
#define iX2 &&p xK val
#define iW2 )lT 3*3*
#define iV2 )iI lI
#define iU2 ,l1 xL1 3
#define iT2 ContainsOtherCandidates
#define iS2 tree.n81
#define iR2 tB1)eO
#define iQ2 begin(),
#define iP2 cond_add
#define iO2 cond_mul
#define iN2 cond_and
#define iM2 ,nW3
#define iL2 bool eI1
#define iK2 ,tree
#define iJ2 leaf1
#define iI2 Optimize()
#define iH2 costree
#define iG2 sintree
#define iF2 leaf_count
#define iE2 =GetParam(
#define iD2 sub_params
#define iC2 xC3.
#define iB2 printf(
#define iA2 cbrt_count
#define i92 sqrt_count
#define i82 l81);nE lD
#define i72 p1 eQ ifp1
#define i62 pcall_tree
#define i52 after_powi
#define i42 ;size_t
#define i32 xC3))n72
#define i22 i12 xC3(
#define i12 ){nW3
#define i02 nT i12 tmp=
#define tZ2 (lS));nE lD
#define tY2 tB1 false)
#define tX2 GetHash().
#define tW2 info=info;
#define tV2 ;nW3
#define tU2 yR t43
#define tT2 cEqual e81
#define tS2 cLog,nW
#define tR2 0x12},{{3,
#define tQ2 )lD4
#define tP2 );x5 l4::
#define tO2 },{{1,
#define tN2 ),0},{
#define tM2 .first.
#define tL2 tree c62
#define tK2 eN cond nC
#define tJ2 (eU==cN3&&
#define tI2 (half&63)-1;
#define tH2 tree eN
#define tG2 tree nC
#define tF2 MakeNEqual
#define tE2 Dump(std::
#define tD2 isInteger(
#define tC2 lM1 r;r eN
#define tB2 Comparison
#define tA2 needs_flip
#define t92 value]
#define t82 ~size_t(0)
#define t72 xD1 xU+1);
#define t62 const std::eV
#define t52 const char*
#define t42 Rule&rule,
#define t32 ,const e2&
#define t22 const lO1&
#define t12 >::res,b8<
#define t02 lI3}
#define eZ2 (size_t
#define eY2 for eZ2
#define eX2 eQ tree);
#define eW2 mul_item
#define eV2 innersub
#define eU2 cbrt_cost
#define eT2 best_cost
#define eS2 condition
#define eR2 nominator
#define eQ2 per_item
#define eP2 item_type
#define eO2 first2
#define eN2 l3 18,1,
#define eM2 ,cGreater
#define eL2 tU 1},0,
#define eK2 l81)));x5
#define eJ2 .what nM1
#define eI2 nC==cNot||
#define eH2 ,(long double)
#define eG2 if(tG2==
#define eF2 xK known)
#define eE2 .min.val
#define eD2 Decision
#define eC2 not_tree
#define eB2 group_by
#define eA2 .Become(xZ
#define e92 (xG3
#define e82 val,xG3);
#define e72 val=t53
#define e62 xG3=
#define e52 ->second
#define e42 cI1 lM1&
#define e32 targetpos
#define e22 ParamSpec
#define e12 rhs.hash2;}
#define e02 rhs.hash1
#define cZ2 struct
#define cY2 Forget()
#define cX2 &&cond eC))
#define cW2 source_tree
#define cV2 .IsImmed(
#define cU2 nB OPCODE
#define cT2 1)cV2
#define cS2 );t4=!t4;}
#define cR2 <tW,xF3>
#define cQ2 p1_evenness
#define cP2 isNegative(
#define cO2 xK cH3(
#define cN2 CopyOnWrite
#define cM2 ,std::cout)
#define cL2 neg_set
#define cK2 );synth.yS
#define cJ2 cNop,cNop}}
#define cI2 cTanh,cNop,
#define cH2 NewHash
#define cG2 >cZ2 cE<
#define cF2 matches
#define cE2 cI1 void*)&
#define cD2 l5 0,1,
#define cC2 cTan,nW
#define cB2 cCos,nW
#define cA2 tU 2},0,xL1
#define c92 (std::move(
#define c82 }else{xB=new
#define c72 ;DumpTree(
#define c62 .SetParam(
#define c52 );y73
#define c42 IsImmed()||
#define c32 +=1 iI n61;
#define c22 negated
#define c12 Specializer
#define c02 params
#define yZ2 coshtree
#define yY2 sinhtree
#define yX2 best_score
#define yW2 mulvalue
#define yV2 pow_item
#define yU2 subgroup
#define yT2 =lW1;bool iQ
#define yS2 PowiResult
#define yR2 0));yP3
#define yQ2 maxValue
#define yP2 minValue
#define yO2 fp_min(y6,
#define yN2 div_tree
#define yM2 pow_tree
#define yL2 nW3 nZ
#define yK2 preserve
#define yJ2 nW3(0.5)
#define yI2 PullResult()
#define yH2 dup_or_fetch
#define yG2 test_order
#define yF2 nJ3,
#define yE2 .param_count
#define yD2 shift(index)
#define yC2 rulenumber
#define yB2 ;i7.Remember(
#define yA2 ;}data;data.
#define y92 (tree))tL1
#define y82 l3 2,1,
#define y72 cLog,y82
#define y62 cTanh,nW
#define y52 cSinh,nW
#define y42 cInv,lU 1,
#define y32 constraints=
#define y22 GetDepth()
#define y12 factor_immed
#define y02 changes
#define xZ2 eQ cond l8
#define xY2 ;tree c0
#define xX2 l8 0));
#define xW2 c0 mul);
#define xV2 exp_diff
#define xU2 ExponentInfo
#define xT2 lower_bound(
#define xS2 factor
#define xR2 is_logical
#define xQ2 newrel_and
#define xP2 tP[c eE
#define xO2 res_stackpos
#define xN2 half_pos
#define xM2 ifdata.ofs
#define xL2 >>1)):(
#define xK2 CodeTreeData
#define xJ2 var_trees
#define xI2 parent_opcode
#define xH2 log2_exponent
#define xG2 dup_fetch_pos
#define xF2 TopLevel)
#define xE2 *)&*start_at;
#define xD2 {e2 start_at;
#define xC2 IsNever cV lD
#define xB2 cSin,nW
#define xA2 Value_EvenInt
#define x92 )){data xD
#define x82 MakeFalse,{l4
#define x72 AddCollection
#define x62 ConditionType
#define x52 SpecialOpcode
#define x42 =i e52.
#define x32 .IsDefined(
#define x22 .match_tree
#define x12 xK cH3 p
#define x02 i4 push_back(
#define nZ2 )lB3
#define nY2 xZ1 c4++a){
#define nX2 assimilated
#define nW2 denominator
#define nV2 fraction
#define nU2 ,eW,synth);
#define nT2 DUP_BOTH();
#define nS2 template lL
#define nR2 0x80000000u
#define nQ2 template lY
#define nP2 if(synth.Find(
#define nO2 IsDescendantOf
#define nN2 nM2.erase(cs_it);
#define nM2 TreeCounts
#define nL2 1)?(poly^(
#define nK2 iM2(1))){
#define nJ2 IsImmed()){if(
#define nI2 iM2(-1)))xO
#define nH2 bool t4 eV3
#define nG2 SetOpcode(
#define nF2 found_log2
#define nE2 div_params
#define nD2 xK val)
#define nC2 immed_sum
#define nB2 :sim.Eat(1,
#define nA2 l73++IP]
#define n92 OPCODE(opcode)
#define n82 ;sim.Push(
#define n72 break;xC3*=
#define n62 FactorStack xH
#define n52 IsAlways cV lD
#define n42 282870 x0
#define n32 c03 nW
#define n22 cNot,nW
#define n12 replacing_slot
#define n02 RefParams
#define lZ2 if_always[
#define lY2 WhatDoWhenCase
#define lX2 exponent_immed
#define lW2 new_base_immed
#define lV2 base_immed
#define lU2 {case IsAlways:
#define lT2 ||op1==
#define lS2 nX DelParams()
#define lR2 data[a eZ3
#define lQ2 .AddItem(atree
#define lP2 if(newrel_or==
#define lO2 .AddOperation(
#define lN2 DUP_ONE(apos);
#define lM2 flipped
#define lL2 .UseGetNeeded(
#define lK2 e8 2,131,
#define lJ2 (IfData&ifdata
#define lI2 minimum_need
#define lH2 )nJ3 t03
#define lG2 [xU-1-offset].
#define lF2 l73 a
#define lE2 yC Immed.size());
#define lD2 1 yC i4 size()
#define lC2 GetOpcode())
#define lB2 {if(GetOpcode()
#define lA2 *start_at){xB=(
#define l92 const lM1
#define l82 OptimizedUsing
#define l72 Var_or_Funcno
#define l62 l72;
#define l52 GetParams(
#define l42 crc32_t
#define l32 signed_chain
#define l22 MinusInf
#define l12 return true;
#define l02 n_immeds
#define iZ1 stack.size()
#define iY1 FindClone(xQ
#define iX1 l73 IP]
#define iW1 needs_rehash
#define iV1 AnyWhere_Rec
#define iU1 ~xH3(0)
#define iT1 41,42,43,44,
#define iS1 p1_logical_b
#define iR1 p0_logical_b
#define iQ1 p1_logical_a
#define iP1 p0_logical_a
#define iO1 2*2*2)lT 3*
#define iN1 ,cGreaterOrEq
#define iM1 eN tG2);
#define iL1 func(val);nK1
#define iK1 *const func)
#define iJ1 TreeCountItem
#define iI1 ;if(op==
#define iH1 else if(
#define iG1 synth.DoDup(
#define iF1 cache_needed
#define iE1 e8 2,1,e8 2,
#define iD1 treelist
#define iC1 has_bad_balance
#define iB1 cP3 xS2
#define iA1 double)xG3
#define i91 fp_abs(i03))
#define i81 fp_abs(min.val)
#define i71 ,xZ l81));nE lD
#define i61 cNEqual
#define i51 Oneness_NotOne|
#define i41 Value_IsInteger
#define i31 Constness_Const
#define i21 DumpHashesFrom(
#define i11 l82(
#define i01 ));xB1 eQ yA l8
#define tZ1 reltype
#define tY1 SequenceOpcodes
#define tX1 sep_list[
#define tW1 eB yU2
#define tV1 goto fail;}
#define tU1 l1 0x4 nM
#define tT1 template<
#define tS1 n02);
#define tR1 TreeCountType xH
#define tQ1 cB set(fp_floor)
#define tP1 >(nW3(1),
#define tO1 nW3(0.0)){y0
#define tN1 );tree.swap(tmp);
#define tM1 lA4 cV2))
#define tL1 goto redo;
#define tK1 n92);
#define tJ1 <<tree.tX2
#define tI1 ,l5 2,1,
#define tH1 std::cout<<"POP "
#define tG1 divgroup
#define tF1 );eB xI3);
#define tE1 stack[iZ1-
#define tD1 <xH3>&ByteCode,size_t&IP,size_t limit,size_t y4
#define tC1 stack.push_back(
#define tB1 .Rehash(
#define tA1 ParsePowiMuli
#define t91 synth.PushImmed(
#define t81 MaxChildDepth
#define t71 std::pair<It,It>
#define t61 m xK known
#define t51 ,cPow,
#define t41 lJ 2},0,0x4 tO2
#define t31 inline xH3
#define t21 lJ 2},0,xL1 1,
#define t11 }inline
#define t01 lJ 1},0,
#define eZ1 Sign_Negative
#define eY1 Value_Logical
#define eX1 new_factor_immed
#define eW1 grammar_rules[*r]
#define eV1 nQ2 void
#define eU1 synth.xI 1
#define eT1 ,i7 nU2
#define eS1 lB 0x4 tO2
#define eR1 t51 lB
#define eQ1 eY2 a=
#define eP1 occurance_pos
#define eO1 exponent_hash
#define eN1 exponent_list
#define eM1 CollectionSet xH
#define eL1 CollectMulGroup(
#define eK1 source_set
#define eJ1 xG3,yI3
#define eI1 operator
#define eH1 FindAndDup(tree);
#define eG1 xZ cT2)&&
#define eF1 back().thenbranch
#define eE1 ParamSpec_Extract
#define eD1 retry_anyparams_3
#define eC1 retry_anyparams_2
#define eB1 e7(),std xV3<
#define eA1 needlist_cached_t
#define e91 CodeTreeImmed xH(
#define e81 ,l0 2,
#define e71 long value
#define e61 .min.known
#define e51 public e7,public std xV3<
#define e41 (*xB)[a].start_at
#define e31 eY2 b=0;b<
#define e21 by_float_exponent
#define e11 ;xG3 tB1
#define e01 fp_equal e92
#define cZ1 new_exp
#define cY1 end()&&i->first==
#define cX1 return BecomeZero;
#define cW1 return BecomeOne;
#define cV1 if(lR.size()<=n2)
#define cU1 addgroup
#define cT1 found_log2by
#define cS1 nC==cB3)
#define cR1 if(keep_powi
#define cQ1 l72)
#define cP1 branch1_backup
#define cO1 l21 a-->0;)if(
#define cN1 branch2_backup
#define cM1 exponent_map
#define cL1 plain_set
#define cK1 .l83 a);
#define cJ1 NeedList
#define cI1 (const
#define cH1 cI1 nW3&
#define cG1 LightWeight(
#define cF1 );xQ.SetParamsMove(
#define cE1 ,PowiCache&i7,
#define cD1 if(value
#define cC1 nQ2 yZ
#define cB1 nQ2 static
#define cA1 cB val=nW3(1);
#define c91 .GetParamCount()==
#define c81 should_regenerate=true;
#define c71 should_regenerate,
#define c61 Collection
#define c51 RelationshipResult
#define c41 Subdivide_Combine(
#define c31 e93 yR
#define c21 rhs c31 hash1
#define c11 eY2 a y7
#define c01 best_sep_factor
#define yZ1 iH1!xC3
#define yY1 needlist_cached
#define yX1 eQ3 bool pad
#define yW1 changed=true;
#define yV1 l83 a);}
#define yU1 MakesInteger(
#define yT1 const nW3&value
#define yS1 best_sep_cost
#define yR1 MultiplicationRange
#define yQ1 pihalf_limits
#define yP1 yO 2,cMul);lD
#define yO1 ;synth.StackTopIs(
#define yN1 yO 2,cPow);lD
#define yM1 n_stacked
#define yL1 cH2.hash1
#define yK1 AnyParams_Rec
#define yJ1 79,122,123,160,161,
#define yI1 );tree.l83
#define yH1 Become(value l8 0))
#define yG1 PositionalParams,0}
#define yF1 always_sincostan
#define yE1 Recheck_RefCount_Div
#define yD1 Recheck_RefCount_Mul
#define yC1 xI3.
#define yB1 xI3;xI3 eN
#define yA1 MultiplyAndMakeLong(
#define y91 cMul);iG3);tmp
#define y81 ByteCodeSynth xH&synth)
#define y71 nW3(0)
#define y61 ;range.multiply(
#define y51 ;m xK template set_if<
#define y41 nQ2 bool
#define y31 ;nQ2
#define y21 }y73 case
#define y11 eZ2 a=0;a<
#define y01 l93 FPoptimizer_Optimize
#define xZ1 ;for y11
#define xY1 Rehash()xY2 p1)
#define xX1 (p0 e61&&p0 eE2>=nW3(0.0))
#define xW1 eZ2 a=c4 a-->0;)
#define xV1 if(synth.FindAndDup(
#define xU1 SynthesizeParam(
#define xT1 grammar_func
#define xS1 tree cV2)cV
#define xR1 252421 x0 24830,
#define xQ1 c1 529654 x0
#define xP1 l2 0,2,165888 x0
#define xO1 cCos,y82
#define xN1 cIf,lB 0x4},{{3,
#define xM1 l1 tA3
#define xL1 0x0},{{
#define xK1 Modulo_Radians},
#define xJ1 GetImmed()
#define xI1 PositionType
#define xH1 CollectionResult
#define xG1 const_offset
#define xF1 inline TriTruthValue
#define xE1 stacktop_desired
#define xD1 SetStackTop(
#define xC1 FPoptimizer_ByteCode
#define xB1 changed_if
#define xA1 l64 covers_plus1
#define x91 (xH3
#define x81 ;for x91 a=0;a<y1;++a)
#define x71 ;std::cout<<
#define x61 y71)
#define x51 xL leaf2 l8
#define x41 cond_type
#define x31 fphash_value_t
#define x21 Recheck_RefCount_RDiv
#define x11 SwapLastTwoInStack();
#define x01 fPExponentIsTooLarge(
#define nZ1 CollectMulGroup_Item(
#define nY1 pair<nW3,yI3>
#define nX1 Rehash()xY2 r);}
#define nW1 nN xD1 xU-1);
#define nV1 covers_full_cycle
#define nU1 lA4.xJ1;
#define nT1 AssembleSequence(
#define nS1 inverse_nominator
#define nR1 252180 x0 281854,
#define nQ1 {DataP slot_holder(y3[
#define nP1 <<std::dec<<")";}
#define nO1 },{l4::MakeNotP0,l4::
#define nN1 :return p eE2
#define nM1 !=xM)if(TestCase(
#define nL1 &&IsLogicalValue(
#define nK1 else*this=model;}
#define nJ1 std::pair<T1,T2>&
#define nI1 tT1 xJ3
#define nH1 has_good_balance_found
#define nG1 n_occurrences
#define nF1 found_log2_on_exponent
#define nE1 covers_minus1
#define nD1 needs_resynth
#define nC1 immed_product
#define nB1 yB3 bitmask&
#define nA1 Sign_Positive
#define n91 ::MakeTrue
#define n81 SetParamMove(
#define n71 CodeTreeImmed(nW3(
#define n61 Suboptimal
#define n51 n_as_tanh_param
#define n41 opposite=
#define n31 x31(
#define n21 MatchResultType
#define n11 needs_sincos
#define n01 resulting_exponent
#define lZ1 val):Value(Value::
#define lY1 Unknown:yT3;}
#define lX1 GetLogicalValue(xZ
#define lW1 GetParam(a)
#define lV1 },{l4::MakeNotNotP1,l4::
#define lU1 },{l4::MakeNotNotP0,l4::
#define lT1 cSin,y82
#define lS1 IsImmed()i12
#define lR1 AddFunctionOpcode(
#define lQ1 ;cH2.hash2+=
#define lP1 xH(rule.repl_param_list,
#define lO1 CodeTree
#define lN1 ,cIf,l0 3,
#define lM1 lO1 xH
#define lL1 void FunctionParserBase
#define lK1 SetParams(l52));
#define lJ1 o<<"("<<std::hex<<data.
#define lI1 IfBalanceGood(
#define lH1 n_as_tan_param
#define lG1 changed_exponent
#define lF1 iX2<nW3(
#define lE1 inverse_denominator
#define lD1 retry_positionalparams_2
#define lC1 xH3 index
#define lB1 situation_flags&
#define lA1 518 x0 400412,
#define l91 );n81 0,xG3);l83 1);
#define l81 1).xJ1
#define l71 cN2();
#define l61 recursioncount
#define l51 PlanNtimesCache(
#define l41 >){int mStackPtr=0;
#define l31 FPoptimizer_Grammar
#define l21 GetParamCount();
#define l11 GetPositivityInfo(tree)!=
#define l01 ParamSpec_SubFunctionData
#define iZ )){lM1
#define iY ,2,1);synth.xV if(found[data.
#define iX AddOperation(cInv,1,1);synth.xV}
#define iW ]);lE3
#define iV t82){synth.yS
#define iU PositionalParams_Rec
#define iT .data.subfunc_opcode
#define iS ,{l4::MakeNotP1,l4::
#define iR l73 xM2+
#define iQ needs_cow=GetRefCount()>1;
#define iP DumpTreeWithIndent(*this);
#define iO CalculateResultBoundaries(
#define iN tT1 xH3 Compare>
#define iM edited_powgroup
#define iL has_unknown_max
#define iK has_unknown_min
#define iJ (tree.GetParamCount()
#define iI ;return
#define iH set(fp_ceil nA3
#define iG {eL3 type cZ3 cond_or:
#define iF y11 c4++a)if(remaining[a])
#define iE static const yP3
#define iD tree.l83 a
#define iC synthed_tree
#define iB 7168 x0 401798,
#define iA SelectedParams,0},0,xL1
#define i9 by_exponent
#define i8 collections
#define i7 cache
#define i6 ;x02 nR2
#define i5 )i6);
#define i4 ByteCode.
#define i3 goto ReplaceTreeWithOne;case
#define i2 c0 comp.cL1[a].value);
#define i1 !=xM lD4 lZ2
#define i0 tB1);tH2 iJ2 nC);tree.
#define tZ e21.data
#define tY l33 xK2(
#define tX needs_sinhcosh
#define tW int_exponent_t
#define tV )x71 std::endl;DumpHashes(
#define tU AnyParams,
#define tT tN2 nW3(
#define tS nQ2 nD
#define tR MakeFalse,l4::
#define tQ :goto ReplaceTreeWithZero;case
#define tP relationships
#define tO c72 tree)x71"\n";
#define tN ,ByteCode,IP,limit,y4,stack);
#define tM 408964 x0 24963,
#define tL tU 0}},{ReplaceParams,
#define tK matched_params
#define tJ [n2 eU3=true;lR[n2 eZ3
#define tI l31::Grammar*
#define tH powgroup l8
#define tG }},{ProduceNewTree,2,1,
#define tF n71(
#define tE has_mulgroups_remaining
#define tD const l01
#define tC MatchInfo xH&
#define tB Rehash();iD2.push_back(
#define tA best_factor
#define t9 RootPowerTable xH::RootPowers[
#define t8 MatchPositionSpec_AnyParams xH
#define t7 l93 FPoptimizer_CodeTree
#define t6 n_as_sinh_param
#define t5 n_as_cosh_param
#define t4 is_signed
#define t3 result_positivity
#define t2 xK known eV3
#define t1 biggest_minimum
#define t0 cond_tree
#define eZ else_tree
#define eY then_tree
#define eX lM1&tree)
#define eW sequencing
#define eV string t13
#define eU valueType
#define eT ReplaceTreeWithParam0;
#define eS =iO xZ
#define eR nO nW3(-lU3
#define eQ .AddParam(
#define eP :{AdoptChildrenWithSameOpcode(tree);
#define eO ;lM1
#define eN .nG2
#define eM std xV3<lM1>
#define eL if_stack
#define eK (l52));yC1 Rehash();
#define eJ lD4 IsNever iI Unknown;}
#define eI n_as_sin_param
#define eH n_as_cos_param
#define eG PowiResolver::
#define eF cIf,tB3
#define eE ].relationship
#define eD PACKED_GRAMMAR_ATTRIBUTE;
#define eC .BalanceGood
#define eB AddParamMove(
#define eA back().endif_location
#define e9 x31 key
#define e8 130,1,
#define e7 MatchPositionSpecBase
#define e6 l33 lO1(
#define e5 smallest_maximum
#define e4 ]!=t82&&found[data.
#define e3 factor_needs_rehashing
#define e2 MatchPositionSpecBaseP
#define e1 xJ3 tR1::yK3
#define e0 eZ2 a=l21 a
#define cZ eE1 xH(nQ.param_list,
#define cY ,cM 124024 x0 139399,
#define cX ,cM 142456 x0 141449,
#define cW ,cM 528504 x0 24713,
#define cV lD4 false;
#define cU 163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,243,244,245,246,249,250,251,253,255,256,257,258,259}};}
#define cT ];};extern"C"{
#define cS 27,28,29,30,31,32,33,35,36,
#define cR const ParamSpec_SubFunction
#define cQ const ParamSpec_ParamHolder
#define cP ComparisonSetBase::
#define cO otherhalf
#define cN StackState
#define cM l2 16,2,
#define cL const SequenceOpCode xH
#define cK paramholder_matches[x6]
#define cJ MatchPositionSpec_PositionalParams xH
#define cI l92&tree,std::ostream&o
#define cH paramholder_matches.
#define cG nW3(1.5)*fp_const_pi xH()
#define cF CalculatePowiFactorCost(
#define cE ImmedHashGenerator
#define cD ::map<fphash_t,std::set<std eC3> >
#define cC T1,xJ3 T2>inline iL2()(
#define cB m.min.
#define cA has_nonlogical_values
#define c9 from_logical_context)
#define c8 tU 0}},{ProduceNewTree,
#define c7 eQ1 xY.l21 a-->0;)
#define c6 POWI_CACHE_SIZE
#define c5 ,l7 2,1,
#define c4 tree.l21
#define c3 );void lR1 xH3 eQ3 c12<
#define c2 ,cAdd,
#define c1 ,l2 18,2,
#define c0 .eB
#define yZ static inline lM1
#define yY ++IP;t02 iB3 c73.
#define yX },{l4::xM,l4::Never},{l4::xM,l4::Never}}
#define yW .FoundChild
#define yV BalanceResultType
#define yU {lM1 tmp iF3
#define yT lE3 GetOpcode(),
#define yS DoDup(found[data.
#define yR {return
#define yQ const yR data->
#define yP +=fp_const_twopi xH();
#define yO );sim.Eat(
#define yN (p0 xK cH3 p0 xK val<=fp_const_negativezero xH())
#define yM for y11 l21++a){if(
#define yL static void nK3 nB fphash_t&cH2,
#define yK MatchPositionSpec_AnyWhere
#define yJ if e43 data.match_type==
#define yI ):xL3(0),Opcode(
#define yH void OutFloatHex(std::ostream&o,
#define yG {lM1 tmp,tmp2;tmp2 eN
#define yF AddParam(CodeTreeImmed(
#define yE ,xJ3 lM1::
#define yD AssembleSequence_Subdivide(
#define yC ]=nR2|xH3(
#define yB ();pow eN cLog);tH2 cMul);
#define yA branch2
#define y9 xH3 c;xH3 short l[
#define y8 {nM2.erase(i);t02
#define y7 =0;a<yH3.l21++a)if(
#define y6 fp_sin(min),fp_sin(max))
#define y5 fp_const_twopi xH()l64
#define y4 factor_stack_base
#define y3 data->Params
#define y2 );t02 if(list tM2 xJ1==nW3(
#define y1 nQ yE2
#define y0 tree.ReplaceWithImmed(
#define xZ tree l8
#define xY branch1
#define xX (nI3 r=range.first;r!=range t03;++r){
#define xW =fp_cosh(cB val);m xK val=fp_cosh(m nD2;
#define xV StackTopIs(*this)iI;}
#define xU StackTop
#define xT FPOPT_autoptr
#define xS +=xC3 iI xC3;}nQ2 inline nW3
#define xR int_exponent
#define xQ newnode
#define xP has_highlevel_opcodes
#define xO {if(needs_cow){l71 goto
#define xN ,t01 xL1 1,
#define xM Unchanged
#define xL .IsIdenticalTo(
#define xK .max.
#define xJ sim.AddConst(
#define xI GetStackTop()-
#define xH <nW3>
#define xG cI=std::cout
#define xF best_selected_sep
#define xE cAnd,tL
#define xD ->Recalculate_Hash_NoRecursion();}
#define xC l21++a)if(ApplyGrammar(y33,lA4,
#define xB position
#define xA for y11 c4++a){yP3
#define x9 std xV3<lO1>
#define x8 SetParam(0,iftree l8 0))eO p1;p1 eN
#define x7 TestImmedConstraints e43 constraints iK2)cV
#define x6 paramholder_index
#define x5 l12 case
#define x4 occurance_counts
#define x3 SwapLastTwoInStack(yO
#define x2 )){tree.FixIncompleteHashes();}
#define x1 );cR1 nD3 1,cInv c52 xJ-1 yN1
#define x0 ,{2,
#define nZ >p eS a)l64 p.
#define nY t63 0))iF3 cInv);tmp lC4)iI
#define nX ));xB1 tB1);tH2 op1);tree.
#define nW l0 1,
#define nV ()xY2 p2);tH2 iftree nC);tL1}
#define nU -->0;){l92&powgroup=lW1;if(powgroup
#define nT if(xZ 0)cV2)
#define nS eR1 0x4 nM
#define nR const FPoptimizer_CodeTree::lM1&tree
#define nQ model_tree
#define nP ),rangehalf xH model=rangehalf xH()){if(known
#define nO return yP3(
#define nN ){using l93 FUNCTIONPARSERTYPES;
#define nM },{{2,
#define nL eL2 xL1
#define nK eM&n02
#define nJ AddParam(xZ
#define nI ConstantFolding_LogicCommon(tree,cP
#define nH nI1 Ref>inline void xT<Ref>::
#define nG cOr,tL 16,1,
#define nF ):data(new xK2 xH(
#define nE goto do_return;}
#define nD xK2 xH::xK2(
#define nC .GetOpcode()
#define nB FUNCTIONPARSERTYPES::
#define nA b;}};tT1>cZ2 Comp<nB
#define n9 l72(),Params(),Hash(),Depth(1),i11 0){}
#define n8 SynthesizeByteCode(synth);
#define n7 while(ApplyGrammar(cE2
#define n6 GetIntegerInfo(xZ 0))==IsAlways)nE3
#define n5 xY2 xB1 nZ2
#define n4 template set_if<cGreater>(nW3(
#define n3 DumpParams xH e43 data.param_list,eG3 data yE2,o);
#define n2 restholder_index
#define n1 lM1 xG3;xG3 l84);xG3 eQ
#define n0 (lS l64 fp_nequal(tmp,x61){y0 nW3(1)/tmp);nE}lD
#define lZ :if(ParamComparer xH()(Params[1],Params[0])){std::swap(Params[0],Params[1]);Opcode=
#define lY <xJ3 nW3>
#define lX xH tmp iF3 cPow);iG3);tmp.yF nW3(
#define lW i31,0x0},
#define lV eB pow l8 1));pow.l83 1);pow tB1);iS2 0,pow);goto NowWeAreMulGroup;}
#define lU GroupFunction,0},lW{{
#define lT iM2(1)/nW3(
#define lS xZ 0).xJ1
#define lR restholder_matches
#define lQ yL1|=key;x31 crc=(key>>10)|(key<<(64-10))lQ1((~n31 crc))*3)^1234567;}};
#define lP (xZ 0)cV2)&&xZ cT2)){y0
#define lO xB1;xB1 iM1 xB1 c0 xZ 0));xB1 eQ xY l8
#define lN nQ2 lM1::lO1(
#define lM tL2 0,xZ 0)xX2 tL2 1,CodeTreeImmed(
#define lL lY void ByteCodeSynth xH::lR1 xH3 eQ3 c12<
#define lK cMul,lU 2,
#define lJ cMul,tU
#define lI iO tmp);}case
#define lH :yE3 comp.AddRelationship(atree l8 0),atree l8 1),cP
#define lG cPow,l0 2
#define lF xJ3 nW3>inline iL2()cH1 xM3 nW3&b)yR a
#define lE {yP3 m eS 0));
#define lD break;case
#define lC eV1 lM1::
#define lB yG1,0,
#define lA l1 0x0 nM
#define l9 ?0:1))eO xB1;xB1 iM1 xB1.SetParamsMove(tree.l52));xB1 tB1);tH2
#define l8 .GetParam(
#define l7 cAdd,tL
#define l6 SelectedParams,0},0,0x0 nM
#define l5 lJ 0}},{ReplaceParams,
#define l4 RangeComparisonData
#define l3 yG1},{ProduceNewTree,
#define l2 yG1},{ReplaceParams,
#define l1 cMul,SelectedParams,0},0,
#define l0 lB xL1
#ifdef _MSC_VER
typedef
xH3
int
l42;
#else
#include <stdint.h>
typedef
uint_least32_t
l42;
#endif
l93
crc32{enum{startvalue=0xFFFFFFFFUL,poly=0xEDB88320UL}
;tT1
l42
crc>cZ2
b8{enum{b1=(crc&nL2
crc
xL2
crc>>1),b2=(b1&nL2
b1
xL2
b1>>1),b3=(b2&nL2
b2
xL2
b2>>1),b4=(b3&nL2
b3
xL2
b3>>1),b5=(b4&nL2
b4
xL2
b4>>1),b6=(b5&nL2
b5
xL2
b5>>1),b7=(b6&nL2
b6
xL2
b6>>1),res=(b7&nL2
b7
xL2
b7>>1)}
;}
;inline
l42
update(l42
crc,xH3
b){
#define B4(n) b8<n t12 n+1 t12 n+2 t12 n+3>::res
#define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
static
const
l42
table[256]={R(0x00),R(0x10),R(0x20),R(0x30),R(0x40),R(0x50),R(0x60),R(0x70),R(0x80),R(0x90),R(0xA0),R(0xB0),R(0xC0),R(0xD0),R(0xE0),R(0xF0)}
;
#undef R
#undef B4
return((crc>>8))^table[(crc^b)&0xFF];t11
l42
calc_upd(l42
c,const
xH3
char*buf,size_t
size){l42
value=c;eY2
p=0;p<size;++p)value=update(value,buf[p])iI
value;t11
l42
calc
cI1
xH3
char*buf,size_t
size)yR
calc_upd(startvalue,buf,size);}
}
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH
nI1
Ref>class
xT{eO3
xT():p(0){}
xT(Ref*b):p(b){xK3}
xT
cI1
xT&b):p(b.p){xK3
t11
Ref&eI1*(c31*p;t11
Ref*eI1->(c31
p;}
xT&eI1=(Ref*b){Set(b)iI*this;}
xT&eI1=cI1
xT&b){Set(b.p)iI*this;}
#ifdef FP_SUPPORT_CXX11_MOVE
xT(xT&&b):p(b.p){b.p=0;}
xT&eI1=(xT&&b){if(p!=b.p){cY2;p=b.p;b.p=0;}
return*this;}
#endif
~xT(){cY2
e33
UnsafeSetP(Ref*newp){p=newp
e33
swap(xT<Ref>&b){Ref*tmp=p;p=b.p;b.p=tmp;}
private:inline
static
void
Have(Ref*p2);inline
void
cY2;inline
void
xK3
inline
void
Set(Ref*p2);private:Ref*p;}
;nH
cY2{if(!p
lD4;p->xL3-=1;if(!p->xL3)delete
p;}
nH
Have(Ref*p2){if(p2)++(p2->xL3);}
nH
Birth(){Have(p);}
nH
Set(Ref*p2){Have(p2);cY2;p=p2;}
#endif
#include <utility>
cZ2
Compare2ndRev{nI1
T>inline
iL2()cI1
T&xM3
T&b
c31
a
t03>b
t03;}
}
;cZ2
Compare1st{nI1
cC
const
nJ1
xM3
nJ1
b
c31
a.first<b.first;}
nI1
cC
const
nJ1
a,T1
b
c31
a.first<b;}
nI1
cC
T1
xM3
nJ1
b
c31
a<b.first;}
}
;
#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH
#ifdef _MSC_VER
typedef
xH3
long
long
x31;
#define FPHASH_CONST(x) x##ULL
#else
#include <stdint.h>
typedef
uint_fast64_t
x31;
#define FPHASH_CONST(x) x##ULL
#endif
l93
FUNCTIONPARSERTYPES{cZ2
fphash_t{x31
hash1,hash2;fphash_t():hash1(0),hash2(0){}
fphash_t
cI1
x31&xM3
x31&b):hash1(a),hash2(b){}
iL2==cI1
fphash_t&c21==e02&&hash2==e12
iL2!=cI1
fphash_t&c21!=e02||hash2!=e12
iL2<cI1
fphash_t&c21!=e02?hash1<e02:hash2<e12}
;}
#endif
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
l93
l31{cZ2
Grammar;}
l93
xC1{nQ2
class
ByteCodeSynth;}
t7{nQ2
class
lO1
y31
cZ2
xK2
y31
class
lO1{typedef
xT<xK2
xH>DataP;DataP
data;eO3
lO1();~lO1();cZ2
OpcodeTag{}
;e6
cU2
o,OpcodeTag);cZ2
FuncOpcodeTag{}
;e6
cU2
o,xH3
f,FuncOpcodeTag);cZ2
xN3{}
;e6
const
nW3&v,xN3);
#ifdef FP_SUPPORT_CXX11_MOVE
e6
nW3&&v,xN3);
#endif
cZ2
VarTag{}
;e6
xH3
varno,VarTag);cZ2
CloneTag{}
;e6
t22
b,CloneTag);void
GenerateFrom
cI1
xJ3
FunctionParserBase
xH::Data&data,bool
keep_powi=false);void
GenerateFrom
cI1
xJ3
FunctionParserBase
xH::Data&data,const
x9&xJ2,bool
keep_powi=false);void
SynthesizeByteCode(std
xV3<xH3>&cU3,std
xV3
xH&immed,size_t&stacktop_max);void
SynthesizeByteCode(xC1
y83&synth,bool
MustPopTemps=true
e93
i42
SynthCommonSubExpressions(xC1::y81
const;void
SetParams
cI1
x9&xO3
SetParamsMove(x9&tS1
lO1
GetUniqueRef();
#ifdef FP_SUPPORT_CXX11_MOVE
void
SetParams(x9&&tS1
#endif
void
SetParam
eZ2
which,t22
b);void
n81
size_t
which,lO1&b);void
AddParam
cI1
lO1&param);void
eB
lO1&param);void
AddParams
cI1
x9&xO3
AddParamsMove(x9&xO3
AddParamsMove(x9&n02,size_t
n12);void
DelParam
eZ2
index);void
DelParams();void
Become
cI1
lO1&b);inline
size_t
GetParamCount(c31
l52).size();t11
lO1&GetParam
eZ2
n)yR
l52)[n];t11
t22
GetParam
eZ2
n
c31
l52)[n];t11
void
nG2
cU2
o)eN3
Opcode=o;t11
cU2
GetOpcode()yQ
Opcode;t11
nB
fphash_t
GetHash()yQ
Hash;t11
const
x9&l52
c31
y3;t11
x9&l52)yR
y3;t11
size_t
y22
yQ
Depth;t11
const
nW3&xJ1
yQ
Value;t11
xH3
GetVar()yQ
l62
t11
xH3
GetFuncNo()yQ
l62
t11
bool
IsDefined(c31
GetOpcode()!=nB
cNop;t11
bool
IsImmed(c31
GetOpcode()==nB
cImmed;t11
bool
IsVar(c31
GetOpcode()==nB
l13;t11
xH3
GetRefCount()yQ
xL3
e33
ReplaceWithImmed
cH1
i);void
Rehash(bool
constantfolding=true);void
Sort();inline
void
Mark_Incompletely_Hashed()eN3
Depth=0;t11
bool
Is_Incompletely_Hashed()yQ
Depth==0;t11
const
tI
GetOptimizedUsing()yQ
l82;t11
void
SetOptimizedUsing
cI1
tI
g)eN3
l82=g;}
bool
RecreateInversionsAndNegations(bool
prefer_base2=false);void
FixIncompleteHashes();void
swap(lO1&b){data.swap(b.data);}
bool
IsIdenticalTo
cI1
lO1&b
e93;void
l71}
y31
cZ2
xK2{int
xL3;cU2
Opcode
tV2
Value;xH3
l62
eM
Params;nB
fphash_t
Hash
i42
Depth;const
tI
l82;xK2();xK2
cI1
xK2&b);tY
cU2
o);tY
cU2
o,xH3
f);tY
const
nW3&i);
#ifdef FP_SUPPORT_CXX11_MOVE
tY
nW3&&i);xK2(xK2&&b);
#endif
bool
IsIdenticalTo
cI1
xK2&b
e93;void
Sort();void
Recalculate_Hash_NoRecursion();private:void
eI1=cI1
xK2&b);}
y31
yZ
CodeTreeImmed
cH1
i)xP3
i
yE
xN3());}
#ifdef FP_SUPPORT_CXX11_MOVE
cC1
CodeTreeImmed(nW3&&i)yR
lM1
c92
i)yE
xN3());}
#endif
cC1
CodeTreeOp(cU2
opcode)xP3
opcode
yE
OpcodeTag());}
cC1
CodeTreeFuncOp(cU2
eQ3
xH3
f)xP3
eQ3
f
yE
FuncOpcodeTag());}
cC1
CodeTreeVar
x91
varno)xP3
varno
yE
VarTag());}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
eV1
DumpHashes(xG)y31
void
DumpTree(xG)y31
void
DumpTreeWithIndent(xG,const
std
eC3&indent="\\"
);
#endif
}
#endif
#endif
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH
#include <iostream>
t7{nQ2
class
lO1;}
l93
l31{enum
ImmedConstraint_Value{ValueMask=0x07,Value_AnyNum=0x0,xA2=0x1,Value_OddInt=0x2,i41=0x3,Value_NonInteger=0x4,eY1=0x5
eA3
ImmedConstraint_Sign{SignMask=0x18,Sign_AnySign=0x00,nA1=0x08,eZ1=0x10,Sign_NoIdea=0x18
eA3
ImmedConstraint_Oneness{OnenessMask=0x60,Oneness_Any=0x00,Oneness_One=0x20,Oneness_NotOne=0x40
eA3
ImmedConstraint_Constness{ConstnessMask=0x180,Constness_Any=0x00,i31=0x80,Constness_NotConst=0x100
eA3
Modulo_Mode{Modulo_None=0,Modulo_Radians=1
eA3
Situation_Flags{LogicalContextOnly=0x01,NotForIntegers=0x02,OnlyForIntegers=0x04,OnlyForComplex=0x08,NotForComplex=0x10
eA3
x52{NumConstant,t23,SubFunction
eA3
ParamMatchingType{PositionalParams,SelectedParams,tU
GroupFunction
eA3
RuleType{ProduceNewTree,ReplaceParams}
;
#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif
typedef
std::pair<x52,const
void*>e22
y31
e22
eE1
x91
paramlist,lC1)y31
bool
ParamSpec_Compare
cI1
void*xM3
void*b,x52
type);xH3
ParamSpec_GetDepCode
cI1
e22&b);cZ2
ParamSpec_ParamHolder{lC1:8;xH3
constraints:9;xH3
depcode:15;}
eD
nQ2
cZ2
ParamSpec_NumConstant{nW3
constvalue;xH3
modulo;}
;cZ2
l01{xH3
param_count:2;xH3
param_list:30;cU2
subfunc_opcode:8;ParamMatchingType
match_type:3;xH3
n2:5;}
eD
cZ2
ParamSpec_SubFunction{l01
data;xH3
constraints:9;xH3
depcode:7;}
eD
cZ2
Rule{RuleType
ruletype:2;xH3
situation_flags:5;xH3
repl_param_count:2+9;xH3
repl_param_list:30;l01
match_tree;}
eD
cZ2
Grammar{xH3
rule_count;xH3
short
rule_list[999
cT
extern
const
Rule
grammar_rules[];}
eV1
DumpParam
cI1
e22&p,std::ostream&o=std::cout)y31
void
DumpParams
x91
paramlist,xH3
count,std::ostream&o=std::cout);}
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define CONSTANT_POS_INF HUGE_VAL
#define CONSTANT_NEG_INF (-HUGE_VAL)
l93
FUNCTIONPARSERTYPES{nQ2
inline
nW3
fp_const_pihalf()yR
fp_const_pi
xH()*yJ2;}
nQ2
inline
nW3
fp_const_twopi(i22
fp_const_pi
xH());lA3
fp_const_twoe(i22
fp_const_e
xH());lA3
fp_const_twoeinv(i22
fp_const_einv
xH());lA3
fp_const_negativezero()yR-Epsilon
xH::value;}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#include <iostream>
y01{using
l93
l31;using
t7;using
l93
FUNCTIONPARSERTYPES
y31
class
MatchInfo{eO3
std
xV3<std::pair<bool,eM> >lR;eM
paramholder_matches;std
xV3<xH3>tK;eO3
MatchInfo():lR(),paramholder_matches(),tK(){}
eO3
bool
SaveOrTestRestHolder
x91
n2,iA3&iD1){cV1{lR.eP3
n2+1);lR
tJ=iD1
lB3
if(lR[n2
eU3==false){lR
tJ=iD1
lB3
iA3&found=lR[n2
eZ3;if(iD1.size()!=found.size()cV
for
y11
iD1.size();++a)if(!iD1[a]xL
found[a])cV
l12}
void
SaveRestHolder
x91
n2,eM&iD1){cV1
lR.eP3
n2+1);lR
tJ.swap(iD1);}
bool
SaveOrTestParamHolder
x91
x6,l92&xQ3){if(cH
size()<=x6){cH
xR3
x6+1);cH
eP3
x6);cH
push_back(xQ3
nZ2
if(!cK
x32)){cK=xQ3
lB3
return
xQ3
xL
cK)e33
SaveMatchedParamIndex(lC1){tK.push_back(index);}
l92&GetParamHolderValueIfFound
x91
x6
e93{static
l92
dummytree;if(cH
size()<=x6
lD4
dummytree
iI
cK;}
l92&GetParamHolderValue
x91
x6
c31
cK;}
bool
HasRestHolder
x91
n2
c31
lR.size()>n2&&lR[n2
eU3==true;}
iA3&GetRestHolderValues
x91
n2
e93{static
iA3
empty_result;cV1
return
empty_result
iI
lR[n2
eZ3;}
const
std
xV3<xH3>&GetMatchedParamIndexes(c31
tK
e33
swap(tC
b){lR.swap(b.lR);cH
swap(b.paramholder_matches);tK.swap(b.tK);}
tC
eI1=cI1
tC
b){lR=b.lR;paramholder_matches=b.paramholder_matches;tK=b.tK
iI*this;}
}
;class
e7
iD3
xT<e7>e2;class
e7{eO3
int
xL3;eO3
e7():xL3(0){}
virtual~e7(){}
}
;cZ2
n21{bool
found;e2
specs;n21(bool
f):found(f),specs(){}
n21(bool
f
t32
s):found(f),specs(s){}
}
y31
void
SynthesizeRule
cI1
t42
lM1&tree
iE3)y31
n21
TestParam
cI1
e22&yF2
l92&tree
t32
start_at
iE3)y31
n21
TestParams(tD&nQ,l92&tree
t32
start_at
iE3,bool
xF2
y31
bool
ApplyGrammar
cI1
Grammar&y33,FPoptimizer_CodeTree::lM1&tree,bool
from_logical_context=false)y31
void
ApplyGrammars(FPoptimizer_CodeTree::eX
y31
bool
IsLogisticallyPlausibleParamsMatch(tD&c02,const
eX;}
l93
l31{eV1
DumpMatch
cI1
t42
nR,const
FPoptimizer_Optimize::tC
info,bool
DidMatch,std::ostream&o=std::cout)y31
void
DumpMatch
cI1
t42
nR,const
FPoptimizer_Optimize::tC
info,t52
tQ3,std::ostream&o=std::cout);}
#endif
#include <string>
t62
l31::x52
yX1=false);t62
cU2
yX1=false);
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
using
l93
l31;using
l93
FUNCTIONPARSERTYPES;t62
l31::x52
yX1){
#if 1
t52
p=0;eL3
opcode
cZ3
e63
p="NumConstant"
;lD
t23:p="ParamHolder"
;lD
SubFunction:p="SubFunction"
;y73
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().size()<12)tmp<<' '
iI
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().size()<5)tmp<<' '
iI
tmp.str();
#endif
}
t62
cU2
yX1){
#if 1
t52
p=0;eL3
opcode
cZ3
cAbs:p="cAbs"
;lD
cAcos:p="cAcos"
;lD
cAcosh:p="cAcosh"
;lD
cArg:p="cArg"
;lD
cAsin:p="cAsin"
;lD
cAsinh:p="cAsinh"
;lD
cAtan:p="cAtan"
;lD
cAtan2:p="cAtan2"
;lD
cAtanh:p="cAtanh"
;lD
cCbrt:p="cCbrt"
;lD
cCeil:p="cCeil"
;lD
cConj:p="cConj"
;lD
cCos:p="cCos"
;lD
cCosh:p="cCosh"
;lD
cCot:p="cCot"
;lD
cCsc:p="cCsc"
;lD
cExp:p="cExp"
;lD
cExp2:p="cExp2"
;lD
cFloor:p="cFloor"
;lD
cHypot:p="cHypot"
;lD
cIf:p="cIf"
;lD
cImag:p="cImag"
;lD
cInt:p="cInt"
;lD
cLog:p="cLog"
;lD
cLog2:p="cLog2"
;lD
cLog10:p="cLog10"
;lD
cMax:p="cMax"
;lD
cMin:p="cMin"
;lD
cPolar:p="cPolar"
;lD
cPow:p="cPow"
;lD
cReal:p="cReal"
;lD
cSec:p="cSec"
;lD
cSin:p="cSin"
;lD
cSinh:p="cSinh"
;lD
cSqrt:p="cSqrt"
;lD
cTan:p="cTan"
;lD
cTanh:p="cTanh"
;lD
cTrunc:p="cTrunc"
;lD
cImmed:p="cImmed"
;lD
cJump:p="cJump"
;lD
cNeg:p="cNeg"
;lD
cAdd:p="cAdd"
;lD
cSub:p="cSub"
;lD
cMul:p="cMul"
;lD
cDiv:p="cDiv"
;lD
cMod:p="cMod"
;lD
cEqual:p="cEqual"
;lD
i61:p="cNEqual"
;lD
cLess:p="cLess"
;lD
cLessOrEq:p="cLessOrEq"
;lD
cGreater:p="cGreater"
;lD
cGreaterOrEq:p="cGreaterOrEq"
;lD
cNot:p="cNot"
;lD
cAnd:p="cAnd"
;lD
cOr:p="cOr"
;lD
cDeg:p="cDeg"
;lD
cRad:p="cRad"
;lD
cFCall:p="cFCall"
;lD
cPCall:p="cPCall"
;break;
#ifdef FP_SUPPORT_OPTIMIZER
case
cFetch:p="cFetch"
;lD
cPopNMov:p="cPopNMov"
;lD
cLog2by:p="cLog2by"
;lD
cNop:p="cNop"
;break;
#endif
case
cSinCos:p="cSinCos"
;lD
cSinhCosh:p="cSinhCosh"
;lD
cB3:p="cAbsNot"
;lD
cAbsNotNot:p="cAbsNotNot"
;lD
cAbsAnd:p="cAbsAnd"
;lD
cAbsOr:p="cAbsOr"
;lD
cAbsIf:p="cAbsIf"
;lD
cDup:p="cDup"
;lD
cInv:p="cInv"
;lD
cSqr:p="cSqr"
;lD
cRDiv:p="cRDiv"
;lD
cRSub:p="cRSub"
;lD
cNotNot:p="cNotNot"
;lD
cRSqrt:p="cRSqrt"
;lD
l13:p="VarBegin"
;y73
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().size()<12)tmp<<' '
iI
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().size()<5)tmp<<' '
iI
tmp.str();
#endif
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#ifndef FP_GENERATING_POWI_TABLE
enum{MAX_POWI_BYTECODE_LENGTH=20}
;
#else
enum{MAX_POWI_BYTECODE_LENGTH=999}
;
#endif
enum{MAX_MULI_BYTECODE_LENGTH=3}
;l93
xC1{nQ2
class
ByteCodeSynth{eO3
ByteCodeSynth():ByteCode(),Immed(),cN(),xU(0),StackMax(0){i4
xR3
64);Immed.xR3
8);cN.xR3
16)e33
Pull(std
xV3<xH3>&bc,std
xV3
xH&imm,size_t&StackTop_max){for
x91
a=0;a<i4
size()tO3
lF2]&=~nR2;}
i4
swap(bc);Immed.swap(imm);StackTop_max=StackMax;}
size_t
GetByteCodeSize(c31
i4
size();}
size_t
GetStackTop(c31
xU
e33
PushVar
x91
varno){x02
varno);t72}
void
PushImmed(nW3
immed
nN
x02
cImmed);Immed.push_back(immed);t72}
void
StackTopIs(nR,int
offset=0){if((int)xU>offset){cN
lG2
first=true;cN
lG2
second=tree;}
}
bool
IsStackTop(nR,int
offset=0
c31(int)xU>offset&&cN
lG2
first&&cN
lG2
second
xL
tree);t11
void
EatNParams
x91
eat_count){xU-=eat_count
e33
ProducedNParams
x91
produce_count){xD1
xU+produce_count)e33
DoPopNMov
eZ2
e32,size_t
srcpos
nN
x02
cPopNMov
xS3
e32
xS3
srcpos);xD1
srcpos+1);cN[e32]=cN[srcpos];xD1
e32+1)e33
DoDup
eZ2
xT3
nN
if(xT3==xU-1){x02
cDup);}
else{x02
cFetch
xS3
xT3);}
t72
cN[xU-1]=cN[xT3];}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tT1
int>void
Dump(){std::ostream&o=std::cout;o<<"Stack state now("
<<xU<<"):\n"
xZ1
xU
tO3
o<<a<<": "
;if(cN[a
eU3){nR=cN[a
eZ3;o<<'['<<std::hex<<(void*)(&tree.l52))<<std::dec<<','<<tree.GetRefCount()<<']'
c72
tree,o);}
else
o<<"?"
;o<<"\n"
;}
o<<std::flush;}
#endif
size_t
xU3
nR
e93{eQ1
xU;a-->0;)if(cN[a
eU3&&cN[a
eZ3
xL
tree
tQ2
a
iI
t82;}
bool
Find(nR
c31
xU3
tree)!=t82;}
bool
FindAndDup(nR){size_t
pos=xU3
tree
l64
pos!=t82){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<l14"duplicate at ["
<<pos<<"]: "
c72
tree)x71" -- issuing cDup or cFetch\n"
;
#endif
DoDup(pos
nZ2
return
t43
cZ2
IfData{size_t
ofs;}
;void
SynthIfStep1
lJ2,cU2
op
nW1
xM2=i4
size();x02
op
i5
x02
nR2)e33
SynthIfStep2
lJ2
nW1
iR
lD2+2);iR
2
lE2
xM2=i4
size();x02
cJump
i5
x02
nR2)e33
SynthIfStep3
lJ2
nW1
i4
back()|=nR2;iR
lD2-1);iR
2
lE2
xD1
xU+1)xZ1
xM2
tO3
if(lF2]==cJump&&lF2+1]==(nR2|(xM2-1))){lF2+lD2-1);lF2+2
lE2
cA3(lF2]cZ3
cAbsIf:case
cIf:case
cJump:case
cPopNMov:a+=2;lD
cFCall:case
cPCall:case
cFetch:a+=1;break;yT3
y73}
}
protected:void
xD1
size_t
value){xU=value;if(xU>lC3{StackMax=xU;cN.eP3
lC3;}
}
protected:std
xV3<xH3>ByteCode;std
xV3
xH
Immed;std
xV3<std::pair<bool,FPoptimizer_CodeTree::lM1> >cN
i42
xU
i42
StackMax;private:void
incStackPtr(){if(xU+2>lC3
cN.eP3
StackMax=xU+2);}
tT1
bool
IsIntType,bool
IsComplexType>cZ2
c12{}
;eO3
void
AddOperation
x91
eQ3
xH3
eat_count,xH3
produce_count=1){EatNParams(eat_count);lR1
opcode);ProducedNParams(produce_count)e33
lR1
xH3
eQ3
c12<false,false>c3
false,true>c3
true,false>c3
true,true>);inline
void
lR1
xH3
opcode){lR1
eQ3
c12<bool(nB
IsIntType
xH::xC3),bool(nB
IsComplexType
xH::xC3)>());}
}
y31
cZ2
SequenceOpCode
y31
cZ2
tY1{static
cL
AddSequence;static
cL
MulSequence;}
y31
void
nT1
long
count,cL&eW,y81;}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;l93
xC1{nQ2
cZ2
SequenceOpCode{nW3
basevalue;xH3
op_flip;xH3
op_normal,op_normal_flip;xH3
op_inverse,op_inverse_flip;}
y31
cL
tY1
xH::AddSequence={y71,cNeg
c2
cAdd,cSub,cRSub}
y31
cL
tY1
xH::MulSequence={nW3(1),cInv,cMul,cMul,cDiv,cRDiv}
;
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode ByteCode
#define mImmed Immed
nS2
false,false
l41
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 0
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif
# include "extrasrc/fp_opcode_add.inc"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nS2
true,false
l41
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 0
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif
# include "extrasrc/fp_opcode_add.inc"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nS2
false,true
l41
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 1
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif
# include "extrasrc/fp_opcode_add.inc"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nS2
true,true
l41
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 1
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc99-extensions"
#endif
# include "extrasrc/fp_opcode_add.inc"
#ifdef __clang__
#pragma clang diagnostic pop
#endif
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#endif
#undef findName
#undef mImmed
#undef mByteCode
#undef mData
#undef TryCompilePowi
}
using
l93
xC1;
#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
l93
xC1{
#ifndef FP_GENERATING_POWI_TABLE
extern
const
xH3
char
powi_table[POWI_TABLE_SIZE];const
#endif
xH3
char
powi_table[POWI_TABLE_SIZE]={0,1,1,1,2,1,2,1,xY3
4,1,2,xZ3
2,1,xY3
8,cW3
y03
15,1,16,1,2,1,4,1,2,xZ3
2,1,4,cW3
1,16,1,25,y03
27,5,8,3,2,1,30,1,31,3,32,1,2,1,xY3
8,1,2,y03
39,1,16,137,2,1,4,cW3
xZ3
45,135,4,31,2,5,32,1,2,131,50,1,51,1,8,3,2,1,54,1,55,3,16,1,57,133,4,137,2,135,60,1,61,3,62,133,63,1,iE1
131,iE1
139,lK2
e8
30,1,130,137,2,31,lK2
e8
e8
130,cW3
1,e8
e8
2,1,130,133,iE1
61,130,133,62,139,130,137,e8
lK2
e8
e8
iE1
131,e8
e8
130,131,2,133,lK2
130,141,e8
130,cW3
1,e8
5,135,e8
lK2
e8
lK2
130,133,130,141,130,131,e8
e8
2,131}
;}
static
x83
c6=256;
#define FPO(x)
l93{class
PowiCache{private:int
i7[c6];int
iF1[c6];eO3
PowiCache():i7(),iF1(){i7[1]=1;}
bool
Plan_Add(e71,int
count){cD1>=c6
cV
iF1[t92+=count
iI
i7[t92!=0
e33
lD3
e71){cD1<c6)i7[t92=1
e33
Start
eZ2
value1_pos){for(int
n=2;n<c6;++n)i7[n]=-1;Remember(1,value1_pos);DumpContents();}
int
Find(e71
e93{cD1<c6){if(i7[t92>=0){FPO(iO3(iS3,"* I found %ld from cache (%u,%d)\n",value,(unsigned)cache[value],iP3 value]))iI
i7[t92;}
}
return-1
e33
Remember(e71,size_t
l04){cD1>=c6
lD4;FPO(iO3(iS3,"* Remembering that %ld can be found at %u (%d uses remain)\n",value,(unsigned)l04,iP3 value]));i7[t92=(int)l04
e33
DumpContents
e83
FPO(for(int a=1;a<POWI_CACHE_SIZE;++a)if(cache[a]>=0||iP3 a]>0){iO3(iS3,"== cache: sp=%d, val=%d, needs=%d\n",cache[a],a,iP3 a]);})}
int
UseGetNeeded(e71){cD1>=0&&value<c6
lD4--iF1[t92
iI
0;}
}
y31
size_t
yD
long
count
cE1
cL&eW,y81
y31
void
c41
size_t
apos,long
aval,size_t
bpos,long
bval
cE1
xH3
cumulation_opcode,xH3
cimulation_opcode_flip,y81;void
l51
e71
cE1
int
need_count,int
l61=0){cD1<1
lD4;
#ifdef FP_GENERATING_POWI_TABLE
if(l61>32)throw
false;
#endif
if(i7.Plan_Add(value,need_count
tQ2;long
y23
1;cD1<POWI_TABLE_SIZE){y23
powi_table[t92
yA3&128){half&=127
yA3&64)y23-tI2
FPO(iO3(iS3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,value/half));l51
half
y13
i7.lD3
half)iI;}
iH1
half&64){y23-tI2}
}
else
cD1&1)y23
value&((1<<POWI_WINDOW_SIZE)-1);else
y23
value/2;long
cO=value-half
yA3>cO||half<0)std::swap(half,cO);FPO(iO3(iS3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,otherhalf))yA3==cO){l51
half,i7,2,l61+1);}
else{l51
half
y13
l51
cO>0?cO:-cO
y13}
i7.lD3
value);}
nQ2
size_t
yD
e71
cE1
cL&eW,y81{int
y43=i7.Find(value
l64
y43>=0)yR
y43;}
long
y23
1;cD1<POWI_TABLE_SIZE){y23
powi_table[t92
yA3&128){half&=127
yA3&64)y23-tI2
FPO(iO3(iS3,"* I want %ld, my plan is %ld * %ld\n",value,half,value/half))i42
xN2=yD
half
eT1
if(i7
lL2
half)>0||xN2!=eU1){iG1
xN2)yB2
half,eU1);}
nT1
value/half
nU2
size_t
l04=eU1
yB2
value,l04);i7.DumpContents()iI
l04;}
iH1
half&64){y23-tI2}
}
else
cD1&1)y23
value&((1<<POWI_WINDOW_SIZE)-1);else
y23
value/2;long
cO=value-half
yA3>cO||half<0)std::swap(half,cO);FPO(iO3(iS3,"* I want %ld, my plan is %ld + %ld\n",value,half,value-half))yA3==cO){size_t
xN2=yD
half
eT1
c41
xN2,half,xN2,half,i7,eW.op_normal,eW.op_normal_flip,synth);}
else{long
part1=half;long
part2=cO>0?cO:-cO
i42
part1_pos=yD
part1
eT1
size_t
part2_pos=yD
part2
eT1
FPO(iO3(iS3,"Subdivide(%ld: %ld, %ld)\n",value,half,otherhalf));c41
part1_pos,part1,part2_pos,part2,i7,cO>0?eW.op_normal:eW.op_inverse,cO>0?eW.op_normal_flip:eW.op_inverse_flip,synth);}
size_t
l04=eU1
yB2
value,l04);i7.DumpContents()iI
l04;}
eV1
c41
size_t
apos,long
aval,size_t
bpos,long
bval
cE1
xH3
cumulation_opcode,xH3
cumulation_opcode_flip,y81{int
a_needed=i7
lL2
aval);int
y53=i7
lL2
bval);bool
lM2
eV3
#define DUP_BOTH() do{if(apos<bpos){size_t tmp=apos;apos=bpos;bpos=tmp;lM2=!lM2;}FPO(iO3(iS3,"-> " iY3 iY3"op\n",(unsigned)apos,(unsigned)bpos));iG1 apos);iG1 apos==bpos?eU1:bpos);}while(0)
#define DUP_ONE(p) do{FPO(iO3(iS3,"-> " iY3"op\n",(unsigned)p));iG1 p);}while(0)
if(a_needed>0){if(y53>0){nT2}
cK3
bpos!=eU1)nT2
else{lN2
lM2=!lM2;}
}
}
iH1
y53>0){if(apos!=eU1)nT2
else
DUP_ONE(bpos);}
cK3
apos==bpos&&apos==eU1)lN2
iH1
apos==eU1&&bpos==synth.xI
2){FPO(iO3(iS3,"-> op\n"));lM2=!lM2;}
iH1
apos==synth.xI
2&&bpos==eU1)FPO(iO3(iS3,"-> op\n"));iH1
apos==eU1)DUP_ONE(bpos);iH1
bpos==eU1){lN2
lM2=!lM2;}
else
nT2}
lE3
lM2?cumulation_opcode_flip:cumulation_opcode,2);}
eV1
cG1
long
count,cL&eW,y81{while
cX3<256){int
y23
xC1::powi_table[count]yA3&128){half&=127;cG1
half
nU2
count/=half;}
else
y73
if
cX3==1
lD4;if(!cX3&1)){lE3
cSqr,1);cG1
count/2
nU2}
else{iG1
eU1);cG1
count-1
nU2
lE3
cMul,2);}
}
}
l93
xC1{eV1
nT1
long
count,cL&eW,y81{if
cX3==0)t91
eW.basevalue);else{bool
tA2
eV3
if
cX3<0){tA2=true;count=-count;}
if(false)cG1
count
nU2
iH1
count>1){PowiCache
i7;l51
count,i7,1)i42
xE1
tS3
GetStackTop();i7.Start(eU1);FPO(iO3(iS3,"Calculating result for %ld...\n",count))i42
xO2=yD
count
eT1
size_t
n_excess
tS3
xI
xE1;if(n_excess>0||xO2!=xE1-1){synth.DoPopNMov(xE1-1,xO2);}
}
if(tA2)lE3
eW.op_flip,1);}
}
}
#endif
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH
t7{l93
lH3{iN
cZ2
Comp{}
;tT1>cZ2
Comp<nB
cLess>{tT1
lF<nA
cLessOrEq>{tT1
lF<=nA
cGreater>{tT1
lF>nA
cGreaterOrEq>{tT1
lF>=nA
cEqual>{tT1
lF==nA
i61>{tT1
lF!=b;}
}
;}
nQ2
cZ2
rangehalf{nW3
val;bool
known;rangehalf():val(),known(false){}
rangehalf
cH1
v):val(v),known(true){t11
y63
cH1
v){known=true;val=v;}
y63(nW3(iK1(nW3
nP
lG3
y63(nW3(iK1
cH1
nP
lG3
iN
void
set_if(nW3
v
iM2(iK1(nW3
nP&&lH3::Comp<Compare>()(val,v)lG3
iN
void
set_if
cH1
v
iM2(iK1
cH1
nP&&lH3::Comp<Compare>()(val,v)lG3}
y31
cZ2
range{rangehalf
xH
min,max;range():min(),max(){}
range(nW3
mi
iM2
ma):min(mi),max(ma){}
range(bool
iM2
ma):min(),max(ma){}
range(nW3
mi,bool):min(mi),max(){}
void
set_abs();void
set_neg();}
y31
bool
IsLogicalTrueValue
cI1
yP3&p,bool
abs)y31
bool
IsLogicalFalseValue
cI1
yP3&p,bool
abs);}
#endif
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH
t7{enum
TriTruthValue{IsAlways,IsNever,Unknown}
y31
yP3
iO
const
eX
y31
bool
IsLogicalValue
cI1
eX
y31
TriTruthValue
GetIntegerInfo
cI1
eX
y31
xF1
GetEvennessInfo
cI1
eX{if(!tree
cV2
tQ2
Unknown;yT1=l54;if(nB
isEvenInteger(value
tQ2
IsAlways;if(nB
isOddInteger(value)eJ
nQ2
xF1
GetPositivityInfo
cI1
eX{yP3
p=iO
tree
l64
p
e61&&p
eE2>=nW3(tQ2
IsAlways;if(p
xK
known
lF1)eJ
nQ2
xF1
GetLogicalValue
e42
tree,bool
abs){yP3
p=iO
tree
l64
IsLogicalTrueValue(p,abs
tQ2
IsAlways;if(IsLogicalFalseValue(p,abs)eJ}
#endif
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH
t7{eV1
ConstantFolding(eX;}
#endif
l93{using
l93
FUNCTIONPARSERTYPES;using
t7;cZ2
ComparisonSetBase{enum{eR3=0x1,Eq_Mask=0x2,Le_Mask=0x3,eS3=0x4,eT3=0x5,Ge_Mask=0x6}
;static
int
Swap_Mask(int
m)yR(m&Eq_Mask)|((m&eR3)?eS3:0)|((m&eS3)?eR3:0);}
enum
c51{Ok,BecomeZero,BecomeOne,n61
eA3
x62{cond_or,iN2,iO2,iP2}
;}
y31
cZ2
ComparisonSet:public
ComparisonSetBase{cZ2
tB2{lM1
a
eO
b;int
relationship;tB2():a(),b(),relationship(){}
}
;std
xV3<tB2>tP;cZ2
Item{lM1
value;bool
c22;Item():value(),c22(false){}
}
;std
xV3<Item>cL1;int
xG1;ComparisonSet():tP(),cL1(),xG1(0){}
c51
AddItem
e42
a,bool
c22,x62
type){eY2
c=0;c<cL1.size();++c)if(cL1[c].value
xL
a)){if(c22!=cL1[c].c22)iG
cW1
case
iP2:cL1.erase(cL1.begin()+c);xG1
c32
case
iN2:case
iO2:cX1}
}
return
n61;}
Item
pole;pole.value=a;pole.c22=c22;cL1.push_back(pole)iI
Ok;}
c51
AddRelationship(lM1
a,lM1
b,int
tZ1,x62
type)iG
if(tZ1==7)cW1
lD
iP2:if(tZ1==7){xG1
c32}
lD
iN2:case
iO2:if(tZ1==0)cX1
y73
if(!(a.GetHash()<b.GetHash())){a.swap(b);tZ1=Swap_Mask(tZ1);}
eY2
c=0;c<tP.size();++c){if(tP[c].a
xL
a)&&tP[c].b
xL
b))iG{int
yD3=xP2|tZ1;if(yD3==7)cW1
xP2=yD3;y73
case
iN2:case
iO2:{int
yD3=xP2&tZ1;if(yD3==0)cX1
xP2=yD3;y73
case
iP2:{int
newrel_or=xP2|tZ1;int
xQ2=xP2&tZ1;lP2
5&&xQ2==0){xP2=eT3
iI
n61;}
lP2
7&&xQ2==0){xG1+=1;tP.erase(tP.begin()+c)iI
n61;}
lP2
7&&xQ2==Eq_Mask){xP2=Eq_Mask;xG1
c32}
t02}
return
n61;}
}
tB2
comp;comp.a=a;comp.b=b;comp.relationship=tZ1;tP.push_back(comp)iI
Ok;}
}
;nI1
nW3,xJ3
CondType>bool
ConstantFolding_LogicCommon(lM1&tree,CondType
x41,bool
xR2){bool
should_regenerate
eV3
ComparisonSet
xH
comp
nY2
xJ3
cP
c51
yE3
cP
Ok;l92&atree=lA4;eL3
atree
nC
cZ3
cEqual
lH
Eq_Mask,x41);lD
i61
lH
eT3,x41);lD
cLess
lH
eR3,x41);lD
cLessOrEq
lH
Le_Mask,x41);lD
cGreater
lH
eS3,x41);lD
cGreaterOrEq
lH
Ge_Mask,x41);lD
cNot:yE3
comp
lQ2
l8
0),true,x41);lD
cNotNot:yE3
comp
lQ2
l8
0),false,x41);break;yT3
if(xR2||IsLogicalValue(atree))yE3
comp
lQ2,false,x41);cA3(change){ReplaceTreeWithZero:y0
0)iI
true;ReplaceTreeWithOne:y0
1);x5
cP
Ok:lD
cP
BecomeZero
tQ
cP
BecomeOne:i3
cP
n61:c81
y73}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_LogicCommon: "
tO
#endif
if(xR2){tree.DelParams();}
else{tX3{l92&atree=xZ
a
l64
IsLogicalValue(atree))iD);}
}
for
y11
comp.cL1.size()tO3
if(comp.cL1[a].c22){tC2
cNot);r
i2
r.nX1
iH1!xR2){tC2
cNotNot);r
i2
r.nX1
else
tree
i2}
for
y11
yF3.size()tO3
tC2
cNop);eL3
yF3[a
eE
cZ3
cP
eR3:r
eN
cLess);lD
cP
Eq_Mask:r
eN
cEqual);lD
cP
eS3:r
eN
cGreater);lD
cP
Le_Mask:r
eN
cLessOrEq);lD
cP
eT3:r
eN
i61);lD
cP
Ge_Mask:r
eN
cGreaterOrEq
c52
r
c0
yF3[a].a);r
c0
yF3[a].b);r.nX1
if(comp.xG1!=0)tree.yF
nW3(comp.xG1)));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_LogicCommon: "
tO
#endif
l12}
return
t43
y41
ConstantFolding_AndLogic(iX3(tree.GetOpcode()==cAnd
iM3()==cAbsAnd)iI
nI
iN2,true);}
y41
ConstantFolding_OrLogic(iX3(tree.GetOpcode()==cOr
iM3()==cAbsOr)iI
nI
cond_or,true);}
y41
ConstantFolding_AddLogicItems(iX3(tree.GetOpcode()==cAdd)iI
nI
iP2,false);}
y41
ConstantFolding_MulLogicItems(iX3(tree.GetOpcode()==cMul)iI
nI
iO2,false);}
}
#include <vector>
#include <map>
#include <algorithm>
l93{using
l93
FUNCTIONPARSERTYPES;using
t7;cZ2
CollectionSetBase{enum
xH1{Ok,n61}
;}
y31
cZ2
CollectionSet:public
CollectionSetBase{cZ2
c61{lM1
value
eO
xS2;bool
e3;c61():value(),xS2(),e3(false){}
c61
e42
v,l92&f):value(v),xS2(f),e3(false){}
}
;std::multimap<fphash_t,c61>i8
iD3
xJ3
std::multimap<fphash_t,c61>::yK3
xI1;CollectionSet():i8(){}
xI1
FindIdenticalValueTo
e42
value){fphash_t
hash=value.GetHash();for(xI1
i=i8.xT2
hash);i!=i8.cY1
hash;++i){cD1
xL
i
e52.value
tQ2
i;}
return
i8.end();}
bool
Found
cI1
xI1&b)yR
b!=i8.end();}
xH1
AddCollectionTo
e42
xS2,const
xI1&into_which){c61&c=into_which
e52;if(c.e3)c.xS2
eQ
xS2);else{lM1
add;add
eN
cAdd);add
c0
c.xS2);add
eQ
xS2);c.xS2.swap(add);c.e3=true;}
return
n61;}
xH1
x72
e42
value,l92&xS2){const
fphash_t
hash=value.GetHash();xI1
i=i8.xT2
hash);for(;i!=i8.cY1
hash;++i){if(i
e52.value
xL
value
tQ2
AddCollectionTo(xS2,i);}
i8.yG3,std::make_pair(hash,c61(value,xS2)))iI
Ok;}
xH1
x72
e42
a)yR
x72(a,n71
1)));}
}
y31
cZ2
ConstantExponentCollection{typedef
eM
yI3
iD3
std::nY1
xU2;std
xV3<xU2>data;ConstantExponentCollection():data(){}
void
MoveToSet_Unique
cH1
eJ1&eK1){data.push_back(std::nY1(eJ1()));data.back()t03.swap(eK1)e33
MoveToSet_NonUnique
cH1
eJ1&eK1){xJ3
std
xV3<xU2>::yK3
i=std::xT2
data.iQ2
data.end(),xG3,Compare1st()l64
i!=data.cY1
xG3){i
e52.yG3
e52.end(),eK1.iQ2
eK1.end());}
else{data.yG3,std::nY1
e92,eK1));}
}
bool
iI2{bool
changed
eV3
std::sort(data.iQ2
data.end(),Compare1st());redo:for
y11
data.size();++a
i12
exp_a=data[a
eU3;lL3
exp_a
iM2(1)))lI3
eY2
b=a+1;b<data.size();++b
i12
exp_b=data[b
eU3
tV2
xV2=exp_b-exp_a;if(xV2>=fp_abs(exp_a))break
tV2
exp_diff_still_probable_integer=xV2*nW3(16
l64
tD2
exp_diff_still_probable_integer)&&!(tD2
exp_b)&&!tD2
xV2))){yI3&a_set=lR2;yI3&b_set=data[b
eZ3;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantExponentCollection iteration:\n"
;tE2
cout);
#endif
if(isEvenInteger(exp_b)&&!isEvenInteger(xV2+exp_a
iZ
tmp2;tmp2
l84);tmp2.SetParamsMove(b_set);tmp2
iR2
tmp
iF3
cAbs);tmp
lC4
tM3
b_set.eP3
1);b_set[0
tT3
tmp);}
a_set.insert(a_set.end(),b_set.iQ2
b_set.end());yI3
b_copy=b_set;data.erase(data.begin()+b);MoveToSet_NonUnique(xV2,b_copy);yW1
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantExponentCollection iteration:\n"
;tE2
cout);
#endif
tL1}
}
}
return
changed;}
#ifdef DEBUG_SUBSTITUTIONS
void
tE2
ostream&out){for
y11
data.size()tO3
out.precision(12);out<<data[a
eU3<<": "
;e31
lR2.size();++b){if(b>0)out<<'*'
c72
lR2[b],out);}
out<<std::endl;}
}
#endif
}
y31
static
lM1
nZ1
lM1&value,bool&xP){eL3
value
nC
cZ3
cPow:{lM1
e62
value
l8
1);value.yH1
iI
xG3;}
case
cRSqrt:value.yH1;xP=true
iI
n71-0.5));case
cInv:value.yH1;xP=true
iI
n71-1));yT3
y73
return
n71
1));}
cB1
void
eL1
eM1&mul,l92&tree,l92&xS2,bool&c71
bool&xP){for
y11
c4++a){lM1
value(lA4)eO
xG3(nZ1
value,xP)l64!xS2.c42
xS2.xJ1!=nW3(1.0
iZ
cZ1;cZ1
l84
eB3
xG3
eB3
xS2);cZ1
tB1);xG3.swap(cZ1);}
#if 0 /* FIXME: This does not work */
cD1
nC==cMul){if(1){bool
exponent_is_even=xG3
cV2)&&isEvenInteger
e92.xJ1);e31
value.eX3{bool
tmp=false
eO
val(value
l8
b))eO
exp(nZ1
val,tmp)l64
exponent_is_even||(exp
cV2)&&isEvenInteger(exp.xJ1)iZ
cZ1;cZ1
l84
eB3
xG3);cZ1
c0
exp);cZ1.ConstantFolding(l64!cZ1.c42!isEvenInteger(cZ1.xJ1)){goto
cannot_adopt_mul;}
}
}
}
eL1
mul,value,xG3,c71
xP);}
else
cannot_adopt_mul:
#endif
{if(mul.x72(value,xG3)==CollectionSetBase::n61)c81}
}
}
y41
ConstantFolding_MulGrouping(eX{bool
xP
eV3
bool
should_regenerate
eV3
eM1
mul;eL1
mul
iK2,n71
1)),c71
xP)iD3
std::pair<lM1,eM>eN1
iD3
std::multimap<fphash_t,eN1>cM1;cM1
i9;tY3
eM1::xI1
j=mul.i8.yJ3
j!=mul.i8.end();++j){lM1&value=j
e52.value
eO&e62
j
e52.xS2;if(j
e52.e3)xG3
tB1);const
fphash_t
eO1=xG3.GetHash();xJ3
cM1::yK3
i=i9.xT2
eO1);for(;i!=i9.cY1
eO1;++i)if(i
e52.first
xL
xG3)){if(!xG3.c42!e01.xJ1
iM2(1)))c81
i
e52
t03.push_back(value);goto
skip_b;}
i9.yG3,std::make_pair(eO1,std::make_pair
e92,eM
eZ2(1),value))));skip_b:;}
#ifdef FP_MUL_COMBINE_EXPONENTS
ConstantExponentCollection
xH
e21;tY3
cM1::yK3
j,i=i9.yJ3
i!=i9.end();i=j){j=i;++j;eN1&list=i
e52;if(list
tM2
lS1
e62
list
tM2
xJ1;if(!e92==x61)e21.MoveToSet_Unique
e92,list
t03);i9.erase(i);}
}
if(e21.iI2)c81
#endif
if(should_regenerate){lM1
before=tree;before.l71
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_MulGrouping: "
c72
before)x71"\n"
;
#endif
tree.DelParams();tY3
cM1::yK3
i=i9.yJ3
i!=i9.end();++i){eN1&list=i
e52;
#ifndef FP_MUL_COMBINE_EXPONENTS
if(list
tM2
lS1
e62
list
tM2
xJ1;if
e92==x61
lI3
if(e01
nK2
tree.AddParamsMove(list
t03);t02}
#endif
lM1
mul;mul
l84);mul.SetParamsMove(list
t03);mul
tB1
l64
xP&&list
tM2
nJ2
list
tM2
xJ1==nW3(1)/nW3(3
iZ
cbrt;cbrt
eN
cCbrt);cbrt
xW2
cbrt
tB1)xY2
cbrt
y2
0.5
iZ
sqrt;sqrt
eN
cSqrt);sqrt
xW2
sqrt
tB1)xY2
sqrt
y2-0.5
iZ
rsqrt;rsqrt
eN
cRSqrt);rsqrt
xW2
rsqrt
tB1)xY2
rsqrt
y2-1
iZ
inv;inv
eN
cInv);inv
xW2
inv
tB1)xY2
inv);t02}
lM1
pow;pow
eN
cPow);pow
xW2
pow
c0
list.first);pow
tB1)xY2
pow);}
#ifdef FP_MUL_COMBINE_EXPONENTS
i9.clear()xZ1
tZ.size();++a
i12
e62
tZ[a
eU3;if(e01
nK2
tree.AddParamsMove(tZ[a
eZ3);t02
lM1
mul;mul
l84);mul.SetParamsMove(tZ[a
eZ3);mul
iR2
pow;pow
eN
cPow);pow
xW2
pow.yF
xG3));pow
tB1)xY2
pow);}
#endif
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_MulGrouping: "
tO
#endif
return!tree
xL
before);}
return
t43
y41
ConstantFolding_AddGrouping(eX{bool
should_regenerate
eV3
eM1
add
nY2
if(lA4
nC==cMul)lI3
if(add.x72(lA4)==CollectionSetBase::n61)c81}
e13
remaining
iJ)i42
tE=0
nY2
l92&xI3=lA4;if
lJ3
nC==cMul){e31
yC1
eX3{if
lJ3
l8
b)cV2))lI3
xJ3
eM1::xI1
c=add.FindIdenticalValueTo
lJ3
l8
b)l64
add.Found(c
iZ
tmp
lJ3
yE
CloneTag());tmp.l83
b
tM3
add.AddCollectionTo(tmp,c);c81
goto
done_a;}
}
remaining[a]=true;tE+=1;done_a:;}
}
if(tE>0){if(tE>1){std
xV3<std::pair<lM1,eY3>x4;std::multimap<fphash_t,eY3
eP1;bool
lW3
eV3
for
iF{e31
lA4.eX3{l92&p=lA4
l8
b);const
fphash_t
p_hash=p.GetHash();for(std::multimap<fphash_t,eY3::const_iterator
i=eP1.xT2
p_hash);i!=eP1.cY1
p_hash;++i){if(x4[i
e52
eU3
xL
p)){x4[i
e52
eZ3+=1;lW3=true;goto
found_mulgroup_item_dup;}
}
x4.push_back(std::make_pair(p,size_t(1)));eP1.insert(std::make_pair(p_hash,x4.size()-1));found_mulgroup_item_dup:;}
}
if(lW3){lM1
eB2;{size_t
max=0;eY2
p=0;p<x4.size();++p)if(x4[p
eZ3<=1)x4[p
eZ3=0;else{x4[p
eZ3*=x4[p]tM2
y22;if(x4[p
eZ3>max){eB2=x4[p
eU3;max=x4[p
eZ3;}
}
}
lM1
group_add;group_add
eN
cAdd);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Duplicate across some trees: "
c72
eB2)x71" in "
tO
#endif
for
iF
e31
lA4.eX3
if(eB2
xL
lA4
l8
b)iZ
tmp(lA4
yE
CloneTag());tmp.l83
b
tM3
group_add
c0
tmp);remaining[a]eV3
y73
group_add
iR2
group;group
l84);group
c0
eB2);group
c0
group_add);group
tB1);add.x72(group);c81}
}
for
iF{if(add.x72(lA4)==CollectionSetBase::n61)c81}
}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_AddGrouping: "
tO
#endif
tree.DelParams();tY3
eM1::xI1
j=add.i8.yJ3
j!=add.i8.end();++j){lM1&value=j
e52.value
eO&coeff=j
e52.xS2;if(j
e52.e3)coeff
tB1
l64
coeff.nJ2
n23
coeff.xJ1,x61)lI3
lL3
coeff.xJ1
nK2
tree
c0
value);t02}
lM1
mul;mul
l84);mul
c0
value);mul
c0
coeff);mul
tB1);tree
xW2}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_AddGrouping: "
tO
#endif
l12}
return
t43}
l93{using
l93
FUNCTIONPARSERTYPES;using
t7
y31
bool
ConstantFolding_IfOperations(iX3(tree.GetOpcode()==cIf
iM3()==cAbsIf);for(;;){if(cF3==cNot){tH2
cIf);xZ
0)eA2
0)xX2
xZ
1).swap(xZ
2));}
iH1
xZ
0)cS1{tH2
t33;xZ
0)eA2
0)xX2
xZ
1).swap(xZ
2));}
else
yB3
lX1
0),tG2==t33)lU2
tree
eA2
1));x5
l53
tree
eA2
2));x5
lY1
if(cF3==cIf||cF3==t33{lM1
cond=xZ
0)eO
lX3;lX3
tK2==cIf?cNotNot:cAbsNotNot);lX3
xZ2
1));ConstantFolding(lX3)eO
lY3;lY3
tK2==cIf?cNotNot:cAbsNotNot);lY3
xZ2
2));ConstantFolding(lY3
l64
lX3.c42
lY3
cV2
iZ
eY;eY
tK2);eY
xZ2
1));eY.nJ
1));eY.nJ
2));eY
iR2
eZ;eZ
tK2);eZ
xZ2
2));eZ.nJ
1));eZ.nJ
2));eZ
tB1);tH2
cond
nC);tL2
0,cond
xX2
iS2
1,eY);iS2
2,eZ
nZ2}
if(xZ
1)nC==xZ
2)nC&&(xZ
1)nC==cIf||xZ
1)nC==cAbsIf
iZ&iJ2=xZ
1)eO&leaf2=xZ
2
l64
iJ2
l8
0)x51
0))&&lK3
1))||iJ2
l8
2)x51
2))iZ
eY;eY
iM1
eY.nJ
0));eY
eQ
iJ2
l8
1));eY
eQ
leaf2
l8
1));eY
iR2
eZ;eZ
iM1
eZ.nJ
0));eZ
eQ
iJ2
l8
2));eZ
eQ
leaf2
l8
2));eZ
i0
SetParam(0,iJ2
xX2
iS2
1,eY);iS2
2,eZ
nZ2
if
lK3
1))&&iJ2
l8
2)x51
2)iZ
t0;t0
iM1
t0
c0
xZ
0));t0
eQ
iJ2
xX2
t0
eQ
leaf2
xX2
t0
i0
n81
0,t0);tL2
2,iJ2
l8
2));tL2
1,iJ2
l8
1)nZ2
if
lK3
2))&&iJ2
l8
2)x51
1)iZ
eC2;eC2
eN
leaf2
nC==cIf?cNot:cB3);eC2
eQ
leaf2
xX2
eC2
iR2
t0;t0
iM1
t0
c0
xZ
0));t0
eQ
iJ2
xX2
t0
c0
eC2);t0
i0
n81
0,t0);tL2
2,iJ2
l8
2));tL2
1,iJ2
l8
1)nZ2}
lM1&xY=xZ
1)eO&yA=xZ
2
l64
xY
xL
yA)){tree
eA2
1)nZ2
const
OPCODE
op1=xY
nC;const
OPCODE
op2=yA
nC;yL3
op2){if(xY
c91
1){lM1
lO
0
i01
0
lS2
n5
if(xY
c91
2&&yA
c91
2){if(xY
l8
0)xL
yA
l8
0)iZ
param0=xY
l8
0)eO
lO
1
i01
1
lS2
xY2
param0)n5
if(xY
l8
1)xL
yA
l8
1)iZ
param1=xY
l8
1)eO
lO
0
i01
0
lS2
xY2
xB1)xY2
param1
nZ2}
yL3
yM3
cMul
lT2
cAnd
lT2
cOr
lT2
cAbsAnd
lT2
cAbsOr
lT2
cMin
lT2
cMax){eM
lZ3;c7{eY2
b=yA.l21
b-->0;){if
lM3
yA
l8
b))){if(lZ3
tH3){xY.l71
yA.l71}
lZ3.push_back(xY
eJ3;yA.l83
b);xY
cK1
y73}
}
if(!lZ3
tH3){xY
tB1);yA
iR2
xB1;xB1
iM1
xB1.SetParamsMove(tree.l52
nX
SetParamsMove(lZ3)n5}
}
yL3
yM3
cMul||(op1==cAnd
nL1
yA))||(op1==cOr
nL1
yA))){c7
if
lM3
yA)){xY.l71
xY
cK1
xY
iR2
cN1=yA;yA=tF
op1==yM3
cOr)l9
op1)xY2
cN1)n5}
if((op1==cAnd
lT2
cOr)&&op2==cNotNot){lM1&n03=yA
l8
0);c7
if
lM3
n03)){xY.l71
xY
cK1
xY
iR2
cN1=n03;yA=tF
op1==cOr)l9
op1)xY2
cN1)n5}
if(op2==cAdd||op2==cMul||(op2==cAnd
nL1
xY))||(op2==cOr
nL1
xY))){eQ1
yA.cO1
yA
l8
a)xL
xY)){yA.l71
yA
cK1
yA
iR2
cP1=xY;xY=tF
op2==cAdd||op2==cOr)l9
op2)xY2
cP1)n5}
if((op2==cAnd||op2==cOr)&&op1==cNotNot){lM1&n13=xY
l8
0);eQ1
yA.cO1
yA
l8
a)xL
n13)){yA.l71
yA
cK1
yA
iR2
cP1=n13;xY=tF
op2==cOr)l9
op2)xY2
cP1)n5}
return
t43}
#include <limits>
l93{using
l93
FUNCTIONPARSERTYPES;using
t7
y31
int
maxFPExponent()yR
std::numeric_limits
xH::max_exponent;}
y41
x01
nW3
base
iM2
xG3){if(base<x61
l12
lL3
base,x61||n23
base
iM2(1))cV
return
xG3>=nW3(maxFPExponent
xH())/fp_log2(base);}
y41
ConstantFolding_PowOperations(iX3(tree.GetOpcode()==cPow);nT&&xZ
1).lS1
const_value=t53
lS,xZ
l81);y0
const_value)iI
t43
if(eG1
n23
xZ
l81
nK2
tree
eA2
0)nZ2
nT&&n23
lS
nK2
y0
1)iI
t43
nT&&xZ
1)nC==cMul){bool
y02=false
tV2
lV2=lS
eO
xI3=xZ
1);eQ1
yC1
cO1
xI3
l8
a).lS1
imm=xI3
l8
a).xJ1;{if(x01
lV2,imm))break
tV2
lW2=t53
lV2,imm);lL3
lW2,x61)break;if(!y02){y02=true;yC1
l71}
lV2=lW2;yC1
l83
a
c52}
if(y02){yC1
Rehash();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before pow-mul change: "
tO
#endif
xZ
0).Become(e91
lV2));xZ
1).Become
lJ3);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After pow-mul change: "
tO
#endif
}
}
if(eG1
cF3==cMul
i12
lX2=xZ
l81
tV2
y12=1.0;bool
y02=false
eO&xI3=xZ
0);eQ1
yC1
cO1
xI3
l8
a).lS1
imm=xI3
l8
a).xJ1;{if(x01
imm,lX2))break
tV2
eX1=t53
imm,lX2);lL3
eX1,x61)break;if(!y02){y02=true;yC1
l71}
y12*=eX1;yC1
l83
a
c52}
if(y02){yC1
Rehash()eO
cY3;cY3
eN
cPow);cY3.SetParamsMove(tree.l52));cY3
tY2;tH2
cMul)xY2
cY3);tree
eQ
e91
y12)nZ2}
if(cF3==cPow&&eG1
xZ
0)l8
1).lS1
a=xZ
0)l8
l81
tV2
b=xZ
l81
tV2
c=a*b;if(isEvenInteger(a)&&!isEvenInteger(c
iZ
n33;n33
l94
n33.nJ
0)xX2
n33
tB1);iS2
0,n33);}
else
tL2
0,xZ
0)xX2
tL2
1,e91
c));}
return
t43}
l93{using
l93
FUNCTIONPARSERTYPES;using
t7;cZ2
l4{enum
eD2{MakeFalse=0,MakeTrue=1,tF2=2,n63=3,MakeNotNotP0=4,MakeNotNotP1=5,MakeNotP0=6,MakeNotP1=7,xM=8
eA3
lY2{Never=0,Eq0=1,Eq1=2,yN3=3,yO3=4}
;eD2
if_identical;eD2
lZ2
4];cZ2{eD2
what:4;lY2
when:4;}
iP1,iQ1,iR1,iS1
y31
eD2
Analyze
e42
a,l92&b
e93{if(a
xL
b
tQ2
if_identical;yP3
p0=iO
a);yP3
p1=iO
b
l64
p0
xK
cH3
p1
e61){if(p0
xK
val<p1
eE2&&lZ2
0]i1
0];if(p0
xK
val<=p1
eE2&&lZ2
1]i1
1];}
if(p0
e61&&p1
eF2{if(p0
eE2>p1
xK
val&&lZ2
2]i1
2];if(p0
eE2>=p1
xK
val&&lZ2
3]i1
3];}
if(IsLogicalValue(a)){if(iP1
eJ2
iP1.when,p1
tQ2
iP1.what;if(iR1
eJ2
iR1.when,p1
tQ2
iR1.what;}
if(IsLogicalValue(b)){if(iQ1
eJ2
iQ1.when,p0
tQ2
iQ1.what;if(iS1
eJ2
iS1.when,p0
tQ2
iS1.what;}
return
xM;}
cB1
bool
TestCase(lY2
when,const
yP3&p){if(!p
e61||!p
xK
known
cV
eL3
when
cZ3
Eq0
nN1==nW3(0.0)iX2==p
eE2;case
Eq1
nN1==nW3(1.0)iX2==p
xK
val;case
yN3
nN1>y71
iX2<=nW3(1);case
yO3
nN1>=y71
lF1
1);yT3;}
return
t43}
;l93
RangeComparisonsData{static
const
l4
Data[6]={{l4
n43
tR
xM,l4::tR
xM
lU1
Eq1
lV1
Eq1
nO1
Eq0}
iS
Eq0}
}
,{l4::x82
n53
xM,l4
n53
xM
lU1
Eq0
lV1
Eq0
nO1
Eq1}
iS
Eq1}
}
,{l4::x82
n53
tF2,l4::tR
MakeFalse
nO1
yN3
lV1
yO3
yX,{l4
n43
xM,l4
n53
tR
n63
nO1
yO3
lV1
yN3
yX,{l4::x82::tR
tR
MakeTrue,l4::tF2
lU1
yO3}
iS
yN3
yX,{l4
n43
tR
n63,l4::xM,l4
n91
lU1
yN3}
iS
yO3
yX}
;}
y41
ConstantFolding_Comparison(eX{using
l93
RangeComparisonsData;assert(tree.GetOpcode()>=cEqual&&tree.GetOpcode()<=cGreaterOrEq);eL3
Data[tG2-cEqual].Analyze(xZ
0),xZ
1))cZ3
l4::MakeFalse:y0
0);x5
l4
n91:y0
1
tP2
n63:tH2
cEqual
tP2
tF2:tH2
i61
tP2
MakeNotNotP0:tH2
cNotNot
yI1
1
tP2
MakeNotNotP1:tH2
cNotNot
yI1
0
tP2
MakeNotP0:tH2
cNot
yI1
1
tP2
MakeNotP1:tH2
cNot
yI1
0
tP2
xM:;}
if(xZ
cT2))eL3
cF3
cZ3
cAsin:lM
fp_sin(xZ
eK2
cAcos:lM
fp_cos(xZ
l81)));tH2
tG2==cLess?cGreater:tG2==cLessOrEq?cGreaterOrEq:tG2==cGreater?cLess:tG2==cGreaterOrEq?cLessOrEq:tG2);x5
cAtan:lM
fp_tan(xZ
eK2
cLog:lM
fp_exp(xZ
eK2
cSinh:lM
fp_asinh(xZ
eK2
cTanh:if(fp_less(fp_abs(xZ
l81)nK2
lM
fp_atanh(xZ
l81))nZ2
break;yT3
y73
return
t43}
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;l93{
#ifdef DEBUG_SUBSTITUTIONS
yH
double
d){union{double
d;uint_least64_t
h
yA2
d=d;lJ1
h
nP1
#ifdef FP_SUPPORT_FLOAT_TYPE
yH
float
f){union{float
f;uint_least32_t
h
yA2
f=f;lJ1
h
nP1
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
yH
long
double
ld){union{long
double
ld;cZ2{uint_least64_t
a;xH3
short
b;}
s
yA2
ld=ld;lJ1
s.b<<data.s.a
nP1
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
yH
long
ld){o<<"("
<<std::hex<<ld
nP1
#endif
#endif
}
t7{lN
nF)){}
lN
const
nW3&i
yE
xN3
nF
i
x92
#ifdef FP_SUPPORT_CXX11_MOVE
lN
nW3&&i
yE
xN3
nF
std::move(i)x92
#endif
lN
xH3
v
yE
VarTag
nF
l13,v
x92
lN
cU2
o
yE
OpcodeTag
nF
o
x92
lN
cU2
o,xH3
f
yE
FuncOpcodeTag
nF
o,f
x92
lN
l92&b
yE
CloneTag
nF*b.data)){}
nQ2
lM1::~lO1(){}
lC
ReplaceWithImmed
cH1
i){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Replacing "
c72*this
l64
IsImmed())OutFloatHex(std::cout,xJ1)x71" with const value "
<<i;OutFloatHex(std::cout,i)x71"\n"
;
#endif
data=new
xK2
xH(i);}
nQ2
cZ2
ParamComparer{iL2()e42
a,l92&b
e93{if(a.y22!=b.y22
lD4
a.y22<b.y22
iI
a.GetHash()<b.GetHash();}
}
y31
void
xK2
xH::Sort(){eL3
Opcode
cZ3
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cHypot:case
cEqual:case
i61:std::sort(iZ2
iQ2
iZ2
end(),ParamComparer
xH());lD
cLess
lZ
cGreater;}
lD
cLessOrEq
lZ
cGreaterOrEq;}
lD
cGreater
lZ
cLess;}
lD
cGreaterOrEq
lZ
cLessOrEq;}
break;yT3
y73}
lC
AddParam
e42
param){y3.push_back(param);}
lC
eB
lM1&param){y3.push_back(lM1());y3.back().swap(param);}
lC
SetParam
eZ2
which,l92&b)nQ1
which
iY2
y3[which]=b;}
lC
n81
size_t
which,lM1&b)nQ1
which
iY2
y3[which
tT3
b);}
lC
AddParams
cI1
nK){y3.insert(y3.end(),n02.iQ2
n02.end());}
lC
AddParamsMove(nK){size_t
endpos=y3.size(),added=n02.size();y3.eP3
endpos+added,lM1());eY2
p=0;p<added;++p)y3[endpos+p
tT3
n02[p]);}
lC
AddParamsMove(nK,size_t
n12)nQ1
n12
iY2
l83
n12);AddParamsMove(tS1}
lC
SetParams
cI1
nK){eM
tmp(tS1
y3.swap(tmp);}
lC
SetParamsMove(nK){y3.swap(tS1
n02.clear();}
#ifdef FP_SUPPORT_CXX11_MOVE
lC
SetParams(eM&&n02){SetParamsMove(tS1}
#endif
lC
DelParam
eZ2
index){eM&Params=y3;
#ifdef FP_SUPPORT_CXX11_MOVE
iZ2
erase(iZ2
begin()+index);
#else
Params[index].data=0;eY2
p=index;p+1<y93;++p)Params[p].data.UnsafeSetP(&*Params[p+1
iY2
Params[y93-1].data.UnsafeSetP(0);iZ2
eP3
y93-1);
#endif
}
lC
DelParams(){y3.clear();}
y41
lM1::IsIdenticalTo
e42
b
e93{if(&*data==&*b.data
lD4
true
iI
data->IsIdenticalTo(*b.data);}
y41
xK2
xH::IsIdenticalTo
cI1
xK2
xH&b
e93{if(Hash!=b.Hash
cV
if(Opcode!=tU3
cV
eL3
Opcode
cZ3
cImmed:return
n23
Value,tV3;case
l13:return
l72==b.l62
case
cFCall:case
cPCall:if(l72!=b.l72
cV
break;yT3
y73
if(y93!=b.y93
cV
for
y11
y93
tO3
if(!Params[a]xL
b.Params[a])cV}
l12}
lC
Become
e42
b){if(&b!=this&&&*data!=&*b.data){DataP
tmp=b.data;l71
data.swap(tmp);}
}
lC
cN2(){if(GetRefCount()>1)data=new
xK2
xH(*data);}
nQ2
lM1
lM1::GetUniqueRef(){if(GetRefCount()>1
lD4
lM1(*this,CloneTag())iI*this;}
tS
yI
cNop),Value(),n9
tS
const
xK2&b
yI
tU3),Value(tV3,l72(b.cQ1,Params(b.Params),Hash(b.Hash),Depth(b.Depth),i11
b.l82){}
tS
const
nW3&i
yI
cImmed),Value(i),n9
#ifdef FP_SUPPORT_CXX11_MOVE
tS
xK2
xH&&b
yI
tU3),Value
c92
tV3),l72(b.cQ1,Params
c92
b.Params)),Hash(b.Hash),Depth(b.Depth),i11
b.l82){}
tS
nW3&&i
yI
cImmed),Value
c92
i)),n9
#endif
tS
cU2
o
yI
o),Value(),n9
tS
cU2
o,xH3
f
yI
o),Value(),l72(f),Params(),Hash(),Depth(1),i11
0){}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <iostream>
using
l93
FUNCTIONPARSERTYPES;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
l93{eV1
i21
nR,std
cD&done,std::ostream&o){for
y11
c4++a)i21
lA4,done,o);std::ostringstream
buf
c72
tree,buf);done[tree.GetHash()].insert(buf.str());}
}
#endif
t7{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
eV1
DumpHashes(cI){std
cD
done;i21
tree,done,o);for(std
cD::const_iterator
i=done.yJ3
i!=done.end();++i){const
std::set<std
eC3>&flist=i
e52;if(flist.size()!=1)o<<"ERROR - HASH COLLISION?\n"
;for(std::set<std
eC3>::const_iterator
j=flist.yJ3
j!=flist.end();++j){o<<'['<<std::hex<<i->first.hash1<<','<<i->first.hash2<<']'<<std::dec;o<<": "
<<*j<<"\n"
;}
}
}
eV1
DumpTree(cI){t52
iR3;eL3
tG2
cZ3
cImmed:o<<l54
t73
l13:o<<"Var"
<<(tree.GetVar()-l13)t73
cAdd:iR3"+"
;lD
cMul:iR3"*"
;lD
cAnd:iR3"&"
;lD
cOr:iR3"|"
;lD
cPow:iR3"^"
;break;yT3
iR3;o<<t13
tG2);eG2
cFCall||tG2==cPCall)o<<':'<<tree.GetFuncNo();}
o<<'(';if
iJ<=1&&sep2[1])o<<(sep2+1)<<' '
nY2
if(a>0)o<<' '
c72
lA4,o
l64
a+1<tree.GetParamCount())o<<sep2;}
o<<')';}
eV1
DumpTreeWithIndent(cI,const
std
eC3&indent){o<<'['<<std::hex<<(void*)(&tree.l52))<<std::dec<<','<<tree.GetRefCount()<<']';o<<indent<<'_';eL3
tG2
cZ3
cImmed:o<<"cImmed "
<<l54;o<<'\n'
t73
l13:o<<"VarBegin "
<<(tree.GetVar()-l13);o<<'\n'
iI;yT3
o<<t13
tG2);eG2
cFCall||tG2==cPCall)o<<':'<<tree.GetFuncNo();o<<'\n';}
for
y11
c4++a){std
eC3
ind=indent;eY2
p=0;p<ind.size();p+=2)if(ind[p]=='\\')ind[p]=' ';ind+=(a+1<tree.GetParamCount())?" |"
:" \\"
;DumpTreeWithIndent(lA4,o,ind);}
o<<std::flush;}
#endif
}
#endif
using
l93
l31;using
l93
FUNCTIONPARSERTYPES;
#include <cctype>
l93
l31{xH3
ParamSpec_GetDepCode
cI1
e22&b){eL3
b.first
cZ3
t23:{cQ*s=(cQ*)b
t03
iI
s->depcode;}
case
SubFunction:{cR*s=(cR*)b
t03
iI
s->depcode;}
yT3
y73
return
0;}
eV1
DumpParam
cI1
e22&yF2
std::ostream&o){static
const
char
ParamHolderNames[][2]={"%"
,"&"
,"x"
,"y"
,"z"
,"a"
,"b"
,"c"
}
;xH3
y32
0;e53
e63{const
ParamSpec_NumConstant
xH
yV3
cI1
ParamSpec_NumConstant
xH*lH2;using
l93
FUNCTIONPARSERTYPES;o.precision(12);o<<eG3
constvalue;y73
case
t23:{cQ
yV3(cQ*lH2;o<<ParamHolderNames[eG3
index];y32
eG3
constraints;y73
case
SubFunction:{cR&param
yC3;y32
eG3
constraints;yJ
GroupFunction){if(param
iT==cNeg){o<<"-"
;n3}
iH1
param
iT==cInv){o<<"/"
;n3}
else{std
eC3
opcode=t13(cU2)param
iT).substr(1)xZ1
opcode.size();++a)opcode[a]=(char)std::toupper(opcode[a]);o<<opcode<<"( "
;n3
o<<" )"
;}
}
else{o<<'('<<t13(cU2)param
iT)<<' ';yJ
PositionalParams)o<<'[';yJ
SelectedParams)o<<'{';n3
if
e43
data.n2!=0)o<<" <"
<<eG3
data.n2<<'>';yJ
PositionalParams)o<<"]"
;yJ
SelectedParams)o<<"}"
;o<<')';}
y73
cA3(ImmedConstraint_Value(constraints&ValueMask)cZ3
ValueMask:lD
Value_AnyNum:lD
xA2:o<<"@E"
;lD
Value_OddInt:o<<"@O"
;lD
i41:o<<"@I"
;lD
Value_NonInteger:o<<"@F"
;lD
eY1:o<<"@L"
;yB3
ImmedConstraint_Sign(constraints&SignMask)cZ3
SignMask:lD
Sign_AnySign:lD
nA1:o<<"@P"
;lD
eZ1:o<<"@N"
;yB3
ImmedConstraint_Oneness(constraints&OnenessMask)cZ3
OnenessMask:lD
Oneness_Any:lD
Oneness_One:o<<"@1"
;lD
Oneness_NotOne:o<<"@M"
;yB3
ImmedConstraint_Constness(constraints&ConstnessMask)cZ3
ConstnessMask:lD
i31:if(nJ3.first==t23){cQ
yV3(cQ*lH2;if
e43
index<2)y73
o<<"@C"
;lD
Constness_NotConst:o<<"@V"
;lD
Oneness_Any:y73}
eV1
DumpParams
x91
paramlist,xH3
count,std::ostream&o){for
x91
a=0;a<count
tO3
if(a>0)o<<' ';const
e22&param=eE1
xH(paramlist,a);DumpParam
xH(param,o);xH3
depcode=ParamSpec_GetDepCode(param
l64
depcode!=0)o<<"@D"
<<depcode;}
}
}
#include <algorithm>
using
l93
l31;using
l93
FUNCTIONPARSERTYPES;l93{cQ
plist_p[37]={{2,0,0x0}
x0
0,0x4}
x0
nA1,0x0}
x0
eZ1|Constness_NotConst,0x0}
x0
Sign_NoIdea,0x0}
x0
eY1,0x0}
,{3,Sign_NoIdea,0x0}
,{3,0,0x0}
,{3,eY1,0x0}
,{3,0,0x8}
,{3,Value_OddInt,0x0}
,{3,Value_NonInteger,0x0}
,{3,xA2,0x0}
,{3,nA1,0x0}
,{0,eZ1|lW{0,lW{0,nA1|lW{0,xA2|lW{0,i31,0x1}
,{0,i41|nA1|lW{0,i51
i31,0x1}
,{0,i51
lW{0,Oneness_One|lW{0,eY1|lW{1,lW{1,xA2|lW{1,i51
lW{1,i41|lW{1,nA1|lW{1,eZ1|lW{6,0,0x0}
,{4,0,0x0}
,{4,i41,0x0}
,{4,lW{4,0,0x16}
,{5,0,0x0}
,{5,lW}
y31
cZ2
plist_n_container{static
const
ParamSpec_NumConstant
xH
plist_n[20];}
y31
const
ParamSpec_NumConstant
xH
plist_n_container
xH::plist_n[20]={{nW3(-2
tT-1
tT-0.5
tT-0.25
tT
0
tN2
fp_const_deg_to_rad
xH(tN2
fp_const_einv
xH(tN2
fp_const_log10inv
xH(tT
0.5
tN2
fp_const_log2
xH(tT
1
tN2
fp_const_log2inv
xH(tT
2
tN2
fp_const_log10
xH(tN2
fp_const_e
xH(tN2
fp_const_rad_to_deg
xH(tN2-fp_const_pihalf
xH(),xK1{y71,xK1{fp_const_pihalf
xH(),xK1{fp_const_pi
xH(),xK1}
;cR
plist_s[517]={{{1,15,n73,398,n73,477,n73,15,cNeg,GroupFunction,0}
,i31,0x1
tO2
15,y42
24,y42
465,y42
466,y42
498,cInv,lU
2,327995
c2
l0
2,48276
c2
l6
260151
c2
l6
470171
c2
l6
169126
c2
l6
48418
c2
l6
1328
c2
l6
283962
c2
l6
169275
c2
l6
39202
c2
l6
283964
c2
l6
283973
c2
l6
476619
c2
l6
296998
c2
l6
47
c2
SelectedParams,0}
,0,0x4
nM
161839
c2
l6
25036
c2
l6
35847
c2
l6
60440
c2
l6
30751
c2
l6
183474
c2
l6
259318
c2
l6
270599
c2
l6
60431
c2
l6
259119
c2
l6
332066
c2
l6
7168
c2
l6
197632
c2
l6
291840
c2
l6
283648
c2
l6
238866
c2
l6
239902
c2
l6
31751
c2
l6
244743
c2
l6
384022
c2
SelectedParams,0}
,0,0x4
nM
385262
c2
l6
386086
c2
l6
393254
c2
SelectedParams,0}
,0,0x5
nM
393254
c2
l6
386095
c2
l6
387312
c2
l6
18662
c2
l6
61670
c2
l6
387397
c2
l6
247855
c2
SelectedParams,0}
,0,0x1
nM
342063
c2
l6
297007
c2
l6
15820
c2
l6
393263
c2
l6
393263
c2
SelectedParams,0}
,0,0x5
nM
161847
c2
l6
258103
c2
l6
249073
c2
l6
249076
c2
iA
0,0
c2
nL
0,0
c2
cA2
1,45
c2
nL
1,53
c2
nL
1,54
c2
nL
1,55
c2
nL
1,56
c2
nL
1,26
c2
nL
1,259
c2
eL2
0x16
tO2
272
c2
cA2
1,323
c2
eL2
0x16
tO2
0
c2
nL
1,21
c2
nL
1,447
c2
eL2
0x4
tO2
449
c2
eL2
0x4
tO2
0
c2
eL2
0x4
tO2
0
c2
tU
2}
,0,0x4
tO2
15
c2
nL
1,24
c2
tU
2}
,0,0x0
nM
58392
c2
cA2
0,0
c2
tU
1}
,nA1,0x0
nM
24591
tP3
33807
tP3
48143
tP3
285720
tP3
290840
tP3
305152,lA
312400,lA
39202,lA
122918,lA
421926,lA
429094,lA
443430,lA
317834,lA
329098,lA
7633,lA
7706,lA
7730,lA
38,lA
50587,lA
406528,lA
24583,lA
31751,lA
405511,lA
321551,xM1
327713,lA
322596,lA
90409,lA
335174,lA
327050,lA
493606,lA
496678,lA
503846,lA
516134,lA
7217,lA
333875,lA
336896,lA
524326,lA
509952,lA
286727,lA
89103,lA
92175,lA
296976,tU1
324623,l1
0x14
nM
332815,l1
0x10}
,{{3,7340056,tU1
289092,lA
93200,xM1
337935
iU2,7340060,l1
tR2
7340176,lA
338959
iU2,7340061,xM1
7206,lA
yW3
lA
357414,lA
368678,lA
370745,l1
0x7}
,{{3,7340177,lA
39277,tU1
426398,l1
tR2
40272286,xM1
490910,l1
tR2
40336798,xM1
50600,lA
426462,xM1
490974,xM1
370726,l1
0x6
nM
371750,l1
0x6
nM
428070
iU2,40336862,xM1
38378,lA
50671
iU2,47662080,lA
477184,lA
568320,lA
371727,l1
0x7}
,{{3,15779306,lA
370703,l1
0x7
nM
39277,lA
39279,l1
0x4}
,{{3,15779238,lA
39338,tU1
436262,lA
508966,lA
39409,tU1
296998,tU1
35847,lA
15,tU1
377894,lA
386063,l1
0x1
nM
15,lA
7192,lA
123928,lA
122904,lA
30751,lA
57,lA
7456,lA
15674
iU2,67579935,lA
39237,lA
58768,lA
62924,lA
122880,lA
15760
iU2,64009216,l1
xL1
0,0,t01
xL1
0,0,t21
2,t01
0x4
tO2
2,t41
3,t01
0x4
tO2
3,t41
38
xN
38,t21
14
xN
57
xN
16,t01
0x0
nM
471103,t01
0x1
tO2
303
xN
323,lJ
2}
,0,0x0
nM
471363,t01
0x16
tO2
293,t01
0x4
tO2
294,t41
295
xN
296,t21
400
xN
0
xN
460
xN
465
xN
16,t01
0x1
tO2
57,lJ
2}
,0,0x1
tO2
0,t21
21
xN
15,t01
0x0
nM
24591
xN
24,t21
517,lJ
2}
,0,0x0
nM
46095,lK
46104,lK
15397,lK
287789,lK
66584,lK
404763,lK
62504,lK
15409,lK
39951,lK
24591,lK
33807,lK
50200,lK
62509,lK
50176,lG,178176
eR1
tA3
283648,lG,19456,lG,27648,lG,91136,lG,86016,lG,488448,lG,14342,lG,58375,lG,46147
nS
46151,lG,284679,lG,7183,lG,46159
nS
38993
nS
50262,lG,50249,lG,283808,lG,284835,lG,24822,lG,10240,lG,11264,lG,7170,lG,yW3
lG,17408,lG,164864,lG,237568,lG,242688
eR1
0x14
nM
476160,lG,25607,lG,122895,lG,50252,lG,39374,lG,50183,lG,7192,lG,122911,lG,252979,lG,46155,lG,38919,lG,50268,lG,50269,lG,50253,lG,46191,lG,50296,lG,7563
eR1
0x10
nM
416811,lG,416819,lG,40047,lG,46192
nS
415795,lG,40048
nS
415787,lG,39016
eR1
0x5
nM
39326
nS
39326,lG,39332
eR1
0x5
nM
39333
eR1
0x1
nM
50590
nS
50590,lG,39338
nS
39338,lG,39335
eR1
0x5
nM
15786
nS
146858,lG,39372,lG,39379,lG,39380,lG,39390
nS
50654
nS
50654,lG,24
eR1
0x6
nM
62,lG,24,lG,62
eR1
0x6
nM
43,lG,43
nS
51,lG,51
nS
50270,lG,50176
nS
50271,lG,39159,lG,39183
nS
7168
nS
31744,lG,100352,lG,31746,lG,101400,lG,39409
nS
39411
nS
39411,lG,39420,lG,39420
nS
15,lG,39026
eR1
0x5
nM
39422,lG,16384,lG,62853,lG,15360,lG,15
eR1
0x1
nM
16,lG,7183
eR1
0x1
nM
7172
t51
yG1,nA1,0x0
nM
24591
t51
lU
2,50200
t51
lU
2,63521
t51
lU
2,62500
t51
lU
2,50453
t51
lU
2,62488
t51
lU
1,0,t83
7,t83
194,t83
0,cAcos,t93
cAcosh,t93
cAsin,t93
cAsinh,nW
120,cAsinh,t93
cAtan
e81
306176,cAtan2
e81
yW3
cAtan2,t93
cAtanh,nW
246,cCeil,t93
cCeil,eS1
0,cB2
0,cCos,eS1
7,cB2
92,cB2
93,cB2
120,cB2
236,cB2
255,cB2
214,l23
236,l23
464,l23
0,cCosh,eS1
0,l23
0,cExp,nW
7,cExp,nW
92,cExp,t93
yX3
7,yX3
92,yX3
246,cFloor,t93
cFloor,lB
0x4
nM
309540,cHypot
e81
316708,cHypot
e81
316724,cHypot,l0
3,32513024,xN1
34627584
lN1
31493120,xN1
89213952
lN1
149042176
lN1
246647808
lN1
301234176
lN1
494360576
lN1
498558976
lN1
62933520
lN1
62933520,xN1
62933526
lN1
62933526,xN1
24670208
lN1
579378176
lN1
573578240
lN1
32513024
lN1
566254592
lN1
7900160
lN1
588822528,cIf,nW
120,cInt,nW
246,tS2
0,tS2
7,tS2
31,tS2
194,tS2
363,tS2
15,cLog,lU
1,24,cLog,lU
1,0,cLog10,t93
cLog2
e81
yW3
cMax
e81
35847,cMax
e81
30751,cMax,t93
cMax,eL2
0x4
nM
yW3
cMin
e81
35847,cMin
e81
30751,cMin,t93
cMin,eL2
0x4
nM
24591,cMin,lU
1,0,xB2
7,xB2
92,xB2
93,xB2
120,xB2
149,xB2
231,cSin,lB
0x5
tO2
246,xB2
255,xB2
254,xB2
0,cSin,eS1
273,cSin,lB
0x1
tO2
214,y52
231,cSinh,lB
0x5
tO2
246,y52
254,y52
255,y52
464,y52
0,cSinh,eS1
0,y52
15,cSqrt,lU
1,0,cC2
0,cTan,eS1
116,cTan,eS1
117,cC2
231,cC2
246,cC2
273,cC2
254,cC2
255,cC2
0,y62
0,cTanh,eS1
213,y62
231,y62
246,y62
254,y62
255,y62
0,cTrunc
e81
15384,cSub,lU
2,15384,cDiv,lU
2,476626,cDiv,lU
2,122937,tT2
yW3
n83
tA3
yW3
tT2
31744,n83
0x20
nM
31751,n83
0x24
nM
31751,tT2
122937,i61
e81
yW3
cLess,lB
tA3
41984,cLess,lB
0x4
nM
41984,cLess
e81
7,cLess
e81
yW3
cLessOrEq
e81
296182,cLessOrEq
e81
7168
eM2,lB
tA3
41984
eM2,lB
0x4
nM
41984
eM2
e81
7
eM2
e81
7168
iN1
e81
296182
iN1,t93
n22
245,n22
7,n22
550,n22
553,n22
554,n22
556,n22
31,n22
559,n22
15,n22
560,cNot
e81
7706,n93
yW3
n93
35847,n93
30751,n93
463903,n93
466975,cAnd,iA
0,0,cAnd,nL
2,yW3
eH3
7706,eH3
35847,eH3
463903,eH3
466975,eH3
30751,cOr,iA
1,0,n32
92,n32
131,n32
245,n32
215,n32
246,cDeg,nW
246,cRad
e81
yW3
cAbsAnd,l6
yW3
cAbsOr,iA
1,0,cB3,t93
cAbsNotNot,l0
3,32513024,c13
lB
0x0}
,}
;}
l93
l31{const
Rule
grammar_rules[262]={{ProduceNewTree,17,1,0,{1,0,cAbs,eN2
409,{1,146,cAtan,eN2
403
x0
1324,cAtan2,eN2
405
x0
307201,cAtan2
c1
253174
x0
255224,cAtan2
c1
259324
x0
257274,cAtan2,eN2
152,{1,252,cCeil,y82
486,{1,68,xO1
482,{1,123,xO1
483,{1,125,xO1
151,{1,126,xO1
419,{1,124,xO1
0,{1,403,cCos,l2
2,1,246,{1,252,cCos,l2
18,1,0,{1,400,xO1
301,{1,406,cCosh,l2
2,1,246,{1,252,cCosh,l2
18,1,0,{1,400,cCosh,y82
458,{1,122,cFloor,eN2
150,{1,252,cFloor,tB3
156,{3,7382016,eF
549,{3,8430592,eF
556,{3,8436736,eF
157,{3,42998784,eF
550,{3,42999808,eF
562,{3,43039744,eF
557,{3,49291264,eF
538,{3,49325056,eF
469,{3,1058318,eF
473,{3,1058324,eF
473,{3,9438734,eF
469,{3,9438740,cIf,l2
0,3,32542225,{3,36732434,cIf,l2
0,3,32542231,{3,36732440,cIf,tC3
573,{3,32513026,cIf,tC3
515,{3,455505423,cIf,tC3
515,{3,433506837,cIf,y82
78,{1,256,y72
69,{1,258,y72
404,{1,72,y72
159,{1,147,cLog,l2
0,1,0
x0
487425,cMax,tL
16,1,445
x0
yY3
cMax,tL
0,1,0
x0
483329,cMin,tL
16,1,446
x0
yY3
cMin,c8
0,1,153
x0
24832
t51
tB3
153
x0
25854
t51
tB3
154
x0
130063
t51
xP1
32055
t51
xP1
32056
t51
xP1
32057
t51
l2
0,2,166288
x0
32137
t51
xP1
33082
t51
l2
0,2,7168
x0
12688
t51
l2
0,2,7434
x0
12553
t51
y82
435
x0
46146
t51
y82
436
x0
46154
t51
y82
437
x0
46150
t51
y82
169
x0
83983
t51
y82
168
x0
131106
t51
y82
175
x0
133154
t51
iC3
476160
x0
471055
t51
iC3
274432
x0
273423
t51
iC3
251904
x0
266274
t51
iC3
251904
x0
263186
t51
y82
171,{1,252,lT1
421,{1,68,lT1
151,{1,123,lT1
419,{1,125,lT1
170,{1,126,lT1
482,{1,124,lT1
0,{1,405,lT1
172,{1,252,cSinh,y82
328,{1,404,cSinh,y82
173,{1,252,yZ3
0,{1,408,yZ3
176,{1,410,yZ3
177,{1,252,cTanh,l2
0,1,442
x0
449551,nG3
441
x0
yY3
nG3
167
x0
268549,nG3
180
x0
276749,nG3
181
x0
276500,nH3
190770
x0
189622,nH3
194748
x0
193723,nH3
202943
x0
196795,nH3
59699
x0
298148,nH3
59714
x0
325815,nH3
59724
x0
343224
c2
c8
2,1,337,{1,333
c2
tU
1
tG
336,{1,338
c2
tU
1}
}
,{ReplaceParams,2,1,340
x0
1363
c5
342
x0
1365
c5
463
x0
472524
c5
47
x0
356711
c5
349
x0
200751
c5
360
x0
199727
c5
480
x0
207053
c5
481
x0
208077
c5
417
x0
211144
c5
209
x0
211145
c5
418
x0
215240
c5
212
x0
212329
c5
204
x0
373097
c5
211
x0
372944
c5
217
x0
201944
c5
221
x0
223448
c5
367
x0
508329
c5
219
x0
508126
c5
224
x0
225705
c5
223
x0
225776
c5
365
x0
230825
c5
426
x0
377057
c5
497
x0
377054
c5
497
x0
204201
c5
426
x0
375280
c5
224
x0
375006,l7
2,2,407781
x0
233698,l7
2,2,59763
x0
233842,nG3
372
x0
1397,cD2
96
x0
24705,cD2
97
x0
24708,cD2
444
x0
449551,cD2
443
x0
yY3
cD2
101
x0
102774,cD2
109
x0
107845,cD2
106
x0
104773,l5
0,2,111631
x0
109893,l5
0,2,108559
x0
110917,lJ
0
tG
113
x0
112658,cMul,SelectedParams,0
tG
567,{1,52,lJ
1
tG
568,{1,42,lJ
1}
}
,{ReplaceParams,2,1,467
x0
45516
tI1
356
x0
51555
tI1
468
x0
49612
tI1
357
x0
47459
tI1
429
x0
438699
tI1
432
x0
441774
tI1
486
x0
498726
tI1
494
x0
504870
tI1
382
x0
435579
tI1
497
x0
435709
tI1
426
x0
508287
tI1
414
x0
500092
tI1
499
x0
352744
tI1
345
x0
367092
tI1
381
x0
425318
tI1
478
x0
425460
tI1
47
x0
512501
tI1
505
x0
355817
tI1
47
x0
516598
tI1
507
x0
518182
tI1
508
x0
358896
tI1
351
x0
388605
tI1
511
x0
360939
tI1
503
x0
354788
tI1
514
x0
525350
tI1
510
x0
394343
tI1
386
x0
351347,l5
2,2,363004
x0
361968,l5
16,1,118
x0
1157,l5
16,1,119
x0
1158,l5
16,1,402
x0
411024,l5
16,2,58768
x0
1472,l5
16,2,15760
x0
1474,l5
17,1,0,{1,400,l5
17,1,57,{1,14,lJ
0}
}
,{ProduceNewTree,4
iZ3
41
eI3,l3
4,1,0
x0
5167
eI3,cM
41984
x0
409641
eI3
cW
cEqual
cX
cEqual
cY
cEqual
xQ1
24849
eI3
c1
tM
cEqual
c1
n42
281873
eI3
c1
iB
cEqual
c1
lA1
cEqual,l3
4,1,562
x0
41,i61,l3
4
iZ3
5167,i61,cM
41984
x0
409641,i61
cW
i61
cX
i61
cY
i61
xQ1
24849,i61
c1
tM
i61
c1
n42
281873,i61
c1
iB
i61
c1
lA1
i61
cW
cLess
cX
cLess
cY
cLess,eN2
571
x0
46080,cLess
xQ1
24832,cLess
c1
xR1
cLess
c1
tM
cLess
c1
n42
281856,cLess
c1
nR1
cLess
c1
iB
cLess
c1
lA1
cLess,l3
20,1,562
x0
409641,cLess
cW
cLessOrEq
cX
cLessOrEq
cY
cLessOrEq,eN2
565
x0
409615,cLessOrEq
xQ1
24832,cLessOrEq
c1
xR1
cLessOrEq
c1
tM
cLessOrEq
c1
n42
281856,cLessOrEq
c1
nR1
cLessOrEq
c1
iB
cLessOrEq
c1
lA1
cLessOrEq,l3
20,1,562
x0
409647,cLessOrEq
cW
cGreater
cX
cGreater
cY
cGreater,eN2
539
x0
409615
eM2
xQ1
24832
eM2
c1
xR1
cGreater
c1
tM
cGreater
c1
n42
281856
eM2
c1
nR1
cGreater
c1
iB
cGreater
c1
lA1
cGreater,l3
20
iZ3
409647
eM2
cW
cGreaterOrEq
cX
cGreaterOrEq
cY
cGreaterOrEq,eN2
572
x0
46080
iN1
xQ1
24832
iN1
c1
xR1
cGreaterOrEq
c1
tM
cGreaterOrEq
c1
n42
281856
iN1
c1
nR1
cGreaterOrEq
c1
iB
cGreaterOrEq
c1
lA1
cGreaterOrEq,l3
20
iZ3
409641
iN1,l3
4,1,519,{1,137,cNot,tC3
571,{1,2,cNot,l2
0,1,452
x0
yY3
xE
0,2,537097,{3,547892744,cAnd,c8
16,1,566,{1,5,cAnd,tU
1}
}
,{ReplaceParams,16,1,569
x0
13314,xE
16,1,544
x0
553498,xE
16,1,546
x0
462369,xE
16,1,548
x0
466465,xE
0,1,457
x0
yY3
nG
570
x0
13314,nG
563
x0
8197,nG
541
x0
553498,nG
542
x0
462369,nG
543
x0
466465,nG
564
x0
143365,cOr,c8
4,1,525,{1,137,c03
tC3
572,{1,2,c03
l3
17,1,0,{1,0,c03
eN2
537,{1,256,cAbsNotNot,c8
18,1,531,{1,254,cAbsNotNot,c8
0,1,572,{3,43039744,c13
tB3
571,{3,49325056,c13
tC3
454,{3,32513586,c13
l2
16,3,32542225,{3,36732434,c13
yG1}
,}
;cZ2
grammar_optimize_abslogical_type{y9
9
cT
grammar_optimize_abslogical_type
grammar_optimize_abslogical={9,{34,192,228,238,242,247,254,260,261}
}
;}
cZ2
grammar_optimize_ignore_if_sideeffects_type{y9
59
cT
grammar_optimize_ignore_if_sideeffects_type
grammar_optimize_ignore_if_sideeffects={59,{0,20,21,22,23,24,25,26,cS
iT1
78,yJ1
cU
cZ2
grammar_optimize_nonshortcut_logical_evaluation_type{y9
56
cT
grammar_optimize_nonshortcut_logical_evaluation_type
grammar_optimize_nonshortcut_logical_evaluation={56,{0,25,cS
iT1
78,yJ1
163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,241,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259}
}
;}
cZ2
grammar_optimize_recreate_type{y9
22
cT
grammar_optimize_recreate_type
grammar_optimize_recreate={22,{18,55,56,57,80,81,82,83,84,85,117,118,120,121,130,131,132,133,134,135,136,137}
}
;}
cZ2
grammar_optimize_round1_type{y9
125
cT
grammar_optimize_round1_type
grammar_optimize_round1={125,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,25,cS
37,38,iT1
45,46,47,48,49,50,51,52,53,54,58,59,60,61,62,63,64,65,66,67,68,69,70,71,78,79,80,81,82,83,84,85,86,87,88,93,94,95,96,97,98,99,100,101,117,118,119,120,121,122,123,124,125,126,127,128,129,138,160,161,162,cU
cZ2
grammar_optimize_round2_type{y9
103
cT
grammar_optimize_round2_type
grammar_optimize_round2={103,{0,15,16,17,25,cS
39,40,iT1
45,46,47,48,49,50,51,52,53,54,59,60,72,73,78,79,86,87,88,89,90,91,92,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,119,122,123,124,125,126,127,128,139,159,160,161,162,cU
cZ2
grammar_optimize_round3_type{y9
79
cT
grammar_optimize_round3_type
grammar_optimize_round3={79,{74,75,76,77,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,170,171,172,173,174,175,176,177,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,197,198,199,201,202,203,205,206,207,208,209,210,211,213,214,215,217,218,219,220,221,222,223,225,226,227,229,230,231,232,233,234,235}
}
;}
cZ2
grammar_optimize_round4_type{y9
12
cT
grammar_optimize_round4_type
grammar_optimize_round4={12,{18,55,56,57,130,131,132,133,134,135,136,137}
}
;}
cZ2
grammar_optimize_shortcut_logical_evaluation_type{y9
53
cT
grammar_optimize_shortcut_logical_evaluation_type
grammar_optimize_shortcut_logical_evaluation={53,{0,25,cS
iT1
78,yJ1
cU}
l93
l31{nQ2
e22
eE1
x91
paramlist,lC1){index=(paramlist>>(index*10))&1023;if(index>=57
lD4
e22(SubFunction,cE2
plist_s[index-57]l64
index>=37
lD4
e22(NumConstant,cE2
plist_n_container
xH::plist_n[index-37])iI
e22(t23,cE2
plist_p[index]);}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <stdio.h>
#include <algorithm>
#include <map>
#include <sstream>
using
l93
FUNCTIONPARSERTYPES;using
l93
l31;using
t7;using
y01;l93{nI1
It,xJ3
T,xJ3
Comp>t71
MyEqualRange(It
first,It
last,const
T&val,Comp
comp){size_t
len=last-first;while(len>0){size_t
y23
len/2;It
xD3(first);xD3+=half;if(comp(*xD3,val)){first=xD3;++first;len=len-half-1;}
iH1
comp(val,*xD3)){len=half;}
else{It
left(first);{It&eO2=left;It
last2(xD3)i42
len2=last2-eO2;while(len2>0){size_t
half2=len2/2;It
eD3(eO2);eD3+=half2;if(comp(*eD3,val)){eO2=eD3;++eO2;len2=len2-half2-1;}
else
len2=half2;}
}
first+=len;It
right(++xD3);{It&eO2=right;It&last2=first
i42
len2=last2-eO2;while(len2>0){size_t
half2=len2/2;It
eD3(eO2);eD3+=half2;if(comp(val,*eD3))len2=half2;else{eO2=eD3;++eO2;len2=len2-half2-1;}
}
}
return
t71(left,right);}
}
return
t71(first,first);}
nQ2
cZ2
OpcodeRuleCompare{iL2()e42
tree,xH3
yC2
e93{const
Rule&rule=grammar_rules[yC2]iI
tG2<rule
x22.subfunc_opcode;}
iL2()x91
yC2,const
eX
const{const
Rule&rule=grammar_rules[yC2]iI
rule
x22.subfunc_opcode<tG2;}
}
y31
bool
TestRuleAndApplyIfMatch
cI1
t42
lM1&tree,bool
c9{MatchInfo
xH
info;n21
found(false,e2()l64(rule.lB1
LogicalContextOnly)&&!c9{tV1
if(nB
IsIntType
xH::xC3){if(rule.lB1
NotForIntegers)tV1
cK3
rule.lB1
OnlyForIntegers)tV1
if(nB
IsComplexType
xH::xC3){if(rule.lB1
NotForComplex)tV1
cK3
rule.lB1
OnlyForComplex)tV1
for(;;){
#ifdef DEBUG_SUBSTITUTIONS
#endif
found=TestParams(rule
x22
iK2,found.specs,info,true
l64
found.found)break;if(!&*found.specs){fail:;
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
iK2,info,false);
#endif
return
t43}
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
iK2,info,true);
#endif
SynthesizeRule(rule
iK2,info
nZ2}
y01{y41
ApplyGrammar
cI1
Grammar&y33,lM1&tree,bool
c9{if(tree.GetOptimizedUsing()==&y33){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Already optimized:  "
c72
tree)x71"\n"
<<std::flush;
#endif
return
t43
if(true){bool
changed
eV3
eL3
tG2
cZ3
cNot:case
cNotNot:case
cAnd:case
cOr:for
y11
tree.xC
true))yW1
lD
cIf:case
cAbsIf:if(ApplyGrammar(y33,xZ
0),tG2==cIf))yW1
eQ1
1;a<tree.xC
c9)yW1
break;yT3
for
y11
tree.xC
false))yW1}
if(changed){tree.Mark_Incompletely_Hashed(nZ2}
typedef
const
xH3
short*nI3;std::pair<nI3,nI3>range=MyEqualRange(y33.rule_list,y33.rule_list+y33.rule_count
iK2,OpcodeRuleCompare
xH());std
xV3<xH3
short>rules;rules.xR3
range
t03-range.first);for
xX
if(IsLogisticallyPlausibleParamsMatch(eW1
x22
iK2))rules.push_back(*r);}
range.first=!rules
tH3?&rules[0]:0;range
t03=!rules
tH3?&rules[rules.size()-1]+1:0;if(range.first!=range
t03){
#ifdef DEBUG_SUBSTITUTIONS
if(range.first!=range
t03){std::cout<<"Input ("
<<t13
tG2)<<")["
<<tree.GetParamCount()<<"]"
;if(c9
std::cout<<"(Logical)"
;xH3
first=iU1,prev=iU1;t52
sep=", rules "
;for
xX
if(first==iU1)first=prev=*r;iH1*r==prev+1)prev=*r;else{std::cout<<sep<<first;sep=","
;if(prev!=first)std::cout<<'-'<<prev;first=prev=*r;}
}
if(first!=iU1){std::cout<<sep<<first;if(prev!=first)std::cout<<'-'<<prev;}
std::cout<<": "
c72
tree)x71"\n"
<<std::flush;}
#endif
bool
changed
eV3
for
xX
#ifndef DEBUG_SUBSTITUTIONS
if(!IsLogisticallyPlausibleParamsMatch(eW1
x22
iK2))lI3
#endif
if(TestRuleAndApplyIfMatch(eW1
iK2,c9){yW1
y73}
if(changed){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Changed."
<<std::endl
x71"Output: "
c72
tree)x71"\n"
<<std::flush;
#endif
tree.Mark_Incompletely_Hashed(nZ2}
tree.SetOptimizedUsing(&y33)iI
t43
y41
ApplyGrammar
cI1
void*p,FPoptimizer_CodeTree::eX
yR
ApplyGrammar(*cI1
Grammar*)p
iK2);}
eV1
ApplyGrammars(FPoptimizer_CodeTree::eX{
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round1\n"
;
#endif
n7
grammar_optimize_round1
iK2
x2
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round2\n"
;
#endif
n7
grammar_optimize_round2
iK2
x2
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round3\n"
;
#endif
n7
grammar_optimize_round3
iK2
x2
#ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_nonshortcut_logical_evaluation\n"
;
#endif
n7
grammar_optimize_nonshortcut_logical_evaluation
iK2
x2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_round4\n"
;
#endif
n7
grammar_optimize_round4
iK2
x2
#ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_shortcut_logical_evaluation\n"
;
#endif
n7
grammar_optimize_shortcut_logical_evaluation
iK2
x2
#endif
#ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_ignore_if_sideeffects\n"
;
#endif
n7
grammar_optimize_ignore_if_sideeffects
iK2
x2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iN3"grammar_optimize_abslogical\n"
;
#endif
n7
grammar_optimize_abslogical
iK2
x2
#undef C
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
#include <cstring>
#include <cmath>
#include <memory> /* for auto_ptr */
using
l93
FUNCTIONPARSERTYPES;using
l93
l31;using
t7;using
y01;l93{y41
TestImmedConstraints
x91
bitmask,const
eX{eL3
bitmask&ValueMask
cZ3
Value_AnyNum:case
ValueMask:lD
xA2:if(GetEvennessInfo(tree)!=n52
Value_OddInt:if(GetEvennessInfo(tree)!=xC2
i41:if(GetIntegerInfo(tree)!=n52
Value_NonInteger:if(GetIntegerInfo(tree)!=xC2
eY1:if(!IsLogicalValue(tree)cV
nB1
SignMask
cZ3
Sign_AnySign:lD
nA1:if(l11
n52
eZ1:if(l11
xC2
Sign_NoIdea:if(l11
Unknown
cV
nB1
OnenessMask
cZ3
Oneness_Any:case
OnenessMask:lD
Oneness_One:if(!xS1
if(!n23
fp_abs(l54)iM2(1))cV
lD
Oneness_NotOne:if(!xS1
lL3
fp_abs(l54)iM2(1))cV
nB1
ConstnessMask
cZ3
Constness_Any:lD
i31:if(!xS1
lD
Constness_NotConst:if(xS1
y73
l12}
tT1
xH3
extent,xH3
nbits,xJ3
eP2=xH3
int>cZ2
nbitmap{private:static
const
xH3
bits_in_char=8;static
const
xH3
eQ2=(c23
eP2)*bits_in_char)/nbits;eP2
data[(extent+eQ2-1)/eQ2];eO3
void
inc(lC1,int
by=1){data[pos(index)]+=by*eP2(1<<yD2);t11
void
dec(lC1){inc(index,-1);}
int
get(lC1
c31(data[pos(index)]>>yD2)&mask()e03
pos(lC1)yR
index/eQ2
e03
shift(lC1)yR
nbits*(index%eQ2)e03
mask()yR(1<<nbits)-1
e03
mask(lC1)yR
mask()<<yD2;}
}
;cZ2
cV3{int
SubTrees:8;int
Others:8;int
lI2:8;int
eM3:8;nbitmap<l13,2>SubTreesDetail;cV3(){std::memset(this,0,c23*this));}
cV3
cI1
cV3&b){std::memcpy(this,&b,c23
b));}
cV3&eI1=cI1
cV3&b){std::memcpy(this,&b,c23
b))iI*this;}
}
y31
cV3
CreateNeedList_uncached(tD&c02){cV3
cJ1;for
x91
a=0;a<c02
yE2
tO3
const
e22&nJ3=eE1
xH(c02.param_list,a);e53
SubFunction:{cR&param
yC3;yJ
GroupFunction)++cJ1.eM3;else{++nC3;assert(param.data.subfunc_opcode<VarBegin);cJ1.SubTreesDetail.inc(param
iT);}
++cJ1.lI2;y73
case
e63
case
t23:++nB3;++cJ1.lI2;y73}
return
cJ1;}
nQ2
cV3&CreateNeedList(tD&c02){typedef
std::map<tD*,cV3>eA1;static
eA1
yY1;eA1::yK3
i=yY1.xT2&c02
l64
i!=yY1.cY1&c02
lD4
i
e52
iI
yY1.yG3,std::make_pair(&c02,CreateNeedList_uncached
xH(c02)))e52;}
nQ2
lM1
CalculateGroupFunction
cI1
e22&yF2
const
tC
info){e53
e63{const
ParamSpec_NumConstant
xH
yV3
cI1
ParamSpec_NumConstant
xH*lH2
iI
CodeTreeImmed
e43
constvalue);}
case
t23:{cQ
yV3(cQ*lH2
iI
tK3
GetParamHolderValueIfFound
e43
index);}
case
SubFunction:{cR&param
yC3
eO
xC3;xC3
eN
param
iT);iC2
l52).xR3
eG3
data
yE2);for
x91
a=0;a<eG3
data
yE2
tO3
lM1
tmp(CalculateGroupFunction(eE1
xH
e43
data.param_list,a),info));xC3
c0
tmp);}
xC3
tB1)iI
xC3;}
}
return
lM1();}
}
y01{y41
IsLogisticallyPlausibleParamsMatch(tD&c02,const
eX{cV3
cJ1(CreateNeedList
xH(c02))i42
tD3=c4
if(tD3<size_t(cJ1.lI2))tU2
for
y11
tD3
tO3
xH3
opcode=lA4
nC;eL3
opcode
cZ3
cImmed:if(cJ1.eM3>0)--cJ1.eM3;else--nB3;lD
l13:case
cFCall:case
cPCall:--nB3;break;yT3
assert(opcode<VarBegin);if(nC3>0&&cJ1.SubTreesDetail.get(opcode)>0){--nC3;cJ1.SubTreesDetail.dec(opcode);}
else--nB3;}
}
if(cJ1.eM3>0||nC3>0||nB3>0)tU2
if(c02.match_type!=AnyParams){if(0||nC3<0||nB3<0)tU2}
l12}
nQ2
n21
TestParam
cI1
e22&yF2
l92&tree
t32
start_at
iE3){e53
e63{const
ParamSpec_NumConstant
xH
yV3
cI1
ParamSpec_NumConstant
xH*lH2;if(!xS1
nW3
imm=l54;switch
e43
modulo
cZ3
Modulo_None:lD
Modulo_Radians:imm=fp_mod(imm,y5
imm<x61
imm
yP
if(imm>fp_const_pi
xH())imm-=fp_const_twopi
xH(c52
return
n23
imm,eG3
constvalue);}
case
t23:{cQ
yV3(cQ*lH2;if(!x7
return
tK3
SaveOrTestParamHolder
e43
index
iK2);}
case
SubFunction:{cR&param
yC3;yJ
GroupFunction){if(!x7
lM1
xT1=CalculateGroupFunction(yF2
info);
#ifdef DEBUG_SUBSTITUTIONS
DumpHashes(xT1)x71*cI1
void**)&xT1.xJ1
x71"\n"
x71*cI1
void**)&l54
x71"\n"
;DumpHashes(tree)x71"Comparing "
c72
xT1)x71" and "
c72
tree)x71": "
x71(xT1
xL
tree)?"true"
:"false"
)x71"\n"
;
#endif
return
xT1
xL
tree);}
cK3!&*start_at){if(!x7
if(tG2!=param
iT
cV}
return
TestParams
e43
data
iK2,start_at,info,false);}
}
}
return
t43
nQ2
cZ2
iU
xD2
MatchInfo
xH
info;iU()e73,info(){}
}
y31
class
MatchPositionSpec_PositionalParams:e51
iU
xH>{eO3
l33
MatchPositionSpec_PositionalParams
eZ2
n):eB1
iU
xH>(n){}
}
;cZ2
iV1
xD2
iV1()e73{}
}
;class
yK:e51
iV1>{eO3
xH3
trypos;l33
yK
eZ2
n):eB1
iV1>(n),trypos(0){}
}
y31
n21
TestParam_AnyWhere
cI1
e22&yF2
l92&tree
t32
start_at
iE3,e13&used,bool
xF2{xT<yK>xB;xH3
a;if(&lA2
yK
xE2
a=xB->trypos;goto
retry_anywhere_2;c82
yK
iJ);a=0;}
tE3
c4++a){if(used[a])lI3
retry_anywhere
l34=TestParam(yF2
lA4,tG3);tL3
used[a]=true;if(xF2
tK3
SaveMatchedParamIndex(a);xB->trypos=a
iI
n21(true,&*xB);}
}
retry_anywhere_2
l44
goto
retry_anywhere;}
}
return
t43
nQ2
cZ2
yK1
xD2
MatchInfo
xH
info;e13
used;l33
yK1
eZ2
tD3)e73,info(),used(tD3){}
}
y31
class
MatchPositionSpec_AnyParams:e51
yK1
xH>{eO3
l33
MatchPositionSpec_AnyParams
eZ2
n,size_t
m):eB1
yK1
xH>(n,yK1
xH(m)){}
}
y31
n21
TestParams(tD&nQ,l92&tree
t32
start_at
iE3,bool
xF2{if(nQ.match_type!=AnyParams){if(y1!=tree.GetParamCount()cV}
if(!IsLogisticallyPlausibleParamsMatch(nQ
iK2))tU2
eL3
nQ.match_type
cZ3
PositionalParams:{xT<cJ>xB;xH3
a;if(&lA2
cJ
xE2
a=y1-1;goto
lD1;c82
cJ(y1);a=0;}
tE3
y1
tO3(*xB)[a].tW2
retry_positionalparams
l34=TestParam(cZ
a),lA4,tG3);tL3
t02}
lD1
l44
info
tF3.info;goto
retry_positionalparams;}
tN3--a;goto
lD1;}
info=(*xB)[0].info
iI
t43
if(xF2
for
x91
a=0;a<y1;++a)tK3
SaveMatchedParamIndex(a)iI
n21(true,&*xB);}
case
SelectedParams:case
AnyParams:{xT<t8>xB;e13
used
iJ);std
xV3<xH3>l43(y1);std
xV3<xH3>yG2(y1)x81{const
e22
nJ3=cZ
a);l43[a]=ParamSpec_GetDepCode(nJ3);}
{xH3
b=0
x81
if(l43[a]!=0)yG2[b++]=a
x81
if(l43[a]==0)yG2[b++]=a;}
xH3
a;if(&lA2
t8
xE2
if(y1==0){a=0;goto
retry_anyparams_4;}
a=y1-1;goto
eC1;c82
t8(y1
iK2.GetParamCount());a=0;if(y1!=0){(*xB)[0].tW2(*xB)[0].used=used;}
}
tE3
y1
tO3
tN3(*xB)[a].tW2(*xB)[a].used=used;}
retry_anyparams
l34=TestParam_AnyWhere
xH(cZ
yG2[a])iK2,tG3,used,xF2;tL3
t02}
eC1
l44
info
tF3.info;used
tF3.used;goto
retry_anyparams;}
eD1:tN3--a;goto
eC1;}
info=(*xB)[0].info
iI
t43
retry_anyparams_4:if(nQ.n2!=0){if(!TopLevel||!tK3
HasRestHolder(nQ.n2)){eM
cF2;cF2.reserve
iJ);for
x91
b=0;b<c4++b){if(c33)lI3
cF2.push_back(xZ
b));c33=true;if(xF2
tK3
SaveMatchedParamIndex(b);}
if(!tK3
SaveOrTestRestHolder(nQ.n2,cF2)){goto
eD1;}
}
else{iA3&cF2=tK3
GetRestHolderValues(nQ.n2)xZ1
cF2.size()tO3
bool
found
eV3
for
x91
b=0;b<c4++b){if(c33)lI3
if(cF2[a]xL
xZ
b))){c33=true;if(xF2
tK3
SaveMatchedParamIndex(b);found=true;y73}
if(!found){goto
eD1;}
}
}
}
return
n21(true,y1?&*xB:0);}
case
GroupFunction:y73
return
t43}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
using
t7;using
y01;l93{nQ2
lM1
xU1
const
e22&yF2
tC
info,bool
inner=true){e53
e63{const
ParamSpec_NumConstant
xH
yV3
cI1
ParamSpec_NumConstant
xH*lH2
iI
CodeTreeImmed
e43
constvalue);}
case
t23:{cQ
yV3(cQ*lH2
iI
tK3
GetParamHolderValue
e43
index);}
case
SubFunction:{cR&param
yC3
eO
tree;tH2
param
iT);for
x91
a=0;a<eG3
data
yE2
tO3
lM1
nparam=xU1
eE1
xH
e43
data.param_list,a),info,true)xY2
nparam);}
if
e43
data.n2!=0){eM
trees(tK3
GetRestHolderValues
e43
data.n2));tree.AddParamsMove(trees);if
iJ==1){assert(tree.GetOpcode()==cAdd iM3()==cMul iM3()==cMin iM3()==cMax iM3()==cAnd iM3()==cOr iM3()==cAbsAnd iM3()==cAbsOr);tree
eA2
0));}
else
if
iJ==0){eL3
tG2
cZ3
cAdd:case
cOr:tree=n71
0));lD
cMul:case
cAnd:tree=n71
1));yT3
y73}
}
if(inner)tree
tB1)iI
tree;}
}
return
lM1();}
}
y01{eV1
SynthesizeRule
cI1
t42
lM1&tree
iE3){eL3
rule.ruletype
cZ3
ProduceNewTree:{tree.Become(xU1
eE1
lP1
0),info,false)c52
case
ReplaceParams:yT3{std
xV3<xH3>list=tK3
GetMatchedParamIndexes();std::sort(list.iQ2
list.end());eQ1
list.size();a-->0;)tree.l83
list[a]);for
x91
a=0;a<rule.repl_param_count
tO3
lM1
nparam=xU1
eE1
lP1
a),info,true)xY2
nparam);}
y73}
}
}
#endif
#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
#include <cstring>
using
l93
FUNCTIONPARSERTYPES;using
l93
l31;using
t7;using
y01;l93
l31{eV1
DumpMatch
cI1
t42
l92&tree,const
tC
info,bool
DidMatch,std::ostream&o){DumpMatch(rule
iK2,info,DidMatch?l14"match"
:l14"mismatch"
,o);}
eV1
DumpMatch
cI1
t42
l92&tree,const
tC
info,t52
tQ3,std::ostream&o){static
const
char
ParamHolderNames[][2]={"%"
,"&"
,"x"
,"y"
,"z"
,"a"
,"b"
,"c"
}
;o<<tQ3<<" (rule "
<<(&rule-grammar_rules)<<")"
<<":\n  Pattern    : "
;{e22
tmp;tmp.first=SubFunction;ParamSpec_SubFunction
tmp2;tmp2.data=rule
x22;tmp
t03=cE2
tmp2;DumpParam
xH(tmp,o);}
o<<"\n  Replacement: "
;DumpParams
lP1
rule.repl_param_count,o);o<<"\n"
;o<<"  Tree       : "
c72
tree,o);o<<"\n"
;if(!std::strcmp(tQ3,l14"match"
))DumpHashes(tree,o)xZ1
tK3
cH
size()tO3
if(!tK3
paramholder_matches[a]x32))lI3
o<<"           "
<<ParamHolderNames[a]<<" = "
c72
tK3
paramholder_matches[a],o);o<<"\n"
;}
e31
tK3
lR.size();++b){if(!tK3
lR[b
eU3)continue
xZ1
tK3
lR[b
eZ3.size()tO3
o<<"         <"
<<b<<"> = "
c72
tK3
lR[b
eZ3[a],o);o<<std::endl;}
}
o<<std::flush;}
}
#endif
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;l93{y41
MarkIncompletes(FPoptimizer_CodeTree::eX{if(tree.Is_Incompletely_Hashed())l12
bool
iW1=false
xZ1
c4++a)iW1|=MarkIncompletes(lA4
l64
iW1)tree.Mark_Incompletely_Hashed()iI
iW1;}
eV1
FixIncompletes(FPoptimizer_CodeTree::eX{if(tree.Is_Incompletely_Hashed()){for
y11
c4++a)FixIncompletes(lA4);tree
tB1);}
}
}
t7{lC
Sort()eN3
Sort();}
lC
Rehash(bool
constantfolding){if(constantfolding)ConstantFolding(*this);else
Sort();data
xD
nQ2
cZ2
cE{c43
nW3
c53
yL1=0;
#if 0
long
double
value=Value;e9=crc32::calc(cI1
xH3
char*)&value,c23
value));key^=(key<<24);
#elif 0
union{cZ2{xH3
char
filler1[16]tV2
v;xH3
char
filler2[16];}
buf2;cZ2{xH3
char
filler3[c23
nW3)+16-c23
x31)];e9;}
buf1;}
data;memset(&data,0,c23
data));data.buf2.v=Value;e9=data.buf1.key;
#else
int
xG3
tV2
nV2=std::frexp(Value,&xG3);e9=x91
e92+0x8000)&0xFFFF
l64
nV2<0){nV2=-nV2;key=key^0xFFFF;}
else
key+=0x10000;nV2-=yJ2;key<<=39;key|=n31(nV2+nV2)*nW3(1u<<31))<<8;
#endif
lQ
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nI1
T
cG2
std::complex<T> >{c43
std::complex<T>c53
cE<T>::nK3
cH2,Value.real());nB
fphash_t
temp;cE<T>::nK3
temp,Value.imag());yL1^=temp.hash2;cH2.hash2^=temp.hash1;}
}
;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
tT1
cG2
long>{yL
long
Value){e9=Value;lQ
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
tT1
cG2
GmpInt>{c43
GmpInt
c53
e9=Value.toInt();lQ
#endif
eV1
xK2
xH::Recalculate_Hash_NoRecursion(){fphash_t
cH2(n31
Opcode)<<56,Opcode*iQ3(0x1131462E270012B));Depth=1;eL3
Opcode
cZ3
cImmed:{cE
xH::nK3
cH2,Value
c52
case
l13:{yL1|=n31
cQ1<<48
lQ1((n31
cQ1)*11)^iQ3(0x3A83A83A83A83A0);y73
case
cFCall:case
cPCall:{yL1|=n31
cQ1<<48
lQ1((~n31
cQ1)*7)^3456789;}
yT3{size_t
t81=0
xZ1
y93
tO3
if(l03
y22>t81)t81=l03
y22;yL1+=((l03
tX2
hash1*(a+1))>>12)lQ1
l03
tX2
hash1
lQ1(3)*iQ3(0x9ABCD801357);cH2.hash2*=iQ3(0xECADB912345)lQ1(~l03
tX2
hash2)^4567890;}
Depth+=t81;}
}
if(Hash!=cH2){Hash=cH2;l82=0;}
}
lC
FixIncompleteHashes(){MarkIncompletes(*this);FixIncompletes(*this);}
}
#endif
#include <cmath>
#include <list>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;l93{using
t7
y31
bool
nT1
l92&tree,long
count,const
xC1::SequenceOpCode
xH&eW,xC1
y83&synth,size_t
max_bytecode_grow_length);static
const
cZ2
SinCosTanDataType{OPCODE
whichopcode;OPCODE
inverse_opcode;enum{eR2,nW2,nS1,lE1}
;OPCODE
codes[4];}
SinCosTanData[12]={{cTan,cCot,{cSin,cCos,cCsc,cSec}
}
,{cCot,cCot,{cCos,cSin,cSec,cCsc}
}
,{cCos,cSec,{cSin,cTan,cCsc,cCot}
}
,{cSec,cCos,{cTan,cSin,cCot,cCsc}
}
,{cSin,cCsc,{cCos,cCot,cSec,cTan}
}
,{cCsc,cSin,{cCot,cCos,cTan,cSec}
}
,{cI2{cSinh,cCosh,cJ2,{cSinh,cNop,{cI2
cNop,cCosh}
}
,{cCosh,cNop,{cSinh,cI2
cNop}
}
,{cNop,cTanh,{cCosh,cSinh,cJ2,{cNop,cSinh,{cNop,cTanh,cCosh,cNop}
}
,{cNop,cCosh,{cTanh,cSinh,cJ2}
;}
t7{lC
SynthesizeByteCode(std
xV3<xH3>&ByteCode,std
xV3
xH&Immed,size_t&stacktop_max){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Making bytecode for:\n"
;iP
#endif
while(RecreateInversionsAndNegations()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"One change issued, produced:\n"
;iP
#endif
FixIncompleteHashes();using
y01;using
l93
l31;const
void*g=cE2
grammar_optimize_recreate;while(ApplyGrammar(*cI1
Grammar*)g,*this)){FixIncompleteHashes();}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Actually synthesizing, after recreating inv/neg:\n"
;iP
#endif
xC1
y83
synth;SynthesizeByteCode(synth,false);synth.Pull(ByteCode,Immed,stacktop_max);}
lC
SynthesizeByteCode(xC1
y83&synth,bool
MustPopTemps
e93{xV1*this))yR;}
for
y11
12
tO3
const
SinCosTanDataType&data=SinCosTanData[a];if(data.whichopcode!=cNop)lB2!=data.whichopcode)lI3
lO1
nL3;nL3.lK1
nL3
eN
data.inverse_opcode);nL3
tY2;xV1
nL3)){synth.iX
else
lB2!=cInv)lI3
if(GetParam(0)nC!=data.inverse_opcode)lI3
xV1
GetParam(0))){synth.iX
size_t
found[4];e31
4;++b){lO1
tree;if(data.tR3]==cNop){tH2
cInv);lO1
nM3;nM3.lK1
nM3
eN
data.tR3^2]);nM3
tY2
xY2
nM3);}
else{tree.lK1
tH2
data.tR3]);}
tree
tY2;found[b]tS3
xU3
tree);}
if(found[data.eR2
e4
nW2]!=iV
eR2]cK2
nW2
iW
cDiv
iY
eR2
e4
lE1]!=iV
eR2]cK2
lE1
iW
cMul
iY
nS1
e4
lE1]!=iV
nS1]cK2
lE1
iW
cRDiv
iY
nS1
e4
nW2]!=iV
nS1]cK2
nW2
iW
cMul,2,1);synth.iX
size_t
n_subexpressions_synthesized=SynthCommonSubExpressions(synth);eL3
lC2{case
l13:synth.PushVar(GetVar());lD
cImmed:t91
xJ1);lD
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:lB2==cMul){bool
c63
eV3
yM
lW1
cV2)&&isLongInteger(lW1.xJ1)){e71=makeLongInteger(lW1.xJ1);lO1
tmp(*this,xJ3
lO1::CloneTag());tmp
cK1
tmp
tB1
l64
nT1
tmp,value,xC1::tY1
xH::AddSequence,synth,MAX_MULI_BYTECODE_LENGTH)){c63=true;y73}
}
if(c63)y73
int
yM1=0;e13
done(GetParamCount(),false);lO1
iC;iC
eN
lC2;for(;;){bool
found
eV3
yM
done[a])lI3
if(synth.IsStackTop(lW1)){found=true;done[a]=true;lW1.n8
iC
eQ
lW1
l64++yM1>1){yT
2);iC
tY2
yO1
iC);yM1=yM1-2+1;}
}
}
if(!found)y73
yM
done[a])lI3
lW1.n8
iC
eQ
lW1
l64++yM1>1){yT
2);iC
tY2
yO1
iC);yM1=yM1-2+1;}
}
if(yM1==0){eL3
lC2{case
cAdd:case
cOr:case
cAbsOr:t91
0);lD
cMul:case
cAnd:case
cAbsAnd:t91
1);lD
cMin:case
cMax:t91
0);break;yT3
y73++yM1;}
assert(n_stacked==1);y73
case
cPow:{t22
p0
iE2
0);t22
p1
iE2
1
l64!p1.c42!isLongInteger
eE3)||!nT1
p0,makeLongInteger
eE3),xC1::tY1
xH::MulSequence,synth,MAX_POWI_BYTECODE_LENGTH)){p0.n8
p1.n8
yT
2);y21
cIf:case
cAbsIf:{xJ3
xC1
y83::IfData
ifdata
tJ3(0)i93
SynthIfStep1(ifdata,lC2
tJ3(1)i93
SynthIfStep2(ifdata)tJ3(2)i93
SynthIfStep3(ifdata
c52
case
cFCall:case
cPCall:{for
y11
l21++a)lW1.n8
yT
x91)GetParamCount());lE3
nR2|GetFuncNo(),0,0
c52
yT3{for
y11
l21++a)lW1.n8
yT
x91)GetParamCount()c52}
synth.StackTopIs(*this
l64
MustPopTemps&&n_subexpressions_synthesized>0){size_t
top
tS3
GetStackTop();synth.DoPopNMov(top-1-n_subexpressions_synthesized,top-1);}
}
}
l93{y41
nT1
l92&tree,long
count,const
xC1::SequenceOpCode
xH&eW,xC1
y83&synth,size_t
max_bytecode_grow_length){if
cX3!=0){xC1
y83
backup=synth;tree.n8
size_t
bytecodesize_backup
tS3
GetByteCodeSize();xC1::nT1
count
nU2
size_t
bytecode_grow_amount
tS3
GetByteCodeSize()-bytecodesize_backup;if(bytecode_grow_amount>max_bytecode_grow_length){synth=backup
iI
t43
l12}
else{xC1::nT1
count,eW,synth
nZ2}
}
#endif
#include <cmath>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;l93{using
t7;
#define FactorStack std xV3
const
cZ2
PowiMuliType{xH3
opcode_square;xH3
opcode_cumulate;xH3
opcode_invert;xH3
opcode_half;xH3
opcode_invhalf;}
iseq_powi={cSqr,cMul,cInv,cSqrt,cRSqrt}
,iseq_muli={iU1
c2
cNeg,iU1,iU1}
y31
nW3
tA1
cI1
PowiMuliType&c73,const
std
xV3
tD1,n62&stack
i22
1);while(IP<limit){iB3
c73.opcode_square){if(!tD2
i32
2;yY
opcode_invert){xC3=-xC3;yY
opcode_half){if(xC3>y71&&isEvenInteger(i32
yJ2;yY
opcode_invhalf){if(xC3>y71&&isEvenInteger(i32
nW3(-0.5);++IP;t02
size_t
xG2=IP
tV2
lhs(1);iB3
cFetch){lC1=nA2;if(index<y4||size_t(index-y4)>=iZ1){IP=xG2;y73
lhs=stack[index-y4];goto
yH2;}
iB3
cDup){lhs=xC3;goto
yH2;yH2:tC1
xC3);++IP
tV2
subexponent=tA1(c73
tN
if(IP>=limit||iX1!=c73.opcode_cumulate){IP=xG2;y73++IP;stack.pop_back();xC3+=lhs*subexponent;t02
y73
return
xC3;}
nQ2
nW3
ParsePowiSequence
cI1
std
xV3
tD1){n62
stack;tC1
nW3(1))iI
tA1(iseq_powi
tN}
nQ2
nW3
ParseMuliSequence
cI1
std
xV3
tD1){n62
stack;tC1
nW3(1))iI
tA1(iseq_muli
tN}
nQ2
class
CodeTreeParserData{eO3
l33
CodeTreeParserData(bool
k_powi):stack(),clones(),keep_powi(k_powi){}
void
Eat
eZ2
tD3,OPCODE
opcode){lM1
xQ;xQ
eN
opcode);eM
c02=Pop(tD3
cF1
c02
l64!keep_powi)eL3
opcode
cZ3
cTanh:{lM1
sinh,cosh;sinh
eN
cSinh);sinh
eQ
xQ
xX2
sinh
tB1);cosh
eN
cCosh);cosh
c0
xQ
xX2
cosh
iR2
pow;pow
eN
cPow);pow
c0
cosh);pow.yF
nW3(-1)));pow
tB1);xQ
l84);xQ.n81
0,sinh);xQ
c0
pow
c52
case
cTan:{lM1
sin,cos;sin
eN
cSin);sin
eQ
xQ
xX2
sin
tB1);cos
eN
cCos);cos
c0
xQ
xX2
cos
iR2
pow;pow
eN
cPow);pow
c0
cos);pow.yF
nW3(-1)));pow
tB1);xQ
l84);xQ.n81
0,sin);xQ
c0
pow
c52
case
cPow:{l92&p0=xQ
l8
0);l92&p1=xQ
l8
1
l64
p1
nC==cAdd){eM
xI3(p1.GetParamCount())xZ1
p1.l21++a){lM1
pow;pow
eN
cPow);pow
eQ
p0);pow
eQ
p1
eJ3;pow
tB1);xI3[a
tT3
pow);}
xQ
l84
cF1
xI3);}
y73
yT3
y73
xQ
tB1!keep_powi);iY1,false);
#ifdef DEBUG_SUBSTITUTIONS
tH1<<tD3<<", "
<<t13
opcode)<<"->"
<<t13
xQ
nC)<<": "
iV3
xQ
tV
xQ);
#endif
tC1
e23
EatFunc
eZ2
tD3,OPCODE
eQ3
xH3
funcno
xX3
CodeTreeFuncOp
xH(eQ3
funcno);eM
c02=Pop(tD3
cF1
c02);xQ
tY2;
#ifdef DEBUG_SUBSTITUTIONS
tH1<<tD3<<", "
iV3
xQ
tV
xQ);
#endif
iY1);tC1
e23
AddConst(yT1
xX3
CodeTreeImmed(value);iY1);Push(e23
AddVar
x91
varno
xX3
CodeTreeVar
xH(varno);iY1);Push(e23
SwapLastTwoInStack(){tE1
1
tT3
tE1
2])e33
Dup(){Fetch(iZ1-1)e33
Fetch
eZ2
which){Push(stack[which]);}
nI1
T>void
Push(T
tree){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<iV3
tree
tV
tree);
#endif
tC1
tree)e33
PopNMov
eZ2
target,size_t
source){stack[target]=stack[source];stack.eP3
target+1);}
lM1
yI2{clones.clear()eO
xC3(stack.back());stack.eP3
iZ1-1)iI
xC3;}
eM
Pop
eZ2
n_pop){eM
xC3(n_pop);for
x91
n=0;n<n_pop;++n)xC3[n
tT3
tE1
n_pop+n]);
#ifdef DEBUG_SUBSTITUTIONS
eY2
n=n_pop;n-->0;){tH1
c72
xC3[n]tV
xC3[n]);}
#endif
stack.eP3
iZ1-n_pop)iI
xC3;}
size_t
GetStackTop(c31
iZ1;}
private:void
FindClone(lM1&,bool=true)yR;}
private:eM
stack;std::multimap<fphash_t,lM1>clones;bool
keep_powi;private:CodeTreeParserData
cI1
CodeTreeParserData&);CodeTreeParserData&eI1=cI1
CodeTreeParserData&);}
y31
cZ2
IfInfo{lM1
eS2
eO
thenbranch
i42
endif_location;IfInfo():eS2(),thenbranch(),endif_location(){}
}
;}
t7{lC
GenerateFrom
cI1
xJ3
FunctionParserBase
xH::Data&xE3,bool
keep_powi){eM
xJ2;xJ2.xR3
xE3.mVariablesAmount);for
x91
n=0;n<xE3.mVariablesAmount;++n){xJ2.push_back(CodeTreeVar
xH(n+l13));}
GenerateFrom(xE3,xJ2,keep_powi);}
lC
GenerateFrom
cI1
xJ3
FunctionParserBase
xH::Data&xE3,const
x9&xJ2,bool
keep_powi){const
std
xV3<xH3>&ByteCode=xE3.mByteCode;const
std
xV3
xH&Immed=xE3.mImmed;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"ENTERS GenerateFrom()\n"
;
#endif
CodeTreeParserData
xH
sim(keep_powi);std
xV3<IfInfo
xH>eL;eY2
IP=0,DP=0;;++IP){i52:while(!eL
tH3&&(eL.eA==IP||(IP<i4
size()&&iX1==cJump&&eL.eF1
x32)))){lO1
elsebranch=sim.yI2
n82
eL.back().eS2)n82
eL.eF1)n82
elsebranch
yO
3,cIf);eL.pop_back();}
if(IP>=i4
size())break;xH3
opcode=iX1;if((opcode==cSqr||opcode==cDup||(opcode==cInv&&!IsIntType
xH::xC3)||opcode==cNeg||opcode==cSqrt||opcode==cRSqrt||opcode==cFetch)){size_t
was_ip=IP
tV2
e62
ParsePowiSequence
xH(ByteCode,IP,eL
tH3?i4
size():eL.eA,sim.xI
1);if
e92!=nW3(1.0)){xJ
xG3
yO
2,cPow);goto
i52;}
if(opcode==cDup||opcode==cFetch||opcode==cNeg
i12
xS2=ParseMuliSequence
xH(ByteCode,IP,eL
tH3?i4
size():eL.eA,sim.xI
1
l64
xS2!=nW3(1.0)){xJ
xS2
yO
2,cMul);goto
i52;}
}
IP=was_ip;}
if(n92>=l13){lC1=opcode-l13
n82
xJ2[index]);}
else{eL3
n92
cZ3
cIf:case
cAbsIf:{eL.eP3
eL.size()+1);lO1
res(sim.yI2);eL.back().eS2.swap(res);eL.eA=i4
size();IP+=2;t02
case
cJump:{lO1
res(sim.yI2);eL.eF1.swap(res);eL.eA=l73
IP+1]+1;IP+=2;t02
case
cImmed:xJ
Immed[DP++]);lD
cDup:sim.Dup();lD
cNop:lD
cFCall:{xH3
funcno=nA2;assert(funcno<fpdata.mFuncPtrs.size());xH3
c02=xE3.mFuncPtrs[funcno].mParams;sim.EatFunc(c02,n92,funcno
c52
case
cPCall:{xH3
funcno=nA2;assert(funcno<fpdata.iT3.size());const
FunctionParserBase
xH&p=*xE3.iT3[funcno].mParserPtr;xH3
c02=xE3.iT3[funcno].mParams;x9
paramlist=sim.Pop(c02);lO1
i62;i62.GenerateFrom(*p.mData,paramlist)n82
i62
c52
case
cInv:xJ
1)tI3
cDiv);lD
cNeg
nB2
cNeg);break;xJ
0)tI3
cSub);lD
cSqr:xJ
2
yN1
cSqrt:xJ
yJ2
yN1
cRSqrt:xJ
nW3(-0.5)yN1
cCbrt:xJ
nW3(1)/nW3(3)yN1
cDeg:xJ
fp_const_rad_to_deg
xH()yP1
cRad:xJ
fp_const_deg_to_rad
xH()yP1
cExp:cR1)goto
nO3;xJ
fp_const_e
xH())tI3
cPow);lD
cExp2:cR1)goto
nO3;xJ
2.0)tI3
cPow);lD
cCot
nB2
cTan
x1
cCsc
nB2
cSin
x1
cSec
nB2
cCos
x1
cInt:
#ifndef __x86_64
cR1
nD3
1,cInt
c52
#endif
xJ
yJ2
nN3
yO
1,cFloor);lD
cLog10
nB2
c83
fp_const_log10inv
xH()yP1
cLog2
nB2
c83
fp_const_log2inv
xH()yP1
cLog2by:nP3
c83
fp_const_log2inv
xH()yO
3,cMul);lD
cHypot:xJ
2
yO
2,cPow);tW3
xJ
2
yO
2,cPow
nN3);xJ
yJ2
yN1
cSinCos:sim.Dup(yO
1,cSin);nP3
cCos);lD
cSinhCosh:sim.Dup(yO
1,cSinh);nP3
cCosh);lD
cRSub:tW3
case
cSub:cR1
nD3
2,cSub
c52
xJ-1
yO
2,cMul
nN3);lD
cRDiv:tW3
case
cDiv:cR1||IsIntType
xH::xC3
nD3
2,cDiv
c52
xJ-1
yO
2,cPow
yP1
cAdd:case
cMul:case
cMod:case
cPow:case
cEqual:case
cLess:case
cGreater:case
i61:case
cLessOrEq:case
cGreaterOrEq:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:sim.Eat(2,tK1
lD
cNot:case
cNotNot:case
cB3:case
cAbsNotNot
nB2
tK1
lD
cFetch:sim.Fetch(nA2);lD
cPopNMov:{xH3
stackOffs_target=nA2;xH3
stackOffs_source=nA2;sim.PopNMov(stackOffs_target,stackOffs_source
c52
yT3
nO3:;xH3
funcno=opcode-cAbs;assert(funcno<FUNC_AMOUNT);const
FuncDefinition&func=Functions[funcno];sim.Eat(func.c02,tK1
y73}
}
Become(sim.yI2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Produced tree:\n"
;iP
#endif
}
}
#endif
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
#include <assert.h>
#define FP_MUL_COMBINE_EXPONENTS
l93{using
l93
FUNCTIONPARSERTYPES;using
t7
y31
static
void
AdoptChildrenWithSameOpcode(eX{
#ifdef DEBUG_SUBSTITUTIONS
bool
nX2
eV3
#endif
tX3
if(lA4
nC==tG2){
#ifdef DEBUG_SUBSTITUTIONS
if(!nX2){std::cout<<"Before assimilation: "
tO
nX2=true;}
#endif
tree.AddParamsMove(lA4.GetUniqueRef().l52),a);}
#ifdef DEBUG_SUBSTITUTIONS
if(nX2){std::cout<<"After assimilation:   "
tO}
#endif
}
}
t7{eV1
ConstantFolding(eX{tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
void*c93=0
x71"["
<<(&c93)<<"]Runs ConstantFolding for: "
tO
DumpHashes(tree)x71
std::flush;
#endif
if(false){redo:;tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&c93)<<"]Re-runs ConstantFolding: "
tO
DumpHashes(tree);
#endif
}
if(tG2!=cImmed){yP3
p=iO
tree
l64
p
e61&&p
x12
eE2==p
nD2{y0
p
eE2);nE}
if(false){ReplaceTreeWithOne:y0
nW3(1));goto
do_return;ReplaceTreeWithZero:y0
x61;goto
do_return;ReplaceTreeWithParam0:
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before replace: "
x71
std::hex<<'['
tJ1
hash1<<','
tJ1
hash2<<']'<<std::dec
tO
#endif
tree
eA2
0));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After replace: "
x71
std::hex<<'['
tJ1
hash1<<','
tJ1
hash2<<']'<<std::dec
tO
#endif
tL1
cA3(tG2
cZ3
cImmed:lD
l13:lD
cAnd:case
cAbsAnd
eP
bool
cA
eV3
tX3{if(!cD3
a)))cA=true;cG3
a),tG2==cAbsAnd)cZ3
IsNever
tQ
IsAlways:iD);lD
lY1
cA3
iJ
cZ3
0:i3
1:tH2
tG2==cAnd?cNotNot:cAbsNotNot);tL1
yT3
eG2
cAnd||!cA)if(ConstantFolding_AndLogic
y92
y21
cOr:case
cAbsOr
eP
bool
cA
eV3
tX3{if(!cD3
a)))cA=true;cG3
a),tG2==cAbsOr))lU2
i3
l53
iD);lD
lY1
cA3
iJ
cZ3
0
tQ
1:tH2
tG2==cOr?cNotNot:cAbsNotNot);tL1
yT3
eG2
cOr||!cA)if(ConstantFolding_OrLogic
y92
y21
cNot:case
cB3:{xH3
n41
0;eL3
cF3
cZ3
cEqual:n41
i61;lD
i61:n41
cEqual;lD
cLess:n41
cGreaterOrEq;lD
cGreater:n41
cLessOrEq;lD
cLessOrEq:n41
cGreater;lD
cGreaterOrEq:n41
cLess;lD
cNotNot:n41
cNot;lD
cNot:n41
cNotNot;lD
cB3:n41
cAbsNotNot;lD
cAbsNotNot:n41
cB3;break;yT3
y73
if(opposite){tH2
OPCODE(opposite));tree.SetParamsMove(xZ
0).GetUniqueRef().l52));tL1
cA3(lX1
0)iK2
cS1
cZ3
IsAlways
tQ
l53
i3
lY1
eG2
cNot&&GetPositivityInfo(xZ
0))==IsAlways)tH2
cB3
l64
cF3==cIf||cF3==t33{lM1
iftree=xZ
0);l92&ifp1=iftree
l8
1);l92&ifp2=iftree
l8
2
l64
ifp1
eI2
ifp1
cS1{tree.x8
ifp1
nC==cNot?cNotNot:cAbsNotNot);i72
xX2
p1.xY1
cC3
iM1
cE3
p2.nR3
if(ifp2
eI2
ifp2
cS1{tree.x8
tG2);i72);p1.xY1
cC3
eN
ifp2
nC==cNot?cNotNot:cAbsNotNot);p2
eQ
ifp2
xX2
p2.nR3
y21
cNotNot:case
cAbsNotNot:{if(cD3
0)))nE3
cG3
0),tG2==cAbsNotNot)cZ3
IsNever
tQ
IsAlways:i3
lY1
eG2
cNotNot&&GetPositivityInfo(xZ
0))==IsAlways)tH2
cAbsNotNot
l64
cF3==cIf||cF3==t33{lM1
iftree=xZ
0);l92&ifp1=iftree
l8
1);l92&ifp2=iftree
l8
2
l64
ifp1
eI2
ifp1
cS1{tL2
0,iftree
xX2
tree
eQ
ifp1)cC3
iM1
cE3
p2.nR3
if(ifp2
eI2
ifp2
cS1{tree.x8
tG2);i72);p1.xY1;tree
nQ3
tH2
iftree
nC);tL1}
y21
cIf:case
cAbsIf:{if(ConstantFolding_IfOperations
y92
y73
case
cMul:{NowWeAreMulGroup:;AdoptChildrenWithSameOpcode(tree)tV2
nC1=nW3(1)i42
l02=0;bool
nD1=false
nY2
if(!tM1
continue
tV2
immed=nU1
if(immed==x61
goto
ReplaceTreeWithZero;nC1*=immed;++l02;}
if(l02>1||(l02==1&&n23
nC1
iM2(1))))nD1=true;if(nD1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cMul: Will add new "
iU3
nC1<<"\n"
;
#endif
tX3
if(tM1{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iU3
nU1
std::cout<<"\n"
;
#endif
tZ3(!n23
nC1
iM2(1)))tree
eQ
e91
nC1));cA3
iJ
cZ3
0:i3
1:nE3
yT3
if(ConstantFolding_MulGrouping
y92
if(ConstantFolding_MulLogicItems
y92
y21
cAdd
eP
nW3
nC2=0.0
i42
l02=0;bool
nD1=false
nY2
if(!tM1
continue
tV2
immed=nU1
nC2+=immed;++l02;}
if(l02>1||(l02==1&&nC2==x61)nD1=true;if(nD1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cAdd: Will add new "
iU3
nC2<<"\n"
x71"In: "
tO
#endif
tX3
if(tM1{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iU3
nU1
std::cout<<"\n"
;
#endif
tZ3(!(nC2==nW3(0.0)))tree
eQ
e91
nC2));cA3
iJ
cZ3
0
tQ
1:nE3
yT3
if(ConstantFolding_AddGrouping
y92
if(ConstantFolding_AddLogicItems
y92
y21
cMin
eP
size_t
yK2=0;yP3
e5
nY2
while(a+1<tree.GetParamCount()&&lA4
xL
xZ
a+1)))iD+1);yQ3
max.cH3(!e5
xK
known||(p
nD2<e5
nD2){e5
xK
val=p
xK
val;e5
xK
known=true;yK2=a;}
}
if(e5
eF2
tX3{yQ3
min.cH3
a!=yK2&&p
eE2>=e5
nD2
tZ3
iJ==1){nE3
y21
cMax
eP
size_t
yK2=0;yP3
t1
nY2
while(a+1<tree.GetParamCount()&&lA4
xL
xZ
a+1)))iD+1);yQ3
min.cH3(!t1
e61||p
eE2>t1
eE2)){t1
eE2=p
eE2;t1
e61=true;yK2=a;}
}
if(t1
e61){tX3{yQ3
max.cH3
a!=yK2&&(p
nD2<t1
eE2){iD);}
}
}
if
iJ==1){nE3
y21
cEqual:case
i61:case
cLess:case
cGreater:case
cLessOrEq:case
cGreaterOrEq:if(ConstantFolding_Comparison
y92
lD
cAbs:{yP3
p0
eS
0));if
xX1
nE3
if
yN{tH2
cMul);tree.yF
nW3(1)));goto
NowWeAreMulGroup;}
if(cF3==cMul){l92&p=xZ
0);eM
nS3;eM
cL2
xZ1
p.l21++a){p0=iO
p
eJ3;if
xX1{nS3.push_back(p
eJ3;}
if
yN{cL2.push_back(p
eJ3;}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Abs: mul group has "
<<nS3.size()<<" pos, "
<<cL2.size()<<"neg\n"
;
#endif
if(!nS3
tH3||!cL2
tH3){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-Before: "
c72
tree)x71"\n"
<<std::flush;DumpHashes(tree
cM2;
#endif
lM1
eK3;eK3
l84)xZ1
p.l21++a){p0=iO
p
eJ3;if(xX1||yN){}
else
eK3
eQ
p
eJ3;}
eK3
iR2
nT3;nT3
l94
nT3
c0
eK3);nT3
iR2
yB1
cMul);xI3
c0
nT3);yC1
AddParamsMove(nS3
l64!cL2
tH3){if(cL2.size()%2)yC1
yF
nW3(-1)));yC1
AddParamsMove(cL2);}
tree.Become
lJ3);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-After: "
c72
tree
cM2
x71"\n"
<<std::flush;DumpHashes(tree
cM2;
#endif
goto
NowWeAreMulGroup;}
}
y73
#define HANDLE_UNARY_CONST_FUNC(funcname) nT){y0 funcname(lS));nE
case
cLog:iJ3(fp_log);if(cF3==cPow){lM1
pow=xZ
0
l64
GetPositivityInfo(pow
l8
0))==IsAlways
yS3
yB
tree.lV
if(GetEvennessInfo(pow
l8
1))==IsAlways
yS3()eO
abs;abs
l94
abs
c0
pow
xX2
abs.Rehash
yB
pow.n81
0,abs);tree.lV}
iH1
cF3==cAbs){lM1
pow=xZ
0)l8
0
l64
pow
nC==cPow
yS3()eO
abs;abs
l94
abs
c0
pow
xX2
abs.Rehash
yB
pow.n81
0,abs);tree.lV}
lD
cAcosh:iJ3(fp_acosh);lD
cAsinh:iJ3(fp_asinh);lD
cAtanh:iJ3(fp_atanh);lD
cAcos:iJ3(fp_acos);lD
cAsin:iJ3(fp_asin);lD
cAtan:iJ3(fp_atan);lD
cCosh:iJ3(fp_cosh);lD
cSinh:iJ3(fp_sinh);lD
cTanh:iJ3(fp_tanh);lD
cSin:iJ3(fp_sin);lD
cCos:iJ3(fp_cos);lD
cTan:iJ3(fp_tan);lD
cCeil:l74(fp_ceil);lD
cTrunc:l74(fp_trunc);lD
cFloor:l74(fp_floor);lD
cInt:l74(fp_int);lD
cCbrt:iJ3(fp_cbrt);lD
cSqrt:iJ3(fp_sqrt);lD
cExp:iJ3(fp_exp);lD
cLog2:iJ3(fp_log2);lD
cLog10:iJ3(fp_log10);lD
cLog2by
l24
fp_log2(lS)*xZ
i82
cArg:iJ3(fp_arg);lD
cConj:iJ3(fp_conj);lD
cImag:iJ3(fp_imag);lD
cReal:iJ3(fp_real);lD
cPolar
l24
fp_polar
lB4
cMod
l24
fp_mod
lB4
cAtan2:{yP3
p0
eS
yR2
p1
eS
1));nT&&n23
lS,x61){if(p1
cO2
p1
nD2<x61{y0
fp_const_pi
cJ3
if(p1
e61&&p1
eE2>=tO1
x61;nE}
if(eG1
n23
xZ
l81,x61){if(p0
cO2
p0
nD2<x61{y0-fp_const_pihalf
cJ3
if(p0
e61&&p0
eE2>x61{y0
fp_const_pihalf
cJ3}
if
lP
fp_atan2(lS,xZ
l81));nE
if((p1
e61&&p1
eE2>x61||(p1
cO2
p1
nD2<fp_const_negativezero
xH()iZ
yM2;yM2
eN
cPow);yM2
c0
xZ
1));yM2.yF
nW3(-1)));yM2
iR2
yN2;yN2
l84);yN2
c0
xZ
0));yN2
c0
yM2);yN2
tB1);tH2
cAtan);iS2
0,yN2
yI1
1);y21
cPow:{if(ConstantFolding_PowOperations
y92
y73
case
cDiv:nT&&eG1
xZ
l81!=tO1
lS/xZ
i82
cInv:nT&&lS!=tO1
nW3(1)/lS);nE
lD
cSub
l24
lS-xZ
i82
cNeg:nT){y0-lS);nE
lD
cRad:nT){y0
RadiansToDegrees
tZ2
cDeg:nT){y0
DegreesToRadians
tZ2
cSqr:nT){y0
lS*lS);nE
lD
cExp2:iJ3(fp_exp2);lD
cRSqrt:nT){y0
nW3(1)/fp_sqrt
tZ2
cCot:i02
fp_tan
n0
cSec:i02
fp_cos
n0
cCsc:i02
fp_sin
n0
cHypot
l24
fp_hypot
lB4
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:lD
cPCall:case
cFCall:y73
do_return:;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&c93)<<"]Done ConstantFolding, result: "
tO
DumpHashes(tree);
#endif
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
t7{eV1
yP3::set_abs(nN
bool
has_negative=!min.known||min.val<nW3();bool
has_positive=!cI3||i03
nU3);bool
crosses_axis=has_negative&&has_positive;rangehalf
xH
newmax;if(min.cH3
cI3)newmax.set(cL3
i81,i91
l64
crosses_axis)min.set(nW3());cK3
min.cH3
cI3)min.set(fp_min(i81,i91);iH1
min.known)min.set(i81);else
min.set(i91;}
max=newmax;}
eV1
yP3::set_neg(){std::swap(min,max);min.val=-min.val;i03=-i03;}
y41
IsLogicalTrueValue
cI1
yP3&p,bool
abs){if(nB
IsIntType
xH::xC3){if(p
e61&&p
eE2>=nW3(1))l12
if(!abs&&p
x12
xK
val<=nW3(-1))l12}
cK3
p
e61&&p
eE2>=yJ2)l12
if(!abs&&p
x12
xK
val<=nW3(-0.5))l12}
return
t43
y41
IsLogicalFalseValue
cI1
yP3&p,bool
abs){if(nB
IsIntType
xH::xC3){if(abs
lD4
p
xK
known
lF1
1);else
return
p
e61&&p
x12
eE2
nU3-1)lF1
1);}
cK3
abs
lD4
p
xK
known
lF1
0.5);else
return
p
e61&&p
x12
eE2
nU3-0.5)lF1
0.5);}
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;using
t7;l93{nI1
T>inline
int
isnan_workaround(T
t)yR(t!=t);}
}
t7{nQ2
yP3
iO
const
eX
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
{using
l93
FUNCTIONPARSERTYPES;yP3
tmp=CalculateResultBoundaries_do(tree)x71"Estimated boundaries: "
;if(tmp
e61)std::cout<<tmp
eE2;else
std::cout<<"-inf"
x71" .. "
;if(tmp
eF2
std::cout<<tmp
xK
val;else
std::cout<<"+inf"
x71": "
c72
tree)x71
std::endl
iI
tmp;}
nQ2
yP3
lM1::CalculateResultBoundaries_do
cI1
eX
#endif
{iE
yQ1(-fp_const_pihalf
xH(),fp_const_pihalf
xH());iE
pi_limits(-fp_const_pi
xH(),fp_const_pi
xH());iE
abs_pi_limits(y71,fp_const_pi
xH());iE
plusminus1_limits(nW3(-lU3
using
l93
std;eL3
tG2
cZ3
cImmed:nO
l54
iK2.xJ1);case
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cNot:case
cB3:case
cNotNot:case
cAbsNotNot:case
cEqual:case
i61:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:{nO
y71
iM2(1));}
case
cAbs:lE
m.set_abs(nA3
cLog:lE
nY3
fp_log);m
x03
fp_log
nA3
cLog2:lE
nY3
fp_log2);m
x03
fp_log2
nA3
cLog10:lE
nY3
fp_log10);m
x03
fp_log10
nA3
cAcosh:lE
cB
template
set_if<cGreaterOrEq
tP1
fp_acosh)y51
cGreaterOrEq
tP1
fp_acosh
nA3
cAsinh:lE
cB
set(fp_asinh);m
xK
set(fp_asinh
nA3
cAtanh:lE
cB
n4-1),fp_atanh)y51
cLess
tP1
fp_atanh
nA3
cAcos:lE
nO(t61&&(m
nD2<nW3(1))?fp_acos(m
nD2:y71,(nZ3&&(cB
val)>=nW3(-1))?fp_acos(cB
val):fp_const_pi
xH());}
case
cAsin:lE
cB
n4-1),fp_asin,yQ1
eE2)y51
cLess
tP1
fp_asin,yQ1
nD2
iI
m;}
case
cAtan:lE
cB
set(fp_atan,yQ1
eE2);m
xK
set(fp_atan,yQ1
nD2
iI
m;}
case
cAtan2:{nT&&n23
lS,x61)yR
abs_pi_limits;}
if(eG1
n23
xZ
l81,x61)yR
yQ1;}
return
pi_limits;}
case
cSin:lE
bool
nV1=!nZ3||!t61||nV3-cB
val)>=(y5
nV1)eR
nW3
min=fp_mod(cB
val,y5
min<x61
min
yP
nW3
max=fp_mod
nV3,y5
max<x61
max
yP
if(max<min)max
yP
bool
covers_plus1=(min<=fp_const_pihalf
xH()&&max>=fp_const_pihalf
xH());bool
nE1=(min<=cG&&max>=cG
xA1&&nE1)eR
if(nE1)nO
nW3(-1),cL3
y6
xA1)nO
yO2
nW3(1));nO
yO2
cL3
y6);}
case
cCos:lE
if(nZ3)cB
val+=fp_const_pihalf
xH(l64
t61)m
xK
val+=fp_const_pihalf
xH();bool
nV1=!nZ3||!t61||nV3-cB
val)>=(y5
nV1)eR
nW3
min=fp_mod(cB
val,y5
min<x61
min
yP
nW3
max=fp_mod
nV3,y5
max<x61
max
yP
if(max<min)max
yP
bool
covers_plus1=(min<=fp_const_pihalf
xH()&&max>=fp_const_pihalf
xH());bool
nE1=(min<=cG&&max>=cG
xA1&&nE1)eR
if(nE1)nO
nW3(-1),cL3
y6
xA1)nO
yO2
nW3(1));nO
yO2
cL3
y6);}
case
cTan:{nO);}
case
cCeil:lE
m
xK
iH
cFloor:lE
tQ1
iI
m;}
case
cTrunc:lE
tQ1;m
xK
iH
cInt:lE
tQ1;m
xK
iH
cSinh:lE
cB
set(fp_sinh);m
xK
set(fp_sinh
nA3
cTanh:lE
cB
set(fp_tanh,plusminus1_limits.min);m
xK
set(fp_tanh,plusminus1_limits.max
nA3
cCosh:lE
if(nZ3){if(t61){if(cB
val>=y71&&m
xK
val>=x61{cB
val
xW}
iH1(cB
val)<y71&&m
xK
val>=x61{nW3
tmp
xW
if(tmp>m
nD2
m
xK
val=tmp;cA1}
else{cB
val
xW
std::swap(cB
val,m
nD2;}
}
cK3
cB
val>=x61{m
t2
cB
val=fp_cosh(cB
val);}
else{m
t2
cA1}
}
}
else{nZ3=true;cA1
if(t61){cB
val=fp_cosh(m
nD2;m
t2}
else
m
t2}
return
m;}
case
cIf:case
cAbsIf:{yP3
res1
eS
1));yP3
res2
eS
2)l64!res2
e61)res1
e61
eV3
iH1
res1
e61&&(res2
eE2)<res1
eE2)res1
eE2=res2
eE2;iH1
isnan_workaround(res2
eE2))res1
eE2=res2
eE2;if(!res2
eF2
res1
t2
iH1
res1
cO2
res2
nD2>res1
nD2
res1
xK
val=res2
xK
val;iH1
isnan_workaround(res2
nD2)res1
xK
val=res2
xK
val
iI
res1;}
case
cMin:{bool
iK
eV3
bool
iL
eV3
yR3;xA
m
eS
cM3
nZ3)iK
eW3
e61||(cB
val)<yU3)yU3=cB
val;if(!t61)iL
eW3
xK
known||(m
nD2<xC3
nD2
i13=m
xK
val;}
if(iK)lR3
iL)xC3
t2
return
xC3;}
case
cMax:{bool
iK
eV3
bool
iL
eV3
yR3;xA
m
eS
cM3
nZ3)iK
eW3
e61||cB
val>yU3)yU3=cB
val;if(!t61)iL
eW3
xK
known||m
xK
val>xC3
nD2
i13=m
xK
val;}
if(iK)lR3
iL)xC3
t2
return
xC3;}
case
cAdd:{yR3(y71,x61;xA
item
eS
a)l64
item
e61)yU3+=item
eE2;else
lR3
item
eF2
i13+=item
xK
val;else
xC3
t2
if(!lP3&&!xC3
eF2
y73
if(lP3&&x43&&yU3>xC3
nD2
std::swap(yU3,xC3
nD2
iI
xC3;}
case
cMul:{cZ2
Value{enum
x13{cN3,l22,nX3}
;x13
eU
tV2
value;Value(x13
t):eU(t),value(0){}
Value(nW3
v):eU(cN3),value(v){}
bool
cP2
c31
eU==l22||tJ2
value<x61
e33
eI1*=cI1
Value&rhs){if
tJ2
rhs.eU==cN3)value*=rhs.value;else
eU=(cP2)!=rhs.cP2)?l22:nX3);}
iL2<cI1
Value&rhs
c31(eU==l22&&rhs.eU!=l22)||tJ2(rhs.eU==nX3||(rhs.eU==cN3&&value<rhs.value)));}
}
;cZ2
yR1{Value
yP2,yQ2;yR1():yP2(Value::nX3),yQ2(Value::l22){}
void
multiply(Value
cO3,const
Value&value2){cO3*=value2;if(cO3<yP2)yP2=cO3;if(yQ2<cO3)yQ2=cO3;}
}
;yR3(nW3(lU3
xA
item
eS
cM3
item
e61&&!item
eF2
nO);Value
x23=lP3?Value(iC2
min.lZ1
l22);Value
x33=x43?Value(lQ3
lZ1
nX3);Value
x53=item
e61?Value(item.min.lZ1
l22);Value
x63=item
xK
known?Value(item
xK
lZ1
nX3);yR1
range
y61
x23,x53)y61
x23,x63)y61
x33,x53)y61
x33,x63
l64
range.yP2.eU==Value::cN3)yU3=range.yP2.value;else
lR3
range.yQ2.eU==Value::cN3)i13=range.yQ2.value;else
xC3
t2
if(!lP3&&!xC3
eF2
y73
if(lP3&&x43&&yU3>xC3
nD2
std::swap(yU3,xC3
nD2
iI
xC3;}
case
cMod:{yP3
x
eS
yR2
y
eS
1)l64
y
eF2{if(y
xK
val>=x61{if(!x
e61||(x
eE2)<x61
nO-y
xK
val,y
nD2;else
nO
y71,y
nD2;}
cK3!x
xK
known||(x
nD2>=x61
nO
y
xK
val,-y
nD2;else
nO
y
xK
val,fp_const_negativezero
xH());}
}
else
nO);}
case
cPow:{if(eG1
xZ
l81==x61{nO
nW3(lU3}
nT&&lS==x61{nO
y71,x61;}
nT&&n23
lS
nK2
nO
nW3(lU3}
if(eG1
xZ
l81>y71&&GetEvennessInfo(xZ
1))==IsAlways
i12
e62
xZ
l81;yP3
tmp
eS
yR2
xC3;lP3=true;yU3=0;if(tmp
e61&&tmp
eE2>=x61
iC2
min.e72
tmp.min.e82
iH1
tmp
xK
cH3
tmp
xK
val<=x61
iC2
min.e72
tmp
xK
e82
xC3
t2
if(tmp
e61&&tmp
eF2{x43=true;i13=cL3
fp_abs(tmp
eE2),fp_abs(tmp
nD2);lQ3
e72
lQ3
e82}
return
xC3;}
yP3
p0
eS
yR2
p1
eS
1))lS3
p0_positivity=(p0
e61&&(p0
eE2)>=x61?IsAlways:(p0
cO2
p0
nD2<y71?l53
Unknown)lS3
cQ2=GetEvennessInfo(xZ
1))lS3
t3=Unknown;eL3
p0_positivity)lU2
t3=IsAlways;lD
l53{t3=cQ2;y73
yT3
eL3
cQ2)lU2
t3=IsAlways;lD
l53
lD
Unknown:{if(eG1!tD2
xZ
l81)&&xZ
l81>=x61{t3=IsAlways;}
y73}
cA3(t3)lU2{nW3
min=y71;if(p0
e61&&p1
e61){min=t53
p0
eE2,p1
eE2
l64
p0
eE2<y71&&(!p1
xK
known||p1
xK
val>=x61&&min>=x61
min=y71;}
if(p0
e61&&p0
eE2>=y71&&p0
xK
cH3
p1
eF2{nW3
max=t53
p0
xK
val,p1
nD2;if(min>max)std::swap(min,max);nO
min,max);}
nO
min,false);}
case
l53{nO
false,fp_const_negativezero
xH());}
yT3{y73
y21
cNeg:lE
m.set_neg(nA3
cSub:yG
cNeg)t63
1))iF3
cAdd);iG3);tmp
lC4
iV2
cInv:{lO1
lX-1))iV2
cDiv:yG
cInv)t63
1))iF3
y91
lC4
iV2
cRad:yU
y91.yF
fp_const_rad_to_deg
xH())iV2
cDeg:yU
y91.yF
fp_const_deg_to_rad
xH())iV2
cSqr:{lO1
lX
2))iV2
cExp:yU
cPow);tmp.yF
fp_const_e
xH()));iG3
iV2
cExp2:yU
cPow);tmp.yF
lT3
iG3
iV2
cCbrt:lE
cB
set(fp_cbrt);m
xK
set(fp_cbrt
nA3
cSqrt:lE
if(nZ3)cB
val=(cB
val)<y71?0:fp_sqrt(cB
val
l64
t61)m
xK
val=(m
nD2<y71?0:fp_sqrt(m
nD2
iI
m;}
case
cRSqrt:{lO1
lX-0.5))iV2
cHypot:{lM1
xsqr,ysqr,add,sqrt;xsqr.nJ
0));xsqr.yF
lT3
ysqr.nJ
1));ysqr.yF
lT3
xsqr
eN
cPow);ysqr
eN
cPow);add
c0
xsqr);add
c0
ysqr);add
eN
cAdd);sqrt
c0
add);sqrt
eN
cSqrt)iI
iO
sqrt);}
case
cLog2by:yG
cLog2)t63
0))iF3
cMul);tmp
lC4);tmp.nJ
1)iV2
cCot:yG
cTan)nY
lI
cSec:yG
cCos)nY
lI
cCsc:yG
cSin)nY
iO
tmp);}
lD
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:case
l13:lD
cArg:case
cConj:case
cImag:case
cReal:case
cPolar:lD
cPCall:lD
cFCall:y73
nO);}
nQ2
TriTruthValue
GetIntegerInfo
cI1
eX{eL3
tG2
cZ3
cImmed:return
tD2
l54)?IsAlways:IsNever;case
cFloor:case
cCeil:case
cTrunc:case
cInt:return
IsAlways;case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cEqual:case
i61:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:return
IsAlways;case
cIf:{TriTruthValue
a=GetIntegerInfo(xZ
1))lS3
b=GetIntegerInfo(xZ
2)l64
a==b
lD4
a
iI
Unknown;}
case
cAdd:case
cMul:{tX3
if(GetIntegerInfo(lA4)!=IsAlways
lD4
Unknown
iI
IsAlways;}
yT3
y73
return
Unknown;}
y41
IsLogicalValue
cI1
eX{eL3
tG2
cZ3
cImmed:return
n23
l54,x61||n23
l54
iM2(1));case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cAbsAnd:case
cAbsOr:case
cB3:case
cAbsNotNot:case
cEqual:case
i61:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:x5
cMul:{tX3
if(!cD3
a))cV
l12}
case
cIf:case
cAbsIf:yR
cD3
1))nL1
xZ
2));}
yT3
y73
return
t43}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;
#if defined(__x86_64) || !defined(FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS)
# define CBRT_IS_SLOW
#endif
#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif
l93
xC1{extern
const
xH3
char
powi_table[256];}
l93{using
t7
y31
bool
IsOptimizableUsingPowi(long
immed,long
penalty=0){xC1
y83
synth;synth.PushVar(l13)i42
bytecodesize_backup
tS3
GetByteCodeSize();xC1::nT1
immed,xC1::tY1
xH::MulSequence,synth)i42
bytecode_grow_amount
tS3
GetByteCodeSize()-bytecodesize_backup
iI
bytecode_grow_amount<size_t(MAX_POWI_BYTECODE_LENGTH-penalty);}
eV1
ChangeIntoRootChain(lM1&tree,bool
l63,long
i92,long
iA2){while(iA2>0)yU
cCbrt);lN3);tmp
tB1
tN1--iA2;}
while(i92>0)yU
cSqrt
l64
l63){tmp
eN
cRSqrt);l63
eV3}
lN3);tmp
tB1
tN1--i92;}
if(l63)yU
cInv);lN3
tN1}
}
nQ2
cZ2
RootPowerTable{static
const
nW3
RootPowers[(1+4)*(1+3)];}
y31
const
nW3
t9(1+4)*(1+3)]={nW3(1)lT
i23
i23
2*i23
2*2*2)lT
3)lT
3*2)lT
3*2*2)lT
3*iO1
2*iO1
3
iW2
2
iW2
2*2
iW2
iO1
3*2*iO1
3*3
iW2
3*2
iW2
3*2*2
iW2
3*iO1
3*3*2*2*2*2)}
;cZ2
PowiResolver{static
const
xH3
MaxSep=4;static
x83
i33=5
iD3
int
cP3
iD3
long
xF3
iD3
long
tW;cZ2
yS2{yS2():n_int_sqrt(0),n_int_cbrt(0),sep_list(),n01(0){}
int
n_int_sqrt;int
n_int_cbrt;int
tX1
MaxSep];tW
n01;}
y31
static
yS2
CreatePowiResult(nW3
xG3){yS2
xC3;cP3
tA=FindIntegerFactor
e92
l64
tA==0){
#ifdef DEBUG_POWI
iB2"no factor found for %Lg\n"
,(cQ3);
#endif
return
xC3;}
cR3=yA1
xG3,tA);xF3
eT2=EvaluateFactorCost(tA,0,0,0)+cF
cR3);int
i43=0;int
i53=0;int
x73=0;
#ifdef DEBUG_POWI
iB2"orig = %Lg\n"
,(cQ3);iB2"plain factor = "
iW3"%ld\n"
,(int)tA,(long)eT2);
#endif
for
x91
n_s=0;n_s<MaxSep;++n_s){int
xF=0;xF3
yS1=eT2;cP3
c01=tA;for(int
s=1;s<i33*4;++s){
#ifdef CBRT_IS_SLOW
if(s>=i33)break;
#endif
int
n_sqrt=s%i33;int
n_cbrt=s/i33;if(n_sqrt+n_cbrt>4)continue
tV2
lG1=xG3;lG1-=t9
s];iB1=FindIntegerFactor(lG1
l64
xS2!=0){tW
xR=yA1
lG1,xS2);xF3
cost=EvaluateFactorCost(xS2,i43+n_sqrt,i53+n_cbrt,x73+1)+cF
xR);
#ifdef DEBUG_POWI
iB2"Candidate sep %u (%d*sqrt %d*cbrt)factor = "
iW3"%ld (for %Lg to %ld)\n"
,s,n_sqrt,n_cbrt,xS2,(long)cost
eH2
lG1,(long)xR);
#endif
if(cost<yS1){xF=s;c01=xS2;yS1=cost;}
}
}
if(!xF)break;
#ifdef DEBUG_POWI
iB2"CHOSEN sep %u (%d*sqrt %d*cbrt)factor = "
iW3"%ld, exponent %Lg->%Lg\n"
,xF,xF%i33,xF/i33,c01,yS1
eH2
e92)eH2
e92-t9
xF]));
#endif
iC2
tX1
n_s]=xF;xG3-=t9
xF];i43+=xF%i33;i53+=xF/i33;eT2=yS1;tA=c01;x73+=1;}
cR3=yA1
xG3,tA);
#ifdef DEBUG_POWI
iB2"resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n"
,cR3,(cQ3
eH2
tA);
#endif
while(tA%2==0){++iC2
n_int_sqrt;tA/=2;}
while(tA%3==0){++iC2
n_int_cbrt;tA/=3;}
return
xC3;}
private:static
xF3
cF
tW
xR){static
std::map
cR2
i7;if(xR<0){xF3
cost=22
iI
cost+cF-xR);}
std::map
cR2::yK3
i=i7.xT2
xR
l64
i!=i7.cY1
xR
lD4
i
e52;std::pair
cR2
xC3(xR,0.0);xF3&cost=iC2
second;while(xR>1){int
xS2=0;if(xR<256){xS2=xC1::powi_table[xR];if(xS2&128)xS2&=127;else
xS2=0;if(xS2&64)xS2=-(xS2&63)-1;}
if(xS2){cost+=cF
xS2);xR/=xS2;t02
if(!(xR&1)){xR/=2;cost+=6;}
else{cost+=7;xR-=1;}
}
i7.yG3,xC3)iI
cost;}
cB1
tW
yA1
yT1,iB1)yR
makeLongInteger(value*nW3(xS2));}
cB1
bool
yU1
yT1,iB1
i12
v=value*nW3(xS2)iI
isLongInteger(v);}
cB1
cP3
FindIntegerFactor(yT1){iB1=(2*2*2*2);
#ifdef CBRT_IS_SLOW
#else
xS2*=(3*3*3);
#endif
cP3
xC3=0;if(yU1
value,xS2)){xC3=xS2;while((xS2%2)==0&&yU1
value,xS2/2))xC3=xS2/=2;while((xS2%3)==0&&yU1
value,xS2/3))xC3=xS2/=3;}
#ifdef CBRT_IS_SLOW
if(xC3==0){if(yU1
value,3
tQ2
3;}
#endif
return
xC3;}
static
int
EvaluateFactorCost(int
xS2,int
s,int
c,int
nmuls){x83
x93=6;
#ifdef CBRT_IS_SLOW
x83
eU2=25;
#else
x83
eU2=8;
#endif
int
xC3=s*x93+c*eU2;while(xS2%2==0){xS2/=2;xC3+=x93;}
while(xS2%3==0){xS2/=3;xC3+=eU2;}
xC3+=nmuls
iI
xC3;}
}
;}
t7{y41
lM1::RecreateInversionsAndNegations(bool
prefer_base2){bool
changed=false
xZ1
l21++a)if(lW1.RecreateInversionsAndNegations(prefer_base2))yW1
if(changed){exit_changed:Mark_Incompletely_Hashed(nZ2
eL3
lC2{case
cMul:{eM
nE2
eO
nF2,cT1;if(true){bool
nF1=false
tV2
xH2=0;xA3
i63
0)i73
tH
cT2)){nF1=true;xH2=tH
l81;y73}
if(nF1
i12
immeds=1.0;xA3
cV2)){immeds*=powgroup.xJ1;yV1}
for
e0-->0;){lM1&powgroup=lW1;if(powgroup
i63
0)i73
tH
cT2
iZ&log2=tH
0);log2.l71
log2
eN
xB3
log2.yF
t53
immeds
iM2(1)/xH2)));log2
tB1
c52}
}
}
xA3
i63
cT2)){l92&exp_param=tH
1)tV2
e62
exp_param.xJ1;if(e01
iM2(-1))){l71
nE2.push_back(lW1
xX2
yV1
iH1
xG3<y71&&isInteger
e92
iZ
iM;iM
eN
cPow);iM
eQ
tH
0));iM.yF-xG3));iM
tB1);nE2.push_back(iM)iH3}
iH1
powgroup
i73!nF2
x32)){nF2=tH
0)iH3
iH1
powgroup
nC==cLog2by&&!cT1
x32)){cT1=powgroup
iH3}
if(!nE2
tH3){changed=true
eO
tG1
iI3
cMul);tG1.SetParamsMove(nE2);tG1
iR2
yB1
cMul);yC1
SetParamsMove
eK
if(yC1
IsImmed()&&n23
yC1
xJ1
nK2
nG2
cInv);eB
tG1);}
cK3
yC1
y22>=tG1.y22){nG2
cDiv
tF1
eB
tG1);}
else{nG2
cRDiv);eB
tG1
tF1}
}
}
if(nF2
x32
iZ
yB1
lC2;yC1
SetParamsMove
eK
while(yC1
RecreateInversionsAndNegations(prefer_base2))yC1
FixIncompleteHashes();nG2
xB3
eB
nF2
tF1
yW1}
if(cT1
x32
iZ
yB1
cMul);xI3
c0
cT1
l8
1));yC1
AddParamsMove
eK
while(yC1
RecreateInversionsAndNegations(prefer_base2))yC1
FixIncompleteHashes();DelParams();nG2
xB3
eB
cT1
l8
0)tF1
yW1
y21
cAdd:{eM
iD2;for
e0-->0;)if(cT3
cMul){nH2
yD1:eO&xI3
yT2
eY2
b=yC1
l21
b-->0;){if
lJ3
l8
b).lS1
xS2=xI3
l8
b).xJ1;lL3
xS2
nI2
yD1;}
yC1
l71
yC1
l83
b
cS2
iH1
n23
xS2
iM2(-2)))xO
yD1;}
yC1
l71
yC1
l83
b);yC1
yF
nW3(2))cS2}
}
if(t4){yC1
tB
xI3);yV1}
iH1
cT3
cDiv&&!IsIntType
xH::xC3){nH2
yE1:eO&tG1
yT2
if(tG1
l8
0).nJ2
n23
tG1
l8
0).xJ1
nI2
yE1;}
tG1.l71
tG1.l83
0)iI3
cInv
cS2}
if(t4)xO
yE1;}
tG1.tB
tG1);yV1}
iH1
cT3
cRDiv&&!IsIntType
xH::xC3){nH2
x21:eO&tG1
yT2
if(tG1
l8
1).nJ2
n23
tG1
l8
l81
nI2
x21;}
tG1.l71
tG1.l83
1)iI3
cInv
cS2}
if(t4)xO
x21;}
tG1.tB
tG1);yV1}
if(!iD2
tH3){
#ifdef DEBUG_SUBSTITUTIONS
iB2"Will make a Sub conversion in:\n"
);fflush(stdout);iP
#endif
lM1
yU2;yU2
eN
cAdd);yU2.SetParamsMove(iD2);yU2
iR2
cU1;cU1
eN
cAdd);cU1.SetParamsMove(l52));cU1
tB1
l64
cU1
cV2)&&n23
cU1.xJ1,x61){nG2
cNeg);tW1);}
cK3
cU1.y22==1){nG2
cRSub);tW1
lF3}
iH1
yU2
nC==cAdd){nG2
cSub
lF3
tW1
xX2
eQ1
1;a<yU2.l21++a){lM1
eV2;eV2
eN
cSub);eV2.SetParamsMove(l52));eV2
tY2;eB
eV2);tW1
eJ3;}
}
else{nG2
cSub
lF3
tW1);}
}
#ifdef DEBUG_SUBSTITUTIONS
iB2"After Sub conversion:\n"
);fflush(stdout);iP
#endif
y21
cPow:{l92&p0
iE2
0);l92&p1
iE2
1
l64
p1.nJ2
p1.xJ1!=y71&&!tD2
p1.xJ1)){eG
yS2
r=eG
CreatePowiResult(fp_abs
eE3)l64
r.n01!=0){bool
l32
eV3
if
eE3<y71&&r.tX1
0]==0&&r.n_int_sqrt>0){l32=true;}
#ifdef DEBUG_POWI
iB2"Will resolve powi %Lg as powi(chain(%d,%d),%ld)"
eH2
fp_abs
eE3),r.n_int_sqrt,r.n_int_cbrt,r.n01);for
x91
n=0;n<eG
MaxSep;++n){if(r
lV3==0)break;int
n_sqrt=r
lV3%eG
i33;int
n_cbrt=r
lV3/eG
i33;iB2"*chain(%d,%d)"
,n_sqrt,n_cbrt);}
iB2"\n"
);
#endif
lM1
cW2
iE2
0)eO
yV2=cW2;yV2.l71
ChangeIntoRootChain(yV2,l32,r.n_int_sqrt,r.n_int_cbrt);yV2
iR2
pow;if(r.n01!=1){pow
eN
cPow);pow
c0
yV2);pow.yF
nW3(r.n01)));}
else
pow.swap(yV2)eO
mul;mul
l84);mul
c0
pow);for
x91
n=0;n<eG
MaxSep;++n){if(r
lV3==0)break;int
n_sqrt=r
lV3%eG
i33;int
n_cbrt=r
lV3/eG
i33
eO
eW2=cW2;eW2.l71
ChangeIntoRootChain(eW2,false,n_sqrt,n_cbrt);eW2
tB1);mul
c0
eW2);}
if
eE3<y71&&!l32){mul
tB1);nG2
cInv);n81
0,mul);l83
1);}
else{nG2
cMul);SetParamsMove(mul.l52));}
#ifdef DEBUG_POWI
iP
#endif
yW1
y73}
}
if(GetOpcode()==cPow&&(!p1.c42!isLongInteger
eE3)||!IsOptimizableUsingPowi
xH(makeLongInteger
eE3)))){if(p0
cV2)&&p0.xJ1
nU3
0.0)){if(prefer_base2
i12
yW2=fp_log2(p0.xJ1);lL3
yW2
nK2
l83
0);}
else{n1
e91
yW2));xG3
eQ
p1)e11
l91}
nG2
cExp2);yW1}
else{nW3
yW2=fp_log(p0.xJ1);lL3
yW2
nK2
l83
0);}
else{n1
e91
yW2));xG3
eQ
p1)e11
l91}
nG2
cExp);yW1}
}
iH1
GetPositivityInfo(p0)==IsAlways){if(prefer_base2){lM1
log;log
eN
cLog2);log
eQ
p0);log
tB1);n1
p1);xG3
c0
log)e11);nG2
cExp2
l91
yW1}
else{lM1
log;log
eN
cLog);log
eQ
p0);log
tB1);n1
p1);xG3
c0
log)e11);nG2
cExp
l91
yW1}
}
y21
cDiv:{if(GetParam(0)cV2)&&n23
GetParam(0).xJ1
nK2
nG2
cInv);l83
0);}
y73
yT3
y73
if(changed)goto
exit_changed
iI
changed;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l93
FUNCTIONPARSERTYPES;l93{using
t7;class
iJ1{size_t
nG1
i42
eH
i42
eI
i42
lH1
i42
t5
i42
t6
i42
n51;eO3
iJ1():nG1(0),eH(0),eI(0),lH1(0),t5(0),t6(0),n51(0){}
void
eF3
OPCODE
op){nG1+=1
iI1
cCos)++eH
iI1
cSin)++eI
iI1
cSec)++eH
iI1
cCsc)++eI
iI1
cTan)++lH1
iI1
cCot)++lH1
iI1
cSinh)++t6
iI1
cCosh)++t5
iI1
cTanh)++n51;}
size_t
GetCSEscore
e83
size_t
xC3=nG1
iI
xC3;}
int
NeedsSinCos
e83
bool
yF1=(nG1==(eH+eI+lH1)l64(lH1&&(eI||eH))||(eI&&eH)){if(yF1
lD4
1
iI
2;}
return
0;}
int
NeedsSinhCosh
e83
bool
yF1=(nG1==(t5+t6+n51)l64(n51&&(t6||t5))||(t6&&t5)){if(yF1
lD4
1
iI
2;}
return
0;}
size_t
MinimumDepth
e83
size_t
n_sincos=std::min(eH,eI)i42
n_sinhcosh=std::min(t5,t6
l64
n_sincos==0&&n_sinhcosh==0
lD4
2
iI
1;}
}
y31
class
TreeCountType:public
std::multimap<fphash_t,std::pair<iJ1,lM1> >{}
y31
void
FindTreeCounts(tR1&nM2,l92&tree,OPCODE
xI2,bool
skip_root=false){e1
i=nM2.xT2
tree.GetHash()l64!skip_root){bool
found
eV3
for(;i!=nM2.cY1
tree.GetHash();++i){if(tree
xL
i
e52
t03)){i
e52
tM2
eF3
xI2);found=true;y73}
if(!found){iJ1
count;count.eF3
xI2);nM2.yG3,std::make_pair(tree.GetHash(),std::make_pair
cX3
iK2)));}
}
for
y11
c4++a)FindTreeCounts(nM2,lA4,tG2);}
cZ2
yV{bool
BalanceGood;bool
FoundChild;}
y31
yV
lI1
l92&root,l92&child){if(root
xL
child)){yV
xC3={true,true}
iI
xC3;}
yV
xC3={true,false}
;if(root
nC==cIf||root
nC==t33{yV
cond=lI1
root
l8
0
xW3
yV
xY=lI1
root
l8
1
xW3
yV
yA=lI1
root
l8
2
xW3
if
nF3||xY
yW||yA
yW){xC3
yW=true;}
xC3
eC=((xY
yW==yA
yW)||nF3
cX2&&(cond
eC||(xY
yW&&yA
yW))&&(xY
eC||nF3
cX2&&(yA
eC||nF3
cX2;}
else{bool
iC1
eV3
bool
nH1
eV3
eY2
b=root.GetParamCount(),a=0;a<b
tO3
yV
tmp=lI1
root
l8
a
xW3
if(tmp
yW)xC3
yW=true;if(tmp
eC==false)iC1=true;iH1
tmp
yW)nH1=true;}
if(iC1&&!nH1)xC3
eC
eV3}
return
xC3;}
y41
iT2
e42
i83
l92&tree,const
xC1
y83&synth,const
tR1&nM2){eY2
b=tree.GetParamCount(),a=0;a<b
tO3
l92&leaf=lA4;e1
synth_it;tY3
tR1::const_iterator
i=nM2.yJ3
i!=nM2.end();++i){if(i->first!=leaf.GetHash())lI3
const
iJ1&occ
x42
first
i42
score=occ.GetCSEscore();l92&candidate
x42
second;nP2
candidate))lI3
if(leaf.y22<occ.MinimumDepth())lI3
if(score<2)lI3
if(lI1
i83
leaf)eC==false)continue
lB3
if(iT2(i83
leaf,synth,nM2))l12}
return
t43
y41
nO2
e42
yH3,l92&expr){c11
yH3
l8
a)xL
expr))l12
c11
nO2(yH3
l8
a),expr
tQ2
true
iI
t43
y41
GoodMomentForCSE
e42
yH3,l92&expr){if(yH3
nC==cIf)l12
c11
yH3
l8
a)xL
expr))l12
size_t
iF2=0;c11
nO2(yH3
l8
a),expr))++iF2
iI
iF2!=1;}
}
t7{nQ2
size_t
lM1::SynthCommonSubExpressions(xC1::y81
const{if(GetParamCount()==0
lD4
0
i42
stacktop_before
tS3
GetStackTop();tR1
nM2;FindTreeCounts(nM2,*this,GetOpcode(),true);for(;;){size_t
yX2=0;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Finding a CSE candidate, root is:"
<<std::endl;DumpHashes(*this);
#endif
e1
cs_it(nM2.end());for(e1
j=nM2.yJ3
j!=nM2.end();){e1
i(j++);const
iJ1&occ
x42
first
i42
score=occ.GetCSEscore();l92&tree
x42
second;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Score "
<<score<<":\n"
<<std::flush;DumpTreeWithIndent(tree);
#endif
nP2
tree))y8
if(tree.y22<occ.MinimumDepth())y8
if(score<2)y8
if(lI1*this
iK2)eC==false)y8
if(iT2(*this
iK2,synth,nM2)){t02
if(!GoodMomentForCSE(*this
iK2))y8
score*=tree.y22;if(score>yX2){yX2=score;cs_it=i;}
}
if(yX2<=0){
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"No more CSE candidates.\n"
<<std::flush;
#endif
y73
l92&tree=cs_it
e52
t03;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<l14"Common Subexpression:"
;DumpTree
xH(tree)x71
std::endl;
#endif
#if 0
int
n11=occ.NeedsSinCos();int
tX=occ.NeedsSinhCosh()eO
iG2,iH2,yY2,yZ2;if(n11){iG2
eX2
iG2
eN
cSin);iG2
tB1);iH2
eX2
iH2
eN
cCos);iH2
tB1);nP2
iG2)||synth.Find(iH2)){if(n11==2){nN2
t02
n11=0;}
}
if(tX){yY2
eX2
yY2
eN
cSinh);yY2
tB1);yZ2
eX2
yZ2
eN
cCosh);yZ2
tB1);nP2
yY2)||synth.Find(yZ2)){if(tX==2){nN2
t02
tX=0;}
}
#endif
tree.SynthesizeByteCode(synth,false);nN2
#ifdef DEBUG_SUBSTITUTIONS_CSE
synth.template
Dump<0>()x71"Done with Common Subexpression:"
;DumpTree
xH(tree)x71
std::endl;
#endif
#if 0
if(n11){if(n11==2||tX){synth.eH1}
lE3
cSinCos,1,2)yO1
iG2,1)yO1
iH2,0);}
if(tX){if(n11)synth.eH1
if(tX==2){synth.eH1}
lE3
cSinhCosh,1,2)yO1
yY2,1)yO1
yZ2,0);}
#endif
}
return
synth.xI
stacktop_before;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
nQ2
lL1
xH::iI2{using
t7;cN2()eO
tree;tree.GenerateFrom(*mData);FPoptimizer_Optimize::ApplyGrammars(tree);std
xV3<xH3>cU3;std
xV3
xH
immed
i42
stacktop_max=0;tree.SynthesizeByteCode(cU3,immed,stacktop_max
l64
mData->mStackSize!=stacktop_max){mData->mStackSize=xH3(stacktop_max);
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
mData->mStack.eP3
stacktop_max);
#endif
}
mData->mByteCode.swap(cU3);mData->mImmed.swap(immed);}
#define FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(type) tT1>lL1<type>::iI2{}
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
iK3(MpfrFloat)
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
iK3(GmpInt)
#endif
#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
iK3(std::complex<double>)
#endif
#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
iK3(std::complex<float>)
#endif
#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
iK3(std::complex<long
double>)
#endif
#define FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(type) template lL1<type>::iI2;
#ifndef FP_DISABLE_DOUBLE_TYPE
iL3(double)
#endif
#ifdef FP_SUPPORT_FLOAT_TYPE
iL3(float)
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
iL3(long
double)
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
iL3(long)
#endif
#endif // FP_SUPPORT_OPTIMIZER

#endif
