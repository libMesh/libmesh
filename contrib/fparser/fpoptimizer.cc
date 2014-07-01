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
#define lA4 lO2 tJ1
#define l94 0].info
#define l84 tM3 lQ[b
#define l74 tK info
#define l64 b.Value)
#define l54 b.Opcode
#define l44 )val=iZ2
#define l34 ){std::
#define l24 ,cAtan2
#define l14 tree y11
#define l04 {lT1 xQ=
#define iZ3 "Found "
#define iY3 stackpos
#define iX3 ,lT 2,
#define iW3 .end()
#define iV3 "dup(%u) "
#define iU3 eU{assert
#define iT3 "%d, cost "
#define iS3 "immed "<<
#define iR3 mFuncParsers
#define iQ3 ;DumpTree(
#define iP3 stderr
#define iO3 (constraints&
#define iN3 sep2=" "
#define iM3 FPHASH_CONST
#define iL3 cache_needed[
#define iK3 "PUSH " iQ3
#define iJ3 fprintf
#define iI3 ::cout<<"Applying "
#define iH3 FUNCTIONPARSER_INSTANTIATE_OPTIMIZE
#define iG3 FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE
#define iF3 HANDLE_UNARY_CONST_FUNC
#define iE3 {data->
#define iD3 second;
#define iC3 within,
#define iB3 AddFrom(
#define iA3 {size_t
#define i93 c_count
#define i83 s_count
#define i73 MaxOp
#define i63 2)lS 2*
#define i53 max.val
#define i43 e9 ifp2
#define i33 sim.xL;
#define i23 ].swap(
#define i13 e12 i31
#define i03 ;typedef
#define tZ3 (Value_t
#define tY3 =synth.
#define tX3 codes[b
#define tW3 whydump
#define tV3 nparams
#define tU3 cTan eY
#define tT3 cLog eY
#define tS3 l4 0,1,
#define tR3 0x12 nN
#define tQ3 cHypot,
#define tP3 nR 0,
#define tO3 cAbs nR
#define tN3 xE(i92
#define tM3 info.
#define tL3 *x6)[a].
#define tK3 continue;}
#define tJ3 ;unsigned
#define tI3 false;}
#define tH3 lF1 y5.
#define tG3 cAbsIf)
#define tF3 eV nC==
#define tE3 ));tmp
#define tD3 params
#define tC3 Immeds
#define tB3 {pow.c02
#define tA3 );tmp2.nI
#define t93 }switch(
#define t83 iD cMul);
#define t73 :start_at()
#define t63 NumConstant:
#define t53 FP_GetOpcodeName
#define t43 &&IsLogicalValue(
#define t33 switch(nJ3.first){case
#define t23 ,const cZ&
#define t13 .size()
#define t03 ].second
#define eZ3 ].first
#define eY3 return p
#define eX3 }return
#define eW3 Ne_Mask
#define eV3 Gt_Mask
#define eU3 Lt_Mask
#define eT3 opcode,
#define eS3 public:
#define eR3 yC val=
#define eQ3 yC known)
#define eP3 =fp_pow(
#define eO3 =true;c81
#define eN3 =eX a)yG3
#define eM3 (mulgroup
#define eL3 ;}case
#define eK3 pclone
#define eJ3 ,1,0,{1,
#define eI3 cOr,l6
#define eH3 newpow
#define eG3 change
#define eF3 (count
#define eE3 133,2,
#define eD3 Params
#define eC3 Needs
#define eB3 byteCode
#define eA3 eL e31);
#define e93 nF2 a eQ
#define e83 .SetParamsMove(
#define e73 lZ1 nC==
#define e63 cLog2by
#define e53 c23 n31
#define e43 long i71
#define e33 factor_t
#define e23 value1
#define e13 else{if(
#define e03 i2 p2;p2
#define cZ3 {tree.x4
#define cY3 cAbsNot
#define cX3 stackptr
#define cW3 cLog);xF
#define cV3 l8 0));
#define cU3 opcodes
#define cT3 did_muli
#define cS3 &Value){
#define cR3 yL const
#define cQ3 used[b]
#define cP3 :{n51 r=
#define cO3 size_t n
#define cN3 sizeof(
#define cM3 cAbsIf,
#define cL3 cNotNot,
#define cK3 l4 16,1,
#define cJ3 281856,
#define cI3 cLess,cL
#define cH3 middle2
#define cG3 ::string
#define cF3 GetIntegerInfo(
#define cE3 (tree)!=
#define cD3 (yC val
#define cC3 l4 20,1,
#define cB3 l4 4,1,
#define cA3 ,i6 2,
#define c93 cO2 2,2,
#define c83 cO2 0,2,
#define c73 450998,
#define c63 cExp2 nR
#define c53 result
#define c43 c53 i4
#define c33 c53 eL3
#define c23 c53.
#define c13 c43=false;if(
#define c03 c53 nP
#define yZ3 c53 xK1
#define yY3 xE c53
#define yX3 ;Value_t
#define yW3 ,tree,info
#define yV3 y63 size()
#define yU3 lJ 2},0,
#define yT3 default:
#define yS3 n21;case
#define yR3 Ge0Lt1
#define yQ3 Gt0Le1
#define yP3 .what t3
#define yO3 .Become(
#define yN3 cAdd lK2
#define yM3 y41;++b)
#define yL3 n03 0;b<
#define yK3 ,xB3 l8
#define yJ3 );range xE
#define yI3 tO1 1));
#define yH3 eF2 1)yI3
#define yG3 );if(
#define yF3 tO1 1)))
#define yE3 iterator
#define yD3 ContainsOtherCandidates(
#define yC3 begin();
#define yB3 TreeSet
#define yA3 parent
#define y93 insert(i
#define y83 newrel
#define y73 ,lT 1,0,
#define y63 eD3.
#define y53 break;}
#define y43 synth lE2
#define y33 synth.xG
#define y23 b_needed
#define y13 cachepos
#define y03 half=
#define xZ3 131,4,1,
#define xY3 .n7 synth.
#define xX3 ;}static x32
#define xW3 ,iM,1,lE1+1);
#define xV3 131,8,1,
#define xU3 4,1,2,1,
#define xT3 ::vector
#define xS3 1 y7 n61
#define xR3 FindPos(
#define xQ3 src_pos
#define xP3 reserve(
#define xO3 const std::eS
#define xN3 const char*
#define xM3 yF1 void
#define xL3 treeptr
#define xK3 .resize(
#define xJ3 const tR
#define xI3 tL1 void
#define xH3 ImmedTag
#define xG3 a,const
#define xF3 RefCount
#define xE3 Birth();
#define xD3 typename
#define xC3 template
#define xB3 leaf1
#define xA3 cost_t
#define x93 fpdata
#define x83 middle
#define x73 };enum
#define x63 for iY nY
#define x53 sqrt_cost
#define x43 const int
#define x33 mul_count
#define x23 TreeCountItem
#define x13 if(op==
#define x03 eF2 2)));
#define nZ3 maxValue1
#define nY3 minValue1
#define nX3 maxValue0
#define nW3 minValue0
#define nV3 ValueType
#define nU3 yC n3 0),
#define nT3 max.known
#define nS3 abs_mul
#define nR3 pos_set
#define nQ3 x81);}if(
#define nP3 default_function_handling
#define nO3 {sim.Eat(
#define nN3 subtree
#define nM3 invtree
#define nL3 MakeHash(
#define nK3 yQ false;
#define nJ3 parampair
#define nI3 rulenumit
#define nH3 a<tree.x9
#define nG3 cAnd l3
#define nF3 cAnd,l6
#define nE3 l0 yJ2
#define nD3 l0 2,
#define nC3 y53 switch(
#define nB3 &param=*(
#define nA3 MakeEqual
#define n93 nH1,l5::
#define n83 nH1,{l5::
#define n73 newbase
#define n63 branch1op
#define n53 branch2op
#define n43 overlap
#define n33 truth_b
#define n23 truth_a
#define n13 found_dup
#define n03 size_t b=
#define lZ3 {tM1 lE
#define lY3 rangeutil
#define lX3 Plan_Has(
#define lW3 StackMax)
#define lV3 ;range.xP2
#define lU3 );xQ e83
#define lT3 ;}void
#define lS3 )lT3
#define lR3 c53 xS
#define lQ3 namespace
#define lP3 ::res,b8<
#define lO3 (cond c1
#define lN3 ),child);
#define lM3 ;if(half
#define lL3 tB=!tB;}
#define lK3 )lS 3*3*
#define lJ3 inverted
#define lI3 IsNever:
#define lH3 .known&&
#define lG3 iftree
#define lF3 goto e1
#define lE3 if(i42==
#define lD3 cPow);lC
#define lC3 depcodes
#define lB3 explicit
#define lA3 ,tree nX
#define l93 VarBegin
#define l83 ,o);o<<"\n";
#define l73 ].data);
#define l63 ?0:1))l7
#define l53 ;x0 l5::
#define l43 ))break yX3
#define l33 begin(),
#define l23 cond_add
#define l13 cond_mul
#define l03 cond_and
#define iZ2 func lQ1
#define iY2 bool eQ1
#define iX2 Optimize()
#define iW2 costree
#define iV2 sintree
#define iU2 leaf_count
#define iT2 &&cond eD)
#define iS2 sub_params
#define iR2 printf(
#define iQ2 cbrt_count
#define iP2 sqrt_count
#define iO2 lO2 n3 0),
#define iN2 p1 e9 ifp1
#define iM2 pcall_tree
#define iL2 after_powi
#define iK2 yH1.SubTrees
#define iJ2 yH1.Others
#define iI2 tD3)
#define iH2 grammar
#define iG2 cCos eY
#define iF2 cEqual,
#define iE2 cLog nR
#define iD2 0x12},{{3,
#define iC2 cNeg,lT 1,
#define iB2 x7 lT 2,
#define iA2 },{{1,
#define i92 ),0},{
#define i82 std::move(
#define i72 data;data.
#define i62 iD cond nC
#define i52 tree iD
#define i42 tree nC
#define i32 MakeNEqual
#define i22 .GetHash().
#define i12 ,cMul l3
#define i02 (half&63)-1;
#define tZ2 t13;++
#define tY2 (lR));nD lC
#define tX2 );synth.yT
#define tW2 .second);
#define tV2 tree l8 1)
#define tU2 tree l8 2)
#define tT2 (i42)
#define tS2 Dump(std::
#define tR2 isInteger(
#define tQ2 lT1 r;r iD
#define tP2 Comparison
#define tO2 needs_flip
#define tN2 synth);
#define tM2 iF1 apos==
#define tL2 value]
#define tK2 ~size_t(0)
#define tJ2 xT1 xU+1);
#define tI2 n21 true;}
#define tH2 e9 tree);
#define tG2 mul_item
#define tF2 innersub
#define tE2 cbrt_cost
#define tD2 best_cost
#define tC2 n21 lH
#define tB2 fp_min(y3
#define tA2 IsAlways
#define t92 condition
#define t82 TopLevel)
#define t72 {std::cout<<
#define t62 ,xE1);lC
#define t52 (xB3 l8 1)xF1
#define t42 }eX3
#define t32 per_item
#define t22 item_type
#define t12 first2
#define t02 l4 18,1,
#define eZ2 cGreater,
#define eY2 cIf,t4 3,
#define eX2 ,cCosh nR
#define eW2 ,cPow eY
#define eV2 lJ 1},0,
#define eU2 eZ 1},0,
#define eT2 yB2)));x0
#define eS2 Decision
#define eR2 not_tree
#define eQ2 group_by
#define eP2 (param.data
#define eO2 ,(long double)
#define eN2 ->second
#define eM2 targetpos
#define eL2 ParamSpec
#define eK2 rhs.hash2;}
#define eJ2 rhs.hash1
#define eI2 struct
#define eH2 Forget()
#define eG2 {SetOpcode(
#define eF2 Value_t(
#define eE2 exponent
#define eD2 source_tree
#define eC2 .n_int_sqrt
#define eB2 tO1-1))){xO
#define eA2 nC==cLog2&&
#define e92 <tI,xA3>
#define e82 p1_evenness
#define e72 (list.first
#define e62 e42 IsNever x02
#define e52 e42 tA2;if(
#define e42 ))return
#define e32 ,eE2);
#define e22 isNegative(
#define e12 (c53
#define e02 neg_set
#define cZ2 cNop,cNop}}
#define cY2 cTanh,cNop,
#define cX2 NewHash
#define cW2 >eI2 cC<
#define cV2 matches
#define cU2 info=(*x6)[
#define cT2 *x6 iI,info
#define cS2 .match_tree
#define cR2 eZ2 cL
#define cQ2 cSin eY
#define cP2 cTan nR
#define cO2 ,cPow,l2
#define cN2 cCos nR
#define cM2 negated
#define cL2 Specializer
#define cK2 yF2 c53(
#define cJ2 *)&*start_at;
#define cI2 const lT1
#define cH2 const void*
#define cG2 coshtree
#define cF2 )continue
#define cE2 sinhtree
#define cD2 best_score
#define cC2 mulvalue
#define cB2 pow_item
#define cA2 subgroup
#define c92 IsDefined(
#define c82 nC==cPow&&tQ
#define c72 PowiResult
#define c62 maxValue
#define c52 minValue
#define c42 const Value&
#define c32 cJ1 fp_floor
#define c22 div_tree
#define c12 pow_tree
#define c02 CopyOnWrite(
#define yZ2 preserve
#define yY2 nC==cNot||
#define yX2 sim.xL nV 1,
#define yW2 PullResult()
#define yV2 dup_or_fetch
#define yU2 nominator]
#define yT2 Rehash(false
#define yS2 test_order
#define yR2 }else{x6=new
#define yQ2 nJ3,
#define yP2 .param_count
#define yO2 shift(index)
#define yN2 rulenumber
#define yM2 cLessOrEq,cL
#define yL2 cTanh nR
#define yK2 cSinh nR
#define yJ2 2,7168,
#define yI2 cInv,lT 1,
#define yH2 constraints=
#define yG2 factor_immed
#define yF2 ){Value_t
#define yE2 y01 yF2
#define yD2 changes
#define yC2 );y53
#define yB2 tree l8 iG
#define yA2 DelParam(a);
#define y92 nC1 i3 op1 yV DelParams(
#define y82 nC1 e9 y5 l8
#define y72 e9 leaf2 l8
#define y62 e9 xB3 l8
#define y52 e9 cond l8
#define y42 for(xD3
#define y32 exp_diff
#define y22 ExponentInfo
#define y12 lower_bound(
#define y02 factor
#define xZ2 is_logical
#define xY2 newrel_and
#define xX2 (tU2);
#define xW2 Suboptimal
#define xV2 eW[c tY
#define xU2 res_stackpos
#define xT2 half_pos
#define xS2 ifdata.ofs
#define xR2 >>1)):(
#define xQ2 CodeTreeData
#define xP2 multiply(
#define xO2 var_trees
#define xN2 parent_opcode
#define xM2 log2_exponent
#define xL2 yV swap(tmp);
#define xK2 nG yF2 tmp=
#define xJ2 .max lH3(
#define xI2 dup_fetch_pos
#define xH2 std xT3<bool>
#define xG2 (t71=0;a eM
#define xF2 {cZ start_at;
#define xE2 t71;if(&*start_at){x6=(
#define xD2 (*x6 iI=r.specs;if(r.found){
#define xC2 IsNever cG lC
#define xB2 cH2)&
#define xA2 cSin nR
#define x92 Value_EvenInt
#define x82 MakeFalse,{l5
#define x72 if(xY l8 a)xK
#define x62 ConditionType
#define x52 ;iM.Remember(
#define x42 tree))eK1
#define x32 inline unsigned
#define x22 nA OPCODE
#define x12 }inline
#define x02 n21 Unknown;}
#define nZ2 (IfData&ifdata
#define nY2 (unsigned
#define nX2 iF|nY2)
#define nW2 SpecialOpcode
#define nV2 ByteCodeSynth
#define nU2 AddParamMove(
#define nT2 =i eN2.
#define nS2 fp_max(y3);
#define nR2 assimilated
#define nQ2 denominator
#define nP2 fraction
#define nO2 tree y01 cG
#define nN2 l2 18,2,
#define nM2 .GetDepth()
#define nL2 DUP_BOTH();
#define nK2 xC3 lL
#define nJ2 0x80000000u
#define nI2 -1-offset].
#define nH2 if(synth.Find(
#define nG2 TreeCounts
#define nF2 =GetParam(
#define nE2 bool tB=false;
#define nD2 found_log2
#define nC2 div_params
#define nB2 immed_sum
#define nA2 ByteCode[++IP]
#define n92 nV 2,cPow)
#define n82 );sim.xL nV 2,
#define n72 :sim.Eat(1,
#define n62 OPCODE(opcode)
#define n52 ;sim.Push(
#define n42 FactorStack xE
#define n32 if(&*(*x6 iI){
#define n22 tA2 cG lC
#define n12 282870 nU
#define n02 cNotNot nR
#define lZ2 replacing_slot
#define lY2 RefParams
#define lX2 if_always[
#define lW2 WhatDoWhenCase
#define lV2 exponent_immed
#define lU2 new_base_immed
#define lT2 base_immed
#define lS2 Rehash(iO
#define lR2 2*2)lS 3*
#define lQ2 )tI2
#define lP2 .Rehash()
#define lO2 );m.max.
#define lN2 minimum_need
#define lM2 (size_t a=
#define lL2 ;for lM2
#define lK2 ||op1==
#define lJ2 ):Value(Value::
#define lI2 tV2 nC==
#define lH2 data[a t03
#define lG2 AddCollection(
#define lF2 if(newrel_or==
#define lE2 .AddOperation(
#define lD2 DUP_ONE(apos);
#define lC2 flipped
#define lB2 .UseGetNeeded(
#define lA2 e5 2,131,
#define l92 Immed t13);
#define l82 OptimizedUsing
#define l72 Var_or_Funcno
#define l62 l72;
#define l52 GetParams(
#define l42 crc32_t
#define l32 signed_chain
#define l22 MinusInf
#define l12 n_immeds
#define l02 (tree,std::cout)
#define iZ1 nG2.erase(cs_it);
#define iY1 e42 true;
#define iX1 stack t13
#define iW1 FindClone(xQ
#define iV1 ByteCode[IP]
#define iU1 GetOpcode())
#define iT1 needs_rehash
#define iS1 AnyWhere_Rec
#define iR1 ~unsigned(0)
#define iQ1 41,42,43,44,
#define iP1 yV DelParam(
#define iO1 p1_logical_b
#define iN1 p0_logical_b
#define iM1 p1_logical_a
#define iL1 p0_logical_a
#define iK1 const Rule&rule,
#define iJ1 ,PowiCache&iM,
#define iI1 nT3)
#define iH1 ;p1.lS2 p1
#define iG1 iD i42);
#define iF1 else if(
#define iE1 synth.DoDup(
#define iD1 cache_needed
#define iC1 e5 2,1,e5 2,
#define iB1 treelist
#define iA1 IsDescendantOf(
#define i91 has_bad_balance
#define i81 e33 y02
#define i71 double)eE2
#define i61 {case tA2:
#define i51 fp_abs(i53))
#define i41 fp_abs(min.val)
#define i31 ))break;c53*=
#define i21 {if(GetOpcode()
#define i11 cNEqual
#define i01 },0,0x0},{{
#define tZ1 eZ 2 i01
#define tY1 Oneness_NotOne|
#define tX1 Value_IsInteger
#define tW1 Constness_Const
#define tV1 DumpHashesFrom(
#define tU1 l82(
#define tT1 reltype
#define tS1 const Value_t&i)
#define tR1 const lU1&
#define tQ1 const Value_t&v
#define tP1 SequenceOpcodes
#define tO1 ,eF2
#define tN1 goto fail;}
#define tM1 xC3<
#define tL1 lY2);
#define tK1 )eL mulgroup)
#define tJ1 xC3 set_if<
#define tI1 xC3 lX
#define tH1 TreeCountType xE
#define tG1 >(eF2 1),
#define tF1 yB2);nD lC
#define tE1 eF2 0.0)){xI
#define tD1 n62);
#define tC1 std::cout<<"POP "
#define tB1 stack[iX1-
#define tA1 stack.push_back(
#define t91 synth.PushImmed(
#define t81 MaxChildDepth
#define t71 unsigned a
#define t61 std::pair<It,It>
#define t51 Sign_Negative
#define t41 Value_Logical
#define t31 new_factor_immed
#define t21 if(remaining[a])
#define t11 )nJ3.second
#define t01 tree l8 a)
#define eZ1 occurance_pos
#define eY1 exponent_hash
#define eX1 exponent_list
#define eW1 CollectionSet xE
#define eV1 CollectMulGroup(
#define eU1 source_set
#define eT1 eE2,yB3
#define eS1 *const func)(
#define eR1 tI1 void
#define eQ1 operator
#define eP1 FindAndDup(tree);
#define eO1 back().thenbranch
#define eN1 ParamSpec_Extract
#define eM1 retry_anyparams_3
#define eL1 retry_anyparams_2
#define eK1 goto redo;
#define eJ1 divgroup
#define eI1 ByteCode,size_t&IP,size_t limit,size_t y1
#define eH1 ;cX2.hash2+=
#define eG1 e4(),std xT3<
#define eF1 needlist_cached_t
#define eE1 grammar_rules[*r]
#define eD1 yU3 0x4 iA2
#define eC1 eV2 0x4 iA2
#define eB1 CodeTreeImmed xE(
#define eA1 by_float_exponent
#define e91 fp_equal(eE2
#define e81 new_exp
#define e71 end()&&i->first==
#define e61 return BecomeZero;
#define e51 return BecomeOne;
#define e41 if(lQ t13<=n1)
#define e31 addgroup
#define e21 found_log2by
#define e11 nC==cY3)
#define e01 ParsePowiMuli(
#define cZ1 l72)
#define cY1 eP 529654 nU
#define cX1 ,cNot nR
#define cW1 branch1_backup
#define cV1 branch2_backup
#define cU1 exponent_map
#define cT1 plain_set
#define cS1 rangehalf
#define cR1 LightWeight(
#define cQ1 y33 1
#define cP1 ,iM,eT,tN2
#define cO1 if(value
#define cN1 tI1 c3
#define cM1 .sep_list[
#define cL1 {Value_t eE2=
#define cK1 tI1 static
#define cJ1 :lD yC set(
#define cI1 eF2 0.5))nV 2,
#define cH1 should_regenerate=true;
#define cG1 should_regenerate,
#define cF1 Collection
#define cE1 RelationshipResult
#define cD1 Subdivide_Combine(
#define cC1 long value
#define cB1 )const yQ
#define cA1 rhs cB1 hash1
#define c91 best_sep_factor
#define c81 iF1!c53
#define c71 &&p nP<eF2
#define c61 GetParamCount()
#define c51 needlist_cached
#define c41 eT3 bool pad
#define c31 changed=true;
#define c21 yA2}
#define c11 MakesInteger(
#define c01 ;eE2 lP2;
#define yZ1 const Value_t&value
#define yY1 best_sep_cost
#define yX1 MultiplicationRange
#define yW1 pihalf_limits
#define yV1 p0.max lH3 p0 t9
#define yU1 n_stacked
#define yT1 cX2.hash1
#define yS1 AnyParams_Rec
#define yR1 Become(value l8 0))
#define yQ1 =comp.AddItem(atree
#define yP1 ByteCode[xS2+
#define yO1 ByteCode.push_back(
#define yN1 PositionalParams,0}
#define yM1 always_sincostan
#define yL1 Recheck_RefCount_Div
#define yK1 Recheck_RefCount_Mul
#define yJ1 MultiplyAndMakeLong(
#define yI1 covers_plus1
#define yH1 NeedList
#define yG1 tI1 bool
#define yF1 ;tI1
#define yE1 mulgroup.
#define yD1 .nU2
#define yC1 }y53 case
#define yB1 GetLogicalValue(
#define yA1 if(synth.FindAndDup(
#define y91 SynthesizeParam(
#define y81 grammar_func
#define y71 cOr l3 16,1,
#define y61 252421 nU 24830,
#define y51 Modulo_Radians},
#define y41 .c61
#define y31 yV SetParam(
#define y21 ;synth.StackTopIs(
#define y11 .GetImmed()
#define y01 .IsImmed()
#define xZ1 PositionType
#define xY1 CollectionResult
#define xX1 const_offset
#define xW1 inline TriTruthValue
#define xV1 stacktop_desired
#define xU1 int mStackPtr=0;
#define xT1 SetStackTop(
#define xS1 FPoptimizer_ByteCode
#define xR1 <<tree i22
#define xQ1 public e4,public std xT3<
#define xP1 1)?(poly^(
#define xO1 ,l1 0x0},{{3,
#define xN1 ,l1 tR3
#define xM1 xE())nV 2,cMul);lC
#define xL1 eF2 0)
#define xK1 .min.val
#define xJ1 ,nV2 xE&synth)
#define xI1 ;std::cout<<
#define xH1 )n21 m eL3
#define xG1 xL1)
#define xF1 xK leaf2 l8
#define xE1 cond_type
#define xD1 fphash_value_t
#define xC1 Recheck_RefCount_RDiv
#define xB1 tmp yD1 tree
#define xA1 cMul);tmp.nI 0 tE3.
#define x91 (lR,yB2));nD
#define x81 tree.DelParam(a
#define x71 t01 y11;
#define x61 t01 y01)
#define x51 SetParams(l52)
#define x41 fPExponentIsTooLarge(
#define x31 CollectMulGroup_Item(
#define x21 nL xT1 xU-1);
#define x11 covers_full_cycle
#define x01 AssembleSequence(
#define nZ1 252180 nU 281854,
#define nY1 <<std::dec<<")";}
#define nX1 },{l5::MakeNotP1,l5::
#define nW1 },{l5::MakeNotP0,l5::
#define nV1 },{l5::MakeNotNotP1,l5::
#define nU1 :eY3 xK1
#define nT1 std::pair<T1,T2>&
#define nS1 tM1 xD3
#define nR1 has_good_balance_found
#define nQ1 n_occurrences
#define nP1 found_log2_on_exponent
#define nO1 covers_minus1
#define nN1 needs_resynth
#define nM1 immed_product
#define nL1 ,2,1 xJ if(found[data.
#define nK1 nC3 bitmask&
#define nJ1 Sign_Positive
#define nI1 {DataP slot_holder(xZ[
#define nH1 ::MakeTrue
#define nG1 SetParamMove(
#define nF1 ::pair<Value_t,yB3>
#define nE1 lQ3 FPoptimizer_Optimize
#define nD1 CodeTreeImmed(eF2
#define nC1 changed_if
#define nB1 n_as_tanh_param
#define nA1 0 e32 DelParam(1);
#define n91 (p1 y11
#define n81 opposite=
#define n71 xD1(
#define n61 ByteCode t13
#define n51 MatchResultType
#define n41 needs_sincos
#define n31 resulting_exponent
#define n21 ;return
#define n11 n21 xW2;}
#define n01 Unknown:yT3;}
#define lZ1 GetParam(a)
#define lY1 inverse_nominator]
#define lX1 cMul l3 0,1,
#define lW1 AddFunctionOpcode(
#define lV1 xE(rule.repl_param_list,
#define lU1 CodeTree
#define lT1 lU1 xE
#define lS1 void FunctionParserBase
#define lR1 o<<"("<<std::hex<<data.
#define lQ1 (val);else*this=model;}
#define lP1 IfBalanceGood(
#define lO1 n_as_tan_param
#define lN1 changed_exponent
#define lM1 inverse_denominator
#define lL1 retry_positionalparams_2
#define lK1 unsigned index
#define lJ1 situation_flags&
#define lI1 518 nU 400412,
#define lH1 data.subfunc_opcode
#define lG1 },{l5::MakeNotNotP0,l5::
#define lF1 c02);
#define lE1 recursioncount
#define lD1 PlanNtimesCache(
#define lC1 FPoptimizer_Grammar
#define lB1 AddOperation(cInv,1,1 xJ}
#define lA1 (nS yP2
#define l91 GetPositivityInfo cE3
#define l81 ParamSpec_SubFunctionData
#define l71 tK2){synth.yT
#define l61 PositionalParams_Rec
#define l51 y41;a-->0;)if(
#define l41 tV2 y01
#define l31 lJ 2 i01 1,
#define l21 DumpTreeWithIndent(*this);
#define l11 switch(type){case cond_or:
#define l01 tM1 unsigned Compare>
#define iZ edited_powgroup
#define iY lM2 c61;a
#define iX has_unknown_max
#define iW has_unknown_min
#define iV nO eF2-1)yI3
#define iU static const range xE
#define iT if(keep_powi
#define iS synthed_tree
#define iR 7168 nU 401798,
#define iQ SelectedParams,0 i01
#define iP collections
#define iO yV nU2
#define iN lS2 r);}
#define iM cache
#define iL ]);y43
#define iK ,ByteCode,IP,limit,y1,stack);
#define iJ by_exponent
#define iI )[a].start_at
#define iH ;if(t82 tM3 SaveMatchedParamIndex(
#define iG 1)y11
#define iF ;yO1 nJ2
#define iE )iF);
#define iD .SetOpcode(
#define iC mulgroup;mulgroup iD
#define iB goto ReplaceTreeWithOne;case
#define iA c83 165888 nU
#define i9 eA1.data
#define i8 lB3 xQ2(
#define i7 needs_sinhcosh
#define i6 cAdd l3 0,
#define i5 ,cIf,l0 3,
#define i4 .min.known
#define i3 .Rehash(yV SetOpcode(
#define i2 ;lT1
#define i1 i92 eF2
#define i0 tI1 nB
#define tZ MakeFalse,l5::
#define tY ].relationship
#define tX p0 i4&&p0 xK1>=eF2 0.0))
#define tW p0=CalculateResultBoundaries(
#define tV 408964 nU 24963,
#define tU 528504 nU 24713,
#define tT matched_params
#define tS [n1 eZ3=true;lQ[n1 t03
#define tR lC1::Grammar*
#define tQ powgroup l8
#define tP tK2&&found[data.
#define tO }},{ProduceNewTree,2,1,
#define tN nD1(
#define tM has_mulgroups_remaining
#define tL valueType
#define tK MatchInfo xE&
#define tJ Rehash();iS2.push_back(
#define tI int_exponent_t
#define tH best_factor
#define tG RootPowerTable xE::RootPowers[
#define tF MatchPositionSpec_AnyParams xE
#define tE lQ3 FPoptimizer_CodeTree
#define tD n_as_sinh_param
#define tC n_as_cosh_param
#define tB is_signed
#define tA result_positivity
#define t9 nP<=fp_const_negativezero xE())
#define t8 biggest_minimum
#define t7 const l81
#define t6 124024 nU 139399,
#define t5 142456 nU 141449,
#define t4 lA 0x4},{{
#define t3 !=Unchanged)if(TestCase(
#define t2 cond_tree
#define t1 else_tree
#define t0 then_tree
#define eZ x7 AnyParams,
#define eY ,l4 2,1,
#define eX CalculateResultBoundaries(tree l8
#define eW relationships
#define eV tree l8 0)
#define eU lT1&tree)
#define eT sequencing
#define eS string t53(
#define eR tree y41
#define eQ );bool needs_cow=GetRefCount()>1;
#define eP ,nN2
#define eO std xT3<lT1>
#define eN if_stack
#define eM <nS yP2;++a)
#define eL ;nU2
#define eK (l52));yE1 Rehash();
#define eJ .max.set(fp_ceil xH1
#define eI n_as_sin_param
#define eH n_as_cos_param
#define eG PowiResolver::
#define eF cIf,tS3
#define eE PACKED_GRAMMAR_ATTRIBUTE;
#define eD .BalanceGood
#define eC nU2 cA2
#define eB back().endif_location
#define eA xD1 key
#define e9 .AddParam(
#define e8 nU2 mul);
#define e7 :goto ReplaceTreeWithZero;case
#define e6 ;p2.lS2 p2 yV SetOpcode(lG3 nC);eK1}
#define e5 130,1,
#define e4 MatchPositionSpecBase
#define e3 lB3 lU1(
#define e2 smallest_maximum
#define e1 ReplaceTreeWithParam0;
#define e0 factor_needs_rehashing
#define cZ MatchPositionSpecBaseP
#define cY xD3 tH1::yE3
#define cX fp_cosh cD3);m nP=fp_cosh(m nP);
#define cW {AdoptChildrenWithSameOpcode(tree);
#define cV eN1 xE(nS.param_list,
#define cU 243,244,245,246,249,250,251,253,255,256,257,258,259}};}
#define cT ];};extern"C"{
#define cS 79,122,123,160,161,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,
#define cR 27,28,29,30,31,32,33,35,36,
#define cQ const ParamSpec_SubFunction
#define cP const ParamSpec_ParamHolder
#define cO otherhalf
#define cN StackState
#define cM 0x4 nN
#define cL l2 16,2,
#define cK const SequenceOpCode xE
#define cJ MatchPositionSpec_PositionalParams xE
#define cI cI2&tree,std::ostream&o
#define cH nT3=false;
#define cG )return false;
#define cF eF2 1.5)*fp_const_pi xE()
#define cE !=Unchanged)return lX2
#define cD CalculatePowiFactorCost(
#define cC ImmedHashGenerator
#define cB ::map<fphash_t,std::set<std cG3> >
#define cA nU2 comp.cT1[a].value);
#define c9 T1,xD3 T2>inline iY2()(
#define c8 has_nonlogical_values
#define c7 from_logical_context)
#define c6 AnyParams,0}},{ProduceNewTree,
#define c5 for lM2 xY y41;a-->0;)
#define c4 POWI_CACHE_SIZE
#define c3 static inline lT1
#define c2 ++IP;tK3 if(iV1==cU3.
#define c1 .FoundChild
#define c0 BalanceResultType
#define yZ paramholder_matches
#define yY {lT1 tmp;tmp iD
#define yX ComparisonSetBase::
#define yW lM2 eR;a-->0;)
#define yV );tree.
#define yU .nI 0 tE3 iD cInv);tmp yD1 tmp2)n21
#define yT DoDup(found[data.
#define yS xF3(0),Opcode(
#define yR );void lW1 unsigned eT3 cL2<
#define yQ {return
#define yP const yQ data->
#define yO +=fp_const_twopi xE();
#define yN y43 GetOpcode(),
#define yM for lM2 0;a<c61;++a){if(
#define yL static void nL3 nA fphash_t&cX2,
#define yK MatchPositionSpec_AnyWhere
#define yJ iQ3 tree)xI1"\n";
#define yI );pow iD cLog yV SetOpcode(cMul);
#define yH if eP2.match_type==
#define yG void OutFloatHex(std::ostream&o,
#define yF nO3 1,cInv yC2 xF-1)n92;lC
#define yE AddParam(CodeTreeImmed(
#define yD lT1 tmp,tmp2;tmp2 iD
#define yC m.min.
#define yB cGreaterOrEq,
#define yA ,xD3 lT1::
#define y9 xE model=cS1 xE()){if(known
#define y8 AssembleSequence_Subdivide(
#define y7 ]=nJ2|unsigned(
#define y6 )xI1 std::endl;DumpHashes(
#define y5 branch2
#define y4 unsigned c tJ3 short l[
#define y3 fp_sin(min),fp_sin(max))
#define y2 fp_const_twopi xE()yG3
#define y1 factor_stack_base
#define y0 ,cPow,lA
#define xZ data->eD3
#define xY branch1
#define xX );tK3 if e72 y11==eF2
#define xW (nI3 r=range.first;r!=range.iD3++r){
#define xV {nG2.erase(i);tK3
#define xU StackTop
#define xT FPOPT_autoptr
#define xS +=c53 n21 c53;}tI1 inline Value_t
#define xR int_exponent
#define xQ newnode
#define xP has_highlevel_opcodes
#define xO if(needs_cow){lF1 goto
#define xN },{l5::Unchanged,l5::Never},{l5::Unchanged,l5::Never}}
#define xM size_t a=0;a<yA3 y41;++a)if(
#define xL SwapLastTwoInStack()
#define xK .IsIdenticalTo(
#define xJ )y21*this)n21;}
#define xI tree.ReplaceWithImmed(
#define xH lM2 0;a<eR;++a)
#define xG GetStackTop()-
#define xF sim.AddConst(
#define xE <Value_t>
#define xD ,lJ 1 i01
#define xC cI=std::cout
#define xB best_selected_sep
#define xA ->Recalculate_Hash_NoRecursion();}
#define x9 c61;++a)if(ApplyGrammar(iH2,t01,
#define x8 i12 2,1,
#define x7 ,cAdd,
#define x6 position
#define x5 for xH{range xE
#define x4 SetParam(0,lG3 l8 0))i2 p1;p1 iD
#define x3 std xT3<lU1>
#define x2 TestImmedConstraints(param.constraints,tree)cG
#define x1 paramholder_index
#define x0 return true;case
#define nZ occurance_counts
#define nY -->0;){cI2&powgroup=lZ1;if(powgroup
#define nX )){tree.FixIncompleteHashes();}
#define nW Value_t>p eN3 p.
#define nV ;sim.Eat(
#define nU ,{2,
#define nT const FPoptimizer_CodeTree::lT1&tree
#define nS model_tree
#define nR ,l0 1,
#define nQ ,cAdd l3 2,1,
#define nP .i53
#define nO return range xE(
#define nN },{{2,
#define nM )){lT1
#define nL ){using lQ3 FUNCTIONPARSERTYPES;
#define nK AnyParams,1 i01
#define nJ eO&lY2
#define nI AddParam(tree l8
#define nH ConstantFolding_LogicCommon(tree,yX
#define nG if(eV y01
#define nF nS1 Ref>inline void xT<Ref>::
#define nE ):data(new xQ2 xE(
#define nD goto do_return;}
#define nC .GetOpcode()
#define nB xQ2 xE::xQ2(
#define nA FUNCTIONPARSERTYPES::
#define n9 b;}};tM1>eI2 Comp<nA
#define n8 l72(),eD3(),Hash(),Depth(1),tU1 0){}
#define n7 SynthesizeByteCode(tN2
#define n6 while(ApplyGrammar((xB2
#define n5 cF3 eV)==tA2)lF3
#define n4 ;tree yD1 nC1 lQ2
#define n3 tJ1 cGreater>(eF2
#define n2 DumpParams xE eP2.param_list,param.data yP2,o);
#define n1 restholder_index
#define n0 lT1 eE2;eE2 t83 eE2 e9
#define lZ lR yG3 fp_nequal(tmp,xG1){xI eF2 1)/tmp);nD}lC
#define lY :if(ParamComparer xE()(eD3[1],eD3[0])l34 swap(eD3[0],eD3[1]);Opcode=
#define lX <xD3 Value_t>
#define lW xE tmp;tmp iD cPow);tmp.nI 0 tE3.yE eF2
#define lV tW1,0x0},
#define lU nU2 pow l8 1));pow.DelParam(1);pow.Rehash(yV nG1 0,pow);goto NowWeAreMulGroup;}
#define lT GroupFunction,0},lV{{
#define lS tO1 1)/eF2
#define lR eV y11
#define lQ restholder_matches
#define lP yT1|=key;xD1 crc=(key>>10)|(key<<(64-10))eH1((~n71 crc))*3)^1234567;}};
#define lO nC1;nC1 iG1 nC1 yD1 eV);nC1 e9 xY l8
#define lN tI1 lT1::lU1(
#define lM tree.SetParam(0,eV l8 0)y31 1,CodeTreeImmed(
#define lL lX void nV2 xE::lW1 unsigned eT3 cL2<
#define lK cMul iX3
#define lJ cMul,AnyParams,
#define lI (eV y01&&l41){xI
#define lH CalculateResultBoundaries(tmp)eL3
#define lG :eG3=comp.AddRelationship(atree l8 0),atree l8 1),yX
#define lF cPow,l0 2
#define lE xD3 Value_t>inline iY2()(const Value_t&xG3 Value_t&b)yQ a
#define lD {range xE m=CalculateResultBoundaries(eV);
#define lC break;case
#define lB eR1 lT1::
#define lA yN1,0,
#define l9 l1 0x0 nN
#define l8 .GetParam(
#define l7 i2 nC1;nC1 iG1 nC1 e83 tree.l52));nC1 i3
#define l6 SelectedParams,0},0,0x0 nN
#define l5 RangeComparisonData
#define l4 yN1},{ProduceNewTree,
#define l3 ,AnyParams,0}},{ReplaceParams,
#define l2 yN1},{ReplaceParams,
#define l1 cMul,SelectedParams,0},0,
#define l0 lA 0x0},{{
#ifdef _MSC_VER
typedef
unsigned
int
l42;
#else
#include <stdint.h>
typedef
uint_least32_t
l42;
#endif
lQ3
crc32{enum{startvalue=0xFFFFFFFFUL,poly=0xEDB88320UL}
;tM1
l42
crc>eI2
b8{enum{b1=(crc&xP1
crc
xR2
crc>>1),b2=(b1&xP1
b1
xR2
b1>>1),b3=(b2&xP1
b2
xR2
b2>>1),b4=(b3&xP1
b3
xR2
b3>>1),b5=(b4&xP1
b4
xR2
b4>>1),b6=(b5&xP1
b5
xR2
b5>>1),b7=(b6&xP1
b6
xR2
b6>>1),res=(b7&xP1
b7
xR2
b7>>1)}
;}
;inline
l42
update(l42
crc,unsigned
b){
#define B4(n) b8<n>lP3 n+1>lP3 n+2>lP3 n+3>::res
#define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
static
const
l42
table[256]={R(0x00),R(0x10),R(0x20),R(0x30),R(0x40),R(0x50),R(0x60),R(0x70),R(0x80),R(0x90),R(0xA0),R(0xB0),R(0xC0),R(0xD0),R(0xE0),R(0xF0)}
;
#undef R
#undef B4
return((crc>>8))^table[(crc^b)&0xFF];x12
l42
calc_upd(l42
c,const
unsigned
char*buf,size_t
size){l42
value=c;for(size_t
p=0;p<size;++p)value=update(value,buf[p])n21
value;x12
l42
calc(const
unsigned
char*buf,size_t
size)yQ
calc_upd(startvalue,buf,size);}
}
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH
nS1
Ref>class
xT{eS3
xT():p(0){}
xT(Ref*b):p(b){xE3}
xT(const
xT&b):p(b.p){xE3
x12
Ref&eQ1*(cB1*p;x12
Ref*eQ1->(cB1
p;}
xT&eQ1=(Ref*b){Set(b)n21*this;}
xT&eQ1=(const
xT&b){Set(b.p)n21*this;}
#ifdef FP_SUPPORT_CXX11_MOVE
xT(xT&&b):p(b.p){b.p=0;}
xT&eQ1=(xT&&b){if(p!=b.p){eH2;p=b.p;b.p=0;eX3*this;}
#endif
~xT(){eH2
lT3
UnsafeSetP(Ref*newp){p=newp
lT3
swap(xT<Ref>&b){Ref*tmp=p;p=b.p;b.p=tmp;}
private:inline
static
void
Have(Ref*p2);inline
void
eH2;inline
void
xE3
inline
void
Set(Ref*p2);private:Ref*p;}
;nF
eH2{if(!p)return;p->xF3-=1;if(!p->xF3)delete
p;}
nF
Have(Ref*p2){if(p2)++(p2->xF3);}
nF
Birth(){Have(p);}
nF
Set(Ref*p2){Have(p2);eH2;p=p2;}
#endif
#include <utility>
eI2
Compare2ndRev{nS1
T>inline
iY2()(const
T&xG3
T&b
cB1
a.second>b.iD3}
}
;eI2
Compare1st{nS1
c9
const
nT1
xG3
nT1
b
cB1
a.first<b.first;}
nS1
c9
const
nT1
a,T1
b
cB1
a.first<b;}
nS1
c9
T1
xG3
nT1
b
cB1
a<b.first;}
}
;
#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH
#ifdef _MSC_VER
typedef
unsigned
long
long
xD1;
#define FPHASH_CONST(x) x##ULL
#else
#include <stdint.h>
typedef
uint_fast64_t
xD1;
#define FPHASH_CONST(x) x##ULL
#endif
lQ3
FUNCTIONPARSERTYPES{eI2
fphash_t{xD1
hash1,hash2;fphash_t():hash1(0),hash2(0){}
fphash_t(const
xD1&xG3
xD1&b):hash1(a),hash2(b){}
iY2==(const
fphash_t&cA1==eJ2&&hash2==eK2
iY2!=(const
fphash_t&cA1!=eJ2||hash2!=eK2
iY2<(const
fphash_t&cA1!=eJ2?hash1<eJ2:hash2<eK2}
;}
#endif
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
lQ3
lC1{eI2
Grammar;}
lQ3
xS1{tI1
class
nV2;}
tE{tI1
class
lU1
yF1
eI2
xQ2
yF1
class
lU1{typedef
xT<xQ2
xE>DataP;DataP
data;eS3
lU1();~lU1();eI2
OpcodeTag{}
;e3
x22
o,OpcodeTag);eI2
FuncOpcodeTag{}
;e3
x22
o,unsigned
f,FuncOpcodeTag);eI2
xH3{}
;e3
tQ1,xH3);
#ifdef FP_SUPPORT_CXX11_MOVE
e3
Value_t&&v,xH3);
#endif
eI2
VarTag{}
;e3
unsigned
varno,VarTag);eI2
CloneTag{}
;e3
tR1
b,CloneTag);void
GenerateFrom(const
xD3
FunctionParserBase
xE::Data&data,bool
keep_powi=false);void
GenerateFrom(const
xD3
FunctionParserBase
xE::Data&data,const
x3&xO2,bool
keep_powi=false);void
SynthesizeByteCode(std
xT3<unsigned>&eB3,std
xT3
xE&immed,size_t&stacktop_max);void
SynthesizeByteCode(xS1::nV2
xE&synth,bool
MustPopTemps=true)const;size_t
SynthCommonSubExpressions(xS1::nV2
xE&synth)const;void
SetParams(const
x3&xI3
SetParamsMove(x3&tL1
lU1
GetUniqueRef();
#ifdef FP_SUPPORT_CXX11_MOVE
void
SetParams(x3&&tL1
#endif
void
SetParam(size_t
which,tR1
b);void
nG1
size_t
which,lU1&b);void
AddParam(tR1
param);void
nU2
lU1&param);void
AddParams(const
x3&xI3
AddParamsMove(x3&xI3
AddParamsMove(x3&lY2,size_t
lZ2);void
DelParam(size_t
index);void
DelParams();void
Become(tR1
b);inline
size_t
c61
const
yQ
l52)t13;x12
lU1&GetParam(cO3)yQ
l52)[n];x12
tR1
GetParam(cO3
cB1
l52)[n];x12
void
SetOpcode(x22
o)iE3
Opcode=o;x12
x22
GetOpcode()yP
Opcode;x12
nA
fphash_t
GetHash()yP
Hash;x12
const
x3&l52
cB1
xZ;x12
x3&l52)yQ
xZ;x12
size_t
GetDepth()yP
Depth;x12
const
Value_t&GetImmed()yP
Value;x12
unsigned
GetVar()yP
l62
x12
unsigned
GetFuncNo()yP
l62
x12
bool
c92
cB1
GetOpcode()!=nA
cNop;x12
bool
IsImmed(cB1
GetOpcode()==nA
cImmed;x12
bool
IsVar(cB1
GetOpcode()==nA
l93;x12
unsigned
GetRefCount()yP
xF3
lT3
ReplaceWithImmed(tS1;void
Rehash(bool
constantfolding=true);void
Sort();inline
void
Mark_Incompletely_Hashed()iE3
Depth=0;x12
bool
Is_Incompletely_Hashed()yP
Depth==0;x12
xJ3
GetOptimizedUsing()yP
l82;x12
void
SetOptimizedUsing(xJ3
g)iE3
l82=g;}
bool
RecreateInversionsAndNegations(bool
prefer_base2=false);void
FixIncompleteHashes();void
swap(lU1&b){data.swap(b.data);}
bool
IsIdenticalTo(tR1
b)const;void
lF1}
yF1
eI2
xQ2{int
xF3;x22
Opcode
yX3
Value
tJ3
l62
eO
eD3;nA
fphash_t
Hash;size_t
Depth;xJ3
l82;xQ2();xQ2(const
xQ2&b);i8
x22
o);i8
x22
o,unsigned
f);i8
tS1;
#ifdef FP_SUPPORT_CXX11_MOVE
i8
Value_t&&i);xQ2(xQ2&&b);
#endif
bool
IsIdenticalTo(const
xQ2&b)const;void
Sort();void
Recalculate_Hash_NoRecursion();private:void
eQ1=(const
xQ2&b);}
yF1
c3
CodeTreeImmed(tS1
yQ
lT1(i
yA
xH3());}
#ifdef FP_SUPPORT_CXX11_MOVE
cN1
CodeTreeImmed
tZ3&&i)yQ
lT1(i82
i)yA
xH3());}
#endif
cN1
CodeTreeOp(x22
opcode)yQ
lT1(opcode
yA
OpcodeTag());}
cN1
CodeTreeFuncOp(x22
eT3
unsigned
f)yQ
lT1(eT3
f
yA
FuncOpcodeTag());}
cN1
CodeTreeVar
nY2
varno)yQ
lT1(varno
yA
VarTag());}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
eR1
DumpHashes(xC)xM3
DumpTree(xC)xM3
DumpTreeWithIndent(xC,const
std
cG3&indent="\\"
);
#endif
}
#endif
#endif
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH
#include <iostream>
tE{tI1
class
lU1;}
lQ3
lC1{enum
ImmedConstraint_Value{ValueMask=0x07,Value_AnyNum=0x0,x92=0x1,Value_OddInt=0x2,tX1=0x3,Value_NonInteger=0x4,t41=0x5
x73
ImmedConstraint_Sign{SignMask=0x18,Sign_AnySign=0x00,nJ1=0x08,t51=0x10,Sign_NoIdea=0x18
x73
ImmedConstraint_Oneness{OnenessMask=0x60,Oneness_Any=0x00,Oneness_One=0x20,Oneness_NotOne=0x40
x73
ImmedConstraint_Constness{ConstnessMask=0x180,Constness_Any=0x00,tW1=0x80,Constness_NotConst=0x100
x73
Modulo_Mode{Modulo_None=0,Modulo_Radians=1
x73
Situation_Flags{LogicalContextOnly=0x01,NotForIntegers=0x02,OnlyForIntegers=0x04,OnlyForComplex=0x08,NotForComplex=0x10
x73
nW2{NumConstant,ParamHolder,SubFunction
x73
ParamMatchingType{PositionalParams,SelectedParams,AnyParams,GroupFunction
x73
RuleType{ProduceNewTree,ReplaceParams}
;
#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif
typedef
std::pair<nW2,cH2>eL2
yF1
eL2
eN1
nY2
paramlist,lK1)yF1
bool
ParamSpec_Compare(cH2
a,cH2
b,nW2
type)tJ3
ParamSpec_GetDepCode(const
eL2&b);eI2
ParamSpec_ParamHolder{lK1:8
tJ3
constraints:9
tJ3
depcode:15;}
eE
tI1
eI2
ParamSpec_NumConstant{Value_t
constvalue
tJ3
modulo;}
;eI2
l81{unsigned
param_count:2
tJ3
param_list:30;x22
subfunc_opcode:8;ParamMatchingType
match_type:3
tJ3
n1:5;}
eE
eI2
ParamSpec_SubFunction{l81
data
tJ3
constraints:9
tJ3
depcode:7;}
eE
eI2
Rule{RuleType
ruletype:2
tJ3
situation_flags:5
tJ3
repl_param_count:2+9
tJ3
repl_param_list:30;l81
match_tree;}
eE
eI2
Grammar{unsigned
rule_count
tJ3
short
rule_list[999
cT
extern
const
Rule
grammar_rules[];}
eR1
DumpParam(const
eL2&p,std::ostream&o=std::cout)xM3
DumpParams
nY2
paramlist,unsigned
count,std::ostream&o=std::cout);}
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define CONSTANT_POS_INF HUGE_VAL
#define CONSTANT_NEG_INF (-HUGE_VAL)
lQ3
FUNCTIONPARSERTYPES{tI1
inline
Value_t
fp_const_pihalf()yQ
fp_const_pi
xE()*eF2
0.5);}
tI1
inline
Value_t
fp_const_twopi(cK2
fp_const_pi
xE());lR3
fp_const_twoe(cK2
fp_const_e
xE());lR3
fp_const_twoeinv(cK2
fp_const_einv
xE());lR3
fp_const_negativezero()yQ-Epsilon
xE::value;}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#include <iostream>
nE1{using
lQ3
lC1;using
tE;using
lQ3
FUNCTIONPARSERTYPES
yF1
class
MatchInfo{eS3
std
xT3<std::pair<bool,eO> >lQ;eO
yZ;std
xT3<unsigned>tT;eS3
MatchInfo():lQ(),yZ(),tT(){}
eS3
bool
SaveOrTestRestHolder
nY2
n1,const
eO&iB1){e41{lQ
xK3
n1+1);lQ
tS=iB1
tI2
if(lQ[n1
eZ3==false){lQ
tS=iB1
tI2
const
eO&found=lQ[n1
t03;if(iB1
t13!=found
t13
cG
for
lM2
0;a<iB1
tZ2
a)if(!iB1[a]xK
found[a])cG
return
true
lT3
SaveRestHolder
nY2
n1,eO&iB1){e41
lQ
xK3
n1+1);lQ
tS.swap(iB1);}
bool
SaveOrTestParamHolder
nY2
x1,cI2&xL3){if(yZ
t13<=x1){yZ.xP3
x1+1);yZ
xK3
x1);yZ.push_back(xL3
lQ2
if(!yZ[x1].c92)){yZ[x1]=xL3
tI2
return
xL3
xK
yZ[x1]lS3
SaveMatchedParamIndex(lK1){tT.push_back(index);}
cI2&GetParamHolderValueIfFound
nY2
x1)const{static
cI2
dummytree;if(yZ
t13<=x1)return
dummytree
n21
yZ[x1];}
cI2&GetParamHolderValue
nY2
x1
cB1
yZ[x1];}
bool
HasRestHolder
nY2
n1
cB1
lQ
t13>n1&&lQ[n1
eZ3==true;}
const
eO&GetRestHolderValues
nY2
n1)const{static
const
eO
empty_result;e41
return
empty_result
n21
lQ[n1
t03;}
const
std
xT3<unsigned>&GetMatchedParamIndexes(cB1
tT
lT3
swap(tK
b){lQ.swap(b.lQ);yZ.swap(b.yZ);tT.swap(b.tT);}
tK
eQ1=(const
tK
b){lQ=b.lQ;yZ=b.yZ;tT=b.tT
n21*this;}
}
;class
e4
i03
xT<e4>cZ;class
e4{eS3
int
xF3;eS3
e4():xF3(0){}
virtual~e4(){}
}
;eI2
n51{bool
found;cZ
specs;n51(bool
f):found(f),specs(){}
n51(bool
f
t23
s):found(f),specs(s){}
}
xM3
SynthesizeRule(iK1
lT1&tree,l74)yF1
n51
TestParam(const
eL2&yQ2
cI2&tree
t23
start_at,l74)yF1
n51
TestParams(t7&nS,cI2&tree
t23
start_at,l74,bool
t82
yF1
bool
ApplyGrammar(const
Grammar&iH2,FPoptimizer_CodeTree::lT1&tree,bool
from_logical_context=false)xM3
ApplyGrammars(FPoptimizer_CodeTree::eU
yF1
bool
IsLogisticallyPlausibleParamsMatch(t7&tD3,const
eU;}
lQ3
lC1{eR1
DumpMatch(iK1
nT,const
FPoptimizer_Optimize::l74,bool
DidMatch,std::ostream&o=std::cout)xM3
DumpMatch(iK1
nT,const
FPoptimizer_Optimize::l74,xN3
tW3,std::ostream&o=std::cout);}
#endif
#include <string>
xO3
lC1::nW2
c41=false);xO3
x22
c41=false);
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
using
lQ3
lC1;using
lQ3
FUNCTIONPARSERTYPES;xO3
lC1::nW2
c41){
#if 1
xN3
p=0;switch(opcode){case
t63
p="NumConstant"
;lC
ParamHolder:p="ParamHolder"
;lC
SubFunction:p="SubFunction"
;y53
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str()t13<12)tmp<<' '
n21
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str()t13<5)tmp<<' '
n21
tmp.str();
#endif
}
xO3
x22
c41){
#if 1
xN3
p=0;switch(opcode){case
cAbs:p="cAbs"
;lC
cAcos:p="cAcos"
;lC
cAcosh:p="cAcosh"
;lC
cArg:p="cArg"
;lC
cAsin:p="cAsin"
;lC
cAsinh:p="cAsinh"
;lC
cAtan:p="cAtan"
;lC
cAtan2:p="cAtan2"
;lC
cAtanh:p="cAtanh"
;lC
cCbrt:p="cCbrt"
;lC
cCeil:p="cCeil"
;lC
cConj:p="cConj"
;lC
cCos:p="cCos"
;lC
cCosh:p="cCosh"
;lC
cCot:p="cCot"
;lC
cCsc:p="cCsc"
;lC
cExp:p="cExp"
;lC
cExp2:p="cExp2"
;lC
cFloor:p="cFloor"
;lC
cHypot:p="cHypot"
;lC
cIf:p="cIf"
;lC
cImag:p="cImag"
;lC
cInt:p="cInt"
;lC
cLog:p="cLog"
;lC
cLog2:p="cLog2"
;lC
cLog10:p="cLog10"
;lC
cMax:p="cMax"
;lC
cMin:p="cMin"
;lC
cPolar:p="cPolar"
;lC
cPow:p="cPow"
;lC
cReal:p="cReal"
;lC
cSec:p="cSec"
;lC
cSin:p="cSin"
;lC
cSinh:p="cSinh"
;lC
cSqrt:p="cSqrt"
;lC
cTan:p="cTan"
;lC
cTanh:p="cTanh"
;lC
cTrunc:p="cTrunc"
;lC
cImmed:p="cImmed"
;lC
cJump:p="cJump"
;lC
cNeg:p="cNeg"
;lC
cAdd:p="cAdd"
;lC
cSub:p="cSub"
;lC
cMul:p="cMul"
;lC
cDiv:p="cDiv"
;lC
cMod:p="cMod"
;lC
cEqual:p="cEqual"
;lC
i11:p="cNEqual"
;lC
cLess:p="cLess"
;lC
cLessOrEq:p="cLessOrEq"
;lC
cGreater:p="cGreater"
;lC
cGreaterOrEq:p="cGreaterOrEq"
;lC
cNot:p="cNot"
;lC
cAnd:p="cAnd"
;lC
cOr:p="cOr"
;lC
cDeg:p="cDeg"
;lC
cRad:p="cRad"
;lC
cFCall:p="cFCall"
;lC
cPCall:p="cPCall"
;break;
#ifdef FP_SUPPORT_OPTIMIZER
case
cFetch:p="cFetch"
;lC
cPopNMov:p="cPopNMov"
;lC
e63:p="cLog2by"
;lC
cNop:p="cNop"
;break;
#endif
case
cSinCos:p="cSinCos"
;lC
cSinhCosh:p="cSinhCosh"
;lC
cY3:p="cAbsNot"
;lC
cAbsNotNot:p="cAbsNotNot"
;lC
cAbsAnd:p="cAbsAnd"
;lC
cAbsOr:p="cAbsOr"
;lC
cAbsIf:p="cAbsIf"
;lC
cDup:p="cDup"
;lC
cInv:p="cInv"
;lC
cSqr:p="cSqr"
;lC
cRDiv:p="cRDiv"
;lC
cRSub:p="cRSub"
;lC
cNotNot:p="cNotNot"
;lC
cRSqrt:p="cRSqrt"
;lC
l93:p="VarBegin"
;y53
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str()t13<12)tmp<<' '
n21
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str()t13<5)tmp<<' '
n21
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
;lQ3
xS1{tI1
class
nV2{eS3
nV2():ByteCode(),Immed(),cN(),xU(0),StackMax(0){ByteCode.xP3
64);Immed.xP3
8);cN.xP3
16
lS3
Pull(std
xT3<unsigned>&bc,std
xT3
xE&imm,size_t&StackTop_max){for(t71=0;a<n61;++a){ByteCode[a]&=~nJ2;}
ByteCode.swap(bc);Immed.swap(imm);StackTop_max=StackMax;}
size_t
GetByteCodeSize(cB1
n61;}
size_t
GetStackTop(cB1
xU
lT3
PushVar
nY2
varno){yO1
varno);tJ2}
void
PushImmed
tZ3
immed
nL
yO1
cImmed);Immed.push_back(immed);tJ2}
void
StackTopIs(nT,int
offset=0){if((int)xU>offset){cN[xU
nI2
first=true;cN[xU
nI2
second=tree;}
}
bool
IsStackTop(nT,int
offset=0
cB1(int)xU>offset&&cN[xU
nI2
first&&cN[xU
nI2
second
xK
tree);x12
void
EatNParams
nY2
eat_count){xU-=eat_count
lT3
ProducedNParams
nY2
produce_count){xT1
xU+produce_count
lS3
DoPopNMov(size_t
eM2,size_t
srcpos
nL
yO1
cPopNMov)nX2
eM2)nX2
srcpos);xT1
srcpos+1);cN[eM2]=cN[srcpos];xT1
eM2+1
lS3
DoDup(size_t
xQ3
nL
if(xQ3==xU-1){yO1
cDup);}
else{yO1
cFetch)nX2
xQ3);}
tJ2
cN[xU-1]=cN[xQ3];}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tM1
int>void
Dump(l34
ostream&o=std::cout;o<<"Stack state now("
<<xU<<"):\n"
lL2
0;a<xU;++a){o<<a<<": "
;if(cN[a
eZ3){nT=cN[a
t03;o<<'['<<std::hex<<(void*)(&tree.l52))<<std::dec<<','<<tree.GetRefCount()<<']'
iQ3
tree,o);}
else
o<<"?"
;o<<"\n"
;}
o<<std::flush;}
#endif
size_t
xR3
nT)const{for
lM2
xU;a-->0;)if(cN[a
eZ3&&cN[a
t03
xK
tree
e42
a
n21
tK2;}
bool
Find(nT
cB1
xR3
tree)!=tK2;}
bool
FindAndDup(nT)iA3
pos=xR3
tree
yG3
pos!=tK2){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<iZ3"duplicate at ["
<<pos<<"]: "
iQ3
tree)xI1" -- issuing cDup or cFetch\n"
;
#endif
DoDup(pos
lQ2
return
tI3
eI2
IfData
iA3
ofs;}
;void
SynthIfStep1
nZ2,x22
op
x21
xS2=n61;yO1
op
iE
yO1
nJ2
lS3
SynthIfStep2
nZ2
x21
yP1
xS3+2);yP1
2
y7
l92
xS2=n61;yO1
cJump
iE
yO1
nJ2
lS3
SynthIfStep3
nZ2
x21
ByteCode.back()|=nJ2;yP1
xS3-1);yP1
2
y7
l92
xT1
xU+1)lL2
0;a<xS2;++a){if(ByteCode[a]==cJump&&ByteCode[a+1]==(nJ2|(xS2-1))){ByteCode[a+xS3-1);ByteCode[a+2
y7
l92
t93
ByteCode[a]){case
cAbsIf:case
cIf:case
cJump:case
cPopNMov:a+=2;lC
cFCall:case
cPCall:case
cFetch:a+=1;break;yT3
y53}
}
protected:void
xT1
size_t
value){xU=value;if(xU>lW3{StackMax=xU;cN
xK3
lW3;}
}
protected:std
xT3<unsigned>ByteCode;std
xT3
xE
Immed;std
xT3<std::pair<bool,FPoptimizer_CodeTree::lT1> >cN;size_t
xU;size_t
StackMax;private:void
incStackPtr(){if(xU+2>lW3
cN
xK3
StackMax=xU+2);}
tM1
bool
IsIntType,bool
IsComplexType>eI2
cL2{}
;eS3
void
AddOperation
nY2
eT3
unsigned
eat_count,unsigned
produce_count=1){EatNParams(eat_count);lW1
opcode);ProducedNParams(produce_count
lS3
lW1
unsigned
eT3
cL2<false,false>yR
false,true>yR
true,false>yR
true,true>);inline
void
lW1
unsigned
opcode){lW1
eT3
cL2<bool(nA
IsIntType
xE::c53),bool(nA
IsComplexType
xE::c53)>());}
}
yF1
eI2
SequenceOpCode
yF1
eI2
tP1{static
cK
AddSequence;static
cK
MulSequence;}
xM3
x01
long
count,cK&eT
xJ1;}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;lQ3
xS1{tI1
eI2
SequenceOpCode{Value_t
basevalue
tJ3
op_flip
tJ3
op_normal,op_normal_flip
tJ3
op_inverse,op_inverse_flip;}
yF1
cK
tP1
xE::AddSequence={xL1,cNeg
x7
cAdd,cSub,cRSub}
yF1
cK
tP1
xE::MulSequence={eF2
1),cInv,cMul,cMul,cDiv,cRDiv}
;
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode ByteCode
#define mImmed Immed
nK2
false,false>){xU1
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nK2
true,false>){xU1
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nK2
false,true>){xU1
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nK2
true,true>){xU1
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
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
lQ3
xS1;
#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
lQ3
xS1{
#ifndef FP_GENERATING_POWI_TABLE
extern
const
unsigned
char
powi_table[POWI_TABLE_SIZE];const
#endif
unsigned
char
powi_table[POWI_TABLE_SIZE]={0,1,1,1,2,1,2,1,xU3
4,1,2,xV3
2,1,xU3
8,eE3
xZ3
15,1,16,1,2,1,4,1,2,xV3
2,1,4,eE3
1,16,1,25,xZ3
27,5,8,3,2,1,30,1,31,3,32,1,2,1,xU3
8,1,2,xZ3
39,1,16,137,2,1,4,eE3
xV3
45,135,4,31,2,5,32,1,2,131,50,1,51,1,8,3,2,1,54,1,55,3,16,1,57,133,4,137,2,135,60,1,61,3,62,133,63,1,iC1
131,iC1
139,lA2
e5
30,1,130,137,2,31,lA2
e5
e5
130,eE3
1,e5
e5
2,1,130,133,iC1
61,130,133,62,139,130,137,e5
lA2
e5
e5
iC1
131,e5
e5
130,131,2,133,lA2
130,141,e5
130,eE3
1,e5
5,135,e5
lA2
e5
lA2
130,133,130,141,130,131,e5
e5
2,131}
;}
static
x43
c4=256;
#define FPO(x)
lQ3{class
PowiCache{private:int
iM[c4];int
iD1[c4];eS3
PowiCache():iM(),iD1(){iM[1]=1;}
bool
Plan_Add(cC1,int
count){cO1>=c4
cG
iD1[tL2+=count
n21
iM[tL2!=0
lT3
lX3
cC1){cO1<c4)iM[tL2=1
lT3
Start(size_t
value1_pos){for(int
n=2;n<c4;++n)iM[n]=-1;Remember(1,value1_pos);DumpContents();}
int
Find(cC1)const{cO1<c4){if(iM[tL2>=0){FPO(iJ3(iP3,"* I found %ld from cache (%u,%d)\n",value,(unsigned)cache[value],iL3 value]))n21
iM[tL2;t42-1
lT3
Remember(cC1,size_t
iY3){cO1>=c4)return;FPO(iJ3(iP3,"* Remembering that %ld can be found at %u (%d uses remain)\n",value,(unsigned)iY3,iL3 value]));iM[tL2=(int)iY3
lT3
DumpContents()const{FPO(for(int a=1;a<POWI_CACHE_SIZE;++a)if(cache[a]>=0||iL3 a]>0){iJ3(iP3,"== cache: sp=%d, val=%d, needs=%d\n",cache[a],a,iL3 a]);})}
int
UseGetNeeded(cC1){cO1>=0&&value<c4)return--iD1[tL2
n21
0;}
}
yF1
size_t
y8
long
count
iJ1
cK&eT
xJ1
xM3
cD1
size_t
apos,long
aval,size_t
bpos,long
bval
iJ1
unsigned
cumulation_opcode,unsigned
cimulation_opcode_flip
xJ1;void
lD1
cC1
iJ1
int
need_count,int
lE1=0){cO1<1)return;
#ifdef FP_GENERATING_POWI_TABLE
if(lE1>32)throw
false;
#endif
if(iM.Plan_Add(value,need_count
e42;long
y03
1;cO1<POWI_TABLE_SIZE){y03
powi_table[tL2
lM3&128){half&=127
lM3&64)y03-i02
FPO(iJ3(iP3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,value/half));lD1
half
xW3
iM.lX3
half)n21;}
iF1
half&64){y03-i02}
}
else
cO1&1)y03
value&((1<<POWI_WINDOW_SIZE)-1);else
y03
value/2;long
cO=value-half
lM3>cO||half<0)std::swap(half,cO);FPO(iJ3(iP3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,otherhalf))lM3==cO){lD1
half,iM,2,lE1+1);}
else{lD1
half
xW3
lD1
cO>0?cO:-cO
xW3}
iM.lX3
value);}
tI1
size_t
y8
cC1
iJ1
cK&eT
xJ1{int
y13=iM.Find(value
yG3
y13>=0)yQ
y13;}
long
y03
1;cO1<POWI_TABLE_SIZE){y03
powi_table[tL2
lM3&128){half&=127
lM3&64)y03-i02
FPO(iJ3(iP3,"* I want %ld, my plan is %ld * %ld\n",value,half,value/half));size_t
xT2=y8
half
cP1
if(iM
lB2
half)>0||xT2!=cQ1){iE1
xT2)x52
half,cQ1);}
x01
value/half,eT,tN2
size_t
iY3=cQ1
x52
value,iY3);iM.DumpContents()n21
iY3;}
iF1
half&64){y03-i02}
}
else
cO1&1)y03
value&((1<<POWI_WINDOW_SIZE)-1);else
y03
value/2;long
cO=value-half
lM3>cO||half<0)std::swap(half,cO);FPO(iJ3(iP3,"* I want %ld, my plan is %ld + %ld\n",value,half,value-half))lM3==cO)iA3
xT2=y8
half
cP1
cD1
xT2,half,xT2,half,iM,eT.op_normal,eT.op_normal_flip,tN2}
else{long
part1=half;long
part2=cO>0?cO:-cO;size_t
part1_pos=y8
part1
cP1
size_t
part2_pos=y8
part2
cP1
FPO(iJ3(iP3,"Subdivide(%ld: %ld, %ld)\n",value,half,otherhalf));cD1
part1_pos,part1,part2_pos,part2,iM,cO>0?eT.op_normal:eT.op_inverse,cO>0?eT.op_normal_flip:eT.op_inverse_flip,tN2}
size_t
iY3=cQ1
x52
value,iY3);iM.DumpContents()n21
iY3;}
eR1
cD1
size_t
apos,long
aval,size_t
bpos,long
bval
iJ1
unsigned
cumulation_opcode,unsigned
cumulation_opcode_flip
xJ1{int
a_needed=iM
lB2
aval);int
y23=iM
lB2
bval);bool
lC2=false;
#define DUP_BOTH() do{if(apos<bpos)iA3 tmp=apos;apos=bpos;bpos=tmp;lC2=!lC2;}FPO(iJ3(iP3,"-> " iV3 iV3"op\n",(unsigned)apos,(unsigned)bpos));iE1 apos);iE1 apos==bpos?cQ1:bpos);}while(0)
#define DUP_ONE(p) do{FPO(iJ3(iP3,"-> " iV3"op\n",(unsigned)p));iE1 p);}while(0)
if(a_needed>0){if(y23>0){nL2}
e13
bpos!=cQ1)nL2
else{lD2
lC2=!lC2;}
}
}
iF1
y23>0){if(apos!=cQ1)nL2
else
DUP_ONE(bpos);}
e13
apos==bpos&&apos==cQ1)lD2
tM2
cQ1&&bpos==y33
2){FPO(iJ3(iP3,"-> op\n"));lC2=!lC2;}
tM2
y33
2&&bpos==cQ1)FPO(iJ3(iP3,"-> op\n"));tM2
cQ1)DUP_ONE(bpos);iF1
bpos==cQ1){lD2
lC2=!lC2;}
else
nL2}
y43
lC2?cumulation_opcode_flip:cumulation_opcode,2);}
eR1
cR1
long
count,cK&eT
xJ1{while
eF3<256){int
y03
xS1::powi_table[count]lM3&128){half&=127;cR1
half,eT,tN2
count/=half;}
else
y53
if
eF3==1)return;if(!eF3&1)){y43
cSqr,1);cR1
count/2,eT,tN2}
else{iE1
cQ1);cR1
count-1,eT,tN2
y43
cMul,2);}
}
}
lQ3
xS1{eR1
x01
long
count,cK&eT
xJ1{if
eF3==0)t91
eT.basevalue);else{bool
tO2=false;if
eF3<0){tO2=true;count=-count;}
if(false)cR1
count,eT,tN2
iF1
count>1){PowiCache
iM;lD1
count,iM,1);size_t
xV1
tY3
GetStackTop();iM.Start(cQ1);FPO(iJ3(iP3,"Calculating result for %ld...\n",count));size_t
xU2=y8
count
cP1
size_t
n_excess=y33
xV1;if(n_excess>0||xU2!=xV1-1){synth.DoPopNMov(xV1-1,xU2);}
}
if(tO2)y43
eT.op_flip,1);}
}
}
#endif
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH
tE{lQ3
lY3{l01
eI2
Comp{}
;tM1>eI2
Comp<nA
cLess>lZ3<n9
cLessOrEq>lZ3<=n9
cGreater>lZ3>n9
cGreaterOrEq>lZ3>=n9
cEqual>lZ3==n9
i11>lZ3!=b;}
}
;}
tI1
eI2
cS1{Value_t
val;bool
known;cS1():val(),known(false){}
cS1(tQ1):val(v),known(true){x12
void
set(tQ1){known=true;val=v
lT3
set(eF2
eS1
Value_t),cS1
y9
l44
void
set(eF2
eS1
const
Value_t&),cS1
y9
l44
l01
void
set_if
tZ3
v
tO1
eS1
Value_t),cS1
y9&&lY3::Comp<Compare>()(val,v)l44
l01
void
set_if(tQ1
tO1
eS1
const
Value_t&),cS1
y9&&lY3::Comp<Compare>()(val,v)l44}
yF1
eI2
range{cS1
xE
min,max;range():min(),max(){}
range
tZ3
mi,Value_t
ma):min(mi),max(ma){}
range(bool,Value_t
ma):min(),max(ma){}
range
tZ3
mi,bool):min(mi),max(){}
void
set_abs();void
set_neg();}
yF1
bool
IsLogicalTrueValue(const
range
xE&p,bool
abs)yF1
bool
IsLogicalFalseValue(const
range
xE&p,bool
abs);}
#endif
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH
tE{enum
TriTruthValue{tA2,IsNever,Unknown}
yF1
range
xE
CalculateResultBoundaries(const
eU
yF1
bool
IsLogicalValue(const
eU
yF1
TriTruthValue
cF3
const
eU
yF1
xW1
GetEvennessInfo(const
eU{if(!tree
y01)return
Unknown;yZ1=l14;if(nA
isEvenInteger(value
e52
nA
isOddInteger(value
e62
tI1
xW1
GetPositivityInfo(const
eU{range
xE
p=CalculateResultBoundaries(tree
yG3
p
i4&&p
xK1>=eF2
e52
p.nT3
c71
e62
tI1
xW1
yB1
cI2&tree,bool
abs){range
xE
p=CalculateResultBoundaries(tree
yG3
IsLogicalTrueValue(p,abs
e52
IsLogicalFalseValue(p,abs
e62}
#endif
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH
tE{eR1
ConstantFolding(eU;}
#endif
lQ3{using
lQ3
FUNCTIONPARSERTYPES;using
tE;eI2
ComparisonSetBase{enum{eU3=0x1,Eq_Mask=0x2,Le_Mask=0x3,eV3=0x4,eW3=0x5,Ge_Mask=0x6}
;static
int
Swap_Mask(int
m)yQ(m&Eq_Mask)|((m&eU3)?eV3:0)|((m&eV3)?eU3:0);}
enum
cE1{Ok,BecomeZero,BecomeOne,xW2
x73
x62{cond_or,l03,l13,l23}
;}
yF1
eI2
ComparisonSet:public
ComparisonSetBase{eI2
tP2{lT1
a
i2
b;int
relationship;tP2():a(),b(),relationship(){}
}
;std
xT3<tP2>eW;eI2
Item{lT1
value;bool
cM2;Item():value(),cM2(false){}
}
;std
xT3<Item>cT1;int
xX1;ComparisonSet():eW(),cT1(),xX1(0){}
cE1
AddItem(cI2&a,bool
cM2,x62
type){for(size_t
c=0;c<cT1
tZ2
c)if(cT1[c].value
xK
a)){if(cM2!=cT1[c].cM2){l11
e51
case
l23:cT1.erase(cT1.begin()+c);xX1+=1
n21
xW2;case
l03:case
l13:e61
t42
xW2;}
Item
pole;pole.value=a;pole.cM2=cM2;cT1.push_back(pole)n21
Ok;}
cE1
AddRelationship(lT1
a,lT1
b,int
tT1,x62
type){l11
if(tT1==7)e51
lC
l23:if(tT1==7){xX1+=1
n11
lC
l03:case
l13:if(tT1==0)e61
y53
if(!(a.GetHash()<b.GetHash())){a.swap(b);tT1=Swap_Mask(tT1);}
for(size_t
c=0;c<eW
tZ2
c){if(eW[c].a
xK
a)&&eW[c].b
xK
b)){l11{int
y83=xV2|tT1;if(y83==7)e51
xV2=y83;y53
case
l03:case
l13:{int
y83=xV2&tT1;if(y83==0)e61
xV2=y83;y53
case
l23:{int
newrel_or=xV2|tT1;int
xY2=xV2&tT1;lF2
5&&xY2==0){xV2=eW3
n11
lF2
7&&xY2==0){xX1+=1;eW.erase(eW.begin()+c)n11
lF2
7&&xY2==Eq_Mask){xV2=Eq_Mask;xX1+=1
n11
continue;t42
xW2;}
}
tP2
comp;comp.a=a;comp.b=b;comp.relationship=tT1;eW.push_back(comp)n21
Ok;}
}
;nS1
Value_t,xD3
CondType>bool
ConstantFolding_LogicCommon(lT1&tree,CondType
xE1,bool
xZ2){bool
should_regenerate=false;ComparisonSet
xE
comp;for
xH{xD3
yX
cE1
eG3=yX
Ok;cI2&atree=t01;switch(atree
nC){case
cEqual
lG
Eq_Mask
t62
i11
lG
eW3
t62
cLess
lG
eU3
t62
cLessOrEq
lG
Le_Mask
t62
cGreater
lG
eV3
t62
cGreaterOrEq
lG
Ge_Mask
t62
cNot:eG3
yQ1
l8
0),true
t62
cNotNot:eG3
yQ1
l8
0),false,xE1);break;yT3
if(xZ2||IsLogicalValue(atree))eG3
yQ1,false,xE1);t93
eG3){ReplaceTreeWithZero:xI
0)n21
true;ReplaceTreeWithOne:xI
1);x0
yX
Ok:lC
yX
BecomeZero
e7
yX
BecomeOne:iB
yX
xW2:cH1
y53}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_LogicCommon: "
yJ
#endif
if(xZ2){tree.DelParams();}
else{for
yW{cI2&atree=t01;if(IsLogicalValue(atree))x81);}
}
for
lM2
0;a<comp.cT1
tZ2
a){if(comp.cT1[a].cM2){tQ2
cNot);r.cA
r.iN
iF1!xZ2){tQ2
cNotNot);r.cA
r.iN
else
tree.cA}
for
lM2
0;a<comp.eW
tZ2
a){tQ2
cNop);switch(comp.eW[a
tY){case
yX
eU3:r
iD
cLess);lC
yX
Eq_Mask:r
iD
cEqual);lC
yX
eV3:r
iD
cGreater);lC
yX
Le_Mask:r
iD
cLessOrEq);lC
yX
eW3:r
iD
i11);lC
yX
Ge_Mask:r
iD
cGreaterOrEq
yC2
r
yD1
comp.eW[a].a);r
yD1
comp.eW[a].b);r.iN
if(comp.xX1!=0)tree.yE
eF2
comp.xX1)));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_LogicCommon: "
yJ
#endif
return
true;eX3
tI3
yG1
ConstantFolding_AndLogic(iU3(tree.GetOpcode()==cAnd||tree.GetOpcode()==cAbsAnd)n21
nH
l03,true);}
yG1
ConstantFolding_OrLogic(iU3(tree.GetOpcode()==cOr||tree.GetOpcode()==cAbsOr)n21
nH
cond_or,true);}
yG1
ConstantFolding_AddLogicItems(iU3(tree.GetOpcode()==cAdd)n21
nH
l23,false);}
yG1
ConstantFolding_MulLogicItems(iU3(tree.GetOpcode()==cMul)n21
nH
l13,false);}
}
#include <vector>
#include <map>
#include <algorithm>
lQ3{using
lQ3
FUNCTIONPARSERTYPES;using
tE;eI2
CollectionSetBase{enum
xY1{Ok,xW2}
;}
yF1
eI2
CollectionSet:public
CollectionSetBase{eI2
cF1{lT1
value
i2
y02;bool
e0;cF1():value(),y02(),e0(false){}
cF1(cI2&v,cI2&f):value(v),y02(f),e0(false){}
}
;std::multimap<fphash_t,cF1>iP
i03
xD3
std::multimap<fphash_t,cF1>::yE3
xZ1;CollectionSet():iP(){}
xZ1
FindIdenticalValueTo(cI2&value){fphash_t
hash=value.GetHash();for(xZ1
i=iP.y12
hash);i!=iP.e71
hash;++i){cO1
xK
i
eN2.value
e42
i;eX3
iP
iW3;}
bool
Found(const
xZ1&b)yQ
b!=iP
iW3;}
xY1
AddCollectionTo(cI2&y02,const
xZ1&into_which){cF1&c=into_which
eN2;if(c.e0)c.y02
e9
y02);else{lT1
add;add
iD
cAdd);add
yD1
c.y02);add
e9
y02);c.y02.swap(add);c.e0=true;eX3
xW2;}
xY1
lG2
cI2&value,cI2&y02){const
fphash_t
hash=value.GetHash();xZ1
i=iP.y12
hash);for(;i!=iP.e71
hash;++i){if(i
eN2.value
xK
value
e42
AddCollectionTo(y02,i);}
iP.y93,std::make_pair(hash,cF1(value,y02)))n21
Ok;}
xY1
lG2
cI2&a)yQ
lG2
a,nD1
1)));}
}
yF1
eI2
ConstantExponentCollection{typedef
eO
yB3
i03
std
nF1
y22;std
xT3<y22>data;ConstantExponentCollection():data(){}
void
MoveToSet_Unique(const
Value_t&eT1&eU1){data.push_back(std
nF1(eT1()));data.back().second.swap(eU1
lS3
MoveToSet_NonUnique(const
Value_t&eT1&eU1){xD3
std
xT3<y22>::yE3
i=std::y12
data.l33
data
iW3,eE2,Compare1st()yG3
i!=data.e71
eE2){i
eN2.y93
eN2
iW3,eU1.l33
eU1
iW3);}
else{data.y93,std
nF1(eE2,eU1));}
}
bool
iX2{bool
changed=false;std::sort(data.l33
data
iW3,Compare1st());redo:for
lM2
0;a<data
tZ2
a
yF2
exp_a=data[a
eZ3;if(fp_equal(exp_a
tO1
1))cF2;for(n03
a+1;b<data
tZ2
b
yF2
exp_b=data[b
eZ3
yX3
y32=exp_b-exp_a;if(y32>=fp_abs(exp_a
l43
exp_diff_still_probable_integer=y32*eF2
16
yG3
tR2
exp_diff_still_probable_integer)&&!(tR2
exp_b)&&!tR2
y32))){yB3&a_set=lH2;yB3&b_set=data[b
t03;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantExponentCollection iteration:\n"
;tS2
cout);
#endif
if(isEvenInteger(exp_b)&&!isEvenInteger(y32+exp_a
nM
tmp2;tmp2
t83
tmp2
e83
b_set);tmp2
lP2
i2
tmp;tmp
iD
cAbs);tmp
yD1
tmp2);tmp
lP2;b_set
xK3
1);b_set[0
i23
tmp);}
a_set.insert(a_set
iW3,b_set.l33
b_set
iW3);yB3
b_copy=b_set;data.erase(data.begin()+b);MoveToSet_NonUnique(y32,b_copy);c31
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantExponentCollection iteration:\n"
;tS2
cout);
#endif
eK1}
t42
changed;}
#ifdef DEBUG_SUBSTITUTIONS
void
tS2
ostream&out){for
lM2
0;a<data
tZ2
a){out.precision(12);out<<data[a
eZ3<<": "
;for(yL3
lH2
tZ2
b){if(b>0)out<<'*'
iQ3
lH2[b],out);}
out<<std::endl;}
}
#endif
}
yF1
static
lT1
x31
lT1&value,bool&xP){switch(value
nC){case
cPow:{lT1
eE2=value
l8
1);value.yR1
n21
eE2
eL3
cRSqrt:value.yR1;xP=true
n21
nD1-0.5));case
cInv:value.yR1;xP=true
n21
nD1-1));yT3
y53
return
nD1
1));}
cK1
void
eV1
eW1&mul,cI2&tree,cI2&y02,bool&cG1
bool&xP){for
xH{lT1
value(t01)i2
eE2(x31
value,xP)yG3!y02
y01||y02
y11!=eF2
1.0
nM
e81;e81
t83
e81
e9
eE2);e81
e9
y02);e81
lP2;eE2.swap(e81);}
#if 0 /* FIXME: This does not work */
cO1
nC==cMul){if(1){bool
exponent_is_even=eE2
y01&&isEvenInteger(eE2
y11);for(yL3
value
yM3{bool
tmp=false
i2
val(value
l8
b))i2
exp(x31
val,tmp)yG3
exponent_is_even||(exp
y01&&isEvenInteger(exp
y11)nM
e81;e81
t83
e81
e9
eE2);e81
yD1
exp);e81.ConstantFolding(yG3!e81
y01||!isEvenInteger(e81
y11)){goto
cannot_adopt_mul;}
}
}
}
eV1
mul,value,eE2,cG1
xP);}
else
cannot_adopt_mul:
#endif
{if(mul.lG2
value,eE2)==CollectionSetBase::xW2)cH1}
}
}
yG1
ConstantFolding_MulGrouping(eU{bool
xP=false;bool
should_regenerate=false;eW1
mul;eV1
mul,tree,nD1
1)),cG1
xP)i03
std::pair<lT1,eO>eX1
i03
std::multimap<fphash_t,eX1>cU1;cU1
iJ;y42
eW1::xZ1
j=mul.iP.yC3
j!=mul.iP
iW3;++j){lT1&value=j
eN2.value
i2&eE2=j
eN2.y02;if(j
eN2.e0)eE2
lP2;const
fphash_t
eY1=eE2.GetHash();xD3
cU1::yE3
i=iJ.y12
eY1);for(;i!=iJ.e71
eY1;++i)if(i
eN2.first
xK
eE2)){if(!eE2
y01||!e91
y11
yF3
cH1
i
eN2.second.push_back(value);goto
skip_b;}
iJ.y93,std::make_pair(eY1,std::make_pair(eE2,eO(size_t(1),value))));skip_b:;}
#ifdef FP_MUL_COMBINE_EXPONENTS
ConstantExponentCollection
xE
eA1;y42
cU1::yE3
j,i=iJ.yC3
i!=iJ
iW3;i=j){j=i;++j;eX1&list=i
eN2;if
e72
y01)cL1
list.first
y11;if(!(eE2==xG1)eA1.MoveToSet_Unique(eE2,list
tW2
iJ.erase(i);}
}
if(eA1.iX2)cH1
#endif
if(should_regenerate){lT1
before=tree;before.lF1
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_MulGrouping: "
iQ3
before)xI1"\n"
;
#endif
tree.DelParams();y42
cU1::yE3
i=iJ.yC3
i!=iJ
iW3;++i){eX1&list=i
eN2;
#ifndef FP_MUL_COMBINE_EXPONENTS
if
e72
y01)cL1
list.first
y11;if(eE2==xG1
continue;if(e91
yF3{tree.AddParamsMove(list
tW2
tK3}
#endif
lT1
mul;mul
t83
mul
e83
list
tW2
mul
lP2;if(xP&&list.first
y01){if
e72
y11==eF2
1)/eF2
3
nM
cbrt;cbrt
iD
cCbrt);cbrt.e8
cbrt.lS2
cbrt
xX
0.5
nM
sqrt;sqrt
iD
cSqrt);sqrt.e8
sqrt.lS2
sqrt
xX-0.5
nM
rsqrt;rsqrt
iD
cRSqrt);rsqrt.e8
rsqrt.lS2
rsqrt
xX-1
nM
inv;inv
iD
cInv);inv.e8
inv.lS2
inv);tK3}
lT1
pow;pow
iD
cPow);pow.e8
pow
yD1
list.first);pow.lS2
pow);}
#ifdef FP_MUL_COMBINE_EXPONENTS
iJ.clear()lL2
0;a<i9
tZ2
a)cL1
i9[a
eZ3;if(e91
yF3{tree.AddParamsMove(i9[a]tW2
tK3
lT1
mul;mul
t83
mul
e83
i9[a]tW2
mul
lP2
i2
pow;pow
iD
cPow);pow.e8
pow.yE
eE2));pow.lS2
pow);}
#endif
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_MulGrouping: "
yJ
#endif
return!tree
xK
before);eX3
tI3
yG1
ConstantFolding_AddGrouping(eU{bool
should_regenerate=false;eW1
add;for
xH{if(t01
nC==cMul
cF2;if(add.lG2
t01)==CollectionSetBase::xW2)cH1}
xH2
remaining(eR);size_t
tM=0;for
xH{cI2&mulgroup=t01;if
eM3
nC==cMul){for(yL3
mulgroup
yM3{if
eM3
l8
b)y01
cF2;xD3
eW1::xZ1
c=add.FindIdenticalValueTo
eM3
l8
b)yG3
add.Found(c
nM
tmp
eM3
yA
CloneTag(tE3.DelParam(b);tmp
lP2;add.AddCollectionTo(tmp,c);cH1
goto
done_a;}
}
remaining[a]=true;tM+=1;done_a:;}
}
if(tM>0){if(tM>1){std
xT3<std::pair<lT1,size_t> >nZ;std::multimap<fphash_t,size_t>eZ1;bool
n13=false;for
xH
t21{for(yL3
t01
yM3{cI2&p=t01
l8
b);const
fphash_t
p_hash=p.GetHash();for(std::multimap<fphash_t,size_t>::const_iterator
i=eZ1.y12
p_hash);i!=eZ1.e71
p_hash;++i){if(nZ[i
eN2
eZ3
xK
p)){nZ[i
eN2
t03+=1;n13=true;goto
found_mulgroup_item_dup;}
}
nZ.push_back(std::make_pair(p,size_t(1)));eZ1.insert(std::make_pair(p_hash,nZ
t13-1));found_mulgroup_item_dup:;}
}
if(n13){lT1
eQ2;iA3
max=0;for(size_t
p=0;p<nZ
tZ2
p)if(nZ[p
t03<=1)nZ[p
t03=0;else{nZ[p
t03*=nZ[p
eZ3
nM2;if(nZ[p
t03>max){eQ2=nZ[p
eZ3;max=nZ[p
t03;}
}
}
lT1
group_add;group_add
iD
cAdd);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Duplicate across some trees: "
iQ3
eQ2)xI1" in "
yJ
#endif
for
xH
t21
for(yL3
t01
yM3
if(eQ2
xK
t01
l8
b)nM
tmp(t01
yA
CloneTag(tE3.DelParam(b);tmp
lP2;group_add
yD1
tmp);remaining[a]=false;y53
group_add
lP2
i2
group;group
t83
group
yD1
eQ2);group
yD1
group_add);group
lP2;add.lG2
group);cH1}
}
for
xH
t21{if(add.lG2
t01)==CollectionSetBase::xW2)cH1}
}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_AddGrouping: "
yJ
#endif
tree.DelParams();y42
eW1::xZ1
j=add.iP.yC3
j!=add.iP
iW3;++j){lT1&value=j
eN2.value
i2&coeff=j
eN2.y02;if(j
eN2.e0)coeff
lP2;if(coeff
y01){if(fp_equal(coeff
y11,xG1
cF2;if(fp_equal(coeff
y11
yF3{tree
yD1
value);tK3}
lT1
mul;mul
t83
mul
yD1
value);mul
yD1
coeff);mul.Rehash(yV
e8}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_AddGrouping: "
yJ
#endif
return
true;eX3
tI3}
lQ3{using
lQ3
FUNCTIONPARSERTYPES;using
tE
yF1
bool
ConstantFolding_IfOperations(iU3(tree.GetOpcode()==cIf||tree.GetOpcode()==cAbsIf);for(;;){if(tF3
cNot){i52
cIf);eV
yO3
eV
cV3
tV2.swap
xX2}
iF1
eV
e11{i52
tG3;eV
yO3
eV
cV3
tV2.swap
xX2}
else
nC3
yB1
eV,i42==tG3)i61
tree
yO3
tV2);x0
lI3
tree.Become
xX2
x0
n01
if(tF3
cIf||tF3
tG3{lT1
cond=eV
i2
n23;n23
i62==cIf?cNotNot:cAbsNotNot);n23
y52
1));ConstantFolding(n23)i2
n33;n33
i62==cIf?cNotNot:cAbsNotNot);n33
y52
2));ConstantFolding(n33
yG3
n23
y01||n33.IsImmed(nM
t0;t0
i62);t0
y52
1));t0.nI
1));t0.nI
2));t0
lP2
i2
t1;t1
i62);t1
y52
2));t1.nI
1));t1.nI
2));t1
i3
cond
nC
y31
0,cond
l8
0)yV
nG1
1,t0
yV
nG1
2,t1
lQ2}
if(lI2
tU2
nC&&(lI2
cIf||lI2
cAbsIf
nM&xB3=tV2
i2&leaf2=tU2;if(xB3
l8
0)xF1
0))&&t52
1))||xB3
l8
2)xF1
2))nM
t0;t0
iG1
t0.nI
0));t0
y62
1));t0
y72
1));t0
lP2
i2
t1;t1
iG1
t1.nI
0));t1
y62
2));t1
y72
2));t1
i3
xB3
nC
y31
0
yK3
0)yV
nG1
1,t0
yV
nG1
2,t1
lQ2
if
t52
1))&&xB3
l8
2)xF1
2)nM
t2;t2
iG1
t2
yD1
eV);t2
y62
0));t2
y72
0));t2
i3
xB3
nC
yV
nG1
0,t2
y31
2
yK3
2)y31
1
yK3
1)lQ2
if
t52
2))&&xB3
l8
2)xF1
1)nM
eR2;eR2
iD
leaf2
nC==cIf?cNot:cY3);eR2
y72
0));eR2
lP2
i2
t2;t2
iG1
t2
yD1
eV);t2
y62
0));t2
yD1
eR2);t2
i3
xB3
nC
yV
nG1
0,t2
y31
2
yK3
2)y31
1
yK3
1)lQ2}
lT1&xY=tV2
i2&y5=tU2;if(xY
xK
y5)){tree
yO3
tV2
lQ2
const
OPCODE
op1=xY
nC;const
OPCODE
op2=y5
nC;if(op1==op2){if(xY
y41==1){lT1
lO
0));y82
0));y92)n4
if(xY
y41==2&&y5
y41==2){if(xY
l8
0)xK
y5
l8
0)nM
param0=xY
l8
0)i2
lO
1));y82
1));y92
iO
param0)n4
if(xY
l8
1)xK
y5
l8
1)nM
param1=xY
l8
1)i2
lO
0));y82
0));y92
iO
nC1
iO
param1
lQ2}
if(op1==yN3
cMul
lK2
cAnd
lK2
cOr
lK2
cAbsAnd
lK2
cAbsOr
lK2
cMin
lK2
cMax){eO
n43;c5{for(n03
y5
y41;b-->0;){x72
y5
l8
b))){if(n43.empty()){xY.tH3
lF1}
n43.push_back(xY
l8
a));y5.DelParam(b);xY.yA2
y53}
}
if(!n43.empty()){xY
lP2;y5
lP2
l7
op1
yV
SetParamsMove(n43)n4}
}
if(op1==yN3
cMul||(op1==cAnd
t43
y5))||(op1==cOr
t43
y5))){c5
x72
y5)){xY.lF1
xY.yA2
xY
lP2
i2
cV1=y5;y5=tN
op1==yN3
cOr)l63
op1
iO
cV1)n4}
if((op1==cAnd
lK2
cOr)&&op2==cNotNot){lT1&n53=y5
l8
0);c5
x72
n53)){xY.lF1
xY.yA2
xY
lP2
i2
cV1=n53;y5=tN
op1==cOr)l63
op1
iO
cV1)n4}
if(op2==cAdd||op2==cMul||(op2==cAnd
t43
xY))||(op2==cOr
t43
xY))){for
lM2
y5
l51
y5
l8
a)xK
xY)){y5.tH3
yA2
y5
lP2
i2
cW1=xY;xY=tN
op2==cAdd||op2==cOr)l63
op2
iO
cW1)n4}
if((op2==cAnd||op2==cOr)&&op1==cNotNot){lT1&n63=xY
l8
0)lL2
y5
l51
y5
l8
a)xK
n63)){y5.tH3
yA2
y5
lP2
i2
cW1=n63;xY=tN
op2==cOr)l63
op2
iO
cW1)n4
eX3
tI3}
#include <limits>
lQ3{using
lQ3
FUNCTIONPARSERTYPES;using
tE
yF1
int
maxFPExponent()yQ
std::numeric_limits
xE::max_exponent;}
yG1
x41
Value_t
base,Value_t
eE2){if(base<xG1
return
true;if(fp_equal(base,xG1||fp_equal(base
tO1
1))cG
return
eE2>=eF2
maxFPExponent
xE())/fp_log2(base);}
yG1
ConstantFolding_PowOperations(iU3(tree.GetOpcode()==cPow);nG&&l41
yF2
const_value
eP3
lR,yB2);xI
const_value)n21
tI3
if(l41&&fp_equal(yB2
yF3{tree
yO3
eV
lQ2
nG&&fp_equal(lR
yF3{xI
1)n21
tI3
nG&&lI2
cMul){bool
yD2=false
yX3
lT2=lR
i2
mulgroup=tV2
lL2
mulgroup
l51
mulgroup
l8
a)yE2
imm=mulgroup
l8
a)y11;{if(x41
lT2,imm
l43
lU2
eP3
lT2,imm
yG3
fp_equal(lU2,xG1)break;if(!yD2){yD2=true;yE1
lF1}
lT2=lU2;yE1
yA2
y53}
if(yD2){yE1
Rehash();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before pow-mul change: "
yJ
#endif
eV
yO3
eB1
lT2));tV2
yO3
mulgroup);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After pow-mul change: "
yJ
#endif
}
}
if(l41&&tF3
cMul
yF2
lV2=yB2
yX3
yG2=1.0;bool
yD2=false
i2&mulgroup=eV
lL2
mulgroup
l51
mulgroup
l8
a)yE2
imm=mulgroup
l8
a)y11;{if(x41
imm,lV2
l43
t31
eP3
imm,lV2
yG3
fp_equal(t31,xG1)break;if(!yD2){yD2=true;yE1
lF1}
yG2*=t31;yE1
yA2
y53}
if(yD2){yE1
Rehash()i2
eH3;eH3
iD
cPow);eH3
e83
tree.l52));eH3.yT2
yV
SetOpcode(cMul
iO
eH3
yV
AddParam(eB1
yG2)lQ2}
if(tF3
cPow&&l41&&eV
l8
1)yE2
a=eV
l8
iG
yX3
b=yB2
yX3
c=a*b;if(isEvenInteger(a)&&!isEvenInteger(c
nM
n73;n73
iD
cAbs);n73.nI
0)cV3
n73.Rehash(yV
nG1
0,n73);}
else
tree.SetParam(0,eV
l8
0)y31
1,eB1
c));eX3
tI3}
lQ3{using
lQ3
FUNCTIONPARSERTYPES;using
tE;eI2
l5{enum
eS2{MakeFalse=0,MakeTrue=1,i32=2,nA3=3,MakeNotNotP0=4,MakeNotNotP1=5,MakeNotP0=6,MakeNotP1=7,Unchanged=8
x73
lW2{Never=0,Eq0=1,Eq1=2,yQ3=3,yR3=4}
;eS2
if_identical;eS2
lX2
4];eI2{eS2
what:4;lW2
when:4;}
iL1,iM1,iN1,iO1
yF1
eS2
Analyze(cI2&a,cI2&b)const{if(a
xK
b
e42
if_identical;range
xE
tW
a
yJ3
p1=CalculateResultBoundaries(b
yG3
p0.max
lH3
p1
i4){if(p0
nP<p1
xK1&&lX2
0]cE
0];if(p0
nP<=p1
xK1&&lX2
1]cE
1];}
if(p0
i4&&p1.iI1{if(p0
xK1>p1
nP&&lX2
2]cE
2];if(p0
xK1>=p1
nP&&lX2
3]cE
3];}
if(IsLogicalValue(a)){if(iL1
yP3
iL1.when,p1
e42
iL1.what;if(iN1
yP3
iN1.when,p1
e42
iN1.what;}
if(IsLogicalValue(b)){if(iM1
yP3
iM1.when,p0
e42
iM1.what;if(iO1
yP3
iO1.when,p0
e42
iO1.what;eX3
Unchanged;}
cK1
bool
TestCase(lW2
when,const
range
xE&p){if(!p
i4||!p.nT3
cG
switch(when){case
Eq0
nU1==eF2
0.0)&&p
nP==p
xK1;case
Eq1
nU1==eF2
1.0)&&p
nP==p
nP;case
yQ3
nU1>xL1&&p
nP<=eF2
1);case
yR3
nU1>=xL1
c71
1);yT3;eX3
tI3}
;lQ3
RangeComparisonsData{static
const
l5
Data[6]={{l5
n83
tZ
Unchanged,l5::tZ
Unchanged
lG1
Eq1
nV1
Eq1
nW1
Eq0
nX1
Eq0}
}
,{l5::x82
n93
Unchanged,l5
n93
Unchanged
lG1
Eq0
nV1
Eq0
nW1
Eq1
nX1
Eq1}
}
,{l5::x82
n93
i32,l5::tZ
MakeFalse
nW1
yQ3
nV1
yR3
xN,{l5
n83
Unchanged,l5
n93
tZ
nA3
nW1
yR3
nV1
yQ3
xN,{l5::x82::tZ
tZ
MakeTrue,l5::i32
lG1
yR3
nX1
yQ3
xN,{l5
n83
tZ
nA3,l5::Unchanged,l5
nH1
lG1
yQ3
nX1
yR3
xN}
;}
yG1
ConstantFolding_Comparison(eU{using
lQ3
RangeComparisonsData;assert(tree.GetOpcode()>=cEqual&&tree.GetOpcode()<=cGreaterOrEq);switch(Data[i42-cEqual].Analyze(eV,tV2)){case
l5::MakeFalse:xI
0);x0
l5
nH1:xI
1)l53
nA3:i52
cEqual)l53
i32:i52
i11)l53
MakeNotNotP0:i52
cNotNot
iP1
1)l53
MakeNotNotP1:i52
cNotNot
iP1
0)l53
MakeNotP0:i52
cNot
iP1
1)l53
MakeNotP1:i52
cNot
iP1
0)l53
Unchanged:;}
if(l41)switch(eV
nC){case
cAsin:lM
fp_sin(eT2
cAcos:lM
fp_cos(yB2))yV
SetOpcode(i42==cLess?cGreater:i42==cLessOrEq?cGreaterOrEq:i42==cGreater?cLess:i42==cGreaterOrEq?cLessOrEq:i42);x0
cAtan:lM
fp_tan(eT2
cLog:lM
fp_exp(eT2
cSinh:lM
fp_asinh(eT2
cTanh:if(fp_less(fp_abs(yB2)yF3{lM
fp_atanh(yB2))lQ2
break;yT3
y53
return
tI3}
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;lQ3{
#ifdef DEBUG_SUBSTITUTIONS
yG
double
d){union{double
d;uint_least64_t
h;}
i72
d=d;lR1
h
nY1
#ifdef FP_SUPPORT_FLOAT_TYPE
yG
float
f){union{float
f;uint_least32_t
h;}
i72
f=f;lR1
h
nY1
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
yG
long
double
ld){union{long
double
ld;eI2{uint_least64_t
a
tJ3
short
b;}
s;}
i72
ld=ld;lR1
s.b<<data.s.a
nY1
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
yG
long
ld){o<<"("
<<std::hex<<ld
nY1
#endif
#endif
}
tE{lN
nE)){}
lN
const
Value_t&i
yA
xH3
nE
i)){data
xA
#ifdef FP_SUPPORT_CXX11_MOVE
lN
Value_t&&i
yA
xH3
nE
i82
i))){data
xA
#endif
lN
unsigned
v
yA
VarTag
nE
l93,v)){data
xA
lN
x22
o
yA
OpcodeTag
nE
o)){data
xA
lN
x22
o,unsigned
f
yA
FuncOpcodeTag
nE
o,f)){data
xA
lN
cI2&b
yA
CloneTag
nE*b.data)){}
tI1
lT1::~lU1(){}
lB
ReplaceWithImmed(tS1{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Replacing "
iQ3*this
yG3
IsImmed())OutFloatHex(std::cout,GetImmed())xI1" with const value "
<<i;OutFloatHex(std::cout,i)xI1"\n"
;
#endif
data=new
xQ2
xE(i);}
tI1
eI2
ParamComparer{iY2()(cI2&a,cI2&b)const{if(a
nM2!=b
nM2)return
a
nM2<b
nM2
n21
a.GetHash()<b.GetHash();}
}
xM3
xQ2
xE::Sort(){switch(Opcode){case
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
i11:std::sort(y63
l33
y63
end(),ParamComparer
xE());lC
cLess
lY
cGreater;}
lC
cLessOrEq
lY
cGreaterOrEq;}
lC
cGreater
lY
cLess;}
lC
cGreaterOrEq
lY
cLessOrEq;}
break;yT3
y53}
lB
AddParam(cI2&param){xZ.push_back(param);}
lB
nU2
lT1&param){xZ.push_back(lT1());xZ.back().swap(param);}
lB
SetParam(size_t
which,cI2&b)nI1
which
l73
xZ[which]=b;}
lB
nG1
size_t
which,lT1&b)nI1
which
l73
xZ[which
i23
b);}
lB
AddParams(const
nJ){xZ.insert(xZ
iW3,lY2.l33
lY2
iW3);}
lB
AddParamsMove(nJ)iA3
endpos=xZ
t13,added=lY2
t13;xZ
xK3
endpos+added,lT1());for(size_t
p=0;p<added;++p)xZ[endpos+p
i23
lY2[p]);}
lB
AddParamsMove(nJ,size_t
lZ2)nI1
lZ2
l73
DelParam(lZ2);AddParamsMove(tL1}
lB
SetParams(const
nJ){eO
tmp(tL1
xZ.swap(tmp);}
lB
SetParamsMove(nJ){xZ.swap(tL1
lY2.clear();}
#ifdef FP_SUPPORT_CXX11_MOVE
lB
SetParams(eO&&lY2){SetParamsMove(tL1}
#endif
lB
DelParam(size_t
index){eO&eD3=xZ;
#ifdef FP_SUPPORT_CXX11_MOVE
y63
erase(y63
begin()+index);
#else
eD3[index].data=0;for(size_t
p=index;p+1<eD3
tZ2
p)eD3[p].data.UnsafeSetP(&*eD3[p+1
l73
eD3[yV3-1].data.UnsafeSetP(0);eD3
xK3
yV3-1);
#endif
}
lB
DelParams(){xZ.clear();}
yG1
lT1::IsIdenticalTo(cI2&b)const{if(&*data==&*b.data)return
true
n21
data->IsIdenticalTo(*b.data);}
yG1
xQ2
xE::IsIdenticalTo(const
xQ2
xE&b)const{if(Hash!=b.Hash
cG
if(Opcode!=l54
cG
switch(Opcode){case
cImmed:return
fp_equal(Value,l64;case
l93:return
l72==b.l62
case
cFCall:case
cPCall:if(l72!=b.l72
cG
break;yT3
y53
if(yV3!=b.yV3
cG
for
lM2
0;a<eD3
tZ2
a){if(!eD3[a]xK
b.eD3[a])cG
eX3
true;}
lB
Become(cI2&b){if(&b!=this&&&*data!=&*b.data){DataP
tmp=b.data;lF1
data.swap(tmp);}
}
lB
c02){if(GetRefCount()>1)data=new
xQ2
xE(*data);}
tI1
lT1
lT1::GetUniqueRef(){if(GetRefCount()>1)return
lT1(*this,CloneTag())n21*this;}
i0):yS
cNop),Value(),n8
i0
const
xQ2&b):yS
l54),Value(l64,l72(b.cZ1,eD3(b.eD3),Hash(b.Hash),Depth(b.Depth),tU1
b.l82){}
i0
tS1:yS
cImmed),Value(i),n8
#ifdef FP_SUPPORT_CXX11_MOVE
i0
xQ2
xE&&b):yS
l54),Value(i82
l64),l72(b.cZ1,eD3(i82
b.eD3)),Hash(b.Hash),Depth(b.Depth),tU1
b.l82){}
i0
Value_t&&i):yS
cImmed),Value(i82
i)),n8
#endif
i0
x22
o):yS
o),Value(),n8
i0
x22
o,unsigned
f):yS
o),Value(),l72(f),eD3(),Hash(),Depth(1),tU1
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
lQ3
FUNCTIONPARSERTYPES;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
lQ3{eR1
tV1
nT,std
cB&done,std::ostream&o){for
xH
tV1
t01,done,o);std::ostringstream
buf
iQ3
tree,buf);done[tree.GetHash()].insert(buf.str());}
}
#endif
tE{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
eR1
DumpHashes(cI){std
cB
done;tV1
tree,done,o);for(std
cB::const_iterator
i=done.yC3
i!=done
iW3;++i){const
std::set<std
cG3>&flist=i
eN2;if(flist
t13!=1)o<<"ERROR - HASH COLLISION?\n"
;for(std::set<std
cG3>::const_iterator
j=flist.yC3
j!=flist
iW3;++j){o<<'['<<std::hex<<i->first.hash1<<','<<i->first.hash2<<']'<<std::dec;o<<": "
<<*j<<"\n"
;}
}
}
eR1
DumpTree(cI){xN3
iN3;switch
tT2{case
cImmed:o<<l14
yS3
l93:o<<"Var"
<<(tree.GetVar()-l93)yS3
cAdd:iN3"+"
;lC
cMul:iN3"*"
;lC
cAnd:iN3"&"
;lC
cOr:iN3"|"
;lC
cPow:iN3"^"
;break;yT3
iN3;o<<t53
tT2;lE3
cFCall||i42==cPCall)o<<':'<<tree.GetFuncNo();}
o<<'(';if(eR<=1&&sep2[1])o<<(sep2+1)<<' ';for
xH{if(a>0)o<<' '
iQ3
t01,o
yG3
a+1<eR)o<<sep2;}
o<<')';}
eR1
DumpTreeWithIndent(cI,const
std
cG3&indent){o<<'['<<std::hex<<(void*)(&tree.l52))<<std::dec<<','<<tree.GetRefCount()<<']';o<<indent<<'_';switch
tT2{case
cImmed:o<<"cImmed "
<<l14;o<<'\n'
yS3
l93:o<<"VarBegin "
<<(tree.GetVar()-l93);o<<'\n'
n21;yT3
o<<t53
tT2;lE3
cFCall||i42==cPCall)o<<':'<<tree.GetFuncNo();o<<'\n';}
for
xH{std
cG3
ind=indent;for(size_t
p=0;p<ind
t13;p+=2)if(ind[p]=='\\')ind[p]=' ';ind+=(a+1<eR)?" |"
:" \\"
;DumpTreeWithIndent(t01,o,ind);}
o<<std::flush;}
#endif
}
#endif
using
lQ3
lC1;using
lQ3
FUNCTIONPARSERTYPES;
#include <cctype>
lQ3
lC1{unsigned
ParamSpec_GetDepCode(const
eL2&b){switch(b.first){case
ParamHolder:{cP*s=(cP*)b.second
n21
s->depcode
eL3
SubFunction:{cQ*s=(cQ*)b.second
n21
s->depcode;}
yT3
y53
return
0;}
eR1
DumpParam(const
eL2&yQ2
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
tJ3
yH2
0;t33
t63{const
ParamSpec_NumConstant
xE
nB3
const
ParamSpec_NumConstant
xE*t11;using
lQ3
FUNCTIONPARSERTYPES;o.precision(12);o<<param.constvalue;y53
case
ParamHolder:{cP
nB3
cP*t11;o<<ParamHolderNames[param.index];yH2
param.constraints;y53
case
SubFunction:{cQ
nB3
cQ*t11;yH2
param.constraints;yH
GroupFunction){if(param.lH1==cNeg){o<<"-"
;n2}
iF1
param.lH1==cInv){o<<"/"
;n2}
else{std
cG3
opcode=t53((x22)param.lH1).substr(1)lL2
0;a<opcode
tZ2
a)opcode[a]=(char)std::toupper(opcode[a]);o<<opcode<<"( "
;n2
o<<" )"
;}
}
else{o<<'('<<t53((x22)param.lH1)<<' ';yH
PositionalParams)o<<'[';yH
SelectedParams)o<<'{';n2
if
eP2.n1!=0)o<<" <"
<<param.data.n1<<'>';yH
PositionalParams)o<<"]"
;yH
SelectedParams)o<<"}"
;o<<')';}
y53
t93
ImmedConstraint_Value
iO3
ValueMask)){case
ValueMask:lC
Value_AnyNum:lC
x92:o<<"@E"
;lC
Value_OddInt:o<<"@O"
;lC
tX1:o<<"@I"
;lC
Value_NonInteger:o<<"@F"
;lC
t41:o<<"@L"
;nC3
ImmedConstraint_Sign
iO3
SignMask)){case
SignMask:lC
Sign_AnySign:lC
nJ1:o<<"@P"
;lC
t51:o<<"@N"
;nC3
ImmedConstraint_Oneness
iO3
OnenessMask)){case
OnenessMask:lC
Oneness_Any:lC
Oneness_One:o<<"@1"
;lC
Oneness_NotOne:o<<"@M"
;nC3
ImmedConstraint_Constness
iO3
ConstnessMask)){case
ConstnessMask:lC
tW1:if(nJ3.first==ParamHolder){cP
nB3
cP*t11;if(param.index<2)y53
o<<"@C"
;lC
Constness_NotConst:o<<"@V"
;lC
Oneness_Any:y53}
eR1
DumpParams
nY2
paramlist,unsigned
count,std::ostream&o){for(t71=0;a<count;++a){if(a>0)o<<' ';const
eL2&param=eN1
xE(paramlist,a);DumpParam
xE(param,o)tJ3
depcode=ParamSpec_GetDepCode(param
yG3
depcode!=0)o<<"@D"
<<depcode;}
}
}
#include <algorithm>
using
lQ3
lC1;using
lQ3
FUNCTIONPARSERTYPES;lQ3{cP
plist_p[37]={{2,0,0x0}
nU
0,0x4}
nU
nJ1,0x0}
nU
t51|Constness_NotConst,0x0}
nU
Sign_NoIdea,0x0}
nU
t41,0x0}
,{3,Sign_NoIdea,0x0}
,{3,0,0x0}
,{3,t41,0x0}
,{3,0,0x8}
,{3,Value_OddInt,0x0}
,{3,Value_NonInteger,0x0}
,{3,x92,0x0}
,{3,nJ1,0x0}
,{0,t51|lV{0,lV{0,nJ1|lV{0,x92|lV{0,tW1,0x1}
,{0,tX1|nJ1|lV{0,tY1
tW1,0x1}
,{0,tY1
lV{0,Oneness_One|lV{0,t41|lV{1,lV{1,x92|lV{1,tY1
lV{1,tX1|lV{1,nJ1|lV{1,t51|lV{6,0,0x0}
,{4,0,0x0}
,{4,tX1,0x0}
,{4,lV{4,0,0x16}
,{5,0,0x0}
,{5,lV}
yF1
eI2
plist_n_container{static
const
ParamSpec_NumConstant
xE
plist_n[20];}
yF1
const
ParamSpec_NumConstant
xE
plist_n_container
xE::plist_n[20]={{eF2-2
i1-1
i1-0.5
i1-0.25
i1
0
i92
fp_const_deg_to_rad
tN3
fp_const_einv
tN3
fp_const_log10inv
xE(i1
0.5
i92
fp_const_log2
xE(i1
1
i92
fp_const_log2inv
xE(i1
2
i92
fp_const_log10
tN3
fp_const_e
tN3
fp_const_rad_to_deg
tN3-fp_const_pihalf
xE(),y51{xL1,y51{fp_const_pihalf
xE(),y51{fp_const_pi
xE(),y51}
;cQ
plist_s[517]={{{1,15,iC2
398,iC2
477,iC2
15,cNeg,GroupFunction,0}
,tW1,0x1
iA2
15,yI2
24,yI2
465,yI2
466,yI2
498,cInv
iX3
327995
x7
nD3
48276
x7
l6
260151
x7
l6
470171
x7
l6
169126
x7
l6
48418
x7
l6
1328
x7
l6
283962
x7
l6
169275
x7
l6
39202
x7
l6
283964
x7
l6
283973
x7
l6
476619
x7
l6
296998
x7
l6
47
x7
SelectedParams,0}
,0,cM
161839
x7
l6
25036
x7
l6
35847
x7
l6
60440
x7
l6
30751
x7
l6
183474
x7
l6
259318
x7
l6
270599
x7
l6
60431
x7
l6
259119
x7
l6
332066
x7
l6
7168
x7
l6
197632
x7
l6
291840
x7
l6
283648
x7
l6
238866
x7
l6
239902
x7
l6
31751
x7
l6
244743
x7
l6
384022
x7
SelectedParams,0}
,0,cM
385262
x7
l6
386086
x7
l6
393254
x7
SelectedParams,0}
,0,0x5
nN
393254
x7
l6
386095
x7
l6
387312
x7
l6
18662
x7
l6
61670
x7
l6
387397
x7
l6
247855
x7
SelectedParams,0}
,0,0x1
nN
342063
x7
l6
297007
x7
l6
15820
x7
l6
393263
x7
l6
393263
x7
SelectedParams,0}
,0,0x5
nN
161847
x7
l6
258103
x7
l6
249073
x7
l6
249076
x7
iQ
0,0
x7
nK
0,0
tZ1
1,45
x7
nK
1,53
x7
nK
1,54
x7
nK
1,55
x7
nK
1,56
x7
nK
1,26
x7
nK
1,259
eU2
0x16
iA2
272
tZ1
1,323
eU2
0x16
iA2
0
x7
nK
1,21
x7
nK
1,447
eU2
0x4
iA2
449
eU2
0x4
iA2
0
eU2
0x4
iA2
0
eZ
2}
,0,0x4
iA2
15
x7
nK
1,24
eZ
2}
,0,0x0
nN
58392
tZ1
0,0
eZ
1}
,nJ1,0x0
nN
24591
iB2
33807
iB2
48143
iB2
285720
iB2
290840
iB2
305152,l9
312400,l9
39202,l9
122918,l9
421926,l9
429094,l9
443430,l9
317834,l9
329098,l9
7633,l9
7706,l9
7730,l9
38,l9
50587,l9
406528,l9
24583,l9
31751,l9
405511,l9
321551
xN1
327713,l9
322596,l9
90409,l9
335174,l9
327050,l9
493606,l9
496678,l9
503846,l9
516134,l9
7217,l9
333875,l9
336896,l9
524326,l9
509952,l9
286727,l9
89103,l9
92175,l9
296976,l1
cM
324623,l1
0x14
nN
332815,l1
0x10}
,{{3,7340056,l1
cM
289092,l9
93200
xN1
337935
xO1
7340060,l1
iD2
7340176,l9
338959
xO1
7340061
xN1
7206,l9
7168,l9
357414,l9
368678,l9
370745,l1
0x7}
,{{3,7340177,l9
39277,l1
cM
426398,l1
iD2
40272286
xN1
490910,l1
iD2
40336798
xN1
50600,l9
426462
xN1
490974
xN1
370726,l1
0x6
nN
371750,l1
0x6
nN
428070
xO1
40336862
xN1
38378,l9
50671
xO1
47662080,l9
477184,l9
568320,l9
371727,l1
0x7}
,{{3,15779306,l9
370703,l1
0x7
nN
39277,l9
39279,l1
0x4}
,{{3,15779238,l9
39338,l1
cM
436262,l9
508966,l9
39409,l1
cM
296998,l1
cM
35847,l9
15,l1
cM
377894,l9
386063,l1
0x1
nN
15,l9
7192,l9
123928,l9
122904,l9
30751,l9
57,l9
7456,l9
15674
xO1
67579935,l9
39237,l9
58768,l9
62924,l9
122880,l9
15760
xO1
64009216,l1
0x0}
,{{0,0
xD
0,0,l31
2,eC1
2,eD1
3,eC1
3,eD1
38
xD
1,38,l31
14
xD
1,57
xD
1,16,eV2
0x0
nN
471103,eV2
0x1
iA2
303
xD
1,323,yU3
0x0
nN
471363,eV2
0x16
iA2
293,eC1
294,eD1
295
xD
1,296,l31
400
xD
1,0
xD
1,460
xD
1,465
xD
1,16,eV2
0x1
iA2
57,yU3
0x1
iA2
0,l31
21
xD
1,15,eV2
0x0
nN
24591
xD
1,24,l31
517,yU3
0x0
nN
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
50176,lF,178176
y0
tR3
283648,lF,19456,lF,27648,lF,91136,lF,86016,lF,488448,lF,14342,lF,58375,lF,46147
y0
cM
46151,lF,284679,lF,7183,lF,46159
y0
cM
38993
y0
cM
50262,lF,50249,lF,283808,lF,284835,lF,24822,lF,10240,lF,11264,lF,7170,lF,7168,lF,17408,lF,164864,lF,237568,lF,242688
y0
0x14
nN
476160,lF,25607,lF,122895,lF,50252,lF,39374,lF,50183,lF,7192,lF,122911,lF,252979,lF,46155,lF,38919,lF,50268,lF,50269,lF,50253,lF,46191,lF,50296,lF,7563
y0
0x10
nN
416811,lF,416819,lF,40047,lF,46192
y0
cM
415795,lF,40048
y0
cM
415787,lF,39016
y0
0x5
nN
39326
y0
cM
39326,lF,39332
y0
0x5
nN
39333
y0
0x1
nN
50590
y0
cM
50590,lF,39338
y0
cM
39338,lF,39335
y0
0x5
nN
15786
y0
cM
146858,lF,39372,lF,39379,lF,39380,lF,39390
y0
cM
50654
y0
cM
50654,lF,24
y0
0x6
nN
62,lF,24,lF,62
y0
0x6
nN
43,lF,43
y0
cM
51,lF,51
y0
cM
50270,lF,50176
y0
cM
50271,lF,39159,lF,39183
y0
cM
7168
y0
cM
31744,lF,100352,lF,31746,lF,101400,lF,39409
y0
cM
39411
y0
cM
39411,lF,39420,lF,39420
y0
cM
15,lF,39026
y0
0x5
nN
39422,lF,16384,lF,62853,lF,15360,lF,15
y0
0x1
nN
16,lF,7183
y0
0x1
nN
7172,cPow,yN1,nJ1,0x0
nN
24591,cPow
iX3
50200,cPow
iX3
63521,cPow
iX3
62500,cPow
iX3
50453,cPow
iX3
62488,cPow
y73
tO3
7,tO3
194,tO3
0,cAcos
tP3
cAcosh
tP3
cAsin
tP3
cAsinh
nR
120,cAsinh
tP3
cAtan,nD3
306176
l24,nE3
cAtan2
tP3
cAtanh
nR
246,cCeil
tP3
cCeil,t4
1,0,cN2
0,cCos,t4
1,7,cN2
92,cN2
93,cN2
120,cN2
236,cN2
255,cN2
214
eX2
236
eX2
464
eX2
0,cCosh,t4
1,0
eX2
0,cExp
nR
7,cExp
nR
92,cExp
tP3
c63
7,c63
92,c63
246,cFloor
tP3
cFloor,lA
cM
309540,tQ3
nD3
316708,tQ3
nD3
316724,tQ3
l0
3,32513024,eY2
34627584
i5
31493120,eY2
89213952
i5
149042176
i5
246647808
i5
301234176
i5
494360576
i5
498558976
i5
62933520
i5
62933520,eY2
62933526
i5
62933526,eY2
24670208
i5
579378176
i5
573578240
i5
32513024
i5
566254592
i5
7900160
i5
588822528,cIf
nR
120,cInt
nR
246,iE2
0,iE2
7,iE2
31,iE2
194,iE2
363,iE2
15,cLog,lT
1,24,cLog
y73
cLog10
tP3
cLog2,nE3
cMax,nD3
35847,cMax,nD3
30751,cMax
tP3
cMax,AnyParams,1}
,0,cM
7168,cMin,nD3
35847,cMin,nD3
30751,cMin
tP3
cMin,AnyParams,1}
,0,cM
24591,cMin
y73
xA2
7,xA2
92,xA2
93,xA2
120,xA2
149,xA2
231,cSin,lA
0x5
iA2
246,xA2
255,xA2
254,xA2
0,cSin,t4
1,273,cSin,lA
0x1
iA2
214,yK2
231,cSinh,lA
0x5
iA2
246,yK2
254,yK2
255,yK2
464,yK2
0,cSinh,t4
1,0,yK2
15,cSqrt
y73
cP2
0,cTan,t4
1,116,cTan,t4
1,117,cP2
231,cP2
246,cP2
273,cP2
254,cP2
255,cP2
0,yL2
0,cTanh,t4
1,213,yL2
231,yL2
246,yL2
254,yL2
255,yL2
0,cTrunc,nD3
15384,cSub
iX3
15384,cDiv
iX3
476626,cDiv
iX3
122937,iF2
nE3
iF2
lA
tR3
7168,iF2
nD3
31744,iF2
lA
0x20
nN
31751,iF2
lA
0x24
nN
31751,iF2
nD3
122937,i11,nE3
cLess,lA
tR3
41984,cLess,lA
cM
41984,cLess,nD3
7,cLess,nE3
cLessOrEq,nD3
296182,cLessOrEq,nE3
eZ2
lA
tR3
41984,eZ2
lA
cM
41984,eZ2
nD3
7,eZ2
nE3
yB
nD3
296182,cGreaterOrEq
nR
0
cX1
245
cX1
7
cX1
550
cX1
553
cX1
554
cX1
556
cX1
31
cX1
559
cX1
15
cX1
560,cNot,nD3
7706,nF3
7168,nF3
35847,nF3
30751,nF3
463903,nF3
466975,cAnd,iQ
0,0,cAnd,nK
yJ2
eI3
7706,eI3
35847,eI3
463903,eI3
466975,eI3
30751,cOr,iQ
1,0,n02
92,n02
131,n02
245,n02
215,n02
246,cDeg
nR
246,cRad,nE3
cAbsAnd,l6
7168,cAbsOr,iQ
1,0,cY3
tP3
cAbsNotNot,l0
3,32513024,cM3
lA
0x0}
,}
;}
lQ3
lC1{const
Rule
grammar_rules[262]={{ProduceNewTree,17
eJ3
0,cAbs,t02
409,{1,146,cAtan,t02
403
nU
1324
l24,t02
405
nU
307201
l24
eP
253174
nU
255224
l24
eP
259324
nU
257274
l24,t02
152,{1,252,cCeil
eY
486,{1,68,iG2
482,{1,123,iG2
483,{1,125,iG2
151,{1,126,iG2
419,{1,124,iG2
0,{1,403,cCos,l2
2,1,246,{1,252,cCos,l2
18
eJ3
400,iG2
301,{1,406,cCosh,l2
2,1,246,{1,252,cCosh,l2
18
eJ3
400,cCosh
eY
458,{1,122,cFloor,t02
150,{1,252,cFloor,tS3
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
0,3,32542231,{3,36732440,cIf,cK3
573,{3,32513026,cIf,cK3
515,{3,455505423,cIf,cK3
515,{3,433506837,cIf
eY
78,{1,256,tT3
69,{1,258,tT3
404,{1,72,tT3
159,{1,147,cLog,l2
0,1,0
nU
487425,cMax
l3
16,1,445
nU
c73
cMax
l3
0,1,0
nU
483329,cMin
l3
16,1,446
nU
c73
cMin,c6
0,1,153
nU
24832,cPow,tS3
153
nU
25854,cPow,tS3
154
nU
130063
iA
32055
iA
32056
iA
32057
c83
166288
nU
32137
iA
33082
c83
7168
nU
12688
c83
7434
nU
12553
eW2
435
nU
46146
eW2
436
nU
46154
eW2
437
nU
46150
eW2
169
nU
83983
eW2
168
nU
131106
eW2
175
nU
133154
c93
476160
nU
471055
c93
274432
nU
273423
c93
251904
nU
266274
c93
251904
nU
263186
eW2
171,{1,252,cQ2
421,{1,68,cQ2
151,{1,123,cQ2
419,{1,125,cQ2
170,{1,126,cQ2
482,{1,124,cQ2
0,{1,405,cQ2
172,{1,252,cSinh
eY
328,{1,404,cSinh
eY
173,{1,252,tU3
0,{1,408,tU3
176,{1,410,tU3
177,{1,252,cTanh,l2
0,1,442
nU
449551,i6
1,441
nU
c73
i6
1,167
nU
268549,i6
1,180
nU
276749,i6
1,181
nU
276500
cA3
190770
nU
189622
cA3
194748
nU
193723
cA3
202943
nU
196795
cA3
59699
nU
298148
cA3
59714
nU
325815
cA3
59724
nU
343224
x7
c6
2,1,337,{1,333
eZ
1
tO
336,{1,338
eZ
1}
}
,{ReplaceParams,2,1,340
nU
1363
nQ
342
nU
1365
nQ
463
nU
472524
nQ
47
nU
356711
nQ
349
nU
200751
nQ
360
nU
199727
nQ
480
nU
207053
nQ
481
nU
208077
nQ
417
nU
211144
nQ
209
nU
211145
nQ
418
nU
215240
nQ
212
nU
212329
nQ
204
nU
373097
nQ
211
nU
372944
nQ
217
nU
201944
nQ
221
nU
223448
nQ
367
nU
508329
nQ
219
nU
508126
nQ
224
nU
225705
nQ
223
nU
225776
nQ
365
nU
230825
nQ
426
nU
377057
nQ
497
nU
377054
nQ
497
nU
204201
nQ
426
nU
375280
nQ
224
nU
375006,cAdd
l3
2,2,407781
nU
233698,cAdd
l3
2,2,59763
nU
233842,i6
1,372
nU
1397,lX1
96
nU
24705,lX1
97
nU
24708,lX1
444
nU
449551,lX1
443
nU
c73
lX1
101
nU
102774,lX1
109
nU
107845,lX1
106
nU
104773
i12
0,2,111631
nU
109893
i12
0,2,108559
nU
110917,lJ
0
tO
113
nU
112658,cMul,SelectedParams,0
tO
567,{1,52,lJ
1
tO
568,{1,42,lJ
1}
}
,{ReplaceParams,2,1,467
nU
45516
x8
356
nU
51555
x8
468
nU
49612
x8
357
nU
47459
x8
429
nU
438699
x8
432
nU
441774
x8
486
nU
498726
x8
494
nU
504870
x8
382
nU
435579
x8
497
nU
435709
x8
426
nU
508287
x8
414
nU
500092
x8
499
nU
352744
x8
345
nU
367092
x8
381
nU
425318
x8
478
nU
425460
x8
47
nU
512501
x8
505
nU
355817
x8
47
nU
516598
x8
507
nU
518182
x8
508
nU
358896
x8
351
nU
388605
x8
511
nU
360939
x8
503
nU
354788
x8
514
nU
525350
x8
510
nU
394343
x8
386
nU
351347
i12
2,2,363004
nU
361968
i12
16,1,118
nU
1157
i12
16,1,119
nU
1158
i12
16,1,402
nU
411024
i12
16,2,58768
nU
1472
i12
16,2,15760
nU
1474
i12
17
eJ3
400
i12
17,1,57,{1,14,lJ
0}
}
,{ProduceNewTree,4,1,538
nU
41,iF2
cB3
0
nU
5167,iF2
cL
41984
nU
409641,iF2
cL
tU
iF2
cL
t5
iF2
cL
t6
cEqual
cY1
24849,cEqual
eP
tV
cEqual
eP
n12
281873,cEqual
eP
iR
cEqual
eP
lI1
iF2
cB3
562
nU
41,i11,cB3
538
nU
5167,i11,cL
41984
nU
409641,i11,cL
tU
i11,cL
t5
i11,cL
t6
i11
cY1
24849,i11
eP
tV
i11
eP
n12
281873,i11
eP
iR
i11
eP
lI1
i11,cL
tU
cI3
t5
cI3
t6
cLess,t02
571
nU
46080,cLess
cY1
24832,cLess
eP
y61
cLess
eP
tV
cLess
eP
n12
cJ3
cLess
eP
nZ1
cLess
eP
iR
cLess
eP
lI1
cLess,cC3
562
nU
409641,cI3
tU
yM2
t5
yM2
t6
cLessOrEq,t02
565
nU
409615,cLessOrEq
cY1
24832,cLessOrEq
eP
y61
cLessOrEq
eP
tV
cLessOrEq
eP
n12
cJ3
cLessOrEq
eP
nZ1
cLessOrEq
eP
iR
cLessOrEq
eP
lI1
cLessOrEq,cC3
562
nU
409647,yM2
tU
cR2
t5
cR2
t6
eZ2
t02
539
nU
409615,cGreater
cY1
24832,cGreater
eP
y61
cGreater
eP
tV
cGreater
eP
n12
cJ3
cGreater
eP
nZ1
cGreater
eP
iR
cGreater
eP
lI1
eZ2
cC3
538
nU
409647,cR2
tU
yB
cL
t5
yB
cL
t6
yB
t02
572
nU
46080,yB
nN2
529654
nU
24832,yB
nN2
y61
yB
nN2
tV
yB
nN2
n12
cJ3
yB
nN2
nZ1
yB
nN2
iR
yB
nN2
lI1
yB
cC3
538
nU
409641,yB
cB3
519,{1,137,cNot,cK3
571,{1,2,cNot,l2
0,1,452
nU
c73
nG3
0,2,537097,{3,547892744,cAnd,c6
16,1,566,{1,5,cAnd,AnyParams,1}
}
,{ReplaceParams,16,1,569
nU
13314,nG3
16,1,544
nU
553498,nG3
16,1,546
nU
462369,nG3
16,1,548
nU
466465,nG3
0,1,457
nU
c73
y71
570
nU
13314,y71
563
nU
8197,y71
541
nU
553498,y71
542
nU
462369,y71
543
nU
466465,y71
564
nU
143365,cOr,c6
4,1,525,{1,137,cL3
cK3
572,{1,2,cL3
l4
17
eJ3
0,cL3
t02
537,{1,256,cAbsNotNot,c6
18,1,531,{1,254,cAbsNotNot,c6
0,1,572,{3,43039744,cM3
tS3
571,{3,49325056,cM3
cK3
454,{3,32513586,cM3
l2
16,3,32542225,{3,36732434,cM3
yN1}
,}
;eI2
grammar_optimize_abslogical_type{y4
9
cT
grammar_optimize_abslogical_type
grammar_optimize_abslogical={9,{34,192,228,238,242,247,254,260,261}
}
;}
eI2
grammar_optimize_ignore_if_sideeffects_type{y4
59
cT
grammar_optimize_ignore_if_sideeffects_type
grammar_optimize_ignore_if_sideeffects={59,{0,20,21,22,23,24,25,26,cR
iQ1
78,cS
cU
eI2
grammar_optimize_nonshortcut_logical_evaluation_type{y4
56
cT
grammar_optimize_nonshortcut_logical_evaluation_type
grammar_optimize_nonshortcut_logical_evaluation={56,{0,25,cR
iQ1
78,cS
241,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259}
}
;}
eI2
grammar_optimize_recreate_type{y4
22
cT
grammar_optimize_recreate_type
grammar_optimize_recreate={22,{18,55,56,57,80,81,82,83,84,85,117,118,120,121,130,131,132,133,134,135,136,137}
}
;}
eI2
grammar_optimize_round1_type{y4
125
cT
grammar_optimize_round1_type
grammar_optimize_round1={125,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,25,cR
37,38,iQ1
45,46,47,48,49,50,51,52,53,54,58,59,60,61,62,63,64,65,66,67,68,69,70,71,78,79,80,81,82,83,84,85,86,87,88,93,94,95,96,97,98,99,100,101,117,118,119,120,121,122,123,124,125,126,127,128,129,138,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cU
eI2
grammar_optimize_round2_type{y4
103
cT
grammar_optimize_round2_type
grammar_optimize_round2={103,{0,15,16,17,25,cR
39,40,iQ1
45,46,47,48,49,50,51,52,53,54,59,60,72,73,78,79,86,87,88,89,90,91,92,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,119,122,123,124,125,126,127,128,139,159,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cU
eI2
grammar_optimize_round3_type{y4
79
cT
grammar_optimize_round3_type
grammar_optimize_round3={79,{74,75,76,77,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,170,171,172,173,174,175,176,177,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,197,198,199,201,202,203,205,206,207,208,209,210,211,213,214,215,217,218,219,220,221,222,223,225,226,227,229,230,231,232,233,234,235}
}
;}
eI2
grammar_optimize_round4_type{y4
12
cT
grammar_optimize_round4_type
grammar_optimize_round4={12,{18,55,56,57,130,131,132,133,134,135,136,137}
}
;}
eI2
grammar_optimize_shortcut_logical_evaluation_type{y4
53
cT
grammar_optimize_shortcut_logical_evaluation_type
grammar_optimize_shortcut_logical_evaluation={53,{0,25,cR
iQ1
78,cS
cU}
lQ3
lC1{tI1
eL2
eN1
nY2
paramlist,lK1){index=(paramlist>>(index*10))&1023;if(index>=57)return
eL2(SubFunction,(xB2
plist_s[index-57]yG3
index>=37)return
eL2(NumConstant,(xB2
plist_n_container
xE::plist_n[index-37])n21
eL2(ParamHolder,(xB2
plist_p[index]);}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <stdio.h>
#include <algorithm>
#include <map>
#include <sstream>
using
lQ3
FUNCTIONPARSERTYPES;using
lQ3
lC1;using
tE;using
nE1;lQ3{nS1
It,xD3
T,xD3
Comp>t61
MyEqualRange(It
first,It
last,const
T&val,Comp
comp)iA3
len=last-first;while(len>0)iA3
y03
len/2;It
x83(first);x83+=half;if(comp(*x83,val)){first=x83;++first;len=len-half-1;}
iF1
comp(val,*x83)){len=half;}
else{It
left(first);{It&t12=left;It
last2(x83);size_t
len2=last2-t12;while(len2>0)iA3
half2=len2/2;It
cH3(t12);cH3+=half2;if(comp(*cH3,val)){t12=cH3;++t12;len2=len2-half2-1;}
else
len2=half2;}
}
first+=len;It
right(++x83);{It&t12=right;It&last2=first;size_t
len2=last2-t12;while(len2>0)iA3
half2=len2/2;It
cH3(t12);cH3+=half2;if(comp(val,*cH3))len2=half2;else{t12=cH3;++t12;len2=len2-half2-1;}
t42
t61(left,right);t42
t61(first,first);}
tI1
eI2
OpcodeRuleCompare{iY2()(cI2&tree,unsigned
yN2)const{const
Rule&rule=grammar_rules[yN2]n21
i42<rule
cS2.subfunc_opcode;}
iY2()nY2
yN2,const
eU
const{const
Rule&rule=grammar_rules[yN2]n21
rule
cS2.subfunc_opcode<i42;}
}
yF1
bool
TestRuleAndApplyIfMatch(iK1
lT1&tree,bool
c7{MatchInfo
xE
info;n51
found(false,cZ()yG3(rule.lJ1
LogicalContextOnly)&&!c7{tN1
if(nA
IsIntType
xE::c53){if(rule.lJ1
NotForIntegers)tN1
e13
rule.lJ1
OnlyForIntegers)tN1
if(nA
IsComplexType
xE::c53){if(rule.lJ1
NotForComplex)tN1
e13
rule.lJ1
OnlyForComplex)tN1
for(;;){
#ifdef DEBUG_SUBSTITUTIONS
#endif
found=TestParams(rule
cS2,tree,found.specs,info,true
yG3
found.found)break;if(!&*found.specs){fail:;
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
yW3,false);
#endif
return
tI3}
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
yW3,true);
#endif
SynthesizeRule(rule
yW3
lQ2}
nE1{yG1
ApplyGrammar(const
Grammar&iH2,lT1&tree,bool
c7{if(tree.GetOptimizedUsing()==&iH2){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Already optimized:  "
iQ3
tree)xI1"\n"
<<std::flush;
#endif
return
tI3
if(true){bool
changed=false;switch
tT2{case
cNot:case
cNotNot:case
cAnd:case
cOr:for
lM2
0;nH3
true))c31
lC
cIf:case
cAbsIf:if(ApplyGrammar(iH2,eV,i42==cIf))c31
for
lM2
1;nH3
c7)c31
break;yT3
for
lM2
0;nH3
false))c31}
if(changed){tree.Mark_Incompletely_Hashed(lQ2}
typedef
const
unsigned
short*nI3;std::pair<nI3,nI3>range=MyEqualRange(iH2.rule_list,iH2.rule_list+iH2.rule_count,tree,OpcodeRuleCompare
xE());std
xT3<unsigned
short>rules;rules.xP3
range.second-range.first);for
xW
if(IsLogisticallyPlausibleParamsMatch(eE1
cS2,tree))rules.push_back(*r);}
range.first=!rules.empty()?&rules[0]:0;range.second=!rules.empty()?&rules[rules
t13-1]+1:0;if(range.first!=range.second){
#ifdef DEBUG_SUBSTITUTIONS
if(range.first!=range.second)t72"Input ("
<<t53
tT2<<")["
<<eR<<"]"
;if(c7
std::cout<<"(Logical)"
tJ3
first=iR1,prev=iR1;xN3
sep=", rules "
;for
xW
if(first==iR1)first=prev=*r;iF1*r==prev+1)prev=*r;else
t72
sep<<first;sep=","
;if(prev!=first)std::cout<<'-'<<prev;first=prev=*r;}
}
if(first!=iR1)t72
sep<<first;if(prev!=first)std::cout<<'-'<<prev;}
std::cout<<": "
iQ3
tree)xI1"\n"
<<std::flush;}
#endif
bool
changed=false;for
xW
#ifndef DEBUG_SUBSTITUTIONS
if(!IsLogisticallyPlausibleParamsMatch(eE1
cS2,tree)cF2;
#endif
if(TestRuleAndApplyIfMatch(eE1,tree,c7){c31
y53}
if(changed){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Changed."
<<std::endl
xI1"Output: "
iQ3
tree)xI1"\n"
<<std::flush;
#endif
tree.Mark_Incompletely_Hashed(lQ2}
tree.SetOptimizedUsing(&iH2)n21
tI3
yG1
ApplyGrammar(cH2
p,FPoptimizer_CodeTree::eU
yQ
ApplyGrammar(*(const
Grammar*)p,tree);}
eR1
ApplyGrammars(FPoptimizer_CodeTree::eU{
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round1\n"
;
#endif
n6
grammar_optimize_round1
lA3
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round2\n"
;
#endif
n6
grammar_optimize_round2
lA3
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round3\n"
;
#endif
n6
grammar_optimize_round3
lA3
#ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_nonshortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_nonshortcut_logical_evaluation
lA3
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round4\n"
;
#endif
n6
grammar_optimize_round4
lA3
#ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_shortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_shortcut_logical_evaluation
lA3
#endif
#ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_ignore_if_sideeffects\n"
;
#endif
n6
grammar_optimize_ignore_if_sideeffects
lA3
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_abslogical\n"
;
#endif
n6
grammar_optimize_abslogical
lA3
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
lQ3
FUNCTIONPARSERTYPES;using
lQ3
lC1;using
tE;using
nE1;lQ3{yG1
TestImmedConstraints
nY2
bitmask,const
eU{switch(bitmask&ValueMask){case
Value_AnyNum:case
ValueMask:lC
x92:if(GetEvennessInfo
cE3
n22
Value_OddInt:if(GetEvennessInfo
cE3
xC2
tX1:if(GetIntegerInfo
cE3
n22
Value_NonInteger:if(GetIntegerInfo
cE3
xC2
t41:if(!IsLogicalValue(tree)cG
nK1
SignMask){case
Sign_AnySign:lC
nJ1:if(l91
n22
t51:if(l91
xC2
Sign_NoIdea:if(l91
Unknown
cG
nK1
OnenessMask){case
Oneness_Any:case
OnenessMask:lC
Oneness_One:if(!nO2
if(!fp_equal(fp_abs(l14)tO1
1))cG
lC
Oneness_NotOne:if(!nO2
if(fp_equal(fp_abs(l14)tO1
1))cG
nK1
ConstnessMask){case
Constness_Any:lC
tW1:if(!nO2
lC
Constness_NotConst:if(nO2
y53
return
true;}
tM1
unsigned
extent,unsigned
nbits,xD3
t22=unsigned
int>eI2
nbitmap{private:static
const
unsigned
bits_in_char=8;static
const
unsigned
t32=(cN3
t22)*bits_in_char)/nbits;t22
data[(extent+t32-1)/t32];eS3
void
inc(lK1,int
by=1){data[pos(index)]+=by*t22(1<<yO2);x12
void
dec(lK1){inc(index,-1);}
int
get(lK1
cB1(data[pos(index)]>>yO2)&mask()xX3
pos(lK1)yQ
index/t32
xX3
shift(lK1)yQ
nbits*(index%t32)xX3
mask()yQ(1<<nbits)-1
xX3
mask(lK1)yQ
mask()<<yO2;}
}
;eI2
eC3{int
SubTrees:8;int
Others:8;int
lN2:8;int
tC3:8;nbitmap<l93,2>SubTreesDetail;eC3(l34
memset(this,0,cN3*this));}
eC3(const
eC3&b
l34
memcpy(this,&b,cN3
b));}
eC3&eQ1=(const
eC3&b
l34
memcpy(this,&b,cN3
b))n21*this;}
}
yF1
eC3
CreateNeedList_uncached(t7&iI2{eC3
yH1;for(t71=0;a<tD3
yP2;++a){const
eL2&nJ3=eN1
xE(tD3.param_list,a);t33
SubFunction:{cQ
nB3
cQ*t11;yH
GroupFunction)++yH1.tC3;else{++iK2;assert(param.data.subfunc_opcode<VarBegin);yH1.SubTreesDetail.inc(param.lH1);}
++yH1.lN2;y53
case
t63
case
ParamHolder:++iJ2;++yH1.lN2;break;t42
yH1;}
tI1
eC3&CreateNeedList(t7&iI2{typedef
std::map<t7*,eC3>eF1;static
eF1
c51;eF1::yE3
i=c51.y12&iI2;if(i!=c51.e71&iI2
return
i
eN2
n21
c51.y93,std::make_pair(&tD3,CreateNeedList_uncached
xE(iI2))eN2;}
tI1
lT1
CalculateGroupFunction(const
eL2&yQ2
const
l74){t33
t63{const
ParamSpec_NumConstant
xE
nB3
const
ParamSpec_NumConstant
xE*t11
n21
CodeTreeImmed(param.constvalue)eL3
ParamHolder:{cP
nB3
cP*t11
n21
tM3
GetParamHolderValueIfFound(param.index)eL3
SubFunction:{cQ
nB3
cQ*t11
i2
c53;c53
iD
param.lH1);c23
l52).reserve
eP2
yP2);for(t71=0;a<param.data
yP2;++a){lT1
tmp(CalculateGroupFunction(eN1
xE
eP2.param_list,a),info));c53
yD1
tmp);}
c53
lP2
n21
c53;t42
lT1();}
}
nE1{yG1
IsLogisticallyPlausibleParamsMatch(t7&tD3,const
eU{eC3
yH1(CreateNeedList
xE(iI2);size_t
tV3=eR;if(tV3<size_t(yH1.lN2))nK3}
for
lM2
0;a<tV3;++a){unsigned
opcode=t01
nC;switch(opcode){case
cImmed:if(yH1.tC3>0)--yH1.tC3;else--iJ2;lC
l93:case
cFCall:case
cPCall:--iJ2;break;yT3
assert(opcode<VarBegin);if(iK2>0&&yH1.SubTreesDetail.get(opcode)>0){--iK2;yH1.SubTreesDetail.dec(opcode);}
else--iJ2;}
}
if(yH1.tC3>0||iK2>0||iJ2>0)nK3}
if(tD3.match_type!=AnyParams){if(0||iK2<0||iJ2<0)nK3
t42
true;}
tI1
n51
TestParam(const
eL2&yQ2
cI2&tree
t23
start_at,l74){t33
t63{const
ParamSpec_NumConstant
xE
nB3
const
ParamSpec_NumConstant
xE*t11;if(!nO2
Value_t
imm=l14;switch(param.modulo){case
Modulo_None:lC
Modulo_Radians:imm=fp_mod(imm,y2
imm<xG1
imm
yO
if(imm>fp_const_pi
xE())imm-=fp_const_twopi
xE(yC2
return
fp_equal(imm,param.constvalue)eL3
ParamHolder:{cP
nB3
cP*t11;if(!x2
return
tM3
SaveOrTestParamHolder(param.index,tree)eL3
SubFunction:{cQ
nB3
cQ*t11;yH
GroupFunction){if(!x2
lT1
y81=CalculateGroupFunction(yQ2
info);
#ifdef DEBUG_SUBSTITUTIONS
DumpHashes(y81)xI1*(cH2*)&y81
y11
xI1"\n"
xI1*(cH2*)&l14
xI1"\n"
;DumpHashes(tree)xI1"Comparing "
iQ3
y81)xI1" and "
iQ3
tree)xI1": "
xI1(y81
xK
tree)?"true"
:"false"
)xI1"\n"
;
#endif
return
y81
xK
tree);}
e13!&*start_at){if(!x2
if(i42!=param.lH1
cG
eX3
TestParams
eP2,tree,start_at,info,false);}
t42
tI3
tI1
eI2
l61
xF2
MatchInfo
xE
info;l61()t73,info(){}
}
yF1
class
MatchPositionSpec_PositionalParams:xQ1
l61
xE>{eS3
lB3
MatchPositionSpec_PositionalParams(cO3):eG1
l61
xE>(n){}
}
;eI2
iS1
xF2
iS1()t73{}
}
;class
yK:xQ1
iS1>{eS3
unsigned
trypos;lB3
yK(cO3):eG1
iS1>(n),trypos(0){}
}
yF1
n51
TestParam_AnyWhere(const
eL2&yQ2
cI2&tree
t23
start_at,l74,xH2&used,bool
t82{xT<yK>x6;xE2
yK
cJ2
a=x6->trypos;goto
retry_anywhere_2;yR2
yK(eR);a=0;}
for(;a<eR;++a){if(used[a]cF2;retry_anywhere
cP3
TestParam(yQ2
t01,(cT2);xD2
used[a]=true
iH
a);x6->trypos=a
n21
n51(true,&*x6);}
}
retry_anywhere_2:n32
goto
retry_anywhere;t42
tI3
tI1
eI2
yS1
xF2
MatchInfo
xE
info;xH2
used;lB3
yS1(size_t
tV3)t73,info(),used(tV3){}
}
yF1
class
MatchPositionSpec_AnyParams:xQ1
yS1
xE>{eS3
lB3
MatchPositionSpec_AnyParams(cO3,size_t
m):eG1
yS1
xE>(n,yS1
xE(m)){}
}
yF1
n51
TestParams(t7&nS,cI2&tree
t23
start_at,l74,bool
t82{if(nS.match_type!=AnyParams){if
lA1!=eR
cG}
if(!IsLogisticallyPlausibleParamsMatch(nS,tree))nK3
t93
nS.match_type){case
PositionalParams:{xT<cJ>x6;xE2
cJ
cJ2
a=nS
yP2-1;goto
lL1;yR2
cJ
lA1);a=0;}
for(;a
eM{(tL3
info=info;retry_positionalparams
cP3
TestParam(cV
a),t01,(cT2);xD2
tK3}
lL1:n32
cU2
a].info;goto
retry_positionalparams;}
if(a>0){--a;goto
lL1;}
cU2
l94
n21
tI3
if(t82
for
xG2
tM3
SaveMatchedParamIndex(a)n21
n51(true,&*x6)eL3
SelectedParams:case
AnyParams:{xT<tF>x6;xH2
used(eR);std
xT3<unsigned>lC3
lA1);std
xT3<unsigned>yS2
lA1);for
xG2{const
eL2
nJ3=cV
a);lC3[a]=ParamSpec_GetDepCode(nJ3);}
{unsigned
b=0;for
xG2
if(lC3[a]!=0)yS2[b++]=a;for
xG2
if(lC3[a]==0)yS2[b++]=a;}
xE2
tF
cJ2
if
lA1==0){a=0;goto
retry_anyparams_4;}
a=nS
yP2-1;goto
eL1;yR2
tF
lA1,eR);a=0;if
lA1!=0){(*x6)[l94=info;(*x6)[0].used=used;}
}
for(;a
eM{if(a>0){(tL3
info=info;(tL3
used=used;}
retry_anyparams
cP3
TestParam_AnyWhere
xE(cV
yS2[a]),tree,(cT2,used,t82;xD2
tK3}
eL1:n32
cU2
a].info;used=(tL3
used;goto
retry_anyparams;}
eM1:if(a>0){--a;goto
eL1;}
cU2
l94
n21
tI3
retry_anyparams_4:if(nS.n1!=0){if(!TopLevel||!tM3
HasRestHolder(nS.n1)){eO
cV2;cV2.xP3
eR);for
nY2
b=0;b<eR;++b){if(cQ3
cF2;cV2.push_back(tree
l8
b));cQ3=true
iH
b);}
if(!tM3
SaveOrTestRestHolder(nS.n1,cV2)){goto
eM1;}
}
else{const
eO&cV2=tM3
GetRestHolderValues(nS.n1)lL2
0;a<cV2
tZ2
a){bool
found=false;for
nY2
b=0;b<eR;++b){if(cQ3
cF2;if(cV2[a]xK
tree
l8
b))){cQ3=true
iH
b);found=true;y53}
if(!found){goto
eM1;}
}
t42
n51(true,nS
yP2?&*x6:0)eL3
GroupFunction:y53
return
tI3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
using
tE;using
nE1;lQ3{tI1
lT1
y91
const
eL2&yQ2
l74,bool
inner=true){t33
t63{const
ParamSpec_NumConstant
xE
nB3
const
ParamSpec_NumConstant
xE*t11
n21
CodeTreeImmed(param.constvalue)eL3
ParamHolder:{cP
nB3
cP*t11
n21
tM3
GetParamHolderValue(param.index)eL3
SubFunction:{cQ
nB3
cQ*t11
i2
tree;i52
param.lH1);for(t71=0;a<param.data
yP2;++a){lT1
nparam=y91
eN1
xE
eP2.param_list,a),info,true
iO
nparam);}
if
eP2.n1!=0){eO
trees(tM3
GetRestHolderValues
eP2.n1)yV
AddParamsMove(trees
yG3
eR==1){assert(tree.GetOpcode()==cAdd||tree.GetOpcode()==cMul||tree.GetOpcode()==cMin||tree.GetOpcode()==cMax||tree.GetOpcode()==cAnd||tree.GetOpcode()==cOr||tree.GetOpcode()==cAbsAnd||tree.GetOpcode()==cAbsOr);tree
yO3
eV);}
iF1
eR==0){switch
tT2{case
cAdd:case
cOr:tree=nD1
0));lC
cMul:case
cAnd:tree=nD1
1));yT3
y53}
}
if(inner)tree
lP2
n21
tree;t42
lT1();}
}
nE1{eR1
SynthesizeRule(iK1
lT1&tree,l74){switch(rule.ruletype){case
ProduceNewTree:{tree
yO3
y91
eN1
lV1
0),info,false)yC2
case
ReplaceParams:yT3{std
xT3<unsigned>list=tM3
GetMatchedParamIndexes();std::sort(list.l33
list
iW3)lL2
list
t13;a-->0;)tree.DelParam(list[a]);for(t71=0;a<rule.repl_param_count;++a){lT1
nparam=y91
eN1
lV1
a),info,true
iO
nparam);}
y53}
}
}
#endif
#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
#include <cstring>
using
lQ3
FUNCTIONPARSERTYPES;using
lQ3
lC1;using
tE;using
nE1;lQ3
lC1{eR1
DumpMatch(iK1
cI2&tree,const
l74,bool
DidMatch,std::ostream&o){DumpMatch(rule
yW3,DidMatch?iZ3"match"
:iZ3"mismatch"
,o);}
eR1
DumpMatch(iK1
cI2&tree,const
l74,xN3
tW3,std::ostream&o){static
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
;o<<tW3<<" (rule "
<<(&rule-grammar_rules)<<")"
<<":\n  Pattern    : "
;{eL2
tmp;tmp.first=SubFunction;ParamSpec_SubFunction
tmp2;tmp2.data=rule
cS2;tmp.second=(xB2
tmp2;DumpParam
xE(tmp,o);}
o<<"\n  Replacement: "
;DumpParams
lV1
rule.repl_param_count
l83
o<<"  Tree       : "
iQ3
tree
l83
if(!std::strcmp(tW3,iZ3"match"
))DumpHashes(tree,o)lL2
0;a<tM3
yZ
tZ2
a){if(!tM3
yZ[a].c92)cF2;o<<"           "
<<ParamHolderNames[a]<<" = "
iQ3
tM3
yZ[a]l83}
for(yL3
tM3
lQ
tZ2
b){if(!l84
eZ3
cF2
lL2
0;a<l84
t03
tZ2
a){o<<"         <"
<<b<<"> = "
iQ3
l84
t03[a],o);o<<std::endl;}
}
o<<std::flush;}
}
#endif
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;lQ3{yG1
MarkIncompletes(FPoptimizer_CodeTree::eU{if(tree.Is_Incompletely_Hashed(iY1
bool
iT1=false;for
xH
iT1|=MarkIncompletes(t01
yG3
iT1)tree.Mark_Incompletely_Hashed()n21
iT1;}
eR1
FixIncompletes(FPoptimizer_CodeTree::eU{if(tree.Is_Incompletely_Hashed()){for
xH
FixIncompletes(t01
yV
Rehash();}
}
}
tE{lB
Sort()iE3
Sort();}
lB
Rehash(bool
constantfolding){if(constantfolding)ConstantFolding(*this);else
Sort();data
xA
tI1
eI2
cC{cR3
Value_t
cS3
yT1=0;
#if 0
long
double
value=Value;eA=crc32::calc((const
unsigned
char*)&value,cN3
value));key^=(key<<24);
#elif 0
union{eI2{unsigned
char
filler1[16]yX3
v
tJ3
char
filler2[16];}
buf2;eI2{unsigned
char
filler3[cN3
Value_t)+16-c
N3
xD1)];eA;}
buf1;}
data;memset(&data,0,cN3
data));data.buf2.v=Value;eA=data.buf1.key;
#else
int
eE2
yX3
nP2=std::frexp(Value,&eE2);eA=nY2(eE2+0x8000)&0xFFFF
yG3
nP2<0){nP2=-nP2;key=key^0xFFFF;}
else
key+=0x10000;nP2-=eF2
0.5);key<<=39;key|=n71(nP2+nP2)*eF2
1u<<31))<<8;
#endif
lP
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nS1
T
cW2
std::complex<T> >{cR3
std::complex<T>cS3
cC<T>::nL3
cX2,Value.real());nA
fphash_t
temp;cC<T>::nL3
temp,Value.imag());yT1^=temp.hash2;cX2.hash2^=temp.hash1;}
}
;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
tM1
cW2
long>{yL
long
Value){eA=Value;lP
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
tM1
cW2
GmpInt>{cR3
GmpInt
cS3
eA=Value.toInt();lP
#endif
eR1
xQ2
xE::Recalculate_Hash_NoRecursion(){fphash_t
cX2(n71
Opcode)<<56,Opcode*iM3(0x1131462E270012B));Depth=1;switch(Opcode){case
cImmed:{cC
xE::nL3
cX2,Value
yC2
case
l93:{yT1|=n71
cZ1<<48
eH1((n71
cZ1)*11)^iM3(0x3A83A83A83A83A0);y53
case
cFCall:case
cPCall:{yT1|=n71
cZ1<<48
eH1((~n71
cZ1)*7)^3456789;}
yT3
iA3
t81=0
lL2
0;a<eD3
tZ2
a){if(eD3[a]nM2>t81)t81=eD3[a]nM2;yT1+=((eD3[a]i22
hash1*(a+1))>>12)eH1
eD3[a]i22
hash1
eH1(3)*iM3(0x9ABCD801357);cX2.hash2*=iM3(0xECADB912345)eH1(~eD3[a]i22
hash2)^4567890;}
Depth+=t81;}
}
if(Hash!=cX2){Hash=cX2;l82=0;}
}
lB
FixIncompleteHashes(){MarkIncompletes(*this);FixIncompletes(*this);}
}
#endif
#include <cmath>
#include <list>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;lQ3{using
tE
yF1
bool
x01
cI2&tree,long
count,const
xS1::SequenceOpCode
xE&eT,xS1::nV2
xE&synth,size_t
max_bytecode_grow_length);static
const
eI2
SinCosTanDataType{OPCODE
whichopcode;OPCODE
inverse_opcode;enum{nominator,nQ2,inverse_nominator,lM1}
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
,{cY2{cSinh,cCosh,cZ2,{cSinh,cNop,{cY2
cNop,cCosh}
}
,{cCosh,cNop,{cSinh,cY2
cNop}
}
,{cNop,cTanh,{cCosh,cSinh,cZ2,{cNop,cSinh,{cNop,cTanh,cCosh,cNop}
}
,{cNop,cCosh,{cTanh,cSinh,cZ2}
;}
tE{lB
SynthesizeByteCode(std
xT3<unsigned>&ByteCode,std
xT3
xE&Immed,size_t&stacktop_max){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Making bytecode for:\n"
;l21
#endif
while(RecreateInversionsAndNegations()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"One change issued, produced:\n"
;l21
#endif
FixIncompleteHashes();using
nE1;using
lQ3
lC1;cH2
g=(xB2
grammar_optimize_recreate;while(ApplyGrammar(*(const
Grammar*)g,*this)){FixIncompleteHashes();}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Actually synthesizing, after recreating inv/neg:\n"
;l21
#endif
xS1::nV2
xE
synth;SynthesizeByteCode(synth,false);synth.Pull(ByteCode,Immed,stacktop_max);}
lB
SynthesizeByteCode(xS1::nV2
xE&synth,bool
MustPopTemps)const{yA1*this))yQ;}
for
lM2
0;a<12;++a){const
SinCosTanDataType&data=SinCosTanData[a];if(data.whichopcode!=cNop)i21!=data.whichopcode
cF2;lU1
nM3;nM3.x51);nM3
iD
data.inverse_opcode);nM3.yT2);yA1
nM3)){synth.lB1
else
i21!=cInv
cF2;if(GetParam(0)nC!=data.inverse_opcode
cF2;yA1
GetParam(0))){synth.lB1
size_t
found[4];for(yL3
4;++b){lU1
tree;if(data.tX3]==cNop){i52
cInv);lU1
nN3;nN3.x51);nN3
iD
data.tX3^2]);nN3.yT2
iO
nN3);}
else{tree.x51
yV
SetOpcode(data.tX3]);}
tree.yT2);found[b]tY3
xR3
tree);}
if(found[data.yU2!=tP
nQ2]!=l71
yU2
tX2
nQ2
iL
cDiv
nL1
yU2!=tP
lM1]!=l71
yU2
tX2
lM1
iL
cMul
nL1
lY1!=tP
lM1]!=l71
lY1
tX2
lM1
iL
cRDiv
nL1
lY1!=tP
nQ2]!=l71
lY1
tX2
nQ2
iL
cMul,2,1);synth.lB1
size_t
n_subexpressions_synthesized=SynthCommonSubExpressions(tN2
switch(iU1{case
l93:synth.PushVar(GetVar());lC
cImmed:t91
GetImmed());lC
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:i21==cMul){bool
cT3=false;yM
lZ1
y01&&isLongInteger(lZ1
y11)){cC1=makeLongInteger(lZ1
y11);lU1
tmp(*this,xD3
lU1::CloneTag(tE3.yA2
tmp
lP2;if(x01
tmp,value,xS1::tP1
xE::AddSequence,synth,MAX_MULI_BYTECODE_LENGTH)){cT3=true;y53}
}
if(cT3)y53
int
yU1=0;xH2
done(c61,false);lU1
iS;iS
iD
iU1;for(;;){bool
found=false;yM
done[a]cF2;if(synth.IsStackTop(lZ1)){found=true;done[a]=true;lZ1.n7
iS
e9
lZ1
yG3++yU1>1){yN
2);iS.yT2)y21
iS);yU1=yU1-2+1;}
}
}
if(!found)y53
yM
done[a]cF2;lZ1.n7
iS
e9
lZ1
yG3++yU1>1){yN
2);iS.yT2)y21
iS);yU1=yU1-2+1;}
}
if(yU1==0){switch(iU1{case
cAdd:case
cOr:case
cAbsOr:t91
0);lC
cMul:case
cAnd:case
cAbsAnd:t91
1);lC
cMin:case
cMax:t91
0);break;yT3
y53++yU1;}
assert(n_stacked==1);y53
case
cPow:{tR1
p0
nF2
0);tR1
p1
nF2
1
yG3!p1
y01||!isLongInteger
n91)||!x01
p0,makeLongInteger
n91),xS1::tP1
xE::MulSequence,synth,MAX_POWI_BYTECODE_LENGTH)){p0.n7
p1.n7
yN
2);yC1
cIf:case
cAbsIf:{xD3
xS1::nV2
xE::IfData
ifdata;GetParam(0)xY3
SynthIfStep1(ifdata,iU1;GetParam(1)xY3
SynthIfStep2(ifdata);GetParam(2)xY3
SynthIfStep3(ifdata
yC2
case
cFCall:case
cPCall:{for
lM2
0;a<c61;++a)lZ1.n7
yN
nY2)c61);y43
nJ2|GetFuncNo(),0,0
yC2
yT3{for
lM2
0;a<c61;++a)lZ1.n7
yN
nY2)c61
yC2}
synth.StackTopIs(*this
yG3
MustPopTemps&&n_subexpressions_synthesized>0)iA3
top
tY3
GetStackTop();synth.DoPopNMov(top-1-n_subexpressions_synthesized,top-1);}
}
}
lQ3{yG1
x01
cI2&tree,long
count,const
xS1::SequenceOpCode
xE&eT,xS1::nV2
xE&synth,size_t
max_bytecode_grow_length){if
eF3!=0){xS1::nV2
xE
backup=synth;tree.n7
size_t
bytecodesize_backup
tY3
GetByteCodeSize();xS1::x01
count,eT,tN2
size_t
bytecode_grow_amount
tY3
GetByteCodeSize()-bytecodesize_backup;if(bytecode_grow_amount>max_bytecode_grow_length){synth=backup
n21
false;eX3
true;}
else{xS1::x01
count,eT,synth
lQ2}
}
#endif
#include <cmath>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;lQ3{using
tE;
#define FactorStack std xT3
const
eI2
PowiMuliType{unsigned
opcode_square
tJ3
opcode_cumulate
tJ3
opcode_invert
tJ3
opcode_half
tJ3
opcode_invhalf;}
iseq_powi={cSqr,cMul,cInv,cSqrt,cRSqrt}
,iseq_muli={iR1
x7
cNeg,iR1,iR1}
yF1
Value_t
e01
const
PowiMuliType&cU3,const
std
xT3<unsigned>&eI1,n42&stack
cK2
1);while(IP<limit){if(iV1==cU3.opcode_square){if(!isInteger
i13
2;c2
opcode_invert){c53=-c53;c2
opcode_half){if
e12>xL1&&isEvenInteger
i13
eF2
0.5);c2
opcode_invhalf){if
e12>xL1&&isEvenInteger
i13
eF2-0.5);++IP;tK3
size_t
xI2=IP
yX3
lhs(1
yG3
iV1==cFetch){lK1=nA2;if(index<y1||size_t(index-y1)>=iX1){IP=xI2;y53
lhs=stack[index-y1];goto
yV2;}
if(iV1==cDup){lhs=c53;goto
yV2;yV2:tA1
c53);++IP
yX3
subexponent=e01
cU3
iK
if(IP>=limit||iV1!=cU3.opcode_cumulate){IP=xI2;y53++IP;stack.pop_back();c53+=lhs*subexponent;tK3
y53
return
c53;}
tI1
Value_t
ParsePowiSequence(const
std
xT3<unsigned>&eI1){n42
stack;tA1
eF2
1))n21
e01
iseq_powi
iK}
tI1
Value_t
ParseMuliSequence(const
std
xT3<unsigned>&eI1){n42
stack;tA1
eF2
1))n21
e01
iseq_muli
iK}
tI1
class
CodeTreeParserData{eS3
lB3
CodeTreeParserData(bool
k_powi):stack(),clones(),keep_powi(k_powi){}
void
Eat(size_t
tV3,OPCODE
opcode){lT1
xQ;xQ
iD
opcode);eO
tD3=Pop(tV3
lU3
iI2;if(!keep_powi)switch(opcode){case
cTanh:{lT1
sinh,cosh;sinh
iD
cSinh);sinh
e9
xQ
cV3
sinh
lP2;cosh
iD
cCosh);cosh
yD1
xQ
cV3
cosh
lP2
i2
pow;pow
iD
cPow);pow
yD1
cosh);pow.yE
eF2-1)));pow
lP2;xQ
t83
xQ.nG1
0,sinh);xQ
yD1
pow
yC2
case
cTan:{lT1
sin,cos;sin
iD
cSin);sin
e9
xQ
cV3
sin
lP2;cos
iD
cCos);cos
yD1
xQ
cV3
cos
lP2
i2
pow;pow
iD
cPow);pow
yD1
cos);pow.yE
eF2-1)));pow
lP2;xQ
t83
xQ.nG1
0,sin);xQ
yD1
pow
yC2
case
cPow:{cI2&p0=xQ
l8
0);cI2&p1=xQ
l8
1
yG3
p1
nC==cAdd){eO
mulgroup(p1
y41)lL2
0;a<p1
y41;++a){lT1
pow;pow
iD
cPow);pow
e9
p0);pow
e9
p1
l8
a));pow
lP2;mulgroup[a
i23
pow);}
xQ
iD
cMul
lU3
mulgroup);}
y53
yT3
y53
xQ.Rehash(!keep_powi);iW1,false);
#ifdef DEBUG_SUBSTITUTIONS
tC1<<tV3<<", "
<<t53(opcode)<<"->"
<<t53(xQ
nC)<<": "
iK3
xQ
y6
xQ);
#endif
tA1
xQ
lS3
EatFunc(size_t
tV3,OPCODE
eT3
unsigned
funcno)l04
CodeTreeFuncOp
xE(eT3
funcno);eO
tD3=Pop(tV3
lU3
iI2;xQ.yT2);
#ifdef DEBUG_SUBSTITUTIONS
tC1<<tV3<<", "
iK3
xQ
y6
xQ);
#endif
iW1);tA1
xQ
lS3
AddConst(yZ1)l04
CodeTreeImmed(value);iW1);Push(xQ
lS3
AddVar
nY2
varno)l04
CodeTreeVar
xE(varno);iW1);Push(xQ
lS3
xL{tB1
1
i23
tB1
2]lS3
Dup(){Fetch(iX1-1
lS3
Fetch(size_t
which){Push(stack[which]);}
nS1
T>void
Push(T
tree){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<iK3
tree
y6
tree);
#endif
tA1
tree
lS3
PopNMov(size_t
target,size_t
source){stack[target]=stack[source];stack
xK3
target+1);}
lT1
yW2{clones.clear()i2
c53(stack.back());stack
xK3
iX1-1)n21
c53;}
eO
Pop(size_t
n_pop){eO
c53(n_pop);for
nY2
n=0;n<n_pop;++n)c53[n
i23
tB1
n_pop+n]);
#ifdef DEBUG_SUBSTITUTIONS
for(cO3=n_pop;n-->0;){tC1;DumpTree
e12[n]y6
c53[n]);}
#endif
stack
xK3
iX1-n_pop)n21
c53;}
size_t
GetStackTop(cB1
iX1;}
private:void
FindClone(lT1&,bool=true)yQ;}
private:eO
stack;std::multimap<fphash_t,lT1>clones;bool
keep_powi;private:CodeTreeParserData(const
CodeTreeParserData&);CodeTreeParserData&eQ1=(const
CodeTreeParserData&);}
yF1
eI2
IfInfo{lT1
t92
i2
thenbranch;size_t
endif_location;IfInfo():t92(),thenbranch(),endif_location(){}
}
;}
tE{lB
GenerateFrom(const
xD3
FunctionParserBase
xE::Data&x93,bool
keep_powi){eO
xO2;xO2.xP3
x93.mVariablesAmount);for
nY2
n=0;n<x93.mVariablesAmount;++n){xO2.push_back(CodeTreeVar
xE(n+l93));}
GenerateFrom(x93,xO2,keep_powi);}
lB
GenerateFrom(const
xD3
FunctionParserBase
xE::Data&x93,const
x3&xO2,bool
keep_powi){const
std
xT3<unsigned>&ByteCode=x93.mByteCode;const
std
xT3
xE&Immed=x93.mImmed;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"ENTERS GenerateFrom()\n"
;
#endif
CodeTreeParserData
xE
sim(keep_powi);std
xT3<IfInfo
xE>eN;for(size_t
IP=0,DP=0;;++IP){iL2:while(!eN.empty()&&(eN.eB==IP||(IP<n61&&iV1==cJump&&eN.eO1.c92)))){lU1
elsebranch=sim.yW2
n52
eN.back().t92)n52
eN.eO1)n52
elsebranch)nV
3,cIf);eN.pop_back();}
if(IP>=n61)break
tJ3
opcode=iV1;if((opcode==cSqr||opcode==cDup||(opcode==cInv&&!IsIntType
xE::c53)||opcode==cNeg||opcode==cSqrt||opcode==cRSqrt||opcode==cFetch))iA3
was_ip=IP
yX3
eE2=ParsePowiSequence
xE(ByteCode,IP,eN.empty()?n61:eN.eB,sim.xG
1
yG3
eE2!=eF2
1.0)){xF
eE2)n92;goto
iL2;}
if(opcode==cDup||opcode==cFetch||opcode==cNeg
yF2
y02=ParseMuliSequence
xE(ByteCode,IP,eN.empty()?n61:eN.eB,sim.xG
1
yG3
y02!=eF2
1.0)){xF
y02)nV
2,cMul);goto
iL2;}
}
IP=was_ip;}
if(n62>=l93){lK1=opcode-l93
n52
xO2[index]);}
else{switch(n62){case
cIf:case
cAbsIf:{eN
xK3
eN
t13+1);lU1
res(sim.yW2);eN.back().t92.swap(res);eN.eB=n61;IP+=2;continue
eL3
cJump:{lU1
res(sim.yW2);eN.eO1.swap(res);eN.eB=ByteCode[IP+1]+1;IP+=2;continue
eL3
cImmed:xF
Immed[DP++]);lC
cDup:sim.Dup();lC
cNop:lC
cFCall:{unsigned
funcno=nA2;assert(funcno<fpdata.mFuncPtrs.size())tJ3
tD3=x93.mFuncPtrs[funcno].mParams;sim.EatFunc(tD3,n62,funcno
yC2
case
cPCall:{unsigned
funcno=nA2;assert(funcno<fpdata.iR3.size());const
FunctionParserBase
xE&p=*x93.iR3[funcno].mParserPtr
tJ3
tD3=x93.iR3[funcno].mParams;x3
paramlist=sim.Pop(iI2;lU1
iM2;iM2.GenerateFrom(*p.mData,paramlist)n52
iM2
yC2
case
cInv:xF
1
n82
cDiv);lC
cNeg
n72
cNeg);break;xF
0
n82
cSub);lC
cSqr:xF
2)n92;lC
cSqrt:xF
cI1
lD3
cRSqrt:xF
eF2-0.5))n92;lC
cCbrt:xF
eF2
1)/eF2
3))n92;lC
cDeg:xF
fp_const_rad_to_deg
xM1
cRad:xF
fp_const_deg_to_rad
xM1
cExp:iT)goto
nP3;xF
fp_const_e
xE()n82
lD3
cExp2:iT)goto
nP3;xF
2.0
n82
lD3
cCot
n72
cTan);iT)yF
cCsc
n72
cSin);iT)yF
cSec
n72
cCos);iT)yF
cInt:
#ifndef __x86_64
iT)nO3
1,cInt
yC2
#endif
xF
cI1
cAdd)nV
1,cFloor);lC
cLog10
n72
cW3
fp_const_log10inv
xM1
cLog2
n72
cW3
fp_const_log2inv
xM1
e63:yX2
cW3
fp_const_log2inv
xE())nV
3,cMul);lC
cHypot:xF
2)n92;i33
xF
2)n92
nV
2,cAdd);xF
cI1
lD3
cSinCos:sim.Dup()nV
1,cSin);yX2
cCos);lC
cSinhCosh:sim.Dup()nV
1,cSinh);yX2
cCosh);lC
cRSub:i33
case
cSub:iT)nO3
2,cSub
yC2
xF-1)nV
2,cMul)nV
2,cAdd);lC
cRDiv:i33
case
cDiv:iT||IsIntType
xE::c53)nO3
2,cDiv
yC2
xF-1)n92
nV
2,cMul);lC
cAdd:case
cMul:case
cMod:case
cPow:case
cEqual:case
cLess:case
cGreater:case
i11:case
cLessOrEq:case
cGreaterOrEq:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:sim.Eat(2,tD1
lC
cNot:case
cNotNot:case
cY3:case
cAbsNotNot
n72
tD1
lC
cFetch:sim.Fetch(nA2);lC
cPopNMov:{unsigned
stackOffs_target=nA2
tJ3
stackOffs_source=nA2;sim.PopNMov(stackOffs_target,stackOffs_source
yC2
yT3
nP3:tJ3
funcno=opcode-cAbs;assert(funcno<FUNC_AMOUNT);const
FuncDefinition&func=Functions[funcno]nV
func.tD3,tD1
y53}
}
Become(sim.yW2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Produced tree:\n"
;l21
#endif
}
}
#endif
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
#include <assert.h>
#define FP_MUL_COMBINE_EXPONENTS
lQ3{using
lQ3
FUNCTIONPARSERTYPES;using
tE
yF1
static
void
AdoptChildrenWithSameOpcode(eU{
#ifdef DEBUG_SUBSTITUTIONS
bool
nR2=false;
#endif
for
yW
if(t01
nC==i42){
#ifdef DEBUG_SUBSTITUTIONS
if(!nR2)t72"Before assimilation: "
yJ
nR2=true;}
#endif
tree.AddParamsMove(t01.GetUniqueRef().l52),a);}
#ifdef DEBUG_SUBSTITUTIONS
if(nR2)t72"After assimilation:   "
yJ}
#endif
}
}
tE{eR1
ConstantFolding(eU{tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
void*cX3=0
xI1"["
<<(&cX3)<<"]Runs ConstantFolding for: "
yJ
DumpHashes(tree)xI1
std::flush;
#endif
if(false){redo:;tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cX3)<<"]Re-runs ConstantFolding: "
yJ
DumpHashes(tree);
#endif
}
if(i42!=cImmed){range
xE
p=CalculateResultBoundaries(tree
yG3
p
i4&&p.max
lH3
p
xK1==p
nP){xI
p
xK1);nD}
if(false){ReplaceTreeWithOne:xI
eF2
1));goto
do_return;ReplaceTreeWithZero:xI
xG1;goto
do_return;ReplaceTreeWithParam0:
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before replace: "
xI1
std::hex<<'['
xR1
hash1<<','
xR1
hash2<<']'<<std::dec
yJ
#endif
tree
yO3
eV);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After replace: "
xI1
std::hex<<'['
xR1
hash1<<','
xR1
hash2<<']'<<std::dec
yJ
#endif
eK1}
switch
tT2{case
cImmed:lC
l93:lC
cAnd:case
cAbsAnd:cW
bool
c8=false;for
yW{if(!IsLogicalValue(t01))c8=true;switch(yB1
t01,i42==cAbsAnd)){case
IsNever
e7
tA2:x81);lC
n01
t93
eR){case
0:iB
1:i52
i42==cAnd?cNotNot:cAbsNotNot);eK1
yT3
lE3
cAnd||!c8)if(ConstantFolding_AndLogic(x42
yC1
cOr:case
cAbsOr:cW
bool
c8=false;for
yW{if(!IsLogicalValue(t01))c8=true;switch(yB1
t01,i42==cAbsOr))i61
iB
lI3
x81);lC
n01
t93
eR){case
0
e7
1:i52
i42==cOr?cNotNot:cAbsNotNot);eK1
yT3
lE3
cOr||!c8)if(ConstantFolding_OrLogic(x42
yC1
cNot:case
cY3:{unsigned
n81
0;switch(eV
nC){case
cEqual:n81
i11;lC
i11:n81
cEqual;lC
cLess:n81
cGreaterOrEq;lC
cGreater:n81
cLessOrEq;lC
cLessOrEq:n81
cGreater;lC
cGreaterOrEq:n81
cLess;lC
cNotNot:n81
cNot;lC
cNot:n81
cNotNot;lC
cY3:n81
cAbsNotNot;lC
cAbsNotNot:n81
cY3;break;yT3
y53
if(opposite){i52
OPCODE(opposite)yV
SetParamsMove(eV.GetUniqueRef().l52));eK1
t93
yB1
eV,tree
e11){case
tA2
e7
lI3
iB
n01
lE3
cNot&&GetPositivityInfo(eV)==tA2)i52
cY3
yG3
tF3
cIf||tF3
tG3{lT1
lG3=eV;cI2&ifp1=lG3
l8
1);cI2&ifp2=lG3
l8
2
yG3
ifp1
yY2
ifp1
e11
cZ3
ifp1
nC==cNot?cNotNot:cAbsNotNot);iN2
l8
0))iH1)e03
iG1
p2
i43)e6
if(ifp2
yY2
ifp2
e11
cZ3
i42);iN2)iH1)e03
iD
ifp2
nC==cNot?cNotNot:cAbsNotNot);p2
i43
l8
0))e6
yC1
cNotNot:case
cAbsNotNot:{if(IsLogicalValue(eV))lF3
switch(yB1
eV,i42==cAbsNotNot)){case
IsNever
e7
tA2:iB
n01
lE3
cNotNot&&GetPositivityInfo(eV)==tA2)i52
cAbsNotNot
yG3
tF3
cIf||tF3
tG3{lT1
lG3=eV;cI2&ifp1=lG3
l8
1);cI2&ifp2=lG3
l8
2
yG3
ifp1
yY2
ifp1
e11{tree.SetParam(0,lG3
l8
0)yV
AddParam(ifp1)e03
iG1
p2
i43)e6
if(ifp2
yY2
ifp2
e11
cZ3
i42);iN2)iH1
yV
AddParam(ifp2
yV
SetOpcode(lG3
nC);eK1}
yC1
cIf:case
cAbsIf:{if(ConstantFolding_IfOperations(x42
y53
case
cMul:{NowWeAreMulGroup:;AdoptChildrenWithSameOpcode(tree)yX3
nM1=eF2
1);size_t
l12=0;bool
nN1=false;for
xH{if(!x61
continue
yX3
immed=x71
if(immed==xG1
goto
ReplaceTreeWithZero;nM1*=immed;++l12;}
if(l12>1||(l12==1&&fp_equal(nM1
yF3)nN1=true;if(nN1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cMul: Will add new "
iS3
nM1<<"\n"
;
#endif
for
yW
if(x61{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iS3
x71
std::cout<<"\n"
;
#endif
nQ3!fp_equal(nM1
yF3
tree
e9
eB1
nM1));t93
eR){case
0:iB
1:lF3
yT3
if(ConstantFolding_MulGrouping(x42
if(ConstantFolding_MulLogicItems(x42
yC1
cAdd:cW
Value_t
nB2=0.0;size_t
l12=0;bool
nN1=false;for
xH{if(!x61
continue
yX3
immed=x71
nB2+=immed;++l12;}
if(l12>1||(l12==1&&nB2==xG1)nN1=true;if(nN1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cAdd: Will add new "
iS3
nB2<<"\n"
xI1"In: "
yJ
#endif
for
yW
if(x61{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iS3
x71
std::cout<<"\n"
;
#endif
nQ3!(nB2==eF2
0.0)))tree
e9
eB1
nB2));t93
eR){case
0
e7
1:lF3
yT3
if(ConstantFolding_AddGrouping(x42
if(ConstantFolding_AddLogicItems(x42
yC1
cMin:cW
size_t
yZ2=0;range
xE
e2;for
xH{while(a+1<eR&&t01
xK
tree
l8
a+1)))x81+1);range<nW
max
lH3(!e2.nT3||(p
nP)<e2
nP)){e2
nP=p
nP;e2.nT3=true;yZ2=a;}
}
if(e2.iI1
for
yW{range<nW
min
lH3
a!=yZ2&&p
xK1>=e2
nP)nQ3
eR==1){lF3
yC1
cMax:cW
size_t
yZ2=0;range
xE
t8;for
xH{while(a+1<eR&&t01
xK
tree
l8
a+1)))x81+1);range<nW
min
lH3(!t8
i4||p
xK1>t8
xK1)){t8
xK1=p
xK1;t8
i4=true;yZ2=a;}
}
if(t8
i4){for
yW{range<nW
max
lH3
a!=yZ2&&(p
nP)<t8
xK1){x81);}
}
}
if(eR==1){lF3
yC1
cEqual:case
i11:case
cLess:case
cGreater:case
cLessOrEq:case
cGreaterOrEq:if(ConstantFolding_Comparison(x42
lC
cAbs:{range
xE
tW
eV
yG3
tX
lF3
if(yV1{i52
cMul
yV
yE
eF2
1)));goto
NowWeAreMulGroup;}
if(tF3
cMul){cI2&p=eV;eO
nR3;eO
e02
lL2
0;a<p
y41;++a){tW
p
l8
a)yG3
tX{nR3.push_back(p
l8
a));}
if(yV1{e02.push_back(p
l8
a));}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Abs: mul group has "
<<nR3
t13<<" pos, "
<<e02
t13<<"neg\n"
;
#endif
if(!nR3.empty()||!e02.empty()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-Before: "
iQ3
tree)xI1"\n"
<<std::flush;DumpHashes
l02;
#endif
lT1
eK3;eK3
iD
cMul)lL2
0;a<p
y41;++a){tW
p
l8
a)yG3(tX||(yV1){}
else
eK3
e9
p
l8
a));}
eK3
lP2
i2
nS3;nS3
iD
cAbs);nS3
yD1
eK3);nS3
lP2
i2
iC
cMul);mulgroup
yD1
nS3);yE1
AddParamsMove(nR3
yG3!e02.empty()){if(e02
t13%2)yE1
yE
eF2-1)));yE1
AddParamsMove(e02);}
tree
yO3
mulgroup);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-After: "
;DumpTree
l02
xI1"\n"
<<std::flush;DumpHashes
l02;
#endif
goto
NowWeAreMulGroup;}
}
y53
#define HANDLE_UNARY_CONST_FUNC(funcname) nG){xI funcname(lR));nD
case
cLog:iF3(fp_log);if(tF3
cPow){lT1
pow=eV;if(GetPositivityInfo(pow
l8
0))==tA2)tB3
yI
tree.lU
if(GetEvennessInfo(pow
l8
1))==tA2)tB3)i2
abs;abs
iD
cAbs);abs
yD1
pow
cV3
abs.Rehash(yI
pow.nG1
0,abs
yV
lU}
iF1
tF3
cAbs){lT1
pow=eV
l8
0
yG3
pow
nC==cPow)tB3)i2
abs;abs
iD
cAbs);abs
yD1
pow
cV3
abs.Rehash(yI
pow.nG1
0,abs
yV
lU}
lC
cAcosh:iF3(fp_acosh);lC
cAsinh:iF3(fp_asinh);lC
cAtanh:iF3(fp_atanh);lC
cAcos:iF3(fp_acos);lC
cAsin:iF3(fp_asin);lC
cAtan:iF3(fp_atan);lC
cCosh:iF3(fp_cosh);lC
cSinh:iF3(fp_sinh);lC
cTanh:iF3(fp_tanh);lC
cSin:iF3(fp_sin);lC
cCos:iF3(fp_cos);lC
cTan:iF3(fp_tan);lC
cCeil:if(n5
iF3(fp_ceil);lC
cTrunc:if(n5
iF3(fp_trunc);lC
cFloor:if(n5
iF3(fp_floor);lC
cInt:if(n5
iF3(fp_int);lC
cCbrt:iF3(fp_cbrt);lC
cSqrt:iF3(fp_sqrt);lC
cExp:iF3(fp_exp);lC
cLog2:iF3(fp_log2);lC
cLog10:iF3(fp_log10);lC
e63:if
lI
fp_log2(lR)*tF1
cArg:iF3(fp_arg);lC
cConj:iF3(fp_conj);lC
cImag:iF3(fp_imag);lC
cReal:iF3(fp_real);lC
cPolar:if
lI
fp_polar
x91
lC
cMod:if
lI
fp_mod
x91
lC
cAtan2:{range
xE
tW
eV
yJ3
p1=eX
1));nG&&fp_equal(lR,xG1){if(p1
xJ2
p1
nP)<xG1{xI
fp_const_pi
xE());nD
if(p1
i4&&p1
xK1>=tE1
xG1;nD}
if(l41&&fp_equal(yB2,xG1){if(p0
xJ2
p0
nP)<xG1{xI-fp_const_pihalf
xE());nD
if(p0
i4&&p0
xK1>xG1{xI
fp_const_pihalf
xE());nD}
if
lI
fp_atan2
x91
if((p1
i4&&p1
xK1>xG1||(p1
xJ2
p1
nP)<fp_const_negativezero
xE()nM
c12;c12
iD
cPow);c12
yD1
tV2);c12.yE
eF2-1)));c12
lP2
i2
c22;c22
t83
c22
yD1
eV);c22
yD1
c12);c22
i3
cAtan
yV
nG1
0,c22
iP1
1);yC1
cPow:{if(ConstantFolding_PowOperations(x42
y53
case
cDiv:nG&&l41&&yB2!=tE1
lR/tF1
cInv:nG&&lR!=tE1
eF2
1)/lR);nD
lC
cSub:if
lI
lR-tF1
cNeg:nG){xI-lR);nD
lC
cRad:nG){xI
RadiansToDegrees
tY2
cDeg:nG){xI
DegreesToRadians
tY2
cSqr:nG){xI
lR*lR);nD
lC
cExp2:iF3(fp_exp2);lC
cRSqrt:nG){xI
eF2
1)/fp_sqrt
tY2
cCot:xK2
fp_tan(lZ
cSec:xK2
fp_cos(lZ
cCsc:xK2
fp_sin(lZ
cHypot:if
lI
fp_hypot
x91
lC
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:lC
cPCall:case
cFCall:y53
do_return:;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cX3)<<"]Done ConstantFolding, result: "
yJ
DumpHashes(tree);
#endif
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
tE{eR1
range
xE::set_abs(nL
bool
has_negative=!min.known||min.val<eF2);bool
has_positive=!nT3||i53>eF2);bool
crosses_axis=has_negative&&has_positive;cS1
xE
newmax;if(min
lH3
iI1
newmax.set(fp_max(i41,i51
yG3
crosses_axis)min.set(eF2));e13
min
lH3
iI1
min.set(fp_min(i41,i51);iF1
min.known)min.set(i41);else
min.set(i51;}
max=newmax;}
eR1
range
xE::set_neg(l34
swap(min,max);min.val=-min.val;i53=-i53;}
yG1
IsLogicalTrueValue(const
range
xE&p,bool
abs){if(nA
IsIntType
xE::c53){if(p
i4&&p
xK1>=eF2
1
iY1
if(!abs&&p.max
lH3
p
nP<=eF2-1
iY1}
e13
p
i4&&p
xK1>=eF2
0.5
iY1
if(!abs&&p.max
lH3
p
nP<=eF2-0.5
iY1
eX3
tI3
yG1
IsLogicalFalseValue(const
range
xE&p,bool
abs){if(nA
IsIntType
xE::c53){if(abs)eY3.nT3
c71
1);else
eY3
i4&&p.max
lH3
p
xK1>eF2-1)c71
1);}
e13
abs)eY3.nT3
c71
0.5);else
eY3
i4&&p.max
lH3
p
xK1>eF2-0.5)c71
0.5);}
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;using
tE;lQ3{nS1
T>inline
int
isnan_workaround(T
t)yQ(t!=t);}
}
tE{tI1
range
xE
CalculateResultBoundaries(const
eU
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
{using
lQ3
FUNCTIONPARSERTYPES;range
xE
tmp=CalculateResultBoundaries_do(tree)xI1"Estimated boundaries: "
;if(tmp
i4)std::cout<<tmp
xK1;else
std::cout<<"-inf"
xI1" .. "
;if(tmp.iI1
std::cout<<tmp
nP;else
std::cout<<"+inf"
xI1": "
iQ3
tree)xI1
std::endl
n21
tmp;}
tI1
range
xE
lT1::CalculateResultBoundaries_do(const
eU
#endif
{iU
yW1(-fp_const_pihalf
xE(),fp_const_pihalf
xE());iU
pi_limits(-fp_const_pi
xE(),fp_const_pi
xE());iU
abs_pi_limits(xL1,fp_const_pi
xE());iU
plusminus1_limits(eF2-1)yI3
using
lQ3
std;switch
tT2{case
cImmed:nO
l14,l14);case
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cNot:case
cY3:case
cNotNot:case
cAbsNotNot:case
cEqual:case
i11:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:{nO
xL1
yI3}
case
cAbs:lD
m.set_abs(xH1
cLog:lD
nU3
fp_log
iO2
fp_log
xH1
cLog2:lD
nU3
fp_log2
iO2
fp_log2
xH1
cLog10:lD
nU3
fp_log10
iO2
fp_log10
xH1
cAcosh:lD
yC
tJ1
cGreaterOrEq
tG1
fp_acosh
lA4
cGreaterOrEq
tG1
fp_acosh
xH1
cAsinh
cJ1
fp_asinh
lO2
set(fp_asinh
xH1
cAtanh:lD
yC
n3-1),fp_atanh
lA4
cLess
tG1
fp_atanh
xH1
cAcos:lD
nO(m
xJ2
m
nP)<eF2
1))?fp_acos(m
nP):xL1,(yC
known&&cD3)>=eF2-1))?fp_acos
cD3):fp_const_pi
xE())eL3
cAsin:lD
yC
n3-1),fp_asin,yW1
xK1
lA4
cLess
tG1
fp_asin,yW1
nP
xH1
cAtan
cJ1
fp_atan,yW1
xK1
lO2
set(fp_atan,yW1
nP
xH1
cAtan2:{nG&&fp_equal(lR,xG1)yQ
abs_pi_limits;}
if(l41&&fp_equal(yB2,xG1)yQ
yW1;eX3
pi_limits
eL3
cSin:lD
bool
x11=!yC
known||!m.nT3||(m
nP-yC
val)>=(y2
x11)iV
Value_t
min=fp_mod
cD3,y2
min<xG1
min
yO
Value_t
max=fp_mod(m
nP,y2
max<xG1
max
yO
if(max<min)max
yO
bool
yI1=(min<=fp_const_pihalf
xE()&&max>=fp_const_pihalf
xE());bool
nO1=(min<=cF&&max>=cF
yG3
yI1&&nO1)iV
if(nO1)nO
eF2-1),nS2
if(yI1)nO
tB2
yI3
nO
tB2,nS2}
case
cCos:lD
if(eQ3
yC
val+=fp_const_pihalf
xE(yG3
m.iI1
m
nP+=fp_const_pihalf
xE();bool
x11=!yC
known||!m.nT3||(m
nP-yC
val)>=(y2
x11)iV
Value_t
min=fp_mod
cD3,y2
min<xG1
min
yO
Value_t
max=fp_mod(m
nP,y2
max<xG1
max
yO
if(max<min)max
yO
bool
yI1=(min<=fp_const_pihalf
xE()&&max>=fp_const_pihalf
xE());bool
nO1=(min<=cF&&max>=cF
yG3
yI1&&nO1)iV
if(nO1)nO
eF2-1),nS2
if(yI1)nO
tB2
yI3
nO
tB2,nS2}
case
cTan:{nO)eL3
cCeil:lD
m
eJ
cFloor
c32
xH1
cTrunc
c32);m
eJ
cInt
c32);m
eJ
cSinh
cJ1
fp_sinh
lO2
set(fp_sinh
xH1
cTanh
cJ1
fp_tanh,plusminus1_limits.min
lO2
set(fp_tanh,plusminus1_limits.max
xH1
cCosh:lD
if(eQ3{if(m.iI1{if
cD3>=xL1&&m
nP>=xG1{eR3
cX}
iF1
cD3)<xL1&&m
nP>=xG1{Value_t
tmp=cX
if(tmp>m
nP)m
nP=tmp;eR3
eF2
1);}
else{eR3
cX
std::swap
cD3,m
nP);}
}
else{if
cD3>=xG1{m.cH
eR3
fp_cosh
cD3);}
else{m.cH
eR3
eF2
1);}
}
}
else{yC
known=true;eR3
eF2
1
yG3
m.iI1{eR3
fp_cosh(m
nP);m.cH}
else
m.cH
eX3
m
eL3
cIf:case
cAbsIf:{range
xE
res1=eX
1)yJ3
res2=eX
2)yG3!res2
i4)res1
i4=false;iF1
res1
i4&&(res2
xK1)<res1
xK1)res1
xK1=res2
xK1;iF1
isnan_workaround(res2
xK1))res1
xK1=res2
xK1;if(!res2.iI1
res1.cH
iF1
res1
xJ2
res2
nP)>res1
nP)res1
nP=res2
nP;iF1
isnan_workaround(res2
nP))res1
nP=res2
nP
n21
res1
eL3
cMin:{bool
iW=false;bool
iX=false;range
yY3;x5
m
eN3!eQ3
iW
eO3
i4||cD3)<yZ3)yZ3=yC
val;if(!m.iI1
iX
eO3.nT3||(m
nP)<c03)c03=m
nP;}
if(iW)c13
iX)c23
cH
return
c33
cMax:{bool
iW=false;bool
iX=false;range
yY3;x5
m
eN3!eQ3
iW
eO3
i4||yC
val>yZ3)yZ3=yC
val;if(!m.iI1
iX
eO3.nT3||m
nP>c03)c03=m
nP;}
if(iW)c13
iX)c23
cH
return
c33
cAdd:{range
yY3(xL1,xG1;x5
item
eN3
item
i4)yZ3+=item
xK1;else
c13
item.iI1
c03+=item
nP;else
c23
cH
if(!c43&&!c23
iI1
y53
if
e12
i4&&c23
max
lH3
yZ3>c03)std::swap
e12
xK1,c03)n21
c33
cMul:{eI2
Value{enum
nV3{Finite,l22,PlusInf}
;nV3
tL
yX3
value;Value(nV3
t):tL(t),value(0){}
Value
tZ3
v):tL(Finite),value(v){}
bool
e22
cB1
tL==l22||(tL==Finite&&value<xG1
lT3
eQ1*=(c42
rhs){if(tL==Finite&&rhs.tL==Finite)value*=rhs.value;else
tL=(e22)!=rhs.e22)?l22:PlusInf);}
iY2<(c42
rhs
cB1(tL==l22&&rhs.tL!=l22)||(tL==Finite&&(rhs.tL==PlusInf||(rhs.tL==Finite&&value<rhs.value)));}
}
;eI2
yX1{Value
c52,c62;yX1():c52(Value::PlusInf),c62(Value::l22){}
void
xP2
Value
e23,c42
value2){e23*=value2;if(e23<c52)c52=e23;if(c62<e23)c62=e23;}
}
;range
yY3(yH3
x5
item
eN3!item
i4&&!item.iI1
nO);Value
nW3=c43?Value
e12
xK1
lJ2
l22);Value
nX3=c23
nT3?Value
e12
nP
lJ2
PlusInf);Value
nY3=item
i4?Value(item
xK1
lJ2
l22);Value
nZ3=item.nT3?Value(item
nP
lJ2
PlusInf);yX1
range
lV3
nW3,nY3)lV3
nW3,nZ3)lV3
nX3,nY3)lV3
nX3,nZ3
yG3
range.c52.tL==Value::Finite)yZ3=range.c52.value;else
c13
range.c62.tL==Value::Finite)c03=range.c62.value;else
c23
cH
if(!c43&&!c23
iI1
y53
if
e12
i4&&c23
max
lH3
yZ3>c03)std::swap
e12
xK1,c03)n21
c33
cMod:{range
xE
x=CalculateResultBoundaries(eV
yJ3
y=eX
1)yG3
y.iI1{if(y
nP>=xG1{if(!x
i4||(x
xK1)<xG1
nO-y
nP,y
nP);else
nO
xL1,y
nP);}
e13!x.nT3||(x
nP)>=xG1
nO
y
nP,-y
nP);else
nO
y
nP,fp_const_negativezero
xE());}
}
else
nO)eL3
cPow:{if(l41&&yB2==xG1{nO
yH3}
nG&&lR==xG1{nO
xL1,xG1;}
nG&&fp_equal(lR
yF3{nO
yH3}
if(l41&&yB2>xL1&&GetEvennessInfo(tV2)==tA2)cL1
yB2;range
xE
tmp=CalculateResultBoundaries(eV
yJ3
c53;c43=true;yZ3=0;if(tmp
i4&&tmp
xK1>=xG1
yZ3
eP3
tmp
xK1
e32
iF1
tmp.max
lH3
tmp
nP<=xG1
yZ3
eP3
tmp
nP
e32
c23
cH
if(tmp
i4&&tmp.iI1{c23
nT3=true;c03=fp_max(fp_abs(tmp
xK1),fp_abs(tmp
nP));c03=fp_pow
e12
nP
e32
eX3
c53;}
range
xE
tW
eV
yJ3
p1=eX
1));TriTruthValue
p0_positivity=(p0
i4&&(p0
xK1)>=xG1?tA2:(p0
xJ2
p0
nP)<xL1?lI3
Unknown);TriTruthValue
e82=GetEvennessInfo(tV2);TriTruthValue
tA=Unknown;switch(p0_positivity)i61
tA=tA2;lC
lI3{tA=e82;y53
yT3
switch(e82)i61
tA=tA2;lC
lI3
lC
Unknown:{if(l41&&!tR2
yB2)&&yB2>=xG1{tA=tA2;}
y53}
t93
tA)i61{Value_t
min=xL1;if(p0
i4&&p1
i4){min
eP3
p0
xK1,p1
xK1
yG3
p0
xK1<xL1&&(!p1.nT3||p1
nP>=xG1&&min>=xG1
min=xL1;}
if(p0
i4&&p0
xK1>=xL1&&p0.max
lH3
p1.iI1{Value_t
max
eP3
p0
nP,p1
nP
yG3
min>max)std::swap(min,max);nO
min,max);}
nO
min,false)eL3
lI3{nO
false,fp_const_negativezero
xE());}
yT3{y53
yC1
cNeg:lD
m.set_neg(xH1
cSub:{yD
cNeg
tA3
1
tE3
iD
cAdd);tmp.nI
0
tE3
yD1
tmp2)tC2
cInv:{lU1
lW-1)))tC2
cDiv:{yD
cInv
tA3
1
tE3
iD
xA1
nU2
tmp2)tC2
cRad:yY
xA1
yE
fp_const_rad_to_deg
xE()))tC2
cDeg:yY
xA1
yE
fp_const_deg_to_rad
xE()))tC2
cSqr:{lU1
lW
2)))tC2
cExp:yY
cPow);tmp.yE
fp_const_e
xE()tE3.nI
0))tC2
cExp2:yY
cPow);tmp.yE
x03
tmp.nI
0))tC2
cCbrt
cJ1
fp_cbrt
lO2
set(fp_cbrt
xH1
cSqrt:lD
if(eQ3
eR3
cD3)<xL1?0:fp_sqrt
cD3
yG3
m.iI1
m
nP=(m
nP)<xL1?0:fp_sqrt(m
nP
xH1
cRSqrt:{lU1
lW-0.5)))tC2
cHypot:{lT1
xsqr,ysqr,add,sqrt;xsqr.nI
0));xsqr.yE
x03
ysqr.nI
1));ysqr.yE
x03
xsqr
iD
cPow);ysqr
iD
cPow);add
yD1
xsqr);add
yD1
ysqr);add
iD
cAdd);sqrt
yD1
add);sqrt
iD
cSqrt)n21
CalculateResultBoundaries(sqrt)eL3
e63:{yD
cLog2
tA3
0
tE3
t83
tmp
yD1
tmp2);tmp.nI
1))tC2
cCot:{yD
cTan);tmp2
yU
lH
cSec:{yD
cCos);tmp2
yU
lH
cCsc:{yD
cSin);tmp2
yU
CalculateResultBoundaries(tmp);}
lC
cRDiv:case
cRSub:case
cDup:case
cFetch:case
cPopNMov:case
cSinCos:case
cSinhCosh:case
cNop:case
cJump:case
l93:lC
cArg:case
cConj:case
cImag:case
cReal:case
cPolar:lC
cPCall:lC
cFCall:y53
nO);}
tI1
TriTruthValue
cF3
const
eU{switch
tT2{case
cImmed:return
tR2
l14)?tA2:IsNever;case
cFloor:case
cCeil:case
cTrunc:case
cInt:return
tA2;case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cEqual:case
i11:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:return
tA2;case
cIf:{TriTruthValue
a=cF3
tV2);TriTruthValue
b=GetIntegerInfo
xX2
if(a==b)return
a
x02
case
cAdd:case
cMul:{for
yW
if(cF3
t01)!=tA2)return
Unknown
n21
tA2;}
yT3
y53
return
Unknown;}
yG1
IsLogicalValue(const
eU{switch
tT2{case
cImmed:return
fp_equal(l14,xG1||fp_equal(l14
yI3
case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cAbsAnd:case
cAbsOr:case
cY3:case
cAbsNotNot:case
cEqual:case
i11:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:x0
cMul:{for
yW
if(!IsLogicalValue(t01)cG
return
true
eL3
cIf:case
cAbsIf:yQ
IsLogicalValue(tV2)&&IsLogicalValue
xX2}
yT3
y53
return
tI3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;
#if defined(__x86_64) || !defined(FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS)
# define CBRT_IS_SLOW
#endif
#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif
lQ3
xS1{extern
const
unsigned
char
powi_table[256];}
lQ3{using
tE
yF1
bool
IsOptimizableUsingPowi(long
immed,long
penalty=0){xS1::nV2
xE
synth;synth.PushVar(l93);size_t
bytecodesize_backup
tY3
GetByteCodeSize();xS1::x01
immed,xS1::tP1
xE::MulSequence,tN2
size_t
bytecode_grow_amount
tY3
GetByteCodeSize()-bytecodesize_backup
n21
bytecode_grow_amount<size_t(MAX_POWI_BYTECODE_LENGTH-penalty);}
eR1
ChangeIntoRootChain(lT1&tree,bool
lJ3,long
iP2,long
iQ2){while(iQ2>0)yY
cCbrt);xB1);tmp.Rehash(xL2--iQ2;}
while(iP2>0)yY
cSqrt
yG3
lJ3){tmp
iD
cRSqrt);lJ3=tI3
xB1);tmp.Rehash(xL2--iP2;}
if(lJ3)yY
cInv);xB1
xL2}
}
tI1
eI2
RootPowerTable{static
const
Value_t
RootPowers[(1+4)*(1+3)];}
yF1
const
Value_t
tG(1+4)*(1+3)]={eF2
1)lS
i63
i63
2*i63
2*2*2)lS
3)lS
3*2)lS
3*lR2
2*lR2
2*2*lR2
3
lK3
2
lK3
lR2
3*2*lR2
3*2*2*lR2
3*3
lK3
3*2
lK3
3*lR2
3*3*2*lR2
3*3*2*2*2*2)}
;eI2
PowiResolver{static
const
unsigned
MaxSep=4;static
x43
i73=5
i03
int
e33
i03
long
xA3
i03
long
tI;eI2
c72{c72():n_int_sqrt(0),n_int_cbrt(0),sep_list(),n31(0){}
int
n_int_sqrt;int
n_int_cbrt;int
sep_list[MaxSep];tI
n31;}
yF1
static
c72
CreatePowiResult
tZ3
eE2){c72
c53;e33
tH=FindIntegerFactor(eE2
yG3
tH==0){
#ifdef DEBUG_POWI
iR2"no factor found for %Lg\n"
,(e43);
#endif
return
c53;}
e53=yJ1
eE2,tH);xA3
tD2=EvaluateFactorCost(tH,0,0,0)+cD
e53);int
i83=0;int
i93=0;int
x33=0;
#ifdef DEBUG_POWI
iR2"orig = %Lg\n"
,(e43);iR2"plain factor = "
iT3"%ld\n"
,(int)tH,(long)tD2);
#endif
for
nY2
n_s=0;n_s<MaxSep;++n_s){int
xB=0;xA3
yY1=tD2;e33
c91=tH;for(int
s=1;s<i73*4;++s){
#ifdef CBRT_IS_SLOW
if(s>=i73)break;
#endif
int
n_sqrt=s%i73;int
n_cbrt=s/i73;if(n_sqrt+n_cbrt>4
cF2
yX3
lN1=eE2;lN1-=tG
s];i81=FindIntegerFactor(lN1
yG3
y02!=0){tI
xR=yJ1
lN1,y02);xA3
cost=EvaluateFactorCost(y02,i83+n_sqrt,i93+n_cbrt,x33+1)+cD
xR);
#ifdef DEBUG_POWI
iR2"Candidate sep %u (%d*sqrt %d*cbrt)factor = "
iT3"%ld (for %Lg to %ld)\n"
,s,n_sqrt,n_cbrt,y02,(long)cost
eO2
lN1,(long)xR);
#endif
if(cost<yY1){xB=s;c91=y02;yY1=cost;}
}
}
if(!xB)break;
#ifdef DEBUG_POWI
iR2"CHOSEN sep %u (%d*sqrt %d*cbrt)factor = "
iT3"%ld, exponent %Lg->%Lg\n"
,xB,xB%i73,xB/i73,c91,yY1
eO2(eE2)eO2(eE2-tG
xB]));
#endif
c53
cM1
n_s]=xB;eE2-=tG
xB];i83+=xB%i73;i93+=xB/i73;tD2=yY1;tH=c91;x33+=1;}
e53=yJ1
eE2,tH);
#ifdef DEBUG_POWI
iR2"resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n"
,e53,(e43
eO2
tH);
#endif
while(tH%2==0){++c53
eC2;tH/=2;}
while(tH%3==0){++c23
n_int_cbrt;tH/=3;eX3
c53;}
private:static
xA3
cD
tI
xR){static
std::map
e92
iM;if(xR<0){xA3
cost=22
n21
cost+cD-xR);}
std::map
e92::yE3
i=iM.y12
xR
yG3
i!=iM.e71
xR)return
i
eN2;std::pair
e92
c53(xR,0.0);xA3&cost=c23
iD3
while(xR>1){int
y02=0;if(xR<256){y02=xS1::powi_table[xR];if(y02&128)y02&=127;else
y02=0;if(y02&64)y02=-(y02&63)-1;}
if(y02){cost+=cD
y02);xR/=y02;tK3
if(!(xR&1)){xR/=2;cost+=6;}
else{cost+=7;xR-=1;}
}
iM.y93,c53)n21
cost;}
cK1
tI
yJ1
yZ1,i81)yQ
makeLongInteger(value*eF2
y02));}
cK1
bool
c11
yZ1,i81
yF2
v=value*eF2
y02)n21
isLongInteger(v);}
cK1
e33
FindIntegerFactor(yZ1){i81=(2*2*2*2);
#ifdef CBRT_IS_SLOW
#else
y02*=(3*3*3);
#endif
e33
c53=0;if(c11
value,y02)){c53=y02;while((y02%2)==0&&c11
value,y02/2))c53=y02/=2;while((y02%3)==0&&c11
value,y02/3))c53=y02/=3;}
#ifdef CBRT_IS_SLOW
if
e12==0){if(c11
value,3
e42
3;}
#endif
return
c53;}
static
int
EvaluateFactorCost(int
y02,int
s,int
c,int
nmuls){x43
x53=6;
#ifdef CBRT_IS_SLOW
x43
tE2=25;
#else
x43
tE2=8;
#endif
int
c53=s*x53+c*tE2;while(y02%2==0){y02/=2;c53+=x53;}
while(y02%3==0){y02/=3;c53+=tE2;}
c53+=nmuls
n21
c53;}
}
;}
tE{yG1
lT1::RecreateInversionsAndNegations(bool
prefer_base2){bool
changed=false
lL2
0;a<c61;++a)if(lZ1.RecreateInversionsAndNegations(prefer_base2))c31
if(changed){exit_changed:Mark_Incompletely_Hashed(lQ2
switch(iU1{case
cMul:{eO
nC2
i2
nD2,e21;if(true){bool
nP1=false
yX3
xM2=0;x63
c82
0)eA2
tQ
1)y01){nP1=true;xM2=tQ
iG;y53}
if(nP1
yF2
immeds=1.0;x63
y01){immeds*=powgroup
y11;c21}
for
iY-->0;){lT1&powgroup=lZ1;if(powgroup
c82
0)eA2
tQ
1).IsImmed(nM&log2=tQ
0);log2.lF1
log2
iD
e63);log2.yE
fp_pow(immeds
tO1
1)/xM2)));log2
lP2;y53}
}
}
x63
c82
1)y01){cI2&exp_param=tQ
1)yX3
eE2=exp_param
y11;if(e91
tO1-1))){lF1
nC2.push_back(lZ1
cV3
c21
iF1
eE2<xL1&&tR2
eE2
nM
iZ;iZ
iD
cPow);iZ
e9
tQ
0));iZ.yE-eE2));iZ
lP2;nC2.push_back(iZ);lF1
c21}
iF1
powgroup
eA2!nD2.c92)){nD2=tQ
0);lF1
c21
iF1
powgroup
nC==e63&&!e21.c92)){e21=powgroup;lF1
c21}
if(!nC2.empty()){changed=true
i2
eJ1;eJ1
t83
eJ1
e83
nC2);eJ1
lP2
i2
iC
cMul);yE1
SetParamsMove
eK
if
eM3
y01&&fp_equal
eM3
y11
yF3
eG2
cInv)eL
eJ1);}
e13
yE1
GetDepth()>=eJ1
nM2)eG2
cDiv
tK1
eL
eJ1);}
else
eG2
cRDiv)eL
eJ1
tK1;}
}
}
if(nD2.c92
nM
iC
iU1;yE1
SetParamsMove
eK
while(yE1
RecreateInversionsAndNegations(prefer_base2))yE1
FixIncompleteHashes();SetOpcode(e63)eL
nD2
tK1;c31}
if(e21.c92
nM
iC
cMul);mulgroup
yD1
e21
l8
1));yE1
AddParamsMove
eK
while(yE1
RecreateInversionsAndNegations(prefer_base2))yE1
FixIncompleteHashes();DelParams();SetOpcode(e63)eL
e21
l8
0)tK1;c31
yC1
cAdd:{eO
iS2;for
iY-->0;)if(e73
cMul){nE2
yK1:i2&mulgroup
e93
for(n03
mulgroup
y41;b-->0;){if
eM3
l8
b)yE2
y02=mulgroup
l8
b)y11;if(fp_equal(y02
eB2
yK1;}
yE1
lF1
yE1
DelParam(b);lL3
iF1
fp_equal(y02
tO1-2))){xO
yK1;}
yE1
lF1
yE1
DelParam(b);yE1
yE
x03
lL3}
}
if(tB){yE1
tJ
mulgroup);c21}
iF1
e73
cDiv&&!IsIntType
xE::c53){nE2
yL1:i2&eJ1
e93
if(eJ1
l8
0)y01){if(fp_equal(eJ1
l8
0)y11
eB2
yL1;}
eJ1.lF1
eJ1.DelParam(0);eJ1
iD
cInv);lL3}
if(tB){xO
yL1;}
eJ1.tJ
eJ1);c21}
iF1
e73
cRDiv&&!IsIntType
xE::c53){nE2
xC1:i2&eJ1
e93
if(eJ1
l8
1)y01){if(fp_equal(eJ1
l8
iG
eB2
xC1;}
eJ1.lF1
eJ1.DelParam(1);eJ1
iD
cInv);lL3}
if(tB){xO
xC1;}
eJ1.tJ
eJ1);c21}
if(!iS2.empty()){
#ifdef DEBUG_SUBSTITUTIONS
iR2"Will make a Sub conversion in:\n"
);fflush(stdout);l21
#endif
lT1
cA2;cA2
iD
cAdd);cA2
e83
iS2);cA2
lP2
i2
e31;e31
iD
cAdd);e31
e83
l52));e31
lP2;if(e31
y01&&fp_equal(e31
y11,xG1)eG2
cNeg);eC);}
e13
e31
nM2==1)eG2
cRSub);eC)eA3}
iF1
cA2
nC==cAdd)eG2
cSub)eA3
eC
l8
0))lL2
1;a<cA2
y41;++a){lT1
tF2;tF2
iD
cSub);tF2
e83
l52));tF2.yT2)eL
tF2);eC
l8
a));}
}
else
eG2
cSub)eA3
eC);}
}
#ifdef DEBUG_SUBSTITUTIONS
iR2"After Sub conversion:\n"
);fflush(stdout);l21
#endif
yC1
cPow:{cI2&p0
nF2
0);cI2&p1
nF2
1
yG3
p1
y01){if
n91!=xL1&&!isInteger
n91)){eG
c72
r=eG
CreatePowiResult(fp_abs
n91)yG3
r.n31!=0){bool
l32=false;if
n91<xL1&&r
cM1
0]==0&&r
eC2>0){l32=true;}
#ifdef DEBUG_POWI
iR2"Will resolve powi %Lg as powi(chain(%d,%d),%ld)"
eO2
fp_abs
n91),r
eC2,r.n_int_cbrt,r.n31);for
nY2
n=0;n<eG
MaxSep;++n){if(r
cM1
n]==0)break;int
n_sqrt=r
cM1
n]%eG
i73;int
n_cbrt=r
cM1
n]/eG
i73;iR2"*chain(%d,%d)"
,n_sqrt,n_cbrt);}
iR2"\n"
);
#endif
lT1
eD2
nF2
0)i2
cB2=eD2;cB2.lF1
ChangeIntoRootChain(cB2,l32,r
eC2,r.n_int_cbrt);cB2
lP2
i2
pow;if(r.n31!=1){pow
iD
cPow);pow
yD1
cB2);pow.yE
eF2
r.n31)));}
else
pow.swap(cB2)i2
mul;mul
t83
mul
yD1
pow);for
nY2
n=0;n<eG
MaxSep;++n){if(r
cM1
n]==0)break;int
n_sqrt=r
cM1
n]%eG
i73;int
n_cbrt=r
cM1
n]/eG
i73
i2
tG2=eD2;tG2.lF1
ChangeIntoRootChain(tG2,false,n_sqrt,n_cbrt);tG2
lP2;mul
yD1
tG2);}
if
n91<xL1&&!l32){mul
lP2;SetOpcode(cInv);nG1
0,mul);DelParam(1);}
else
eG2
cMul);SetParamsMove(mul.l52));}
#ifdef DEBUG_POWI
l21
#endif
c31
y53}
}
if(GetOpcode()==cPow&&(!p1
y01||!isLongInteger
n91)||!IsOptimizableUsingPowi
xE(makeLongInteger
n91)))){if(p0
y01&&p0
y11>eF2
0.0)){if(prefer_base2
yF2
cC2=fp_log2(p0
y11
yG3
fp_equal(cC2
yF3{DelParam(0);}
else{n0
eB1
cC2));eE2
e9
p1)c01
nG1
nA1}
SetOpcode(cExp2);c31}
else{Value_t
cC2=fp_log(p0
y11
yG3
fp_equal(cC2
yF3{DelParam(0);}
else{n0
eB1
cC2));eE2
e9
p1)c01
nG1
nA1}
SetOpcode(cExp);c31}
}
iF1
GetPositivityInfo(p0)==tA2){if(prefer_base2){lT1
log;log
iD
cLog2);log
e9
p0);log
lP2;n0
p1);eE2
yD1
log)c01
SetOpcode(cExp2);nG1
nA1
c31}
else{lT1
log;log
iD
cLog);log
e9
p0);log
lP2;n0
p1);eE2
yD1
log)c01
SetOpcode(cExp);nG1
nA1
c31}
}
yC1
cDiv:{if(GetParam(0)y01&&fp_equal(GetParam(0)y11
yF3
eG2
cInv);DelParam(0);}
y53
yT3
y53
if(changed)goto
exit_changed
n21
changed;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lQ3
FUNCTIONPARSERTYPES;lQ3{using
tE;class
x23
iA3
nQ1;size_t
eH;size_t
eI;size_t
lO1;size_t
tC;size_t
tD;size_t
nB1;eS3
x23():nQ1(0),eH(0),eI(0),lO1(0),tC(0),tD(0),nB1(0){}
void
iB3
OPCODE
op){nQ1+=1;x13
cCos)++eH;x13
cSin)++eI;x13
cSec)++eH;x13
cCsc)++eI;x13
cTan)++lO1;x13
cCot)++lO1;x13
cSinh)++tD;x13
cCosh)++tC;x13
cTanh)++nB1;}
size_t
GetCSEscore()const
iA3
c53=nQ1
n21
c53;}
int
NeedsSinCos()const{bool
yM1=(nQ1==(eH+eI+lO1)yG3(lO1&&(eI||eH))||(eI&&eH)){if(yM1)return
1
n21
2;eX3
0;}
int
NeedsSinhCosh()const{bool
yM1=(nQ1==(tC+tD+nB1)yG3(nB1&&(tD||tC))||(tD&&tC)){if(yM1)return
1
n21
2;eX3
0;}
size_t
MinimumDepth()const
iA3
n_sincos=std::min(eH,eI);size_t
n_sinhcosh=std::min(tC,tD
yG3
n_sincos==0&&n_sinhcosh==0)return
2
n21
1;}
}
yF1
class
TreeCountType:public
std::multimap<fphash_t,std::pair<x23,lT1> >{}
xM3
FindTreeCounts(tH1&nG2,cI2&tree,OPCODE
xN2,bool
skip_root=false){cY
i=nG2.y12
tree.GetHash()yG3!skip_root){bool
found=false;for(;i!=nG2.e71
tree.GetHash();++i){if(tree
xK
i
eN2.second)){i
eN2.first.iB3
xN2);found=true;y53}
if(!found){x23
count;count.iB3
xN2);nG2.y93,std::make_pair(tree.GetHash(),std::make_pair
eF3,tree)));}
}
for
xH
FindTreeCounts(nG2,t01,i42);}
eI2
c0{bool
BalanceGood;bool
FoundChild;}
yF1
c0
lP1
cI2&root,cI2&child){if(root
xK
child)){c0
c53={true,true}
n21
c53;}
c0
c53={true,false}
;if(root
nC==cIf||root
nC==tG3{c0
cond=lP1
root
l8
0
lN3
c0
xY=lP1
root
l8
1
lN3
c0
y5=lP1
root
l8
2
lN3
if
lO3||xY
c1||y5
c1){c53
c1=true;}
c53
eD=((xY
c1==y5
c1)||lO3
iT2)&&(cond
eD||(xY
c1&&y5
c1))&&(xY
eD||lO3
iT2)&&(y5
eD||lO3
iT2);}
else{bool
i91=false;bool
nR1=false;for(n03
root
y41,a=0;a<b;++a){c0
tmp=lP1
root
l8
a
lN3
if(tmp
c1)c53
c1=true;if(tmp
eD==false)i91=true;iF1
tmp
c1)nR1=true;}
if(i91&&!nR1)c53
eD=false;eX3
c53;}
yG1
yD3
cI2&iC3
cI2&tree,const
xS1::nV2
xE&synth,const
tH1&nG2){for(n03
eR,a=0;a<b;++a){cI2&leaf=t01;cY
synth_it;y42
tH1::const_iterator
i=nG2.yC3
i!=nG2
iW3;++i){if(i->first!=leaf.GetHash()cF2;const
x23&occ
nT2
first;size_t
score=occ.GetCSEscore();cI2&candidate
nT2
iD3
nH2
candidate)cF2;if(leaf
nM2<occ.MinimumDepth()cF2;if(score<2
cF2;if(lP1
iC3
leaf)eD==false
cF2
tI2
if(yD3
iC3
leaf,synth,nG2
iY1
eX3
tI3
yG1
iA1
cI2&yA3,cI2&expr){for(xM
yA3
l8
a)xK
expr
iY1
for(xM
iA1
yA3
l8
a),expr
e42
true
n21
tI3
yG1
GoodMomentForCSE(cI2&yA3,cI2&expr){if(yA3
nC==cIf)return
true;for(xM
yA3
l8
a)xK
expr
iY1
size_t
iU2=0;for(xM
iA1
yA3
l8
a),expr))++iU2
n21
iU2!=1;}
}
tE{tI1
size_t
lT1::SynthCommonSubExpressions(xS1::nV2
xE&synth)const{if(c61==0)return
0;size_t
stacktop_before
tY3
GetStackTop();tH1
nG2;FindTreeCounts(nG2,*this,GetOpcode(),true);for(;;)iA3
cD2=0;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Finding a CSE candidate, root is:"
<<std::endl;DumpHashes(*this);
#endif
cY
cs_it(nG2
iW3);for(cY
j=nG2.yC3
j!=nG2
iW3;){cY
i(j++);const
x23&occ
nT2
first;size_t
score=occ.GetCSEscore();cI2&tree
nT2
iD3
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Score "
<<score<<":\n"
<<std::flush;DumpTreeWithIndent(tree);
#endif
nH2
tree))xV
if(tree
nM2<occ.MinimumDepth())xV
if(score<2)xV
if(lP1*this,tree)eD==false)xV
if(yD3*this,tree,synth,nG2)){tK3
if(!GoodMomentForCSE(*this,tree))xV
score*=tree
nM2;if(score>cD2){cD2=score;cs_it=i;}
}
if(cD2<=0){
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"No more CSE candidates.\n"
<<std::flush;
#endif
y53
cI2&tree=cs_it
eN2.iD3
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<iZ3"Common Subexpression:"
;DumpTree
xE(tree)xI1
std::endl;
#endif
#if 0
int
n41=occ.NeedsSinCos();int
i7=occ.NeedsSinhCosh()i2
iV2,iW2,cE2,cG2;if(n41){iV2
tH2
iV2
iD
cSin);iV2
lP2;iW2
tH2
iW2
iD
cCos);iW2
lP2;nH2
iV2)||synth.Find(iW2)){if(n41==2){iZ1
tK3
n41=0;}
}
if(i7){cE2
tH2
cE2
iD
cSinh);cE2
lP2;cG2
tH2
cG2
iD
cCosh);cG2
lP2;nH2
cE2)||synth.Find(cG2)){if(i7==2){iZ1
tK3
i7=0;}
}
#endif
tree.SynthesizeByteCode(synth,false);iZ1
#ifdef DEBUG_SUBSTITUTIONS_CSE
synth.xC3
Dump<0>()xI1"Done with Common Subexpression:"
;DumpTree
xE(tree)xI1
std::endl;
#endif
#if 0
if(n41){if(n41==2||i7){synth.eP1}
y43
cSinCos,1,2)y21
iV2,1)y21
iW2,0);}
if(i7){if(n41)synth.eP1
if(i7==2){synth.eP1}
y43
cSinhCosh,1,2)y21
cE2,1)y21
cG2,0);}
#endif
eX3
y33
stacktop_before;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
tI1
lS1
xE::iX2{using
tE;c02)i2
tree;tree.GenerateFrom(*mData);FPoptimizer_Optimize::ApplyGrammars(tree);std
xT3<unsigned>eB3;std
xT3
xE
immed;size_t
stacktop_max=0;tree.SynthesizeByteCode(eB3,immed,stacktop_max
yG3
mData->mStackSize!=stacktop_max){mData->mStackSize=unsigned(stacktop_max);
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
mData->mStack
xK3
stacktop_max);
#endif
}
mData->mByteCode.swap(eB3);mData->mImmed.swap(immed);}
#define FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(type) tM1>lS1<type>::iX2{}
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
iG3(MpfrFloat)
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
iG3(GmpInt)
#endif
#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
iG3(std::complex<double>)
#endif
#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
iG3(std::complex<float>)
#endif
#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
iG3(std::complex<long
double>)
#endif
#define FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(type) xC3 lS1<type>::iX2;
#ifndef FP_DISABLE_DOUBLE_TYPE
iH3(double)
#endif
#ifdef FP_SUPPORT_FLOAT_TYPE
iH3(float)
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
iH3(long
double)
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
iH3(long)
#endif
#endif // FP_SUPPORT_OPTIMIZER

#endif
