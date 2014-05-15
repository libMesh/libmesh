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
#define lA4 ||tree.GetOpcode
#define l94 l74 data
#define l84 param iW
#define l74 param.
#define l64 (l74
#define l54 :if lI
#define l44 b.Value)
#define l34 b.Opcode
#define l24 .end()
#define l14 1,252,
#define l04 {if(n31
#define iZ3 sim.xA1
#define iY3 "Found "
#define iX3 lT 1,0,
#define iW3 stackpos
#define iV3 ,lT 2,
#define iU3 "dup(%u) "
#define iT3 ,cCos lY1
#define iS3 eK{assert
#define iR3 "%d, cost "
#define iQ3 "PUSH " yH2
#define iP3 "immed "<<
#define iO3 mFuncParsers
#define iN3 stderr
#define iM3 sep2=" "
#define iL3 FPHASH_CONST
#define iK3 cache_needed[
#define iJ3 fprintf
#define iI3 ::cout<<"Applying "
#define iH3 FUNCTIONPARSER_INSTANTIATE_OPTIMIZE
#define iG3 FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE
#define iF3 HANDLE_UNARY_CONST_FUNC
#define iE3 within,
#define iD3 c_count
#define iC3 s_count
#define iB3 MaxOp
#define iA3 2)lS 2*
#define i93 NumConstant:
#define i83 switch(xD2.first xI3
#define i73 FP_GetOpcodeName
#define i63 e22&xD2
#define i53 xD3))lZ2
#define i43 {cV1 xD3(
#define i33 .what nW1
#define i23 ;unsigned
#define i13 ].swap(
#define i03 codes[b
#define tZ3 whydump
#define tY3 info.
#define tX3 tY3 SaveMatchedParamIndex(
#define tW3 :if(yT3
#define tV3 for(;a<
#define tU3 nparams
#define tT3 281856,
#define tS3 (count
#define tR3 l4 0,1,
#define tQ3 0x12 nH
#define tP3 iO 1,0,
#define tO3 nQ 0,
#define tN3 cAbs nQ
#define tM3 xC cMul);
#define tL3 xD lT 2,
#define tK3 },{{1,0
#define tJ3 ;case
#define tI3 iF tJ3
#define tH3 lP2 xC2
#define tG3 lP2 cK1
#define tF3 lP2 cJ1
#define tE3 nD1,l5::
#define tD3 nD1 lP2
#define tC3 fp_pow(
#define tB3 cAbsIf)
#define tA3 l21++b)
#define t93 t73 l31
#define t83 );range.xI2
#define t73 t5 a)
#define t63 .second
#define t53 ]t63
#define t43 ].first
#define t33 yM r);}
#define t23 Ne_Mask
#define t13 Gt_Mask
#define t03 Lt_Mask
#define eZ3 FindPos
#define eY3 i6|l51)
#define eX3 exponent
#define eW3 (eX3
#define eV3 resize(
#define eU3 .nJ 0))
#define eT3 ,yD3 l8
#define eS3 .nJ 1));
#define eR3 (*x5 cG3
#define eQ3 public:
#define eP3 {data->
#define eO3 IsNever
#define eN3 pclone
#define eM3 cOr,l6
#define eL3 xK(tG2
#define eK3 newpow
#define eJ3 change
#define eI3 133,2,
#define eH3 Params
#define eG3 Needs
#define eF3 byteCode
#define eE3 AddFrom(
#define eD3 eN cU1);
#define eC3 tF1){if(
#define eB3 cLog2by
#define eA3 long iA1
#define e93 factor_t
#define e83 value1
#define e73 synth.
#define e63 iC2 yQ
#define e53 cV1(1))
#define e43 else{if(
#define e33 xK());nD
#define e23 cO val)<
#define e13 ;x0 l5
#define e03 nN1);}if
#define cZ3 xF1 tL{
#define cY3 goto cV
#define cX3 p2 tM ifp2
#define cW3 iG p2;p2
#define cV3 cAbsNot
#define cU3 (m cO val
#define cT3 eI2))i3
#define cS3 eI2 nC==
#define cR3 stackptr
#define cQ3 cLog);xL
#define cP3 cV1(0.5)
#define cO3 l8 0));
#define cN3 opcodes
#define cM3 ::vector<unsigned>&eR1
#define cL3 ::nL2 xK
#define cK3 did_muli
#define cJ3 &Value){
#define cI3 yC const
#define cH3 used[b]
#define cG3 lS1,info
#define cF3 :{n51 r=
#define cE3 std::yD2
#define cD3 sizeof(
#define cC3 cAbsIf,
#define cB3 cNotNot,
#define cA3 l4 16,1,
#define c93 cLess,cD
#define c83 l4 20,1,
#define c73 l4 4,1,
#define c63 cTan lY1
#define c53 ,tree,info
#define c43 n53 2,2,
#define c33 ,cPow lY1
#define c23 n53 0,2,
#define c13 450998,
#define c03 cLog lY1
#define yZ3 cExp2 nQ
#define yY3 lJ 2},0,
#define yX3 if(param
#define yW3 ::string
#define yV3 t5 1)nC==
#define yU3 );if(
#define yT3 &*(*x5 lS1){
#define yS3 (*x5)[0].info
#define yR3 param=*
#define yQ3 Ge0Lt1
#define yP3 ;for l51
#define yO3 .n7 synth
#define yN3 break;}
#define yM3 default:
#define yL3 Gt0Le1
#define yK3 &&p cO val
#define yJ3 if(p0 cO
#define yI3 range<nU
#define yH3 range nQ3
#define yG3 range xK
#define yF3 cAdd lJ2
#define yE3 (op1==
#define yD3 leaf1
#define yC3 .tW1 n]
#define yB3 )))iF lH
#define yA3 t63);
#define y93 1),e53;
#define y83 ;}static yO1
#define y73 (cond yT cU2
#define y63 e53)
#define y53 iterator
#define y43 begin();
#define y33 TreeSet
#define y23 parent
#define y13 insert(i
#define y03 newrel
#define xZ3 b_needed
#define xY3 cachepos
#define xX3 grammar
#define xW3 half=
#define xV3 131,4,1,
#define xU3 ,i7,1,iZ+1);
#define xT3 131,8,1,
#define xS3 4,1,2,1,
#define xR3 src_pos
#define xQ3 reserve(
#define xP3 n11 void
#define xO3 treeptr
#define xN3 tM1 void
#define xM3 ImmedTag
#define xL3 IsLogicalValue(t5
#define xK3 a,const
#define xJ3 RefCount
#define xI3 ){case
#define xH3 Birth();
#define xG3 typename
#define xF3 mulgroup
#define xE3 template
#define xD3 result
#define xC3 cost_t
#define xB3 .known
#define xA3 fpdata
#define x93 middle
#define x83 7168,
#define x73 };enum
#define x63 sqrt_cost
#define x53 const int
#define x43 mul_count
#define x33 cV1(2)));
#define x23 maxValue1
#define x13 minValue1
#define x03 maxValue0
#define nZ3 minValue0
#define nY3 ValueType
#define nX3 ;TriTruthValue
#define nW3 xD3 iN
#define nV3 xD3 lD3
#define nU3 xD3 e8
#define nT3 xD3 cO known
#define nS3 nW3=false;if(
#define nR3 xD3 cO val
#define nQ3 xK xD3
#define nP3 lR,eL)yS2
#define nO3 abs_mul
#define nN3 l43 size()
#define nM3 t5 0)nC==
#define nL3 l8 a));
#define nK3 pos_set
#define nJ3 p1 yM p1)
#define nI3 default_function_handling
#define nH3 y6 2,cAdd
#define nG3 synth lE2
#define nF3 ,eJ,synth);
#define nE3 ;if(half
#define nD3 ;}void
#define nC3 xQ)nD3
#define nB3 cMul);xQ.
#define nA3 subtree
#define n93 invtree
#define n83 MakeHash(
#define n73 rulenumit
#define n63 cAnd l3
#define n53 ,cPow,l2
#define n43 cAnd,l6
#define n33 cEqual,
#define n23 n33 lA
#define n13 t31},{{3,
#define n03 MakeEqual
#define lZ3 newbase
#define lY3 branch1op
#define lX3 branch2op
#define lW3 l8 a)xF
#define lV3 overlap
#define lU3 truth_b
#define lT3 truth_a
#define lS3 found_dup
#define lR3 ContainsOtherCandidates
#define lQ3 cR1 lR1
#define lP3 rangeutil
#define lO3 Plan_Has(
#define lN3 StackMax)
#define lM3 const nV2
#define lL3 namespace
#define lK3 ::res,b8<
#define lJ3 ){sim.Eat(
#define lI3 ,bool abs)
#define lH3 Rehash();
#define lG3 ByteCode[
#define lF3 inverted
#define lE3 eO3:
#define lD3 ;}case
#define lC3 fp_mod(m
#define lB3 return p
#define lA3 iftree
#define l93 depcodes
#define l83 explicit
#define l73 cCosh nQ
#define l63 t31 nH
#define l53 VarBegin
#define l43 eH3.
#define l33 ].data);
#define l23 eL)));x0
#define l13 nC==cNot||
#define l03 ),child);
#define iZ2 cV1(-1)))xB
#define iY2 ;break;
#define iX2 )iY2}
#define iW2 ;std::cout<<
#define iV2 ?0:1))l7
#define iU2 .size();++
#define iT2 (yD3 l8
#define iS2 for(xG3
#define iR2 begin(),
#define iQ2 cond_add
#define iP2 cond_mul
#define iO2 cond_and
#define iN2 func lM1
#define iM2 const eH
#define iL2 bool eH1
#define iK2 (xF3
#define iJ2 Optimize()
#define iI2 costree
#define iH2 sintree
#define iG2 leaf_count
#define iF2 sub_params
#define iE2 .n_int_sqrt
#define iD2 );eY=!eY;}
#define iC2 xD3.
#define iB2 printf(
#define iA2 swap(tmp);
#define i92 cbrt_count
#define i82 sqrt_count
#define i72 eX3);
#define i62 t5 1)tF1&&
#define i52 PlusInf
#define i42 Finite
#define i32 (e63
#define i22 c71 n3 0),
#define i12 fp_abs(max.val))
#define i02 fp_abs(yQ)
#define tZ2 p1 tM ifp1
#define tY2 pcall_tree
#define tX2 after_powi
#define tW2 GetHash().
#define tV2 else{x5=new
#define tU2 cO1.SubTrees
#define tT2 cO1.Others
#define tS2 cO1.Immeds
#define tR2 )xD2 t63
#define tQ2 yJ false;}
#define tP2 eO3)cL
#define tO2 cEqual t41
#define tN2 cLog nQ
#define tM2 TreeCountItem
#define tL2 ;if(op==
#define tK2 ,o);o<<"\n";
#define tJ2 ,(long double)
#define tI2 cNeg,lT 1,
#define tH2 ),Value
#define tG2 ),0},{
#define tF2 xC cond nC
#define tE2 (half&63)-1;
#define tD2 tB2 eO3 xQ2
#define tC2 tB2 IsAlways;if(
#define tB2 ))return
#define tA2 ;pow.lH3
#define t92 ,cMul l3
#define t82 tree xC
#define t72 tR2;
#define t62 {std::cout<<
#define t52 tree nC
#define t42 MakeNEqual
#define t32 Rehash()iG
#define t22 Dump(std::
#define t12 isInteger(
#define t02 Comparison
#define eZ2 needs_flip
#define eY2 value]
#define eX2 l51 opcode
#define eW2 ~size_t(0)
#define eV2 xL1 xU+1);
#define eU2 const std::eI
#define eT2 const char*
#define eS2 Rule&rule,
#define eR2 ,const cT&
#define eQ2 tM tree);
#define eP2 mul_item
#define eO2 innersub
#define eN2 cbrt_cost
#define eM2 best_cost
#define eL2 >(cV1(1),
#define eK2 switch(lZ1
#define eJ2 }switch eI2
#define eI2 (tree
#define eH2 condition
#define eG2 per_item
#define eF2 item_type
#define eE2 first2
#define eD2 l4 18,1,
#define eC2 cIf,iO 3,
#define eB2 lJ 1},0,
#define eA2 tH 1},0,
#define e92 Decision
#define e82 not_tree
#define e72 Become(t5
#define e62 group_by
#define e52 eX3=
#define e42 ->second
#define e32 targetpos
#define e22 ParamSpec
#define e12 ,t51 0x1 nH
#define e02 ,unsigned
#define cZ2 rhs.hash2;}
#define cY2 rhs.hash1
#define cX2 tree nV
#define cW2 struct
#define cV2 Forget()
#define cU2 &&cond e5))
#define cT2 source_tree
#define cS2 nC==cLog2&&
#define cR2 <tN,xC3>
#define cQ2 CodeTree lW
#define cP2 p1_evenness
#define cO2 isNegative(
#define cN2 nS){cV1 tmp=
#define cM2 cO known&&(
#define cL2 neg_set
#define cK2 cNop,cNop}}
#define cJ2 cTanh,cNop,
#define cI2 >cW2 c6<
#define cH2 matches
#define cG2 IsAlways)cL
#define cF2 .match_tree
#define cE2 cR1 void*)&
#define cD2 cGreaterOrEq xH
#define cC2 ,cGreaterOrEq
#define cB2 cGreater,cD
#define cA2 cTan nQ
#define c92 cCos nQ
#define c82 (std::move(
#define c72 }data;data.
#define c62 +=1 iF nB1;
#define c52 negated
#define c42 Specializer
#define c32 params
#define c22 coshtree
#define c12 sinhtree
#define c02 best_score
#define yZ2 mulvalue
#define yY2 pow_item
#define yX2 subgroup
#define yW2 )lS 3*3*
#define yV2 if(list.first
#define yU2 nC==cPow&&tF
#define yT2 PowiResult
#define yS2 );nD lC
#define yR2 0));yG3
#define yQ2 maxValue
#define yP2 minValue
#define yO2 fp_min(yL,
#define yN2 div_tree
#define yM2 pow_tree
#define yL2 cV1(0.0)){xO
#define yK2 preserve
#define yJ2 PullResult()
#define yI2 );sim.nW 2,
#define yH2 ;DumpTree(
#define yG2 dup_or_fetch
#define yF2 nominator]
#define yE2 test_order
#define yD2 vector<bool>
#define yC2 ):start_at()
#define yB2 .param_count
#define yA2 minimum_need
#define y92 shift(index)
#define y82 rulenumber
#define y72 middle2
#define y62 cTanh nQ
#define y52 cSinh nQ
#define y42 cInv,lT 1,
#define y32 constraints=
#define y22 factor_immed
#define y12 changes
#define y02 tM leaf2 l8
#define xZ2 tM yD3 l8
#define xY2 tM cond l8
#define xX2 exp_diff
#define xW2 ExponentInfo
#define xV2 lower_bound(
#define xU2 factor
#define xT2 is_logical
#define xS2 newrel_and
#define xR2 ;i7.Remember(
#define xQ2 iF Unknown;}
#define xP2 res_stackpos
#define xO2 half_pos
#define xN2 ifdata.ofs
#define xM2 >>1)):(
#define xL2 CodeTreeData
#define xK2 ;sim.Push(
#define xJ2 =GetParam
#define xI2 multiply(
#define xH2 (IfData&ifdata
#define xG2 i4 push_back(
#define xF2 i4 size()
#define xE2 var_trees
#define xD2 parampair
#define xC2 MakeFalse
#define xB2 return false;}
#define xA2 :sim.Eat(1,
#define x92 const lR1
#define x82 e73 PushImmed(
#define x72 );e73 DoDup iQ
#define x62 )iW2
#define x52 );m cO n3 0),
#define x42 ,2,1 xM if iQ
#define x32 )eN xF3);
#define x22 );SetParamMove(
#define x12 ;tree c51
#define x02 .DelParam(
#define nZ2 ;nA1 tM y7 l8
#define nY2 )){data x7
#define nX2 iF true;}
#define nW2 nA OPCODE
#define nV2 CodeTree&
#define nU2 parent_opcode
#define nT2 log2_exponent
#define nS2 dup_fetch_pos
#define nR2 {cT start_at;
#define nQ2 cSin nQ
#define nP2 Value_EvenInt
#define nO2 AddCollection
#define nN2 ConditionType
#define nM2 SpecialOpcode
#define nL2 ByteCodeSynth
#define nK2 AddParamMove(
#define nJ2 .IsDefined(
#define nI2 fp_max(yL);
#define nH2 cO known&&p
#define nG2 assimilated
#define nF2 denominator
#define nE2 fraction
#define nD2 .GetDepth()
#define nC2 DUP_BOTH();
#define nB2 xE3 lL
#define nA2 0x80000000u
#define n92 if(e73 Find(
#define n82 IsDescendantOf
#define n72 TreeCounts
#define n62 bool eY=false;
#define n52 SetOpcode(
#define n42 found_log2
#define n32 div_params
#define n22 immed_sum
#define n12 lG3++IP]
#define n02 OPCODE(opcode)
#define lZ2 break;xD3*=
#define lY2 FactorStack xK
#define lX2 Rehash(false);
#define lW2 a;if(tZ1){x5=(
#define lV2 cLessOrEq,
#define lU2 282870 nP
#define lT2 cNotNot nQ
#define lS2 cNot nQ
#define lR2 replacing_slot
#define lQ2 RefParams
#define lP2 ,{l5::
#define lO2 if_always[
#define lN2 WhatDoWhenCase
#define lM2 exponent_immed
#define lL2 new_base_immed
#define lK2 base_immed
#define lJ2 ||op1==
#define lI2 data[a t53
#define lH2 if(newrel_or==
#define lG2 TopLevel)
#define lF2 ):Value(Value::
#define lE2 .AddOperation(
#define lD2 DUP_ONE(apos);
#define lC2 flipped
#define lB2 .UseGetNeeded(
#define lA2 e0 2,131,
#define l92 [xU-1-offset].
#define l82 lG3 a
#define l72 y8 Immed.size());
#define l62 1 y8 xF2
#define l52 OptimizedUsing
#define l42 Var_or_Funcno
#define l32 l42;
#define l22 GetParams(
#define l12 crc32_t
#define l02 signed_chain
#define iZ1 MinusInf
#define iY1 n_immeds
#define iX1 stack.size()
#define iW1 FindClone(xQ
#define iV1 lG3 IP]
#define iU1 GetOpcode())
#define iT1 needs_rehash
#define iS1 AnyWhere_Rec
#define iR1 ~unsigned(0)
#define iQ1 41,42,43,44,
#define iP1 p1_logical_b
#define iO1 p0_logical_b
#define iN1 p1_logical_a
#define iM1 p0_logical_a
#define iL1 {if(GetOpcode()
#define iK1 ,PowiCache&i7,
#define iJ1 xC t52);
#define iI1 *const func)
#define iH1 else if(
#define iG1 e73 DoDup(
#define iF1 cache_needed
#define iE1 e0 2,1,e0 2,
#define iD1 treelist
#define iC1 has_bad_balance
#define iB1 e93 xU2
#define iA1 double)eX3
#define i91 {case IsAlways:
#define i81 cNEqual
#define i71 tH 2},0,0x0},{{
#define i61 Oneness_NotOne|
#define i51 Value_IsInteger
#define i41 Constness_Const
#define i31 DumpHashesFrom(
#define i21 l52(
#define i11 reltype
#define i01 SequenceOpcodes
#define tZ1 &*start_at
#define tY1 (*x5 lS1=r.specs;if(r.found){
#define tX1 )return false
#define tW1 sep_list[
#define tV1 (*x5)[a].info
#define tU1 if(remaining[a])
#define tT1 grammar_rules[*r]
#define tS1 (size_t
#define tR1 for tS1 b=0;b<
#define tQ1 ;e73 StackTopIs(
#define tP1 goto fail;}
#define tO1 l1 0x4 nH
#define tN1 xE3<
#define tM1 lQ2);
#define tL1 ;pow xC cPow);pow
#define tK1 xE3 set_if<
#define tJ1 xE3 lX
#define tI1 TreeCountType xK
#define tH1 2*2*2)lS 3
#define tG1 c71 set(fp_floor
#define tF1 .IsImmed()
#define tE1 a)tF1)
#define tD1 n02);
#define tC1 std::cout<<"POP "
#define tB1 stack[iX1-
#define tA1 stack.push_back(
#define t91 MaxChildDepth
#define t81 repl_param_list,
#define t71 std::pair<It,It>
#define t61 );tree x02
#define t51 cPow,lA
#define t41 ,l0 2,
#define t31 ,l1 0x12
#define t21 Sign_Negative
#define t11 Value_Logical
#define t01 new_factor_immed
#define eZ1 eI2,std::cout);
#define eY1 ));nA1 tO op1);eP1
#define eX1 }break lD3
#define eW1 <<tree.tW2
#define eV1 tJ1 void
#define eU1 cO known)
#define eT1 e73 xI 1
#define eS1 ,i7 nF3
#define eR1 ByteCode,size_t&IP,size_t limit,size_t y4
#define eQ1 tree.SetParam(
#define eP1 tree.DelParams()
#define eO1 occurance_pos
#define eN1 exponent_hash
#define eM1 exponent_list
#define eL1 CollectionSet xK
#define eK1 CollectMulGroup(
#define eJ1 source_set
#define eI1 eX3,y33
#define eH1 operator
#define eG1 FindAndDup eI2);
#define eF1 back().thenbranch
#define eE1 retry_anyparams_3
#define eD1 retry_anyparams_2
#define eC1 needlist_cached_t
#define eB1 yY3 0x4},{{1,
#define eA1 eB2 0x4},{{1,
#define e91 CodeTreeImmed xK(
#define e81 lL3 FPoptimizer_Optimize
#define e71 GetParamCount()==
#define e61 by_float_exponent
#define e51 fp_equal eW3
#define e41 new_exp
#define e31 n72.erase(cs_it);
#define e21 tB2 true;
#define e11 ){CodeTree xK xQ
#define e01 ,l1 0x0},{{3,
#define cZ1 end()&&i->first==
#define cY1 return BecomeZero;
#define cX1 return BecomeOne;
#define cW1 if(lQ.size()<=n1)
#define cV1 Value_t
#define cU1 addgroup
#define cT1 found_log2by
#define cS1 =comp.AddItem(atree
#define cR1 (const
#define cQ1 cR1 cV1&
#define cP1 nC==cV3)
#define cO1 NeedList
#define cN1 if(keep_powi
#define cM1 ParsePowiMuli(
#define cL1 l42)
#define cK1 MakeNotP1,l5::
#define cJ1 MakeNotP0,l5::
#define cI1 branch1_backup
#define cH1 divgroup
#define cG1 l21 a-->0;)if(
#define cF1 branch2_backup
#define cE1 exponent_map
#define cD1 plain_set
#define cC1 rangehalf
#define cB1 LightWeight(
#define cA1 if(value
#define c91 tJ1 yW
#define c81 tJ1 static
#define c71 :lD m.min.
#define c61 xF3.
#define c51 .nK2
#define c41 should_regenerate=true;
#define c31 should_regenerate,
#define c21 Collection
#define c11 CodeTree xK r;r xC
#define c01 RelationshipResult
#define yZ1 Subdivide_Combine(
#define yY1 long value
#define yX1 )iF m lD3
#define yW1 public cZ,public std::vector<
#define yV1 ):cZ(),std::vector<
#define yU1 )const yJ
#define yT1 rhs yU1 hash1
#define yS1 for tS1 a y2
#define yR1 best_sep_factor
#define yQ1 iH1!xD3
#define yP1 needlist_cached
#define yO1 inline unsigned
#define yN1 opcode,bool pad
#define yM1 changed=true;
#define yL1 DelParam(a);}
#define yK1 MakesInteger(
#define yJ1 const cV1&value
#define yI1 best_sep_cost
#define yH1 MultiplicationRange
#define yG1 pihalf_limits
#define yF1 ;NewHash.hash2+=
#define yE1 y6 2,cPow);lC
#define yD1 n_stacked
#define yC1 NewHash.hash1
#define yB1 AnyParams_Rec
#define yA1 cGreaterOrEq,
#define y91 continue;
#define y81 Become(value l8 0))
#define y71 PositionalParams,0}
#define y61 always_sincostan
#define y51 Recheck_RefCount_Div
#define y41 Recheck_RefCount_Mul
#define y31 xF3;xF3 xC
#define y21 MultiplyAndMakeLong(
#define y11 cV1(0)
#define y01 covers_plus1
#define xZ1 if(e73 FindAndDup(
#define xY1 SynthesizeParam(
#define xX1 grammar_func
#define xW1 cOr l3 16,1,
#define xV1 252421 nP 24830,
#define xU1 xH 529654 nP
#define xT1 Modulo_Radians},
#define xS1 PositionType
#define xR1 CollectionResult
#define xQ1 tJ1 bool
#define xP1 const_offset
#define xO1 inline TriTruthValue
#define xN1 stacktop_desired
#define xM1 int mStackPtr=0;
#define xL1 SetStackTop(
#define xK1 FPoptimizer_ByteCode
#define xJ1 1)?(poly^(
#define xI1 ,nL2 xK&synth)
#define xH1 y11)
#define xG1 xF leaf2 l8
#define xF1 ;for tS1 a=
#define xE1 cond_type
#define xD1 fphash_value_t
#define xC1 Recheck_RefCount_RDiv
#define xB1 cMul);tmp eU3;tmp.
#define xA1 SwapLastTwoInStack();
#define x91 ParamSpec_Extract xK(
#define x81 fPExponentIsTooLarge(
#define x71 CollectMulGroup_Item(
#define x61 pair<cV1,y33>
#define x51 nL xL1 xU-1);
#define x41 covers_full_cycle
#define x31 xK()y6 2,cMul);lC
#define x21 AssembleSequence(
#define x11 252180 nP 281854,
#define x01 {DataP slot_holder(y3[
#define nZ1 <<std::dec<<")";}
#define nY1 c23 165888 nP
#define nX1 }lP2 MakeNotNotP0,l5::
#define nW1 !=xA)if(TestCase(
#define nV1 &&IsLogicalValue(
#define nU1 AddFunctionOpcode
#define nT1 }inline
#define nS1 std::pair<T1,T2>&
#define nR1 tN1 xG3
#define nQ1 has_good_balance_found
#define nP1 n_occurrences
#define nO1 found_log2_on_exponent
#define nN1 tree x02 a
#define nM1 tmp c51 tree);
#define nL1 ,cIf,l0 3,
#define nK1 CodeTree xK tmp;tmp xC
#define nJ1 covers_minus1
#define nI1 needs_resynth
#define nH1 immed_product
#define nG1 yP3 a=0;a<xY;++a)
#define nF1 Sign_Positive
#define nE1 {l5::MakeNotNotP1,l5::
#define nD1 ::MakeTrue
#define nC1 CodeTreeImmed(cV1(
#define nB1 Suboptimal
#define nA1 changed_if
#define n91 n_as_tanh_param
#define n81 );p2 yM p2);t82 lA3 nC);i3}
#define n71 opposite=
#define n61 xD1(
#define n51 MatchResultType
#define n41 yJ CodeTree xK(
#define n31 needs_sincos
#define n21 resulting_exponent
#define n11 ;tJ1
#define n01 Unknown:yM3;}
#define lZ1 GetLogicalValue(t5
#define lY1 ,l4 2,1,
#define lX1 GetParam(a)
#define lW1 inverse_nominator]
#define lV1 cMul l3 0,1,
#define lU1 cSin lY1
#define lT1 GetParamCount())
#define lS1 )[a].start_at
#define lR1 CodeTree xK&
#define lQ1 void FunctionParserBase
#define lP1 xC cLog);t82 cMul);
#define lO1 SetParams(l22));
#define lN1 o<<"("<<std::hex<<data.
#define lM1 (val);else*this=model;}
#define lL1 IfBalanceGood(
#define lK1 n_as_tan_param
#define lJ1 0,i72 DelParam(1);
#define lI1 changed_exponent
#define lH1 yK3<cV1(
#define lG1 inverse_denominator
#define lF1 retry_positionalparams_2
#define lE1 tF1){cV1
#define lD1 );eX3.lH3
#define lC1 )){CodeTree xK
#define lB1 tF1 tX1;
#define lA1 situation_flags&
#define l91 518 nP 400412,
#define l81 CopyOnWrite();
#define l71 PlanNtimesCache(
#define l61 FPoptimizer_Grammar
#define l51 (unsigned
#define l41 l51 index
#define l31 .GetImmed()
#define l21 GetParamCount();
#define l11 AddOperation(cInv,1,1 xM}
#define l01 GetPositivityInfo eI2)!=
#define iZ recursioncount
#define iY ParamSpec_SubFunctionData
#define iX PositionalParams_Rec
#define iW .data.subfunc_opcode
#define iV lG3 xN2+
#define iU DumpTreeWithIndent(*this);
#define iT CalculateResultBoundaries(
#define iS tN1 unsigned Compare>
#define iR >=xH1
#define iQ (found[data.
#define iP tJ1 nB
#define iO lA 0x4},{{
#define iN .min xB3
#define iM yY3 0x0},{{1,
#define iL edited_powgroup
#define iK has_unknown_max
#define iJ has_unknown_min
#define iI {switch(type xI3 cond_or:
#define iH (a);bool needs_cow=GetRefCount()>1;
#define iG ;CodeTree xK
#define iF ;return
#define iE set(fp_ceil yX1
#define iD nM cV1(-y93
#define iC static const yG3
#define iB synthed_tree
#define iA 7168 nP 401798,
#define i9 SelectedParams,0},0,0x0},{{
#define i8 collections
#define i7 cache
#define i6 ;xG2 nA2
#define i5 )i6);
#define i4 ByteCode.
#define i3 goto redo;
#define i2 goto ReplaceTreeWithOne tJ3
#define i1 ]);nG3
#define i0 !=xA)return lO2
#define tZ e61.data
#define tY l83 xL2(
#define tX needs_sinhcosh
#define tW cAdd l3 0,
#define tV tG2 cV1(
#define tU xC2,l5::
#define tT by_exponent
#define tS x62 std::endl;DumpHashes(
#define tR ,ByteCode,IP,limit,y4,stack);
#define tQ 408964 nP 24963,
#define tP 528503 nP 24713,
#define tO .lH3 t82
#define tN int_exponent_t
#define tM .AddParam(
#define tL 0;a<tree.l21++a)
#define tK matched_params
#define tJ [n1 t43=true;lQ[n1 t53
#define tI l61::Grammar*
#define tH xD AnyParams,
#define tG relationships
#define tF powgroup l8
#define tE eW2&&found[data.
#define tD }},{ProduceNewTree,2,1,
#define tC nC1(
#define tB has_mulgroups_remaining
#define tA const iY
#define t9 MatchInfo xK&
#define t8 lH3 iF2.push_back(
#define t7 best_factor
#define t6 RootPowerTable xK::RootPowers[
#define t5 tree l8
#define t4 :{AdoptChildrenWithSameOpcode eI2);
#define t3 :goto ReplaceTreeWithZero tJ3
#define t2 MatchPositionSpec_AnyParams xK
#define t1 lL3 FPoptimizer_CodeTree
#define t0 n_as_sinh_param
#define eZ n_as_cosh_param
#define eY is_signed
#define eX result_positivity
#define eW biggest_minimum
#define eV 122999 nP 139399,
#define eU 142455 nP 141449,
#define eT cond_tree
#define eS else_tree
#define eR then_tree
#define eQ valueType
#define eP (p0 iN&&p0.yQ>=cV1(0.0))
#define eO }yV2 l31==cV1(
#define eN ;nK2
#define eM =iT t5
#define eL t5 1)l31
#define eK lR1 tree)
#define eJ sequencing
#define eI string i73(
#define eH std::vector<CodeTree xK>
#define eG if_stack
#define eF (l22));c61 lH3
#define eE n_as_sin_param
#define eD n_as_cos_param
#define eC PowiResolver::
#define eB ;tree.SetParamMove(
#define eA cIf,tR3
#define e9 ].relationship
#define e8 cO known=false;
#define e7 const ParamSpec_SubFunction
#define e6 const ParamSpec_ParamHolder
#define e5 .BalanceGood
#define e4 nK2 yX2
#define e3 back().endif_location
#define e2 xD1 key
#define e1 nK2 mul);
#define e0 130,1,
#define cZ MatchPositionSpecBase
#define cY l83 CodeTree(
#define cX smallest_maximum
#define cW }PACKED_GRAMMAR_ATTRIBUTE;
#define cV ReplaceTreeWithParam0;
#define cU factor_needs_rehashing
#define cT MatchPositionSpecBaseP
#define cS xG3 tI1::y53
#define cR x91 nN.param_list,
#define cQ cV1(1.5)*fp_const_pi xK()
#define cP cO known&&p0 cO val<=fp_const_negativezero xK())
#define cO .max.
#define cN yH2 tree x62"\n";
#define cM return false iY2}switch(bitmask&
#define cL return false;lC
#define cK tS1 a=tree.l21 a-->0;)
#define cJ 243,244,245,246,249,250,251,253,255,256,257,258,259}};}
#define cI ];};extern"C"{
#define cH 79,122,123,160,161,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,
#define cG 27,28,29,30,31,32,33,35,36,
#define cF otherhalf
#define cE StackState
#define cD l2 16,2,
#define cC const SequenceOpCode xK
#define cB paramholder_matches[x1]
#define cA MatchPositionSpec_PositionalParams xK
#define c9 x92 tree,std::ostream&o
#define c8 paramholder_matches.
#define c7 CalculatePowiFactorCost(
#define c6 ImmedHashGenerator
#define c5 ::map<fphash_t,std::set<std yW3> >
#define c4 ComparisonSetBase::
#define c3 nK2 comp.cD1[a].value);
#define c2 T1,xG3 T2>inline iL2()(
#define c1 has_nonlogical_values
#define c0 from_logical_context)
#define yZ AnyParams,0}},{ProduceNewTree,
#define yY for tS1 a=y1.l21 a-->0;)
#define yX POWI_CACHE_SIZE
#define yW static inline CodeTree xK
#define yV ++IP;y91}if(iV1==cN3.
#define yU }lP2 xA,l5::Never}lP2 xA,l5::Never}}
#define yT .FoundChild
#define yS BalanceResultType
#define yR );void nU1 eX2,c42<
#define yQ min.val
#define yP xJ3(0),Opcode(
#define yO CodeTree xK tmp,tmp2;tmp2 xC
#define yN +=fp_const_twopi xK();
#define yM .Rehash()x12
#define yL fp_sin(min),fp_sin(max))
#define yK fp_const_twopi xK()yU3
#define yJ {return
#define yI const yJ data->
#define yH lE2 GetOpcode(),
#define yG for tS1 a=0;a<l21++a){if(
#define yF MatchPositionSpec_AnyWhere
#define yE yX3.data.match_type==
#define yD void OutFloatHex(std::ostream&o,
#define yC {static void n83 nA fphash_t&NewHash,
#define yB AddParam(CodeTreeImmed(
#define yA ,xG3 CodeTree xK::
#define y9 AssembleSequence_Subdivide(
#define y8 ]=nA2|unsigned(
#define y7 branch2
#define y6 );sim.Eat(
#define y5 unsigned c i23 short l[
#define y4 factor_stack_base
#define y3 data->eH3
#define y2 =0;a<y23.l21++a)if(
#define y1 branch1
#define y0 .SetParam(0,lA3 l8 0))iG p1;p1 xC
#define xZ (n73 r=range.first;r!=range t63;++r){
#define xY nN yB2
#define xX !=eW2){e73 DoDup iQ
#define xW =fp_cosh(m.yQ);m cO val=fp_cosh cU3);
#define xV {n72.erase(i);y91}
#define xU StackTop
#define xT FPOPT_autoptr
#define xS +=xD3 iF xD3;}tJ1 inline cV1
#define xR int_exponent
#define xQ newnode
#define xP has_highlevel_opcodes
#define xO tree.ReplaceWithImmed(
#define xN ,t51 0x4 nH
#define xM )tQ1*this)iF;}
#define xL sim.AddConst(
#define xK <cV1>
#define xJ xK model=cC1 xK()){if(known
#define xI GetStackTop()-
#define xH ,l2 18,2,
#define xG eB2 0x0},{{
#define xF .IsIdenticalTo(
#define xE t92 2,1,
#define xD ,cAdd,
#define xC .n52
#define xB {if(needs_cow){l81 goto
#define xA Unchanged
#define x9 c9=std::cout
#define x8 best_selected_sep
#define x7 ->Recalculate_Hash_NoRecursion();}
#define x6 l21++a)if(ApplyGrammar(xX3,t73,
#define x5 position
#define x4 for tS1 a=tL{yG3
#define x3 std::vector<CodeTree>
#define x2 TestImmedConstraints l64 constraints,tree)tX1;
#define x1 paramholder_index
#define x0 return true tJ3
#define nZ occurance_counts
#define nY ;tmp2 eU3;tmp xC cInv);tmp c51 tmp2)iF
#define nX -->0;){x92 powgroup=lX1;if(powgroup
#define nW SwapLastTwoInStack(y6
#define nV )){tree.FixIncompleteHashes();}
#define nU cV1>p eM a)yU3 p.
#define nT );cN1 lJ3 1,cInv iX2 xL-1 yE1
#define nS if(t5 0)tF1
#define nR ,cAdd l3 2,1,
#define nQ ,l0 1,
#define nP ,{2,
#define nO const FPoptimizer_CodeTree::lR1 tree
#define nN model_tree
#define nM return yG3(
#define nL ){using lL3 FUNCTIONPARSERTYPES;
#define nK eH&lQ2
#define nJ AddParam(t5
#define nI ConstantFolding_LogicCommon eI2,c4
#define nH },{{2,
#define nG nR1 Ref>inline void xT<Ref>::
#define nF AnyParams,1},0,0x0},{{
#define nE ):data(new xL2 xK(
#define nD goto do_return;}
#define nC .GetOpcode()
#define nB xL2 xK::xL2(
#define nA FUNCTIONPARSERTYPES::
#define n9 b;}};tN1>cW2 Comp<nA
#define n8 l42(),eH3(),Hash(),Depth(1),i21 0){}
#define n7 SynthesizeByteCode(synth);
#define n6 while(ApplyGrammar(cE2
#define n5 GetIntegerInfo(t5 0))==IsAlways)cY3
#define n4 x12 nA1)nX2
#define n3 tK1 cGreater>(cV1(
#define n2 DumpParams xK l64 data.param_list,l94 yB2,o);
#define n1 restholder_index
#define n0 CodeTree xK eX3;eX3 tM3 eX3 tM
#define lZ lR yU3 fp_nequal(tmp,xH1){xO cV1(1)/tmp);nD}lC
#define lY :if(ParamComparer xK()(eH3[1],eH3[0])){std::swap(eH3[0],eH3[1]);Opcode=
#define lX <xG3 cV1>
#define lW xK tmp;tmp xC cPow);tmp eU3;tmp.yB cV1(
#define lV i41,0x0},
#define lU nK2 pow l8 1));pow x02 1);pow.Rehash()eB 0,pow);goto NowWeAreMulGroup;}
#define lT GroupFunction,0},lV{{
#define lS ,cV1(1)/cV1(
#define lR t5 0)l31
#define lQ restholder_matches
#define lP yC1|=key;xD1 crc=(key>>10)|(key<<(64-10))yF1((~n61 crc))*3)^1234567;}};
#define lO nA1;nA1 iJ1 nA1 c51 t5 0));nA1 tM y1 l8
#define lN tJ1 CodeTree xK::CodeTree(
#define lM eQ1 0,t5 0)cO3 eQ1 1,CodeTreeImmed(
#define lL lX void nL2 xK::nU1 eX2,c42<
#define lK cMul iV3
#define lJ cMul,AnyParams,
#define lI (t5 0)tF1&&t5 1)tF1){xO
#define lH iT tmp)lD3
#define lG :eJ3=comp.AddRelationship(atree l8 0),atree l8 1),c4
#define lF cPow,l0 2
#define lE xG3 cV1>inline iL2()cQ1 xK3 cV1&b)yJ a
#define lD {yG3 m eM 0));
#define lC break tJ3
#define lB eV1 CodeTree xK::
#define lA y71,0,
#define l9 l1 0x0 nH
#define l8 .GetParam(
#define l7 iG nA1;nA1 iJ1 nA1.SetParamsMove eI2.l22));nA1 tO
#define l6 SelectedParams,0},0,0x0 nH
#define l5 RangeComparisonData
#define l4 y71},{ProduceNewTree,
#define l3 ,AnyParams,0}},{ReplaceParams,
#define l2 y71},{ReplaceParams,
#define l1 cMul,SelectedParams,0},0,
#define l0 lA 0x0},{{
#ifdef _MSC_VER
typedef
unsigned
int
l12;
#else
#include <stdint.h>
typedef
uint_least32_t
l12;
#endif
lL3
crc32{enum{startvalue=0xFFFFFFFFUL,poly=0xEDB88320UL}
;tN1
l12
crc>cW2
b8{enum{b1=(crc&xJ1
crc
xM2
crc>>1),b2=(b1&xJ1
b1
xM2
b1>>1),b3=(b2&xJ1
b2
xM2
b2>>1),b4=(b3&xJ1
b3
xM2
b3>>1),b5=(b4&xJ1
b4
xM2
b4>>1),b6=(b5&xJ1
b5
xM2
b5>>1),b7=(b6&xJ1
b6
xM2
b6>>1),res=(b7&xJ1
b7
xM2
b7>>1)}
;}
;inline
l12
update(l12
crc
e02
b){
#define B4(n) b8<n>lK3 n+1>lK3 n+2>lK3 n+3>::res
#define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
static
const
l12
table[256]={R(0x00),R(0x10),R(0x20),R(0x30),R(0x40),R(0x50),R(0x60),R(0x70),R(0x80),R(0x90),R(0xA0),R(0xB0),R(0xC0),R(0xD0),R(0xE0),R(0xF0)}
;
#undef R
#undef B4
return((crc>>8))^table[(crc^b)&0xFF];nT1
l12
calc_upd(l12
c,const
unsigned
char*buf,size_t
size){l12
value=c;for
tS1
p=0;p<size;++p)value=update(value,buf[p])iF
value;nT1
l12
calc
cR1
unsigned
char*buf,size_t
size)yJ
calc_upd(startvalue,buf,size);}
}
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH
nR1
Ref>class
xT{eQ3
xT():p(0){}
xT(Ref*b):p(b){xH3}
xT
cR1
xT&b):p(b.p){xH3
nT1
Ref&eH1*(yU1*p;nT1
Ref*eH1->(yU1
p;}
xT&eH1=(Ref*b){Set(b)iF*this;}
xT&eH1=cR1
xT&b){Set(b.p)iF*this;}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
xT(xT&&b):p(b.p){b.p=0;}
xT&eH1=(xT&&b){if(p!=b.p){cV2;p=b.p;b.p=0;}
return*this;}
#endif
~xT(){cV2
nD3
UnsafeSetP(Ref*newp){p=newp
nD3
swap(xT<Ref>&b){Ref*tmp=p;p=b.p;b.p=tmp;}
private:inline
static
void
Have(Ref*p2);inline
void
cV2;inline
void
xH3
inline
void
Set(Ref*p2);private:Ref*p;}
;nG
cV2{if(!p)return;p->xJ3-=1;if(!p->xJ3)delete
p;}
nG
Have(Ref*p2){if(p2)++(p2->xJ3);}
nG
Birth(){Have(p);}
nG
Set(Ref*p2){Have(p2);cV2;p=p2;}
#endif
#include <utility>
cW2
Compare2ndRev{nR1
T>inline
iL2()cR1
T&xK3
T&b
yU1
a
t63>b
t63;}
}
;cW2
Compare1st{nR1
c2
const
nS1
xK3
nS1
b
yU1
a.first<b.first;}
nR1
c2
const
nS1
a,T1
b
yU1
a.first<b;}
nR1
c2
T1
xK3
nS1
b
yU1
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
lL3
FUNCTIONPARSERTYPES{cW2
fphash_t{xD1
hash1,hash2;fphash_t():hash1(0),hash2(0){}
fphash_t
cR1
xD1&xK3
xD1&b):hash1(a),hash2(b){}
iL2==cR1
fphash_t&yT1==cY2&&hash2==cZ2
iL2!=cR1
fphash_t&yT1!=cY2||hash2!=cZ2
iL2<cR1
fphash_t&yT1!=cY2?hash1<cY2:hash2<cZ2}
;}
#endif
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
lL3
l61{cW2
Grammar;}
lL3
xK1{tJ1
class
nL2;}
t1{tJ1
class
CodeTree
n11
cW2
xL2
n11
class
CodeTree{typedef
xT<xL2
xK>DataP;DataP
data;eQ3
CodeTree();~CodeTree();cW2
OpcodeTag{}
;cY
nW2
o,OpcodeTag);cW2
FuncOpcodeTag{}
;cY
nW2
o
e02
f,FuncOpcodeTag);cW2
xM3{}
;cY
const
cV1&v,xM3);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
cY
cV1&&v,xM3);
#endif
cW2
VarTag{}
;cY
unsigned
varno,VarTag);cW2
CloneTag{}
;cY
lM3
b,CloneTag);void
GenerateFrom
cR1
xG3
FunctionParserBase
xK::Data&data,bool
keep_powi=false);void
GenerateFrom
cR1
xG3
FunctionParserBase
xK::Data&data,const
x3&xE2,bool
keep_powi=false);void
SynthesizeByteCode(std::vector<unsigned>&eF3,std::vector
xK&immed,size_t&stacktop_max);void
SynthesizeByteCode(xK1
cL3&synth,bool
MustPopTemps=true)const;size_t
SynthCommonSubExpressions(xK1
cL3&synth)const;void
SetParams
cR1
x3&xN3
SetParamsMove(x3&tM1
CodeTree
GetUniqueRef();
#ifdef __GXX_EXPERIMENTAL_CXX0X__
void
SetParams(x3&&tM1
#endif
void
SetParam
tS1
which,lM3
b);void
SetParamMove
tS1
which,nV2
b);void
AddParam
cR1
nV2
param);void
nK2
nV2
param);void
AddParams
cR1
x3&xN3
AddParamsMove(x3&xN3
AddParamsMove(x3&lQ2,size_t
lR2);void
DelParam
tS1
index);void
DelParams();void
Become
cR1
nV2
b);inline
size_t
GetParamCount(yU1
l22).size();nT1
nV2
GetParam
tS1
n)yJ
l22)[n];nT1
lM3
GetParam
tS1
n
yU1
l22)[n];nT1
void
n52
nW2
o)eP3
Opcode=o;nT1
nW2
GetOpcode()yI
Opcode;nT1
nA
fphash_t
GetHash()yI
Hash;nT1
const
x3&l22
yU1
y3;nT1
x3&l22)yJ
y3;nT1
size_t
GetDepth()yI
Depth;nT1
const
cV1&GetImmed()yI
Value;nT1
unsigned
GetVar()yI
l32
nT1
unsigned
GetFuncNo()yI
l32
nT1
bool
IsDefined(yU1
GetOpcode()!=nA
cNop;nT1
bool
IsImmed(yU1
GetOpcode()==nA
cImmed;nT1
bool
IsVar(yU1
GetOpcode()==nA
l53;nT1
unsigned
GetRefCount()yI
xJ3
nD3
ReplaceWithImmed
cQ1
i);void
Rehash(bool
constantfolding=true);void
Sort();inline
void
Mark_Incompletely_Hashed()eP3
Depth=0;nT1
bool
Is_Incompletely_Hashed()yI
Depth==0;nT1
const
tI
GetOptimizedUsing()yI
l52;nT1
void
SetOptimizedUsing
cR1
tI
g)eP3
l52=g;}
bool
RecreateInversionsAndNegations(bool
prefer_base2=false);void
FixIncompleteHashes();void
swap(nV2
b){data.swap(b.data);}
bool
IsIdenticalTo
cR1
nV2
b)const;void
l81}
n11
cW2
xL2{int
xJ3;nW2
Opcode;cV1
Value
i23
l32
eH
eH3;nA
fphash_t
Hash;size_t
Depth;const
tI
l52;xL2();xL2
cR1
xL2&b);tY
nW2
o);tY
nW2
o
e02
f);tY
const
cV1&i);
#ifdef __GXX_EXPERIMENTAL_CXX0X__
tY
cV1&&i);xL2(xL2&&b);
#endif
bool
IsIdenticalTo
cR1
xL2&b)const;void
Sort();void
Recalculate_Hash_NoRecursion();private:void
eH1=cR1
xL2&b);}
n11
yW
CodeTreeImmed
cQ1
i)n41
i
yA
xM3());}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
c91
CodeTreeImmed(cV1&&i)n41
std::move(i)yA
xM3());}
#endif
c91
CodeTreeOp(nW2
opcode)n41
opcode
yA
OpcodeTag());}
c91
CodeTreeFuncOp(nW2
opcode
e02
f)n41
opcode,f
yA
FuncOpcodeTag());}
c91
CodeTreeVar
l51
varno)n41
varno
yA
VarTag());}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
eV1
DumpHashes(x9)xP3
DumpTree(x9)xP3
DumpTreeWithIndent(x9,const
std
yW3&indent="\\"
);
#endif
}
#endif
#endif
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH
#include <iostream>
t1{tJ1
class
CodeTree;}
lL3
l61{enum
ImmedConstraint_Value{ValueMask=0x07,Value_AnyNum=0x0,nP2=0x1,Value_OddInt=0x2,i51=0x3,Value_NonInteger=0x4,t11=0x5
x73
ImmedConstraint_Sign{SignMask=0x18,Sign_AnySign=0x00,nF1=0x08,t21=0x10,Sign_NoIdea=0x18
x73
ImmedConstraint_Oneness{OnenessMask=0x60,Oneness_Any=0x00,Oneness_One=0x20,Oneness_NotOne=0x40
x73
ImmedConstraint_Constness{ConstnessMask=0x180,Constness_Any=0x00,i41=0x80,Constness_NotConst=0x100
x73
Modulo_Mode{Modulo_None=0,Modulo_Radians=1
x73
Situation_Flags{LogicalContextOnly=0x01,NotForIntegers=0x02,OnlyForIntegers=0x04,OnlyForComplex=0x08,NotForComplex=0x10
x73
nM2{NumConstant,ParamHolder,SubFunction
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
std::pair<nM2,const
void*>e22
n11
e22
ParamSpec_Extract
l51
paramlist
e02
index)n11
bool
ParamSpec_Compare
cR1
void*xK3
void*b,nM2
type)i23
ParamSpec_GetDepCode
cR1
e22&b);cW2
ParamSpec_ParamHolder{unsigned
index:8
i23
constraints:9
i23
depcode:15;cW
tJ1
cW2
ParamSpec_NumConstant{cV1
constvalue
i23
modulo;}
;cW2
iY{unsigned
param_count:2
i23
param_list:30;nW2
subfunc_opcode:8;ParamMatchingType
match_type:3
i23
n1:5;cW
cW2
ParamSpec_SubFunction{iY
data
i23
constraints:9
i23
depcode:7;cW
cW2
Rule{RuleType
ruletype:2
i23
situation_flags:5
i23
repl_param_count:2+9
i23
repl_param_list:30;iY
match_tree;cW
cW2
Grammar{unsigned
rule_count
i23
short
rule_list[999
cI
extern
const
Rule
grammar_rules[];}
eV1
DumpParam
cR1
e22&p,std::ostream&o=std::cout)xP3
DumpParams
l51
paramlist
e02
count,std::ostream&o=std::cout);}
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define CONSTANT_POS_INF HUGE_VAL
#define CONSTANT_NEG_INF (-HUGE_VAL)
lL3
FUNCTIONPARSERTYPES{tJ1
inline
cV1
fp_const_pihalf()yJ
fp_const_pi
xK()*cP3;}
tJ1
inline
cV1
fp_const_twopi()i43
fp_const_pi
xK());xD3
xS
fp_const_twoe()i43
fp_const_e
xK());xD3
xS
fp_const_twoeinv()i43
fp_const_einv
xK());xD3
xS
fp_const_negativezero()yJ-Epsilon
xK::value;}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#include <iostream>
e81{using
lL3
l61;using
t1;using
lL3
FUNCTIONPARSERTYPES
n11
class
MatchInfo{eQ3
std::vector<std::pair<bool,eH> >lQ;eH
paramholder_matches;std::vector<unsigned>tK;eQ3
MatchInfo():lQ(),paramholder_matches(),tK(){}
eQ3
bool
SaveOrTestRestHolder
l51
n1,iM2&iD1){cW1{lQ.eV3
n1+1);lQ
tJ=iD1
nX2
if(lQ[n1
t43==false){lQ
tJ=iD1
nX2
iM2&found=lQ[n1
t53;if(iD1.size()!=found.size()tX1
xF1
0;a<iD1
iU2
a)if(!iD1[a]xF
found[a])tX1
nX2
void
SaveRestHolder
l51
n1,eH&iD1){cW1
lQ.eV3
n1+1);lQ
tJ.swap(iD1);}
bool
SaveOrTestParamHolder
l51
x1,x92
xO3){if(c8
size()<=x1){c8
xQ3
x1+1);c8
eV3
x1);c8
push_back(xO3)nX2
if(!cB
nJ2)){cB=xO3
nX2
return
xO3
xF
cB)nD3
SaveMatchedParamIndex
l41){tK.push_back(index);}
x92
GetParamHolderValueIfFound
l51
x1)const{static
const
CodeTree
xK
dummytree;if(c8
size()<=x1)return
dummytree
iF
cB;}
x92
GetParamHolderValue
l51
x1
yU1
cB;}
bool
HasRestHolder
l51
n1
yU1
lQ.size()>n1&&lQ[n1
t43==true;}
iM2&GetRestHolderValues
l51
n1)const{static
iM2
empty_result;cW1
return
empty_result
iF
lQ[n1
t53;}
const
std::vector<unsigned>&GetMatchedParamIndexes(yU1
tK
nD3
swap(t9
b){lQ.swap(b.lQ);c8
swap(b.paramholder_matches);tK.swap(b.tK);}
t9
eH1=cR1
t9
b){lQ=b.lQ;paramholder_matches=b.paramholder_matches;tK=b.tK
iF*this;}
}
;class
cZ;typedef
xT<cZ>cT;class
cZ{eQ3
int
xJ3;eQ3
cZ():xJ3(0){}
virtual~cZ(){}
}
;cW2
n51{bool
found;cT
specs;n51(bool
f):found(f),specs(){}
n51(bool
f
eR2
s):found(f),specs(s){}
}
xP3
SynthesizeRule
cR1
eS2
lR1
tree,t9
info)n11
n51
TestParam
cR1
i63,x92
tree
eR2
start_at,t9
info)n11
n51
TestParams(tA&nN,x92
tree
eR2
start_at,t9
info,bool
lG2
n11
bool
ApplyGrammar
cR1
Grammar&xX3,FPoptimizer_CodeTree::lR1
tree,bool
from_logical_context=false)xP3
ApplyGrammars(FPoptimizer_CodeTree::eK
n11
bool
IsLogisticallyPlausibleParamsMatch(tA&c32,const
eK;}
lL3
l61{eV1
DumpMatch
cR1
eS2
nO,const
FPoptimizer_Optimize::t9
info,bool
DidMatch,std::ostream&o=std::cout)xP3
DumpMatch
cR1
eS2
nO,const
FPoptimizer_Optimize::t9
info,eT2
tZ3,std::ostream&o=std::cout);}
#endif
#include <string>
eU2
l61::nM2
yN1=false);eU2
nW2
yN1=false);
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
using
lL3
l61;using
lL3
FUNCTIONPARSERTYPES;eU2
l61::nM2
yN1){
#if 1
eT2
p=0;switch(opcode
xI3
i93
p="NumConstant"
;lC
ParamHolder:p="ParamHolder"
;lC
SubFunction:p="SubFunction"
iY2}
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().size()<12)tmp<<' '
iF
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().size()<5)tmp<<' '
iF
tmp.str();
#endif
}
eU2
nW2
yN1){
#if 1
eT2
p=0;switch(opcode
xI3
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
i81:p="cNEqual"
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
iY2
#ifdef FP_SUPPORT_OPTIMIZER
case
cFetch:p="cFetch"
;lC
cPopNMov:p="cPopNMov"
;lC
eB3:p="cLog2by"
;lC
cNop:p="cNop"
iY2
#endif
case
cSinCos:p="cSinCos"
;lC
cSinhCosh:p="cSinhCosh"
;lC
cV3:p="cAbsNot"
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
l53:p="VarBegin"
iY2}
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str().size()<12)tmp<<' '
iF
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str().size()<5)tmp<<' '
iF
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
;lL3
xK1{tJ1
class
nL2{eQ3
nL2():ByteCode(),Immed(),cE(),xU(0),StackMax(0){i4
xQ3
64);Immed.xQ3
8);cE.xQ3
16)nD3
Pull(std::vector<unsigned>&bc,std::vector
xK&imm,size_t&StackTop_max){for
l51
a=0;a<xF2;++a){l82]&=~nA2;}
i4
swap(bc);Immed.swap(imm);StackTop_max=StackMax;}
size_t
GetByteCodeSize(yU1
xF2;}
size_t
GetStackTop(yU1
xU
nD3
PushVar
l51
varno){xG2
varno);eV2}
void
PushImmed(cV1
immed
nL
xG2
cImmed);Immed.push_back(immed);eV2}
void
StackTopIs(nO,int
offset=0){if((int)xU>offset){cE
l92
first=true;cE
l92
second=tree;}
}
bool
IsStackTop(nO,int
offset=0
yU1(int)xU>offset&&cE
l92
first&&cE
l92
second
xF
tree);nT1
void
EatNParams
l51
eat_count){xU-=eat_count
nD3
ProducedNParams
l51
produce_count){xL1
xU+produce_count)nD3
DoPopNMov
tS1
e32,size_t
srcpos
nL
xG2
cPopNMov)eY3
e32)eY3
srcpos);xL1
srcpos+1);cE[e32]=cE[srcpos];xL1
e32+1)nD3
DoDup
tS1
xR3
nL
if(xR3==xU-1){xG2
cDup);}
else{xG2
cFetch)eY3
xR3);}
eV2
cE[xU-1]=cE[xR3];}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tN1
int>void
Dump(){std::ostream&o=std::cout;o<<"Stack state now("
<<xU<<"):\n"
xF1
0;a<xU;++a){o<<a<<": "
;if(cE[a
t43){nO=cE[a
t53;o<<'['<<std::hex<<(void*)(&tree.l22))<<std::dec<<','<<tree.GetRefCount()<<']'
yH2
tree,o);}
else
o<<"?"
;o<<"\n"
;}
o<<std::flush;}
#endif
size_t
eZ3(nO)const{for
tS1
a=xU;a-->0;)if(cE[a
t43&&cE[a
t53
xF
tree
tB2
a
iF
eW2;}
bool
Find(nO
yU1
eZ3
eI2)!=eW2;}
bool
FindAndDup(nO){size_t
pos=eZ3
eI2
yU3
pos!=eW2){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<iY3"duplicate at ["
<<pos<<"]: "
yH2
tree
x62" -- issuing cDup or cFetch\n"
;
#endif
DoDup(pos)nX2
xB2
cW2
IfData{size_t
ofs;}
;void
SynthIfStep1
xH2,nW2
op
x51
xN2=xF2;xG2
op
i5
xG2
nA2)nD3
SynthIfStep2
xH2
x51
iV
l62+2);iV
2
l72
xN2=xF2;xG2
cJump
i5
xG2
nA2)nD3
SynthIfStep3
xH2
x51
i4
back()|=nA2;iV
l62-1);iV
2
l72
xL1
xU+1)xF1
0;a<xN2;++a){if(l82]==cJump&&l82+1]==(nA2|(xN2-1))){l82+l62-1);l82+2
l72}
switch(l82]xI3
cAbsIf:case
cIf:case
cJump:case
cPopNMov:a+=2;lC
cFCall:case
cPCall:case
cFetch:a+=1
iY2
yM3
yN3}
}
protected:void
xL1
size_t
value){xU=value;if(xU>lN3{StackMax=xU;cE.eV3
lN3;}
}
protected:std::vector<unsigned>ByteCode;std::vector
xK
Immed;std::vector<std::pair<bool,FPoptimizer_CodeTree::CodeTree
xK> >cE;size_t
xU;size_t
StackMax;private:void
incStackPtr(){if(xU+2>lN3
cE.eV3
StackMax=xU+2);}
tN1
bool
IsIntType,bool
IsComplexType>cW2
c42{}
;eQ3
void
AddOperation
eX2
e02
eat_count
e02
produce_count=1){EatNParams(eat_count);nU1(opcode);ProducedNParams(produce_count)nD3
nU1
eX2,c42<false,false>yR
false,true>yR
true,false>yR
true,true>);inline
void
nU1
eX2){nU1(opcode,c42<bool(nA
IsIntType
xK::xD3),bool(nA
IsComplexType
xK::xD3)>());}
}
n11
cW2
SequenceOpCode
n11
cW2
i01{static
cC
AddSequence;static
cC
MulSequence;}
xP3
x21
long
count,cC&eJ
xI1;}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lL3
FUNCTIONPARSERTYPES;lL3
xK1{tJ1
cW2
SequenceOpCode{cV1
basevalue
i23
op_flip
i23
op_normal,op_normal_flip
i23
op_inverse,op_inverse_flip;}
n11
cC
i01
xK::AddSequence={y11,cNeg
xD
cAdd,cSub,cRSub}
n11
cC
i01
xK::MulSequence={cV1(1),cInv,cMul,cMul,cDiv,cRDiv}
;
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode ByteCode
#define mImmed Immed
nB2
false,false>){xM1
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nB2
true,false>){xM1
# define FP_FLOAT_VERSION 0
# define FP_COMPLEX_VERSION 0
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nB2
false,true>){xM1
# define FP_FLOAT_VERSION 1
# define FP_COMPLEX_VERSION 1
# include "extrasrc/fp_opcode_add.inc"
# undef FP_COMPLEX_VERSION
# undef FP_FLOAT_VERSION
}
nB2
true,true>){xM1
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
lL3
xK1;
#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
lL3
xK1{
#ifndef FP_GENERATING_POWI_TABLE
extern
const
unsigned
char
powi_table[POWI_TABLE_SIZE];const
#endif
unsigned
char
powi_table[POWI_TABLE_SIZE]={0,1,1,1,2,1,2,1,xS3
4,1,2,xT3
2,1,xS3
8,eI3
xV3
15,1,16,1,2,1,4,1,2,xT3
2,1,4,eI3
1,16,1,25,xV3
27,5,8,3,2,1,30,1,31,3,32,1,2,1,xS3
8,1,2,xV3
39,1,16,137,2,1,4,eI3
xT3
45,135,4,31,2,5,32,1,2,131,50,1,51,1,8,3,2,1,54,1,55,3,16,1,57,133,4,137,2,135,60,1,61,3,62,133,63,1,iE1
131,iE1
139,lA2
e0
30,1,130,137,2,31,lA2
e0
e0
130,eI3
1,e0
e0
2,1,130,133,iE1
61,130,133,62,139,130,137,e0
lA2
e0
e0
iE1
131,e0
e0
130,131,2,133,lA2
130,141,e0
130,eI3
1,e0
5,135,e0
lA2
e0
lA2
130,133,130,141,130,131,e0
e0
2,131}
;}
static
x53
yX=256;
#define FPO(x)
lL3{class
PowiCache{private:int
i7[yX];int
iF1[yX];eQ3
PowiCache():i7(),iF1(){i7[1]=1;}
bool
Plan_Add(yY1,int
count){cA1>=yX
tX1;iF1[eY2+=count
iF
i7[eY2!=0
nD3
lO3
yY1){cA1<yX)i7[eY2=1
nD3
Start
tS1
value1_pos){for(int
n=2;n<yX;++n)i7[n]=-1;Remember(1,value1_pos);DumpContents();}
int
Find(yY1)const{cA1<yX){if(i7[eY2>=0){FPO(iJ3(iN3,"* I found %ld from cache (%u,%d)\n",value,(unsigned)cache[value],iK3 value]))iF
i7[eY2;}
}
return-1
nD3
Remember(yY1,size_t
iW3){cA1>=yX)return;FPO(iJ3(iN3,"* Remembering that %ld can be found at %u (%d uses remain)\n",value,(unsigned)iW3,iK3 value]));i7[eY2=(int)iW3
nD3
DumpContents()const{FPO(for(int a=1;a<POWI_CACHE_SIZE;++a)if(cache[a]>=0||iK3 a]>0){iJ3(iN3,"== cache: sp=%d, val=%d, needs=%d\n",cache[a],a,iK3 a]);})}
int
UseGetNeeded(yY1){cA1>=0&&value<yX)return--iF1[eY2
iF
0;}
}
n11
size_t
y9
long
count
iK1
cC&eJ
xI1
xP3
yZ1
size_t
apos,long
aval,size_t
bpos,long
bval
iK1
unsigned
cumulation_opcode
e02
cimulation_opcode_flip
xI1;void
l71
yY1
iK1
int
need_count,int
iZ=0){cA1<1)return;
#ifdef FP_GENERATING_POWI_TABLE
if(iZ>32)throw
false;
#endif
if(i7.Plan_Add(value,need_count
tB2;long
xW3
1;cA1<POWI_TABLE_SIZE){xW3
powi_table[eY2
nE3&128){half&=127
nE3&64)xW3-tE2
FPO(iJ3(iN3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,value/half));l71
half
xU3
i7.lO3
half)iF;}
iH1
half&64){xW3-tE2}
}
else
cA1&1)xW3
value&((1<<POWI_WINDOW_SIZE)-1);else
xW3
value/2;long
cF=value-half
nE3>cF||half<0)std::swap(half,cF);FPO(iJ3(iN3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,otherhalf))nE3==cF){l71
half,i7,2,iZ+1);}
else{l71
half
xU3
l71
cF>0?cF:-cF
xU3}
i7.lO3
value);}
tJ1
size_t
y9
yY1
iK1
cC&eJ
xI1{int
xY3=i7.Find(value
yU3
xY3>=0)yJ
xY3;}
long
xW3
1;cA1<POWI_TABLE_SIZE){xW3
powi_table[eY2
nE3&128){half&=127
nE3&64)xW3-tE2
FPO(iJ3(iN3,"* I want %ld, my plan is %ld * %ld\n",value,half,value/half));size_t
xO2=y9
half
eS1
if(i7
lB2
half)>0||xO2!=eT1){iG1
xO2)xR2
half,eT1);}
x21
value/half
nF3
size_t
iW3=eT1
xR2
value,iW3);i7.DumpContents()iF
iW3;}
iH1
half&64){xW3-tE2}
}
else
cA1&1)xW3
value&((1<<POWI_WINDOW_SIZE)-1);else
xW3
value/2;long
cF=value-half
nE3>cF||half<0)std::swap(half,cF);FPO(iJ3(iN3,"* I want %ld, my plan is %ld + %ld\n",value,half,value-half))nE3==cF){size_t
xO2=y9
half
eS1
yZ1
xO2,half,xO2,half,i7,eJ.op_normal,eJ.op_normal_flip,synth);}
else{long
part1=half;long
part2=cF>0?cF:-cF;size_t
part1_pos=y9
part1
eS1
size_t
part2_pos=y9
part2
eS1
FPO(iJ3(iN3,"Subdivide(%ld: %ld, %ld)\n",value,half,otherhalf));yZ1
part1_pos,part1,part2_pos,part2,i7,cF>0?eJ.op_normal:eJ.op_inverse,cF>0?eJ.op_normal_flip:eJ.op_inverse_flip,synth);}
size_t
iW3=eT1
xR2
value,iW3);i7.DumpContents()iF
iW3;}
eV1
yZ1
size_t
apos,long
aval,size_t
bpos,long
bval
iK1
unsigned
cumulation_opcode
e02
cumulation_opcode_flip
xI1{int
a_needed=i7
lB2
aval);int
xZ3=i7
lB2
bval);bool
lC2=false;
#define DUP_BOTH() do{if(apos<bpos){size_t tmp=apos;apos=bpos;bpos=tmp;lC2=!lC2;}FPO(iJ3(iN3,"-> " iU3 iU3"op\n",(unsigned)apos,(unsigned)bpos));iG1 apos);iG1 apos==bpos?eT1:bpos);}while(0)
#define DUP_ONE(p) do{FPO(iJ3(iN3,"-> " iU3"op\n",(unsigned)p));iG1 p);}while(0)
if(a_needed>0){if(xZ3>0){nC2}
e43
bpos!=eT1)nC2
else{lD2
lC2=!lC2;}
}
}
iH1
xZ3>0){if(apos!=eT1)nC2
else
DUP_ONE(bpos);}
e43
apos==bpos&&apos==eT1)lD2
iH1
apos==eT1&&bpos==e73
xI
2){FPO(iJ3(iN3,"-> op\n"));lC2=!lC2;}
iH1
apos==e73
xI
2&&bpos==eT1)FPO(iJ3(iN3,"-> op\n"));iH1
apos==eT1)DUP_ONE(bpos);iH1
bpos==eT1){lD2
lC2=!lC2;}
else
nC2}
nG3
lC2?cumulation_opcode_flip:cumulation_opcode,2);}
eV1
cB1
long
count,cC&eJ
xI1{while
tS3<256){int
xW3
xK1::powi_table[count]nE3&128){half&=127;cB1
half
nF3
count/=half;}
else
yN3
if
tS3==1)return;if(!tS3&1)){nG3
cSqr,1);cB1
count/2
nF3}
else{iG1
eT1);cB1
count-1
nF3
nG3
cMul,2);}
}
}
lL3
xK1{eV1
x21
long
count,cC&eJ
xI1{if
tS3==0)x82
eJ.basevalue);else{bool
eZ2=false;if
tS3<0){eZ2=true;count=-count;}
if(false)cB1
count
nF3
iH1
count>1){PowiCache
i7;l71
count,i7,1);size_t
xN1=e73
GetStackTop();i7.Start(eT1);FPO(iJ3(iN3,"Calculating result for %ld...\n",count));size_t
xP2=y9
count
eS1
size_t
n_excess=e73
xI
xN1;if(n_excess>0||xP2!=xN1-1){e73
DoPopNMov(xN1-1,xP2);}
}
if(eZ2)nG3
eJ.op_flip,1);}
}
}
#endif
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH
t1{lL3
lP3{iS
cW2
Comp{}
;tN1>cW2
Comp<nA
cLess>{tN1
lE<n9
cLessOrEq>{tN1
lE<=n9
cGreater>{tN1
lE>n9
cGreaterOrEq>{tN1
lE>=n9
cEqual>{tN1
lE==n9
i81>{tN1
lE!=b;}
}
;}
tJ1
cW2
cC1{cV1
val;bool
known;cC1():val(),known(false){}
cC1
cQ1
v):val(v),known(true){nT1
void
set
cQ1
v){known=true;val=v
nD3
set(cV1(iI1(cV1),cC1
xJ)val=iN2
void
set(cV1(iI1
cQ1),cC1
xJ)val=iN2
iS
void
set_if(cV1
v,cV1(iI1(cV1),cC1
xJ&&lP3::Comp<Compare>()(val,v))val=iN2
iS
void
set_if
cQ1
v,cV1(iI1
cQ1),cC1
xJ&&lP3::Comp<Compare>()(val,v))val=iN2}
n11
cW2
range{cC1
xK
min,max;range():min(),max(){}
range(cV1
mi,cV1
ma):min(mi),max(ma){}
range(bool,cV1
ma):min(),max(ma){}
range(cV1
mi,bool):min(mi),max(){}
void
set_abs();void
set_neg();}
n11
bool
IsLogicalTrueValue
cR1
yG3&p
lI3
n11
bool
IsLogicalFalseValue
cR1
yG3&p
lI3;}
#endif
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH
t1{enum
TriTruthValue{IsAlways,eO3,Unknown}
n11
yG3
iT
const
eK
n11
bool
IsLogicalValue
cR1
eK
n11
TriTruthValue
GetIntegerInfo
cR1
eK
n11
xO1
GetEvennessInfo
cR1
eK{if(!tree
tF1)return
Unknown;yJ1=tree
l31;if(nA
isEvenInteger(value
tC2
nA
isOddInteger(value
tD2
tJ1
xO1
GetPositivityInfo
cR1
eK{yG3
p=iT
tree
yU3
p
iN&&p.yQ>=cV1(tC2
p
cO
known
lH1
tD2
tJ1
xO1
GetLogicalValue
lQ3
tree
lI3{yG3
p=iT
tree
yU3
IsLogicalTrueValue(p,abs
tC2
IsLogicalFalseValue(p,abs
tD2}
#endif
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH
t1{eV1
ConstantFolding(eK;}
#endif
lL3{using
lL3
FUNCTIONPARSERTYPES;using
t1;cW2
ComparisonSetBase{enum{t03=0x1,Eq_Mask=0x2,Le_Mask=0x3,t13=0x4,t23=0x5,Ge_Mask=0x6}
;static
int
Swap_Mask(int
m)yJ(m&Eq_Mask)|((m&t03)?t13:0)|((m&t13)?t03:0);}
enum
c01{Ok,BecomeZero,BecomeOne,nB1
x73
nN2{cond_or,iO2,iP2,iQ2}
;}
n11
cW2
ComparisonSet:public
ComparisonSetBase{cW2
t02{CodeTree
xK
a
iG
b;int
relationship;t02():a(),b(),relationship(){}
}
;std::vector<t02>tG;cW2
Item{CodeTree
xK
value;bool
c52;Item():value(),c52(false){}
}
;std::vector<Item>cD1;int
xP1;ComparisonSet():tG(),cD1(),xP1(0){}
c01
AddItem
lQ3
a,bool
c52,nN2
type){for
tS1
c=0;c<cD1
iU2
c)if(cD1[c].value
xF
a)){if(c52!=cD1[c].c52)iI
cX1
case
iQ2:cD1.erase(cD1.begin()+c);xP1
c62
case
iO2:case
iP2:cY1}
}
return
nB1;}
Item
pole;pole.value=a;pole.c52=c52;cD1.push_back(pole)iF
Ok;}
c01
AddRelationship(CodeTree
xK
a,CodeTree
xK
b,int
i11,nN2
type)iI
if(i11==7)cX1
lC
iQ2:if(i11==7){xP1
c62}
lC
iO2:case
iP2:if(i11==0)cY1
yN3
if(!(a.GetHash()<b.GetHash())){a.swap(b);i11=Swap_Mask(i11);}
for
tS1
c=0;c<tG
iU2
c){if(tG[c].a
xF
a)&&tG[c].b
xF
b))iI{int
y03=tG[c
e9|i11;if(y03==7)cX1
tG[c
e9=y03
iY2}
case
iO2:case
iP2:{int
y03=tG[c
e9&i11;if(y03==0)cY1
tG[c
e9=y03
iY2}
case
iQ2:{int
newrel_or=tG[c
e9|i11;int
xS2=tG[c
e9&i11;lH2
5&&xS2==0){tG[c
e9=t23
iF
nB1;}
lH2
7&&xS2==0){xP1+=1;tG.erase(tG.begin()+c)iF
nB1;}
lH2
7&&xS2==Eq_Mask){tG[c
e9=Eq_Mask;xP1
c62}
y91}
}
return
nB1;}
}
t02
comp;comp.a=a;comp.b=b;comp.relationship=i11;tG.push_back(comp)iF
Ok;}
}
;nR1
cV1,xG3
CondType>bool
ConstantFolding_LogicCommon(lR1
tree,CondType
xE1,bool
xT2){bool
should_regenerate=false;ComparisonSet
xK
comp
cZ3
xG3
c4
c01
eJ3=c4
Ok;x92
atree=t73;switch(atree
nC
xI3
cEqual
lG
Eq_Mask,xE1);lC
i81
lG
t23,xE1);lC
cLess
lG
t03,xE1);lC
cLessOrEq
lG
Le_Mask,xE1);lC
cGreater
lG
t13,xE1);lC
cGreaterOrEq
lG
Ge_Mask,xE1);lC
cNot:eJ3
cS1
l8
0),true,xE1);lC
cNotNot:eJ3
cS1
l8
0),false,xE1)iY2
yM3
if(xT2||IsLogicalValue(atree))eJ3
cS1,false,xE1);}
switch(eJ3){ReplaceTreeWithZero:xO
0)iF
true;ReplaceTreeWithOne:xO
1);x0
c4
Ok:lC
c4
BecomeZero
t3
c4
BecomeOne:i2
c4
nB1:c41
yN3}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_LogicCommon: "
cN
#endif
if(xT2){eP1;}
else{for
cK{x92
atree=t5
a
yU3
IsLogicalValue(atree))nN1);}
}
for
tS1
a=0;a<comp.cD1
iU2
a){if(comp.cD1[a].c52){c11
cNot);r.c3
r
t33
iH1!xT2){c11
cNotNot);r.c3
r
t33
else
tree.c3}
for
tS1
a=0;a<comp.tG
iU2
a){c11
cNop);switch(comp.tG[a
e9
xI3
c4
t03:r
xC
cLess);lC
c4
Eq_Mask:r
xC
cEqual);lC
c4
t13:r
xC
cGreater);lC
c4
Le_Mask:r
xC
cLessOrEq);lC
c4
t23:r
xC
i81);lC
c4
Ge_Mask:r
xC
cGreaterOrEq
iX2
r
c51
comp.tG[a].a);r
c51
comp.tG[a].b);r
t33
if(comp.xP1!=0)tree.yB
cV1(comp.xP1)));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_LogicCommon: "
cN
#endif
return
true;}
xB2
xQ1
ConstantFolding_AndLogic(iS3(tree.GetOpcode()==cAnd
lA4()==cAbsAnd)iF
nI
iO2,true);}
xQ1
ConstantFolding_OrLogic(iS3(tree.GetOpcode()==cOr
lA4()==cAbsOr)iF
nI
cond_or,true);}
xQ1
ConstantFolding_AddLogicItems(iS3(tree.GetOpcode()==cAdd)iF
nI
iQ2,false);}
xQ1
ConstantFolding_MulLogicItems(iS3(tree.GetOpcode()==cMul)iF
nI
iP2,false);}
}
#include <vector>
#include <map>
#include <algorithm>
lL3{using
lL3
FUNCTIONPARSERTYPES;using
t1;cW2
CollectionSetBase{enum
xR1{Ok,nB1}
;}
n11
cW2
CollectionSet:public
CollectionSetBase{cW2
c21{CodeTree
xK
value
iG
xU2;bool
cU;c21():value(),xU2(),cU(false){}
c21
lQ3
v,x92
f):value(v),xU2(f),cU(false){}
}
;std::multimap<fphash_t,c21>i8;typedef
xG3
std::multimap<fphash_t,c21>::y53
xS1;CollectionSet():i8(){}
xS1
FindIdenticalValueTo
lQ3
value){fphash_t
hash=value.GetHash();for(xS1
i=i8.xV2
hash);i!=i8.cZ1
hash;++i){cA1
xF
i
e42.value
tB2
i;}
return
i8
l24;}
bool
Found
cR1
xS1&b)yJ
b!=i8
l24;}
xR1
AddCollectionTo
lQ3
xU2,const
xS1&into_which){c21&c=into_which
e42;if(c.cU)c.xU2
tM
xU2);else{CodeTree
xK
add;add
xC
cAdd);add
c51
c.xU2);add
tM
xU2);c.xU2.swap(add);c.cU=true;}
return
nB1;}
xR1
nO2
lQ3
value,x92
xU2){const
fphash_t
hash=value.GetHash();xS1
i=i8.xV2
hash);for(;i!=i8.cZ1
hash;++i){if(i
e42.value
xF
value
tB2
AddCollectionTo(xU2,i);}
i8.y13,std::make_pair(hash,c21(value,xU2)))iF
Ok;}
xR1
nO2
lQ3
a)yJ
nO2(a,nC1
1)));}
}
n11
cW2
ConstantExponentCollection{typedef
eH
y33;typedef
std::x61
xW2;std::vector<xW2>data;ConstantExponentCollection():data(){}
void
MoveToSet_Unique
cQ1
eI1&eJ1){data.push_back(std::x61(eI1()));data.back()t63.swap(eJ1)nD3
MoveToSet_NonUnique
cQ1
eI1&eJ1){xG3
std::vector<xW2>::y53
i=std::xV2
data.iR2
data
l24,eX3,Compare1st()yU3
i!=data.cZ1
eX3){i
e42.y13
e42
l24,eJ1.iR2
eJ1
l24);}
else{data.y13,std::x61
eW3,eJ1));}
}
bool
iJ2{bool
changed=false;std::sort(data.iR2
data
l24,Compare1st());redo:for
tS1
a=0;a<data
iU2
a){cV1
exp_a=data[a
t43;if(fp_equal(exp_a,y63
y91
for
tS1
b=a+1;b<data
iU2
b){cV1
exp_b=data[b
t43;cV1
xX2=exp_b-exp_a;if(xX2>=fp_abs(exp_a))break;cV1
exp_diff_still_probable_integer=xX2*cV1(16
yU3
t12
exp_diff_still_probable_integer)&&!(t12
exp_b)&&!t12
xX2))){y33&a_set=lI2;y33&b_set=data[b
t53;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantExponentCollection iteration:\n"
;t22
cout);
#endif
if(isEvenInteger(exp_b)&&!isEvenInteger(xX2+exp_a
lC1
tmp2;tmp2
tM3
tmp2.SetParamsMove(b_set);tmp2.t32
tmp;tmp
xC
cAbs);tmp
c51
tmp2);tmp.lH3
b_set.eV3
1);b_set[0].iA2}
a_set.insert(a_set
l24,b_set.iR2
b_set
l24);y33
b_copy=b_set;data.erase(data.begin()+b);MoveToSet_NonUnique(xX2,b_copy);yM1
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantExponentCollection iteration:\n"
;t22
cout);
#endif
i3}
}
}
return
changed;}
#ifdef DEBUG_SUBSTITUTIONS
void
t22
ostream&out){for
tS1
a=0;a<data
iU2
a){out.precision(12);out<<data[a
t43<<": "
;tR1
lI2
iU2
b){if(b>0)out<<'*'
yH2
lI2[b],out);}
out<<std::endl;}
}
#endif
}
n11
static
CodeTree
xK
x71
lR1
value,bool&xP){switch(value
nC
xI3
cPow:{CodeTree
xK
e52
value
l8
1);value.y81
iF
eX3
lD3
cRSqrt:value.y81;xP=true
iF
nC1-0.5))tJ3
cInv:value.y81;xP=true
iF
nC1-1));yM3
yN3
return
nC1
1));}
c81
void
eK1
eL1&mul,x92
tree,x92
xU2,bool&c31
bool&xP){for
tS1
a=tL{CodeTree
xK
value(t73)iG
eX3(x71
value,xP)yU3!xU2
tF1||xU2
l31!=cV1(1.0
lC1
e41;e41
tM3
e41
tM
i72
e41
tM
xU2);e41.lH3
eX3.swap(e41);}
#if 0 /* FIXME: This does not work */
cA1
nC==cMul){if(1){bool
exponent_is_even=eX3
tF1&&isEvenInteger
eW3
l31);tR1
value.tA3{bool
tmp=false
iG
val(value
l8
b))iG
exp(x71
val,tmp)yU3
exponent_is_even||(exp
tF1&&isEvenInteger(exp
l31)lC1
e41;e41
tM3
e41
tM
i72
e41
c51
exp);e41.ConstantFolding(yU3!e41
tF1||!isEvenInteger(e41
l31)){goto
cannot_adopt_mul;}
}
}
}
eK1
mul,value,eX3,c31
xP);}
else
cannot_adopt_mul:
#endif
{if(mul.nO2(value,eX3)==CollectionSetBase::nB1)c41}
}
}
xQ1
ConstantFolding_MulGrouping(eK{bool
xP=false;bool
should_regenerate=false;eL1
mul;eK1
mul,tree,nC1
1)),c31
xP);typedef
std::pair<CodeTree
xK,eH>eM1;typedef
std::multimap<fphash_t,eM1>cE1;cE1
tT;iS2
eL1::xS1
j=mul.i8.y43
j!=mul.i8
l24;++j){lR1
value=j
e42.value
iG&e52
j
e42.xU2;if(j
e42.cU)eX3.lH3
const
fphash_t
eN1=eX3.GetHash();xG3
cE1::y53
i=tT.xV2
eN1);for(;i!=tT.cZ1
eN1;++i)if(i
e42.first
xF
eX3)){if(!eX3
tF1||!e51
l31,y63
c41
i
e42
t63.push_back(value);goto
skip_b;}
tT.y13,std::make_pair(eN1,std::make_pair
eW3,eH
tS1(1),value))));skip_b:;}
#ifdef FP_MUL_COMBINE_EXPONENTS
ConstantExponentCollection
xK
e61;iS2
cE1::y53
j,i=tT.y43
i!=tT
l24;i=j){j=i;++j;eM1&list=i
e42;yV2
lE1
e52
list.first
l31;if(!eW3==xH1)e61.MoveToSet_Unique
eW3,list
yA3
tT.erase(i);}
}
if(e61.iJ2)c41
#endif
if(should_regenerate){CodeTree
xK
before=tree;before.l81
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_MulGrouping: "
yH2
before
x62"\n"
;
#endif
eP1;iS2
cE1::y53
i=tT.y43
i!=tT
l24;++i){eM1&list=i
e42;
#ifndef FP_MUL_COMBINE_EXPONENTS
yV2
lE1
e52
list.first
l31;if
eW3==xH1
y91
if(e51,y63{tree.AddParamsMove(list
yA3
y91}
}
#endif
CodeTree
xK
mul;mul
tM3
mul.SetParamsMove(list
yA3
mul.lH3
if(xP&&list.first
tF1){yV2
l31==cV1(1)/cV1(3
lC1
cbrt;cbrt
xC
cCbrt);cbrt.e1
cbrt
yM
cbrt);y91
eO
0.5
lC1
sqrt;sqrt
xC
cSqrt);sqrt.e1
sqrt
yM
sqrt);y91
eO-0.5
lC1
rsqrt;rsqrt
xC
cRSqrt);rsqrt.e1
rsqrt
yM
rsqrt);y91
eO-1
lC1
inv;inv
xC
cInv);inv.e1
inv
yM
inv);y91}
}
CodeTree
xK
pow
tL1.e1
pow
c51
list.first);pow
yM
pow);}
#ifdef FP_MUL_COMBINE_EXPONENTS
tT.clear()xF1
0;a<tZ
iU2
a){cV1
e52
tZ[a
t43;if(e51,y63{tree.AddParamsMove(tZ[a]yA3
y91}
CodeTree
xK
mul;mul
tM3
mul.SetParamsMove(tZ[a]yA3
mul.t32
pow
tL1.e1
pow.yB
eX3));pow
yM
pow);}
#endif
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_MulGrouping: "
cN
#endif
return!tree
xF
before);}
xB2
xQ1
ConstantFolding_AddGrouping(eK{bool
should_regenerate=false;eL1
add
cZ3
if(t73
nC==cMul)y91
if(add.nO2(t73)==CollectionSetBase::nB1)c41}
cE3
remaining
eI2.lT1;size_t
tB=0
cZ3
x92
xF3=t73;if
iK2
nC==cMul){tR1
c61
tA3{if
iK2
l8
b)tF1)y91
xG3
eL1::xS1
c=add.FindIdenticalValueTo
iK2
l8
b)yU3
add.Found(c
lC1
tmp
iK2
yA
CloneTag());tmp
x02
b);tmp.lH3
add.AddCollectionTo(tmp,c);c41
goto
done_a;}
}
remaining[a]=true;tB+=1;done_a:;}
}
if(tB>0){if(tB>1){std::vector<std::pair<CodeTree
xK,size_t> >nZ;std::multimap<fphash_t,size_t>eO1;bool
lS3=false
xF1
tL
tU1{tR1
t73.tA3{x92
p=t73
l8
b);const
fphash_t
p_hash=p.GetHash();for(std::multimap<fphash_t,size_t>::const_iterator
i=eO1.xV2
p_hash);i!=eO1.cZ1
p_hash;++i){if(nZ[i
e42
t43
xF
p)){nZ[i
e42
t53+=1;lS3=true;goto
found_mulgroup_item_dup;}
}
nZ.push_back(std::make_pair(p,size_t(1)));eO1.insert(std::make_pair(p_hash,nZ.size()-1));found_mulgroup_item_dup:;}
}
if(lS3){CodeTree
xK
e62;{size_t
max=0;for
tS1
p=0;p<nZ
iU2
p)if(nZ[p
t53<=1)nZ[p
t53=0;else{nZ[p
t53*=nZ[p
t43
nD2;if(nZ[p
t53>max){e62=nZ[p
t43;max=nZ[p
t53;}
}
}
CodeTree
xK
group_add;group_add
xC
cAdd);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Duplicate across some trees: "
yH2
e62
x62" in "
cN
#endif
for
tS1
a=tL
tU1
tR1
t73.tA3
if(e62
xF
t73
l8
b)lC1
tmp(t73
yA
CloneTag());tmp
x02
b);tmp.lH3
group_add
c51
tmp);remaining[a]=false
iY2}
group_add.t32
group;group
tM3
group
c51
e62);group
c51
group_add);group.lH3
add.nO2(group);c41}
}
for
tS1
a=tL
tU1{if(add.nO2(t73)==CollectionSetBase::nB1)c41}
}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_AddGrouping: "
cN
#endif
eP1;iS2
eL1::xS1
j=add.i8.y43
j!=add.i8
l24;++j){lR1
value=j
e42.value
iG&coeff=j
e42.xU2;if(j
e42.cU)coeff.lH3
if(coeff
eC3
fp_equal(coeff
l31,xH1)y91
if(fp_equal(coeff
l31,y63{tree
c51
value);y91}
}
CodeTree
xK
mul;mul
tM3
mul
c51
value);mul
c51
coeff);mul
yM
mul);}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_AddGrouping: "
cN
#endif
return
true;}
xB2}
lL3{using
lL3
FUNCTIONPARSERTYPES;using
t1
n11
bool
ConstantFolding_IfOperations(iS3(tree.GetOpcode()==cIf
lA4()==cAbsIf);for(;;){if(nM3
cNot){t82
cIf);t5
0).e72
0)cO3
t5
1).swap(t5
2));}
iH1
t5
0)cP1{t82
tB3;t5
0).e72
0)cO3
t5
1).swap(t5
2));}
else
yN3
eK2
0),t52==tB3)i91
tree.e72
1));x0
lE3
tree.e72
2));x0
n01
if(nM3
cIf||nM3
tB3{CodeTree
xK
cond=t5
0)iG
lT3;lT3
tF2==cIf?cNotNot:cAbsNotNot);lT3
xY2
1));ConstantFolding(lT3)iG
lU3;lU3
tF2==cIf?cNotNot:cAbsNotNot);lU3
xY2
2));ConstantFolding(lU3
yU3
lT3
tF1||lU3.IsImmed(lC1
eR;eR
tF2);eR
xY2
1));eR
eS3
eR.nJ
2));eR.t32
eS;eS
tF2);eS
xY2
2));eS
eS3
eS.nJ
2));eS
tO
cond
nC);eQ1
0,cond
l8
0))eB
1,eR)eB
2,eS)nX2}
if(yV3
t5
2)nC&&(yV3
cIf||yV3
cAbsIf
lC1&yD3=t5
1)iG&leaf2=t5
2);if
iT2
0)xG1
0))&&iT2
1)xG1
1))||yD3
l8
2)xG1
2))lC1
eR;eR
iJ1
eR
eU3;eR
xZ2
1));eR
y02
1));eR.t32
eS;eS
iJ1
eS
eU3;eS
xZ2
2));eS
y02
2));eS
tO
yD3
nC);eQ1
0
eT3
0))eB
1,eR)eB
2,eS)nX2
if
iT2
1)xG1
1))&&yD3
l8
2)xG1
2)lC1
eT;eT
iJ1
eT
c51
t5
0));eT
xZ2
0));eT
y02
0));eT
tO
yD3
nC)eB
0,eT);eQ1
2
eT3
2));eQ1
1
eT3
1))nX2
if
iT2
1)xG1
2))&&yD3
l8
2)xG1
1)lC1
e82;e82
xC
leaf2
nC==cIf?cNot:cV3);e82
y02
0));e82.t32
eT;eT
iJ1
eT
c51
t5
0));eT
xZ2
0));eT
c51
e82);eT
tO
yD3
nC)eB
0,eT);eQ1
2
eT3
2));eQ1
1
eT3
1))nX2}
lR1
y1=t5
1)iG&y7=t5
2
yU3
y1
xF
y7)){tree.e72
1))nX2
const
OPCODE
op1=y1
nC;const
OPCODE
op2=y7
nC;if
yE3
op2){if(y1.e71
1){CodeTree
xK
lO
0))nZ2
0
eY1
n4
if(y1.e71
2&&y7.e71
2){if(y1
l8
0)xF
y7
l8
0)lC1
param0=y1
l8
0)iG
lO
1))nZ2
1
eY1
x12
param0)n4
if(y1
l8
1)xF
y7
l8
1)lC1
param1=y1
l8
1)iG
lO
0))nZ2
0
eY1
x12
nA1)x12
param1)nX2}
if
yE3
yF3
cMul
lJ2
cAnd
lJ2
cOr
lJ2
cAbsAnd
lJ2
cAbsOr
lJ2
cMin
lJ2
cMax){eH
lV3;yY{for
tS1
b=y7.l21
b-->0;){if(y1
lW3
y7
l8
b))){if(lV3.empty()){y1.l81
y7.l81}
lV3.push_back(y1
nL3
y7
x02
b);y1
x02
a
iX2}
}
if(!lV3.empty()){y1.lH3
y7.Rehash()l7
op1);tree.SetParamsMove(lV3)n4}
}
if
yE3
yF3
cMul||yE3
cAnd
nV1
y7))||yE3
cOr
nV1
y7))){yY
if(y1
lW3
y7)){y1.l81
y1
x02
a);y1.t32
cF1=y7;y7=tC
op1==yF3
cOr)iV2
op1)x12
cF1)n4}
if(yE3
cAnd
lJ2
cOr)&&op2==cNotNot){lR1
lX3=y7
l8
0);yY
if(y1
lW3
lX3)){y1.l81
y1
x02
a);y1.t32
cF1=lX3;y7=tC
op1==cOr)iV2
op1)x12
cF1)n4}
if(op2==cAdd||op2==cMul||(op2==cAnd
nV1
y1))||(op2==cOr
nV1
y1))){for
tS1
a=y7.cG1
y7
lW3
y1)){y7.l81
y7
x02
a);y7.t32
cI1=y1;y1=tC
op2==cAdd||op2==cOr)iV2
op2)x12
cI1)n4}
if((op2==cAnd||op2==cOr)&&op1==cNotNot){lR1
lY3=y1
l8
0)xF1
y7.cG1
y7
lW3
lY3)){y7.l81
y7
x02
a);y7.t32
cI1=lY3;y1=tC
op2==cOr)iV2
op2)x12
cI1)n4}
xB2}
#include <limits>
lL3{using
lL3
FUNCTIONPARSERTYPES;using
t1
n11
int
maxFPExponent()yJ
std::numeric_limits
xK::max_exponent;}
xQ1
x81
cV1
base,cV1
eX3){if(base<xH1
return
true;if(fp_equal(base,xH1||fp_equal(base,e53
tX1
iF
eX3>=cV1(maxFPExponent
xK())/fp_log2(base);}
xQ1
ConstantFolding_PowOperations(iS3(tree.GetOpcode()==cPow);nS&&t5
1)lE1
const_value=tC3
lR,eL);xO
const_value)iF
false;}
if(i62
fp_equal(eL,y63{tree.e72
0))nX2
nS&&fp_equal(lR,y63{xO
1)iF
false;}
nS&&yV3
cMul){bool
y12=false;cV1
lK2=lR
iG
xF3=t5
1)xF1
c61
cG1
xF3
l8
a)lE1
imm=xF3
l8
a)l31;{if(x81
lK2,imm))break;cV1
lL2=tC3
lK2,imm
yU3
fp_equal(lL2,xH1)break;if(!y12){y12=true;c61
l81}
lK2=lL2;c61
DelParam(a
iX2}
if(y12){c61
lH3
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before pow-mul change: "
cN
#endif
t5
0).Become(e91
lK2));t5
1).Become
iK2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After pow-mul change: "
cN
#endif
}
}
if(i62
nM3
cMul){cV1
lM2=eL;cV1
y22=1.0;bool
y12=false
iG&xF3=t5
0)xF1
c61
cG1
xF3
l8
a)lE1
imm=xF3
l8
a)l31;{if(x81
imm,lM2))break;cV1
t01=tC3
imm,lM2
yU3
fp_equal(t01,xH1)break;if(!y12){y12=true;c61
l81}
y22*=t01;c61
DelParam(a
iX2}
if(y12){c61
t32
eK3;eK3
xC
cPow);eK3.SetParamsMove
eI2.l22));eK3.lX2
t82
cMul)x12
eK3);tree
tM
e91
y22))nX2}
if(nM3
cPow&&i62
t5
0)l8
1)lE1
a=t5
0)l8
1)l31;cV1
b=eL;cV1
c=a*b;if(isEvenInteger(a)&&!isEvenInteger(c
lC1
lZ3;lZ3
xC
cAbs);lZ3.nJ
0)cO3
lZ3.Rehash()eB
0,lZ3);}
else
eQ1
0,t5
0)cO3
eQ1
1,e91
c));}
xB2}
lL3{using
lL3
FUNCTIONPARSERTYPES;using
t1;cW2
l5{enum
e92{xC2=0,MakeTrue=1,t42=2,n03=3,MakeNotNotP0=4,MakeNotNotP1=5,MakeNotP0=6,MakeNotP1=7,xA=8
x73
lN2{Never=0,Eq0=1,Eq1=2,yL3=3,yQ3=4}
;e92
if_identical;e92
lO2
4];cW2{e92
what:4;lN2
when:4;}
iM1,iN1,iO1,iP1
n11
e92
Analyze
lQ3
a,x92
b)const{if(a
xF
b
tB2
if_identical;yG3
p0=iT
a);yG3
p1=iT
b);yJ3
known&&p1
iN){yJ3
val<p1.yQ&&lO2
0]i0
0];yJ3
val<=p1.yQ&&lO2
1]i0
1];}
if(p0
iN&&p1
eU1{if(p0.yQ>p1
cO
val&&lO2
2]i0
2];if(p0.yQ>=p1
cO
val&&lO2
3]i0
3];}
if(IsLogicalValue(a)){if(iM1
i33
iM1.when,p1
tB2
iM1.what;if(iO1
i33
iO1.when,p1
tB2
iO1.what;}
if(IsLogicalValue(b)){if(iN1
i33
iN1.when,p0
tB2
iN1.what;if(iP1
i33
iP1.when,p0
tB2
iP1.what;}
return
xA;}
c81
bool
TestCase(lN2
when,const
yG3&p){if(!p
iN||!p
eU1
return
false;switch(when
xI3
Eq0:lB3.yQ==cV1(0.0)yK3==p.yQ
tJ3
Eq1:lB3.yQ==cV1(1.0)yK3==p
cO
val
tJ3
yL3:lB3.yQ>y11
yK3<=cV1(1)tJ3
yQ3:lB3.yQ>=y11
lH1
1);yM3;}
xB2}
;lL3
RangeComparisonsData{static
const
l5
Data[6]={{l5
tD3
tU
xA,l5::tU
xA
nX1
Eq1}
,nE1
Eq1}
tF3
Eq0}
tG3
Eq0}
}
tH3,{l5
tE3
xA,l5
tE3
xA
nX1
Eq0}
,nE1
Eq0}
tF3
Eq1}
tG3
Eq1}
}
tH3,{l5
tE3
t42,l5::tU
xC2}
tF3
yL3}
,nE1
yQ3
yU,{l5
tD3
xA,l5
tE3
tU
n03}
tF3
yQ3}
,nE1
yL3
yU
tH3
lP2
tU
tU
MakeTrue,l5::t42
nX1
yQ3}
tG3
yL3
yU,{l5
tD3
tU
n03,l5::xA,l5
nD1
nX1
yL3}
tG3
yQ3
yU}
;}
xQ1
ConstantFolding_Comparison(eK{using
lL3
RangeComparisonsData;assert(tree.GetOpcode()>=cEqual&&tree.GetOpcode()<=cGreaterOrEq);switch(Data[t52-cEqual].Analyze(t5
0),t5
1))xI3
l5::xC2:xO
0)e13
nD1:xO
1)e13::n03:t82
cEqual)e13::t42:t82
i81)e13::MakeNotNotP0:t82
cNotNot
t61
1)e13::MakeNotNotP1:t82
cNotNot
t61
0)e13::MakeNotP0:t82
cNot
t61
1)e13::MakeNotP1:t82
cNot
t61
0)e13::xA:;}
if(t5
1)tF1)switch(t5
0)nC
xI3
cAsin:lM
fp_sin(l23
cAcos:lM
fp_cos(eL)));t82
t52==cLess?cGreater:t52==cLessOrEq?cGreaterOrEq:t52==cGreater?cLess:t52==cGreaterOrEq?cLessOrEq:t52);x0
cAtan:lM
fp_tan(l23
cLog:lM
fp_exp(l23
cSinh:lM
fp_asinh(l23
cTanh:if(fp_less(fp_abs(eL),y63{lM
fp_atanh(eL)))nX2
break;yM3
yN3
xB2}
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lL3
FUNCTIONPARSERTYPES;lL3{
#ifdef DEBUG_SUBSTITUTIONS
yD
double
d){union{double
d;uint_least64_t
h;c72
d=d;lN1
h
nZ1
#ifdef FP_SUPPORT_FLOAT_TYPE
yD
float
f){union{float
f;uint_least32_t
h;c72
f=f;lN1
h
nZ1
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
yD
long
double
ld){union{long
double
ld;cW2{uint_least64_t
a
i23
short
b;}
s;c72
ld=ld;lN1
s.b<<data.s.a
nZ1
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
yD
long
ld){o<<"("
<<std::hex<<ld
nZ1
#endif
#endif
}
t1{lN
nE)){}
lN
const
cV1&i
yA
xM3
nE
i
nY2
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lN
cV1&&i
yA
xM3
nE
std::move(i)nY2
#endif
lN
unsigned
v
yA
VarTag
nE
l53,v
nY2
lN
nW2
o
yA
OpcodeTag
nE
o
nY2
lN
nW2
o
e02
f
yA
FuncOpcodeTag
nE
o,f
nY2
lN
x92
b
yA
CloneTag
nE*b.data)){}
tJ1
CodeTree
xK::~CodeTree(){}
lB
ReplaceWithImmed
cQ1
i){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Replacing "
yH2*this
yU3
IsImmed())OutFloatHex(std::cout,GetImmed()x62" with const value "
<<i;OutFloatHex(std::cout,i
x62"\n"
;
#endif
data=new
xL2
xK(i);}
tJ1
cW2
ParamComparer{iL2()lQ3
a,x92
b)const{if(a
nD2!=b
nD2)return
a
nD2<b
nD2
iF
a.GetHash()<b.GetHash();}
}
xP3
xL2
xK::Sort(){switch(Opcode
xI3
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
i81:std::sort(l43
iR2
l43
end(),ParamComparer
xK());lC
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
break;yM3
yN3}
lB
AddParam
lQ3
param){y3.push_back(param);}
lB
nK2
lR1
param){y3.push_back(CodeTree
xK());y3.back().swap(param);}
lB
SetParam
tS1
which,x92
b)x01
which
l33
y3[which]=b;}
lB
SetParamMove
tS1
which,lR1
b)x01
which
l33
y3[which
i13
b);}
lB
AddParams
cR1
nK){y3.insert(y3
l24,lQ2.iR2
lQ2
l24);}
lB
AddParamsMove(nK){size_t
endpos=y3.size(),added=lQ2.size();y3.eV3
endpos+added,CodeTree
xK());for
tS1
p=0;p<added;++p)y3[endpos+p
i13
lQ2[p]);}
lB
AddParamsMove(nK,size_t
lR2)x01
lR2
l33
DelParam(lR2);AddParamsMove(tM1}
lB
SetParams
cR1
nK){eH
tmp(tM1
y3.iA2}
lB
SetParamsMove(nK){y3.swap(tM1
lQ2.clear();}
#ifdef __GXX_EXPERIMENTAL_CXX0X__
lB
SetParams(eH&&lQ2){SetParamsMove(tM1}
#endif
lB
DelParam
tS1
index){eH&eH3=y3;
#ifdef __GXX_EXPERIMENTAL_CXX0X__
l43
erase(l43
begin()+index);
#else
eH3[index].data=0;for
tS1
p=index;p+1<eH3
iU2
p)eH3[p].data.UnsafeSetP(&*eH3[p+1
l33
eH3[nN3-1].data.UnsafeSetP(0);l43
eV3
nN3-1);
#endif
}
lB
DelParams(){y3.clear();}
xQ1
CodeTree
xK::IsIdenticalTo
lQ3
b)const{if(&*data==&*b.data)return
true
iF
data->IsIdenticalTo(*b.data);}
xQ1
xL2
xK::IsIdenticalTo
cR1
xL2
xK&b)const{if(Hash!=b.Hash
tX1;if(Opcode!=l34
tX1;switch(Opcode
xI3
cImmed:return
fp_equal(Value,l44
tJ3
l53:return
l42==b.l32
case
cFCall:case
cPCall:if(l42!=b.cL1
return
false
iY2
yM3
yN3
if(nN3!=b.nN3
tX1
xF1
0;a<eH3
iU2
a){if(!eH3[a]xF
b.eH3[a])tX1;}
return
true;}
lB
Become
lQ3
b){if(&b!=this&&&*data!=&*b.data){DataP
tmp=b.data;l81
data.iA2}
}
lB
CopyOnWrite(){if(GetRefCount()>1)data=new
xL2
xK(*data);}
tJ1
CodeTree
xK
CodeTree
xK::GetUniqueRef(){if(GetRefCount()>1)return
CodeTree
xK(*this,CloneTag())iF*this;}
iP):yP
cNop
tH2(),n8
iP
const
xL2&b):yP
l34
tH2(l44,l42(b.cL1,eH3(b.eH3),Hash(b.Hash),Depth(b.Depth),i21
b.l52){}
iP
const
cV1&i):yP
cImmed
tH2(i),n8
#ifdef __GXX_EXPERIMENTAL_CXX0X__
iP
xL2
xK&&b):yP
l34
tH2
c82
l44),l42(b.cL1,eH3
c82
b.eH3)),Hash(b.Hash),Depth(b.Depth),i21
b.l52){}
iP
cV1&&i):yP
cImmed
tH2
c82
i)),n8
#endif
iP
nW2
o):yP
o
tH2(),n8
iP
nW2
o
e02
f):yP
o
tH2(),l42(f),eH3(),Hash(),Depth(1),i21
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
lL3
FUNCTIONPARSERTYPES;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
lL3{eV1
i31
nO,std
c5&done,std::ostream&o){for
tS1
a=tL
i31
t73,done,o);std::ostringstream
buf
yH2
tree,buf);done[tree.GetHash()].insert(buf.str());}
}
#endif
t1{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
eV1
DumpHashes(c9){std
c5
done;i31
tree,done,o);for(std
c5::const_iterator
i=done.y43
i!=done
l24;++i){const
std::set<std
yW3>&flist=i
e42;if(flist.size()!=1)o<<"ERROR - HASH COLLISION?\n"
;for(std::set<std
yW3>::const_iterator
j=flist.y43
j!=flist
l24;++j){o<<'['<<std::hex<<i->first.hash1<<','<<i->first.hash2<<']'<<std::dec;o<<": "
<<*j<<"\n"
;}
}
}
eV1
DumpTree(c9){eT2
iM3;switch
eI2
nC
xI3
cImmed:o<<tree
l31
tI3
l53:o<<"Var"
<<eI2.GetVar()-l53)tI3
cAdd:iM3"+"
;lC
cMul:iM3"*"
;lC
cAnd:iM3"&"
;lC
cOr:iM3"|"
;lC
cPow:iM3"^"
iY2
yM3
iM3;o<<i73
eI2
nC);if
cS3
cFCall||t52==cPCall)o<<':'<<tree.GetFuncNo();}
o<<'(';if
eI2.GetParamCount()<=1&&sep2[1])o<<(sep2+1)<<' '
cZ3
if(a>0)o<<' '
yH2
t73,o
yU3
a+1<tree.lT1
o<<sep2;}
o<<')';}
eV1
DumpTreeWithIndent(c9,const
std
yW3&indent){o<<'['<<std::hex<<(void*)(&tree.l22))<<std::dec<<','<<tree.GetRefCount()<<']';o<<indent<<'_';switch
eI2
nC
xI3
cImmed:o<<"cImmed "
<<tree
l31;o<<'\n'
tI3
l53:o<<"VarBegin "
<<eI2.GetVar()-l53);o<<'\n'
iF;yM3
o<<i73
eI2
nC);if
cS3
cFCall||t52==cPCall)o<<':'<<tree.GetFuncNo();o<<'\n';}
for
tS1
a=tL{std
yW3
ind=indent;for
tS1
p=0;p<ind.size();p+=2)if(ind[p]=='\\')ind[p]=' ';ind+=(a+1<tree.lT1?" |"
:" \\"
;DumpTreeWithIndent(t73,o,ind);}
o<<std::flush;}
#endif
}
#endif
using
lL3
l61;using
lL3
FUNCTIONPARSERTYPES;
#include <cctype>
lL3
l61{unsigned
ParamSpec_GetDepCode
cR1
e22&b){switch(b.first
xI3
ParamHolder:{e6*s=(e6*)b
t63
iF
s->depcode
lD3
SubFunction:{e7*s=(e7*)b
t63
iF
s->depcode;}
yM3
yN3
return
0;}
eV1
DumpParam
cR1
i63,std::ostream&o){static
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
i23
y32
0;i83
i93{const
ParamSpec_NumConstant
xK&yR3
cR1
ParamSpec_NumConstant
xK*t72
using
lL3
FUNCTIONPARSERTYPES;o.precision(12);o<<l74
constvalue
iY2}
case
ParamHolder:{e6&yR3(e6*t72
o<<ParamHolderNames[l74
index];y32
l74
constraints
iY2}
case
SubFunction:{e7&yR3(e7*t72
y32
l74
constraints;yE
GroupFunction){yX3
iW==cNeg){o<<"-"
;n2}
iH1
l84==cInv){o<<"/"
;n2}
else{std
yW3
opcode=i73((nW2)l84).substr(1)xF1
0;a<opcode
iU2
a)opcode[a]=(char)std::toupper(opcode[a]);o<<opcode<<"( "
;n2
o<<" )"
;}
}
else{o<<'('<<i73((nW2)l84)<<' ';yE
PositionalParams)o<<'[';yE
SelectedParams)o<<'{';n2
yX3.data.n1!=0)o<<" <"
<<l94.n1<<'>';yE
PositionalParams)o<<"]"
;yE
SelectedParams)o<<"}"
;o<<')';}
yN3}
switch(ImmedConstraint_Value(constraints&ValueMask)xI3
ValueMask:lC
Value_AnyNum:lC
nP2:o<<"@E"
;lC
Value_OddInt:o<<"@O"
;lC
i51:o<<"@I"
;lC
Value_NonInteger:o<<"@F"
;lC
t11:o<<"@L"
iY2}
switch(ImmedConstraint_Sign(constraints&SignMask)xI3
SignMask:lC
Sign_AnySign:lC
nF1:o<<"@P"
;lC
t21:o<<"@N"
iY2}
switch(ImmedConstraint_Oneness(constraints&OnenessMask)xI3
OnenessMask:lC
Oneness_Any:lC
Oneness_One:o<<"@1"
;lC
Oneness_NotOne:o<<"@M"
iY2}
switch(ImmedConstraint_Constness(constraints&ConstnessMask)xI3
ConstnessMask:lC
i41:if(xD2.first==ParamHolder){e6&yR3(e6*t72
yX3.index<2)yN3
o<<"@C"
;lC
Constness_NotConst:o<<"@V"
;lC
Oneness_Any:yN3}
eV1
DumpParams
l51
paramlist
e02
count,std::ostream&o){for
l51
a=0;a<count;++a){if(a>0)o<<' ';const
e22&param=x91
paramlist,a);DumpParam
xK(param,o)i23
depcode=ParamSpec_GetDepCode(param
yU3
depcode!=0)o<<"@D"
<<depcode;}
}
}
#include <algorithm>
using
lL3
l61;using
lL3
FUNCTIONPARSERTYPES;lL3{e6
plist_p[37]={{2,0,0x0}
nP
0,0x4}
nP
nF1,0x0}
nP
t21|Constness_NotConst,0x0}
nP
Sign_NoIdea,0x0}
nP
t11,0x0}
,{3,Sign_NoIdea,0x0}
,{3,0,0x0}
,{3,t11,0x0}
,{3,0,0x8}
,{3,Value_OddInt,0x0}
,{3,Value_NonInteger,0x0}
,{3,nP2,0x0}
,{3,nF1,0x0}
,{0,t21|lV{0,lV{0,nF1|lV{0,nP2|lV{0,i41,0x1}
,{0,i51|nF1|lV{0,i61
i41,0x1}
,{0,i61
lV{0,Oneness_One|lV{0,t11|lV{1,lV{1,nP2|lV{1,i61
lV{1,i51|lV{1,nF1|lV{1,t21|lV{6,0,0x0}
,{4,0,0x0}
,{4,i51,0x0}
,{4,lV{4,0,0x16}
,{5,0,0x0}
,{5,lV}
n11
cW2
plist_n_container{static
const
ParamSpec_NumConstant
xK
plist_n[20];}
n11
const
ParamSpec_NumConstant
xK
plist_n_container
xK::plist_n[20]={{cV1(-2
tV-1
tV-0.5
tV-0.25
tV
0
tG2
fp_const_deg_to_rad
eL3
fp_const_einv
eL3
fp_const_log10inv
xK(tV
0.5
tG2
fp_const_log2
xK(tV
1
tG2
fp_const_log2inv
xK(tV
2
tG2
fp_const_log10
eL3
fp_const_e
eL3
fp_const_rad_to_deg
eL3-fp_const_pihalf
xK(),xT1{y11,xT1{fp_const_pihalf
xK(),xT1{fp_const_pi
xK(),xT1}
;e7
plist_s[517]={{{1,15,tI2
398,tI2
477,tI2
15,cNeg,GroupFunction,0}
,i41,0x1}
,{{1,15,y42
24,y42
465,y42
466,y42
498,cInv
iV3
327995
xD
l0
2,48276
xD
l6
260151
xD
l6
470171
xD
l6
169126
xD
l6
48418
xD
l6
1328
xD
l6
283962
xD
l6
169275
xD
l6
39202
xD
l6
283964
xD
l6
283973
xD
l6
476619
xD
l6
296998
xD
l6
47
xD
SelectedParams,0}
,0,0x4
nH
161839
xD
l6
25036
xD
l6
35847
xD
l6
60440
xD
l6
30751
xD
l6
270599
xD
l6
60431
xD
l6
259119
xD
l6
183474
xD
l6
332066
xD
l6
7168
xD
l6
197632
xD
l6
291840
xD
l6
283648
xD
l6
238866
xD
l6
239902
xD
l6
31751
xD
l6
244743
xD
l6
384022
xD
SelectedParams,0}
,0,0x4
nH
385262
xD
l6
386086
xD
l6
393254
xD
SelectedParams,0}
,0,0x5
nH
393254
xD
l6
386095
xD
l6
387312
xD
l6
18662
xD
l6
61670
xD
l6
387397
xD
l6
247855
xD
SelectedParams,0}
,0,0x1
nH
342063
xD
l6
297007
xD
l6
15820
xD
l6
393263
xD
l6
393263
xD
SelectedParams,0}
,0,0x5
nH
161847
xD
l6
258103
xD
l6
249073
xD
l6
249076
xD
i9
0,0
xD
nF
0,0
i71
1,45
xD
nF
1,53
xD
nF
1,54
xD
nF
1,55
xD
nF
1,56
xD
nF
1,26
xD
nF
1,259
eA2
0x16}
,{{1,253
xD
nF
1,272
i71
1,323
eA2
0x16
tK3
xD
nF
1,21
xD
nF
1,447
eA2
0x4}
,{{1,449
eA2
0x4
tK3
eA2
0x4
tK3
tH
2}
,0,0x4}
,{{1,15
xD
nF
1,24
tH
2}
,0,0x0
nH
58392
i71
0,0
tH
1}
,nF1,0x0
nH
24591
tL3
33807
tL3
48143
tL3
285720
tL3
290840
tL3
305152,l9
312400,l9
39202,l9
121894,l9
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
l63
327713,l9
322596,l9
88361,l9
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
90127,l9
131087,l9
296976,tO1
324623,l1
0x14
nH
332815,l1
0x10}
,{{3,7340056,tO1
289092,l9
92176
l63
337935
e01
7340060
n13
7340176,l9
338959
e01
7340061
l63
7206,l9
x83
l9
357414,l9
368678,l9
370745,l1
0x7}
,{{3,7340177,l9
39277,tO1
426398
n13
40272286
l63
490910
n13
40336798
l63
50600,l9
426462
l63
490974
l63
370726,l1
0x6
nH
371750,l1
0x6
nH
428070
e01
40336862
l63
38378,l9
50671
e01
47662080,l9
477184,l9
568320,l9
371727,l1
0x7}
,{{3,15779306,l9
370703,l1
0x7
nH
39277,l9
39279,l1
0x4}
,{{3,15779238,l9
39338,tO1
436262,l9
508966,l9
39409,tO1
296998,tO1
35847,l9
15,tO1
377894,l9
386063,l1
0x1
nH
15,l9
7192,l9
122904,l9
121880,l9
30751,l9
57,l9
7456,l9
15674
e01
67579935,l9
39237,l9
58768,l9
62924,l9
121856,l9
15760
e01
64009216,l1
0x0}
,{{0,0,xG
0,0,iM
2,eA1
2,eB1
3,eA1
3,eB1
38,xG
1,38,iM
14,xG
1,57,xG
1,16,eB2
0x0
nH
471103,eB2
0x1}
,{{1,303,xG
1,323,yY3
0x0
nH
471363,eB2
0x16}
,{{1,293,eA1
294,eB1
295,xG
1,296,iM
400,xG
1,0,xG
1,460,xG
1,465,xG
1,16,eB2
0x1}
,{{1,57,yY3
0x1
tK3,iM
21,xG
1,15,eB2
0x0
nH
24591,xG
1,24,iM
517,yY3
0x0
nH
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
50176,lF,178176,t51
tQ3
283648,lF,19456,lF,27648,lF,89088,lF,86016,lF,488448,lF,14342,lF,58375,lF,46147
xN
46151,lF,284679,lF,7183,lF,46159
xN
38993
xN
50265,lF,50249,lF,283808,lF,284835,lF,24822,lF,10240,lF,11264,lF,7170,lF,x83
lF,17408,lF,164864,lF,237568,lF,242688,t51
0x14
nH
476160,lF,25607,lF,121871,lF,50252,lF,39374,lF,50183,lF,7192,lF,121887,lF,252979,lF,46155,lF,38919,lF,50267,lF,50268,lF,50253,lF,46190,lF,50295,lF,7563,t51
0x10
nH
416811,lF,416819,lF,40046,lF,46191
xN
415795,lF,40047
xN
415787,lF,39015,t51
0x5
nH
39326
xN
39326,lF,39332,t51
0x5
nH
39333
e12
50590
xN
50590,lF,39338
xN
39338,lF,39335,t51
0x5
nH
15786
xN
146858,lF,39372,lF,39379,lF,39380,lF,39390
xN
50654
xN
50654,lF,24,t51
0x6
nH
62,lF,24,lF,62,t51
0x6
nH
43,lF,43
xN
51,lF,51
xN
50269,lF,50176
xN
50270,lF,39159,lF,39183
xN
7168
xN
31744,lF,99328,lF,31746,lF,100376,lF,39409
xN
39411
xN
39411,lF,39420,lF,39420
xN
15,lF,39025,t51
0x5
nH
39422,lF,16384,lF,62853,lF,15360,lF,15
e12
16,lF,7183
e12
7172,cPow,y71,nF1,0x0
nH
24591,cPow
iV3
50200,cPow
iV3
63521,cPow
iV3
62500,cPow
iV3
50453,cPow
iV3
62488,cPow,iX3
tN3
7,tN3
194,tN3
0,cAcos
tO3
cAcosh
tO3
cAsin
tO3
cAsinh
nQ
119,cAsinh
tO3
cAtan
t41
306176,cAtan2
t41
x83
cAtan2
tO3
cAtanh
nQ
246,cCeil
tO3
cCeil,tP3
c92
0,cCos,iO
1,7,c92
91,c92
92,c92
119,c92
236,c92
255,c92
214,l73
236,l73
464,l73
0,cCosh,tP3
l73
0,cExp
nQ
7,cExp
nQ
91,cExp
tO3
yZ3
7,yZ3
91,yZ3
246,cFloor
tO3
cFloor,lA
0x4
nH
309540,cHypot
t41
316708,cHypot
t41
316724,cHypot,l0
3,32513024,eC2
34627584
nL1
31493120,eC2
89213952
nL1
149042176
nL1
246647808
nL1
301234176
nL1
494360576
nL1
498558976
nL1
62933520
nL1
62933520,eC2
62933526
nL1
62933526,eC2
24670208
nL1
579378176
nL1
573578240
nL1
32513024
nL1
566254592
nL1
7900160
nL1
588822528,cIf
nQ
119,cInt
nQ
246,tN2
0,tN2
7,tN2
31,tN2
194,tN2
363,tN2
15,cLog,lT
1,24,cLog,iX3
cLog10
tO3
cLog2
t41
x83
cMax
t41
35847,cMax
t41
30751,cMax
tO3
cMax,AnyParams,1}
,0,0x4
nH
x83
cMin
t41
35847,cMin
t41
30751,cMin
tO3
cMin,AnyParams,1}
,0,0x4
nH
24591,cMin,iX3
nQ2
7,nQ2
91,nQ2
92,nQ2
119,nQ2
149,nQ2
231,cSin,lA
0x5}
,{{1,246,nQ2
255,nQ2
254,nQ2
0,cSin,iO
1,273,cSin,lA
0x1}
,{{1,214,y52
231,cSinh,lA
0x5}
,{{1,246,y52
254,y52
255,y52
464,y52
0,cSinh,tP3
y52
15,cSqrt,iX3
cA2
0,cTan,iO
1,115,cTan,iO
1,116,cA2
231,cA2
246,cA2
273,cA2
254,cA2
255,cA2
0,y62
0,cTanh,iO
1,213,y62
231,y62
246,y62
254,y62
255,y62
0,cTrunc
t41
15384,cSub
iV3
15384,cDiv
iV3
476626,cDiv
iV3
121913,tO2
x83
n23
tQ3
x83
tO2
31744,n23
0x20
nH
31751,n23
0x24
nH
31751,tO2
121913,i81
t41
x83
cLess,lA
tQ3
41984,cLess,lA
0x4
nH
41984,cLess
t41
7,cLess
t41
x83
cLessOrEq
t41
296182,cLessOrEq
t41
x83
cGreater,lA
tQ3
41984,cGreater,lA
0x4
nH
41984,cGreater
t41
7,cGreater
t41
x83
yA1
l0
2,296182
cC2
tO3
lS2
245,lS2
7,lS2
550,lS2
553,lS2
554,lS2
556,lS2
31,lS2
559,lS2
15,lS2
560,cNot
t41
7706,n43
x83
n43
35847,n43
30751,n43
463903,n43
466975,cAnd,i9
0,0,cAnd,nF
2,x83
eM3
7706,eM3
35847,eM3
463903,eM3
466975,eM3
30751,cOr,i9
1,0,lT2
91,lT2
131,lT2
245,lT2
215,lT2
246,cDeg
nQ
246,cRad
t41
x83
cAbsAnd,l6
x83
cAbsOr,i9
1,0,cV3
tO3
cAbsNotNot,l0
3,32513024,cC3
lA
0x0}
,}
;}
lL3
l61{const
Rule
grammar_rules[262]={{ProduceNewTree,17,1,0,{1,0,cAbs,eD2
409,{1,146,cAtan,eD2
403
nP
1324,cAtan2,eD2
405
nP
307201,cAtan2
xH
253174
nP
255224,cAtan2
xH
259324
nP
257274,cAtan2,eD2
152,{l14
cCeil
lY1
486,{1,68
iT3
482,{1,122
iT3
483,{1,124
iT3
151,{1,125
iT3
419,{1,123
iT3
0,{1,403,cCos,l2
2,1,246,{l14
cCos,l2
18,1,0,{1,400
iT3
301,{1,406,cCosh,l2
2,1,246,{l14
cCosh,l2
18,1,0,{1,400,cCosh
lY1
458,{1,121,cFloor,eD2
150,{l14
cFloor,tR3
156,{3,7382016,eA
549,{3,8430592,eA
556,{3,8436736,eA
157,{3,42998784,eA
550,{3,42999808,eA
562,{3,43039744,eA
557,{3,49291264,eA
538,{3,49325056,eA
469,{3,1058318,eA
473,{3,1058324,eA
473,{3,9438734,eA
469,{3,9438740,cIf,l2
0,3,32542225,{3,36732434,cIf,l2
0,3,32542231,{3,36732440,cIf,cA3
573,{3,32513026,cIf,cA3
515,{3,455505423,cIf,cA3
515,{3,433506837,cIf
lY1
78,{1,256,c03
69,{1,258,c03
404,{1,72,c03
159,{1,147,cLog,l2
0,1,0
nP
487425,cMax
l3
16,1,445
nP
c13
cMax
l3
0,1,0
nP
483329,cMin
l3
16,1,446
nP
c13
cMin,yZ
0,1,153
nP
24832,cPow,tR3
153
nP
25854,cPow,tR3
154
nP
129039
nY1
32055
nY1
32056
nY1
32057
c23
166288
nP
32137
nY1
33082
c23
7168
nP
12688
c23
7434
nP
12553
c33
435
nP
46146
c33
436
nP
46154
c33
437
nP
46150
c33
169
nP
83983
c33
168
nP
130082
c33
175
nP
133154
c43
476160
nP
471055
c43
274432
nP
273423
c43
251904
nP
266274
c43
251904
nP
263186
c33
171,{l14
lU1
421,{1,68,lU1
151,{1,122,lU1
419,{1,124,lU1
170,{1,125,lU1
482,{1,123,lU1
0,{1,405,lU1
172,{l14
cSinh
lY1
328,{1,404,cSinh
lY1
173,{l14
c63
0,{1,408,c63
176,{1,410,c63
177,{l14
cTanh,l2
0,1,442
nP
449551,tW
1,441
nP
c13
tW
1,167
nP
268549,tW
1,181
nP
276749,tW
1,180
nP
276500,tW
2,190770
nP
189622,tW
2,194748
nP
193723,tW
2,202943
nP
196795,tW
2,59699
nP
298148,tW
2,59714
nP
325815,tW
2,59724
nP
343224
xD
yZ
2,1,337,{1,333
tH
1
tD
336,{1,338
tH
1}
}
,{ReplaceParams,2,1,340
nP
1363
nR
342
nP
1365
nR
463
nP
472524
nR
47
nP
356711
nR
349
nP
200751
nR
360
nP
199727
nR
480
nP
207053
nR
481
nP
208077
nR
417
nP
211144
nR
209
nP
211145
nR
418
nP
215240
nR
212
nP
212329
nR
204
nP
373097
nR
211
nP
372944
nR
217
nP
201944
nR
221
nP
223448
nR
367
nP
508329
nR
219
nP
508126
nR
224
nP
225705
nR
223
nP
225776
nR
365
nP
230825
nR
426
nP
377057
nR
497
nP
377054
nR
497
nP
204201
nR
426
nP
375280
nR
224
nP
375006,cAdd
l3
2,2,407781
nP
233698,cAdd
l3
2,2,59763
nP
233842,tW
1,372
nP
1397,lV1
95
nP
24705,lV1
96
nP
24708,lV1
444
nP
449551,lV1
443
nP
c13
lV1
100
nP
101750,lV1
108
nP
106821,lV1
105
nP
103749
t92
0,2,110607
nP
108869
t92
0,2,107535
nP
109893,lJ
0
tD
112
nP
111634,cMul,SelectedParams,0
tD
567,{1,52,lJ
1
tD
568,{1,42,lJ
1}
}
,{ReplaceParams,2,1,467
nP
45516
xE
356
nP
51555
xE
468
nP
49612
xE
357
nP
47459
xE
429
nP
438699
xE
432
nP
441774
xE
486
nP
498726
xE
494
nP
504870
xE
382
nP
435579
xE
497
nP
435709
xE
426
nP
508287
xE
414
nP
500092
xE
499
nP
352744
xE
345
nP
367092
xE
381
nP
425318
xE
478
nP
425460
xE
47
nP
512501
xE
505
nP
355817
xE
47
nP
516598
xE
507
nP
518182
xE
508
nP
358896
xE
351
nP
388605
xE
511
nP
360939
xE
503
nP
354788
xE
514
nP
525350
xE
510
nP
394342
xE
386
nP
351346
t92
2,2,363004
nP
361968
t92
16,1,117
nP
1157
t92
16,1,118
nP
1158
t92
16,1,402
nP
411024
t92
16,2,58768
nP
1472
t92
16,2,15760
nP
1474
t92
17,1,0,{1,400
t92
17,1,57,{1,14,lJ
0}
}
,{ProduceNewTree,4,1,538
nP
41,n33
c73
0
nP
5167,n33
cD
41984
nP
409641,n33
cD
tP
n33
cD
eU
n33
cD
eV
cEqual
xU1
24849,cEqual
xH
tQ
cEqual
xH
lU2
281873,cEqual
xH
iA
cEqual
xH
l91
n33
c73
562
nP
41,i81,c73
538
nP
5167,i81,cD
41984
nP
409641,i81,cD
tP
i81,cD
eU
i81,cD
eV
i81
xU1
24849,i81
xH
tQ
i81
xH
lU2
281873,i81
xH
iA
i81
xH
l91
i81,cD
tP
c93
eU
c93
eV
cLess,eD2
571
nP
46080,cLess
xU1
24832,cLess
xH
xV1
cLess
xH
tQ
cLess
xH
lU2
tT3
cLess
xH
x11
cLess
xH
iA
cLess
xH
l91
cLess,c83
562
nP
409641,c93
tP
lV2
cD
eU
lV2
cD
eV
lV2
eD2
565
nP
409615,cLessOrEq
xU1
24832,cLessOrEq
xH
xV1
cLessOrEq
xH
tQ
cLessOrEq
xH
lU2
tT3
cLessOrEq
xH
x11
cLessOrEq
xH
iA
cLessOrEq
xH
l91
lV2
c83
562
nP
409647,lV2
cD
tP
cB2
eU
cB2
eV
cGreater,eD2
539
nP
409615,cGreater
xU1
24832,cGreater
xH
xV1
cGreater
xH
tQ
cGreater
xH
lU2
tT3
cGreater
xH
x11
cGreater
xH
iA
cGreater
xH
l91
cGreater,c83
538
nP
409647,cB2
tP
yA1
cD
eU
yA1
cD
eV
yA1
eD2
572
nP
46080
cC2
xU1
24832
cC2
xH
xV1
cD2
tQ
cD2
lU2
281856
cC2
xH
x11
cD2
iA
cD2
l91
yA1
c83
538
nP
409641,yA1
c73
519,{1,137,cNot,cA3
571,{1,2,cNot,l2
0,1,452
nP
c13
n63
0,2,537097,{3,547892744,cAnd,yZ
16,1,566,{1,5,cAnd,AnyParams,1}
}
,{ReplaceParams,16,1,569
nP
13314,n63
16,1,544
nP
553498,n63
16,1,546
nP
462369,n63
16,1,548
nP
466465,n63
0,1,457
nP
c13
xW1
570
nP
13314,xW1
563
nP
8197,xW1
541
nP
553498,xW1
542
nP
462369,xW1
543
nP
466465,xW1
564
nP
143365,cOr,yZ
4,1,525,{1,137,cB3
cA3
572,{1,2,cB3
l4
17,1,0,{1,0,cB3
eD2
537,{1,256,cAbsNotNot,yZ
18,1,531,{1,254,cAbsNotNot,yZ
0,1,572,{3,43039744,cC3
tR3
571,{3,49325056,cC3
cA3
454,{3,32513586,cC3
l2
16,3,32542225,{3,36732434,cC3
y71}
,}
;cW2
grammar_optimize_abslogical_type{y5
9
cI
grammar_optimize_abslogical_type
grammar_optimize_abslogical={9,{34,192,228,238,242,247,254,260,261}
}
;}
cW2
grammar_optimize_ignore_if_sideeffects_type{y5
59
cI
grammar_optimize_ignore_if_sideeffects_type
grammar_optimize_ignore_if_sideeffects={59,{0,20,21,22,23,24,25,26,cG
iQ1
78,cH
cJ
cW2
grammar_optimize_nonshortcut_logical_evaluation_type{y5
56
cI
grammar_optimize_nonshortcut_logical_evaluation_type
grammar_optimize_nonshortcut_logical_evaluation={56,{0,25,cG
iQ1
78,cH
241,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259}
}
;}
cW2
grammar_optimize_recreate_type{y5
22
cI
grammar_optimize_recreate_type
grammar_optimize_recreate={22,{18,55,56,57,80,81,82,83,84,85,117,118,120,121,130,131,132,133,134,135,136,137}
}
;}
cW2
grammar_optimize_round1_type{y5
125
cI
grammar_optimize_round1_type
grammar_optimize_round1={125,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,25,cG
37,38,iQ1
45,46,47,48,49,50,51,52,53,54,58,59,60,61,62,63,64,65,66,67,68,69,70,71,78,79,80,81,82,83,84,85,86,87,88,93,94,95,96,97,98,99,100,101,117,118,119,120,121,122,123,124,125,126,127,128,129,138,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cJ
cW2
grammar_optimize_round2_type{y5
103
cI
grammar_optimize_round2_type
grammar_optimize_round2={103,{0,15,16,17,25,cG
39,40,iQ1
45,46,47,48,49,50,51,52,53,54,59,60,72,73,78,79,86,87,88,89,90,91,92,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,119,122,123,124,125,126,127,128,139,159,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cJ
cW2
grammar_optimize_round3_type{y5
79
cI
grammar_optimize_round3_type
grammar_optimize_round3={79,{74,75,76,77,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,170,171,172,173,174,175,176,177,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,197,198,199,201,202,203,205,206,207,208,209,210,211,213,214,215,217,218,219,220,221,222,223,225,226,227,229,230,231,232,233,234,235}
}
;}
cW2
grammar_optimize_round4_type{y5
12
cI
grammar_optimize_round4_type
grammar_optimize_round4={12,{18,55,56,57,130,131,132,133,134,135,136,137}
}
;}
cW2
grammar_optimize_shortcut_logical_evaluation_type{y5
53
cI
grammar_optimize_shortcut_logical_evaluation_type
grammar_optimize_shortcut_logical_evaluation={53,{0,25,cG
iQ1
78,cH
cJ}
lL3
l61{tJ1
e22
ParamSpec_Extract
l51
paramlist
e02
index){index=(paramlist>>(index*10))&1023;if(index>=57)return
e22(SubFunction,cE2
plist_s[index-57]yU3
index>=37)return
e22(NumConstant,cE2
plist_n_container
xK::plist_n[index-37])iF
e22(ParamHolder,cE2
plist_p[index]);}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <stdio.h>
#include <algorithm>
#include <map>
#include <sstream>
using
lL3
FUNCTIONPARSERTYPES;using
lL3
l61;using
t1;using
e81;lL3{nR1
It,xG3
T,xG3
Comp>t71
MyEqualRange(It
first,It
last,const
T&val,Comp
comp){size_t
len=last-first;while(len>0){size_t
xW3
len/2;It
x93(first);x93+=half;if(comp(*x93,val)){first=x93;++first;len=len-half-1;}
iH1
comp(val,*x93)){len=half;}
else{It
left(first);{It&eE2=left;It
last2(x93);size_t
len2=last2-eE2;while(len2>0){size_t
half2=len2/2;It
y72(eE2);y72+=half2;if(comp(*y72,val)){eE2=y72;++eE2;len2=len2-half2-1;}
else
len2=half2;}
}
first+=len;It
right(++x93);{It&eE2=right;It&last2=first;size_t
len2=last2-eE2;while(len2>0){size_t
half2=len2/2;It
y72(eE2);y72+=half2;if(comp(val,*y72))len2=half2;else{eE2=y72;++eE2;len2=len2-half2-1;}
}
}
return
t71(left,right);}
}
return
t71(first,first);}
tJ1
cW2
OpcodeRuleCompare{iL2()lQ3
tree
e02
y82)const{const
Rule&rule=grammar_rules[y82]iF
t52<rule
cF2.subfunc_opcode;}
iL2()l51
y82,const
eK
const{const
Rule&rule=grammar_rules[y82]iF
rule
cF2.subfunc_opcode<t52;}
}
n11
bool
TestRuleAndApplyIfMatch
cR1
eS2
lR1
tree,bool
c0{MatchInfo
xK
info;n51
found(false,cT()yU3(rule.lA1
LogicalContextOnly)&&!c0{tP1
if(nA
IsIntType
xK::xD3){if(rule.lA1
NotForIntegers)tP1
e43
rule.lA1
OnlyForIntegers)tP1
if(nA
IsComplexType
xK::xD3){if(rule.lA1
NotForComplex)tP1
e43
rule.lA1
OnlyForComplex)tP1
for(;;){
#ifdef DEBUG_SUBSTITUTIONS
#endif
found=TestParams(rule
cF2,tree,found.specs,info,true
yU3
found.found)break;if(!&*found.specs){fail:;
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
c53,false);
#endif
xB2}
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch(rule
c53,true);
#endif
SynthesizeRule(rule
c53)nX2}
e81{xQ1
ApplyGrammar
cR1
Grammar&xX3,lR1
tree,bool
c0{if
eI2.GetOptimizedUsing()==&xX3){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Already optimized:  "
yH2
tree
x62"\n"
<<std::flush;
#endif
xB2
if(true){bool
changed=false;switch
eI2
nC
xI3
cNot:case
cNotNot:case
cAnd:case
cOr:for
tS1
a=0;a<tree.x6
true))yM1
lC
cIf:case
cAbsIf:if(ApplyGrammar(xX3,t5
0),t52==cIf))changed=true
xF1
1;a<tree.x6
c0)yM1
break;yM3
for
tS1
a=0;a<tree.x6
false))yM1}
if(changed){tree.Mark_Incompletely_Hashed()nX2}
typedef
const
unsigned
short*n73;std::pair<n73,n73>range=MyEqualRange(xX3.rule_list,xX3.rule_list+xX3.rule_count,tree,OpcodeRuleCompare
xK());std::vector<unsigned
short>rules;rules.xQ3
range
t63-range.first);for
xZ
if(IsLogisticallyPlausibleParamsMatch(tT1
cF2,tree))rules.push_back(*r);}
range.first=!rules.empty()?&rules[0]:0;range
t63=!rules.empty()?&rules[rules.size()-1]+1:0;if(range.first!=range
t63){
#ifdef DEBUG_SUBSTITUTIONS
if(range.first!=range
t63)t62"Input ("
<<i73
eI2
nC)<<")["
<<tree.GetParamCount()<<"]"
;if(c0
std::cout<<"(Logical)"
i23
first=iR1,prev=iR1;eT2
sep=", rules "
;for
xZ
if(first==iR1)first=prev=*r;iH1*r==prev+1)prev=*r;else
t62
sep<<first;sep=","
;if(prev!=first)std::cout<<'-'<<prev;first=prev=*r;}
}
if(first!=iR1)t62
sep<<first;if(prev!=first)std::cout<<'-'<<prev;}
std::cout<<": "
yH2
tree
x62"\n"
<<std::flush;}
#endif
bool
changed=false;for
xZ
#ifndef DEBUG_SUBSTITUTIONS
if(!IsLogisticallyPlausibleParamsMatch(tT1
cF2,tree))y91
#endif
if(TestRuleAndApplyIfMatch(tT1,tree,c0){yM1
yN3}
if(changed){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Changed."
<<std::endl
iW2"Output: "
yH2
tree
x62"\n"
<<std::flush;
#endif
tree.Mark_Incompletely_Hashed()nX2}
tree.SetOptimizedUsing(&xX3)iF
false;}
xQ1
ApplyGrammar
cR1
void*p,FPoptimizer_CodeTree::eK
yJ
ApplyGrammar(*cR1
Grammar*)p,tree);}
eV1
ApplyGrammars(FPoptimizer_CodeTree::eK{
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round1\n"
;
#endif
n6
grammar_optimize_round1,cX2
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round2\n"
;
#endif
n6
grammar_optimize_round2,cX2
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round3\n"
;
#endif
n6
grammar_optimize_round3,cX2
#ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_nonshortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_nonshortcut_logical_evaluation,cX2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_round4\n"
;
#endif
n6
grammar_optimize_round4,cX2
#ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_shortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_shortcut_logical_evaluation,cX2
#endif
#ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_ignore_if_sideeffects\n"
;
#endif
n6
grammar_optimize_ignore_if_sideeffects,cX2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
iI3"grammar_optimize_abslogical\n"
;
#endif
n6
grammar_optimize_abslogical,cX2
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
lL3
FUNCTIONPARSERTYPES;using
lL3
l61;using
t1;using
e81;lL3{xQ1
TestImmedConstraints
l51
bitmask,const
eK{switch(bitmask&ValueMask
xI3
Value_AnyNum:case
ValueMask:lC
nP2:if(GetEvennessInfo
eI2)!=cG2
Value_OddInt:if(GetEvennessInfo
eI2)!=tP2
i51:if(GetIntegerInfo
eI2)!=cG2
Value_NonInteger:if(GetIntegerInfo
eI2)!=tP2
t11:if(!IsLogicalValue
eI2))cM
SignMask
xI3
Sign_AnySign:lC
nF1:if(l01
cG2
t21:if(l01
tP2
Sign_NoIdea:if(l01
Unknown)cM
OnenessMask
xI3
Oneness_Any:case
OnenessMask:lC
Oneness_One:if(!tree
lB1
if(!fp_equal(fp_abs
eI2
l31),y63
cL
Oneness_NotOne:if(!tree
lB1
if(fp_equal(fp_abs
eI2
l31),y63
cM
ConstnessMask
xI3
Constness_Any:lC
i41:if(!tree
tF1)cL
Constness_NotConst:if
eI2
lB1
yN3
return
true;}
tN1
unsigned
extent
e02
nbits,xG3
eF2=unsigned
int>cW2
nbitmap{private:static
const
unsigned
bits_in_char=8;static
const
unsigned
eG2=(cD3
eF2)*bits_in_char)/nbits;eF2
data[(extent+eG2-1)/eG2];eQ3
void
inc
l41,int
by=1){data[pos(index)]+=by*eF2(1<<y92);nT1
void
dec
l41){inc(index,-1);}
int
get
l41
yU1(data[pos(index)]>>y92)&mask()y83
pos
l41)yJ
index/eG2
y83
shift
l41)yJ
nbits*(index%eG2)y83
mask()yJ(1<<nbits)-1
y83
mask
l41)yJ
mask()<<y92;}
}
;cW2
eG3{int
SubTrees:8;int
Others:8;int
yA2:8;int
Immeds:8;nbitmap<l53,2>SubTreesDetail;eG3(){std::memset(this,0,cD3*this));}
eG3
cR1
eG3&b){std::memcpy(this,&b,cD3
b));}
eG3&eH1=cR1
eG3&b){std::memcpy(this,&b,cD3
b))iF*this;}
}
n11
eG3
CreateNeedList_uncached(tA&c32){eG3
cO1
yP3
a=0;a<c32
yB2;++a){const
i63=x91
c32.param_list,a);i83
SubFunction:{e7&yR3(e7*t72
yE
GroupFunction)++tS2;else{++tU2;assert(param.data.subfunc_opcode<VarBegin);cO1.SubTreesDetail.inc(l84);}
++cO1.yA2
iY2}
case
i93
case
ParamHolder:++tT2;++cO1.yA2
iY2}
}
return
cO1;}
tJ1
eG3&CreateNeedList(tA&c32){typedef
std::map<tA*,eG3>eC1;static
eC1
yP1;eC1::y53
i=yP1.xV2&c32
yU3
i!=yP1.cZ1&c32)return
i
e42
iF
yP1.y13,std::make_pair(&c32,CreateNeedList_uncached
xK(c32)))e42;}
tJ1
CodeTree
xK
CalculateGroupFunction
cR1
i63,const
t9
info){i83
i93{const
ParamSpec_NumConstant
xK&yR3
cR1
ParamSpec_NumConstant
xK*tR2
iF
CodeTreeImmed
l64
constvalue)lD3
ParamHolder:{e6&yR3(e6*tR2
iF
tY3
GetParamHolderValueIfFound
l64
index)lD3
SubFunction:{e7&yR3(e7*tR2
iG
xD3;xD3
xC
l84);iC2
l22).xQ3
l94
yB2)yP3
a=0;a<l94
yB2;++a){CodeTree
xK
tmp(CalculateGroupFunction(x91
l94.param_list,a),info));xD3
c51
tmp);}
iC2
Rehash()iF
xD3;}
}
return
CodeTree
xK();}
}
e81{xQ1
IsLogisticallyPlausibleParamsMatch(tA&c32,const
eK{eG3
cO1(CreateNeedList
xK(c32));size_t
tU3=tree.l21
if(tU3<size_t(cO1.yA2))tQ2
for
tS1
a=0;a<tU3;++a){unsigned
opcode=t73
nC;switch(opcode
xI3
cImmed:if(tS2>0)--tS2;else--tT2;lC
l53:case
cFCall:case
cPCall:--tT2
iY2
yM3
assert(opcode<VarBegin);if(tU2>0&&cO1.SubTreesDetail.get(opcode)>0){--tU2;cO1.SubTreesDetail.dec(opcode);}
else--tT2;}
}
if(tS2>0||tU2>0||tT2>0)tQ2
if(c32.match_type!=AnyParams){if(0||tU2<0||tT2<0)tQ2}
return
true;}
tJ1
n51
TestParam
cR1
i63,x92
tree
eR2
start_at,t9
info){i83
i93{const
ParamSpec_NumConstant
xK&yR3
cR1
ParamSpec_NumConstant
xK*t72
if(!tree
lB1
cV1
imm=tree
l31;switch
l64
modulo
xI3
Modulo_None:lC
Modulo_Radians:imm=fp_mod(imm,yK
imm<xH1
imm
yN
if(imm>fp_const_pi
xK())imm-=fp_const_twopi
xK(iX2
return
fp_equal(imm,l74
constvalue)lD3
ParamHolder:{e6&yR3(e6*t72
if(!x2
return
tY3
SaveOrTestParamHolder
l64
index,tree)lD3
SubFunction:{e7&yR3(e7*t72
yE
GroupFunction){if(!x2
CodeTree
xK
xX1=CalculateGroupFunction(xD2,info);
#ifdef DEBUG_SUBSTITUTIONS
DumpHashes(xX1
x62*cR1
void**)&xX1
l31
iW2"\n"
iW2*cR1
void**)&tree
l31
iW2"\n"
;DumpHashes
eI2
x62"Comparing "
yH2
xX1
x62" and "
yH2
tree
x62": "
iW2(xX1
xF
tree)?"true"
:"false"
x62"\n"
;
#endif
return
xX1
xF
tree);}
e43!tZ1){if(!x2
if
eI2
nC!=l84
tX1;}
return
TestParams
l64
data,tree,start_at,info,false);}
}
}
xB2
tJ1
cW2
iX
nR2
MatchInfo
xK
info;iX(yC2,info(){}
}
n11
class
MatchPositionSpec_PositionalParams:yW1
iX
xK>{eQ3
l83
MatchPositionSpec_PositionalParams
tS1
n
yV1
iX
xK>(n){}
}
;cW2
iS1
nR2
iS1(yC2{}
}
;class
yF:yW1
iS1>{eQ3
unsigned
trypos;l83
yF
tS1
n
yV1
iS1>(n),trypos(0){}
}
n11
n51
TestParam_AnyWhere
cR1
i63,x92
tree
eR2
start_at,t9
info,cE3&used,bool
lG2{xT<yF>x5
i23
lW2
yF*)tZ1;a=x5->trypos;goto
retry_anywhere_2;}
tV2
yF
eI2.lT1;a=0;}
tV3
tree.l21++a){if(used[a])y91
retry_anywhere
cF3
TestParam(xD2,t73,eR3);tY1
used[a]=true;if(lG2
tX3
a);x5->trypos=a
iF
n51(true,&*x5);}
}
retry_anywhere_2
tW3
goto
retry_anywhere;}
}
xB2
tJ1
cW2
yB1
nR2
MatchInfo
xK
info;cE3
used;l83
yB1
tS1
tU3
yC2,info(),used(tU3){}
}
n11
class
MatchPositionSpec_AnyParams:yW1
yB1
xK>{eQ3
l83
MatchPositionSpec_AnyParams
tS1
n,size_t
m
yV1
yB1
xK>(n,yB1
xK(m)){}
}
n11
n51
TestParams(tA&nN,x92
tree
eR2
start_at,t9
info,bool
lG2{if(nN.match_type!=AnyParams){if(xY!=tree.lT1
xB2
if(!IsLogisticallyPlausibleParamsMatch(nN,tree))tQ2
switch(nN.match_type
xI3
PositionalParams:{xT<cA>x5
i23
lW2
cA*)tZ1;a=xY-1;goto
lF1;}
tV2
cA(xY);a=0;}
tV3
xY;++a){tV1=info;retry_positionalparams
cF3
TestParam(cR
a),t73,eR3);tY1
y91}
}
lF1
tW3
info=tV1;goto
retry_positionalparams;}
if(a>0){--a;goto
lF1;}
info=yS3
iF
false;}
if(lG2
for
l51
a=0;a<xY;++a)tX3
a)iF
n51(true,&*x5)lD3
SelectedParams:case
AnyParams:{xT<t2>x5;cE3
used
eI2.lT1;std::vector<unsigned>l93(xY);std::vector<unsigned>yE2(xY)nG1{const
e22
xD2=cR
a);l93[a]=ParamSpec_GetDepCode(xD2);}
{unsigned
b=0
nG1
if(l93[a]!=0)yE2[b++]=a
nG1
if(l93[a]==0)yE2[b++]=a;}
unsigned
lW2
t2*)tZ1;if(xY==0){a=0;goto
retry_anyparams_4;}
a=xY-1;goto
eD1;}
tV2
t2(xY,tree.lT1;a=0;if(xY!=0){yS3=info;(*x5)[0].used=used;}
}
tV3
xY;++a){if(a>0){tV1=info;(*x5)[a].used=used;}
retry_anyparams
cF3
TestParam_AnyWhere
xK(cR
yE2[a]),tree,eR3,used,lG2;tY1
y91}
}
eD1
tW3
info=tV1;used=(*x5)[a].used;goto
retry_anyparams;}
eE1:if(a>0){--a;goto
eD1;}
info=yS3
iF
false;}
retry_anyparams_4:if(nN.n1!=0){if(!TopLevel||!tY3
HasRestHolder(nN.n1)){eH
cH2;cH2.reserve
eI2.lT1
yP3
b=0;b<tree.tA3{if(cH3)y91
cH2.push_back(t5
b));cH3=true;if(lG2
tX3
b);}
if(!tY3
SaveOrTestRestHolder(nN.n1,cH2)){goto
eE1;}
}
else{iM2&cH2=tY3
GetRestHolderValues(nN.n1)xF1
0;a<cH2
iU2
a){bool
found=false
yP3
b=0;b<tree.tA3{if(cH3)y91
if(cH2[a]xF
t5
b))){cH3=true;if(lG2
tX3
b);found=true
iY2}
}
if(!found){goto
eE1;}
}
}
}
return
n51(true,xY?&*x5:0)lD3
GroupFunction:yN3
xB2}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
using
t1;using
e81;lL3{tJ1
CodeTree
xK
xY1
const
i63,t9
info,bool
inner=true){i83
i93{const
ParamSpec_NumConstant
xK&yR3
cR1
ParamSpec_NumConstant
xK*tR2
iF
CodeTreeImmed
l64
constvalue)lD3
ParamHolder:{e6&yR3(e6*tR2
iF
tY3
GetParamHolderValue
l64
index)lD3
SubFunction:{e7&yR3(e7*tR2
iG
tree;t82
l84)yP3
a=0;a<l94
yB2;++a){CodeTree
xK
nparam=xY1
x91
l94.param_list,a),info,true)x12
nparam);}
yX3.data.n1!=0){eH
trees(tY3
GetRestHolderValues
l64
data.n1));tree.AddParamsMove(trees);if
eI2.e71
1){assert(tree.GetOpcode()==cAdd lA4()==cMul lA4()==cMin lA4()==cMax lA4()==cAnd lA4()==cOr lA4()==cAbsAnd lA4()==cAbsOr);tree.e72
0));}
iH1
tree.e71
0){switch
eI2
nC
xI3
cAdd:case
cOr:tree=nC1
0));lC
cMul:case
cAnd:tree=nC1
1));yM3
yN3}
}
if(inner)tree.Rehash()iF
tree;}
}
return
CodeTree
xK();}
}
e81{eV1
SynthesizeRule
cR1
eS2
lR1
tree,t9
info){switch(rule.ruletype
xI3
ProduceNewTree:{tree.Become(xY1
x91
rule.t81
0),info,false)iX2
case
ReplaceParams:yM3{std::vector<unsigned>list=tY3
GetMatchedParamIndexes();std::sort(list.iR2
list
l24)xF1
list.size();a-->0;)tree
x02
list[a])yP3
a=0;a<rule.repl_param_count;++a){CodeTree
xK
nparam=xY1
x91
rule.t81
a),info,true)x12
nparam);}
yN3}
}
}
#endif
#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
#include <cstring>
using
lL3
FUNCTIONPARSERTYPES;using
lL3
l61;using
t1;using
e81;lL3
l61{eV1
DumpMatch
cR1
eS2
x92
tree,const
t9
info,bool
DidMatch,std::ostream&o){DumpMatch(rule
c53,DidMatch?iY3"match"
:iY3"mismatch"
,o);}
eV1
DumpMatch
cR1
eS2
x92
tree,const
t9
info,eT2
tZ3,std::ostream&o){static
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
;o<<tZ3<<" (rule "
<<(&rule-grammar_rules)<<")"
<<":\n  Pattern    : "
;{e22
tmp;tmp.first=SubFunction;ParamSpec_SubFunction
tmp2;tmp2.data=rule
cF2;tmp
t63=cE2
tmp2;DumpParam
xK(tmp,o);}
o<<"\n  Replacement: "
;DumpParams
xK(rule.t81
rule.repl_param_count
tK2
o<<"  Tree       : "
yH2
tree
tK2
if(!std::strcmp(tZ3,iY3"match"
))DumpHashes
eI2,o)xF1
0;a<tY3
c8
size();++a){if(!tY3
paramholder_matches[a]nJ2))y91
o<<"           "
<<ParamHolderNames[a]<<" = "
yH2
tY3
paramholder_matches[a]tK2}
tR1
tY3
lQ
iU2
b){if(!tY3
lQ[b
t43)continue
xF1
0;a<tY3
lQ[b
t53
iU2
a){o<<"         <"
<<b<<"> = "
yH2
tY3
lQ[b
t53[a],o);o<<std::endl;}
}
o<<std::flush;}
}
#endif
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
lL3
FUNCTIONPARSERTYPES;lL3{xQ1
MarkIncompletes(FPoptimizer_CodeTree::eK{if
eI2.Is_Incompletely_Hashed(e21
bool
iT1=false
xF1
tL
iT1|=MarkIncompletes(t73
yU3
iT1)tree.Mark_Incompletely_Hashed()iF
iT1;}
eV1
FixIncompletes(FPoptimizer_CodeTree::eK{if
eI2.Is_Incompletely_Hashed()){for
tS1
a=tL
FixIncompletes(t73);tree.lH3}
}
}
t1{lB
Sort()eP3
Sort();}
lB
Rehash(bool
constantfolding){if(constantfolding)ConstantFolding(*this);else
Sort();data
x7
tJ1
cW2
c6
cI3
cV1
cJ3
yC1=0;
#if 0
long
double
value=Value;e2=crc32::calc(cR1
unsigned
char*)&value,cD3
value));key^=(key<<24);
#elif 0
union{cW2{unsigned
char
filler1[16];cV1
v
i23
char
filler2[16];}
buf2;cW2{unsigned
char
filler3[cD3
cV1)+16-cD3
xD1)];e2;}
buf1;}
data;memset(&data,0,cD3
data));data.buf2.v=Value;e2=data.buf1.key;
#else
int
eX3;cV1
nE2=std::frexp(Value,&i72
e2=l51
eW3+0x8000)&0xFFFF
yU3
nE2<0){nE2=-nE2;key=key^0xFFFF;}
else
key+=0x10000;nE2-=cP3;key<<=39;key|=n61(nE2+nE2)*cV1(1u<<31))<<8;
#endif
lP
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nR1
T
cI2
std::complex<T> >cI3
std::complex<T>cJ3
c6<T>::n83
NewHash,Value.real());nA
fphash_t
temp;c6<T>::n83
temp,Value.imag());yC1^=temp.hash2;NewHash.hash2^=temp.hash1;}
}
;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
tN1
cI2
long>yC
long
Value){e2=Value;lP
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
tN1
cI2
GmpInt>cI3
GmpInt
cJ3
e2=Value.toInt();lP
#endif
eV1
xL2
xK::Recalculate_Hash_NoRecursion(){fphash_t
NewHash(n61
Opcode)<<56,Opcode*iL3(0x1131462E270012B));Depth=1;switch(Opcode
xI3
cImmed:{c6
xK::n83
NewHash,Value
iX2
case
l53:{yC1|=n61
cL1<<48
yF1((n61
cL1)*11)^iL3(0x3A83A83A83A83A0)iY2}
case
cFCall:case
cPCall:{yC1|=n61
cL1<<48
yF1((~n61
cL1)*7)^3456789;}
yM3{size_t
t91=0
xF1
0;a<eH3
iU2
a){if(eH3[a]nD2>t91)t91=eH3[a]nD2;yC1+=((eH3[a].tW2
hash1*(a+1))>>12)yF1
eH3[a].tW2
hash1
yF1(3)*iL3(0x9ABCD801357);NewHash.hash2*=iL3(0xECADB912345)yF1(~eH3[a].tW2
hash2)^4567890;}
Depth+=t91;}
}
if(Hash!=NewHash){Hash=NewHash;l52=0;}
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
lL3
FUNCTIONPARSERTYPES;lL3{using
t1
n11
bool
x21
x92
tree,long
count,const
xK1::SequenceOpCode
xK&eJ,xK1
cL3&synth,size_t
max_bytecode_grow_length);static
const
cW2
SinCosTanDataType{OPCODE
whichopcode;OPCODE
inverse_opcode;enum{nominator,nF2,inverse_nominator,lG1}
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
,{cJ2{cSinh,cCosh,cK2,{cSinh,cNop,{cJ2
cNop,cCosh}
}
,{cCosh,cNop,{cSinh,cJ2
cNop}
}
,{cNop,cTanh,{cCosh,cSinh,cK2,{cNop,cSinh,{cNop,cTanh,cCosh,cNop}
}
,{cNop,cCosh,{cTanh,cSinh,cK2}
;}
t1{lB
SynthesizeByteCode(std::vector<unsigned>&ByteCode,std::vector
xK&Immed,size_t&stacktop_max){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Making bytecode for:\n"
;iU
#endif
while(RecreateInversionsAndNegations()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"One change issued, produced:\n"
;iU
#endif
FixIncompleteHashes();using
e81;using
lL3
l61;const
void*g=cE2
grammar_optimize_recreate;while(ApplyGrammar(*cR1
Grammar*)g,*this)){FixIncompleteHashes();}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Actually synthesizing, after recreating inv/neg:\n"
;iU
#endif
xK1
cL3
synth;SynthesizeByteCode(synth,false);e73
Pull(ByteCode,Immed,stacktop_max);}
lB
SynthesizeByteCode(xK1
cL3&synth,bool
MustPopTemps)const{xZ1*this))yJ;}
for
tS1
a=0;a<12;++a){const
SinCosTanDataType&data=SinCosTanData[a];if(data.whichopcode!=cNop)iL1!=data.whichopcode)y91
CodeTree
n93;n93.lO1
n93
xC
data.inverse_opcode);n93.lX2
xZ1
n93)){e73
l11
else
iL1!=cInv)y91
if(GetParam(0)nC!=data.inverse_opcode)y91
xZ1
GetParam(0))){e73
l11
size_t
found[4];tR1
4;++b){CodeTree
tree;if(data.i03]==cNop){t82
cInv);CodeTree
nA3;nA3.lO1
nA3
xC
data.i03^2]);nA3.lX2
tree
c51
nA3);}
else{tree.lO1
t82
data.i03]);}
tree.lX2
found[b]=e73
eZ3
eI2);}
if
iQ
yF2!=tE
nF2]xX
yF2
x72
nF2
i1
cDiv
x42
yF2!=tE
lG1]xX
yF2
x72
lG1
i1
cMul
x42
lW1!=tE
lG1]xX
lW1
x72
lG1
i1
cRDiv
x42
lW1!=tE
nF2]xX
lW1
x72
nF2
i1
cMul,2,1);e73
l11
size_t
n_subexpressions_synthesized=SynthCommonSubExpressions(synth);switch(iU1{case
l53:e73
PushVar(GetVar());lC
cImmed:x82
GetImmed());lC
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:iL1==cMul){bool
cK3=false;yG
lX1
tF1&&isLongInteger(lX1
l31)){yY1=makeLongInteger(lX1
l31);CodeTree
tmp(*this,xG3
CodeTree::CloneTag());tmp
x02
a);tmp.lH3
if(x21
tmp,value,xK1::i01
xK::AddSequence,synth,MAX_MULI_BYTECODE_LENGTH)){cK3=true
iY2}
}
}
if(cK3)yN3
int
yD1=0;cE3
done(GetParamCount(),false);CodeTree
iB;iB
xC
iU1;for(;;){bool
found=false;yG
done[a])y91
if(e73
IsStackTop(lX1)){found=true;done[a]=true;lX1.n7
iB
tM
lX1
yU3++yD1>1){synth
yH
2);iB.Rehash(false)tQ1
iB);yD1=yD1-2+1;}
}
}
if(!found)yN3
yG
done[a])y91
lX1.n7
iB
tM
lX1
yU3++yD1>1){synth
yH
2);iB.Rehash(false)tQ1
iB);yD1=yD1-2+1;}
}
if(yD1==0){switch(iU1{case
cAdd:case
cOr:case
cAbsOr:x82
0);lC
cMul:case
cAnd:case
cAbsAnd:x82
1);lC
cMin:case
cMax:x82
0)iY2
yM3
yN3++yD1;}
assert(n_stacked==1)iY2}
case
cPow:{lM3
p0
xJ2(0);lM3
p1
xJ2(1
yU3!p1
tF1||!isLongInteger(p1
l31)||!x21
p0,makeLongInteger(p1
l31),xK1::i01
xK::MulSequence,synth,MAX_POWI_BYTECODE_LENGTH)){p0.n7
p1
yO3
yH
2);eX1
cIf:case
cAbsIf:{xG3
xK1
cL3::IfData
ifdata;GetParam(0)yO3.SynthIfStep1(ifdata,iU1;GetParam(1)yO3.SynthIfStep2(ifdata);GetParam(2)yO3.SynthIfStep3(ifdata
iX2
case
cFCall:case
cPCall:{for
tS1
a=0;a<l21++a)lX1
yO3
yH
l51)lT1;nG3
nA2|GetFuncNo(),0,0
iX2
yM3{for
tS1
a=0;a<l21++a)lX1
yO3
yH
l51)lT1
iY2}
}
e73
StackTopIs(*this
yU3
MustPopTemps&&n_subexpressions_synthesized>0){size_t
top=e73
GetStackTop();e73
DoPopNMov(top-1-n_subexpressions_synthesized,top-1);}
}
}
lL3{xQ1
x21
x92
tree,long
count,const
xK1::SequenceOpCode
xK&eJ,xK1
cL3&synth,size_t
max_bytecode_grow_length){if
tS3!=0){xK1
cL3
backup=synth;tree.n7
size_t
bytecodesize_backup=e73
GetByteCodeSize();xK1::x21
count
nF3
size_t
bytecode_grow_amount=e73
GetByteCodeSize()-bytecodesize_backup;if(bytecode_grow_amount>max_bytecode_grow_length){synth=backup
iF
false;}
return
true;}
else{xK1::x21
count,eJ,synth)nX2}
}
#endif
#include <cmath>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
lL3
FUNCTIONPARSERTYPES;lL3{using
t1;
#define FactorStack std::vector
const
cW2
PowiMuliType{unsigned
opcode_square
i23
opcode_cumulate
i23
opcode_invert
i23
opcode_half
i23
opcode_invhalf;}
iseq_powi={cSqr,cMul,cInv,cSqrt,cRSqrt}
,iseq_muli={iR1
xD
cNeg,iR1,iR1}
n11
cV1
cM1
const
PowiMuliType&cN3,const
std
cM3,lY2&stack)i43
1);while(IP<limit){if(iV1==cN3.opcode_square){if(!t12
i53
2;yV
opcode_invert){xD3=-xD3;yV
opcode_half){if(xD3>y11&&isEvenInteger(i53
cP3;yV
opcode_invhalf){if(xD3>y11&&isEvenInteger(i53
cV1(-0.5);++IP;y91}
size_t
nS2=IP;cV1
lhs(1
yU3
iV1==cFetch){unsigned
index=n12;if(index<y4||size_t(index-y4)>=iX1){IP=nS2
iY2}
lhs=stack[index-y4];goto
yG2;}
if(iV1==cDup){lhs=xD3;goto
yG2;yG2:tA1
xD3);++IP;cV1
subexponent=cM1
cN3
tR
if(IP>=limit||iV1!=cN3.opcode_cumulate){IP=nS2
iY2}
++IP;stack.pop_back();xD3+=lhs*subexponent;y91}
yN3
return
xD3;}
tJ1
cV1
ParsePowiSequence
cR1
std
cM3){lY2
stack;tA1
e53
iF
cM1
iseq_powi
tR}
tJ1
cV1
ParseMuliSequence
cR1
std
cM3){lY2
stack;tA1
e53
iF
cM1
iseq_muli
tR}
tJ1
class
CodeTreeParserData{eQ3
l83
CodeTreeParserData(bool
k_powi):stack(),clones(),keep_powi(k_powi){}
void
Eat
tS1
tU3,OPCODE
opcode
e11;xQ
xC
opcode);eH
c32=Pop(tU3);xQ.SetParamsMove(c32
yU3!keep_powi)switch(opcode
xI3
cTanh:{CodeTree
xK
sinh,cosh;sinh
xC
cSinh);sinh
tM
xQ
cO3
sinh.lH3
cosh
xC
cCosh);cosh
c51
xQ
cO3
cosh.t32
pow
tL1
c51
cosh);pow.yB
cV1(-1)))tA2
xQ
xC
nB3
SetParamMove(0,sinh);xQ
c51
pow
iX2
case
cTan:{CodeTree
xK
sin,cos;sin
xC
cSin);sin
tM
xQ
cO3
sin.lH3
cos
xC
cCos);cos
c51
xQ
cO3
cos.t32
pow
tL1
c51
cos);pow.yB
cV1(-1)))tA2
xQ
xC
nB3
SetParamMove(0,sin);xQ
c51
pow
iX2
case
cPow:{x92
p0=xQ
l8
0);x92
p1=xQ
l8
1
yU3
p1
nC==cAdd){eH
xF3(p1.lT1
xF1
0;a<p1.l21++a){CodeTree
xK
pow
tL1
tM
p0);pow
tM
p1
l8
a))tA2
xF3[a
i13
pow);}
xQ
xC
nB3
SetParamsMove
iK2);}
yN3
yM3
yN3
xQ.Rehash(!keep_powi);iW1,false);
#ifdef DEBUG_SUBSTITUTIONS
tC1<<tU3<<", "
<<i73(opcode)<<"->"
<<i73(xQ
nC)<<": "
iQ3
xQ
tS
xQ);
#endif
tA1
nC3
EatFunc
tS1
tU3,OPCODE
opcode
e02
funcno
e11=CodeTreeFuncOp
xK(opcode,funcno);eH
c32=Pop(tU3);xQ.SetParamsMove(c32);xQ.lX2
#ifdef DEBUG_SUBSTITUTIONS
tC1<<tU3<<", "
iQ3
xQ
tS
xQ);
#endif
iW1);tA1
nC3
AddConst(yJ1
e11=CodeTreeImmed(value);iW1);Push(nC3
AddVar
l51
varno
e11=CodeTreeVar
xK(varno);iW1);Push(nC3
SwapLastTwoInStack(){tB1
1
i13
tB1
2])nD3
Dup(){Fetch(iX1-1)nD3
Fetch
tS1
which){Push(stack[which]);}
nR1
T>void
Push(T
tree){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<iQ3
tree
tS
tree);
#endif
tA1
tree)nD3
PopNMov
tS1
target,size_t
source){stack[target]=stack[source];stack.eV3
target+1);}
CodeTree
xK
yJ2{clones.clear()iG
xD3(stack.back());stack.eV3
iX1-1)iF
xD3;}
eH
Pop
tS1
n_pop){eH
xD3(n_pop)yP3
n=0;n<n_pop;++n)xD3[n
i13
tB1
n_pop+n]);
#ifdef DEBUG_SUBSTITUTIONS
for
tS1
n=n_pop;n-->0;){tC1
yH2
xD3[n]tS
xD3[n]);}
#endif
stack.eV3
iX1-n_pop)iF
xD3;}
size_t
GetStackTop(yU1
iX1;}
private:void
FindClone(lR1,bool=true)yJ;}
private:eH
stack;std::multimap<fphash_t,CodeTree
xK>clones;bool
keep_powi;private:CodeTreeParserData
cR1
CodeTreeParserData&);CodeTreeParserData&eH1=cR1
CodeTreeParserData&);}
n11
cW2
IfInfo{CodeTree
xK
eH2
iG
thenbranch;size_t
endif_location;IfInfo():eH2(),thenbranch(),endif_location(){}
}
;}
t1{lB
GenerateFrom
cR1
xG3
FunctionParserBase
xK::Data&xA3,bool
keep_powi){eH
xE2;xE2.xQ3
xA3.mVariablesAmount)yP3
n=0;n<xA3.mVariablesAmount;++n){xE2.push_back(CodeTreeVar
xK(n+l53));}
GenerateFrom(xA3,xE2,keep_powi);}
lB
GenerateFrom
cR1
xG3
FunctionParserBase
xK::Data&xA3,const
x3&xE2,bool
keep_powi){const
std::vector<unsigned>&ByteCode=xA3.mByteCode;const
std::vector
xK&Immed=xA3.mImmed;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"ENTERS GenerateFrom()\n"
;
#endif
CodeTreeParserData
xK
sim(keep_powi);std::vector<IfInfo
xK>eG;for
tS1
IP=0,DP=0;;++IP){tX2:while(!eG.empty()&&(eG.e3==IP||(IP<xF2&&iV1==cJump&&eG.eF1
nJ2)))){CodeTree
elsebranch=sim.yJ2
xK2
eG.back().eH2)xK2
eG.eF1)xK2
elsebranch
y6
3,cIf);eG.pop_back();}
if(IP>=xF2)break
i23
opcode=iV1;if((opcode==cSqr||opcode==cDup||(opcode==cInv&&!IsIntType
xK::xD3)||opcode==cNeg||opcode==cSqrt||opcode==cRSqrt||opcode==cFetch)){size_t
was_ip=IP;cV1
e52
ParsePowiSequence
xK(ByteCode,IP,eG.empty()?xF2:eG.e3,sim.xI
1
yU3
eX3!=cV1(1.0)){xL
eX3
y6
2,cPow);goto
tX2;}
if(opcode==cDup||opcode==cFetch||opcode==cNeg){cV1
xU2=ParseMuliSequence
xK(ByteCode,IP,eG.empty()?xF2:eG.e3,sim.xI
1
yU3
xU2!=cV1(1.0)){xL
xU2
y6
2,cMul);goto
tX2;}
}
IP=was_ip;}
if(n02>=l53){unsigned
index=opcode-l53
xK2
xE2[index]);}
else{switch(n02
xI3
cIf:case
cAbsIf:{eG.eV3
eG.size()+1);CodeTree
res(sim.yJ2);eG.back().eH2.swap(res);eG.e3=xF2;IP+=2;y91}
case
cJump:{CodeTree
res(sim.yJ2);eG.eF1.swap(res);eG.e3=lG3
IP+1]+1;IP+=2;y91}
case
cImmed:xL
Immed[DP++]);lC
cDup:sim.Dup();lC
cNop:lC
cFCall:{unsigned
funcno=n12;assert(funcno<fpdata.mFuncPtrs.size())i23
c32=xA3.mFuncPtrs[funcno].mParams;sim.EatFunc(c32,n02,funcno
iX2
case
cPCall:{unsigned
funcno=n12;assert(funcno<fpdata.iO3.size());const
FunctionParserBase
xK&p=*xA3.iO3[funcno].mParserPtr
i23
c32=xA3.iO3[funcno].mParams;x3
paramlist=sim.Pop(c32);CodeTree
tY2;tY2.GenerateFrom(*p.mData,paramlist)xK2
tY2
iX2
case
cInv:xL
1
yI2
cDiv);lC
cNeg
xA2
cNeg)iY2
xL
0
yI2
cSub);lC
cSqr:xL
2
yE1
cSqrt:xL
cP3
yE1
cRSqrt:xL
cV1(-0.5)yE1
cCbrt:xL
cV1(1)/cV1(3)yE1
cDeg:xL
fp_const_rad_to_deg
x31
cRad:xL
fp_const_deg_to_rad
x31
cExp:cN1)goto
nI3;xL
fp_const_e
xK()yI2
cPow);lC
cExp2:cN1)goto
nI3;xL
2.0
yI2
cPow);lC
cCot
xA2
cTan
nT
cCsc
xA2
cSin
nT
cSec
xA2
cCos
nT
cInt:
#ifndef __x86_64
cN1
lJ3
1,cInt
iX2
#endif
xL
cP3
nH3
y6
1,cFloor);lC
cLog10
xA2
cQ3
fp_const_log10inv
x31
cLog2
xA2
cQ3
fp_const_log2inv
x31
eB3:sim.nW
1,cQ3
fp_const_log2inv
xK()y6
3,cMul);lC
cHypot:xL
2
y6
2,cPow);iZ3
xL
2
y6
2,cPow
nH3);xL
cP3
yE1
cSinCos:sim.Dup(y6
1,cSin);sim.nW
1,cCos);lC
cSinhCosh:sim.Dup(y6
1,cSinh);sim.nW
1,cCosh);lC
cRSub:iZ3
case
cSub:cN1
lJ3
2,cSub
iX2
xL-1
y6
2,cMul
nH3);lC
cRDiv:iZ3
case
cDiv:cN1||IsIntType
xK::xD3
lJ3
2,cDiv
iX2
xL-1
y6
2,cPow
y6
2,cMul);lC
cAdd:case
cMul:case
cMod:case
cPow:case
cEqual:case
cLess:case
cGreater:case
i81:case
cLessOrEq:case
cGreaterOrEq:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:sim.Eat(2,tD1
lC
cNot:case
cNotNot:case
cV3:case
cAbsNotNot
xA2
tD1
lC
cFetch:sim.Fetch(n12);lC
cPopNMov:{unsigned
stackOffs_target=n12
i23
stackOffs_source=n12;sim.PopNMov(stackOffs_target,stackOffs_source
iX2
yM3
nI3:i23
funcno=opcode-cAbs;assert(funcno<FUNC_AMOUNT);const
FuncDefinition&func=Functions[funcno];sim.Eat(func.c32,tD1
yN3}
}
Become(sim.yJ2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Produced tree:\n"
;iU
#endif
}
}
#endif
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
#include <assert.h>
#define FP_MUL_COMBINE_EXPONENTS
lL3{using
lL3
FUNCTIONPARSERTYPES;using
t1
n11
static
void
AdoptChildrenWithSameOpcode(eK{
#ifdef DEBUG_SUBSTITUTIONS
bool
nG2=false;
#endif
for
cK
if(t73
nC==t52){
#ifdef DEBUG_SUBSTITUTIONS
if(!nG2)t62"Before assimilation: "
cN
nG2=true;}
#endif
tree.AddParamsMove(t73.GetUniqueRef().l22),a);}
#ifdef DEBUG_SUBSTITUTIONS
if(nG2)t62"After assimilation:   "
cN}
#endif
}
}
t1{eV1
ConstantFolding(eK{tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
void*cR3=0
iW2"["
<<(&cR3)<<"]Runs ConstantFolding for: "
cN
DumpHashes
eI2
x62
std::flush;
#endif
if(false){redo:;tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cR3)<<"]Re-runs ConstantFolding: "
cN
DumpHashes
eI2);
#endif
}
if
eI2
nC!=cImmed){yG3
p=iT
tree
yU3
p
iN&&p
nH2.yQ==p
cO
val){xO
p.yQ);nD}
if(false){ReplaceTreeWithOne:xO
e53;goto
do_return;ReplaceTreeWithZero:xO
xH1;goto
do_return;ReplaceTreeWithParam0:
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before replace: "
iW2
std::hex<<'['
eW1
hash1<<','
eW1
hash2<<']'<<std::dec
cN
#endif
tree.e72
0));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After replace: "
iW2
std::hex<<'['
eW1
hash1<<','
eW1
hash2<<']'<<std::dec
cN
#endif
i3
eJ2
nC
xI3
cImmed:lC
l53:lC
cAnd:case
cAbsAnd
t4
bool
c1=false;for
cK{if(!xL3
a)))c1=true;eK2
a),t52==cAbsAnd)xI3
eO3
t3
IsAlways:nN1);lC
n01
eJ2.lT1{case
0:i2
1:t82
t52==cAnd?cNotNot:cAbsNotNot);i3
yM3
if
cS3
cAnd||!c1)if(ConstantFolding_AndLogic
cT3
eX1
cOr:case
cAbsOr
t4
bool
c1=false;for
cK{if(!xL3
a)))c1=true;eK2
a),t52==cAbsOr))i91
i2
lE3
nN1);lC
n01
eJ2.lT1{case
0
t3
1:t82
t52==cOr?cNotNot:cAbsNotNot);i3
yM3
if
cS3
cOr||!c1)if(ConstantFolding_OrLogic
cT3
eX1
cNot:case
cV3:{unsigned
n71
0;switch(t5
0)nC
xI3
cEqual:n71
i81;lC
i81:n71
cEqual;lC
cLess:n71
cGreaterOrEq;lC
cGreater:n71
cLessOrEq;lC
cLessOrEq:n71
cGreater;lC
cGreaterOrEq:n71
cLess;lC
cNotNot:n71
cNot;lC
cNot:n71
cNotNot;lC
cV3:n71
cAbsNotNot;lC
cAbsNotNot:n71
cV3
iY2
yM3
yN3
if(opposite){t82
OPCODE(opposite));tree.SetParamsMove(t5
0).GetUniqueRef().l22));i3}
eK2
0),tree
cP1
xI3
IsAlways
t3
lE3
i2
n01
if
cS3
cNot&&GetPositivityInfo(t5
0))==IsAlways)t82
cV3
yU3
nM3
cIf||nM3
tB3{CodeTree
xK
lA3=t5
0);x92
ifp1=lA3
l8
1);x92
ifp2=lA3
l8
2
yU3
ifp1
l13
ifp1
cP1{tree
y0
ifp1
nC==cNot?cNotNot:cAbsNotNot);tZ2
cO3
nJ3
cW3
iJ1
cX3
n81
if(ifp2
l13
ifp2
cP1{tree
y0
t52);tZ2);nJ3
cW3
xC
ifp2
nC==cNot?cNotNot:cAbsNotNot);cX3
l8
0)n81
eX1
cNotNot:case
cAbsNotNot:{if(xL3
0)))cY3
eK2
0),t52==cAbsNotNot)xI3
eO3
t3
IsAlways:i2
n01
if
cS3
cNotNot&&GetPositivityInfo(t5
0))==IsAlways)t82
cAbsNotNot
yU3
nM3
cIf||nM3
tB3{CodeTree
xK
lA3=t5
0);x92
ifp1=lA3
l8
1);x92
ifp2=lA3
l8
2
yU3
ifp1
l13
ifp1
cP1{eQ1
0,lA3
cO3
tree
tM
ifp1)cW3
iJ1
cX3
n81
if(ifp2
l13
ifp2
cP1{tree
y0
t52);tZ2);nJ3;tree
tM
ifp2);t82
lA3
nC);i3}
eX1
cIf:case
cAbsIf:{if(ConstantFolding_IfOperations
cT3
break
lD3
cMul:{NowWeAreMulGroup:;AdoptChildrenWithSameOpcode
eI2);cV1
nH1=cV1(1);size_t
iY1=0;bool
nI1=false
cZ3
if(!t5
tE1
y91
cV1
immed=t93;if(immed==xH1
goto
ReplaceTreeWithZero;nH1*=immed;++iY1;}
if(iY1>1||(iY1==1&&fp_equal(nH1,y63)nI1=true;if(nI1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cMul: Will add new "
iP3
nH1<<"\n"
;
#endif
for
cK
if(t5
tE1{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iP3
t93
iW2"\n"
;
#endif
e03(!fp_equal(nH1,y63
tree
tM
e91
nH1));eJ2.lT1{case
0:i2
1:cY3
yM3
if(ConstantFolding_MulGrouping
cT3
if(ConstantFolding_MulLogicItems
cT3
eX1
cAdd
t4
cV1
n22=0.0;size_t
iY1=0;bool
nI1=false
cZ3
if(!t5
tE1
y91
cV1
immed=t93;n22+=immed;++iY1;}
if(iY1>1||(iY1==1&&n22==xH1)nI1=true;if(nI1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cAdd: Will add new "
iP3
n22<<"\n"
iW2"In: "
cN
#endif
for
cK
if(t5
tE1{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
iP3
t93
iW2"\n"
;
#endif
e03(!(n22==cV1(0.0)))tree
tM
e91
n22));eJ2.lT1{case
0
t3
1:cY3
yM3
if(ConstantFolding_AddGrouping
cT3
if(ConstantFolding_AddLogicItems
cT3
eX1
cMin
t4
size_t
yK2=0;yG3
cX
cZ3
while(a+1<tree.GetParamCount()&&t73
xF
t5
a+1)))nN1+1);yI3
max
xB3&&(!cX
cO
known||(p
e23
cX
cO
val)){cX
cO
val=p
cO
val;cX
cO
known=true;yK2=a;}
}
if(cX
eU1
for
cK{yI3
min
xB3&&a!=yK2&&p.yQ>=cX
cO
val)e03
eI2.e71
1){cY3
eX1
cMax
t4
size_t
yK2=0;yG3
eW
cZ3
while(a+1<tree.GetParamCount()&&t73
xF
t5
a+1)))nN1+1);yI3
min
xB3&&(!eW
iN||p.yQ>eW.yQ)){eW.yQ=p.yQ;eW
iN=true;yK2=a;}
}
if(eW
iN){for
cK{yI3
max
xB3&&a!=yK2&&(p
e23
eW.yQ){nN1);}
}
}
if
eI2.e71
1){cY3
eX1
cEqual:case
i81:case
cLess:case
cGreater:case
cLessOrEq:case
cGreaterOrEq:if(ConstantFolding_Comparison
cT3
lC
cAbs:{yG3
p0
eM
0));if
eP
cY3
if(p0
cP{t82
cMul);tree.yB
y63;goto
NowWeAreMulGroup;}
if(nM3
cMul){x92
p=t5
0);eH
nK3;eH
cL2
xF1
0;a<p.l21++a){p0=iT
p
nL3
if
eP{nK3.push_back(p
nL3}
if(p0
cP{cL2.push_back(p
nL3}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Abs: mul group has "
<<nK3.size()<<" pos, "
<<cL2.size()<<"neg\n"
;
#endif
if(!nK3.empty()||!cL2.empty()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-Before: "
yH2
tree
x62"\n"
<<std::flush;DumpHashes
eZ1
#endif
CodeTree
xK
eN3;eN3
xC
cMul)xF1
0;a<p.l21++a){p0=iT
p
nL3
if(eP||(p0
cP){}
else
eN3
tM
p
nL3}
eN3.t32
nO3;nO3
xC
cAbs);nO3
c51
eN3);nO3.t32
y31
cMul);xF3
c51
nO3);c61
AddParamsMove(nK3
yU3!cL2.empty()){if(cL2.size()%2)c61
yB
cV1(-1)));c61
AddParamsMove(cL2);}
tree.Become
iK2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-After: "
;DumpTree
eZ1
std::cout<<"\n"
<<std::flush;DumpHashes
eZ1
#endif
goto
NowWeAreMulGroup;}
}
yN3
#define HANDLE_UNARY_CONST_FUNC(funcname) nS){xO funcname(lR));nD
case
cLog:iF3(fp_log);if(nM3
cPow){CodeTree
xK
pow=t5
0
yU3
GetPositivityInfo(pow
l8
0))==IsAlways){pow.l81
pow
lP1
tree.lU
if(GetEvennessInfo(pow
l8
1))==IsAlways){pow.CopyOnWrite()iG
abs;abs
xC
cAbs);abs
c51
pow
cO3
abs.lH3
pow
lP1
pow.SetParamMove(0,abs);tree.lU}
iH1
nM3
cAbs){CodeTree
xK
pow=t5
0)l8
0
yU3
pow
nC==cPow){pow.CopyOnWrite()iG
abs;abs
xC
cAbs);abs
c51
pow
cO3
abs.lH3
pow
lP1
pow.SetParamMove(0,abs);tree.lU}
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
eB3
l54
fp_log2(lR)*eL
yS2
cArg:iF3(fp_arg);lC
cConj:iF3(fp_conj);lC
cImag:iF3(fp_imag);lC
cReal:iF3(fp_real);lC
cPolar
l54
fp_polar(nP3
cMod
l54
fp_mod(nP3
cAtan2:{yG3
p0
eM
yR2
p1
eM
1));nS&&fp_equal(lR,xH1){if(p1
cM2
p1
e23
xH1{xO
fp_const_pi
e33
if(p1
iN&&p1.yQ>=yL2
xH1;nD}
if(i62
fp_equal(eL,xH1){if(p0
cM2
p0
e23
xH1{xO-fp_const_pihalf
e33
if(p0
iN&&p0.yQ>xH1{xO
fp_const_pihalf
e33}
if
lI
fp_atan2(lR,eL));nD
if((p1
iN&&p1.yQ>xH1||(p1
cM2
p1
e23
fp_const_negativezero
xK()lC1
yM2;yM2
xC
cPow);yM2
c51
t5
1));yM2.yB
cV1(-1)));yM2.t32
yN2;yN2
tM3
yN2
c51
t5
0));yN2
c51
yM2);yN2
tO
cAtan)eB
0,yN2
t61
1);eX1
cPow:{if(ConstantFolding_PowOperations
cT3
break
lD3
cDiv:nS&&i62
eL!=yL2
lR/eL
yS2
cInv:nS&&lR!=yL2
cV1(1)/lR
yS2
cSub
l54
lR-eL
yS2
cNeg:nS){xO-lR
yS2
cRad:nS){xO
RadiansToDegrees(lR)yS2
cDeg:nS){xO
DegreesToRadians(lR)yS2
cSqr:nS){xO
lR*lR
yS2
cExp2:iF3(fp_exp2);lC
cRSqrt:nS){xO
cV1(1)/fp_sqrt(lR)yS2
cCot:cN2
fp_tan(lZ
cSec:cN2
fp_cos(lZ
cCsc:cN2
fp_sin(lZ
cHypot
l54
fp_hypot(nP3
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
cFCall:yN3
do_return:;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&cR3)<<"]Done ConstantFolding, result: "
cN
DumpHashes
eI2);
#endif
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
t1{eV1
yG3::set_abs(nL
bool
has_negative=!min
xB3||yQ<cV1();bool
has_positive=!max
xB3||max.val>cV1();bool
crosses_axis=has_negative&&has_positive;cC1
xK
newmax;if(min
xB3&&max
xB3)newmax.set(fp_max(i02,i12
yU3
crosses_axis)min.set(cV1());e43
min
xB3&&max
xB3)min.set(fp_min(i02,i12);iH1
min
xB3)min.set(i02);else
min.set(i12;}
max=newmax;}
eV1
yG3::set_neg(){std::swap(min,max);yQ=-yQ;max.val=-max.val;}
xQ1
IsLogicalTrueValue
cR1
yG3&p
lI3{if(nA
IsIntType
xK::xD3){if(p
iN&&p.yQ>=cV1(1
e21
if(!abs&&p
nH2
cO
val<=cV1(-1
e21}
e43
p
iN&&p.yQ>=cV1(0.5
e21
if(!abs&&p
nH2
cO
val<=cV1(-0.5
e21}
xB2
xQ1
IsLogicalFalseValue
cR1
yG3&p
lI3{if(nA
IsIntType
xK::xD3){if(abs)lB3
cO
known
lH1
1);else
lB3
iN&&p
nH2.yQ>cV1(-1)lH1
1);}
e43
abs)lB3
cO
known
lH1
0.5);else
lB3
iN&&p
nH2.yQ>cV1(-0.5)lH1
0.5);}
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lL3
FUNCTIONPARSERTYPES;using
t1;lL3{nR1
T>inline
int
isnan_workaround(T
t)yJ(t!=t);}
}
t1{tJ1
yG3
iT
const
eK
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
{using
lL3
FUNCTIONPARSERTYPES;yG3
tmp=CalculateResultBoundaries_do
eI2
x62"Estimated boundaries: "
;if(tmp
iN)std::cout<<tmp.yQ;else
std::cout<<"-inf"
iW2" .. "
;if(tmp
eU1
std::cout<<tmp
cO
val;else
std::cout<<"+inf"
iW2": "
yH2
tree
x62
std::endl
iF
tmp;}
tJ1
yG3
CodeTree
xK::CalculateResultBoundaries_do
cR1
eK
#endif
{iC
yG1(-fp_const_pihalf
xK(),fp_const_pihalf
xK());iC
pi_limits(-fp_const_pi
xK(),fp_const_pi
xK());iC
abs_pi_limits(y11,fp_const_pi
xK());iC
plusminus1_limits(cV1(-y93
using
lL3
std;switch
eI2
nC
xI3
cImmed:nM
tree
l31,tree
l31)tJ3
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cNot:case
cV3:case
cNotNot:case
cAbsNotNot:case
cEqual:case
i81:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:{nM
y11,e53
lD3
cAbs:lD
m.set_abs(yX1
cLog
i22
fp_log
x52
fp_log
yX1
cLog2
i22
fp_log2
x52
fp_log2
yX1
cLog10
i22
fp_log10
x52
fp_log10
yX1
cAcosh
c71
tK1
cGreaterOrEq
eL2
fp_acosh);m
cO
tK1
cGreaterOrEq
eL2
fp_acosh
yX1
cAsinh
c71
set(fp_asinh);m
cO
set(fp_asinh
yX1
cAtanh
c71
n3-1),fp_atanh);m
cO
tK1
cLess
eL2
fp_atanh
yX1
cAcos:lD
nM(m
cM2
m
e23
e53?fp_acos
cU3):y11,(m
iN&&(m.yQ)>=cV1(-1))?fp_acos(m.yQ):fp_const_pi
xK())lD3
cAsin
c71
n3-1),fp_asin,yG1.yQ);m
cO
tK1
cLess
eL2
fp_asin,yG1
cO
val
yX1
cAtan
c71
set(fp_atan,yG1.yQ);m
cO
set(fp_atan,yG1
cO
val
yX1
cAtan2:{nS&&fp_equal(lR,xH1)yJ
abs_pi_limits;}
if(i62
fp_equal(eL,xH1)yJ
yG1;}
return
pi_limits
lD3
cSin:lD
bool
x41=!m
iN||!m
cO
known||cU3-m.yQ)>=(yK
x41)iD
cV1
min=lC3.yQ,yK
min<xH1
min
yN
cV1
max=lC3
cO
val,yK
max<xH1
max
yN
if(max<min)max
yN
bool
y01=(min<=fp_const_pihalf
xK()&&max>=fp_const_pihalf
xK());bool
nJ1=(min<=cQ&&max>=cQ
yU3
y01&&nJ1)iD
if(nJ1)nM
cV1(-1),nI2
if(y01)nM
yO2
e53;nM
yO2
nI2}
case
cCos:lD
if(m
iN)m.yQ+=fp_const_pihalf
xK(yU3
m
eU1
m
cO
val+=fp_const_pihalf
xK();bool
x41=!m
iN||!m
cO
known||cU3-m.yQ)>=(yK
x41)iD
cV1
min=lC3.yQ,yK
min<xH1
min
yN
cV1
max=lC3
cO
val,yK
max<xH1
max
yN
if(max<min)max
yN
bool
y01=(min<=fp_const_pihalf
xK()&&max>=fp_const_pihalf
xK());bool
nJ1=(min<=cQ&&max>=cQ
yU3
y01&&nJ1)iD
if(nJ1)nM
cV1(-1),nI2
if(y01)nM
yO2
e53;nM
yO2
nI2}
case
cTan:{nM)lD3
cCeil:lD
m
cO
iE
cFloor
tG1
yX1
cTrunc
tG1);m
cO
iE
cInt
tG1);m
cO
iE
cSinh
c71
set(fp_sinh);m
cO
set(fp_sinh
yX1
cTanh
c71
set(fp_tanh,plusminus1_limits.min);m
cO
set(fp_tanh,plusminus1_limits.max
yX1
cCosh:lD
if(m
iN){if(m
eU1{if(m.yQ>=y11&&m
cO
val
iR{m.yQ
xW}
iH1(m.yQ)<y11&&m
cO
val
iR{cV1
tmp
xW
if(tmp>m
cO
val)m
cO
val=tmp;m.yQ=cV1(1);}
else{m.yQ
xW
std::swap(m.yQ,m
cO
val);}
}
e43
m.yQ
iR{m
e8
m.yQ=fp_cosh(m.yQ);}
else{m
e8
m.yQ=cV1(1);}
}
}
else{m
iN=true;m.yQ=cV1(1
yU3
m
eU1{m.yQ=fp_cosh
cU3);m
e8}
else
m
e8}
return
m
lD3
cIf:case
cAbsIf:{yG3
res1
eM
1));yG3
res2
eM
2)yU3!res2
iN)res1
iN=false;iH1
res1
iN&&(res2.yQ)<res1.yQ)res1.yQ=res2.yQ;iH1
isnan_workaround(res2.yQ))res1.yQ=res2.yQ;if(!res2
eU1
res1
e8
iH1
res1
cM2
res2
cO
val)>res1
cO
val)res1
cO
val=res2
cO
val;iH1
isnan_workaround(res2
cO
val))res1
cO
val=res2
cO
val
iF
res1
lD3
cMin:{bool
iJ=false;bool
iK=false;yH3;x4
m
eM
a)yU3!m
iN)iJ=true;yQ1
iN||(m.yQ)<e63)e63=m.yQ;if(!m
eU1
iK=true;yQ1
cO
known||cU3)<nR3)nR3=m
cO
val;}
if(iJ)nS3
iK)nU3
return
nV3
cMax:{bool
iJ=false;bool
iK=false;yH3;x4
m
eM
a)yU3!m
iN)iJ=true;yQ1
iN||m.yQ>e63)e63=m.yQ;if(!m
eU1
iK=true;yQ1
cO
known||m
cO
val>nR3)nR3=m
cO
val;}
if(iJ)nS3
iK)nU3
return
nV3
cAdd:{yH3(y11,xH1;x4
item
eM
a)yU3
item
iN)e63+=item.yQ;else
nS3
item
eU1
nR3+=item
cO
val;else
nU3
if(!nW3&&!xD3
eU1
yN3
if(nW3&&nT3&&e63>nR3)std::swap
i32,nR3)iF
nV3
cMul:{cW2
Value{enum
nY3{i42,iZ1,i52}
;nY3
eQ;cV1
value;Value(nY3
t):eQ(t),value(0){}
Value(cV1
v):eQ(i42),value(v){}
bool
cO2
yU1
eQ==iZ1||(eQ==i42&&value<xH1
nD3
eH1*=cR1
Value&rhs){if(eQ==i42&&rhs.eQ==i42)value*=rhs.value;else
eQ=(cO2)!=rhs.cO2)?iZ1:i52);}
iL2<cR1
Value&rhs
yU1(eQ==iZ1&&rhs.eQ!=iZ1)||(eQ==i42&&(rhs.eQ==i52||(rhs.eQ==i42&&value<rhs.value)));}
}
;cW2
yH1{Value
yP2,yQ2;yH1():yP2(Value::i52),yQ2(Value::iZ1){}
void
xI2
Value
e83,const
Value&value2){e83*=value2;if(e83<yP2)yP2=e83;if(yQ2<e83)yQ2=e83;}
}
;yH3(cV1(y93
x4
item
eM
a)yU3!item
iN&&!item
eU1
nM);Value
nZ3=nW3?Value
i32
lF2
iZ1);Value
x03=nT3?Value(nR3
lF2
i52);Value
x13=item
iN?Value(item.yQ
lF2
iZ1);Value
x23=item
cO
known?Value(item
cO
val
lF2
i52);yH1
range;range.xI2
nZ3,x13
t83
nZ3,x23
t83
x03,x13
t83
x03,x23
yU3
range.yP2.eQ==Value::i42)e63=range.yP2.value;else
nS3
range.yQ2.eQ==Value::i42)nR3=range.yQ2.value;else
nU3
if(!nW3&&!xD3
eU1
yN3
if(nW3&&nT3&&e63>nR3)std::swap
i32,nR3)iF
nV3
cMod:{yG3
x
eM
yR2
y
eM
1)yU3
y
eU1{if(y
cO
val
iR{if(!x
iN||(x.yQ)<xH1
nM-y
cO
val,y
cO
val);else
nM
y11,y
cO
val);}
e43!x
cO
known||(x
cO
val)iR
nM
y
cO
val,-y
cO
val);else
nM
y
cO
val,fp_const_negativezero
xK());}
}
else
nM)lD3
cPow:{if(i62
eL==xH1{nM
cV1(y93}
nS&&lR==xH1{nM
y11,xH1;}
nS&&fp_equal(lR,y63{nM
cV1(y93}
if(i62
eL>y11&&GetEvennessInfo(t5
1))==IsAlways){cV1
e52
eL;yG3
tmp
eM
yR2
xD3;nW3=true;e63=0;if(tmp
iN&&tmp.yQ
iR
e63=tC3
tmp.yQ,i72
iH1
tmp
cO
known&&tmp
cO
val<=xH1
e63=tC3
tmp
cO
val,i72
nU3
if(tmp
iN&&tmp
eU1{nT3=true;nR3=fp_max(fp_abs(tmp.yQ),fp_abs(tmp
cO
val));nR3=tC3
nR3,i72}
return
xD3;}
yG3
p0
eM
yR2
p1
eM
1))nX3
p0_positivity=(p0
iN&&(p0.yQ)iR?IsAlways:(p0
cM2
p0
e23
y11?lE3
Unknown)nX3
cP2=GetEvennessInfo(t5
1))nX3
eX=Unknown;switch(p0_positivity)i91
eX=IsAlways;lC
lE3{eX=cP2
iY2}
yM3
switch(cP2)i91
eX=IsAlways;lC
lE3
lC
Unknown:{if(i62!t12
eL)&&eL
iR{eX=IsAlways;}
yN3}
}
switch(eX)i91{cV1
min=y11;if(p0
iN&&p1
iN){min=tC3
p0.yQ,p1.yQ
yU3
p0.yQ<y11&&(!p1
cO
known||p1
cO
val
iR&&min
iR
min=y11;}
if(p0
iN&&p0.yQ>=y11&&p0
cO
known&&p1
eU1{cV1
max=tC3
p0
cO
val,p1
cO
val
yU3
min>max)std::swap(min,max);nM
min,max);}
nM
min,false)lD3
lE3{nM
false,fp_const_negativezero
xK());}
yM3{yN3
eX1
cNeg:lD
m.set_neg(yX1
cSub:{yO
cNeg);tmp2
eS3
tmp
xC
cAdd);tmp
eU3;tmp
c51
tmp2)iF
lH
cInv:{cQ2-1
yB3
cDiv:{yO
cInv);tmp2
eS3
tmp
xC
xB1
nK2
tmp2)iF
lH
cRad:{nK1
xB1
yB
fp_const_rad_to_deg
xK(yB3
cDeg:{nK1
xB1
yB
fp_const_deg_to_rad
xK(yB3
cSqr:{cQ2
2
yB3
cExp:{nK1
cPow);tmp.yB
fp_const_e
xK()));tmp
eU3
iF
lH
cExp2:{nK1
cPow);tmp.yB
x33
tmp
eU3
iF
lH
cCbrt
c71
set(fp_cbrt);m
cO
set(fp_cbrt
yX1
cSqrt:lD
if(m
iN)m.yQ=(m.yQ)<y11?0:fp_sqrt(m.yQ
yU3
m
eU1
m
cO
val=cU3)<y11?0:fp_sqrt
cU3
yX1
cRSqrt:{cQ2-0.5
yB3
cHypot:{CodeTree
xK
xsqr,ysqr,add,sqrt;xsqr
eU3;xsqr.yB
x33
ysqr
eS3
ysqr.yB
x33
xsqr
xC
cPow);ysqr
xC
cPow);add
c51
xsqr);add
c51
ysqr);add
xC
cAdd);sqrt
c51
add);sqrt
xC
cSqrt)iF
iT
sqrt)lD3
eB3:{yO
cLog2);tmp2
eU3;tmp
tM3
tmp
c51
tmp2);tmp.nJ
1))iF
lH
cCot:{yO
cTan)nY
lH
cSec:{yO
cCos)nY
lH
cCsc:{yO
cSin)nY
iT
tmp);}
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
l53:lC
cArg:case
cConj:case
cImag:case
cReal:case
cPolar:lC
cPCall:lC
cFCall:yN3
nM);}
tJ1
TriTruthValue
GetIntegerInfo
cR1
eK{switch
eI2
nC
xI3
cImmed:return
isInteger
eI2
l31)?IsAlways:eO3
tJ3
cFloor:case
cCeil:case
cTrunc:case
cInt:return
IsAlways
tJ3
cAnd:case
cOr:case
cNot:case
cNotNot:case
cEqual:case
i81:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:return
IsAlways
tJ3
cIf:{TriTruthValue
a=GetIntegerInfo(t5
1))nX3
b=GetIntegerInfo(t5
2)yU3
a==b)return
a
xQ2
case
cAdd:case
cMul:{for
cK
if(GetIntegerInfo(t73)!=IsAlways)return
Unknown
iF
IsAlways;}
yM3
yN3
return
Unknown;}
xQ1
IsLogicalValue
cR1
eK{switch
eI2
nC
xI3
cImmed:return
fp_equal
eI2
l31,xH1||fp_equal
eI2
l31,e53
tJ3
cAnd:case
cOr:case
cNot:case
cNotNot:case
cAbsAnd:case
cAbsOr:case
cV3:case
cAbsNotNot:case
cEqual:case
i81:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:x0
cMul:{for
cK
if(!xL3
a))tX1
nX2
case
cIf:case
cAbsIf:yJ
xL3
1))nV1
t5
2));}
yM3
yN3
xB2}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lL3
FUNCTIONPARSERTYPES;
#if defined(__x86_64) || !defined(FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS)
# define CBRT_IS_SLOW
#endif
#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif
lL3
xK1{extern
const
unsigned
char
powi_table[256];}
lL3{using
t1
n11
bool
IsOptimizableUsingPowi(long
immed,long
penalty=0){xK1
cL3
synth;e73
PushVar(l53);size_t
bytecodesize_backup=e73
GetByteCodeSize();xK1::x21
immed,xK1::i01
xK::MulSequence,synth);size_t
bytecode_grow_amount=e73
GetByteCodeSize()-bytecodesize_backup
iF
bytecode_grow_amount<size_t(MAX_POWI_BYTECODE_LENGTH-penalty);}
eV1
ChangeIntoRootChain(lR1
tree,bool
lF3,long
i82,long
i92){while(i92>0){nK1
cCbrt);nM1
tmp.lH3
tree.iA2--i92;}
while(i82>0){nK1
cSqrt
yU3
lF3){tmp
xC
cRSqrt);lF3=false;}
nM1
tmp.lH3
tree.iA2--i82;}
if(lF3){nK1
cInv);nM1
tree.iA2}
}
tJ1
cW2
RootPowerTable{static
const
cV1
RootPowers[(1+4)*(1+3)];}
n11
const
cV1
t6(1+4)*(1+3)]={cV1(1)lS
iA3
iA3
2*iA3
tH1)lS
3*2)lS
3*2*2)lS
3*tH1*2*tH1*3
yW2
2
yW2
2*2
yW2
tH1*3*2*tH1*3*3
yW2
3*2
yW2
3*2*2
yW2
3*tH1*3*3*2*2*2*2)}
;cW2
PowiResolver{static
const
unsigned
MaxSep=4;static
x53
iB3=5;typedef
int
e93;typedef
long
xC3;typedef
long
tN;cW2
yT2{yT2():n_int_sqrt(0),n_int_cbrt(0),sep_list(),n21(0){}
int
n_int_sqrt;int
n_int_cbrt;int
tW1
MaxSep];tN
n21;}
n11
static
yT2
CreatePowiResult(cV1
eX3){yT2
xD3;e93
t7=FindIntegerFactor(i72
if(t7==0){
#ifdef DEBUG_POWI
iB2"no factor found for %Lg\n"
,(eA3);
#endif
return
xD3;}
iC2
n21=y21
eX3,t7);xC3
eM2=EvaluateFactorCost(t7,0,0,0)+c7
iC2
n21);int
iC3=0;int
iD3=0;int
x43=0;
#ifdef DEBUG_POWI
iB2"orig = %Lg\n"
,(eA3);iB2"plain factor = "
iR3"%ld\n"
,(int)t7,(long)eM2);
#endif
for
l51
n_s=0;n_s<MaxSep;++n_s){int
x8=0;xC3
yI1=eM2;e93
yR1=t7;for(int
s=1;s<iB3*4;++s){
#ifdef CBRT_IS_SLOW
if(s>=iB3)break;
#endif
int
n_sqrt=s%iB3;int
n_cbrt=s/iB3;if(n_sqrt+n_cbrt>4)y91
cV1
lI1=eX3;lI1-=t6
s];iB1=FindIntegerFactor(lI1
yU3
xU2!=0){tN
xR=y21
lI1,xU2);xC3
cost=EvaluateFactorCost(xU2,iC3+n_sqrt,iD3+n_cbrt,x43+1)+c7
xR);
#ifdef DEBUG_POWI
iB2"Candidate sep %u (%d*sqrt %d*cbrt)factor = "
iR3"%ld (for %Lg to %ld)\n"
,s,n_sqrt,n_cbrt,xU2,(long)cost
tJ2
lI1,(long)xR);
#endif
if(cost<yI1){x8=s;yR1=xU2;yI1=cost;}
}
}
if(!x8)break;
#ifdef DEBUG_POWI
iB2"CHOSEN sep %u (%d*sqrt %d*cbrt)factor = "
iR3"%ld, exponent %Lg->%Lg\n"
,x8,x8%iB3,x8/iB3,yR1,yI1
tJ2
eW3)tJ2
eW3-t6
x8]));
#endif
iC2
tW1
n_s]=x8;eX3-=t6
x8];iC3+=x8%iB3;iD3+=x8/iB3;eM2=yI1;t7=yR1;x43+=1;}
iC2
n21=y21
eX3,t7);
#ifdef DEBUG_POWI
iB2"resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n"
,iC2
n21,(eA3
tJ2
t7);
#endif
while(t7%2==0){++iC2
n_int_sqrt;t7/=2;}
while(t7%3==0){++iC2
n_int_cbrt;t7/=3;}
return
xD3;}
private:static
xC3
c7
tN
xR){static
std::map
cR2
i7;if(xR<0){xC3
cost=22
iF
cost+c7-xR);}
std::map
cR2::y53
i=i7.xV2
xR
yU3
i!=i7.cZ1
xR)return
i
e42;std::pair
cR2
xD3(xR,0.0);xC3&cost=iC2
second;while(xR>1){int
xU2=0;if(xR<256){xU2=xK1::powi_table[xR];if(xU2&128)xU2&=127;else
xU2=0;if(xU2&64)xU2=-(xU2&63)-1;}
if(xU2){cost+=c7
xU2);xR/=xU2;y91}
if(!(xR&1)){xR/=2;cost+=6;}
else{cost+=7;xR-=1;}
}
i7.y13,xD3)iF
cost;}
c81
tN
y21
yJ1,iB1)yJ
makeLongInteger(value*cV1(xU2));}
c81
bool
yK1
yJ1,iB1){cV1
v=value*cV1(xU2)iF
isLongInteger(v);}
c81
e93
FindIntegerFactor(yJ1){iB1=(2*2*2*2);
#ifdef CBRT_IS_SLOW
#else
xU2*=(3*3*3);
#endif
e93
xD3=0;if(yK1
value,xU2)){xD3=xU2;while((xU2%2)==0&&yK1
value,xU2/2))xD3=xU2/=2;while((xU2%3)==0&&yK1
value,xU2/3))xD3=xU2/=3;}
#ifdef CBRT_IS_SLOW
if(xD3==0){if(yK1
value,3
tB2
3;}
#endif
return
xD3;}
static
int
EvaluateFactorCost(int
xU2,int
s,int
c,int
nmuls){x53
x63=6;
#ifdef CBRT_IS_SLOW
x53
eN2=25;
#else
x53
eN2=8;
#endif
int
xD3=s*x63+c*eN2;while(xU2%2==0){xU2/=2;xD3+=x63;}
while(xU2%3==0){xU2/=3;xD3+=eN2;}
xD3+=nmuls
iF
xD3;}
}
;}
t1{xQ1
CodeTree
xK::RecreateInversionsAndNegations(bool
prefer_base2){bool
changed=false
xF1
0;a<l21++a)if(lX1.RecreateInversionsAndNegations(prefer_base2))yM1
if(changed){exit_changed:Mark_Incompletely_Hashed()nX2
switch(iU1{case
cMul:{eH
n32
iG
n42,cT1;if(true){bool
nO1=false;cV1
nT2=0
xF1
l21
a
nX
yU2
0)cS2
tF
1)tF1){nO1=true;nT2=tF
1)l31
iY2}
}
if(nO1){cV1
immeds=1.0
xF1
l21
a
nX
tF1){immeds*=powgroup
l31;yL1}
for
tS1
a=l21
a-->0;){lR1
powgroup=lX1;if(powgroup
yU2
0)cS2
tF
1).IsImmed(lC1&log2=tF
0);log2.l81
log2
xC
eB3);log2.yB
tC3
immeds,cV1(1)/nT2)));log2.Rehash(iX2}
}
}
for
tS1
a=l21
a
nX
yU2
1)tF1){x92
exp_param=tF
1);cV1
e52
exp_param
l31;if(e51,cV1(-1))){l81
n32.push_back(lX1
cO3
yL1
iH1
eX3<y11&&t12
eX3
lC1
iL;iL
xC
cPow);iL
tM
tF
0));iL.yB-eX3));iL.lH3
n32.push_back(iL);l81
yL1}
iH1
powgroup
cS2!n42
nJ2)){n42=tF
0);l81
yL1
iH1
powgroup
nC==eB3&&!cT1
nJ2)){cT1=powgroup;l81
yL1}
if(!n32.empty()){changed=true
iG
cH1;cH1
tM3
cH1.SetParamsMove(n32);cH1.t32
y31
cMul);c61
SetParamsMove
eF
if(c61
IsImmed()&&fp_equal
iK2
l31,y63{n52
cInv)eN
cH1);}
e43
c61
GetDepth()>=cH1
nD2){n52
cDiv)eN
xF3)eN
cH1);}
else{n52
cRDiv)eN
cH1
x32}
}
}
if(n42
nJ2
lC1
y31
iU1;c61
SetParamsMove
eF
while(c61
RecreateInversionsAndNegations(prefer_base2))c61
FixIncompleteHashes();n52
eB3)eN
n42
x32
yM1}
if(cT1
nJ2
lC1
y31
cMul);xF3
c51
cT1
l8
1));c61
AddParamsMove
eF
while(c61
RecreateInversionsAndNegations(prefer_base2))c61
FixIncompleteHashes();DelParams();n52
eB3)eN
cT1
l8
0)x32
yM1
eX1
cAdd:{eH
iF2
xF1
cG1
lX1
nC==cMul){n62
y41:iG&xF3
xJ2
iH
for
tS1
b=c61
l21
b-->0;){if
iK2
l8
b)lE1
xU2=xF3
l8
b)l31;if(fp_equal(xU2,iZ2
y41;}
c61
l81
c61
DelParam(b
iD2
iH1
fp_equal(xU2,cV1(-2)))xB
y41;}
c61
l81
c61
DelParam(b);c61
yB
cV1(2))iD2}
}
if(eY){c61
t8
xF3);yL1}
iH1
lX1
nC==cDiv&&!IsIntType
xK::xD3){n62
y51:iG&cH1
xJ2
iH
if(cH1
l8
0)eC3
fp_equal(cH1
l8
0)l31,iZ2
y51;}
cH1.l81
cH1
x02
0);cH1
xC
cInv
iD2}
if(eY)xB
y51;}
cH1.t8
cH1);yL1}
iH1
lX1
nC==cRDiv&&!IsIntType
xK::xD3){n62
xC1:iG&cH1
xJ2
iH
if(cH1
l8
1)eC3
fp_equal(cH1
l8
1)l31,iZ2
xC1;}
cH1.l81
cH1
x02
1);cH1
xC
cInv
iD2}
if(eY)xB
xC1;}
cH1.t8
cH1);yL1}
if(!iF2.empty()){
#ifdef DEBUG_SUBSTITUTIONS
iB2"Will make a Sub conversion in:\n"
);fflush(stdout);iU
#endif
CodeTree
xK
yX2;yX2
xC
cAdd);yX2.SetParamsMove(iF2);yX2.t32
cU1;cU1
xC
cAdd);cU1.SetParamsMove(l22));cU1.lH3
if(cU1
tF1&&fp_equal(cU1
l31,xH1){n52
cNeg);e4);}
e43
cU1
nD2==1){n52
cRSub);e4)eD3}
iH1
yX2
nC==cAdd){n52
cSub)eD3
e4
l8
0))xF1
1;a<yX2.l21++a){CodeTree
xK
eO2;eO2
xC
cSub);eO2.SetParamsMove(l22));eO2.Rehash(false)eN
eO2);e4
nL3}
}
else{n52
cSub)eD3
e4);}
}
#ifdef DEBUG_SUBSTITUTIONS
iB2"After Sub conversion:\n"
);fflush(stdout);iU
#endif
eX1
cPow:{x92
p0
xJ2(0);x92
p1
xJ2(1
yU3
p1
eC3
p1
l31!=y11&&!t12
p1
l31)){eC
yT2
r=eC
CreatePowiResult(fp_abs(p1
l31)yU3
r.n21!=0){bool
l02=false;if(p1
l31<y11&&r.tW1
0]==0&&r
iE2>0){l02=true;}
#ifdef DEBUG_POWI
iB2"Will resolve powi %Lg as powi(chain(%d,%d),%ld)"
tJ2
fp_abs(p1
l31),r
iE2,r.n_int_cbrt,r.n21)yP3
n=0;n<eC
MaxSep;++n){if(r
yC3==0)break;int
n_sqrt=r
yC3%eC
iB3;int
n_cbrt=r
yC3/eC
iB3;iB2"*chain(%d,%d)"
,n_sqrt,n_cbrt);}
iB2"\n"
);
#endif
CodeTree
xK
cT2
xJ2(0)iG
yY2=cT2;yY2.l81
ChangeIntoRootChain(yY2,l02,r
iE2,r.n_int_cbrt);yY2.t32
pow;if(r.n21!=1){pow
xC
cPow);pow
c51
yY2);pow.yB
cV1(r.n21)));}
else
pow.swap(yY2)iG
mul;mul
tM3
mul
c51
pow)yP3
n=0;n<eC
MaxSep;++n){if(r
yC3==0)break;int
n_sqrt=r
yC3%eC
iB3;int
n_cbrt=r
yC3/eC
iB3
iG
eP2=cT2;eP2.l81
ChangeIntoRootChain(eP2,false,n_sqrt,n_cbrt);eP2.lH3
mul
c51
eP2);}
if(p1
l31<y11&&!l02){mul.lH3
n52
cInv
x22
0,mul);DelParam(1);}
else{n52
cMul);SetParamsMove(mul.l22));}
#ifdef DEBUG_POWI
iU
#endif
yM1
yN3}
}
if(GetOpcode()==cPow&&(!p1
tF1||!isLongInteger(p1
l31)||!IsOptimizableUsingPowi
xK(makeLongInteger(p1
l31)))){if(p0
tF1&&p0
l31>cV1(0.0)){if(prefer_base2){cV1
yZ2=fp_log2(p0
l31
yU3
fp_equal(yZ2,y63{DelParam(0);}
else{n0
e91
yZ2));eX3
tM
p1
lD1
SetParamMove(lJ1}
n52
cExp2);yM1}
else{cV1
yZ2=fp_log(p0
l31
yU3
fp_equal(yZ2,y63{DelParam(0);}
else{n0
e91
yZ2));eX3
tM
p1
lD1
SetParamMove(lJ1}
n52
cExp);yM1}
}
iH1
GetPositivityInfo(p0)==IsAlways){if(prefer_base2){CodeTree
xK
log;log
xC
cLog2);log
tM
p0);log.lH3
n0
p1);eX3
c51
log
lD1
n52
cExp2
x22
lJ1
yM1}
else{CodeTree
xK
log;log
xC
cLog);log
tM
p0);log.lH3
n0
p1);eX3
c51
log
lD1
n52
cExp
x22
lJ1
yM1}
}
eX1
cDiv:{if(GetParam(0)tF1&&fp_equal(GetParam(0)l31,y63{n52
cInv);DelParam(0);}
yN3
yM3
yN3
if(changed)goto
exit_changed
iF
changed;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
lL3
FUNCTIONPARSERTYPES;lL3{using
t1;class
tM2{size_t
nP1;size_t
eD;size_t
eE;size_t
lK1;size_t
eZ;size_t
t0;size_t
n91;eQ3
tM2():nP1(0),eD(0),eE(0),lK1(0),eZ(0),t0(0),n91(0){}
void
eE3
OPCODE
op){nP1+=1
tL2
cCos)++eD
tL2
cSin)++eE
tL2
cSec)++eD
tL2
cCsc)++eE
tL2
cTan)++lK1
tL2
cCot)++lK1
tL2
cSinh)++t0
tL2
cCosh)++eZ
tL2
cTanh)++n91;}
size_t
GetCSEscore()const{size_t
xD3=nP1
iF
xD3;}
int
NeedsSinCos()const{bool
y61=(nP1==(eD+eE+lK1)yU3(lK1&&(eE||eD))||(eE&&eD)){if(y61)return
1
iF
2;}
return
0;}
int
NeedsSinhCosh()const{bool
y61=(nP1==(eZ+t0+n91)yU3(n91&&(t0||eZ))||(t0&&eZ)){if(y61)return
1
iF
2;}
return
0;}
size_t
MinimumDepth()const{size_t
n_sincos=std::min(eD,eE);size_t
n_sinhcosh=std::min(eZ,t0
yU3
n_sincos==0&&n_sinhcosh==0)return
2
iF
1;}
}
n11
class
TreeCountType:public
std::multimap<fphash_t,std::pair<tM2,CodeTree
xK> >{}
xP3
FindTreeCounts(tI1&n72,x92
tree,OPCODE
nU2,bool
skip_root=false){cS
i=n72.xV2
tree.GetHash()yU3!skip_root){bool
found=false;for(;i!=n72.cZ1
tree.GetHash();++i){if
eI2
xF
i
e42
t63)){i
e42.first.eE3
nU2);found=true
iY2}
}
if(!found){tM2
count;count.eE3
nU2);n72.y13,std::make_pair
eI2.GetHash(),std::make_pair
tS3,tree)));}
}
for
tS1
a=tL
FindTreeCounts(n72,t73,t52);}
cW2
yS{bool
BalanceGood;bool
FoundChild;}
n11
yS
lL1
x92
root,x92
child){if(root
xF
child)){yS
xD3={true,true}
iF
xD3;}
yS
xD3={true,false}
;if(root
nC==cIf||root
nC==tB3{yS
cond=lL1
root
l8
0
l03
yS
y1=lL1
root
l8
1
l03
yS
y7=lL1
root
l8
2
l03
if(cond
yT||y1
yT||y7
yT){xD3
yT=true;}
xD3
e5=((y1
yT==y7
yT)||y73&&(cond
e5||(y1
yT&&y7
yT))&&(y1
e5||y73&&(y7
e5||y73;}
else{bool
iC1=false;bool
nQ1=false;for
tS1
b=root.GetParamCount(),a=0;a<b;++a){yS
tmp=lL1
root
l8
a
l03
if(tmp
yT)xD3
yT=true;if(tmp
e5==false)iC1=true;iH1
tmp
yT)nQ1=true;}
if(iC1&&!nQ1)xD3
e5=false;}
return
xD3;}
xQ1
lR3
lQ3
iE3
x92
tree,const
xK1
cL3&synth,const
tI1&n72){for
tS1
b=tree.GetParamCount(),a=0;a<b;++a){x92
leaf=t73;cS
synth_it;iS2
tI1::const_iterator
i=n72.y43
i!=n72
l24;++i){if(i->first!=leaf.GetHash())y91
const
tM2&occ=i
e42.first;size_t
score=occ.GetCSEscore();x92
candidate=i
e42
t63;n92
candidate))y91
if(leaf
nD2<occ.MinimumDepth())y91
if(score<2)y91
if(lL1
iE3
leaf)e5==false)continue
nX2
if(lR3(iE3
leaf,synth,n72
e21}
xB2
xQ1
n82
lQ3
y23,x92
expr){yS1
y23
lW3
expr
e21
yS1
n82(y23
l8
a),expr
tB2
true
iF
false;}
xQ1
GoodMomentForCSE
lQ3
y23,x92
expr){if(y23
nC==cIf)return
true;yS1
y23
lW3
expr
e21
size_t
iG2=0;yS1
n82(y23
l8
a),expr))++iG2
iF
iG2!=1;}
}
t1{tJ1
size_t
CodeTree
xK::SynthCommonSubExpressions(xK1
cL3&synth)const{if(e71
0)return
0;size_t
stacktop_before=e73
GetStackTop();tI1
n72;FindTreeCounts(n72,*this,GetOpcode(),true);for(;;){size_t
c02=0;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Finding a CSE candidate, root is:"
<<std::endl;DumpHashes(*this);
#endif
cS
cs_it(n72
l24);for(cS
j=n72.y43
j!=n72
l24;){cS
i(j++);const
tM2&occ=i
e42.first;size_t
score=occ.GetCSEscore();x92
tree=i
e42
t63;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Score "
<<score<<":\n"
<<std::flush;DumpTreeWithIndent
eI2);
#endif
n92
tree))xV
if
eI2
nD2<occ.MinimumDepth())xV
if(score<2)xV
if(lL1*this,tree)e5==false)xV
if(lR3(*this,tree,synth,n72)){y91}
if(!GoodMomentForCSE(*this,tree))xV
score*=tree
nD2;if(score>c02){c02=score;cs_it=i;}
}
if(c02<=0){
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"No more CSE candidates.\n"
<<std::flush;
#endif
yN3
x92
tree=cs_it
e42
t63;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<iY3"Common Subexpression:"
;DumpTree
xK
eI2
x62
std::endl;
#endif
#if 0
int
n31=occ.NeedsSinCos();int
tX=occ.NeedsSinhCosh()iG
iH2,iI2,c12,c22;if(n31){iH2
eQ2
iH2
xC
cSin);iH2.lH3
iI2
eQ2
iI2
xC
cCos);iI2.lH3
n92
iH2)||e73
Find(iI2))l04==2){e31
y91}
n31=0;}
}
if(tX){c12
eQ2
c12
xC
cSinh);c12.lH3
c22
eQ2
c22
xC
cCosh);c22.lH3
n92
c12)||e73
Find(c22)){if(tX==2){e31
y91}
tX=0;}
}
#endif
tree.SynthesizeByteCode(synth,false);e31
#ifdef DEBUG_SUBSTITUTIONS_CSE
e73
xE3
Dump<0>(x62"Done with Common Subexpression:"
;DumpTree
xK
eI2
x62
std::endl;
#endif
#if 0
if(n31)l04==2||tX){e73
eG1}
nG3
cSinCos,1,2)tQ1
iH2,1)tQ1
iI2,0);}
if(tX)l04)e73
eG1
if(tX==2){e73
eG1}
nG3
cSinhCosh,1,2)tQ1
c12,1)tQ1
c22,0);}
#endif
}
return
e73
xI
stacktop_before;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
tJ1
lQ1
xK::iJ2{using
t1;CopyOnWrite()iG
tree;tree.GenerateFrom(*mData);FPoptimizer_Optimize::ApplyGrammars
eI2);std::vector<unsigned>eF3;std::vector
xK
immed;size_t
stacktop_max=0;tree.SynthesizeByteCode(eF3,immed,stacktop_max
yU3
mData->mStackSize!=stacktop_max){mData->mStackSize=unsigned(stacktop_max);
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
mData->mStack.eV3
stacktop_max);
#endif
}
mData->mByteCode.swap(eF3);mData->mImmed.swap(immed);}
#define FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(type) tN1>lQ1<type>::iJ2{}
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
#define FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(type) xE3 lQ1<type>::iJ2;
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
