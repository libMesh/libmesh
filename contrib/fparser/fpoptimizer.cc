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
#define iC3 tree.x5
#define iB3 eB1(-1)
#define iA3 if(n5 tK3
#define i93 tree e6
#define i83 i73"\n";
#define i73 ,o);o<<
#define i63 tree nF
#define i53 :if lJ
#define i43 .empty()
#define i33 "Found "
#define i23 stackpos
#define i13 lS));nG lD
#define i03 "dup(%u) "
#define tZ3 "%d, cost "
#define tY3 "PUSH " eL3
#define tX3 ::cout<<tO3
#define tW3 "immed "<<
#define tV3 mFuncParsers
#define tU3 cS2{assert
#define tT3 stderr
#define tS3 sep2=" "
#define tR3 FPHASH_CONST
#define tQ3 cache_needed[
#define tP3 fprintf
#define tO3 "Applying "
#define tN3 ||tree.GetOpcode
#define tM3 FUNCTIONPARSER_INSTANTIATE_OPTIMIZE
#define tL3 FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE
#define tK3 HANDLE_UNARY_CONST_FUNC
#define tJ3 second;
#define tI3 within,
#define tH3 t0 iK1)
#define tG3 c_count
#define tF3 s_count
#define tE3 2)lT 2*
#define tD3 n61 2))
#define tC3 ;tmp yX
#define tB3 eR3 y01
#define tA3 );m cC2
#define t93 max.val
#define t83 sim.nY1
#define t73 (list t03
#define t63 ].swap(
#define t53 codes[b
#define t43 whydump
#define t33 params
#define t23 b.Value)
#define t13 b.Opcode
#define t03 .first
#define eZ3 .size()
#define eY3 ].second
#define eX3 e6;++b)
#define eW3 Immeds
#define eV3 yK lU 2
#define eU3 ?cNotNot:cAbsNotNot)
#define eT3 nF==cIf?
#define eS3 ,nY3 l9
#define eR3 );m.max.
#define eQ3 },{{1,
#define eP3 AddFrom(
#define eO3 info.SaveMatchedParamIndex(
#define eN3 i93);
#define eM3 NumConstant:
#define eL3 ;DumpTree(
#define eK3 for(;a<
#define eJ3 nparams
#define eI3 0x12 nO
#define eH3 cHypot,
#define eG3 t9 1,0,
#define eF3 nW 0,
#define eE3 cAbs,nW
#define eD3 fp_pow(
#define eC3 false;}
#define eB3 cAbsIf)
#define eA3 yX cond
#define e93 tree yX
#define e83 (exponent
#define e73 FP_GetOpcodeName(
#define e63 =false;
#define e53 (count
#define e43 ]t03
#define e33 Ne_Mask
#define e23 Gt_Mask
#define e13 Lt_Mask
#define e03 {data->
#define cZ3 opcode,
#define cY3 public:
#define cX3 result
#define cW3 cX3 xQ
#define cV3 cX3 xF
#define cU3 cX3 e3
#define cT3 1),eB1(1));
#define cS3 cX3 cU2
#define cR3 cX3))lP2
#define cQ3 cX3(
#define cP3 ++a)if(
#define cO3 pclone
#define cN3 cOr,l6
#define cM3 newpow
#define cL3 change
#define cK3 133,2,
#define cJ3 Params
#define cI3 size_t
#define cH3 Needs
#define cG3 byteCode
#define cF3 child)
#define cE3 lS1 nF==
#define cD3 cLog2by
#define cC3 ,iC,1,iX+1);
#define cB3 ,e62 7168
#define cA3 nF==cLog2&&
#define c93 nF==cPow
#define c83 long tX1
#define c73 factor_t
#define c63 min=x71;
#define c53 value1
#define c43 a));if(!
#define c33 fp_mod(
#define c23 yG known
#define c13 yG n3 0)
#define c03 else{if(
#define yZ3 iJ);}if(
#define yY3 p2 cL ifp2
#define yX3 switch(lT1
#define yW3 xE p2;p2
#define yV3 stackptr
#define yU3 cLog);xV
#define yT3 ==cOr)i22
#define yS3 IsLogicalValue(y1
#define yR3 l9 0));
#define yQ3 opcodes
#define yP3 did_muli
#define yO3 &Value){
#define yN3 yL const
#define yM3 used[b]
#define yL3 :if(&*lD1){
#define yK3 :{lZ1 r=
#define yJ3 ){switch(
#define yI3 lD1=r.specs;if(r.found){
#define yH3 lD1,info
#define yG3 [0].info
#define yF3 sizeof(
#define yE3 l3 16,1,
#define yD3 cLess,cJ
#define yC3 450998,
#define yB3 cExp2,nW
#define yA3 lK 2},0,
#define y93 yX cMul);
#define y83 default:
#define y73 Ge0Lt1
#define y63 Gt0Le1
#define y53 ,y6 0x5 nO
#define y43 default_function_handling
#define y33 for yW1=
#define y23 cAnd,l6
#define y13 cL tree);
#define y03 cAdd l52
#define xZ3 y1 0)nF==
#define xY3 ;goto redo;
#define xX3 lC3 size()
#define xW3 eB1(1)))
#define xV3 iterator
#define xU3 begin();
#define xT3 TreeSet
#define xS3 parent
#define xR3 insert(i
#define xQ3 newrel
#define xP3 b_needed
#define xO3 cachepos
#define xN3 =value
#define xM3 131,4,1,
#define xL3 131,8,1,
#define xK3 4,1,2,1,
#define xJ3 1 y8 n01
#define xI3 eE1 i11
#define xH3 FindPos(
#define xG3 src_pos
#define xF3 reserve(
#define xE3 tree.GetHash()
#define xD3 iW1 tree
#define xC3 eY void
#define xB3 treeptr
#define xA3 .resize(
#define x93 xT1 xE>&
#define x83 tE1 void
#define x73 ImmedTag
#define x63 eY class
#define x53 a,const
#define x43 RefCount
#define x33 Birth();
#define x23 mulgroup
#define x13 unsigned
#define x03 IsAlways
#define nZ3 start_at
#define nY3 leaf1
#define nX3 cost_t
#define nW3 iftree
#define nV3 fpdata
#define nU3 middle
#define nT3 ifdata
#define nS3 synth.
#define nR3 sqrt_cost
#define nQ3 const int
#define nP3 mul_count
#define nO3 maxValue1
#define nN3 minValue1
#define nM3 maxValue0
#define nL3 minValue0
#define nK3 ValueType
#define nJ3 max.known
#define nI3 e0);nG lD
#define nH3 abs_mul
#define nG3 l9 a));
#define nF3 pos_set
#define nE3 p1.lR2 p1
#define nD3 [funcno].
#define nC3 ParamHolder
#define nB3 ){case
#define nA3 eE1[++IP]
#define n93 sim.x3 1,
#define n83 {sim.Eat(
#define n73 eE1[IP]==
#define n63 subtree
#define n53 invtree
#define n43 MakeHash(
#define n33 middle2
#define n23 (constraints&
#define n13 std::string
#define n03 (cI3 n
#define lZ3 yU false;
#define lY3 rulenumit
#define lX3 {tH1 lG
#define lW3 =-tX2 63)-1;
#define lV3 ,l02&iC,
#define lU3 l7 0,2,
#define lT3 l7 0,1,
#define lS3 lB 0x4 nO
#define lR3 cNeg,lU 1
#define lQ3 MakeEqual
#define lP3 n91,l4::
#define lO3 n91,{l4::
#define lN3 newbase
#define lM3 branch1op
#define lL3 branch2op
#define lK3 fp_equal(
#define lJ3 if(lK3
#define lI3 ContainsOtherCandidates
#define lH3 l9 a)xI
#define lG3 overlap
#define lF3 truth_b
#define lE3 truth_a
#define lD3 found_dup
#define lC3 cJ3.
#define lB3 cQ1 xE&
#define lA3 );nZ l4
#define l93 nQ r;r yX
#define l83 rangeutil
#define l73 Plan_Has(
#define l63 ))break;eB1
#define l53 StackMax)
#define l43 eX true;}
#define l33 const x32
#define l23 xB2 hash1
#define l13 namespace
#define l03 ::res,b8<
#define iZ2 inverted
#define iY2 IsNever:
#define iX2 .known&&
#define iW2 );}eY bool
#define iV2 const std::t4
#define iU2 const char*
#define iT2 ,const eB&
#define iS2 param.
#define iR2 &param=*
#define iQ2 }switch(
#define iP2 depcodes
#define iO2 explicit
#define iN2 ,tree x8
#define iM2 cPow,l2 0,2,
#define iL2 cPow,xO1
#define iK2 cCosh,nW
#define iJ2 VarBegin
#define iI2 ].data);
#define iH2 e0)));nZ
#define iG2 };enum
#define iF2 iE2 IsNever xF2
#define iE2 ))return
#define iD2 ;}void
#define iC2 .what nN1
#define iB2 nU){eB1 tmp=
#define iA2 .min.val
#define i92 i82 iA2
#define i82 )cX3
#define i72 )iD2
#define i62 t82 switch(
#define i52 eZ3;++
#define i42 cN1.SubTrees
#define i32 cN1.Others
#define i22 ?0:1))l8
#define i12 begin(),
#define i02 cond_add
#define tZ2 cond_mul
#define tY2 cond_and
#define tX2 (half&
#define tW2 bool eD1
#define tV2 Optimize()
#define tU2 costree
#define tT2 sintree
#define tS2 leaf_count
#define tR2 ;tG=!tG;}}
#define tQ2 sub_params
#define tP2 printf(
#define tO2 cbrt_count
#define tN2 sqrt_count
#define tM2 ,(long double)
#define tL2 exponent);
#define tK2 PlusInf
#define tJ2 Finite
#define tI2 p1 cL ifp1
#define tH2 xL 2,cAdd)
#define tG2 pcall_tree
#define tF2 after_powi
#define tE2 .n8 nS3
#define tD2 GetHash().
#define tC2 t33)
#define tB2 grammar
#define tA2 0x12},{{3,
#define t92 ),0},{
#define t82 break;}
#define t72 data;data.
#define t62 MakeNEqual
#define t52 Dump(std::
#define t42 ;tmp y41 tmp2)
#define t32 isInteger(
#define t22 Comparison
#define t12 needs_flip
#define t02 l72 apos==
#define eZ2 .Remember(
#define eY2 value]
#define eX2 for(typename
#define eW2 ::ByteCodeSynth yB
#define eV2 std::vector<bool>
#define eU2 )eX m t11
#define eT2 y33 0;a<
#define eS2 );p2.lR2 p2 yH yZ
#define eR2 ~cI3(0)
#define eQ2 xB1 xS+1);
#define eP2 Rule&rule,
#define eO2 SetOpcode(
#define eN2 mul_item
#define eM2 innersub
#define eL2 cbrt_cost
#define eK2 best_cost
#define eJ2 2)lT 3*3*
#define eI2 >(eB1(1),
#define eH2 condition
#define eG2 }}return
#define eF2 per_item
#define eE2 item_type
#define eD2 first2
#define eC2 l3 18,1,
#define eB2 ,cLog,nW
#define eA2 info.lR[b].
#define e92 (iS2 data
#define e82 ,cEqual
#define e72 cIf,t9 3,
#define e62 l0 2,
#define e52 lK 1},0,
#define e42 tB 1},0,
#define e32 Decision
#define e22 not_tree
#define e12 Become(y1
#define e02 group_by
#define cZ2 exponent=
#define cY2 ->second
#define cX2 else{x9=new
#define cW2 {data xC lO
#define cV2 cT2 t82}
#define cU2 .min.known
#define cT2 DelParam(a);
#define cS2 xE&tree)
#define cR2 targetpos
#define cQ2 ParamSpec
#define cP2 rhs.hash2;}
#define cO2 rhs.hash1
#define cN2 struct
#define cM2 Forget()
#define cL2 &&cond eQ))
#define cK2 source_tree
#define cJ2 <tM,nX3>
#define cI2 (std::move(
#define cH2 );tmp.x0 0))
#define cG2 );tmp2.x0
#define cF2 (p0 iA2
#define cE2 p1_evenness
#define cD2 isNegative(
#define cC2 .max.n3 0),
#define cB2 ,std::cout)
#define cA2 neg_set
#define c92 0.5))xY1;lD
#define c82 );sim.x3 2,
#define c72 cNop,cNop}}
#define c62 cTanh,cNop,
#define c52 >cN2 cE<
#define c42 matches
#define c32 .match_tree
#define c22 (rule,tree,info
#define c12 cQ1 void*)&
#define c02 l5 0,1,
#define yZ2 cTan,nW
#define yY2 cCos,nW
#define yX2 e1[c i0
#define yW2 +=1 eX n51;
#define yV2 negated
#define yU2 Specializer
#define yT2 coshtree
#define yS2 sinhtree
#define yR2 best_score
#define yQ2 mulvalue
#define yP2 pow_item
#define yO2 subgroup
#define yN2 PowiResult
#define yM2 maxValue
#define yL2 minValue
#define yK2 fp_min(yE,
#define yJ2 )&*nZ3;
#define yI2 range<eB1 nX
#define yH2 (cX3 iA2
#define yG2 ;range yB
#define yF2 fp_max(yE)
#define yE2 div_tree
#define yD2 pow_tree
#define yC2 eB1(0.0)){xB
#define yB2 preserve
#define yA2 PullResult()
#define y92 dup_or_fetch
#define y82 nominator]
#define y72 Rehash(false
#define y62 test_order
#define y52 parampair
#define y42 y52,
#define y32 .param_count
#define y22 shift(index)
#define y12 {std::cout<<
#define y02 rulenumber
#define xZ2 cLessOrEq,cJ
#define xY2 cTan l3 2,1,
#define xX2 cLog l3 2,1,
#define xW2 cEqual,e62
#define xV2 cTanh,nW
#define xU2 cSinh,nW
#define xT2 cInv,lU 1,
#define xS2 constraints=
#define xR2 factor_immed
#define xQ2 changes
#define xP2 n41 cL y5 l9
#define xO2 cL leaf2 l9
#define xN2 cL nY3 l9
#define xM2 cL cond l9
#define xL2 exp_diff
#define xK2 ExponentInfo
#define xJ2 lower_bound(
#define xI2 factor
#define xH2 is_logical
#define xG2 newrel_and
#define xF2 eX Unknown;}
#define xE2 res_stackpos
#define xD2 half_pos
#define xC2 )const yU
#define xB2 rhs xC2
#define xA2 eQ1;if(&*nZ3){x9=(
#define x92 >>1)):(
#define x82 CodeTreeData
#define x72 x23)
#define x62 multiply(
#define x52 var_trees
#define x42 nE OPCODE
#define x32 CodeTree&
#define x22 parent_opcode
#define x12 {eB nZ3;
#define x02 .UseGetNeeded(
#define nZ2 nS3 AddOperation(
#define nY2 ,t5,synth);
#define nX2 =GetParam(
#define nW2 ,cPow,l2 2,2,
#define nV2 changed=true;
#define nU2 log2_exponent
#define nT2 yH swap(tmp);
#define nS2 lS,e0));nG lD
#define nR2 cAbsNot
#define nQ2 dup_fetch_pos
#define nP2 IsNever cI lD
#define nO2 cSin,nW
#define nN2 Value_EvenInt
#define nM2 MakeFalse,{l4
#define nL2 AddCollection
#define nK2 ConditionType
#define nJ2 DUP_ONE(apos)
#define nI2 (x13
#define nH2 cW|nI2)
#define nG2 SpecialOpcode
#define nF2 TreeCountItem
#define nE2 ;if(op==
#define nD2 l01 yI++a){
#define nC2 assimilated
#define nB2 .SetParamsMove(
#define nA2 );nS3 nQ1
#define n92 denominator
#define n82 fraction
#define n72 std::vector<xE>
#define n62 .GetDepth()
#define n52 xI leaf2 l9
#define n42 DUP_BOTH();
#define n32 template lM
#define n22 -1-offset].
#define n12 IsDescendantOf
#define n02 TreeCounts
#define lZ2 bool tG e63
#define lY2 found_log2
#define lX2 div_params
#define lW2 immed_sum
#define lV2 :sim.Eat(1,
#define lU2 ;sim.Push(
#define lT2 minimum_need
#define lS2 .Rehash();
#define lR2 Rehash tA
#define lQ2 OPCODE(opcode)
#define lP2 break;cX3*=
#define lO2 FactorStack yB
#define lN2 x03 cI lD
#define lM2 282870 xG
#define lL2 cNotNot,nW
#define lK2 cNot,nW
#define lJ2 replacing_slot
#define lI2 RefParams
#define lH2 if_always[
#define lG2 WhatDoWhenCase
#define lF2 exponent_immed
#define lE2 base_immed
#define lD2 .IsDefined()
#define lC2 {case x03:
#define lB2 TopLevel)
#define lA2 )l43
#define l92 .IsImmed()
#define l82 1)l92&&
#define l72 else if(
#define l62 nS3 Find(
#define l52 ||op1==
#define l42 data[a eY3
#define l32 iE2 x03;if(
#define l22 lR2 r);}
#define l12 if(newrel_or==
#define l02 PowiCache
#define iZ1 eH 2,131,
#define iY1 Immed eZ3);
#define iX1 eE1.push_back(
#define iW1 const xE&
#define iV1 OptimizedUsing
#define iU1 Var_or_Funcno
#define iT1 iU1;
#define iS1 crc32_t
#define iR1 signed_chain
#define iQ1 MinusInf
#define iP1 n_immeds
#define iO1 stack eZ3
#define iN1 std::cout<<"POP "
#define iM1 FindClone(xO
#define iL1 ,cPow l3
#define iK1 divgroup
#define iJ1 GetOpcode())
#define iI1 needs_rehash
#define iH1 AnyWhere_Rec
#define iG1 ~x13(0)
#define iF1 41,42,43,44,
#define iE1 yH DelParam(
#define iD1 p1_logical_b
#define iC1 p0_logical_b
#define iB1 p1_logical_a
#define iA1 p0_logical_a
#define i91 yX i63);
#define i81 func(val);nM1
#define i71 *const func)
#define i61 nS3 DoDup(
#define i51 cache_needed
#define i41 {if(GetOpcode()
#define i31 2*2*2)lT 3
#define i21 eH 2,1,eH 2,
#define i11 [nT3.ofs+
#define i01 treelist
#define tZ1 has_bad_balance
#define tY1 c73 xI2
#define tX1 double)exponent
#define tW1 fp_abs(t93))
#define tV1 fp_abs(min.val)
#define tU1 cNEqual
#define tT1 tB 2},0,0x0},{{
#define tS1 Oneness_NotOne|
#define tR1 Value_IsInteger
#define tQ1 Constness_Const
#define tP1 DumpHashesFrom(
#define tO1 iV1(
#define tN1 reltype
#define tM1 SequenceOpcodes
#define tL1 goto fail;}
#define tK1 CollectionSet yB
#define tJ1 .GetImmed()
#define tI1 l1 0x4 nO
#define tH1 template<
#define tG1 ,l5 2,1,
#define tF1 GetParams()
#define tE1 lI2);
#define tD1 nS3 PushImmed(
#define tC1 );nS3 StackTopIs(
#define tB1 n02.erase(cs_it);
#define tA1 iE2 true
#define t91 yG set(fp_floor)
#define t81 ,l1 0x0},{{3,
#define t71 ;pow yX cPow);pow
#define t61 ;for(eQ1=0;a<y7;++a)
#define t51 ,cGreater
#define t41 yG val
#define t31 =comp.AddItem(atree
#define t21 ;flipped=!flipped;}
#define t11 ;}case
#define t01 <<tree.tD2
#define eZ1 TreeCountType yB
#define eY1 Value(Value::
#define eX1 y1 a)tJ1
#define eW1 y1 a)l92)
#define eV1 lQ2);
#define eU1 stack[iO1-
#define eT1 stack.push_back(
#define eS1 MaxChildDepth
#define eR1 (*x9)[a].info
#define eQ1 x13 a
#define eP1 std::pair<It,It>
#define eO1 Sign_Negative
#define eN1 Value_Logical
#define eM1 new_factor_immed
#define eL1 e6;a-->0;)if(
#define eK1 occurance_pos
#define eJ1 exponent_hash
#define eI1 exponent_list
#define eH1 CollectMulGroup(
#define eG1 source_set
#define eF1 exponent,xT3
#define eE1 ByteCode
#define eD1 operator
#define eC1 AddParamMove(
#define eB1 Value_t
#define eA1 FindAndDup(tree);
#define e91 back().thenbranch
#define e81 ParamSpec_Extract
#define e71 retry_anyparams_3
#define e61 retry_anyparams_2
#define e51 )y52.tJ3
#define e41 nS3 xH 1
#define e31 ,iC nY2
#define e21 eG(),std::vector<
#define e11 needlist_cached_t
#define e01 grammar_rules[*r]
#define cZ1 yA3 0x4 eQ3
#define cY1 e52 0x4 eQ3
#define cX1 CodeTreeImmed yB(
#define cW1 ;NewHash.hash2+=
#define cV1 for(cI3 b=0;b<
#define cU1 by_float_exponent
#define cT1 ;exponent
#define cS1 lK3 exponent
#define cR1 new_exp
#define cQ1 (const
#define cP1 cQ1 eB1&
#define cO1 l13 FPoptimizer_Optimize
#define cN1 NeedList
#define cM1 xT1 x13>&eE1,cI3&IP,cI3 limit,cI3 y3
#define cL1 end()&&i->first==
#define cK1 return BecomeZero;
#define cJ1 return BecomeOne;
#define cI1 if(lR eZ3<=n1)
#define cH1 addgroup
#define cG1 found_log2by
#define cF1 ParsePowiMuli(
#define cE1 iU1)
#define cD1 yF 529654 xG
#define cC1 new_base_immed
#define cB1 branch1_backup
#define cA1 branch2_backup
#define c91 exponent_map
#define c81 plain_set
#define c71 LightWeight(
#define c61 if(value
#define c51 .sep_list[
#define c41 nC cPow)n7 eB1
#define c31 should_regenerate=true;
#define c21 should_regenerate,
#define c11 Collection
#define c01 RelationshipResult
#define yZ1 Subdivide_Combine(
#define yY1 long value
#define yX1 ;n41 cM op1 yH DelParams
#define yW1 (cI3 a
#define yV1 for yW1 yC
#define yU1 best_sep_factor
#define yT1 l72!cX3
#define yS1 &&p xF<eB1(
#define yR1 tree l92 cI
#define yQ1 GetParamCount()
#define yP1 needlist_cached
#define yO1 inline x13
#define yN1 252421 xG 24830
#define yM1 cZ3 bool pad
#define yL1 cT2}
#define yK1 MakesInteger(
#define yJ1 const eB1&value
#define yI1 best_sep_cost
#define yH1 MultiplicationRange
#define yG1 pihalf_limits
#define yF1 n_stacked
#define yE1 NewHash.hash1
#define yD1 AnyParams_Rec
#define yC1 continue;
#define yB1 Become(value l9 0))
#define yA1 PositionalParams,0}
#define y91 always_sincostan
#define y81 Recheck_RefCount_Div
#define y71 Recheck_RefCount_Mul
#define y61 yB())xL 2,cMul);lD
#define y51 x23.
#define y41 .eC1
#define y31 x23;x23 yX
#define y21 MultiplyAndMakeLong(
#define y11 covers_plus1
#define y01 template set_if<
#define xZ1 }break t11
#define xY1 xL 2,cPow)
#define xX1 if(nS3 FindAndDup(
#define xW1 SynthesizeParam(
#define xV1 ByteCodeSynth yB&synth)
#define xU1 cU2&&
#define xT1 const std::vector<
#define xS1 xR1 std::endl;DumpHashes(
#define xR1 ;std::cout<<
#define xQ1 grammar_func
#define xP1 252180 xG 281854
#define xO1 l2 0,2,165888 xG
#define xN1 cCos l3 2,1,
#define xM1 (p1 tJ1
#define xL1 public eG,public std::vector<
#define xK1 l1 eI3
#define xJ1 Modulo_Radians},
#define xI1 yH SetParam(
#define xH1 PositionType
#define xG1 CollectionResult
#define xF1 const_offset
#define xE1 inline TriTruthValue
#define xD1 stacktop_desired
#define xC1 int mStackPtr=0;
#define xB1 SetStackTop(
#define xA1 }inline
#define x91 FPoptimizer_ByteCode
#define x81 1)?(poly^(
#define x71 eB1(0)
#define x61 },{l4::MakeNotP1,l4::
#define x51 },{l4::MakeNotP0,l4::
#define x41 x71)
#define x31 i63==
#define x21 cond_type
#define x11 fphash_value_t
#define x01 Recheck_RefCount_RDiv
#define nZ1 tmp y41 tree
#define nY1 SwapLastTwoInStack();
#define nX1 SetParams(tF1
#define nW1 fPExponentIsTooLarge(
#define nV1 CollectMulGroup_Item(
#define nU1 pair<eB1,xT3>
#define nT1 nN xB1 xS-1);
#define nS1 covers_full_cycle
#define nR1 AssembleSequence(
#define nQ1 DoDup(found[data.
#define nP1 <<std::dec<<")";}
#define nO1 :return p iA2
#define nN1 !=xK)if(TestCase(
#define nM1 else*this=model;}
#define nL1 std::pair<T1,T2>&
#define nK1 tH1 typename
#define nJ1 has_good_balance_found
#define nI1 n_occurrences
#define nH1 found_log2_on_exponent
#define nG1 covers_minus1
#define nF1 needs_resynth
#define nE1 immed_product
#define nD1 i62 bitmask&
#define nC1 518 xG 400412
#define nB1 Sign_Positive
#define nA1 {DataP slot_holder(y2[
#define n91 ::MakeTrue
#define n81 (nY3 l9 1)n52
#define n71 SetParamMove(
#define n61 CodeTreeImmed(eB1(
#define n51 Suboptimal
#define n41 changed_if
#define n31 n_as_tanh_param
#define n21 opposite=
#define n11 x11(
#define n01 eE1 eZ3
#define lZ1 MatchResultType
#define lY1 needs_sincos
#define lX1 resulting_exponent
#define lW1 Unknown:y83;}
#define lV1 matched_params
#define lU1 ,cIf,l0 3,
#define lT1 GetLogicalValue(y1
#define lS1 GetParam(a)
#define lR1 inverse_nominator]
#define lQ1 cSin l3 2,1,
#define lP1 },{l4::MakeNotNotP1,l4::
#define lO1 },{l4::MakeNotNotP0,l4::
#define lN1 IsImmed()){eB1
#define lM1 AddFunctionOpcode(
#define lL1 void FunctionParserBase
#define lK1 o<<"("<<std::hex<<data.
#define lJ1 IfBalanceGood(
#define lI1 n_as_tan_param
#define lH1 changed_exponent
#define lG1 inverse_denominator
#define lF1 yB(rule.repl_param_list,
#define lE1 retry_positionalparams_2
#define lD1 (*x9)[a].nZ3
#define lC1 x13 index
#define lB1 situation_flags&
#define lA1 ,cLessOrEq yF
#define l91 7168 xG 401798
#define l81 data.subfunc_opcode
#define l71 ,typename xE::
#define l61 PlanNtimesCache(
#define l51 FPoptimizer_Grammar
#define l41 static inline xE
#define l31 yX cLog yH eO2 cMul
#define l21 GetPositivityInfo(tree)!=
#define l11 y33 yQ1;a
#define l01 ;eT2
#define iZ );n71 0,tL2 DelParam(1);
#define iY CopyOnWrite();
#define iX recursioncount
#define iW ParamSpec_SubFunctionData
#define iV (p0 xU1 p0 iA2>=eB1(0.0))
#define iU yW1=yI a-->0;)
#define iT PositionalParams_Rec
#define iS DumpTreeWithIndent(*this);
#define iR switch(type nB3 cond_or:
#define iQ CalculateResultBoundaries(
#define iP tH1 x13 Compare>
#define iO yA3 0x0 eQ3
#define iN edited_powgroup
#define iM has_unknown_max
#define iL has_unknown_min
#define iK static const range yB
#define iJ tree.DelParam(a
#define iI if(keep_powi
#define iH synthed_tree
#define iG 408964 xG 24963
#define iF 528504 xG 24713
#define iE SelectedParams,0},0,0x0},{{
#define iD collections
#define iC cache
#define iB goto ReplaceTreeWithOne;case
#define iA ,2,1);nS3 xU if(found[data.
#define i9 AddOperation(cInv,1,1);nS3 xU}
#define i8 xD3,std::ostream&o
#define i7 ]);nZ2
#define i6 !=xK)return lH2
#define i5 cU1.data
#define i4 iO2 x82(
#define i3 needs_sinhcosh
#define i2 t92 eB1(
#define i1 MakeFalse,l4::
#define i0 ].relationship
#define tZ ,eE1,IP,limit,y3,stack);
#define tY 124024 xG 139399
#define tX 142456 xG 141449
#define tW AnyParams,0}},{ReplaceParams,
#define tV [n1 e43=true;lR[n1 eY3
#define tU l51::Grammar*
#define tT powgroup l9
#define tS eR2&&found[data.
#define tR }},{ProduceNewTree,2,1,
#define tQ n61(
#define tP has_mulgroups_remaining
#define tO by_exponent
#define tN Rehash();tQ2.push_back(
#define tM int_exponent_t
#define tL RootPowerTable yB::RootPowers[
#define tK MatchPositionSpec_AnyParams yB
#define tJ l13 FPoptimizer_CodeTree
#define tI n_as_sinh_param
#define tH n_as_cosh_param
#define tG is_signed
#define tF tF1);y51 Rehash();
#define tE result_positivity
#define tD biggest_minimum
#define tC (tree))goto redo;
#define tB yK AnyParams,
#define tA (yH eC1
#define t9 lB 0x4},{{
#define t8 cond_tree
#define t7 else_tree
#define t6 then_tree
#define t5 sequencing
#define t4 string e73
#define t3 ,l7 2,1,
#define t2 best_factor
#define t1 nR eB1(-cT3
#define t0 ;eC1
#define eZ const iW
#define eY template lY
#define eX ;return
#define eW if_stack
#define eV .max.set(fp_ceil eU2
#define eU n_as_sin_param
#define eT n_as_cos_param
#define eS PowiResolver::
#define eR cIf l3 0,1,
#define eQ .BalanceGood
#define eP eC1 yO2
#define eO {if(needs_cow){iY goto
#define eN );bool needs_cow=GetRefCount()>1;
#define eM valueType
#define eL back().endif_location
#define eK x11 key
#define eJ yW1=0;a<yI cP3 remaining[a])
#define eI eC1 mul);
#define eH 130,1,
#define eG MatchPositionSpecBase
#define eF iO2 CodeTree(
#define eE smallest_maximum
#define eD }PACKED_GRAMMAR_ATTRIBUTE;
#define eC factor_needs_rehashing
#define eB MatchPositionSpecBaseP
#define eA typename eZ1::xV3
#define e9 ,cGreaterOrEq
#define e8 fp_cosh(t41);m xF=fp_cosh(m xF);
#define e7 goto ReplaceTreeWithParam0;
#define e6 .yQ1
#define e5 e81 yB(nS.param_list,
#define e4 e3&&p0 xF<=fp_const_negativezero yB())
#define e3 .nJ3
#define e2 =iQ y1
#define e1 relationships
#define e0 y1 1)tJ1
#define cZ 243,244,245,246,249,250,251,253,255,256,257,258,259}};}
#define cY ];};extern"C"{
#define cX i8=std::cout
#define cW );iX1 0x80000000u
#define cV for(lY3 r=range t03;r!=range.tJ3++r){
#define cU 79,122,123,160,161,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,
#define cT 27,28,29,30,31,32,33,35,36,
#define cS const ParamSpec_SubFunction
#define cR const ParamSpec_ParamHolder
#define cQ }if t73 tJ1==eB1(
#define cP otherhalf
#define cO :{AdoptChildrenWithSameOpcode(tree);
#define cN :goto ReplaceTreeWithZero;case
#define cM .Rehash(yH eO2
#define cL .AddParam(
#define cK StackState
#define cJ l2 16,2,
#define cI )return false;
#define cH nJ3 e63
#define cG eB1(1.5)*fp_const_pi yB()
#define cF CalculatePowiFactorCost(
#define cE ImmedHashGenerator
#define cD ,y6 0x4 nO
#define cC ::map<fphash_t,std::set<n13> >
#define cB eC1 comp.c81[a].value);
#define cA T1,typename T2>inline tW2()(
#define c9 has_nonlogical_values
#define c8 from_logical_context)
#define c7 AnyParams,0}},{ProduceNewTree,
#define c6 y33 y0 e6;a-->0;)
#define c5 POWI_CACHE_SIZE
#define c4 ++IP;yC1}if(n73 yQ3.
#define c3 },{l4::xK,l4::Never},{l4::xK,l4::Never}}
#define c2 .FoundChild
#define c1 BalanceResultType
#define c0 x43(0),Opcode(
#define yZ eO2 nW3 nF)xY3}
#define yY !=eR2){nS3 nQ1
#define yX .eO2
#define yW nZ2 GetOpcode(),
#define yV );void lM1 x13 cZ3 yU2<
#define yU {return
#define yT const yU data->
#define yS +=fp_const_twopi yB();
#define yR eT2 yQ1;++a){if(
#define yQ MatchPositionSpec_AnyWhere
#define yP if e92.match_type==
#define yO ComparisonSetBase::
#define yN eL3 tree)xR1"\n";
#define yM void OutFloatHex(std::ostream&o,
#define yL {static void n43 nE fphash_t&NewHash,
#define yK ,cAdd,
#define yJ paramholder_matches
#define yI i93;
#define yH );tree.
#define yG m.min.
#define yF ,l2 18,2,
#define yE fp_sin(min),fp_sin(max))
#define yD fp_const_twopi yB());if(
#define yC =0;a<xS3 e6;cP3
#define yB <eB1>
#define yA MatchPositionSpec_PositionalParams yB
#define y9 AssembleSequence_Subdivide(
#define y8 ]=0x80000000u|x13(
#define y7 nS y32
#define y6 cPow,lB
#define y5 branch2
#define y4 x13 c;x13 short l[
#define y3 factor_stack_base
#define y2 data->cJ3
#define y1 tree l9
#define y0 branch1
#define xZ MatchInfo yB&
#define xY 0))tC3 cInv)t42
#define xX const SequenceOpCode yB
#define xW ,iB3))eO
#define xV sim.AddConst(
#define xU StackTopIs(*this)eX;}
#define xT {n02.erase(i);yC1}
#define xS StackTop
#define xR FPOPT_autoptr
#define xQ +=cX3 eX cX3;}eY inline eB1
#define xP int_exponent
#define xO newnode
#define xN e52 0x0},{{
#define xM has_highlevel_opcodes
#define xL ;sim.Eat(
#define xK Unchanged
#define xJ best_selected_sep
#define xI .IsIdenticalTo(
#define xH GetStackTop()-
#define xG ,{2,
#define xF .t93
#define xE CodeTree yB
#define xD cAnd,tW
#define xC ->Recalculate_Hash_NoRecursion();}
#define xB tree.ReplaceWithImmed(
#define xA yQ1;cP3 ApplyGrammar(tB2,y1 a),
#define x9 position
#define x8 )){tree.FixIncompleteHashes();}
#define x7 eT2 yI++a){range yB
#define x6 std::vector<CodeTree>
#define x5 SetParam(0,nW3 yR3 xE p1;p1 yX
#define x4 TestImmedConstraints(iS2 constraints,tree)cI
#define x3 SwapLastTwoInStack()xL
#define x2 n83 1,cInv);t82 xV-1)xY1;lD
#define x1 paramholder_index
#define x0 AddParam(y1
#define nZ return true;case
#define nY occurance_counts
#define nX >p e2 a));if(p.
#define nW l0 1,
#define nV -->0;){iW1 powgroup=lS1;if(powgroup
#define nU if(y1 0)l92
#define nT const FPoptimizer_CodeTree::xE&tree
#define nS model_tree
#define nR return range yB(
#define nQ ){xE
#define nP ),rangehalf yB model=rangehalf yB()){if(known
#define nO },{{2,
#define nN ){using l13 FUNCTIONPARSERTYPES;
#define nM AnyParams,1},0,0x0},{{
#define nL n72&lI2
#define nK ConstantFolding_LogicCommon(tree,yO
#define nJ nK1 Ref>inline void xR<Ref>::
#define nI cOr,tW 16,1,
#define nH ):data(new x82 yB(
#define nG goto do_return;}
#define nF .GetOpcode()
#define nE FUNCTIONPARSERTYPES::
#define nD x82 yB::x82(
#define nC xE tmp tC3
#define nB b;}};tH1>cN2 Comp<nE
#define nA xE tmp,tmp2;tmp2 yX
#define n9 iU1(),cJ3(),Hash(),Depth(1),tO1 0){}
#define n8 SynthesizeByteCode(synth);
#define n7 ;tmp.x0 0));tmp cL CodeTreeImmed(
#define n6 while(ApplyGrammar(c12
#define n5 GetIntegerInfo(y1 0))==x03)e7
#define n4 ;tree y41 n41 lA2
#define n3 y01 cGreater>(eB1(
#define n2 DumpParams yB e92.param_list,iS2 data y32,o);
#define n1 restholder_index
#define n0 (lS);if(fp_nequal(tmp,x41){xB eB1(1)/tmp);nG}lD
#define lZ :if(ParamComparer yB()(cJ3[1],cJ3[0])){std::swap(cJ3[0],cJ3[1]);Opcode=
#define lY <typename eB1>
#define lX xE exponent cT1 yX cMul)cT1 cL
#define lW tQ1,0x0},
#define lV eC1 pow l9 1));pow.DelParam(1);pow.Rehash(yH n71 0,pow);goto NowWeAreMulGroup;}
#define lU GroupFunction,0},lW{{
#define lT ,eB1(1)/eB1(
#define lS y1 0)tJ1
#define lR restholder_matches
#define lQ yE1|=key;x11 crc=(key>>10)|(key<<(64-10))cW1((~n11 crc))*3)^1234567;}};
#define lP n41;n41 i91 n41 y41 y1 0));n41 cL y0 l9
#define lO eY xE::CodeTree(
#define lN tree.SetParam(0,y1 0)l9 0)xI1 1,CodeTreeImmed(
#define lM lY void ByteCodeSynth yB::lM1 x13 cZ3 yU2<
#define lL cMul,lU 2,
#define lK cMul,AnyParams,
#define lJ (y1 0)l92&&y1 1)l92){xB
#define lI :cL3=comp.AddRelationship(atree l9 0),atree l9 1),yO
#define lH cPow,l0 2
#define lG typename eB1>inline tW2()cP1 x53 eB1&b)yU a
#define lF eX iQ tmp)t11
#define lE {range yB m e2 0));
#define lD break;case
#define lC xC3 xE::
#define lB yA1,0,
#define lA l1 0x0 nO
#define l9 .GetParam(
#define l8 ;xE n41;n41 i91 n41 nB2 tree.tF1);n41 cM
#define l7 cAdd,tW
#define l6 SelectedParams,0},0,0x0 nO
#define l5 lK 0}},{ReplaceParams,
#define l4 RangeComparisonData
#define l3 ,yA1},{ProduceNewTree,
#define l2 yA1},{ReplaceParams,
#define l1 cMul,SelectedParams,0},0,
#define l0 lB 0x0},{{
#ifdef _MSC_VER
typedef
x13
int
iS1;
#else
#include <stdint.h>
typedef
uint_least32_t
iS1;
#endif
l13
crc32{enum{startvalue=0xFFFFFFFFUL,poly=0xEDB88320UL}
;tH1
iS1
crc>cN2
b8{enum{b1=(crc&x81
crc
x92
crc>>1),b2=(b1&x81
b1
x92
b1>>1),b3=(b2&x81
b2
x92
b2>>1),b4=(b3&x81
b3
x92
b3>>1),b5=(b4&x81
b4
x92
b4>>1),b6=(b5&x81
b5
x92
b5>>1),b7=(b6&x81
b6
x92
b6>>1),res=(b7&x81
b7
x92
b7>>1)}
;}
;inline
iS1
update(iS1
crc,x13
b){
#define B4(n) b8<n>l03 n+1>l03 n+2>l03 n+3>::res
#define R(n) B4(n),B4(n+4),B4(n+8),B4(n+12)
static
const
iS1
table[256]={R(0x00),R(0x10),R(0x20),R(0x30),R(0x40),R(0x50),R(0x60),R(0x70),R(0x80),R(0x90),R(0xA0),R(0xB0),R(0xC0),R(0xD0),R(0xE0),R(0xF0)}
;
#undef R
#undef B4
return((crc>>8))^table[(crc^b)&0xFF];xA1
iS1
calc_upd(iS1
c,const
x13
char*buf,cI3
size){iS1
value=c;for(cI3
p=0;p<size;++p)value=update(value,buf[p])eX
value;xA1
iS1
calc
cQ1
x13
char*buf,cI3
size)yU
calc_upd(startvalue,buf,size);}
}
#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH
nK1
Ref>class
xR{cY3
xR():p(0){}
xR(Ref*b):p(b){x33}
xR
cQ1
xR&b):p(b.p){x33
xA1
Ref&eD1*(xC2*p;xA1
Ref*eD1->(xC2
p;}
xR&eD1=(Ref*b){Set(b)eX*this;}
xR&eD1=cQ1
xR&b){Set(b.p)eX*this;}
#ifdef FP_SUPPORT_CXX11_MOVE
xR(xR&&b):p(b.p){b.p=0;}
xR&eD1=(xR&&b){if(p!=b.p){cM2;p=b.p;b.p=0;}
return*this;}
#endif
~xR(){cM2
iD2
UnsafeSetP(Ref*newp){p=newp
iD2
swap(xR<Ref>&b){Ref*tmp=p;p=b.p;b.p=tmp;}
private:inline
static
void
Have(Ref*p2);inline
void
cM2;inline
void
x33
inline
void
Set(Ref*p2);private:Ref*p;}
;nJ
cM2{if(!p)return;p->x43-=1;if(!p->x43)delete
p;}
nJ
Have(Ref*p2){if(p2)++(p2->x43);}
nJ
Birth(){Have(p);}
nJ
Set(Ref*p2){Have(p2);cM2;p=p2;}
#endif
#include <utility>
cN2
Compare2ndRev{nK1
T>inline
tW2()cQ1
T&x53
T&b
xC2
a.second>b.tJ3}
}
;cN2
Compare1st{nK1
cA
const
nL1
x53
nL1
b
xC2
a
t03<b
t03;}
nK1
cA
const
nL1
a,T1
b
xC2
a
t03<b;}
nK1
cA
T1
x53
nL1
b
xC2
a<b
t03;}
}
;
#ifndef FPoptimizerHashHH
#define FPoptimizerHashHH
#ifdef _MSC_VER
typedef
x13
long
long
x11;
#define FPHASH_CONST(x) x##ULL
#else
#include <stdint.h>
typedef
uint_fast64_t
x11;
#define FPHASH_CONST(x) x##ULL
#endif
l13
FUNCTIONPARSERTYPES{cN2
fphash_t{x11
hash1,hash2;fphash_t():hash1(0),hash2(0){}
fphash_t
cQ1
x11&x53
x11&b):hash1(a),hash2(b){}
tW2==cQ1
fphash_t&l23==cO2&&hash2==cP2
tW2!=cQ1
fphash_t&l23!=cO2||hash2!=cP2
tW2<cQ1
fphash_t&l23!=cO2?hash1<cO2:hash2<cP2}
;}
#endif
#ifndef FPOptimizer_CodeTreeHH
#define FPOptimizer_CodeTreeHH
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
l13
l51{cN2
Grammar;}
l13
x91{x63
ByteCodeSynth;}
tJ{x63
CodeTree;eY
cN2
x82;x63
CodeTree{typedef
xR<x82
yB>DataP;DataP
data;cY3
CodeTree();~CodeTree();cN2
OpcodeTag{}
;eF
x42
o,OpcodeTag);cN2
FuncOpcodeTag{}
;eF
x42
o,x13
f,FuncOpcodeTag);cN2
x73{}
;eF
const
eB1&v,x73);
#ifdef FP_SUPPORT_CXX11_MOVE
eF
eB1&&v,x73);
#endif
cN2
VarTag{}
;eF
x13
varno,VarTag);cN2
CloneTag{}
;eF
l33
b,CloneTag);void
GenerateFrom
cQ1
typename
FunctionParserBase
yB::Data&data,bool
keep_powi=false);void
GenerateFrom
cQ1
typename
FunctionParserBase
yB::Data&data,const
x6&x52,bool
keep_powi=false);void
SynthesizeByteCode(std::vector<x13>&cG3,std::vector
yB&immed,cI3&stacktop_max);void
SynthesizeByteCode(x91
eW2&synth,bool
MustPopTemps=true)const;cI3
SynthCommonSubExpressions(x91::xV1
const;void
SetParams
cQ1
x6&x83
SetParamsMove(x6&tE1
CodeTree
GetUniqueRef();
#ifdef FP_SUPPORT_CXX11_MOVE
void
SetParams(x6&&tE1
#endif
void
SetParam(cI3
which,l33
b);void
n71
cI3
which,x32
b);void
AddParam
cQ1
x32
param);void
eC1
x32
param);void
AddParams
cQ1
x6&x83
AddParamsMove(x6&x83
AddParamsMove(x6&lI2,cI3
lJ2);void
DelParam(cI3
index);void
DelParams();void
Become
cQ1
x32
b);inline
cI3
yQ1
const
yU
tF1
eZ3;xA1
x32
GetParam
n03)yU
tF1[n];xA1
l33
GetParam
n03
xC2
tF1[n];xA1
void
eO2
x42
o)e03
Opcode=o;xA1
x42
GetOpcode()yT
Opcode;xA1
nE
fphash_t
GetHash()yT
Hash;xA1
const
x6&tF1
const
yU
y2;xA1
x6&tF1
yU
y2;xA1
cI3
GetDepth()yT
Depth;xA1
const
eB1&GetImmed()yT
Value;xA1
x13
GetVar()yT
iT1
xA1
x13
GetFuncNo()yT
iT1
xA1
bool
IsDefined(xC2
GetOpcode()!=nE
cNop;xA1
bool
IsImmed(xC2
GetOpcode()==nE
cImmed;xA1
bool
IsVar(xC2
GetOpcode()==nE
iJ2;xA1
x13
GetRefCount()yT
x43
iD2
ReplaceWithImmed
cP1
i);void
Rehash(bool
constantfolding=true);void
Sort();inline
void
Mark_Incompletely_Hashed()e03
Depth=0;xA1
bool
Is_Incompletely_Hashed()yT
Depth==0;xA1
const
tU
GetOptimizedUsing()yT
iV1;xA1
void
SetOptimizedUsing
cQ1
tU
g)e03
iV1=g;}
bool
RecreateInversionsAndNegations(bool
prefer_base2=false);void
FixIncompleteHashes();void
swap(x32
b){data.swap(b.data);}
bool
IsIdenticalTo
cQ1
x32
b)const;void
iY}
;eY
cN2
x82{int
x43;x42
Opcode;eB1
Value;x13
iT1
n72
cJ3;nE
fphash_t
Hash;cI3
Depth;const
tU
iV1;x82();x82
cQ1
x82&b);i4
x42
o);i4
x42
o,x13
f);i4
const
eB1&i);
#ifdef FP_SUPPORT_CXX11_MOVE
i4
eB1&&i);x82(x82&&b);
#endif
bool
IsIdenticalTo
cQ1
x82&b)const;void
Sort();void
Recalculate_Hash_NoRecursion();private:void
eD1=cQ1
x82&b);}
;eY
l41
CodeTreeImmed
cP1
i)yU
xE(i
l71
x73());}
#ifdef FP_SUPPORT_CXX11_MOVE
eY
l41
CodeTreeImmed(eB1&&i)yU
xE
cI2
i)l71
x73());}
#endif
eY
l41
CodeTreeOp(x42
opcode)yU
xE(opcode
l71
OpcodeTag());}
eY
l41
CodeTreeFuncOp(x42
cZ3
x13
f)yU
xE(cZ3
f
l71
FuncOpcodeTag());}
eY
l41
CodeTreeVar
nI2
varno)yU
xE(varno
l71
VarTag());}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
xC3
DumpHashes(cX);xC3
DumpTree(cX);xC3
DumpTreeWithIndent(cX,const
n13&indent="\\"
);
#endif
}
#endif
#endif
#ifndef FPOptimizer_GrammarHH
#define FPOptimizer_GrammarHH
#include <iostream>
tJ{x63
CodeTree;}
l13
l51{enum
ImmedConstraint_Value{ValueMask=0x07,Value_AnyNum=0x0,nN2=0x1,Value_OddInt=0x2,tR1=0x3,Value_NonInteger=0x4,eN1=0x5
iG2
ImmedConstraint_Sign{SignMask=0x18,Sign_AnySign=0x00,nB1=0x08,eO1=0x10,Sign_NoIdea=0x18
iG2
ImmedConstraint_Oneness{OnenessMask=0x60,Oneness_Any=0x00,Oneness_One=0x20,Oneness_NotOne=0x40
iG2
ImmedConstraint_Constness{ConstnessMask=0x180,Constness_Any=0x00,tQ1=0x80,Constness_NotConst=0x100
iG2
Modulo_Mode{Modulo_None=0,Modulo_Radians=1
iG2
Situation_Flags{LogicalContextOnly=0x01,NotForIntegers=0x02,OnlyForIntegers=0x04,OnlyForComplex=0x08,NotForComplex=0x10
iG2
nG2{NumConstant,nC3,SubFunction
iG2
ParamMatchingType{PositionalParams,SelectedParams,AnyParams,GroupFunction
iG2
RuleType{ProduceNewTree,ReplaceParams}
;
#ifdef __GNUC__
# define PACKED_GRAMMAR_ATTRIBUTE __attribute__((packed))
#else
# define PACKED_GRAMMAR_ATTRIBUTE
#endif
typedef
std::pair<nG2,const
void*>cQ2;eY
cQ2
e81
nI2
paramlist,lC1);eY
bool
ParamSpec_Compare
cQ1
void*x53
void*b,nG2
type);x13
ParamSpec_GetDepCode
cQ1
cQ2&b);cN2
ParamSpec_ParamHolder{lC1:8;x13
constraints:9;x13
depcode:15;eD
eY
cN2
ParamSpec_NumConstant{eB1
constvalue;x13
modulo;}
;cN2
iW{x13
param_count:2;x13
param_list:30;x42
subfunc_opcode:8;ParamMatchingType
match_type:3;x13
n1:5;eD
cN2
ParamSpec_SubFunction{iW
data;x13
constraints:9;x13
depcode:7;eD
cN2
Rule{RuleType
ruletype:2;x13
situation_flags:5;x13
repl_param_count:2+9;x13
repl_param_list:30;iW
match_tree;eD
cN2
Grammar{x13
rule_count;x13
short
rule_list[999
cY
extern
const
Rule
grammar_rules[];}
xC3
DumpParam
cQ1
cQ2&p,std::ostream&o=std::cout);xC3
DumpParams
nI2
paramlist,x13
count,std::ostream&o=std::cout);}
#endif
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#define CONSTANT_POS_INF HUGE_VAL
#define CONSTANT_NEG_INF (-HUGE_VAL)
l13
FUNCTIONPARSERTYPES{eY
inline
eB1
fp_const_pihalf()yU
fp_const_pi
yB()*eB1(0.5);}
eY
inline
eB1
fp_const_twopi(){eB1
cQ3
fp_const_pi
yB());cW3
fp_const_twoe(){eB1
cQ3
fp_const_e
yB());cW3
fp_const_twoeinv(){eB1
cQ3
fp_const_einv
yB());cW3
fp_const_negativezero()yU-Epsilon
yB::value;}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <vector>
#include <utility>
#include <iostream>
cO1{using
l13
l51;using
tJ;using
l13
FUNCTIONPARSERTYPES;x63
MatchInfo{cY3
std::vector<std::pair<bool,n72> >lR;n72
yJ;std::vector<x13>lV1;cY3
MatchInfo():lR(),yJ(),lV1(){}
cY3
bool
SaveOrTestRestHolder
nI2
n1,x93
i01){cI1{lR
xA3
n1+1);lR
tV=i01
l43
if(lR[n1
e43==false){lR
tV=i01
l43
x93
found=lR[n1
eY3;if(i01
eZ3!=found
eZ3
cI
eT2
i01
i52
a)if(!i01[a]xI
found[a])cI
return
true
iD2
SaveRestHolder
nI2
n1,n72&i01){cI1
lR
xA3
n1+1);lR
tV.swap(i01);}
bool
SaveOrTestParamHolder
nI2
x1,iW1
xB3){if(yJ
eZ3<=x1){yJ.xF3
x1+1);yJ
xA3
x1);yJ.push_back(xB3
lA2
if(!yJ[x1]lD2){yJ[x1]=xB3
l43
return
xB3
xI
yJ[x1]i72
SaveMatchedParamIndex(lC1){lV1.push_back(index);}
iW1
GetParamHolderValueIfFound
nI2
x1)const{static
const
xE
dummytree;if(yJ
eZ3<=x1)return
dummytree
eX
yJ[x1];}
iW1
GetParamHolderValue
nI2
x1
xC2
yJ[x1];}
bool
HasRestHolder
nI2
n1
xC2
lR
eZ3>n1&&lR[n1
e43==true;}
x93
GetRestHolderValues
nI2
n1)const{static
xT1
xE>empty_result;cI1
return
empty_result
eX
lR[n1
eY3;}
xT1
x13>&GetMatchedParamIndexes(xC2
lV1
iD2
swap(xZ
b){lR.swap(b.lR);yJ.swap(b.yJ);lV1.swap(b.lV1);}
xZ
eD1=cQ1
xZ
b){lR=b.lR;yJ=b.yJ;lV1=b.lV1
eX*this;}
}
;class
eG;typedef
xR<eG>eB;class
eG{cY3
int
x43;cY3
eG():x43(0){}
virtual~eG(){}
}
;cN2
lZ1{bool
found;eB
specs;lZ1(bool
f):found(f),specs(){}
lZ1(bool
f
iT2
s):found(f),specs(s){}
}
;xC3
SynthesizeRule
cQ1
eP2
xE&tree,xZ
info);eY
lZ1
TestParam
cQ1
cQ2&y42
xD3
iT2
nZ3,xZ
info);eY
lZ1
TestParams(eZ&nS,xD3
iT2
nZ3,xZ
info,bool
lB2;eY
bool
ApplyGrammar
cQ1
Grammar&tB2,FPoptimizer_CodeTree::xE&tree,bool
from_logical_context=false);xC3
ApplyGrammars(FPoptimizer_CodeTree::cS2;eY
bool
IsLogisticallyPlausibleParamsMatch(eZ&t33,xD3);}
l13
l51{xC3
DumpMatch
cQ1
eP2
nT,const
FPoptimizer_Optimize::xZ
info,bool
DidMatch,std::ostream&o=std::cout);xC3
DumpMatch
cQ1
eP2
nT,const
FPoptimizer_Optimize::xZ
info,iU2
t43,std::ostream&o=std::cout);}
#endif
#include <string>
iV2
l51::nG2
yM1=false);iV2
x42
yM1=false);
#include <string>
#include <sstream>
#include <assert.h>
#include <iostream>
using
l13
l51;using
l13
FUNCTIONPARSERTYPES;iV2
l51::nG2
yM1){
#if 1
iU2
p=0;switch(opcode
nB3
eM3
p="NumConstant"
;lD
nC3:p="ParamHolder"
;lD
SubFunction:p="SubFunction"
;t82
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str()eZ3<12)tmp<<' '
eX
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str()eZ3<5)tmp<<' '
eX
tmp.str();
#endif
}
iV2
x42
yM1){
#if 1
iU2
p=0;switch(opcode
nB3
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
tU1:p="cNEqual"
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
cD3:p="cLog2by"
;lD
cNop:p="cNop"
;break;
#endif
case
cSinCos:p="cSinCos"
;lD
cSinhCosh:p="cSinhCosh"
;lD
nR2:p="cAbsNot"
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
iJ2:p="VarBegin"
;t82
std::ostringstream
tmp;assert(p);tmp<<p;if(pad)while(tmp.str()eZ3<12)tmp<<' '
eX
tmp.str();
#else
std::ostringstream
tmp;tmp<<opcode;if(pad)while(tmp.str()eZ3<5)tmp<<' '
eX
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
;l13
x91{x63
ByteCodeSynth{cY3
ByteCodeSynth():eE1(),Immed(),cK(),xS(0),StackMax(0){eE1.xF3
64);Immed.xF3
8);cK.xF3
16
i72
Pull(std::vector<x13>&bc,std::vector
yB&imm,cI3&StackTop_max){for(eQ1=0;a<n01;++a){eE1[a]&=~0x80000000u;}
eE1.swap(bc);Immed.swap(imm);StackTop_max=StackMax;}
cI3
GetByteCodeSize(xC2
n01;}
cI3
GetStackTop(xC2
xS
iD2
PushVar
nI2
varno){iX1
varno);eQ2}
void
PushImmed(eB1
immed
nN
iX1
cImmed);Immed.push_back(immed);eQ2}
void
StackTopIs(nT,int
offset=0){if((int)xS>offset){cK[xS
n22
first=true;cK[xS
n22
second=tree;}
}
bool
IsStackTop(nT,int
offset=0
xC2(int)xS>offset&&cK[xS
n22
first&&cK[xS
n22
second
xI
tree);xA1
void
EatNParams
nI2
eat_count){xS-=eat_count
iD2
ProducedNParams
nI2
produce_count){xB1
xS+produce_count
i72
DoPopNMov(cI3
cR2,cI3
srcpos
nN
iX1
cPopNMov
nH2
cR2
nH2
srcpos);xB1
srcpos+1);cK[cR2]=cK[srcpos];xB1
cR2+1
i72
DoDup(cI3
xG3
nN
if(xG3==xS-1){iX1
cDup);}
else{iX1
cFetch
nH2
xG3);}
eQ2
cK[xS-1]=cK[xG3];}
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
tH1
int>void
Dump(){std::ostream&o=std::cout;o<<"Stack state now("
<<xS<<"):\n"
l01
xS;++a){o<<a<<": "
;if(cK[a
e43){nT=cK[a
eY3;o<<'['<<std::hex<<(void*)(&tree.tF1)<<std::dec<<','<<tree.GetRefCount()<<']'
eL3
tree,o);}
else
o<<"?"
;o<<"\n"
;}
o<<std::flush;}
#endif
cI3
xH3
nT)const{y33
xS;a-->0;)if(cK[a
e43&&cK[a
eY3
xI
tree
iE2
a
eX
eR2;}
bool
Find(nT
xC2
xH3
tree)!=eR2;}
bool
FindAndDup(nT){cI3
pos=xH3
tree);if(pos!=eR2){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<i33"duplicate at ["
<<pos<<"]: "
eL3
tree)xR1" -- issuing cDup or cFetch\n"
;
#endif
DoDup(pos
lA2
return
eC3
cN2
IfData{cI3
ofs;}
;void
SynthIfStep1(IfData&nT3,x42
op
nT1
nT3.ofs=n01;iX1
op
cW
cW
i72
SynthIfStep2(IfData&nT3
nT1
xI3
xJ3+2);xI3
2
y8
iY1
nT3.ofs=n01;iX1
cJump
cW
cW
i72
SynthIfStep3(IfData&nT3
nT1
eE1.back()|=0x80000000u;xI3
xJ3-1);xI3
2
y8
iY1
xB1
xS+1)l01
nT3.ofs;++a){if(eE1[a]==cJump&&eE1[a+1]==(0x80000000u|(nT3.ofs-1))){eE1[a+xJ3-1);eE1[a+2
y8
iY1
iQ2
eE1[a]nB3
cAbsIf:case
cIf:case
cJump:case
cPopNMov:a+=2;lD
cFCall:case
cPCall:case
cFetch:a+=1;break;y83
t82}
}
protected:void
xB1
cI3
value){xS
xN3;if(xS>l53{StackMax=xS;cK
xA3
l53;}
}
protected:std::vector<x13>eE1;std::vector
yB
Immed;std::vector<std::pair<bool,FPoptimizer_CodeTree::xE> >cK;cI3
xS;cI3
StackMax;private:void
incStackPtr(){if(xS+2>l53
cK
xA3
StackMax=xS+2);}
tH1
bool
IsIntType,bool
IsComplexType>cN2
yU2{}
;cY3
void
AddOperation
nI2
cZ3
x13
eat_count,x13
produce_count=1){EatNParams(eat_count);lM1
opcode);ProducedNParams(produce_count
i72
lM1
x13
cZ3
yU2<false,false>yV
false,true>yV
true,false>yV
true,true>);inline
void
lM1
x13
opcode){lM1
cZ3
yU2<bool(nE
IsIntType
yB::cX3),bool(nE
IsComplexType
yB::cX3)>());}
}
;eY
cN2
SequenceOpCode;eY
cN2
tM1{static
xX
AddSequence;static
xX
MulSequence;}
;xC3
nR1
long
count,xX&t5,xV1;}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l13
FUNCTIONPARSERTYPES;l13
x91{eY
cN2
SequenceOpCode{eB1
basevalue;x13
op_flip;x13
op_normal,op_normal_flip;x13
op_inverse,op_inverse_flip;}
;eY
xX
tM1
yB::AddSequence={x71,cNeg
yK
cAdd,cSub,cRSub}
;eY
xX
tM1
yB::MulSequence={eB1(1),cInv,cMul,cMul,cDiv,cRDiv}
;
#define findName(a,b,c) "var"
#define TryCompilePowi(o) false
#define mData this
#define mByteCode eE1
#define mImmed Immed
n32
false,false>){xC1
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
n32
true,false>){xC1
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
n32
false,true>){xC1
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
n32
true,true>){xC1
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
l13
x91;
#define POWI_TABLE_SIZE 256
#define POWI_WINDOW_SIZE 3
l13
x91{
#ifndef FP_GENERATING_POWI_TABLE
extern
const
x13
char
powi_table[POWI_TABLE_SIZE];const
#endif
x13
char
powi_table[POWI_TABLE_SIZE]={0,1,1,1,2,1,2,1,xK3
4,1,2,xL3
2,1,xK3
8,cK3
xM3
15,1,16,1,2,1,4,1,2,xL3
2,1,4,cK3
1,16,1,25,xM3
27,5,8,3,2,1,30,1,31,3,32,1,2,1,xK3
8,1,2,xM3
39,1,16,137,2,1,4,cK3
xL3
45,135,4,31,2,5,32,1,2,131,50,1,51,1,8,3,2,1,54,1,55,3,16,1,57,133,4,137,2,135,60,1,61,3,62,133,63,1,i21
131,i21
139,iZ1
eH
30,1,130,137,2,31,iZ1
eH
eH
130,cK3
1,eH
eH
2,1,130,133,i21
61,130,133,62,139,130,137,eH
iZ1
eH
eH
i21
131,eH
eH
130,131,2,133,iZ1
130,141,eH
130,cK3
1,eH
5,135,eH
iZ1
eH
iZ1
130,133,130,141,130,131,eH
eH
2,131}
;}
static
nQ3
c5=256;
#define FPO(x)
l13{class
l02{private:int
iC[c5];int
i51[c5];cY3
l02():iC(),i51(){iC[1]=1;}
bool
Plan_Add(yY1,int
count){c61>=c5
cI
i51[eY2+=count
eX
iC[eY2!=0
iD2
l73
yY1){c61<c5)iC[eY2=1
iD2
Start(cI3
value1_pos){for(int
n=2;n<c5;++n)iC[n]=-1;Remember(1,value1_pos);DumpContents();}
int
Find(yY1)const{c61<c5){if(iC[eY2>=0){FPO(tP3(tT3,"* I found %ld from cache (%u,%d)\n",value,(unsigned)cache[value],tQ3 value]))eX
iC[eY2;eG2-1
iD2
Remember(yY1,cI3
i23){c61>=c5)return;FPO(tP3(tT3,"* Remembering that %ld can be found at %u (%d uses remain)\n",value,(unsigned)i23,tQ3 value]));iC[eY2=(int)i23
iD2
DumpContents()const{FPO(for(int a=1;a<POWI_CACHE_SIZE;++a)if(cache[a]>=0||tQ3 a]>0){tP3(tT3,"== cache: sp=%d, val=%d, needs=%d\n",cache[a],a,tQ3 a]);})}
int
UseGetNeeded(yY1){c61>=0&&value<c5)return--i51[eY2
eX
0;}
}
;eY
cI3
y9
long
count
lV3
xX&t5,xV1;xC3
yZ1
cI3
apos,long
aval,cI3
bpos,long
bval
lV3
x13
cumulation_opcode,x13
cimulation_opcode_flip,xV1;void
l61
yY1
lV3
int
need_count,int
iX=0){c61<1)return;
#ifdef FP_GENERATING_POWI_TABLE
if(iX>32)throw
false;
#endif
if(iC.Plan_Add(value,need_count
iE2;long
half=1;c61<POWI_TABLE_SIZE){half=powi_table[eY2;if
tX2
128){half&=127;if
tX2
64)half
lW3
FPO(tP3(tT3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,value/half));l61
half
cC3
iC.l73
half)eX;}
l72
half&64){half
lW3}
}
else
c61&1)half
xN3&((1<<POWI_WINDOW_SIZE)-1);else
half
xN3/2;long
cP
xN3-half;if(half>cP||half<0)std::swap(half,cP);FPO(tP3(tT3,"value=%ld, half=%ld, otherhalf=%ld\n",value,half,otherhalf));if(half==cP){l61
half,iC,2,iX+1);}
else{l61
half
cC3
l61
cP>0?cP:-cP
cC3}
iC.l73
value);}
eY
cI3
y9
yY1
lV3
xX&t5,xV1{int
xO3=iC.Find(value);if(xO3>=0)yU
xO3;}
long
half=1;c61<POWI_TABLE_SIZE){half=powi_table[eY2;if
tX2
128){half&=127;if
tX2
64)half
lW3
FPO(tP3(tT3,"* I want %ld, my plan is %ld * %ld\n",value,half,value/half));cI3
xD2=y9
half
e31
if(iC
x02
half)>0||xD2!=e41){i61
xD2);iC
eZ2
half,e41);}
nR1
value/half
nY2
cI3
i23=e41;iC
eZ2
value,i23);iC.DumpContents()eX
i23;}
l72
half&64){half
lW3}
}
else
c61&1)half
xN3&((1<<POWI_WINDOW_SIZE)-1);else
half
xN3/2;long
cP
xN3-half;if(half>cP||half<0)std::swap(half,cP);FPO(tP3(tT3,"* I want %ld, my plan is %ld + %ld\n",value,half,value-half));if(half==cP){cI3
xD2=y9
half
e31
yZ1
xD2,half,xD2,half,iC,t5.op_normal,t5.op_normal_flip,synth);}
else{long
part1=half;long
part2=cP>0?cP:-cP;cI3
part1_pos=y9
part1
e31
cI3
part2_pos=y9
part2
e31
FPO(tP3(tT3,"Subdivide(%ld: %ld, %ld)\n",value,half,otherhalf));yZ1
part1_pos,part1,part2_pos,part2,iC,cP>0?t5.op_normal:t5.op_inverse,cP>0?t5.op_normal_flip:t5.op_inverse_flip,synth);}
cI3
i23=e41;iC
eZ2
value,i23);iC.DumpContents()eX
i23;}
xC3
yZ1
cI3
apos,long
aval,cI3
bpos,long
bval
lV3
x13
cumulation_opcode,x13
cumulation_opcode_flip,xV1{int
a_needed=iC
x02
aval);int
xP3=iC
x02
bval);bool
flipped
e63
#define DUP_BOTH() do{if(apos<bpos){cI3 tmp=apos;apos=bpos;bpos=tmp t21 FPO(tP3(tT3,"-> " i03 i03"op\n",(unsigned)apos,(unsigned)bpos));i61 apos);i61 apos==bpos?e41:bpos);}while(0)
#define DUP_ONE(p) do{FPO(tP3(tT3,"-> " i03"op\n",(unsigned)p));i61 p);}while(0)
if(a_needed>0){if(xP3>0){n42}
c03
bpos!=e41)n42
else{nJ2
t21}
}
l72
xP3>0){if(apos!=e41)n42
else
DUP_ONE(bpos);}
c03
apos==bpos&&apos==e41)nJ2;t02
e41&&bpos==nS3
xH
2){FPO(tP3(tT3,"-> op\n"))t21
t02
nS3
xH
2&&bpos==e41)FPO(tP3(tT3,"-> op\n"));t02
e41)DUP_ONE(bpos);l72
bpos==e41){nJ2
t21
else
n42}
nZ2
flipped?cumulation_opcode_flip:cumulation_opcode,2);}
xC3
c71
long
count,xX&t5,xV1{while
e53<256){int
half=x91::powi_table[count];if
tX2
128){half&=127;c71
half
nY2
count/=half;}
else
t82
if
e53==1)return;if(!e53&1)){nZ2
cSqr,1);c71
count/2
nY2}
else{i61
e41);c71
count-1
nY2
nZ2
cMul,2);}
}
}
l13
x91{xC3
nR1
long
count,xX&t5,xV1{if
e53==0)tD1
t5.basevalue);else{bool
t12
e63
if
e53<0){t12=true;count=-count;}
if(false)c71
count
nY2
l72
count>1){l02
iC;l61
count,iC,1);cI3
xD1=nS3
GetStackTop();iC.Start(e41);FPO(tP3(tT3,"Calculating result for %ld...\n",count));cI3
xE2=y9
count
e31
cI3
n_excess=nS3
xH
xD1;if(n_excess>0||xE2!=xD1-1){nS3
DoPopNMov(xD1-1,xE2);}
}
if(t12)nZ2
t5.op_flip,1);}
}
}
#endif
#ifndef FPOptimizer_ValueRangeHH
#define FPOptimizer_ValueRangeHH
tJ{l13
l83{iP
cN2
Comp{}
;tH1>cN2
Comp<nE
cLess>lX3<nB
cLessOrEq>lX3<=nB
cGreater>lX3>nB
cGreaterOrEq>lX3>=nB
cEqual>lX3==nB
tU1>lX3!=b;}
}
;}
eY
cN2
rangehalf{eB1
val;bool
known;rangehalf():val(),known(false){}
rangehalf
cP1
v):val(v),known(true){xA1
void
set
cP1
v){known=true;val=v
iD2
set(eB1(i71(eB1
nP)val=i81
void
set(eB1(i71
cP1
nP)val=i81
iP
void
set_if(eB1
v,eB1(i71(eB1
nP&&l83::Comp<Compare>()(val,v))val=i81
iP
void
set_if
cP1
v,eB1(i71
cP1
nP&&l83::Comp<Compare>()(val,v))val=i81}
;eY
cN2
range{rangehalf
yB
min,max;range():min(),max(){}
range(eB1
mi,eB1
ma):min(mi),max(ma){}
range(bool,eB1
ma):min(),max(ma){}
range(eB1
mi,bool):min(mi),max(){}
void
set_abs();void
set_neg();}
;eY
bool
IsLogicalTrueValue
cQ1
range
yB&p,bool
abs);eY
bool
IsLogicalFalseValue
cQ1
range
yB&p,bool
abs);}
#endif
#ifndef FPOptimizer_RangeEstimationHH
#define FPOptimizer_RangeEstimationHH
tJ{enum
TriTruthValue{x03,IsNever,Unknown}
;eY
range
yB
iQ
xD3);eY
bool
IsLogicalValue
cQ1
cS2;eY
TriTruthValue
GetIntegerInfo
cQ1
cS2;eY
xE1
GetEvennessInfo
cQ1
cS2{if(!tree
l92)return
Unknown;yJ1=tree
tJ1;if(nE
isEvenInteger(value
l32
nE
isOddInteger(value
iF2
eY
xE1
GetPositivityInfo
cQ1
cS2{range
yB
p=iQ
tree);if(p
xU1
p
iA2>=eB1(l32
p
e3
yS1
iF2
eY
xE1
GetLogicalValue
lB3
tree,bool
abs){range
yB
p=iQ
tree);if(IsLogicalTrueValue(p,abs
l32
IsLogicalFalseValue(p,abs
iF2}
#endif
#ifndef FPOptimizer_ConstantFoldingHH
#define FPOptimizer_ConstantFoldingHH
tJ{xC3
ConstantFolding(cS2;}
#endif
l13{using
l13
FUNCTIONPARSERTYPES;using
tJ;cN2
ComparisonSetBase{enum{e13=0x1,Eq_Mask=0x2,Le_Mask=0x3,e23=0x4,e33=0x5,Ge_Mask=0x6}
;static
int
Swap_Mask(int
m)yU(m&Eq_Mask)|((m&e13)?e23:0)|((m&e23)?e13:0);}
enum
c01{Ok,BecomeZero,BecomeOne,n51
iG2
nK2{cond_or,tY2,tZ2,i02}
;}
;eY
cN2
ComparisonSet:public
ComparisonSetBase{cN2
t22{xE
a;xE
b;int
relationship;t22():a(),b(),relationship(){}
}
;std::vector<t22>e1;cN2
Item{xE
value;bool
yV2;Item():value(),yV2(false){}
}
;std::vector<Item>c81;int
xF1;ComparisonSet():e1(),c81(),xF1(0){}
c01
AddItem
lB3
a,bool
yV2,nK2
type){for(cI3
c=0;c<c81
i52
c)if(c81[c].value
xI
a)){if(yV2!=c81[c].yV2){iR
cJ1
case
i02:c81.erase(c81.begin()+c);xF1
yW2
case
tY2:case
tZ2:cK1
eG2
n51;}
Item
pole;pole.value=a;pole.yV2=yV2;c81.push_back(pole)eX
Ok;}
c01
AddRelationship(xE
a,xE
b,int
tN1,nK2
type){iR
if(tN1==7)cJ1
lD
i02:if(tN1==7){xF1
yW2}
lD
tY2:case
tZ2:if(tN1==0)cK1
t82
if(!(a.GetHash()<b.GetHash())){a.swap(b);tN1=Swap_Mask(tN1);}
for(cI3
c=0;c<e1
i52
c){if(e1[c].a
xI
a)&&e1[c].b
xI
b)){iR{int
xQ3=yX2|tN1;if(xQ3==7)cJ1
yX2=xQ3;break
t11
tY2:case
tZ2:{int
xQ3=yX2&tN1;if(xQ3==0)cK1
yX2=xQ3;break
t11
i02:{int
newrel_or=yX2|tN1;int
xG2=yX2&tN1;l12
5&&xG2==0){yX2=e33
eX
n51;}
l12
7&&xG2==0){xF1+=1;e1.erase(e1.begin()+c)eX
n51;}
l12
7&&xG2==Eq_Mask){yX2=Eq_Mask;xF1
yW2}
yC1
eG2
n51;}
}
t22
comp;comp.a=a;comp.b=b;comp.relationship=tN1;e1.push_back(comp)eX
Ok;}
}
;nK1
eB1,typename
CondType>bool
ConstantFolding_LogicCommon(xE&tree,CondType
x21,bool
xH2){bool
should_regenerate
e63
ComparisonSet
yB
comp
nD2
typename
yO
c01
cL3=yO
Ok;iW1
atree=y1
a);switch(atree
nF
nB3
cEqual
lI
Eq_Mask,x21);lD
tU1
lI
e33,x21);lD
cLess
lI
e13,x21);lD
cLessOrEq
lI
Le_Mask,x21);lD
cGreater
lI
e23,x21);lD
cGreaterOrEq
lI
Ge_Mask,x21);lD
cNot:cL3
t31
l9
0),true,x21);lD
cNotNot:cL3
t31
l9
0),false,x21);break;y83
if(xH2||IsLogicalValue(atree))cL3
t31,false,x21);iQ2
cL3){ReplaceTreeWithZero:xB
0)eX
true;ReplaceTreeWithOne:xB
1);nZ
yO
Ok:lD
yO
BecomeZero
cN
yO
BecomeOne:iB
yO
n51:c31
t82}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_LogicCommon: "
yN
#endif
if(xH2){tree.DelParams();}
else{for
iU{iW1
atree=y1
a);if(IsLogicalValue(atree))iJ);}
}
eT2
comp.c81
i52
a){if(comp.c81[a].yV2
l93
cNot);r.cB
r.l22
l72!xH2
l93
cNotNot);r.cB
r.l22
else
tree.cB}
eT2
comp.e1
i52
a
l93
cNop);switch(comp.e1[a
i0
nB3
yO
e13:r
yX
cLess);lD
yO
Eq_Mask:r
yX
cEqual);lD
yO
e23:r
yX
cGreater);lD
yO
Le_Mask:r
yX
cLessOrEq);lD
yO
e33:r
yX
tU1);lD
yO
Ge_Mask:r
yX
cGreaterOrEq);t82
r
y41
comp.e1[a].a);r
y41
comp.e1[a].b);r.l22
if(comp.xF1!=0)tree
cL
n61
comp.xF1)));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_LogicCommon: "
yN
#endif
return
true;}
return
eC3
eY
bool
ConstantFolding_AndLogic(tU3(tree.GetOpcode()==cAnd
tN3()==cAbsAnd)eX
nK
tY2,true
iW2
ConstantFolding_OrLogic(tU3(tree.GetOpcode()==cOr
tN3()==cAbsOr)eX
nK
cond_or,true
iW2
ConstantFolding_AddLogicItems(tU3(tree.GetOpcode()==cAdd)eX
nK
i02,false
iW2
ConstantFolding_MulLogicItems(tU3(tree.GetOpcode()==cMul)eX
nK
tZ2,false);}
}
#include <vector>
#include <map>
#include <algorithm>
l13{using
l13
FUNCTIONPARSERTYPES;using
tJ;cN2
CollectionSetBase{enum
xG1{Ok,n51}
;}
;eY
cN2
CollectionSet:public
CollectionSetBase{cN2
c11{xE
value;xE
xI2;bool
eC;c11():value(),xI2(),eC(false){}
c11
lB3
v,iW1
f):value(v),xI2(f),eC(false){}
}
;std::multimap<fphash_t,c11>iD;typedef
typename
std::multimap<fphash_t,c11>::xV3
xH1;CollectionSet():iD(){}
xH1
FindIdenticalValueTo
lB3
value){fphash_t
hash
xN3.GetHash();for(xH1
i=iD.xJ2
hash);i!=iD.cL1
hash;++i){c61
xI
i
cY2.value
iE2
i;}
return
iD.end();}
bool
Found
cQ1
xH1&b)yU
b!=iD.end();}
xG1
AddCollectionTo
lB3
xI2,const
xH1&into_which){c11&c=into_which
cY2;if(c.eC)c.xI2
cL
xI2);else{xE
add;add
yX
cAdd);add
y41
c.xI2);add
cL
xI2);c.xI2.swap(add);c.eC=true;}
return
n51;}
xG1
nL2
lB3
value,iW1
xI2){const
fphash_t
hash
xN3.GetHash();xH1
i=iD.xJ2
hash);for(;i!=iD.cL1
hash;++i){if(i
cY2.value
xI
value
iE2
AddCollectionTo(xI2,i);}
iD.xR3,std::make_pair(hash,c11(value,xI2)))eX
Ok;}
xG1
nL2
lB3
a)yU
nL2(a,n61
1)));}
}
;eY
cN2
ConstantExponentCollection{typedef
n72
xT3;typedef
std::nU1
xK2;std::vector<xK2>data;ConstantExponentCollection():data(){}
void
MoveToSet_Unique
cP1
eF1&eG1){data.push_back(std::nU1(eF1()));data.back().second.swap(eG1
i72
MoveToSet_NonUnique
cP1
eF1&eG1){typename
std::vector<xK2>::xV3
i=std::xJ2
data.i12
data.end(),exponent,Compare1st());if(i!=data.cL1
exponent){i
cY2.xR3
cY2.end(),eG1.i12
eG1.end());}
else{data.xR3,std::nU1
e83,eG1));}
}
bool
tV2{bool
changed
e63
std::sort(data.i12
data.end(),Compare1st());redo:eT2
data
i52
a){eB1
exp_a=data[a
e43;lJ3
exp_a,xW3
yC1
for(cI3
b=a+1;b<data
i52
b){eB1
exp_b=data[b
e43;eB1
xL2=exp_b-exp_a;if(xL2>=fp_abs(exp_a
l63
exp_diff_still_probable_integer=xL2*eB1(16);if(t32
exp_diff_still_probable_integer)&&!(t32
exp_b)&&!t32
xL2))){xT3&a_set=l42;xT3&b_set=data[b
eY3;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantExponentCollection iteration:\n"
;t52
cout);
#endif
if(isEvenInteger(exp_b)&&!isEvenInteger(xL2+exp_a)nQ
tmp2;tmp2
y93
tmp2
nB2
b_set);tmp2
lS2
nC
cAbs)t42;tmp
lS2
b_set
xA3
1);b_set[0
t63
tmp);}
a_set.insert(a_set.end(),b_set.i12
b_set.end());xT3
b_copy=b_set;data.erase(data.begin()+b);MoveToSet_NonUnique(xL2,b_copy);nV2
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantExponentCollection iteration:\n"
;t52
cout);
#endif
goto
redo;}
eG2
changed;}
#ifdef DEBUG_SUBSTITUTIONS
void
t52
ostream&out){eT2
data
i52
a){out.precision(12);out<<data[a
e43<<": "
;cV1
l42
i52
b){if(b>0)out<<'*'
eL3
l42[b],out);}
out<<std::endl;}
}
#endif
}
;eY
static
xE
nV1
xE&value,bool&xM
yJ3
value
nF
nB3
cPow:{xE
cZ2
value
l9
1);value.yB1
eX
exponent
t11
cRSqrt:value.yB1;xM=true
eX
n61-0.5));case
cInv:value.yB1;xM=true
eX
n61-1));y83
t82
return
n61
1));}
eY
static
void
eH1
tK1&mul,xD3,iW1
xI2,bool&c21
bool&xM){eT2
yI++a
nQ
value(y1
a));xE
exponent(nV1
value,xM));if(!xI2
l92||xI2
tJ1!=eB1(1.0)nQ
cR1;cR1
y93
cR1
cL
tL2
cR1
cL
xI2);cR1.Rehash()cT1.swap(cR1);}
#if 0 /* FIXME: This does not work */
c61
nF==cMul){if(1){bool
exponent_is_even=exponent
l92&&isEvenInteger
e83
tJ1);cV1
value
eX3{bool
tmp
e63
xE
val(value
l9
b));xE
exp(nV1
val,tmp));if(exponent_is_even||(exp
l92&&isEvenInteger(exp
tJ1))nQ
cR1;cR1
y93
cR1
cL
tL2
cR1
y41
exp);cR1.ConstantFolding();if(!cR1
l92||!isEvenInteger(cR1
tJ1)){goto
cannot_adopt_mul;}
}
}
}
eH1
mul,value,exponent,c21
xM);}
else
cannot_adopt_mul:
#endif
{if(mul.nL2(value,exponent)==CollectionSetBase::n51)c31}
}
}
eY
bool
ConstantFolding_MulGrouping(cS2{bool
xM
e63
bool
should_regenerate
e63
tK1
mul;eH1
mul,tree,n61
1)),c21
xM);typedef
std::pair<xE,n72>eI1;typedef
std::multimap<fphash_t,eI1>c91;c91
tO;eX2
tK1::xH1
j=mul.iD.xU3
j!=mul.iD.end();++j
nQ&value=j
cY2.value;xE&cZ2
j
cY2.xI2;if(j
cY2.eC)exponent
lS2
const
fphash_t
eJ1=exponent.GetHash();typename
c91::xV3
i=tO.xJ2
eJ1);for(;i!=tO.cL1
eJ1;++i)if(i
cY2
t03
xI
exponent)){if(!exponent
l92||!cS1
tJ1,xW3
c31
i
cY2.second.push_back(value);goto
skip_b;}
tO.xR3,std::make_pair(eJ1,std::make_pair
e83,n72(cI3(1),value))));skip_b:;}
#ifdef FP_MUL_COMBINE_EXPONENTS
ConstantExponentCollection
yB
cU1;eX2
c91::xV3
j,i=tO.xU3
i!=tO.end();i=j){j=i;++j;eI1&list=i
cY2;if
t73.lN1
cZ2
list
t03
tJ1;if(!e83==x41)cU1.MoveToSet_Unique
e83,list.second);tO.erase(i);}
}
if(cU1.tV2)c31
#endif
if(should_regenerate
nQ
before=tree;before.iY
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_MulGrouping: "
eL3
before)xR1"\n"
;
#endif
tree.DelParams();eX2
c91::xV3
i=tO.xU3
i!=tO.end();++i){eI1&list=i
cY2;
#ifndef FP_MUL_COMBINE_EXPONENTS
if
t73.lN1
cZ2
list
t03
tJ1;if
e83==x41
yC1
if(cS1,xW3{tree.AddParamsMove(list.second);yC1}
}
#endif
xE
mul;mul
y93
mul
nB2
list.second);mul
lS2
if(xM&&list
t03
l92){if
t73
tJ1==eB1(1)/eB1(3)nQ
cbrt;cbrt
yX
cCbrt);cbrt.eI
cbrt.lR2
cbrt);yC1
cQ
0.5)nQ
sqrt;sqrt
yX
cSqrt);sqrt.eI
sqrt.lR2
sqrt);yC1
cQ-0.5)nQ
rsqrt;rsqrt
yX
cRSqrt);rsqrt.eI
rsqrt.lR2
rsqrt);yC1
cQ-1)nQ
inv;inv
yX
cInv);inv.eI
inv.lR2
inv);yC1}
}
xE
pow
t71.eI
pow
y41
list
t03);pow.lR2
pow);}
#ifdef FP_MUL_COMBINE_EXPONENTS
tO.clear()l01
i5
i52
a){eB1
cZ2
i5[a
e43;if(cS1,xW3{tree.AddParamsMove(i5[a
eY3);yC1}
xE
mul;mul
y93
mul
nB2
i5[a
eY3);mul
lS2
xE
pow
t71.eI
pow
cL
CodeTreeImmed
e83));pow.lR2
pow);}
#endif
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_MulGrouping: "
yN
#endif
return!tree
xI
before);}
return
eC3
eY
bool
ConstantFolding_AddGrouping(cS2{bool
should_regenerate
e63
tK1
add
nD2
if(y1
a)nF==cMul)yC1
if(add.nL2(y1
a))==CollectionSetBase::n51)c31}
eV2
remaining(eN3
cI3
tP=0
nD2
iW1
x23=y1
a);if(x23
nF==cMul){cV1
x23
eX3{if(x23
l9
b)l92)yC1
typename
tK1::xH1
c=add.FindIdenticalValueTo(x23
l9
b));if(add.Found(c)nQ
tmp(x23
l71
CloneTag());tmp.DelParam(b);tmp
lS2
add.AddCollectionTo(tmp,c);c31
goto
done_a;}
}
remaining[a]=true;tP+=1;done_a:;}
}
if(tP>0){if(tP>1){std::vector<std::pair<xE,cI3> >nY;std::multimap<fphash_t,cI3>eK1;bool
lD3
e63
for
eJ{cV1
y1
a)eX3{iW1
p=y1
a)l9
b);const
fphash_t
p_hash=p.GetHash();for(std::multimap<fphash_t,cI3>::const_iterator
i=eK1.xJ2
p_hash);i!=eK1.cL1
p_hash;++i){if(nY[i
cY2
e43
xI
p)){nY[i
cY2
eY3+=1;lD3=true;goto
found_mulgroup_item_dup;}
}
nY.push_back(std::make_pair(p,cI3(1)));eK1.insert(std::make_pair(p_hash,nY
eZ3-1));found_mulgroup_item_dup:;}
}
if(lD3
nQ
e02;{cI3
max=0;for(cI3
p=0;p<nY
i52
p)if(nY[p
eY3<=1)nY[p
eY3=0;else{nY[p
eY3*=nY[p
e43
n62;if(nY[p
eY3>max){e02=nY[p
e43;max=nY[p
eY3;}
}
}
xE
group_add;group_add
yX
cAdd);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Duplicate across some trees: "
eL3
e02)xR1" in "
yN
#endif
for
eJ
cV1
y1
a)eX3
if(e02
xI
y1
a)l9
b))nQ
tmp(y1
a)l71
CloneTag());tmp.DelParam(b);tmp
lS2
group_add
y41
tmp);remaining[a]e63
t82
group_add
lS2
xE
group;group
y93
group
y41
e02);group
y41
group_add);group
lS2
add.nL2(group);c31}
}
for
eJ{if(add.nL2(y1
a))==CollectionSetBase::n51)c31}
}
if(should_regenerate){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before ConstantFolding_AddGrouping: "
yN
#endif
tree.DelParams();eX2
tK1::xH1
j=add.iD.xU3
j!=add.iD.end();++j
nQ&value=j
cY2.value;xE&coeff=j
cY2.xI2;if(j
cY2.eC)coeff
lS2
if(coeff
l92){lJ3
coeff
tJ1,x41)yC1
lJ3
coeff
tJ1,xW3{tree
y41
value);yC1}
}
xE
mul;mul
y93
mul
y41
value);mul
y41
coeff);mul.Rehash(yH
eI}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After ConstantFolding_AddGrouping: "
yN
#endif
return
true;}
return
eC3}
l13{using
l13
FUNCTIONPARSERTYPES;using
tJ;eY
bool
ConstantFolding_IfOperations(tU3(tree.GetOpcode()==cIf
tN3()==cAbsIf);for(;;){if(xZ3
cNot){e93
cIf);y1
0).e12
0)yR3
y1
1).swap(y1
2));}
l72
xZ3
nR2){e93
eB3;y1
0).e12
0)yR3
y1
1).swap(y1
2));}
else
i62
lT1
0),x31
eB3)lC2
tree.e12
1));nZ
iY2
tree.e12
2));nZ
lW1
if(xZ3
cIf||xZ3
cAbsIf
nQ
cond=y1
0);xE
lE3;lE3
eA3
eT3
cNotNot:cAbsNotNot);lE3
xM2
1));ConstantFolding(lE3);xE
lF3;lF3
eA3
eT3
cNotNot:cAbsNotNot);lF3
xM2
2));ConstantFolding(lF3);if(lE3
l92||lF3
l92
nQ
t6;t6
eA3
nF);t6
xM2
1));t6.x0
1));t6.x0
2));t6
lS2
xE
t7;t7
eA3
nF);t7
xM2
2));t7.x0
1));t7.x0
2));t7
cM
cond
nF
xI1
0,cond
l9
0)yH
n71
1,t6
yH
n71
2,t7
lA2}
if(y1
1)nF==y1
2)nF&&(y1
1)nF==cIf||y1
1)nF==eB3
nQ&nY3=y1
1);xE&leaf2=y1
2);if(nY3
l9
0)n52
0))&&n81
1))||nY3
l9
2)n52
2)))nQ
t6;t6
i91
t6.x0
0));t6
xN2
1));t6
xO2
1));t6
lS2
xE
t7;t7
i91
t7.x0
0));t7
xN2
2));t7
xO2
2));t7
cM
nY3
nF
xI1
0
eS3
0)yH
n71
1,t6
yH
n71
2,t7
lA2
if
n81
1))&&nY3
l9
2)n52
2))nQ
t8;t8
i91
t8
y41
y1
0));t8
xN2
0));t8
xO2
0));t8
cM
nY3
nF
yH
n71
0,t8
xI1
2
eS3
2)xI1
1
eS3
1)lA2
if
n81
2))&&nY3
l9
2)n52
1))nQ
e22;e22
yX
leaf2
eT3
cNot:nR2);e22
xO2
0));e22
lS2
xE
t8;t8
i91
t8
y41
y1
0));t8
xN2
0));t8
y41
e22);t8
cM
nY3
nF
yH
n71
0,t8
xI1
2
eS3
2)xI1
1
eS3
1)lA2}
xE&y0=y1
1);xE&y5=y1
2);if(y0
xI
y5)){tree.e12
1)lA2
const
OPCODE
op1=y0
nF;const
OPCODE
op2=y5
nF;if(op1==op2){if(y0
e6==1
nQ
lP
0));xP2
0))yX1()n4
if(y0
e6==2&&y5
e6==2){if(y0
l9
0)xI
y5
l9
0))nQ
param0=y0
l9
0);xE
lP
1));xP2
1))yX1
tA
param0)n4
if(y0
l9
1)xI
y5
l9
1))nQ
param1=y0
l9
1);xE
lP
0));xP2
0))yX1
tA
n41
yH
eC1
param1
lA2}
if(op1==y03
cMul
l52
cAnd
l52
cOr
l52
cAbsAnd
l52
cAbsOr
l52
cMin
l52
cMax){n72
lG3;c6{for(cI3
b=y5
e6;b-->0;){if(y0
lH3
y5
l9
b))){if(lG3
i43){y0.iY
y5.iY}
lG3.push_back(y0
nG3
y5.DelParam(b);y0.cV2}
if(!lG3
i43){y0
lS2
y5.Rehash()l8
op1
yH
SetParamsMove(lG3)n4}
}
if(op1==y03
cMul||(op1==cAnd&&IsLogicalValue(y5))||(op1==cOr&&IsLogicalValue(y5))){c6
if(y0
lH3
y5)){y0.iY
y0.cT2
y0
lS2
xE
cA1=y5;y5=tQ
op1==y03
cOr)i22
op1
yH
eC1
cA1)n4}
if((op1==cAnd
l52
cOr)&&op2==cNotNot
nQ&lL3=y5
l9
0);c6
if(y0
lH3
lL3)){y0.iY
y0.cT2
y0
lS2
xE
cA1=lL3;y5=tQ
op1
yT3
op1
yH
eC1
cA1)n4}
if(op2==cAdd||op2==cMul||(op2==cAnd&&IsLogicalValue(y0))||(op2==cOr&&IsLogicalValue(y0))){y33
y5
eL1
y5
lH3
y0)){y5.iY
y5.cT2
y5
lS2
xE
cB1=y0;y0=tQ
op2==cAdd||op2
yT3
op2
yH
eC1
cB1)n4}
if((op2==cAnd||op2==cOr)&&op1==cNotNot
nQ&lM3=y0
l9
0);y33
y5
eL1
y5
lH3
lM3)){y5.iY
y5.cT2
y5
lS2
xE
cB1=lM3;y0=tQ
op2
yT3
op2
yH
eC1
cB1)n4}
return
eC3}
#include <limits>
l13{using
l13
FUNCTIONPARSERTYPES;using
tJ;eY
int
maxFPExponent()yU
std::numeric_limits
yB::max_exponent;}
eY
bool
nW1
eB1
base,eB1
exponent){if(base<x41
return
true;lJ3
base,x41||lK3
base,eB1(1))cI
return
exponent>=eB1(maxFPExponent
yB())/fp_log2(base
iW2
ConstantFolding_PowOperations(tU3(tree.GetOpcode()==cPow);nU&&y1
1).lN1
const_value=eD3
lS,e0);xB
const_value)eX
eC3
if(y1
l82
lK3
e0,xW3{tree.e12
0)lA2
nU&&lK3
lS,xW3{xB
1)eX
eC3
nU&&y1
1)nF==cMul){bool
xQ2
e63
eB1
lE2=lS;xE
x23=y1
1);y33
x23
eL1
x23
l9
a).lN1
imm=x23
l9
a)tJ1;{if(nW1
lE2,imm
l63
cC1=eD3
lE2,imm);lJ3
cC1,x41||lK3
cC1,xW3
break;if(!xQ2){xQ2=true;y51
iY}
lE2=cC1;y51
cV2
if(xQ2){y51
Rehash();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before pow-mul change: "
yN
#endif
y1
0).Become(cX1
lE2));y1
1).Become(x72;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After pow-mul change: "
yN
#endif
}
}
if(y1
l82
xZ3
cMul){eB1
lF2=e0;eB1
xR2=1.0;bool
xQ2
e63
xE&x23=y1
0);y33
x23
eL1
x23
l9
a).lN1
imm=x23
l9
a)tJ1;{if(nW1
imm,lF2
l63
eM1=eD3
imm,lF2);lJ3
eM1,x41)break;if(!xQ2){xQ2=true;y51
iY}
xR2*=eM1;y51
cV2
if(xQ2){y51
Rehash();xE
cM3;cM3
yX
cPow);cM3
nB2
tree.tF1);cM3.y72
yH
eO2
cMul
yH
eC1
cM3
yH
AddParam(cX1
xR2)lA2}
if(xZ3
cPow&&y1
l82
y1
0)l9
1).lN1
a=y1
0)l9
1)tJ1;eB1
b=e0;eB1
c=a*b;if(isEvenInteger(a)&&!isEvenInteger(c)nQ
lN3;lN3
yX
cAbs);lN3.x0
0)yR3
lN3.Rehash(yH
n71
0,lN3);}
else
tree.SetParam(0,y1
0)l9
0)xI1
1,cX1
c));}
return
eC3}
l13{using
l13
FUNCTIONPARSERTYPES;using
tJ;cN2
l4{enum
e32{MakeFalse=0,MakeTrue=1,t62=2,lQ3=3,MakeNotNotP0=4,MakeNotNotP1=5,MakeNotP0=6,MakeNotP1=7,xK=8
iG2
lG2{Never=0,Eq0=1,Eq1=2,y63=3,y73=4}
;e32
if_identical;e32
lH2
4];cN2{e32
what:4;lG2
when:4;}
iA1,iB1,iC1,iD1;eY
e32
Analyze
lB3
a,iW1
b)const{if(a
xI
b
iE2
if_identical
yG2
p0=iQ
a)yG2
p1=iQ
b);if(p0
e3&&p1
cU2){if(p0
xF<p1
iA2&&lH2
0]i6
0];if(p0
xF<=p1
iA2&&lH2
1]i6
1];}
if(p0
xU1
p1
e3){if
cF2>p1
xF&&lH2
2]i6
2];if
cF2>=p1
xF&&lH2
3]i6
3];}
if(IsLogicalValue(a)){if(iA1
iC2
iA1.when,p1
iE2
iA1.what;if(iC1
iC2
iC1.when,p1
iE2
iC1.what;}
if(IsLogicalValue(b)){if(iB1
iC2
iB1.when,p0
iE2
iB1.what;if(iD1
iC2
iD1.when,p0
iE2
iD1.what;}
return
xK;}
eY
static
bool
TestCase(lG2
when,const
range
yB&p){if(!p
cU2||!p
e3
cI
switch(when
nB3
Eq0
nO1==eB1(0.0)&&p
xF==p
iA2;case
Eq1
nO1==eB1(1.0)&&p
xF==p
xF;case
y63
nO1>x71&&p
xF<=eB1(1);case
y73
nO1>=x71
yS1
1);y83;}
return
eC3}
;l13
RangeComparisonsData{static
const
l4
Data[6]={{l4
lO3
i1
xK,l4::i1
xK
lO1
Eq1
lP1
Eq1
x51
Eq0
x61
Eq0}
}
,{l4::nM2
lP3
xK,l4
lP3
xK
lO1
Eq0
lP1
Eq0
x51
Eq1
x61
Eq1}
}
,{l4::nM2
lP3
t62,l4::i1
MakeFalse
x51
y63
lP1
y73
c3,{l4
lO3
xK,l4
lP3
i1
lQ3
x51
y73
lP1
y63
c3,{l4::nM2::i1
i1
MakeTrue,l4::t62
lO1
y73
x61
y63
c3,{l4
lO3
i1
lQ3,l4::xK,l4
n91
lO1
y63
x61
y73
c3}
;}
eY
bool
ConstantFolding_Comparison(cS2{using
l13
RangeComparisonsData;assert(tree.GetOpcode()>=cEqual&&tree.GetOpcode()<=cGreaterOrEq);switch(Data[i63-cEqual].Analyze(y1
0),y1
1))nB3
l4::MakeFalse:xB
0
lA3
n91:xB
1
lA3::lQ3:e93
cEqual
lA3::t62:e93
tU1
lA3::MakeNotNotP0:e93
cNotNot
iE1
1
lA3::MakeNotNotP1:e93
cNotNot
iE1
0
lA3::MakeNotP0:e93
cNot
iE1
1
lA3::MakeNotP1:e93
cNot
iE1
0
lA3::xK:;}
if(y1
1)l92)switch(y1
0)nF
nB3
cAsin:lN
fp_sin(iH2
cAcos:lN
fp_cos(e0))yH
eO2
x31
cLess?cGreater:x31
cLessOrEq?cGreaterOrEq:x31
cGreater?cLess:x31
cGreaterOrEq?cLessOrEq:i63);nZ
cAtan:lN
fp_tan(iH2
cLog:lN
fp_exp(iH2
cSinh:lN
fp_asinh(iH2
cTanh:if(fp_less(fp_abs(e0),xW3{lN
fp_atanh(e0))lA2
break;y83
t82
return
eC3}
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
l13
FUNCTIONPARSERTYPES;l13{
#ifdef DEBUG_SUBSTITUTIONS
yM
double
d){union{double
d;uint_least64_t
h;}
t72
d=d;lK1
h
nP1
#ifdef FP_SUPPORT_FLOAT_TYPE
yM
float
f){union{float
f;uint_least32_t
h;}
t72
f=f;lK1
h
nP1
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
yM
long
double
ld){union{long
double
ld;cN2{uint_least64_t
a;x13
short
b;}
s;}
t72
ld=ld;lK1
s.b<<data.s.a
nP1
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
yM
long
ld){o<<"("
<<std::hex<<ld
nP1
#endif
#endif
}
tJ{lO
nH)){}
lO
const
eB1&i
l71
x73
nH
i)){data
xC
#ifdef FP_SUPPORT_CXX11_MOVE
lO
eB1&&i
l71
x73
nH
std::move(i))){data
xC
#endif
lO
x13
v
l71
VarTag
nH
iJ2,v))cW2
x42
o
l71
OpcodeTag
nH
o))cW2
x42
o,x13
f
l71
FuncOpcodeTag
nH
o,f))cW2
iW1
b
l71
CloneTag
nH*b.data)){}
eY
xE::~CodeTree(){}
lC
ReplaceWithImmed
cP1
i){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Replacing "
eL3*this);if(IsImmed())OutFloatHex(std::cout,GetImmed())xR1" with const value "
<<i;OutFloatHex(std::cout,i)xR1"\n"
;
#endif
data=new
x82
yB(i);}
eY
cN2
ParamComparer{tW2()lB3
a,iW1
b)const{if(a
n62!=b
n62)return
a
n62<b
n62
eX
a.GetHash()<b.GetHash();}
}
;xC3
x82
yB::Sort(yJ3
Opcode
nB3
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
tU1:std::sort(lC3
i12
lC3
end(),ParamComparer
yB());lD
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
break;y83
t82}
lC
AddParam
lB3
param){y2.push_back(param);}
lC
eC1
xE&param){y2.push_back(xE());y2.back().swap(param);}
lC
SetParam(cI3
which,iW1
b)nA1
which
iI2
y2[which]=b;}
lC
n71
cI3
which,xE&b)nA1
which
iI2
y2[which
t63
b);}
lC
AddParams
cQ1
nL){y2.insert(y2.end(),lI2.i12
lI2.end());}
lC
AddParamsMove(nL){cI3
endpos=y2
eZ3,added=lI2
eZ3;y2
xA3
endpos+added,xE());for(cI3
p=0;p<added;++p)y2[endpos+p
t63
lI2[p]);}
lC
AddParamsMove(nL,cI3
lJ2)nA1
lJ2
iI2
DelParam(lJ2);AddParamsMove(tE1}
lC
SetParams
cQ1
nL){n72
tmp(tE1
y2.swap(tmp);}
lC
SetParamsMove(nL){y2.swap(tE1
lI2.clear();}
#ifdef FP_SUPPORT_CXX11_MOVE
lC
SetParams(n72&&lI2){SetParamsMove(tE1}
#endif
lC
DelParam(cI3
index){n72&cJ3=y2;
#ifdef FP_SUPPORT_CXX11_MOVE
lC3
erase(lC3
begin()+index);
#else
cJ3[index].data=0;for(cI3
p=index;p+1<cJ3
i52
p)cJ3[p].data.UnsafeSetP(&*cJ3[p+1
iI2
cJ3[xX3-1].data.UnsafeSetP(0);lC3
resize(xX3-1);
#endif
}
lC
DelParams(){y2.clear(iW2
xE::IsIdenticalTo
lB3
b)const{if(&*data==&*b.data)return
true
eX
data->IsIdenticalTo(*b.data
iW2
x82
yB::IsIdenticalTo
cQ1
x82
yB&b)const{if(Hash!=b.Hash
cI
if(Opcode!=t13
cI
switch(Opcode
nB3
cImmed:return
lK3
Value,t23;case
iJ2:return
iU1==b.iT1
case
cFCall:case
cPCall:if(iU1!=b.iU1
cI
break;y83
t82
if(xX3!=b.xX3
cI
eT2
cJ3
i52
a){if(!cJ3[a]xI
b.cJ3[a])cI}
return
true;}
lC
Become
lB3
b){if(&b!=this&&&*data!=&*b.data){DataP
tmp=b.data;iY
data.swap(tmp);}
}
lC
CopyOnWrite(){if(GetRefCount()>1)data=new
x82
yB(*data);}
eY
xE
xE::GetUniqueRef(){if(GetRefCount()>1)return
xE(*this,CloneTag())eX*this;}
eY
nD):c0
cNop),Value(),n9
eY
nD
const
x82&b):c0
t13),Value(t23,iU1(b.cE1,cJ3(b.cJ3),Hash(b.Hash),Depth(b.Depth),tO1
b.iV1){}
eY
nD
const
eB1&i):c0
cImmed),Value(i),n9
#ifdef FP_SUPPORT_CXX11_MOVE
eY
nD
x82
yB&&b):c0
t13),Value
cI2
t23),iU1(b.cE1,cJ3
cI2
b.cJ3)),Hash(b.Hash),Depth(b.Depth),tO1
b.iV1){}
eY
nD
eB1&&i):c0
cImmed),Value
cI2
i)),n9
#endif
eY
nD
x42
o):c0
o),Value(),n9
eY
nD
x42
o,x13
f):c0
o),Value(),iU1(f),cJ3(),Hash(),Depth(1),tO1
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
l13
FUNCTIONPARSERTYPES;
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
l13{xC3
tP1
nT,std
cC&done,std::ostream&o){eT2
yI++a)tP1
y1
a),done,o);std::ostringstream
buf
eL3
tree,buf);done[xE3].insert(buf.str());}
}
#endif
tJ{
#ifdef FUNCTIONPARSER_SUPPORT_DEBUGGING
xC3
DumpHashes(i8){std
cC
done;tP1
tree,done,o);for(std
cC::const_iterator
i=done.xU3
i!=done.end();++i){const
std::set<n13>&flist=i
cY2;if(flist
eZ3!=1)o<<"ERROR - HASH COLLISION?\n"
;for(std::set<n13>::const_iterator
j=flist.xU3
j!=flist.end();++j){o<<'['<<std::hex<<i->first.hash1<<','<<i->first.hash2<<']'<<std::dec;o<<": "
<<*j<<"\n"
;}
}
}
xC3
DumpTree(i8){iU2
tS3;switch(i63
nB3
cImmed:o<<tree
tJ1
eX;case
iJ2:o<<"Var"
<<(tree.GetVar()-iJ2)eX;case
cAdd:tS3"+"
;lD
cMul:tS3"*"
;lD
cAnd:tS3"&"
;lD
cOr:tS3"|"
;lD
cPow:tS3"^"
;break;y83
tS3;o<<e73
i63);if(x31
cFCall||x31
cPCall)o<<':'<<tree.GetFuncNo();}
o<<'(';if(i93<=1&&sep2[1])o<<(sep2+1)<<' '
nD2
if(a>0)o<<' '
eL3
y1
a),o);if(a+1<i93)o<<sep2;}
o<<')';}
xC3
DumpTreeWithIndent(i8,const
n13&indent){o<<'['<<std::hex<<(void*)(&tree.tF1)<<std::dec<<','<<tree.GetRefCount()<<']';o<<indent<<'_';switch(i63
nB3
cImmed:o<<"cImmed "
<<tree
tJ1;o<<'\n'
eX;case
iJ2:o<<"VarBegin "
<<(tree.GetVar()-iJ2);o<<'\n'
eX;y83
o<<e73
i63);if(x31
cFCall||x31
cPCall)o<<':'<<tree.GetFuncNo();o<<'\n';}
eT2
yI++a){n13
ind=indent;for(cI3
p=0;p<ind
eZ3;p+=2)if(ind[p]=='\\')ind[p]=' ';ind+=(a+1<i93)?" |"
:" \\"
;DumpTreeWithIndent(y1
a),o,ind);}
o<<std::flush;}
#endif
}
#endif
using
l13
l51;using
l13
FUNCTIONPARSERTYPES;
#include <cctype>
l13
l51{x13
ParamSpec_GetDepCode
cQ1
cQ2&b
yJ3
b
t03
nB3
nC3:{cR*s=(cR*)b.second
eX
s->depcode
t11
SubFunction:{cS*s=(cS*)b.second
eX
s->depcode;}
y83
t82
return
0;}
xC3
DumpParam
cQ1
cQ2&y42
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
;x13
xS2
0;switch(y52
t03
nB3
eM3{const
ParamSpec_NumConstant
yB
iR2
cQ1
ParamSpec_NumConstant
yB*e51
using
l13
FUNCTIONPARSERTYPES;o.precision(12);o<<iS2
constvalue;break
t11
nC3:{cR
iR2(cR*e51
o<<ParamHolderNames[iS2
index];xS2
iS2
constraints;break
t11
SubFunction:{cS
iR2(cS*e51
xS2
iS2
constraints;yP
GroupFunction){if(iS2
l81==cNeg){o<<"-"
;n2}
l72
iS2
l81==cInv){o<<"/"
;n2}
else{n13
opcode=e73(x42)iS2
l81).substr(1)l01
opcode
i52
a)opcode[a]=(char)std::toupper(opcode[a]);o<<opcode<<"( "
;n2
o<<" )"
;}
}
else{o<<'('<<e73(x42)iS2
l81)<<' ';yP
PositionalParams)o<<'[';yP
SelectedParams)o<<'{';n2
if
e92.n1!=0)o<<" <"
<<iS2
data.n1<<'>';yP
PositionalParams)o<<"]"
;yP
SelectedParams)o<<"}"
;o<<')';}
t82
iQ2
ImmedConstraint_Value
n23
ValueMask)nB3
ValueMask:lD
Value_AnyNum:lD
nN2:o<<"@E"
;lD
Value_OddInt:o<<"@O"
;lD
tR1:o<<"@I"
;lD
Value_NonInteger:o<<"@F"
;lD
eN1:o<<"@L"
;i62
ImmedConstraint_Sign
n23
SignMask)nB3
SignMask:lD
Sign_AnySign:lD
nB1:o<<"@P"
;lD
eO1:o<<"@N"
;i62
ImmedConstraint_Oneness
n23
OnenessMask)nB3
OnenessMask:lD
Oneness_Any:lD
Oneness_One:o<<"@1"
;lD
Oneness_NotOne:o<<"@M"
;i62
ImmedConstraint_Constness
n23
ConstnessMask)nB3
ConstnessMask:lD
tQ1:if(y52
t03==nC3){cR
iR2(cR*e51
if(iS2
index<2)t82
o<<"@C"
;lD
Constness_NotConst:o<<"@V"
;lD
Oneness_Any:t82}
xC3
DumpParams
nI2
paramlist,x13
count,std::ostream&o){for(eQ1=0;a<count;++a){if(a>0)o<<' ';const
cQ2&param=e81
yB(paramlist,a);DumpParam
yB(param,o);x13
depcode=ParamSpec_GetDepCode(param);if(depcode!=0)o<<"@D"
<<depcode;}
}
}
#include <algorithm>
using
l13
l51;using
l13
FUNCTIONPARSERTYPES;l13{cR
plist_p[37]={{2,0,0x0}
xG
0,0x4}
xG
nB1,0x0}
xG
eO1|Constness_NotConst,0x0}
xG
Sign_NoIdea,0x0}
xG
eN1,0x0}
,{3,Sign_NoIdea,0x0}
,{3,0,0x0}
,{3,eN1,0x0}
,{3,0,0x8}
,{3,Value_OddInt,0x0}
,{3,Value_NonInteger,0x0}
,{3,nN2,0x0}
,{3,nB1,0x0}
,{0,eO1|lW{0,lW{0,nB1|lW{0,nN2|lW{0,tQ1,0x1}
,{0,tR1|nB1|lW{0,tS1
tQ1,0x1}
,{0,tS1
lW{0,Oneness_One|lW{0,eN1|lW{1,lW{1,nN2|lW{1,tS1
lW{1,tR1|lW{1,nB1|lW{1,eO1|lW{6,0,0x0}
,{4,0,0x0}
,{4,tR1,0x0}
,{4,lW{4,0,0x16}
,{5,0,0x0}
,{5,lW}
;eY
cN2
plist_n_container{static
const
ParamSpec_NumConstant
yB
plist_n[20];}
;eY
const
ParamSpec_NumConstant
yB
plist_n_container
yB::plist_n[20]={{eB1(-2
i2-1
i2-0.5
i2-0.25
i2
0
t92
fp_const_deg_to_rad
yB(t92
fp_const_einv
yB(t92
fp_const_log10inv
yB(i2
0.5
t92
fp_const_log2
yB(i2
1
t92
fp_const_log2inv
yB(i2
2
t92
fp_const_log10
yB(t92
fp_const_e
yB(t92
fp_const_rad_to_deg
yB(t92-fp_const_pihalf
yB(),xJ1{x71,xJ1{fp_const_pihalf
yB(),xJ1{fp_const_pi
yB(),xJ1}
;cS
plist_s[517]={{{1,15,lR3,398,lR3,477,lR3,15,cNeg,GroupFunction,0}
,tQ1,0x1
eQ3
15,xT2
24,xT2
465,xT2
466,xT2
498,cInv,lU
2,327995
yK
e62
48276
yK
l6
260151
yK
l6
470171
yK
l6
169126
yK
l6
48418
yK
l6
1328
yK
l6
283962
yK
l6
169275
yK
l6
39202
yK
l6
283964
yK
l6
283973
yK
l6
476619
yK
l6
296998
yK
l6
47
yK
SelectedParams,0}
,0,0x4
nO
161839
yK
l6
25036
yK
l6
35847
yK
l6
60440
yK
l6
30751
yK
l6
183474
yK
l6
259318
yK
l6
270599
yK
l6
60431
yK
l6
259119
yK
l6
332066
yK
l6
7168
yK
l6
197632
yK
l6
291840
yK
l6
283648
yK
l6
238866
yK
l6
239902
yK
l6
31751
yK
l6
244743
yK
l6
384022
yK
SelectedParams,0}
,0,0x4
nO
385262
yK
l6
386086
yK
l6
393254
yK
SelectedParams,0}
,0,0x5
nO
393254
yK
l6
386095
yK
l6
387312
yK
l6
18662
yK
l6
61670
yK
l6
387397
yK
l6
247855
yK
SelectedParams,0}
,0,0x1
nO
342063
yK
l6
297007
yK
l6
15820
yK
l6
393263
yK
l6
393263
yK
SelectedParams,0}
,0,0x5
nO
161847
yK
l6
258103
yK
l6
249073
yK
l6
249076
yK
iE
0,0
yK
nM
0,0
tT1
1,45
yK
nM
1,53
yK
nM
1,54
yK
nM
1,55
yK
nM
1,56
yK
nM
1,26
yK
nM
1,259
e42
0x16
eQ3
272
tT1
1,323
e42
0x16
eQ3
0
yK
nM
1,21
yK
nM
1,447
e42
0x4
eQ3
449
e42
0x4
eQ3
0
e42
0x4
eQ3
0
tB
2}
,0,0x4
eQ3
15
yK
nM
1,24
tB
2}
,0,0x0
nO
58392
tT1
0,0
tB
1}
,nB1,0x0
nO
24591
eV3,33807
eV3,48143
eV3,285720
eV3,290840
eV3,305152,lA
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
321551,xK1
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
296976,tI1
324623,l1
0x14
nO
332815,l1
0x10}
,{{3,7340056,tI1
289092,lA
93200,xK1
337935
t81
7340060,l1
tA2
7340176,lA
338959
t81
7340061,xK1
7206,lA
7168,lA
357414,lA
368678,lA
370745,l1
0x7}
,{{3,7340177,lA
39277,tI1
426398,l1
tA2
40272286,xK1
490910,l1
tA2
40336798,xK1
50600,lA
426462,xK1
490974,xK1
370726,l1
0x6
nO
371750,l1
0x6
nO
428070
t81
40336862,xK1
38378,lA
50671
t81
47662080,lA
477184,lA
568320,lA
371727,l1
0x7}
,{{3,15779306,lA
370703,l1
0x7
nO
39277,lA
39279,l1
0x4}
,{{3,15779238,lA
39338,tI1
436262,lA
508966,lA
39409,tI1
296998,tI1
35847,lA
15,tI1
377894,lA
386063,l1
0x1
nO
15,lA
7192,lA
123928,lA
122904,lA
30751,lA
57,lA
7456,lA
15674
t81
67579935,lA
39237,lA
58768,lA
62924,lA
122880,lA
15760
t81
64009216,l1
0x0}
,{{0,0,xN
0,0,iO
2,cY1
2,cZ1
3,cY1
3,cZ1
38,xN
1,38,iO
14,xN
1,57,xN
1,16,e52
0x0
nO
471103,e52
0x1
eQ3
303,xN
1,323,yA3
0x0
nO
471363,e52
0x16
eQ3
293,cY1
294,cZ1
295,xN
1,296,iO
400,xN
1,0,xN
1,460,xN
1,465,xN
1,16,e52
0x1
eQ3
57,yA3
0x1
eQ3
0,iO
21,xN
1,15,e52
0x0
nO
24591,xN
1,24,iO
517,yA3
0x0
nO
46095,lL
46104,lL
15397,lL
287789,lL
66584,lL
404763,lL
62504,lL
15409,lL
39951,lL
24591,lL
33807,lL
50200,lL
62509,lL
50176,lH,178176,y6
eI3
283648,lH,19456,lH,27648,lH,91136,lH,86016,lH,488448,lH,14342,lH,58375,lH,46147
cD
46151,lH,284679,lH,7183,lH,46159
cD
38993
cD
50262,lH,50249,lH,283808,lH,284835,lH,24822,lH,10240,lH,11264,lH,7170,lH,7168,lH,17408,lH,164864,lH,237568,lH,242688,y6
0x14
nO
476160,lH,25607,lH,122895,lH,50252,lH,39374,lH,50183,lH,7192,lH,122911,lH,252979,lH,46155,lH,38919,lH,50268,lH,50269,lH,50253,lH,46191,lH,50296,lH,7563,y6
0x10
nO
416811,lH,416819,lH,40047,lH,46192
cD
415795,lH,40048
cD
415787,lH,39016
y53
39326
cD
39326,lH,39332
y53
39333,y6
0x1
nO
50590
cD
50590,lH,39338
cD
39338,lH,39335
y53
15786
cD
146858,lH,39372,lH,39379,lH,39380,lH,39390
cD
50654
cD
50654,lH,24,y6
0x6
nO
62,lH,24,lH,62,y6
0x6
nO
43,lH,43
cD
51,lH,51
cD
50270,lH,50176
cD
50271,lH,39159,lH,39183
cD
7168
cD
31744,lH,100352,lH,31746,lH,101400,lH,39409
cD
39411
cD
39411,lH,39420,lH,39420
cD
15,lH,39026
y53
39422,lH,16384,lH,62853,lH,15360,lH,15,y6
0x1
nO
16,lH,7183,y6
0x1
nO
7172,cPow,yA1,nB1,0x0
nO
24591,cPow,lU
2,50200,cPow,lU
2,63521,cPow,lU
2,62500,cPow,lU
2,50453,cPow,lU
2,62488,cPow,lU
1,0,eE3
7,eE3
194,eE3
0,cAcos,eF3
cAcosh,eF3
cAsin,eF3
cAsinh,nW
120,cAsinh,eF3
cAtan,e62
306176,cAtan2
cB3,cAtan2,eF3
cAtanh,nW
246,cCeil,eF3
cCeil,eG3
yY2
0,cCos,t9
1,7,yY2
92,yY2
93,yY2
120,yY2
236,yY2
255,yY2
214,iK2
236,iK2
464,iK2
0,cCosh,eG3
iK2
0,cExp,nW
7,cExp,nW
92,cExp,eF3
yB3
7,yB3
92,yB3
246,cFloor,eF3
cFloor,lS3
309540,eH3
e62
316708,eH3
e62
316724,eH3
l0
3,32513024,e72
34627584
lU1
31493120,e72
89213952
lU1
149042176
lU1
246647808
lU1
301234176
lU1
494360576
lU1
498558976
lU1
62933520
lU1
62933520,e72
62933526
lU1
62933526,e72
24670208
lU1
579378176
lU1
573578240
lU1
32513024
lU1
566254592
lU1
7900160
lU1
588822528,cIf,nW
120,cInt,nW
246
eB2
0
eB2
7
eB2
31
eB2
194
eB2
363
eB2
15,cLog,lU
1,24,cLog,lU
1,0,cLog10,eF3
cLog2
cB3,cMax,e62
35847,cMax,e62
30751,cMax,eF3
cMax,AnyParams,1}
,0,0x4
nO
7168,cMin,e62
35847,cMin,e62
30751,cMin,eF3
cMin,AnyParams,1}
,0,0x4
nO
24591,cMin,lU
1,0,nO2
7,nO2
92,nO2
93,nO2
120,nO2
149,nO2
231,cSin,lB
0x5
eQ3
246,nO2
255,nO2
254,nO2
0,cSin,t9
1,273,cSin,lB
0x1
eQ3
214,xU2
231,cSinh,lB
0x5
eQ3
246,xU2
254,xU2
255,xU2
464,xU2
0,cSinh,eG3
xU2
15,cSqrt,lU
1,0,yZ2
0,cTan,t9
1,116,cTan,t9
1,117,yZ2
231,yZ2
246,yZ2
273,yZ2
254,yZ2
255,yZ2
0,xV2
0,cTanh,t9
1,213,xV2
231,xV2
246,xV2
254,xV2
255,xV2
0,cTrunc,e62
15384,cSub,lU
2,15384,cDiv,lU
2,476626,cDiv,lU
2,122937,xW2
7168
e82,lB
eI3
7168,xW2
31744
e82,lB
0x20
nO
31751
e82,lB
0x24
nO
31751,xW2
122937,tU1
cB3,cLess,lB
eI3
41984,cLess,lS3
41984,cLess,e62
7,cLess
cB3,cLessOrEq,e62
296182,cLessOrEq
cB3
t51,lB
eI3
41984
t51,lS3
41984
t51,e62
7
t51
cB3
e9,e62
296182
e9,eF3
lK2
245,lK2
7,lK2
550,lK2
553,lK2
554,lK2
556,lK2
31,lK2
559,lK2
15,lK2
560,cNot,e62
7706,y23
7168,y23
35847,y23
30751,y23
463903,y23
466975,cAnd,iE
0,0,cAnd,nM
2,7168,cN3
7706,cN3
35847,cN3
463903,cN3
466975,cN3
30751,cOr,iE
1,0,lL2
92,lL2
131,lL2
245,lL2
215,lL2
246,cDeg,nW
246,cRad
cB3,cAbsAnd,l6
7168,cAbsOr,iE
1,0,nR2,eF3
cAbsNotNot,l0
3,32513024,cAbsIf,lB
0x0}
,}
;}
l13
l51{const
Rule
grammar_rules[262]={{ProduceNewTree,17,1,0,{1,0,cAbs
eC2
409,{1,146,cAtan
eC2
403
xG
1324,cAtan2
eC2
405
xG
307201,cAtan2
yF
253174
xG
255224,cAtan2
yF
259324
xG
257274,cAtan2
eC2
152,{1,252,cCeil
l3
2,1,486,{1,68,xN1
482,{1,123,xN1
483,{1,125,xN1
151,{1,126,xN1
419,{1,124,xN1
0,{1,403,cCos,l2
2,1,246,{1,252,cCos,l2
18,1,0,{1,400,xN1
301,{1,406,cCosh,l2
2,1,246,{1,252,cCosh,l2
18,1,0,{1,400,cCosh
l3
2,1,458,{1,122,cFloor
eC2
150,{1,252,cFloor
l3
0,1,156,{3,7382016,eR
549,{3,8430592,eR
556,{3,8436736,eR
157,{3,42998784,eR
550,{3,42999808,eR
562,{3,43039744,eR
557,{3,49291264,eR
538,{3,49325056,eR
469,{3,1058318,eR
473,{3,1058324,eR
473,{3,9438734,eR
469,{3,9438740,cIf,l2
0,3,32542225,{3,36732434,cIf,l2
0,3,32542231,{3,36732440,cIf
yE3
573,{3,32513026,cIf
yE3
515,{3,455505423,cIf
yE3
515,{3,433506837,cIf
l3
2,1,78,{1,256,xX2
69,{1,258,xX2
404,{1,72,xX2
159,{1,147,cLog,l2
0,1,0
xG
487425,cMax,tW
16,1,445
xG
yC3
cMax,tW
0,1,0
xG
483329,cMin,tW
16,1,446
xG
yC3
cMin,c7
0,1,153
xG
24832
iL1
0,1,153
xG
25854
iL1
0,1,154
xG
130063,iL2
32055,iL2
32056,iL2
32057,iM2
166288
xG
32137,iL2
33082,iM2
7168
xG
12688,iM2
7434
xG
12553
iL1
2,1,435
xG
46146
iL1
2,1,436
xG
46154
iL1
2,1,437
xG
46150
iL1
2,1,169
xG
83983
iL1
2,1,168
xG
131106
iL1
2,1,175
xG
133154
nW2
476160
xG
471055
nW2
274432
xG
273423
nW2
251904
xG
266274
nW2
251904
xG
263186
iL1
2,1,171,{1,252,lQ1
421,{1,68,lQ1
151,{1,123,lQ1
419,{1,125,lQ1
170,{1,126,lQ1
482,{1,124,lQ1
0,{1,405,lQ1
172,{1,252,cSinh
l3
2,1,328,{1,404,cSinh
l3
2,1,173,{1,252,xY2
0,{1,408,xY2
176,{1,410,xY2
177,{1,252,cTanh,l2
0,1,442
xG
449551,lT3
441
xG
yC3
lT3
167
xG
268549,lT3
180
xG
276749,lT3
181
xG
276500,lU3
190770
xG
189622,lU3
194748
xG
193723,lU3
202943
xG
196795,lU3
59699
xG
298148,lU3
59714
xG
325815,lU3
59724
xG
343224
yK
c7
2,1,337,{1,333
tB
1
tR
336,{1,338
tB
1}
}
,{ReplaceParams,2,1,340
xG
1363
t3
342
xG
1365
t3
463
xG
472524
t3
47
xG
356711
t3
349
xG
200751
t3
360
xG
199727
t3
480
xG
207053
t3
481
xG
208077
t3
417
xG
211144
t3
209
xG
211145
t3
418
xG
215240
t3
212
xG
212329
t3
204
xG
373097
t3
211
xG
372944
t3
217
xG
201944
t3
221
xG
223448
t3
367
xG
508329
t3
219
xG
508126
t3
224
xG
225705
t3
223
xG
225776
t3
365
xG
230825
t3
426
xG
377057
t3
497
xG
377054
t3
497
xG
204201
t3
426
xG
375280
t3
224
xG
375006,l7
2,2,407781
xG
233698,l7
2,2,59763
xG
233842,lT3
372
xG
1397,c02
96
xG
24705,c02
97
xG
24708,c02
444
xG
449551,c02
443
xG
yC3
c02
101
xG
102774,c02
109
xG
107845,c02
106
xG
104773,l5
0,2,111631
xG
109893,l5
0,2,108559
xG
110917,lK
0
tR
113
xG
112658,cMul,SelectedParams,0
tR
567,{1,52,lK
1
tR
568,{1,42,lK
1}
}
,{ReplaceParams,2,1,467
xG
45516
tG1
356
xG
51555
tG1
468
xG
49612
tG1
357
xG
47459
tG1
429
xG
438699
tG1
432
xG
441774
tG1
486
xG
498726
tG1
494
xG
504870
tG1
382
xG
435579
tG1
497
xG
435709
tG1
426
xG
508287
tG1
414
xG
500092
tG1
499
xG
352744
tG1
345
xG
367092
tG1
381
xG
425318
tG1
478
xG
425460
tG1
47
xG
512501
tG1
505
xG
355817
tG1
47
xG
516598
tG1
507
xG
518182
tG1
508
xG
358896
tG1
351
xG
388605
tG1
511
xG
360939
tG1
503
xG
354788
tG1
514
xG
525350
tG1
510
xG
394343
tG1
386
xG
351347,l5
2,2,363004
xG
361968,l5
16,1,118
xG
1157,l5
16,1,119
xG
1158,l5
16,1,402
xG
411024,l5
16,2,58768
xG
1472,l5
16,2,15760
xG
1474,l5
17,1,0,{1,400,l5
17,1,57,{1,14,lK
0}
}
,{ProduceNewTree,4,1,538
xG
41
e82
l3
4,1,0
xG
5167
e82,cJ
41984
xG
409641
e82,cJ
iF
e82,cJ
tX
e82,cJ
tY
e82
cD1
24849
e82
yF
iG
e82
yF
lM2
281873
e82
yF
l91
e82
yF
nC1
e82
l3
4,1,562
xG
41,tU1
l3
4,1,538
xG
5167,tU1,cJ
41984
xG
409641,tU1,cJ
iF,tU1,cJ
tX,tU1,cJ
tY,tU1
cD1
24849,tU1
yF
iG,tU1
yF
lM2
281873,tU1
yF
l91,tU1
yF
nC1,tU1,cJ
iF,yD3
tX,yD3
tY,cLess
eC2
571
xG
46080,cLess
cD1
24832,cLess
yF
yN1,cLess
yF
iG,cLess
yF
lM2
281856,cLess
yF
xP1,cLess
yF
l91,cLess
yF
nC1,cLess
l3
20,1,562
xG
409641,yD3
iF,xZ2
tX,xZ2
tY,cLessOrEq
eC2
565
xG
409615
lA1
529654
xG
24832
lA1
yN1
lA1
iG
lA1
lM2
281856
lA1
xP1
lA1
l91
lA1
nC1,cLessOrEq
l3
20,1,562
xG
409647,xZ2
iF
t51,cJ
tX
t51,cJ
tY
t51
eC2
539
xG
409615
t51
cD1
24832
t51
yF
yN1
t51
yF
iG
t51
yF
lM2
281856
t51
yF
xP1
t51
yF
l91
t51
yF
nC1
t51
l3
20,1,538
xG
409647
t51,cJ
iF
e9,cJ
tX
e9,cJ
tY
e9
eC2
572
xG
46080
e9
cD1
24832
e9
yF
yN1
e9
yF
iG
e9
yF
lM2
281856
e9
yF
xP1
e9
yF
l91
e9
yF
nC1
e9
l3
20,1,538
xG
409641
e9
l3
4,1,519,{1,137,cNot
yE3
571,{1,2,cNot,l2
0,1,452
xG
yC3
xD
0,2,537097,{3,547892744,cAnd,c7
16,1,566,{1,5,cAnd,AnyParams,1}
}
,{ReplaceParams,16,1,569
xG
13314,xD
16,1,544
xG
553498,xD
16,1,546
xG
462369,xD
16,1,548
xG
466465,xD
0,1,457
xG
yC3
nI
570
xG
13314,nI
563
xG
8197,nI
541
xG
553498,nI
542
xG
462369,nI
543
xG
466465,nI
564
xG
143365,cOr,c7
4,1,525,{1,137,cNotNot
yE3
572,{1,2,cNotNot
l3
17,1,0,{1,0,cNotNot
eC2
537,{1,256,cAbsNotNot,c7
18,1,531,{1,254,cAbsNotNot,c7
0,1,572,{3,43039744,cAbsIf
l3
0,1,571,{3,49325056,cAbsIf
yE3
454,{3,32513586,cAbsIf,l2
16,3,32542225,{3,36732434,cAbsIf,yA1}
,}
;cN2
grammar_optimize_abslogical_type{y4
9
cY
grammar_optimize_abslogical_type
grammar_optimize_abslogical={9,{34,192,228,238,242,247,254,260,261}
}
;}
cN2
grammar_optimize_ignore_if_sideeffects_type{y4
59
cY
grammar_optimize_ignore_if_sideeffects_type
grammar_optimize_ignore_if_sideeffects={59,{0,20,21,22,23,24,25,26,cT
iF1
78,cU
cZ
cN2
grammar_optimize_nonshortcut_logical_evaluation_type{y4
56
cY
grammar_optimize_nonshortcut_logical_evaluation_type
grammar_optimize_nonshortcut_logical_evaluation={56,{0,25,cT
iF1
78,cU
241,243,244,245,246,248,249,250,251,252,253,255,256,257,258,259}
}
;}
cN2
grammar_optimize_recreate_type{y4
22
cY
grammar_optimize_recreate_type
grammar_optimize_recreate={22,{18,55,56,57,80,81,82,83,84,85,117,118,120,121,130,131,132,133,134,135,136,137}
}
;}
cN2
grammar_optimize_round1_type{y4
125
cY
grammar_optimize_round1_type
grammar_optimize_round1={125,{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,25,cT
37,38,iF1
45,46,47,48,49,50,51,52,53,54,58,59,60,61,62,63,64,65,66,67,68,69,70,71,78,79,80,81,82,83,84,85,86,87,88,93,94,95,96,97,98,99,100,101,117,118,119,120,121,122,123,124,125,126,127,128,129,138,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cZ
cN2
grammar_optimize_round2_type{y4
103
cY
grammar_optimize_round2_type
grammar_optimize_round2={103,{0,15,16,17,25,cT
39,40,iF1
45,46,47,48,49,50,51,52,53,54,59,60,72,73,78,79,86,87,88,89,90,91,92,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,119,122,123,124,125,126,127,128,139,159,160,161,162,163,164,165,166,167,168,169,178,179,180,200,204,212,216,224,236,237,239,240,cZ
cN2
grammar_optimize_round3_type{y4
79
cY
grammar_optimize_round3_type
grammar_optimize_round3={79,{74,75,76,77,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,170,171,172,173,174,175,176,177,181,182,183,184,185,186,187,188,189,190,191,193,194,195,196,197,198,199,201,202,203,205,206,207,208,209,210,211,213,214,215,217,218,219,220,221,222,223,225,226,227,229,230,231,232,233,234,235}
}
;}
cN2
grammar_optimize_round4_type{y4
12
cY
grammar_optimize_round4_type
grammar_optimize_round4={12,{18,55,56,57,130,131,132,133,134,135,136,137}
}
;}
cN2
grammar_optimize_shortcut_logical_evaluation_type{y4
53
cY
grammar_optimize_shortcut_logical_evaluation_type
grammar_optimize_shortcut_logical_evaluation={53,{0,25,cT
iF1
78,cU
cZ}
l13
l51{eY
cQ2
e81
nI2
paramlist,lC1){index=(paramlist>>(index*10))&1023;if(index>=57)return
cQ2(SubFunction,c12
plist_s[index-57]);if(index>=37)return
cQ2(NumConstant,c12
plist_n_container
yB::plist_n[index-37])eX
cQ2(nC3,c12
plist_p[index]);}
}
#ifdef FP_SUPPORT_OPTIMIZER
#include <stdio.h>
#include <algorithm>
#include <map>
#include <sstream>
using
l13
FUNCTIONPARSERTYPES;using
l13
l51;using
tJ;using
cO1;l13{nK1
It,typename
T,typename
Comp>eP1
MyEqualRange(It
first,It
last,const
T&val,Comp
comp){cI3
len=last-first;while(len>0){cI3
half=len/2;It
nU3(first);nU3+=half;if(comp(*nU3,val)){first=nU3;++first;len=len-half-1;}
l72
comp(val,*nU3)){len=half;}
else{It
left(first);{It&eD2=left;It
last2(nU3);cI3
len2=last2-eD2;while(len2>0){cI3
half2=len2/2;It
n33(eD2);n33+=half2;if(comp(*n33,val)){eD2=n33;++eD2;len2=len2-half2-1;}
else
len2=half2;}
}
first+=len;It
right(++nU3);{It&eD2=right;It&last2=first;cI3
len2=last2-eD2;while(len2>0){cI3
half2=len2/2;It
n33(eD2);n33+=half2;if(comp(val,*n33))len2=half2;else{eD2=n33;++eD2;len2=len2-half2-1;}
eG2
eP1(left,right);eG2
eP1(first,first);}
eY
cN2
OpcodeRuleCompare{tW2()lB3
tree,x13
y02)const{const
Rule&rule=grammar_rules[y02]eX
i63<rule
c32.subfunc_opcode;}
tW2()nI2
y02,xD3)const{const
Rule&rule=grammar_rules[y02]eX
rule
c32.subfunc_opcode<i63;}
}
;eY
bool
TestRuleAndApplyIfMatch
cQ1
eP2
xE&tree,bool
c8{MatchInfo
yB
info;lZ1
found(false,eB());if((rule.lB1
LogicalContextOnly)&&!c8{tL1
if(nE
IsIntType
yB::cX3){if(rule.lB1
NotForIntegers)tL1
c03
rule.lB1
OnlyForIntegers)tL1
if(nE
IsComplexType
yB::cX3){if(rule.lB1
NotForComplex)tL1
c03
rule.lB1
OnlyForComplex)tL1
for(;;){
#ifdef DEBUG_SUBSTITUTIONS
#endif
found=TestParams(rule
c32,tree,found.specs,info,true);if(found.found)break;if(!&*found.specs){fail:;
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch
c22,false);
#endif
return
eC3}
#ifdef DEBUG_SUBSTITUTIONS
DumpMatch
c22,true);
#endif
SynthesizeRule
c22
lA2}
cO1{eY
bool
ApplyGrammar
cQ1
Grammar&tB2,xE&tree,bool
c8{if(tree.GetOptimizedUsing()==&tB2){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Already optimized:  "
eL3
tree)xR1"\n"
<<std::flush;
#endif
return
eC3
if(true){bool
changed
e63
switch(i63
nB3
cNot:case
cNotNot:case
cAnd:case
cOr:eT2
tree.xA
true))nV2
lD
cIf:case
cAbsIf:if(ApplyGrammar(tB2,y1
0),x31
cIf))nV2
y33
1;a<tree.xA
c8)nV2
break;y83
eT2
tree.xA
false))nV2}
if(changed){tree.Mark_Incompletely_Hashed(lA2}
typedef
const
x13
short*lY3;std::pair<lY3,lY3>range=MyEqualRange(tB2.rule_list,tB2.rule_list+tB2.rule_count,tree,OpcodeRuleCompare
yB());std::vector<x13
short>rules;rules.xF3
range.second-range
t03);cV
if(IsLogisticallyPlausibleParamsMatch(e01
c32,tree))rules.push_back(*r);}
range
t03=!rules
i43?&rules[0]:0;range.second=!rules
i43?&rules[rules
eZ3-1]+1:0;if(range
t03!=range.second){
#ifdef DEBUG_SUBSTITUTIONS
if(range
t03!=range.second)y12"Input ("
<<e73
i63)<<")["
<<i93<<"]"
;if(c8
std::cout<<"(Logical)"
;x13
first=iG1,prev=iG1;iU2
sep=", rules "
;cV
if(first==iG1)first=prev=*r;l72*r==prev+1)prev=*r;else
y12
sep<<first;sep=","
;if(prev!=first)std::cout<<'-'<<prev;first=prev=*r;}
}
if(first!=iG1)y12
sep<<first;if(prev!=first)std::cout<<'-'<<prev;}
std::cout<<": "
eL3
tree)xR1"\n"
<<std::flush;}
#endif
bool
changed
e63
cV
#ifndef DEBUG_SUBSTITUTIONS
if(!IsLogisticallyPlausibleParamsMatch(e01
c32,tree))yC1
#endif
if(TestRuleAndApplyIfMatch(e01,tree,c8){nV2
t82}
if(changed){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Changed."
<<std::endl
xR1"Output: "
eL3
tree)xR1"\n"
<<std::flush;
#endif
tree.Mark_Incompletely_Hashed(lA2}
tree.SetOptimizedUsing(&tB2)eX
eC3
eY
bool
ApplyGrammar
cQ1
void*p,FPoptimizer_CodeTree::cS2
yU
ApplyGrammar(*cQ1
Grammar*)p,tree);}
xC3
ApplyGrammars(FPoptimizer_CodeTree::cS2{
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_round1\n"
;
#endif
n6
grammar_optimize_round1
iN2
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_round2\n"
;
#endif
n6
grammar_optimize_round2
iN2
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_round3\n"
;
#endif
n6
grammar_optimize_round3
iN2
#ifndef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_nonshortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_nonshortcut_logical_evaluation
iN2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_round4\n"
;
#endif
n6
grammar_optimize_round4
iN2
#ifdef FP_ENABLE_SHORTCUT_LOGICAL_EVALUATION
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_shortcut_logical_evaluation\n"
;
#endif
n6
grammar_optimize_shortcut_logical_evaluation
iN2
#endif
#ifdef FP_ENABLE_IGNORE_IF_SIDEEFFECTS
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_ignore_if_sideeffects\n"
;
#endif
n6
grammar_optimize_ignore_if_sideeffects
iN2
#endif
#ifdef DEBUG_SUBSTITUTIONS
std
tX3"grammar_optimize_abslogical\n"
;
#endif
n6
grammar_optimize_abslogical
iN2
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
l13
FUNCTIONPARSERTYPES;using
l13
l51;using
tJ;using
cO1;l13{eY
bool
TestImmedConstraints
nI2
bitmask,xD3
yJ3
bitmask&ValueMask
nB3
Value_AnyNum:case
ValueMask:lD
nN2:if(GetEvennessInfo(tree)!=lN2
Value_OddInt:if(GetEvennessInfo(tree)!=nP2
tR1:if(GetIntegerInfo(tree)!=lN2
Value_NonInteger:if(GetIntegerInfo(tree)!=nP2
eN1:if(!IsLogicalValue(tree)cI
nD1
SignMask
nB3
Sign_AnySign:lD
nB1:if(l21
lN2
eO1:if(l21
nP2
Sign_NoIdea:if(l21
Unknown
cI
nD1
OnenessMask
nB3
Oneness_Any:case
OnenessMask:lD
Oneness_One:if(!yR1
if(!lK3
fp_abs(tree
tJ1),eB1(1))cI
lD
Oneness_NotOne:if(!yR1
lJ3
fp_abs(tree
tJ1),eB1(1))cI
nD1
ConstnessMask
nB3
Constness_Any:lD
tQ1:if(!yR1
lD
Constness_NotConst:if(yR1
t82
return
true;}
tH1
x13
extent,x13
nbits,typename
eE2=x13
int>cN2
nbitmap{private:static
const
x13
bits_in_char=8;static
const
x13
eF2=(yF3
eE2)*bits_in_char)/nbits;eE2
data[(extent+eF2-1)/eF2];cY3
void
inc(lC1,int
by=1){data[pos(index)]+=by*eE2(1<<y22);xA1
void
dec(lC1){inc(index,-1);}
int
get(lC1
xC2(data[pos(index)]>>y22)&mask();}
static
yO1
pos(lC1)yU
index/eF2;}
static
yO1
shift(lC1)yU
nbits*(index%eF2);}
static
yO1
mask()yU(1<<nbits)-1;}
static
yO1
mask(lC1)yU
mask()<<y22;}
}
;cN2
cH3{int
SubTrees:8;int
Others:8;int
lT2:8;int
eW3:8;nbitmap<iJ2,2>SubTreesDetail;cH3(){std::memset(this,0,yF3*this));}
cH3
cQ1
cH3&b){std::memcpy(this,&b,yF3
b));}
cH3&eD1=cQ1
cH3&b){std::memcpy(this,&b,yF3
b))eX*this;}
}
;eY
cH3
CreateNeedList_uncached(eZ&tC2{cH3
cN1;for(eQ1=0;a<t33
y32;++a){const
cQ2&y52=e81
yB(t33.param_list,a);switch(y52
t03
nB3
SubFunction:{cS
iR2(cS*e51
yP
GroupFunction)++cN1.eW3;else{++i42;assert(param.data.subfunc_opcode<VarBegin);cN1.SubTreesDetail.inc(iS2
l81);}
++cN1.lT2;break
t11
eM3
case
nC3:++i32;++cN1.lT2;break;eG2
cN1;}
eY
cH3&CreateNeedList(eZ&tC2{typedef
std::map<eZ*,cH3>e11;static
e11
yP1;e11::xV3
i=yP1.xJ2&tC2;if(i!=yP1.cL1&tC2
return
i
cY2
eX
yP1.xR3,std::make_pair(&t33,CreateNeedList_uncached
yB(tC2))cY2;}
eY
xE
CalculateGroupFunction
cQ1
cQ2&y42
const
xZ
info
yJ3
y52
t03
nB3
eM3{const
ParamSpec_NumConstant
yB
iR2
cQ1
ParamSpec_NumConstant
yB*)y52.second
eX
CodeTreeImmed(iS2
constvalue)t11
nC3:{cR
iR2(cR*)y52.second
eX
info.GetParamHolderValueIfFound(iS2
index)t11
SubFunction:{cS
iR2(cS*e51
xE
cX3;cX3
yX
iS2
l81);cX3.tF1.reserve
e92
y32);for(eQ1=0;a<iS2
data
y32;++a
nQ
tmp(CalculateGroupFunction(e81
yB
e92.param_list,a),info));cX3
y41
tmp);}
cX3.Rehash()eX
cX3;eG2
xE();}
}
cO1{eY
bool
IsLogisticallyPlausibleParamsMatch(eZ&t33,xD3){cH3
cN1(CreateNeedList
yB(tC2);cI3
eJ3=yI
if(eJ3<cI3(cN1.lT2))lZ3}
eT2
eJ3;++a){x13
opcode=y1
a)nF;switch(opcode
nB3
cImmed:if(cN1.eW3>0)--cN1.eW3;else--i32;lD
iJ2:case
cFCall:case
cPCall:--i32;break;y83
assert(opcode<VarBegin);if(i42>0&&cN1.SubTreesDetail.get(opcode)>0){--i42;cN1.SubTreesDetail.dec(opcode);}
else--i32;}
}
if(cN1.eW3>0||i42>0||i32>0)lZ3}
if(t33.match_type!=AnyParams){if(0||i42<0||i32<0)lZ3
eG2
true;}
eY
lZ1
TestParam
cQ1
cQ2&y42
xD3
iT2
nZ3,xZ
info
yJ3
y52
t03
nB3
eM3{const
ParamSpec_NumConstant
yB
iR2
cQ1
ParamSpec_NumConstant
yB*e51
if(!yR1
eB1
imm=tree
tJ1;switch(iS2
modulo
nB3
Modulo_None:lD
Modulo_Radians:imm=c33
imm,yD
imm<x41
imm
yS
if(imm>fp_const_pi
yB())imm-=fp_const_twopi
yB();t82
return
lK3
imm,iS2
constvalue)t11
nC3:{cR
iR2(cR*e51
if(!x4
return
info.SaveOrTestParamHolder(iS2
index,tree)t11
SubFunction:{cS
iR2(cS*e51
yP
GroupFunction){if(!x4
xE
xQ1=CalculateGroupFunction(y42
info);
#ifdef DEBUG_SUBSTITUTIONS
DumpHashes(xQ1)xR1*cQ1
void**)&xQ1
tJ1
xR1"\n"
xR1*cQ1
void**)&tree
tJ1
xR1"\n"
;DumpHashes(tree)xR1"Comparing "
eL3
xQ1)xR1" and "
eL3
tree)xR1": "
xR1(xQ1
xI
tree)?"true"
:"false"
)xR1"\n"
;
#endif
return
xQ1
xI
tree);}
c03!&*nZ3){if(!x4
if(i63!=iS2
l81
cI}
return
TestParams
e92,tree,nZ3,info,false);}
eG2
eC3
eY
cN2
iT
x12
MatchInfo
yB
info;iT():nZ3(),info(){}
}
;x63
MatchPositionSpec_PositionalParams:xL1
iT
yB>{cY3
iO2
MatchPositionSpec_PositionalParams
n03):e21
iT
yB>(n){}
}
;cN2
iH1
x12
iH1():nZ3(){}
}
;class
yQ:xL1
iH1>{cY3
x13
trypos;iO2
yQ
n03):e21
iH1>(n),trypos(0){}
}
;eY
lZ1
TestParam_AnyWhere
cQ1
cQ2&y42
xD3
iT2
nZ3,xZ
info,eV2&used,bool
lB2{xR<yQ>x9;xA2
yQ*yJ2
a=x9->trypos;goto
retry_anywhere_2;}
cX2
yQ(eN3
a=0;}
eK3
yI++a){if(used[a])yC1
retry_anywhere
yK3
TestParam(y42
y1
a),yH3);yI3
used[a]=true;if(lB2
eO3
a);x9->trypos=a
eX
lZ1(true,&*x9);}
}
retry_anywhere_2
yL3
goto
retry_anywhere;eG2
eC3
eY
cN2
yD1
x12
MatchInfo
yB
info;eV2
used;iO2
yD1(cI3
eJ3):nZ3(),info(),used(eJ3){}
}
;x63
MatchPositionSpec_AnyParams:xL1
yD1
yB>{cY3
iO2
MatchPositionSpec_AnyParams
n03,cI3
m):e21
yD1
yB>(n,yD1
yB(m)){}
}
;eY
lZ1
TestParams(eZ&nS,xD3
iT2
nZ3,xZ
info,bool
lB2{if(nS.match_type!=AnyParams){if(y7!=i93
cI}
if(!IsLogisticallyPlausibleParamsMatch(nS,tree))lZ3
iQ2
nS.match_type
nB3
PositionalParams:{xR<yA>x9;xA2
yA*yJ2
a=y7-1;goto
lE1;}
cX2
yA(y7);a=0;}
eK3
y7;++a){eR1=info;retry_positionalparams
yK3
TestParam(e5
a),y1
a),yH3);yI3
yC1}
}
lE1
yL3
info=eR1;goto
retry_positionalparams;}
if(a>0){--a;goto
lE1;}
info=(*x9)yG3
eX
eC3
if(lB2
for(eQ1=0;a<y7;++a)eO3
a)eX
lZ1(true,&*x9)t11
SelectedParams:case
AnyParams:{xR<tK>x9;eV2
used(eN3
std::vector<x13>iP2(y7);std::vector<x13>y62(y7)t61{const
cQ2
y52=e5
a);iP2[a]=ParamSpec_GetDepCode(y52);}
{x13
b=0
t61
if(iP2[a]!=0)y62[b++]=a
t61
if(iP2[a]==0)y62[b++]=a;}
xA2
tK*yJ2
if(y7==0){a=0;goto
retry_anyparams_4;}
a=y7-1;goto
e61;}
cX2
tK(y7,eN3
a=0;if(y7!=0){(*x9)yG3=info;(*x9)[0].used=used;}
}
eK3
y7;++a){if(a>0){eR1=info;(*x9)[a].used=used;}
retry_anyparams
yK3
TestParam_AnyWhere
yB(e5
y62[a]),tree,yH3,used,lB2;yI3
yC1}
}
e61
yL3
info=eR1;used=(*x9)[a].used;goto
retry_anyparams;}
e71:if(a>0){--a;goto
e61;}
info=(*x9)yG3
eX
eC3
retry_anyparams_4:if(nS.n1!=0){if(!TopLevel||!info.HasRestHolder(nS.n1)){n72
c42;c42.xF3
eN3
for
nI2
b=0;b<yI++b){if(yM3)yC1
c42.push_back(y1
b));yM3=true;if(lB2
eO3
b);}
if(!info.SaveOrTestRestHolder(nS.n1,c42)){goto
e71;}
}
else{x93
c42=info.GetRestHolderValues(nS.n1)l01
c42
i52
a){bool
found
e63
for
nI2
b=0;b<yI++b){if(yM3)yC1
if(c42[a]xI
y1
b))){yM3=true;if(lB2
eO3
b);found=true;t82}
if(!found){goto
e71;}
}
eG2
lZ1(true,y7?&*x9:0)t11
GroupFunction:t82
return
eC3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
#include <algorithm>
#include <assert.h>
using
tJ;using
cO1;l13{eY
xE
xW1
const
cQ2&y42
xZ
info,bool
inner=true
yJ3
y52
t03
nB3
eM3{const
ParamSpec_NumConstant
yB
iR2
cQ1
ParamSpec_NumConstant
yB*)y52.second
eX
CodeTreeImmed(iS2
constvalue)t11
nC3:{cR
iR2(cR*)y52.second
eX
info.GetParamHolderValue(iS2
index)t11
SubFunction:{cS
iR2(cS*e51
xE
tree;e93
iS2
l81);for(eQ1=0;a<iS2
data
y32;++a
nQ
nparam=xW1
e81
yB
e92.param_list,a),info,true
yH
eC1
nparam);}
if
e92.n1!=0){n72
trees(info.GetRestHolderValues
e92.n1)yH
AddParamsMove(trees);if(i93==1){assert(tree.GetOpcode()==cAdd tN3()==cMul tN3()==cMin tN3()==cMax tN3()==cAnd tN3()==cOr tN3()==cAbsAnd tN3()==cAbsOr);tree.e12
0));}
l72
i93==0
yJ3
i63
nB3
cAdd:case
cOr:tree=n61
0));lD
cMul:case
cAnd:tree=n61
1));y83
t82}
}
if(inner)tree.Rehash()eX
tree;eG2
xE();}
}
cO1{xC3
SynthesizeRule
cQ1
eP2
xE&tree,xZ
info
yJ3
rule.ruletype
nB3
ProduceNewTree:{tree.Become(xW1
e81
lF1
0),info,false));break
t11
ReplaceParams:y83{std::vector<x13>list=info.GetMatchedParamIndexes();std::sort(list.i12
list.end());y33
list
eZ3;a-->0;)tree.DelParam(list[a]);for(eQ1=0;a<rule.repl_param_count;++a
nQ
nparam=xW1
e81
lF1
a),info,true
yH
eC1
nparam);}
t82}
}
}
#endif
#ifdef DEBUG_SUBSTITUTIONS
#include <sstream>
#include <cstring>
using
l13
FUNCTIONPARSERTYPES;using
l13
l51;using
tJ;using
cO1;l13
l51{xC3
DumpMatch
cQ1
eP2
xD3,const
xZ
info,bool
DidMatch,std::ostream&o){DumpMatch
c22,DidMatch?i33"match"
:i33"mismatch"
,o);}
xC3
DumpMatch
cQ1
eP2
xD3,const
xZ
info,iU2
t43,std::ostream&o){static
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
;o<<t43<<" (rule "
<<(&rule-grammar_rules)<<")"
<<":\n  Pattern    : "
;{cQ2
tmp;tmp
t03=SubFunction;ParamSpec_SubFunction
tmp2;tmp2.data=rule
c32;tmp.second=c12
tmp2;DumpParam
yB(tmp,o);}
o<<"\n  Replacement: "
;DumpParams
lF1
rule.repl_param_count
i83
o<<"  Tree       : "
eL3
tree
i83
if(!std::strcmp(t43,i33"match"
))DumpHashes(tree,o)l01
info.yJ
i52
a){if(!info.yJ[a]lD2)yC1
o<<"           "
<<ParamHolderNames[a]<<" = "
eL3
info.yJ[a]i83}
cV1
info.lR
i52
b){if(!eA2
first)continue
l01
eA2
second
i52
a){o<<"         <"
<<b<<"> = "
eL3
eA2
second[a]i73
std::endl;}
}
o<<std::flush;}
}
#endif
#include <list>
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
using
l13
FUNCTIONPARSERTYPES;l13{eY
bool
MarkIncompletes(FPoptimizer_CodeTree::cS2{if(tree.Is_Incompletely_Hashed(tA1;bool
iI1=false
l01
yI++a)iI1|=MarkIncompletes(y1
a));if(iI1)tree.Mark_Incompletely_Hashed()eX
iI1;}
xC3
FixIncompletes(FPoptimizer_CodeTree::cS2{if(tree.Is_Incompletely_Hashed()){eT2
yI++a)FixIncompletes(y1
a)yH
Rehash();}
}
}
tJ{lC
Sort()e03
Sort();}
lC
Rehash(bool
constantfolding){if(constantfolding)ConstantFolding(*this);else
Sort();data
xC
eY
cN2
cE
yN3
eB1
yO3
yE1=0;
#if 0
long
double
value=Value;eK=crc32::calc(cQ1
x13
char*)&value,yF3
value));key^=(key<<24);
#elif 0
union{cN2{x13
char
filler1[16];eB1
v;x13
char
filler2[16];}
buf2;cN2{x13
char
filler3[yF3
eB1)+16-yF3
x11)];eK;}
buf1;}
data;memset(&data,0,yF3
data));data.buf2.v=Value;eK=data.buf1.key;
#else
int
exponent;eB1
n82=std::frexp(Value,&tL2
eK=nI2
e83+0x8000)&0xFFFF);if(n82<0){n82=-n82;key=key^0xFFFF;}
else
key+=0x10000;n82-=eB1(0.5);key<<=39;key|=n11(n82+n82)*eB1(1u<<31))<<8;
#endif
lQ
#ifdef FP_SUPPORT_COMPLEX_NUMBERS
nK1
T
c52
std::complex<T> >yN3
std::complex<T>yO3
cE<T>::n43
NewHash,Value.real());nE
fphash_t
temp;cE<T>::n43
temp,Value.imag());yE1^=temp.hash2;NewHash.hash2^=temp.hash1;}
}
;
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
tH1
c52
long>yL
long
Value){eK=Value;lQ
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
tH1
c52
GmpInt>yN3
GmpInt
yO3
eK=Value.toInt();lQ
#endif
xC3
x82
yB::Recalculate_Hash_NoRecursion(){fphash_t
NewHash(n11
Opcode)<<56,Opcode*tR3(0x1131462E270012B));Depth=1;switch(Opcode
nB3
cImmed:{cE
yB::n43
NewHash,Value);break
t11
iJ2:{yE1|=n11
cE1<<48
cW1((n11
cE1)*11)^tR3(0x3A83A83A83A83A0);break
t11
cFCall:case
cPCall:{yE1|=n11
cE1<<48
cW1((~n11
cE1)*7)^3456789;}
y83{cI3
eS1=0
l01
cJ3
i52
a){if(cJ3[a]n62>eS1)eS1=cJ3[a]n62;yE1+=((cJ3[a].tD2
hash1*(a+1))>>12)cW1
cJ3[a].tD2
hash1
cW1(3)*tR3(0x9ABCD801357);NewHash.hash2*=tR3(0xECADB912345)cW1(~cJ3[a].tD2
hash2)^4567890;}
Depth+=eS1;}
}
if(Hash!=NewHash){Hash=NewHash;iV1=0;}
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
l13
FUNCTIONPARSERTYPES;l13{using
tJ;eY
bool
nR1
xD3,long
count,const
x91::SequenceOpCode
yB&t5,x91
eW2&synth,cI3
max_bytecode_grow_length);static
const
cN2
SinCosTanDataType{OPCODE
whichopcode;OPCODE
inverse_opcode;enum{nominator,n92,inverse_nominator,lG1}
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
,{c62{cSinh,cCosh,c72,{cSinh,cNop,{c62
cNop,cCosh}
}
,{cCosh,cNop,{cSinh,c62
cNop}
}
,{cNop,cTanh,{cCosh,cSinh,c72,{cNop,cSinh,{cNop,cTanh,cCosh,cNop}
}
,{cNop,cCosh,{cTanh,cSinh,c72}
;}
tJ{lC
SynthesizeByteCode(std::vector<x13>&eE1,std::vector
yB&Immed,cI3&stacktop_max){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Making bytecode for:\n"
;iS
#endif
while(RecreateInversionsAndNegations()){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"One change issued, produced:\n"
;iS
#endif
FixIncompleteHashes();using
cO1;using
l13
l51;const
void*g=c12
grammar_optimize_recreate;while(ApplyGrammar(*cQ1
Grammar*)g,*this)){FixIncompleteHashes();}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Actually synthesizing, after recreating inv/neg:\n"
;iS
#endif
x91
eW2
synth;SynthesizeByteCode(synth,false);nS3
Pull(eE1,Immed,stacktop_max);}
lC
SynthesizeByteCode(x91
eW2&synth,bool
MustPopTemps)const{xX1*this))yU;}
eT2
12;++a){const
SinCosTanDataType&data=SinCosTanData[a];if(data.whichopcode!=cNop)i41!=data.whichopcode)yC1
CodeTree
n53;n53.nX1);n53
yX
data.inverse_opcode);n53.y72);xX1
n53)){nS3
i9
else
i41!=cInv)yC1
if(GetParam(0)nF!=data.inverse_opcode)yC1
xX1
GetParam(0))){nS3
i9
cI3
found[4];cV1
4;++b){CodeTree
tree;if(data.t53]==cNop){e93
cInv);CodeTree
n63;n63.nX1);n63
yX
data.t53^2]);n63.y72
yH
eC1
n63);}
else{tree.nX1
yH
eO2
data.t53]);}
tree.y72);found[b]=nS3
xH3
tree);}
if(found[data.y82!=tS
n92]yY
y82
nA2
n92
i7
cDiv
iA
y82!=tS
lG1]yY
y82
nA2
lG1
i7
cMul
iA
lR1!=tS
lG1]yY
lR1
nA2
lG1
i7
cRDiv
iA
lR1!=tS
n92]yY
lR1
nA2
n92
i7
cMul,2,1);nS3
i9
cI3
n_subexpressions_synthesized=SynthCommonSubExpressions(synth);switch(iJ1{case
iJ2:nS3
PushVar(GetVar());lD
cImmed:tD1
GetImmed());lD
cAdd:case
cMul:case
cMin:case
cMax:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:i41==cMul){bool
yP3
e63
yR
lS1
l92&&isLongInteger(lS1
tJ1)){yY1=makeLongInteger(lS1
tJ1);CodeTree
tmp(*this,typename
CodeTree::CloneTag());tmp.cT2
tmp
lS2
if(nR1
tmp,value,x91::tM1
yB::AddSequence,synth,MAX_MULI_BYTECODE_LENGTH)){yP3=true;t82}
}
if(yP3)t82
int
yF1=0;eV2
done(yQ1,false);CodeTree
iH;iH
yX
iJ1;for(;;){bool
found
e63
yR
done[a])yC1
if(nS3
IsStackTop(lS1)){found=true;done[a]=true;lS1.n8
iH
cL
lS1);if(++yF1>1){yW
2);iH.y72
tC1
iH);yF1=yF1-2+1;}
}
}
if(!found)t82
yR
done[a])yC1
lS1.n8
iH
cL
lS1);if(++yF1>1){yW
2);iH.y72
tC1
iH);yF1=yF1-2+1;}
}
if(yF1==0
yJ3
iJ1{case
cAdd:case
cOr:case
cAbsOr:tD1
0);lD
cMul:case
cAnd:case
cAbsAnd:tD1
1);lD
cMin:case
cMax:tD1
0);break;y83
t82++yF1;}
assert(n_stacked==1);break
t11
cPow:{l33
p0
nX2
0);l33
p1
nX2
1);if(!p1
l92||!isLongInteger
xM1)||!nR1
p0,makeLongInteger
xM1),x91::tM1
yB::MulSequence,synth,MAX_POWI_BYTECODE_LENGTH)){p0.n8
p1.n8
yW
2);xZ1
cIf:case
cAbsIf:{typename
x91
eW2::IfData
nT3;GetParam(0)tE2
SynthIfStep1(nT3,iJ1;GetParam(1)tE2
SynthIfStep2(nT3);GetParam(2)tE2
SynthIfStep3(nT3);break
t11
cFCall:case
cPCall:{eT2
yQ1;++a)lS1.n8
yW
nI2)yQ1);nZ2
0x80000000u|GetFuncNo(),0,0);t82
y83{eT2
yQ1;++a)lS1.n8
yW
nI2)yQ1);t82}
nS3
StackTopIs(*this);if(MustPopTemps&&n_subexpressions_synthesized>0){cI3
top=nS3
GetStackTop();nS3
DoPopNMov(top-1-n_subexpressions_synthesized,top-1);}
}
}
l13{eY
bool
nR1
xD3,long
count,const
x91::SequenceOpCode
yB&t5,x91
eW2&synth,cI3
max_bytecode_grow_length){if
e53!=0){x91
eW2
backup=synth;tree.n8
cI3
bytecodesize_backup=nS3
GetByteCodeSize();x91::nR1
count
nY2
cI3
bytecode_grow_amount=nS3
GetByteCodeSize()-bytecodesize_backup;if(bytecode_grow_amount>max_bytecode_grow_length){synth=backup
eX
eC3
return
true;}
else{x91::nR1
count,t5,synth
lA2}
}
#endif
#include <cmath>
#include <cassert>
#ifdef FP_SUPPORT_OPTIMIZER
using
l13
FUNCTIONPARSERTYPES;l13{using
tJ;
#define FactorStack std::vector
const
cN2
PowiMuliType{x13
opcode_square;x13
opcode_cumulate;x13
opcode_invert;x13
opcode_half;x13
opcode_invhalf;}
iseq_powi={cSqr,cMul,cInv,cSqrt,cRSqrt}
,iseq_muli={iG1
yK
cNeg,iG1,iG1}
;eY
eB1
cF1
const
PowiMuliType&yQ3,cM1,lO2&stack){eB1
cQ3
1);while(IP<limit){if(n73
yQ3.opcode_square){if(!t32
cR3
2;c4
opcode_invert){cX3=-cX3;c4
opcode_half){if(cX3>x71&&isEvenInteger(cR3
eB1(0.5);c4
opcode_invhalf){if(cX3>x71&&isEvenInteger(cR3
eB1(-0.5);++IP;yC1}
cI3
nQ2=IP;eB1
lhs(1);if(n73
cFetch){lC1=nA3;if(index<y3||cI3(index-y3)>=iO1){IP=nQ2;t82
lhs=stack[index-y3];goto
y92;}
if(n73
cDup){lhs=cX3;goto
y92;y92:eT1
cX3);++IP;eB1
subexponent=cF1
yQ3
tZ
if(IP>=limit||eE1[IP]!=yQ3.opcode_cumulate){IP=nQ2;t82++IP;stack.pop_back();cX3+=lhs*subexponent;yC1}
t82
return
cX3;}
eY
eB1
ParsePowiSequence(cM1){lO2
stack;eT1
eB1(1))eX
cF1
iseq_powi
tZ}
eY
eB1
ParseMuliSequence(cM1){lO2
stack;eT1
eB1(1))eX
cF1
iseq_muli
tZ}
x63
CodeTreeParserData{cY3
iO2
CodeTreeParserData(bool
k_powi):stack(),clones(),keep_powi(k_powi){}
void
Eat(cI3
eJ3,OPCODE
opcode
nQ
xO;xO
yX
opcode);n72
t33=Pop(eJ3);xO
nB2
tC2;if(!keep_powi)switch(opcode
nB3
cTanh:{xE
sinh,cosh;sinh
yX
cSinh);sinh
cL
xO
yR3
sinh
lS2
cosh
yX
cCosh);cosh
y41
xO
yR3
cosh
lS2
xE
pow
t71
y41
cosh);pow
cL
n61-1)));pow
lS2
xO
y93
xO.n71
0,sinh);xO
y41
pow);break
t11
cTan:{xE
sin,cos;sin
yX
cSin);sin
cL
xO
yR3
sin
lS2
cos
yX
cCos);cos
y41
xO
yR3
cos
lS2
xE
pow
t71
y41
cos);pow
cL
n61-1)));pow
lS2
xO
y93
xO.n71
0,sin);xO
y41
pow);break
t11
cPow:{iW1
p0=xO
l9
0);iW1
p1=xO
l9
1);if(p1
nF==cAdd){n72
x23(p1
e6)l01
p1
e6;++a
nQ
pow
t71
cL
p0);pow
cL
p1
nG3
pow
lS2
x23[a
t63
pow);}
xO
y93
xO
nB2
x72;}
t82
y83
t82
xO.Rehash(!keep_powi);iM1,false);
#ifdef DEBUG_SUBSTITUTIONS
iN1<<eJ3<<", "
<<e73
opcode)<<"->"
<<e73
xO
nF)<<": "
tY3
xO)xS1
xO);
#endif
eT1
xO
i72
EatFunc(cI3
eJ3,OPCODE
cZ3
x13
funcno
nQ
xO=CodeTreeFuncOp
yB(cZ3
funcno);n72
t33=Pop(eJ3);xO
nB2
tC2;xO.y72);
#ifdef DEBUG_SUBSTITUTIONS
iN1<<eJ3<<", "
tY3
xO)xS1
xO);
#endif
iM1);eT1
xO
i72
AddConst(yJ1
nQ
xO=CodeTreeImmed(value);iM1);Push(xO
i72
AddVar
nI2
varno
nQ
xO=CodeTreeVar
yB(varno);iM1);Push(xO
i72
SwapLastTwoInStack(){eU1
1
t63
eU1
2]i72
Dup(){Fetch(iO1-1
i72
Fetch(cI3
which){Push(stack[which]);}
nK1
T>void
Push(T
tree){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<tY3
tree)xS1
tree);
#endif
eT1
tree
i72
PopNMov(cI3
target,cI3
source){stack[target]=stack[source];stack
xA3
target+1);}
xE
yA2{clones.clear();xE
cQ3
stack.back());stack
xA3
iO1-1)eX
cX3;}
n72
Pop(cI3
n_pop){n72
cQ3
n_pop);for
nI2
n=0;n<n_pop;++n
i82[n
t63
eU1
n_pop+n]);
#ifdef DEBUG_SUBSTITUTIONS
for
n03=n_pop;n-->0;){iN1
eL3
cX3[n])xS1
cX3[n]);}
#endif
stack
xA3
iO1-n_pop)eX
cX3;}
cI3
GetStackTop(xC2
iO1;}
private:void
FindClone(xE&,bool=true)yU;}
private:n72
stack;std::multimap<fphash_t,xE>clones;bool
keep_powi;private:CodeTreeParserData
cQ1
CodeTreeParserData&);CodeTreeParserData&eD1=cQ1
CodeTreeParserData&);}
;eY
cN2
IfInfo{xE
eH2;xE
thenbranch;cI3
endif_location;IfInfo():eH2(),thenbranch(),endif_location(){}
}
;}
tJ{lC
GenerateFrom
cQ1
typename
FunctionParserBase
yB::Data&nV3,bool
keep_powi){n72
x52;x52.xF3
nV3.mVariablesAmount);for
nI2
n=0;n<nV3.mVariablesAmount;++n){x52.push_back(CodeTreeVar
yB(n+iJ2));}
GenerateFrom(nV3,x52,keep_powi);}
lC
GenerateFrom
cQ1
typename
FunctionParserBase
yB::Data&nV3,const
x6&x52,bool
keep_powi){xT1
x13>&eE1=nV3.mByteCode;const
std::vector
yB&Immed=nV3.mImmed;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"ENTERS GenerateFrom()\n"
;
#endif
CodeTreeParserData
yB
sim(keep_powi);std::vector<IfInfo
yB>eW;for(cI3
IP=0,DP=0;;++IP){tF2:while(!eW
i43&&(eW.eL==IP||(IP<n01&&n73
cJump&&eW.e91
lD2))){CodeTree
elsebranch=sim.yA2
lU2
eW.back().eH2)lU2
eW.e91)lU2
elsebranch)xL
3,cIf);eW.pop_back();}
if(IP>=n01)break;x13
opcode=eE1[IP];if((opcode==cSqr||opcode==cDup||(opcode==cInv&&!IsIntType
yB::cX3)||opcode==cNeg||opcode==cSqrt||opcode==cRSqrt||opcode==cFetch)){cI3
was_ip=IP;eB1
cZ2
ParsePowiSequence
yB(eE1,IP,eW
i43?n01:eW.eL,sim.xH
1);if
e83!=eB1(1.0)){xV
exponent)xY1;goto
tF2;}
if(opcode==cDup||opcode==cFetch||opcode==cNeg){eB1
xI2=ParseMuliSequence
yB(eE1,IP,eW
i43?n01:eW.eL,sim.xH
1);if(xI2!=eB1(1.0)){xV
xI2)xL
2,cMul);goto
tF2;}
}
IP=was_ip;}
if(lQ2>=iJ2){lC1=opcode-iJ2
lU2
x52[index]);}
else{switch(lQ2
nB3
cIf:case
cAbsIf:{eW
xA3
eW
eZ3+1);CodeTree
res(sim.yA2);eW.back().eH2.swap(res);eW.eL=n01;IP+=2;yC1}
case
cJump:{CodeTree
res(sim.yA2);eW.e91.swap(res);eW.eL=eE1[IP+1]+1;IP+=2;yC1}
case
cImmed:xV
Immed[DP++]);lD
cDup:sim.Dup();lD
cNop:lD
cFCall:{x13
funcno=nA3;assert(funcno<fpdata.mFuncPtrs.size());x13
t33=nV3.mFuncPtrs
nD3
mParams;sim.EatFunc(t33,lQ2,funcno);break
t11
cPCall:{x13
funcno=nA3;assert(funcno<fpdata.tV3.size());const
FunctionParserBase
yB&p=*nV3.tV3
nD3
mParserPtr;x13
t33=nV3.tV3
nD3
mParams;x6
paramlist=sim.Pop(tC2;CodeTree
tG2;tG2.GenerateFrom(*p.mData,paramlist)lU2
tG2);break
t11
cInv:xV
1
c82
cDiv);lD
cNeg
lV2
cNeg);break;xV
0
c82
cSub);lD
cSqr:xV
2)xY1;lD
cSqrt:xV
eB1(c92
cRSqrt:xV
eB1(-c92
cCbrt:xV
eB1(1)/eB1(3))xY1;lD
cDeg:xV
fp_const_rad_to_deg
y61
cRad:xV
fp_const_deg_to_rad
y61
cExp:iI)goto
y43;xV
fp_const_e
yB()c82
cPow);lD
cExp2:iI)goto
y43;xV
2.0
c82
cPow);lD
cCot
lV2
cTan);iI)x2
cCsc
lV2
cSin);iI)x2
cSec
lV2
cCos);iI)x2
cInt:
#ifndef __x86_64
iI)n83
1,cInt);t82
#endif
xV
eB1(0.5))tH2
xL
1,cFloor);lD
cLog10
lV2
yU3
fp_const_log10inv
y61
cLog2
lV2
yU3
fp_const_log2inv
y61
cD3:n93
yU3
fp_const_log2inv
yB())xL
3,cMul);lD
cHypot:xV
2)xY1;t83
xV
2)xY1
tH2;xV
eB1(c92
cSinCos:sim.Dup()xL
1,cSin);n93
cCos);lD
cSinhCosh:sim.Dup()xL
1,cSinh);n93
cCosh);lD
cRSub:t83
case
cSub:iI)n83
2,cSub);t82
xV-1)xL
2,cMul)tH2;lD
cRDiv:t83
case
cDiv:iI||IsIntType
yB::cX3)n83
2,cDiv);t82
xV-1)xY1
xL
2,cMul);lD
cAdd:case
cMul:case
cMod:case
cPow:case
cEqual:case
cLess:case
cGreater:case
tU1:case
cLessOrEq:case
cGreaterOrEq:case
cAnd:case
cOr:case
cAbsAnd:case
cAbsOr:sim.Eat(2,eV1
lD
cNot:case
cNotNot:case
nR2:case
cAbsNotNot
lV2
eV1
lD
cFetch:sim.Fetch(nA3);lD
cPopNMov:{x13
stackOffs_target=nA3;x13
stackOffs_source=nA3;sim.PopNMov(stackOffs_target,stackOffs_source);t82
y83
y43:;x13
funcno=opcode-cAbs;assert(funcno<FUNC_AMOUNT);const
FuncDefinition&func=Functions[funcno]xL
func.t33,eV1
t82}
}
Become(sim.yA2);
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Produced tree:\n"
;iS
#endif
}
}
#endif
#include <algorithm>
#ifdef FP_SUPPORT_OPTIMIZER
#include <assert.h>
#define FP_MUL_COMBINE_EXPONENTS
l13{using
l13
FUNCTIONPARSERTYPES;using
tJ;eY
static
void
AdoptChildrenWithSameOpcode(cS2{
#ifdef DEBUG_SUBSTITUTIONS
bool
nC2
e63
#endif
for
iU
if(y1
a)nF==i63){
#ifdef DEBUG_SUBSTITUTIONS
if(!nC2)y12"Before assimilation: "
yN
nC2=true;}
#endif
tree.AddParamsMove(y1
a).GetUniqueRef().tF1,a);}
#ifdef DEBUG_SUBSTITUTIONS
if(nC2)y12"After assimilation:   "
yN}
#endif
}
}
tJ{xC3
ConstantFolding(cS2{tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
void*yV3=0
xR1"["
<<(&yV3)<<"]Runs ConstantFolding for: "
yN
DumpHashes(tree)xR1
std::flush;
#endif
if(false){redo:;tree.Sort();
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&yV3)<<"]Re-runs ConstantFolding: "
yN
DumpHashes(tree);
#endif
}
if(i63!=cImmed){range
yB
p=iQ
tree);if(p
xU1
p
e3&&p
iA2==p
xF){xB
p
iA2);nG}
if(false){ReplaceTreeWithOne:xB
eB1(1));goto
do_return;ReplaceTreeWithZero:xB
x41;goto
do_return;ReplaceTreeWithParam0:
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Before replace: "
xR1
std::hex<<'['
t01
hash1<<','
t01
hash2<<']'<<std::dec
yN
#endif
tree.e12
0));
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"After replace: "
xR1
std::hex<<'['
t01
hash1<<','
t01
hash2<<']'<<std::dec
yN
#endif
goto
redo;iQ2
i63
nB3
cImmed:lD
iJ2:lD
cAnd:case
cAbsAnd
cO
bool
c9
e63
for
iU{if(!yS3
a)))c9=true;yX3
a),x31
cAbsAnd)nB3
IsNever
cN
x03:iJ);lD
lW1
iQ2
i93
nB3
0:iB
1:e93
x31
cAnd
eU3
xY3
y83
if(x31
cAnd||!c9)if(ConstantFolding_AndLogic
tC
xZ1
cOr:case
cAbsOr
cO
bool
c9
e63
for
iU{if(!yS3
a)))c9=true;yX3
a),x31
cAbsOr))lC2
iB
iY2
iJ);lD
lW1
iQ2
i93
nB3
0
cN
1:e93
x31
cOr
eU3
xY3
y83
if(x31
cOr||!c9)if(ConstantFolding_OrLogic
tC
xZ1
cNot:case
nR2:{x13
n21
0;switch(y1
0)nF
nB3
cEqual:n21
tU1;lD
tU1:n21
cEqual;lD
cLess:n21
cGreaterOrEq;lD
cGreater:n21
cLessOrEq;lD
cLessOrEq:n21
cGreater;lD
cGreaterOrEq:n21
cLess;lD
cNotNot:n21
cNot;lD
cNot:n21
cNotNot;lD
nR2:n21
cAbsNotNot;lD
cAbsNotNot:n21
nR2;break;y83
t82
if(opposite){e93
OPCODE(opposite)yH
SetParamsMove(y1
0).GetUniqueRef().tF1)xY3
iQ2
lT1
0),x31
nR2)nB3
x03
cN
iY2
iB
lW1
if(x31
cNot&&GetPositivityInfo(y1
0))==x03)e93
nR2);if(xZ3
cIf||xZ3
cAbsIf
nQ
nW3=y1
0);iW1
ifp1=nW3
l9
1);iW1
ifp2=nW3
l9
2);if(ifp1
nF==cNot||ifp1
nF==nR2){iC3
ifp1
nF==cNot
eU3;tI2
yR3
nE3);yW3
i91
yY3
eS2
if(ifp2
nF==cNot||ifp2
nF==nR2){iC3
i63);tI2);nE3);yW3
yX
ifp2
nF==cNot
eU3;yY3
l9
0)eS2
xZ1
cNotNot:case
cAbsNotNot:{if(yS3
0)))e7
yX3
0),x31
cAbsNotNot)nB3
IsNever
cN
x03:iB
lW1
if(x31
cNotNot&&GetPositivityInfo(y1
0))==x03)e93
cAbsNotNot);if(xZ3
cIf||xZ3
cAbsIf
nQ
nW3=y1
0);iW1
ifp1=nW3
l9
1);iW1
ifp2=nW3
l9
2);if(ifp1
nF==cNot||ifp1
nF==nR2){tree.SetParam(0,nW3
l9
0)yH
AddParam(ifp1);yW3
i91
yY3
eS2
if(ifp2
nF==cNot||ifp2
nF==nR2){iC3
i63);tI2);nE3
yH
AddParam(ifp2
yH
yZ
xZ1
cIf:case
cAbsIf:{if(ConstantFolding_IfOperations
tC
break
t11
cMul:{NowWeAreMulGroup:;AdoptChildrenWithSameOpcode(tree);eB1
nE1=eB1(1);cI3
iP1=0;bool
nF1=false
nD2
if(!eW1
yC1
eB1
immed=eX1;if(immed==x41
goto
ReplaceTreeWithZero;nE1*=immed;++iP1;}
if(iP1>1||(iP1==1&&lK3
nE1,xW3)nF1=true;if(nF1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cMul: Will add new "
tW3
nE1<<"\n"
;
#endif
for
iU
if(eW1{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
tW3
eX1
xR1"\n"
;
#endif
yZ3!lK3
nE1,xW3
tree
cL
cX1
nE1));iQ2
i93
nB3
0:iB
1:e7
y83
if(ConstantFolding_MulGrouping
tC
if(ConstantFolding_MulLogicItems
tC
xZ1
cAdd
cO
eB1
lW2=0.0;cI3
iP1=0;bool
nF1=false
nD2
if(!eW1
yC1
eB1
immed=eX1;lW2+=immed;++iP1;}
if(iP1>1||(iP1==1&&lW2==x41)nF1=true;if(nF1){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"cAdd: Will add new "
tW3
lW2<<"\n"
xR1"In: "
yN
#endif
for
iU
if(eW1{
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<" - For that, deleting "
tW3
eX1
xR1"\n"
;
#endif
yZ3!(lW2==eB1(0.0)))tree
cL
cX1
lW2));iQ2
i93
nB3
0
cN
1:e7
y83
if(ConstantFolding_AddGrouping
tC
if(ConstantFolding_AddLogicItems
tC
xZ1
cMin
cO
cI3
yB2=0
yG2
eE
nD2
while(a+1<i93&&y1
a)xI
y1
a+1)))iJ+1);yI2
max
iX2(!eE
e3||(p
xF)<eE
xF)){eE
xF=p
xF;eE
e3=true;yB2=a;}
}
if(eE
e3)for
iU{yI2
min
iX2
a!=yB2&&p
iA2>=eE
xF)yZ3
i93==1){e7
xZ1
cMax
cO
cI3
yB2=0
yG2
tD
nD2
while(a+1<i93&&y1
a)xI
y1
a+1)))iJ+1);yI2
min
iX2(!tD
cU2||p
iA2>tD
iA2)){tD
iA2=p
iA2;tD
cU2=true;yB2=a;}
}
if(tD
cU2){for
iU{yI2
max
iX2
a!=yB2&&(p
xF)<tD
iA2){iJ);}
}
}
if(i93==1){e7
xZ1
cEqual:case
tU1:case
cLess:case
cGreater:case
cLessOrEq:case
cGreaterOrEq:if(ConstantFolding_Comparison
tC
lD
cAbs:{range
yB
p0
e2
0));if
iV
e7
if(p0
e4{e93
cMul
yH
AddParam(n61
1)));goto
NowWeAreMulGroup;}
if(xZ3
cMul){iW1
p=y1
0);n72
nF3;n72
cA2
l01
p
e6;++a){p0=iQ
p
nG3
if
iV{nF3.push_back(p
nG3}
if(p0
e4{cA2.push_back(p
nG3}
}
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"Abs: mul group has "
<<nF3
eZ3<<" pos, "
<<cA2
eZ3<<"neg\n"
;
#endif
if(!nF3
i43||!cA2
i43){
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-Before: "
eL3
tree)xR1"\n"
<<std::flush;DumpHashes(tree
cB2;
#endif
xE
cO3;cO3
yX
cMul)l01
p
e6;++a){p0=iQ
p
nG3
if(iV||(p0
e4){}
else
cO3
cL
p
nG3}
cO3
lS2
xE
nH3;nH3
yX
cAbs);nH3
y41
cO3);nH3
lS2
xE
y31
cMul);x23
y41
nH3);y51
AddParamsMove(nF3);if(!cA2
i43){if(cA2
eZ3%2)x23
cL
n61-1)));y51
AddParamsMove(cA2);}
tree.Become(x72;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"AbsReplace-After: "
eL3
tree
cB2
xR1"\n"
<<std::flush;DumpHashes(tree
cB2;
#endif
goto
NowWeAreMulGroup;}
}
t82
#define HANDLE_UNARY_CONST_FUNC(funcname) nU){xB funcname(lS));nG
case
cLog:tK3(fp_log);if(xZ3
cPow
nQ
pow=y1
0);if(GetPositivityInfo(pow
l9
0))==x03){pow.iY
pow
l31
yH
lV
if(GetEvennessInfo(pow
l9
1))==x03){pow.iY
xE
abs;abs
yX
cAbs);abs
y41
pow
yR3
abs
lS2
pow
l31);pow.n71
0,abs
yH
lV}
l72
xZ3
cAbs
nQ
pow=y1
0)l9
0);if(pow
c93){pow.iY
xE
abs;abs
yX
cAbs);abs
y41
pow
yR3
abs
lS2
pow
l31);pow.n71
0,abs
yH
lV}
lD
cAcosh:tK3(fp_acosh);lD
cAsinh:tK3(fp_asinh);lD
cAtanh:tK3(fp_atanh);lD
cAcos:tK3(fp_acos);lD
cAsin:tK3(fp_asin);lD
cAtan:tK3(fp_atan);lD
cCosh:tK3(fp_cosh);lD
cSinh:tK3(fp_sinh);lD
cTanh:tK3(fp_tanh);lD
cSin:tK3(fp_sin);lD
cCos:tK3(fp_cos);lD
cTan:tK3(fp_tan);lD
cCeil:iA3(fp_ceil);lD
cTrunc:iA3(fp_trunc);lD
cFloor:iA3(fp_floor);lD
cInt:iA3(fp_int);lD
cCbrt:tK3(fp_cbrt);lD
cSqrt:tK3(fp_sqrt);lD
cExp:tK3(fp_exp);lD
cLog2:tK3(fp_log2);lD
cLog10:tK3(fp_log10);lD
cD3
i53
fp_log2(lS)*nI3
cArg:tK3(fp_arg);lD
cConj:tK3(fp_conj);lD
cImag:tK3(fp_imag);lD
cReal:tK3(fp_real);lD
cPolar
i53
fp_polar(nS2
cMod
i53
c33
nS2
cAtan2:{range
yB
p0
e2
0))yG2
p1
e2
1));nU&&lK3
lS,x41){if(p1
e3&&(p1
xF)<x41{xB
fp_const_pi
yB());nG
if(p1
xU1
p1
iA2>=yC2
x41;nG}
if(y1
l82
lK3
e0,x41){if(p0
e3&&(p0
xF)<x41{xB-fp_const_pihalf
yB());nG
if(p0
xU1
p0
iA2>x41{xB
fp_const_pihalf
yB());nG}
if
lJ
fp_atan2(lS,e0));nG
if((p1
xU1
p1
iA2>x41||(p1
e3&&(p1
xF)<fp_const_negativezero
yB())nQ
yD2;yD2
yX
cPow);yD2
y41
y1
1));yD2
cL
n61-1)));yD2
lS2
xE
yE2;yE2
y93
yE2
y41
y1
0));yE2
y41
yD2);yE2
cM
cAtan
yH
n71
0,yE2
iE1
1);xZ1
cPow:{if(ConstantFolding_PowOperations
tC
break
t11
cDiv:nU&&y1
l82
e0!=yC2
lS/nI3
cInv:nU&&lS!=yC2
eB1(1)/lS);nG
lD
cSub
i53
lS-nI3
cNeg:nU){xB-lS);nG
lD
cRad:nU){xB
RadiansToDegrees(i13
cDeg:nU){xB
DegreesToRadians(i13
cSqr:nU){xB
lS*lS);nG
lD
cExp2:tK3(fp_exp2);lD
cRSqrt:nU){xB
eB1(1)/fp_sqrt(i13
cCot:iB2
fp_tan
n0
cSec:iB2
fp_cos
n0
cCsc:iB2
fp_sin
n0
cHypot
i53
fp_hypot(nS2
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
cFCall:t82
do_return:;
#ifdef DEBUG_SUBSTITUTIONS
std::cout<<"["
<<(&yV3)<<"]Done ConstantFolding, result: "
yN
DumpHashes(tree);
#endif
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
tJ{xC3
range
yB::set_abs(nN
bool
has_negative=!min.known||min.val<eB1();bool
has_positive=!nJ3||t93>eB1();bool
crosses_axis=has_negative&&has_positive;rangehalf
yB
newmax;if(min
iX2
nJ3)newmax.set(fp_max(tV1,tW1);if(crosses_axis)min.set(eB1());c03
min
iX2
nJ3)min.set(fp_min(tV1,tW1);l72
min.known)min.set(tV1);else
min.set(tW1;}
max=newmax;}
xC3
range
yB::set_neg(){std::swap(min,max);min.val=-min.val;t93=-t93;}
eY
bool
IsLogicalTrueValue
cQ1
range
yB&p,bool
abs){if(nE
IsIntType
yB::cX3){if(p
xU1
p
iA2>=eB1(1
tA1;if(!abs&&p
e3&&p
xF<=eB1(-1
tA1;}
c03
p
xU1
p
iA2>=eB1(0.5
tA1;if(!abs&&p
e3&&p
xF<=eB1(-0.5
tA1;}
return
eC3
eY
bool
IsLogicalFalseValue
cQ1
range
yB&p,bool
abs){if(nE
IsIntType
yB::cX3){if(abs)return
p
e3
yS1
1);else
return
p
xU1
p
e3&&p
iA2>iB3
yS1
1);}
c03
abs)return
p
e3
yS1
0.5);else
return
p
xU1
p
e3&&p
iA2>eB1(-0.5)yS1
0.5);}
}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l13
FUNCTIONPARSERTYPES;using
tJ;l13{nK1
T>inline
int
isnan_workaround(T
t)yU(t!=t);}
}
tJ{eY
range
yB
iQ
xD3)
#ifdef DEBUG_SUBSTITUTIONS_extra_verbose
{using
l13
FUNCTIONPARSERTYPES
yG2
tmp=CalculateResultBoundaries_do(tree)xR1"Estimated boundaries: "
;if(tmp
cU2)std::cout<<tmp
iA2;else
std::cout<<"-inf"
xR1" .. "
;if(tmp
e3)std::cout<<tmp
xF;else
std::cout<<"+inf"
xR1": "
eL3
tree)xR1
std::endl
eX
tmp;}
eY
range
yB
CalculateResultBoundaries_do
cQ1
cS2
#endif
{iK
yG1(-fp_const_pihalf
yB(),fp_const_pihalf
yB());iK
pi_limits(-fp_const_pi
yB(),fp_const_pi
yB());iK
abs_pi_limits(x71,fp_const_pi
yB());iK
plusminus1_limits(eB1(-cT3
using
l13
std;switch(i63
nB3
cImmed:nR
tree
tJ1,tree
tJ1);case
cAnd:case
cAbsAnd:case
cOr:case
cAbsOr:case
cNot:case
nR2:case
cNotNot:case
cAbsNotNot:case
cEqual:case
tU1:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:{nR
x71,eB1(1))t11
cAbs:lE
m.set_abs(eU2
cLog:lE
c13,fp_log
tA3
fp_log
eU2
cLog2:lE
c13,fp_log2
tA3
fp_log2
eU2
cLog10:lE
c13,fp_log10
tA3
fp_log10
eU2
cAcosh:lE
yG
y01
cGreaterOrEq
eI2
fp_acosh
tB3
cGreaterOrEq
eI2
fp_acosh
eU2
cAsinh:lE
yG
set(fp_asinh
eR3
set(fp_asinh
eU2
cAtanh:lE
yG
n3-1),fp_atanh
tB3
cLess
eI2
fp_atanh
eU2
cAcos:lE
nR(m
e3&&(m
xF)<eB1(1))?fp_acos(m
xF):x71,(c23&&(t41)>=iB3)?fp_acos(t41):fp_const_pi
yB())t11
cAsin:lE
yG
n3-1),fp_asin,yG1
iA2
tB3
cLess
eI2
fp_asin,yG1
xF
eU2
cAtan:lE
yG
set(fp_atan,yG1
iA2
eR3
set(fp_atan,yG1
xF
eU2
cAtan2:{nU&&lK3
lS,x41)yU
abs_pi_limits;}
if(y1
l82
lK3
e0,x41)yU
yG1;}
return
pi_limits
t11
cSin:lE
bool
nS1=!c23||!m
e3||(m
xF-t41)>=(yD
nS1)t1
eB1
min=c33
t41,yD
min<x41
min
yS
eB1
max=c33
m
xF,yD
max<x41
max
yS
if(max<min)max
yS
bool
y11=(min<=fp_const_pihalf
yB()&&max>=fp_const_pihalf
yB());bool
nG1=(min<=cG&&max>=cG);if(y11&&nG1)t1
if(nG1)nR
iB3,yF2;if(y11)nR
yK2
eB1(1));nR
yK2
yF2
t11
cCos:lE
if(c23)t41+=fp_const_pihalf
yB();if(m
e3)m
xF+=fp_const_pihalf
yB();bool
nS1=!c23||!m
e3||(m
xF-t41)>=(yD
nS1)t1
eB1
min=c33
t41,yD
min<x41
min
yS
eB1
max=c33
m
xF,yD
max<x41
max
yS
if(max<min)max
yS
bool
y11=(min<=fp_const_pihalf
yB()&&max>=fp_const_pihalf
yB());bool
nG1=(min<=cG&&max>=cG);if(y11&&nG1)t1
if(nG1)nR
iB3,yF2;if(y11)nR
yK2
eB1(1));nR
yK2
yF2
t11
cTan:{nR)t11
cCeil:lE
m
eV
cFloor:lE
t91
eX
m
t11
cTrunc:lE
t91;m
eV
cInt:lE
t91;m
eV
cSinh:lE
yG
set(fp_sinh
eR3
set(fp_sinh
eU2
cTanh:lE
yG
set(fp_tanh,plusminus1_limits.min
eR3
set(fp_tanh,plusminus1_limits.max
eU2
cCosh:lE
if(c23){if(m
e3){if(t41>=x71&&m
xF>=x41{t41=e8}
l72(t41)<x71&&m
xF>=x41{eB1
tmp=e8
if(tmp>m
xF)m
xF=tmp;t41=eB1(1);}
else{t41=e8
std::swap(t41,m
xF);}
}
c03
t41>=x41{m.cH
t41=fp_cosh(t41);}
else{m.cH
t41=eB1(1);}
}
}
else{c23=true;t41=eB1(1);if(m
e3){t41=fp_cosh(m
xF);m.cH}
else
m.cH}
return
m
t11
cIf:case
cAbsIf:{range
yB
res1
e2
1))yG2
res2
e2
2));if(!res2
cU2)res1
cU2
e63
l72
res1
xU1(res2
iA2)<res1
iA2)res1
iA2=res2
iA2;l72
isnan_workaround(res2
iA2))res1
iA2=res2
iA2;if(!res2
e3)res1.cH
l72
res1
e3&&(res2
xF)>res1
xF)res1
xF=res2
xF;l72
isnan_workaround(res2
xF))res1
xF=res2
xF
eX
res1
t11
cMin:{bool
iL
e63
bool
iM=false
yG2
cX3;x7
m
e2
c43
c23)iL=true;yT1
cU2||(t41)<cX3
iA2
i92=t41;if(!m
e3)iM=true;yT1
e3||(m
xF)<cV3
i82
xF=m
xF;}
if(iL
i82
cU2
e63
if(iM
i82.cH
return
cX3
t11
cMax:{bool
iL
e63
bool
iM=false
yG2
cX3;x7
m
e2
c43
c23)iL=true;yT1
cU2||t41>cX3
iA2
i92=t41;if(!m
e3)iM=true;yT1
e3||m
xF>cV3
i82
xF=m
xF;}
if(iL
i82
cU2
e63
if(iM
i82.cH
return
cX3
t11
cAdd:{range
yB
cQ3
x71,x41;x7
item
e2
a));if(item
cU2
i92+=item
iA2;else
cS3
e63
if(item
e3
i82
xF+=item
xF;else
cX3.cH
if(!cX3
xU1!cU3)t82
if(cX3
xU1
cU3&&cX3
iA2>cV3)std::swap
yH2,cV3)eX
cX3
t11
cMul:{cN2
Value{enum
nK3{tJ2,iQ1,tK2}
;nK3
eM;eB1
value;Value(nK3
t):eM(t),value(0){}
Value(eB1
v):eM(tJ2),value(v){}
bool
cD2
xC2
eM==iQ1||(eM==tJ2&&value<x41
iD2
eD1*=cQ1
Value&rhs){if(eM==tJ2&&rhs.eM==tJ2)value*=rhs.value;else
eM=(cD2)!=rhs.cD2)?iQ1:tK2);}
tW2<cQ1
Value&xB2(eM==iQ1&&rhs.eM!=iQ1)||(eM==tJ2&&(rhs.eM==tK2||(rhs.eM==tJ2&&value<rhs.value)));}
}
;cN2
yH1{Value
yL2,yM2;yH1():yL2(Value::tK2),yM2(Value::iQ1){}
void
x62
Value
c53,const
Value&value2){c53*=value2;if(c53<yL2)yL2=c53;if(yM2<c53)yM2=c53;}
}
yG2
cQ3
eB1(cT3
x7
item
e2
c43
item
xU1!item
e3)nR);Value
nL3=cS3?Value
yH2):eY1
iQ1);Value
nM3=cU3?Value(cV3):eY1
tK2);Value
nN3=item
cU2?Value(item
iA2):eY1
iQ1);Value
nO3=item
e3?Value(item
xF):eY1
tK2);yH1
range;range.x62
nL3,nN3);range.x62
nL3,nO3);range.x62
nM3,nN3);range.x62
nM3,nO3);if(range.yL2.eM==Value::tJ2
i92=range.yL2.value;else
cS3
e63
if(range.yM2.eM==Value::tJ2
i82
xF=range.yM2.value;else
cX3.cH
if(!cX3
xU1!cU3)t82
if(cX3
xU1
cU3&&cX3
iA2>cV3)std::swap
yH2,cV3)eX
cX3
t11
cMod:{range
yB
x
e2
0))yG2
y
e2
1));if(y
e3){if(y
xF>=x41{if(!x
cU2||(x
iA2)<x41
nR-y
xF,y
xF);else
nR
x71,y
xF);}
c03!x
e3||(x
xF)>=x41
nR
y
xF,-y
xF);else
nR
y
xF,fp_const_negativezero
yB());}
}
else
nR)t11
cPow:{if(y1
l82
e0==x41{nR
eB1(cT3}
nU&&lS==x41{nR
x71,x41;}
nU&&lK3
lS,xW3{nR
eB1(cT3}
if(y1
l82
e0>x71&&GetEvennessInfo(y1
1))==x03){eB1
cZ2
e0
yG2
tmp
e2
0))yG2
cX3;cS3=true;cX3
iA2=0;if(tmp
xU1
tmp
iA2>=x41
cX3
iA2=eD3
tmp
iA2,tL2
l72
tmp
e3&&tmp
xF<=x41
cX3
iA2=eD3
tmp
xF,tL2
cX3.cH
if(tmp
xU1
tmp
e3){cU3=true;cV3=fp_max(fp_abs(tmp
iA2),fp_abs(tmp
xF));cV3=eD3
cV3,tL2}
return
cX3;}
range
yB
p0
e2
0))yG2
p1
e2
1));TriTruthValue
p0_positivity=(p0
xU1
cF2)>=x41?x03:(p0
e3&&(p0
xF)<x71?iY2
Unknown);TriTruthValue
cE2=GetEvennessInfo(y1
1));TriTruthValue
tE=Unknown;switch(p0_positivity)lC2
tE=x03;lD
iY2{tE=cE2;t82
y83
switch(cE2)lC2
tE=x03;lD
iY2
lD
Unknown:{if(y1
l82!t32
e0)&&e0>=x41{tE=x03;}
t82}
iQ2
tE)lC2{eB1
c63
if(p0
xU1
p1
cU2){min=fp_pow
cF2,p1
iA2);if
cF2<x71&&(!p1
e3||p1
xF>=x41&&min>=x41
c63
if
cF2==x71||p1
iA2<x41
c63}
if(p0
xU1
p0
iA2>=x71&&p0
e3&&p1
e3){eB1
max=eD3
p0
xF,p1
xF);if(min>max)std::swap(min,max);nR
min,max);}
nR
min,false)t11
iY2{nR
false,fp_const_negativezero
yB());}
y83{t82
xZ1
cNeg:lE
m.set_neg(eU2
cSub:{nA
cNeg
cG2
1))tC3
cAdd
cH2
t42
lF
cInv:{c41(-1)))lF
cDiv:{nA
cInv
cG2
1))tC3
cMul
cH2
t42
lF
cRad:{nC
cMul)n7
fp_const_rad_to_deg
yB()))lF
cDeg:{nC
cMul)n7
fp_const_deg_to_rad
yB()))lF
cSqr:{c41(2)))lF
cExp:{nC
cPow);tmp
cL
CodeTreeImmed(fp_const_e
yB())cH2
lF
cExp2:{nC
cPow);tmp
cL
tD3
cH2
lF
cCbrt:lE
yG
set(fp_cbrt
eR3
set(fp_cbrt
eU2
cSqrt:lE
if(c23)t41=(t41)<x71?0:fp_sqrt(t41);if(m
e3)m
xF=(m
xF)<x71?0:fp_sqrt(m
xF
eU2
cRSqrt:{c41(-0.5)))lF
cHypot:{xE
xsqr,ysqr,add,sqrt;xsqr.x0
0));xsqr
cL
tD3);ysqr.x0
1));ysqr
cL
tD3);xsqr
yX
cPow);ysqr
yX
cPow);add
y41
xsqr);add
y41
ysqr);add
yX
cAdd);sqrt
y41
add);sqrt
yX
cSqrt)eX
iQ
sqrt)t11
cD3:{nA
cLog2
cG2
0))tC3
cMul)t42;tmp.x0
1))lF
cCot:{nA
cTan
cG2
xY
lF
cSec:{nA
cCos
cG2
xY
lF
cCsc:{nA
cSin
cG2
xY
eX
iQ
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
iJ2:lD
cArg:case
cConj:case
cImag:case
cReal:case
cPolar:lD
cPCall:lD
cFCall:t82
nR);}
eY
TriTruthValue
GetIntegerInfo
cQ1
cS2{switch(i63
nB3
cImmed:return
t32
tree
tJ1)?x03:IsNever;case
cFloor:case
cCeil:case
cTrunc:case
cInt:return
x03;case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cEqual:case
tU1:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:return
x03;case
cIf:{TriTruthValue
a=GetIntegerInfo(y1
1));TriTruthValue
b=GetIntegerInfo(y1
2));if(a==b)return
a
eX
Unknown
t11
cAdd:case
cMul:{for
iU
if(GetIntegerInfo(y1
a))!=x03)return
Unknown
eX
x03;}
y83
t82
return
Unknown;}
eY
bool
IsLogicalValue
cQ1
cS2{switch(i63
nB3
cImmed:return
lK3
tree
tJ1,x41||lK3
tree
tJ1,eB1(1));case
cAnd:case
cOr:case
cNot:case
cNotNot:case
cAbsAnd:case
cAbsOr:case
nR2:case
cAbsNotNot:case
cEqual:case
tU1:case
cLess:case
cLessOrEq:case
cGreater:case
cGreaterOrEq:nZ
cMul:{for
iU
if(!yS3
a))cI
return
true
t11
cIf:case
cAbsIf:yU
yS3
1))&&yS3
2));}
y83
t82
return
eC3}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l13
FUNCTIONPARSERTYPES;
#if defined(__x86_64) || !defined(FP_SUPPORT_CPLUSPLUS11_MATH_FUNCS)
# define CBRT_IS_SLOW
#endif
#if defined(DEBUG_POWI) || defined(DEBUG_SUBSTITUTIONS)
#include <cstdio>
#endif
l13
x91{extern
const
x13
char
powi_table[256];}
l13{using
tJ;eY
bool
IsOptimizableUsingPowi(long
immed,long
penalty=0){x91
eW2
synth;nS3
PushVar(iJ2);cI3
bytecodesize_backup=nS3
GetByteCodeSize();x91::nR1
immed,x91::tM1
yB::MulSequence,synth);cI3
bytecode_grow_amount=nS3
GetByteCodeSize()-bytecodesize_backup
eX
bytecode_grow_amount<cI3(MAX_POWI_BYTECODE_LENGTH-penalty);}
xC3
ChangeIntoRootChain(xE&tree,bool
iZ2,long
tN2,long
tO2){while(tO2>0){nC
cCbrt);nZ1);tmp.Rehash(nT2--tO2;}
while(tN2>0){nC
cSqrt);if(iZ2){tmp
yX
cRSqrt);iZ2
e63}
nZ1);tmp.Rehash(nT2--tN2;}
if(iZ2){nC
cInv);nZ1
nT2}
}
eY
cN2
RootPowerTable{static
const
eB1
RootPowers[(1+4)*(1+3)];}
;eY
const
eB1
tL(1+4)*(1+3)]={eB1(1)lT
tE3
tE3
2*tE3
i31)lT
3*2)lT
3*2*2)lT
3*i31*2*i31*3)lT
3*3*eJ2
2*eJ2
i31*3*2*i31*3*3)lT
3*3*3*eJ2
3*2*eJ2
3*i31*3*3*2*2*2*2)}
;cN2
PowiResolver{static
const
x13
MaxSep=4;static
nQ3
MaxOp=5;typedef
int
c73;typedef
long
nX3;typedef
long
tM;cN2
yN2{yN2():n_int_sqrt(0),n_int_cbrt(0),sep_list(),lX1(0){}
int
n_int_sqrt;int
n_int_cbrt;int
sep_list[MaxSep];tM
lX1;}
;eY
static
yN2
CreatePowiResult(eB1
exponent){yN2
cX3;c73
t2=FindIntegerFactor(tL2
if(t2==0){
#ifdef DEBUG_POWI
tP2"no factor found for %Lg\n"
,(c83);
#endif
return
cX3;}
cX3.lX1=y21
exponent,t2);nX3
eK2=EvaluateFactorCost(t2,0,0,0)+cF
cX3.lX1);int
tF3=0;int
tG3=0;int
nP3=0;
#ifdef DEBUG_POWI
tP2"orig = %Lg\n"
,(c83);tP2"plain factor = "
tZ3"%ld\n"
,(int)t2,(long)eK2);
#endif
for
nI2
n_s=0;n_s<MaxSep;++n_s){int
xJ=0;nX3
yI1=eK2;c73
yU1=t2;for(int
s=1;s<MaxOp*4;++s){
#ifdef CBRT_IS_SLOW
if(s>=MaxOp)break;
#endif
int
n_sqrt=s%MaxOp;int
n_cbrt=s/MaxOp;if(n_sqrt+n_cbrt>4)yC1
eB1
lH1=exponent;lH1-=tL
s];tY1=FindIntegerFactor(lH1);if(xI2!=0){tM
xP=y21
lH1,xI2);nX3
cost=EvaluateFactorCost(xI2,tF3+n_sqrt,tG3+n_cbrt,nP3+1)+cF
xP);
#ifdef DEBUG_POWI
tP2"Candidate sep %u (%d*sqrt %d*cbrt)factor = "
tZ3"%ld (for %Lg to %ld)\n"
,s,n_sqrt,n_cbrt,xI2,(long)cost
tM2
lH1,(long)xP);
#endif
if(cost<yI1){xJ=s;yU1=xI2;yI1=cost;}
}
}
if(!xJ)break;
#ifdef DEBUG_POWI
tP2"CHOSEN sep %u (%d*sqrt %d*cbrt)factor = "
tZ3"%ld, exponent %Lg->%Lg\n"
,xJ,xJ%MaxOp,xJ/MaxOp,yU1,yI1
tM2
e83)tM2
e83-tL
xJ]));
#endif
cX3
c51
n_s]=xJ
cT1-=tL
xJ];tF3+=xJ%MaxOp;tG3+=xJ/MaxOp;eK2=yI1;t2=yU1;nP3+=1;}
cX3.lX1=y21
exponent,t2);
#ifdef DEBUG_POWI
tP2"resulting exponent is %ld (from exponent=%Lg, best_factor=%Lg)\n"
,cX3.lX1,(c83
tM2
t2);
#endif
while(t2%2==0){++cX3.n_int_sqrt;t2/=2;}
while(t2%3==0){++cX3.n_int_cbrt;t2/=3;}
return
cX3;}
private:static
nX3
cF
tM
xP){static
std::map
cJ2
iC;if(xP<0){nX3
cost=22
eX
cost+cF-xP);}
std::map
cJ2::xV3
i=iC.xJ2
xP);if(i!=iC.cL1
xP)return
i
cY2;std::pair
cJ2
cQ3
xP,0.0);nX3&cost=cX3.tJ3
while(xP>1){int
xI2=0;if(xP<256){xI2=x91::powi_table[xP];if(xI2&128)xI2&=127;else
xI2=0;if(xI2&64)xI2=-(xI2&63)-1;}
if(xI2){cost+=cF
xI2);xP/=xI2;yC1}
if(!(xP&1)){xP/=2;cost+=6;}
else{cost+=7;xP-=1;}
}
iC.xR3,cX3)eX
cost;}
eY
static
tM
y21
yJ1,tY1)yU
makeLongInteger(value*eB1(xI2));}
eY
static
bool
yK1
yJ1,tY1){eB1
v
xN3*eB1(xI2)eX
isLongInteger(v);}
eY
static
c73
FindIntegerFactor(yJ1){tY1=(2*2*2*2);
#ifdef CBRT_IS_SLOW
#else
xI2*=(3*3*3);
#endif
c73
cX3=0;if(yK1
value,xI2)){cX3=xI2;while((xI2%2)==0&&yK1
value,xI2/2)i82=xI2/=2;while((xI2%3)==0&&yK1
value,xI2/3)i82=xI2/=3;}
#ifdef CBRT_IS_SLOW
if(cX3==0){if(yK1
value,3
iE2
3;}
#endif
return
cX3;}
static
int
EvaluateFactorCost(int
xI2,int
s,int
c,int
nmuls){nQ3
nR3=6;
#ifdef CBRT_IS_SLOW
nQ3
eL2=25;
#else
nQ3
eL2=8;
#endif
int
cX3=s*nR3+c*eL2;while(xI2%2==0){xI2/=2;cX3+=nR3;}
while(xI2%3==0){xI2/=3;cX3+=eL2;}
cX3+=nmuls
eX
cX3;}
}
;}
tJ{eY
bool
xE::RecreateInversionsAndNegations(bool
prefer_base2){bool
changed=false
l01
yQ1;cP3
lS1.RecreateInversionsAndNegations(prefer_base2))nV2
if(changed){exit_changed:Mark_Incompletely_Hashed(lA2
switch(iJ1{case
cMul:{n72
lX2;xE
lY2,cG1;if(true){bool
nH1
e63
eB1
nU2=0;l11
nV
c93&&tT
0)cA3
tT
1)l92){nH1=true;nU2=tT
1)tJ1;t82}
if(nH1){eB1
immeds=1.0;l11
nV
l92){immeds*=powgroup
tJ1;yL1}
l11-->0;nQ&powgroup=lS1;if(powgroup
c93&&tT
0)cA3
tT
1)l92
nQ&log2=tT
0);log2.iY
log2
yX
cD3);log2
cL
CodeTreeImmed(eD3
immeds,eB1(1)/nU2)));log2
lS2
t82}
}
}
l11
nV
c93&&tT
1)l92){iW1
exp_param=tT
1);eB1
cZ2
exp_param
tJ1;if(cS1,iB3)){iY
lX2.push_back(lS1
yR3
yL1
l72
exponent<x71&&t32
exponent)nQ
iN;iN
yX
cPow);iN
cL
tT
0));iN
cL
CodeTreeImmed(-exponent));iN
lS2
lX2.push_back(iN);iY
yL1}
l72
powgroup
cA3!lY2
lD2){lY2=tT
0);iY
yL1
l72
powgroup
nF==cD3&&!cG1
lD2){cG1=powgroup;iY
yL1}
if(!lX2
i43){nV2
xE
iK1;iK1
y93
iK1
nB2
lX2);iK1
lS2
xE
y31
cMul);y51
SetParamsMove(tF
if(y51
IsImmed()&&lK3
y51
GetImmed(),xW3{eO2
cInv)tH3;}
c03
y51
GetDepth()>=iK1
n62){eO2
cDiv)t0
x72
tH3;}
else{eO2
cRDiv)tH3
t0
x72;}
}
}
if(lY2
lD2
nQ
y31
iJ1;y51
SetParamsMove(tF
while(y51
RecreateInversionsAndNegations(prefer_base2))y51
FixIncompleteHashes();eO2
cD3)t0
lY2)t0
x72;nV2}
if(cG1
lD2
nQ
y31
cMul);x23
y41
cG1
l9
1));y51
AddParamsMove(tF
while(y51
RecreateInversionsAndNegations(prefer_base2))y51
FixIncompleteHashes();DelParams();eO2
cD3)t0
cG1
l9
0))t0
x72;nV2
xZ1
cAdd:{n72
tQ2;l11-->0;)if(cE3
cMul){lZ2
y71:;xE&x23
nX2
a
eN
for(cI3
b=x23
e6;b-->0;){if(x23
l9
b).lN1
xI2=x23
l9
b)tJ1;lJ3
xI2
xW
y71;}
y51
iY
y51
DelParam(b);tG=!tG;}
l72
lK3
xI2,eB1(-2)))eO
y71;}
y51
iY
y51
DelParam(b);x23
cL
tD3)tR2}
if(tG){y51
tN
x72;yL1}
l72
cE3
cDiv&&!IsIntType
yB::cX3){lZ2
y81:;xE&iK1
nX2
a
eN
if(iK1
l9
0)l92){lJ3
iK1
l9
0)tJ1
xW
y81;}
iK1.iY
iK1.DelParam(0);iK1
yX
cInv)tR2
if(tG)eO
y81;}
iK1.tN
iK1);yL1}
l72
cE3
cRDiv&&!IsIntType
yB::cX3){lZ2
x01:;xE&iK1
nX2
a
eN
if(iK1
l9
1)l92){lJ3
iK1
l9
1)tJ1
xW
x01;}
iK1.iY
iK1.DelParam(1);iK1
yX
cInv)tR2
if(tG)eO
x01;}
iK1.tN
iK1);yL1}
if(!tQ2
i43){
#ifdef DEBUG_SUBSTITUTIONS
tP2"Will make a Sub conversion in:\n"
);fflush(stdout);iS
#endif
xE
yO2;yO2
yX
cAdd);yO2
nB2
tQ2);yO2
lS2
xE
cH1;cH1
yX
cAdd);cH1
nB2
tF1);cH1
lS2
if(cH1
l92&&lK3
cH1
tJ1,x41){eO2
cNeg);eP);}
c03
cH1
n62==1){eO2
cRSub);eP)t0
cH1);}
l72
yO2
nF==cAdd){eO2
cSub)t0
cH1);eP
yR3
y33
1;a<yO2
e6;++a
nQ
eM2;eM2
yX
cSub);eM2
nB2
tF1);eM2.y72)t0
eM2);eP
nG3}
}
else{eO2
cSub)t0
cH1);eP);}
}
#ifdef DEBUG_SUBSTITUTIONS
tP2"After Sub conversion:\n"
);fflush(stdout);iS
#endif
xZ1
cPow:{iW1
p0
nX2
0);iW1
p1
nX2
1);if(p1
l92){if
xM1!=x71&&!isInteger
xM1)){eS
yN2
r=eS
CreatePowiResult(fp_abs
xM1));if(r.lX1!=0){bool
iR1
e63
if
xM1<x71&&r
c51
0]==0&&r.n_int_sqrt>0){iR1=true;}
#ifdef DEBUG_POWI
tP2"Will resolve powi %Lg as powi(chain(%d,%d),%ld)"
tM2
fp_abs
xM1),r.n_int_sqrt,r.n_int_cbrt,r.lX1);for
nI2
n=0;n<eS
MaxSep;++n){if(r
c51
n]==0)break;int
n_sqrt=r
c51
n]%eS
MaxOp;int
n_cbrt=r
c51
n]/eS
MaxOp;tP2"*chain(%d,%d)"
,n_sqrt,n_cbrt);}
tP2"\n"
);
#endif
xE
cK2
nX2
0);xE
yP2=cK2;yP2.iY
ChangeIntoRootChain(yP2,iR1,r.n_int_sqrt,r.n_int_cbrt);yP2
lS2
xE
pow;if(r.lX1!=1){pow
yX
cPow);pow
y41
yP2);pow
cL
n61
r.lX1)));}
else
pow.swap(yP2);xE
mul;mul
y93
mul
y41
pow);for
nI2
n=0;n<eS
MaxSep;++n){if(r
c51
n]==0)break;int
n_sqrt=r
c51
n]%eS
MaxOp;int
n_cbrt=r
c51
n]/eS
MaxOp;xE
eN2=cK2;eN2.iY
ChangeIntoRootChain(eN2,false,n_sqrt,n_cbrt);eN2
lS2
mul
y41
eN2);}
if
xM1<x71&&!iR1){mul
lS2
eO2
cInv);n71
0,mul);DelParam(1);}
else{eO2
cMul);SetParamsMove(mul.tF1);}
#ifdef DEBUG_POWI
iS
#endif
nV2
t82}
}
if(GetOpcode()==cPow&&(!p1
l92||!isLongInteger
xM1)||!IsOptimizableUsingPowi
yB(makeLongInteger
xM1)))){if(p0
l92&&p0
tJ1>eB1(0.0)){if(prefer_base2){eB1
yQ2=fp_log2(p0
tJ1);lJ3
yQ2,xW3{DelParam(0);}
else{lX
cX1
yQ2))cT1
cL
p1)cT1.Rehash(iZ}
eO2
cExp2);nV2}
else{eB1
yQ2=fp_log(p0
tJ1);lJ3
yQ2,xW3{DelParam(0);}
else{lX
cX1
yQ2))cT1
cL
p1)cT1.Rehash(iZ}
eO2
cExp);nV2}
}
l72
GetPositivityInfo(p0)==x03){if(prefer_base2
nQ
log;log
yX
cLog2);log
cL
p0);log
lS2
lX
p1)cT1
y41
log)cT1
lS2
eO2
cExp2
iZ
nV2}
else{xE
log;log
yX
cLog);log
cL
p0);log
lS2
lX
p1)cT1
y41
log)cT1
lS2
eO2
cExp
iZ
nV2}
}
xZ1
cDiv:{if(GetParam(0)l92&&lK3
GetParam(0)tJ1,xW3{eO2
cInv);DelParam(0);}
t82
y83
t82
if(changed)goto
exit_changed
eX
changed;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
using
l13
FUNCTIONPARSERTYPES;l13{using
tJ;class
nF2{cI3
nI1;cI3
eT;cI3
eU;cI3
lI1;cI3
tH;cI3
tI;cI3
n31;cY3
nF2():nI1(0),eT(0),eU(0),lI1(0),tH(0),tI(0),n31(0){}
void
eP3
OPCODE
op){nI1+=1
nE2
cCos)++eT
nE2
cSin)++eU
nE2
cSec)++eT
nE2
cCsc)++eU
nE2
cTan)++lI1
nE2
cCot)++lI1
nE2
cSinh)++tI
nE2
cCosh)++tH
nE2
cTanh)++n31;}
cI3
GetCSEscore()const{cI3
cX3=nI1
eX
cX3;}
int
NeedsSinCos()const{bool
y91=(nI1==(eT+eU+lI1));if((lI1&&(eU||eT))||(eU&&eT)){if(y91)return
1
eX
2;}
return
0;}
int
NeedsSinhCosh()const{bool
y91=(nI1==(tH+tI+n31));if((n31&&(tI||tH))||(tI&&tH)){if(y91)return
1
eX
2;}
return
0;}
cI3
MinimumDepth()const{cI3
n_sincos=std::min(eT,eU);cI3
n_sinhcosh=std::min(tH,tI);if(n_sincos==0&&n_sinhcosh==0)return
2
eX
1;}
}
;x63
TreeCountType:public
std::multimap<fphash_t,std::pair<nF2,xE> >{}
;xC3
FindTreeCounts(eZ1&n02,xD3,OPCODE
x22,bool
skip_root=false){eA
i=n02.xJ2
xE3);if(!skip_root){bool
found
e63
for(;i!=n02.cL1
xE3;++i){if(tree
xI
i
cY2.second)){i
cY2
t03.eP3
x22);found=true;t82}
if(!found){nF2
count;count.eP3
x22);n02.xR3,std::make_pair(xE3,std::make_pair
e53,tree)));}
}
eT2
yI++a)FindTreeCounts(n02,y1
a),i63);}
cN2
c1{bool
BalanceGood;bool
FoundChild;}
;eY
c1
lJ1
iW1
root,iW1
cF3{if(root
xI
cF3){c1
cX3={true,true}
eX
cX3;}
c1
cX3={true,false}
;if(root
nF==cIf||root
nF==eB3{c1
cond=lJ1
root
l9
0),cF3;c1
y0=lJ1
root
l9
1),cF3;c1
y5=lJ1
root
l9
2),cF3;if(cond
c2||y0
c2||y5
c2){cX3
c2=true;}
cX3
eQ=((y0
c2==y5
c2)||(cond
c2
cL2&&(cond
eQ||(y0
c2&&y5
c2))&&(y0
eQ||(cond
c2
cL2&&(y5
eQ||(cond
c2
cL2;}
else{bool
tZ1
e63
bool
nJ1
e63
for(cI3
b=root
e6,a=0;a<b;++a){c1
tmp=lJ1
root
l9
a),cF3;if(tmp
c2
i82
c2=true;if(tmp
eQ==false)tZ1=true;l72
tmp
c2)nJ1=true;}
if(tZ1&&!nJ1
i82
eQ
e63}
return
cX3;}
eY
bool
lI3
lB3
tI3
xD3,const
x91
eW2&synth,const
eZ1&n02){for(cI3
b=i93,a=0;a<b;++a){iW1
leaf=y1
a);eA
synth_it;eX2
eZ1::const_iterator
i=n02.xU3
i!=n02.end();++i){if(i->first!=leaf.GetHash())yC1
const
nF2&occ=i
cY2
t03;cI3
score=occ.GetCSEscore();iW1
candidate=i
cY2.tJ3
if(l62
candidate))yC1
if(leaf
n62<occ.MinimumDepth())yC1
if(score<2)yC1
if(lJ1
tI3
leaf)eQ==false)continue
l43
if(lI3(tI3
leaf,synth,n02
tA1;}
return
eC3
eY
bool
n12
lB3
xS3,iW1
expr){yV1
xS3
lH3
expr
tA1;yV1
n12(xS3
l9
a),expr
tA1
eX
eC3
eY
bool
GoodMomentForCSE
lB3
xS3,iW1
expr){if(xS3
nF==cIf)return
true;yV1
xS3
lH3
expr
tA1;cI3
tS2=0;yV1
n12(xS3
l9
a),expr))++tS2
eX
tS2!=1;}
}
tJ{eY
cI3
xE::SynthCommonSubExpressions(x91::xV1
const{if(yQ1==0)return
0;cI3
stacktop_before=nS3
GetStackTop();eZ1
n02;FindTreeCounts(n02,*this,GetOpcode(),true);for(;;){cI3
yR2=0;
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Finding a CSE candidate, root is:"
<<std::endl;DumpHashes(*this);
#endif
eA
cs_it(n02.end());for(eA
j=n02.xU3
j!=n02.end();){eA
i(j++);const
nF2&occ=i
cY2
t03;cI3
score=occ.GetCSEscore();xD3=i
cY2.tJ3
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"Score "
<<score<<":\n"
<<std::flush;DumpTreeWithIndent(tree);
#endif
if(l62
tree))xT
if(tree
n62<occ.MinimumDepth())xT
if(score<2)xT
if(lJ1*this,tree)eQ==false)xT
if(lI3(*this,tree,synth,n02)){yC1}
if(!GoodMomentForCSE(*this,tree))xT
score*=tree
n62;if(score>yR2){yR2=score;cs_it=i;}
}
if(yR2<=0){
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<"No more CSE candidates.\n"
<<std::flush;
#endif
t82
xD3=cs_it
cY2.tJ3
#ifdef DEBUG_SUBSTITUTIONS_CSE
std::cout<<i33"Common Subexpression:"
;DumpTree
yB(tree)xR1
std::endl;
#endif
#if 0
int
lY1=occ.NeedsSinCos();int
i3=occ.NeedsSinhCosh();xE
tT2,tU2,yS2,yT2;if(lY1){tT2
y13
tT2
yX
cSin);tT2
lS2
tU2
y13
tU2
yX
cCos);tU2
lS2
if(l62
tT2)||l62
tU2)){if(lY1==2){tB1
yC1}
lY1=0;}
}
if(i3){yS2
y13
yS2
yX
cSinh);yS2
lS2
yT2
y13
yT2
yX
cCosh);yT2
lS2
if(l62
yS2)||l62
yT2)){if(i3==2){tB1
yC1}
i3=0;}
}
#endif
tree.SynthesizeByteCode(synth,false);tB1
#ifdef DEBUG_SUBSTITUTIONS_CSE
nS3
template
Dump<0>()xR1"Done with Common Subexpression:"
;DumpTree
yB(tree)xR1
std::endl;
#endif
#if 0
if(lY1){if(lY1==2||i3){nS3
eA1}
nZ2
cSinCos,1,2
tC1
tT2,1
tC1
tU2,0);}
if(i3){if(lY1)nS3
eA1
if(i3==2){nS3
eA1}
nZ2
cSinhCosh,1,2
tC1
yS2,1
tC1
yT2,0);}
#endif
}
return
nS3
xH
stacktop_before;}
}
#endif
#ifdef FP_SUPPORT_OPTIMIZER
eY
lL1
yB::tV2{using
tJ;iY
xE
tree;tree.GenerateFrom(*mData);FPoptimizer_Optimize::ApplyGrammars(tree);std::vector<x13>cG3;std::vector
yB
immed;cI3
stacktop_max=0;tree.SynthesizeByteCode(cG3,immed,stacktop_max);if(mData->mStackSize!=stacktop_max){mData->mStackSize=x13(stacktop_max);
#if !defined(FP_USE_THREAD_SAFE_EVAL) && \
    !defined(FP_USE_THREAD_SAFE_EVAL_WITH_ALLOCA)
mData->mStack
xA3
stacktop_max);
#endif
}
mData->mByteCode.swap(cG3);mData->mImmed.swap(immed);}
#define FUNCTIONPARSER_INSTANTIATE_EMPTY_OPTIMIZE(type) tH1>lL1<type>::tV2{}
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
tL3(MpfrFloat)
#endif
#ifdef FP_SUPPORT_GMP_INT_TYPE
tL3(GmpInt)
#endif
#ifdef FP_SUPPORT_COMPLEX_DOUBLE_TYPE
tL3(std::complex<double>)
#endif
#ifdef FP_SUPPORT_COMPLEX_FLOAT_TYPE
tL3(std::complex<float>)
#endif
#ifdef FP_SUPPORT_COMPLEX_LONG_DOUBLE_TYPE
tL3(std::complex<long
double>)
#endif
#define FUNCTIONPARSER_INSTANTIATE_OPTIMIZE(type) template lL1<type>::tV2;
#ifndef FP_DISABLE_DOUBLE_TYPE
tM3(double)
#endif
#ifdef FP_SUPPORT_FLOAT_TYPE
tM3(float)
#endif
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
tM3(long
double)
#endif
#ifdef FP_SUPPORT_LONG_INT_TYPE
tM3(long)
#endif
#endif // FP_SUPPORT_OPTIMIZER

#endif
