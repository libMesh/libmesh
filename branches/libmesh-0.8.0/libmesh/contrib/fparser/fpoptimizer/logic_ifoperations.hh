#include "codetree.hh"

namespace
{
    using namespace FUNCTIONPARSERTYPES;
    using namespace FPoptimizer_CodeTree;

    /**************************************/
    /* IF OPERATIONS                      */
    /**************************************/

    template<typename Value_t>
    bool ConstantFolding_IfOperations(CodeTree<Value_t>& tree)
    {
        assert(tree.GetOpcode() == cIf || tree.GetOpcode() == cAbsIf);

        // If the If() condition begins with a cNot,
        // remove the cNot and swap the branches.
        for(;;)
        {
            if(tree.GetParam(0).GetOpcode() == cNot) // TEST 20/ifnot
            {
                tree.SetOpcode(cIf);
                tree.GetParam(0).Become( tree.GetParam(0).GetParam(0) );
                tree.GetParam(1).swap(tree.GetParam(2));
            }
            else if(tree.GetParam(0).GetOpcode() == cAbsNot) // TEST 20/ifabsnot
            {
                tree.SetOpcode(cAbsIf);
                tree.GetParam(0).Become( tree.GetParam(0).GetParam(0) );
                tree.GetParam(1).swap(tree.GetParam(2));
            }
            else break;
        }

        // If the sub-expression evaluates to approx. zero, yield param3.
        // If the sub-expression evaluates to approx. nonzero, yield param2.
        switch(GetLogicalValue(tree.GetParam(0), tree.GetOpcode()==cAbsIf))
        { // TEST 20/ifconst
            case IsAlways: // true
                tree.Become(tree.GetParam(1));
                return true; // rerun optimization (opcode changed)
            case IsNever: // false
                tree.Become(tree.GetParam(2));
                return true; // rerun optimization (opcode changed)
            case Unknown: default: ;
        }

        if(tree.GetParam(0).GetOpcode() == cIf
        || tree.GetParam(0).GetOpcode() == cAbsIf)
        {
            // TEST 20/ififconst
            //     if(if(x, a,b), c,d)
            //  -> if(x, if(a, c,d), if(b, c,d))
            // when either a or b is constantly true/false
            CodeTree<Value_t> cond = tree.GetParam(0);
            CodeTree<Value_t> truth_a;
            truth_a.SetOpcode(cond.GetOpcode() == cIf ? cNotNot : cAbsNotNot);
            truth_a.AddParam(cond.GetParam(1));
            ConstantFolding(truth_a);
            CodeTree<Value_t> truth_b;
            truth_b.SetOpcode(cond.GetOpcode() == cIf ? cNotNot : cAbsNotNot);
            truth_b.AddParam(cond.GetParam(2));
            ConstantFolding(truth_b);
            if(truth_a.IsImmed() || truth_b.IsImmed())
            {
                CodeTree<Value_t> then_tree;
                then_tree.SetOpcode(cond.GetOpcode());
                then_tree.AddParam(cond.GetParam(1));
                then_tree.AddParam(tree.GetParam(1));
                then_tree.AddParam(tree.GetParam(2));
                then_tree.Rehash();
                CodeTree<Value_t> else_tree;
                else_tree.SetOpcode(cond.GetOpcode());
                else_tree.AddParam(cond.GetParam(2));
                else_tree.AddParam(tree.GetParam(1));
                else_tree.AddParam(tree.GetParam(2));
                else_tree.Rehash();
                tree.SetOpcode(cond.GetOpcode());
                tree.SetParam(0, cond.GetParam(0));
                tree.SetParamMove(1, then_tree);
                tree.SetParamMove(2, else_tree);
                return true; // rerun cIf optimization
            }
        }
        if(tree.GetParam(1).GetOpcode() == tree.GetParam(2).GetOpcode()
        && (tree.GetParam(1).GetOpcode() == cIf
         || tree.GetParam(1).GetOpcode() == cAbsIf))
        {
            CodeTree<Value_t>& leaf1 = tree.GetParam(1);
            CodeTree<Value_t>& leaf2 = tree.GetParam(2);
            if(leaf1.GetParam(0).IsIdenticalTo(leaf2.GetParam(0))
            && (leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(1))
             || leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(2))))
            {
            // TEST 20/ifmerge
            //     if(x, if(y,a,b), if(y,c,d))
            // ->  if(y, if(x,a,c), if(x,b,d))
            // when either a,c are identical or b,d are identical
                CodeTree<Value_t> then_tree;
                then_tree.SetOpcode(tree.GetOpcode());
                then_tree.AddParam(tree.GetParam(0));
                then_tree.AddParam(leaf1.GetParam(1));
                then_tree.AddParam(leaf2.GetParam(1));
                then_tree.Rehash();
                CodeTree<Value_t> else_tree;
                else_tree.SetOpcode(tree.GetOpcode());
                else_tree.AddParam(tree.GetParam(0));
                else_tree.AddParam(leaf1.GetParam(2));
                else_tree.AddParam(leaf2.GetParam(2));
                else_tree.Rehash();
                tree.SetOpcode(leaf1.GetOpcode());
                tree.SetParam(0, leaf1.GetParam(0));
                tree.SetParamMove(1, then_tree);
                tree.SetParamMove(2, else_tree);
                return true; // rerun cIf optimization
            // cIf [x (cIf [y a z]) (cIf [y z b])] : (cXor x y) z (cIf[x a b])
            // ^ if only we had cXor opcode.
            }
            if(leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(1))
            && leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(2)))
            {
                // TEST 20/ifmerge2
                //    if(x, if(y,a,b), if(z,a,b))
                // -> if( if(x, y,z), a,b)
                CodeTree<Value_t> cond_tree;
                cond_tree.SetOpcode(tree.GetOpcode());
                cond_tree.AddParamMove(tree.GetParam(0));
                cond_tree.AddParam(leaf1.GetParam(0));
                cond_tree.AddParam(leaf2.GetParam(0));
                cond_tree.Rehash();
                tree.SetOpcode(leaf1.GetOpcode());
                tree.SetParamMove(0, cond_tree);
                tree.SetParam(2, leaf1.GetParam(2));
                tree.SetParam(1, leaf1.GetParam(1));
                return true; // rerun cIf optimization
            }
            if(leaf1.GetParam(1).IsIdenticalTo(leaf2.GetParam(2))
            && leaf1.GetParam(2).IsIdenticalTo(leaf2.GetParam(1)))
            {
                // TEST 20/ifmerge2b
                //    if(x, if(y,a,b), if(z,b,a))
                // -> if( if(x, y,!z), a,b)
                CodeTree<Value_t> not_tree;
                not_tree.SetOpcode(leaf2.GetOpcode() == cIf ? cNot : cAbsNot);
                not_tree.AddParam(leaf2.GetParam(0));
                not_tree.Rehash();
                CodeTree<Value_t> cond_tree;
                cond_tree.SetOpcode(tree.GetOpcode());
                cond_tree.AddParamMove(tree.GetParam(0));
                cond_tree.AddParam(leaf1.GetParam(0));
                cond_tree.AddParamMove(not_tree);
                cond_tree.Rehash();
                tree.SetOpcode(leaf1.GetOpcode());
                tree.SetParamMove(0, cond_tree);
                tree.SetParam(2, leaf1.GetParam(2));
                tree.SetParam(1, leaf1.GetParam(1));
                return true; // rerun cIf optimization
            }
        }

        CodeTree<Value_t>& branch1 = tree.GetParam(1);
        CodeTree<Value_t>& branch2 = tree.GetParam(2);

        if(branch1.IsIdenticalTo(branch2))
        {
            // TEST 20/ifnop
            // If both branches of an If() are identical, the test becomes unnecessary
            // NOTE: Possible side-effect on condition removed
            tree.Become(tree.GetParam(1));
            return true; // rerun optimization (opcode changed)
        }

        const OPCODE op1 = branch1.GetOpcode();
        const OPCODE op2 = branch2.GetOpcode();
        if(op1 == op2)
        {
            // TEST 20/if_extract_sin
            // TEST 20/if_extract_abs
            // If both branches apply the same unary function to different values,
            // extract the function. E.g. if(x,sin(a),sin(b)) -> sin(if(x,a,b))
            if(branch1.GetParamCount() == 1)
            {
                CodeTree<Value_t> changed_if;
                changed_if.SetOpcode(tree.GetOpcode());
                changed_if.AddParamMove(tree.GetParam(0));
                changed_if.AddParam(branch1.GetParam(0));
                changed_if.AddParam(branch2.GetParam(0));
                changed_if.Rehash();
                tree.SetOpcode(op1);
                tree.DelParams();
                tree.AddParamMove(changed_if);
                return true; // rerun optimization (opcode changed)
            }
            // TEST 20/if_extract_div
            // If both branches apply the same binary function to a set of parameters
            // where only one parameter differs, extract the function.
            // E.g. if(x, y/2, z/2) --> if(x, y,z)/2 (integer mode)
            if(branch1.GetParamCount() == 2
            && branch2.GetParamCount() == 2)
            {
                if(branch1.GetParam(0).IsIdenticalTo(branch2.GetParam(0)))
                {
                    CodeTree<Value_t> param0 = branch1.GetParam(0);
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.AddParamMove(tree.GetParam(0));
                    changed_if.AddParam(branch1.GetParam(1));
                    changed_if.AddParam(branch2.GetParam(1));
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.DelParams();
                    tree.AddParamMove(param0);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
                if(branch1.GetParam(1).IsIdenticalTo(branch2.GetParam(1)))
                {
                    CodeTree<Value_t> param1 = branch1.GetParam(1);
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.AddParamMove(tree.GetParam(0));
                    changed_if.AddParam(branch1.GetParam(0));
                    changed_if.AddParam(branch2.GetParam(0));
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.DelParams();
                    tree.AddParamMove(changed_if);
                    tree.AddParamMove(param1);
                    return true; // rerun optimization (opcode changed)
                }
            }
            // TEST 20/if_extract_add
            // TEST 20/if_extract_mul
            // TEST 20/if_extract_min
            if(op1 == cAdd    || op1 == cMul
            || op1 == cAnd    || op1 == cOr
            || op1 == cAbsAnd || op1 == cAbsOr
            || op1 == cMin    || op1 == cMax)
            {
                // If the two commutative groups contain one
                // or more identical values, extract them.
                std::vector<CodeTree<Value_t> > overlap;
                for(size_t a=branch1.GetParamCount(); a-- > 0; )
                {
                    for(size_t b=branch2.GetParamCount(); b-- > 0; )
                    {
                        if(branch1.GetParam(a).IsIdenticalTo(branch2.GetParam(b)))
                        {
                            if(overlap.empty()) { branch1.CopyOnWrite(); branch2.CopyOnWrite(); }
                            overlap.push_back(branch1.GetParam(a));
                            branch2.DelParam(b);
                            branch1.DelParam(a);
                            break;
                        }
                    }
                }
                if(!overlap.empty())
                {
                    branch1.Rehash();
                    branch2.Rehash();
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.SetParamsMove(overlap);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
            }
        }
        // if(x, y+z, y) -> if(x, z,0)+y
        if(op1 == cAdd
        || op1 == cMul
        || (op1 == cAnd && IsLogicalValue(branch2))
        || (op1 == cOr  && IsLogicalValue(branch2))
          )
        {
            // TEST 20/if_extract_add1
            // TEST 20/if_extract_mul1
            // TEST 20/if_extract_and1_l
            // TEST 20/if_extract_and1_nl
            // TEST 20/if_extract_or1_l
            // TEST 20/if_extract_or1_nl
            for(size_t a=branch1.GetParamCount(); a-- > 0; )
                if(branch1.GetParam(a).IsIdenticalTo(branch2))
                {
                    branch1.CopyOnWrite();
                    branch1.DelParam(a);
                    branch1.Rehash();
                    CodeTree<Value_t> branch2_backup = branch2;
                    branch2 = CodeTreeImmed( Value_t( (op1==cAdd||op1==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.AddParamMove(branch2_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, y&z, !!y) -> if(x, z, 1) & y
        if((op1 == cAnd || op1 == cOr) && op2 == cNotNot)
        {
            CodeTree<Value_t>& branch2op = branch2.GetParam(0);
            for(size_t a=branch1.GetParamCount(); a-- > 0; )
                if(branch1.GetParam(a).IsIdenticalTo(branch2op))
                {
                    branch1.CopyOnWrite();
                    branch1.DelParam(a);
                    branch1.Rehash();
                    CodeTree<Value_t> branch2_backup = branch2op;
                    branch2 = CodeTreeImmed( Value_t( (op1==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op1);
                    tree.AddParamMove(branch2_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, y, y+z) -> if(x, 0,z)+y
        if(op2 == cAdd
        || op2 == cMul
        || (op2 == cAnd && IsLogicalValue(branch1))
        || (op2 == cOr  && IsLogicalValue(branch1))
          )
        {
            // TEST 20/if_extract_add2
            // TEST 20/if_extract_mul2
            // TEST 20/if_extract_and2_l
            // TEST 20/if_extract_and2_nl
            // TEST 20/if_extract_or2_l
            // TEST 20/if_extract_or2_nl
            for(size_t a=branch2.GetParamCount(); a-- > 0; )
                if(branch2.GetParam(a).IsIdenticalTo(branch1))
                {
                    branch2.CopyOnWrite();
                    branch2.DelParam(a);
                    branch2.Rehash();
                    CodeTree<Value_t> branch1_backup = branch1;
                    branch1 = CodeTreeImmed( Value_t( (op2==cAdd||op2==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op2);
                    tree.AddParamMove(branch1_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        // if(x, !!y, y&z) -> if(x, 1, z) & y
        if((op2 == cAnd || op2 == cOr) && op1 == cNotNot)
        {
            CodeTree<Value_t>& branch1op = branch1.GetParam(0);
            for(size_t a=branch2.GetParamCount(); a-- > 0; )
                if(branch2.GetParam(a).IsIdenticalTo(branch1op))
                {
                    branch2.CopyOnWrite();
                    branch2.DelParam(a);
                    branch2.Rehash();
                    CodeTree<Value_t> branch1_backup = branch1op;
                    branch1 = CodeTreeImmed( Value_t( (op2==cOr) ? 0 : 1 ) );
                    CodeTree<Value_t> changed_if;
                    changed_if.SetOpcode(tree.GetOpcode());
                    changed_if.SetParamsMove(tree.GetParams());
                    changed_if.Rehash();
                    tree.SetOpcode(op2);
                    tree.AddParamMove(branch1_backup);
                    tree.AddParamMove(changed_if);
                    return true; // rerun optimization (opcode changed)
                }
        }
        return false; // No changes
    }
}
