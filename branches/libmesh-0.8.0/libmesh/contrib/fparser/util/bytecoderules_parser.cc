static const char* const kVersionNumber = "1.0.0.0";

#include <ctype.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <map>
#include <set>

#include "../lib/crc32.hh"

//#define USE_CONTINUATIONS

namespace
{
    std::string trim(const std::string& s);

    const unsigned ShortLabelLength = 3u;

    struct Operation
    {
        enum { Opcode, ImmedFunc, OpcodeFunc } type;
        std::string result;
    };
    struct Match
    {
        enum { FixedOpcode, Immed, AnyOpcode } type;
        std::string name;         // opcode name such as cInt, or holder name such as x or X
        std::string condition;    // condition that applies when Immed or AnyOpcode

        bool operator==(const Match& b) const
        {
            return type==b.type
                && name==b.name
                && condition==b.condition;
        }

        std::vector<Operation> operations;
        bool has_operations;
        unsigned defined_on_line;
    };
    struct Node
    {
        Match opcode;
        std::vector<Node*> predecessors;
    };

    Node global_head;
    std::set<std::string> PossiblyUnusedLabelList;
    unsigned DefaultLabelCounter;

    std::string trim(const std::string& s)
    {
        size_t rm_begin = 0;
        while(rm_begin < s.size() && s[rm_begin] == ' ') ++rm_begin;
        size_t rm_end = 0;
        while(rm_end < s.size()-rm_begin && s[s.size()-1-rm_end] == ' ') ++rm_end;
        std::string result (s, rm_begin, s.size()-rm_end-rm_begin );
        do {
            std::string::size_type spos = result.find(' ');
            if(spos == result.npos) break;
            result.erase(spos, 1);
        } while(true);
        return result;
    }

    std::string Indent(size_t ntabs, size_t nspaces=0)
    {
        return std::string(ntabs, '\t') + std::string(nspaces, ' ');
    }

    std::string Bexpr(size_t pos)
    {
        if(pos == 0) return "opcode";
        std::ostringstream tmp;
        tmp << "ByteCodePtr[" << -int(pos-1) << "]";
        return tmp.str();
    }

    std::string Iexpr(size_t pos)
    {
        std::ostringstream tmp;
        tmp << "ImmedPtr[" << -int(pos) << "]";
        return tmp.str();
    }

    std::string BexprName(size_t pos)
    {
        if(pos == 0) return "opcode";
        std::ostringstream tmp;
        tmp << "op_" << pos;
        return tmp.str();
    }

    bool HasHandlingFor(const std::string& opcode)
    {
        if(!(opcode[0] == 'c' && isupper(opcode[1])))
            return true;

        if(opcode == "cImmed" || opcode == "cDup") return false;

        for(size_t b=0; b<global_head.predecessors.size(); ++b)
            if(global_head.predecessors[b]->opcode.type == Match::FixedOpcode
            && global_head.predecessors[b]->opcode.name == opcode)
                return true;

        for(size_t b=0; b<global_head.predecessors.size(); ++b)
            if(global_head.predecessors[b]->opcode.type == Match::AnyOpcode)
                return true;
        return false;
    }
    bool HasTailCallFor(const std::string& opcode)
    {
        for(size_t b=0; b<global_head.predecessors.size(); ++b)
            if(global_head.predecessors[b]->opcode.type == Match::FixedOpcode
            && global_head.predecessors[b]->opcode.name == opcode)
                return true;
        return false;
    }

    struct LineSeqs
    {
        void OutChain(std::ostream& out,
                      const std::vector<std::string>& chain,
                      size_t indent)
        {
            std::string codehash, prevlabel;
            for(size_t a=chain.size(); a-- > 0; )
            {
                codehash += chain[a];
                std::string label = GenLabel(codehash);
                Chains::iterator i = code.lower_bound(label);
                if(i != code.end() && i->first == label)
                {
                    std::string got_chain;
                    for(Chains::iterator j = i; j != code.end(); )
                    {
                        got_chain.insert(0, j->second.code);
                        if(j->second.nextlabel.empty()) break;
                        j = code.find(j->second.nextlabel);
                    }
                    /*
                    std::cerr << "expected: <" << codehash << ">\n"
                                 "got:      <" << got_chain << ">\n",
                    */
                    assert(got_chain == codehash);

                    // nothing to do
                    i->second.n_uses += 1;
                    if(!prevlabel.empty())
                        code.find(prevlabel)->second.n_uses -= 1;
                }
                else
                {
                    ChainItem item;
                    item.code         = chain[a];
                    item.nextlabel    = prevlabel;
                    item.n_uses       = 1;
                    code.insert(i, std::make_pair(label, item));
                }
                prevlabel = label;
            }
            if(!prevlabel.empty())
            {
                heads.push_back(prevlabel);
                out << Indent(indent) << "goto " << ChangeLabel(prevlabel) << ";\n";
            }
        }
        void Flush(std::ostream& out)
        {
            std::set<std::string> done;
            std::vector<std::string> remain;

            for(size_t a=0; a<heads.size(); ++a)
            {
                if(done.find(heads[a]) != done.end())
                    continue;

                size_t mini = 0;
                remain.push_back(heads[a]);
                while(!remain.empty())
                {
                    Chains::const_iterator i = code.find(remain.back());
                    remain.pop_back();

                    const ChainItem& item = i->second;
                    done.insert(i->first);

                    if(item.n_uses > mini)
                    {
                        std::string l = ChangeLabel(i->first);
                        out << l << ": ";
                    }
                    else
                        out << std::string(ShortLabelLength+2, ' ');
                    //out << " /*" << item.n_uses << "*/ ";
                    out << item.code;
                    if(!item.nextlabel.empty())
                    {
                        if(done.find(item.nextlabel) != done.end())
                        {
                            std::string l = ChangeLabel(item.nextlabel);
                            out << " goto " << l << ';';
                        }
                        else
                            remain.push_back(item.nextlabel);
                    }
                    else if(item.code.compare(0,5,"goto ") != false)
                        out << " return;";
                    out << "\n";

                    mini = 1;
                }
            }
        }
        void clear()
        {
            heads.clear();
            label_trans.clear();
            code.clear();
        }
    private:
        std::string GenLabel(unsigned crc, unsigned len)
        {
            /* If you get duplicate label errors from the resulting code
             * and the only thing you have done is add more optimization
             * rules, just uncomment one of these characters to enable
             * its use in the label name. If they are all already in use,
             * you will need to increase the ShortLabellength setting by one.
             */
            static const char table[] =
                /*"0123456789"*/
                /*"ABCDEFGHIJKLMNOPQRSTUVWXYZ"*/
                "abcdefghijklmnopq"
                /*"rstuvwxyz_"*/;
            char result[16] = {0};
            int o=15;
            while(true)
            {
                result[--o] = table[crc % (sizeof(table)-1)];
                crc /= (sizeof(table)-1);
                if(crc == 0 && (16-o) > 2) break;
            }
            result[--o] = 'L';
            std::string result2(result+o);
            result2.resize(len,' ');
            return result2;
        }
        std::string GenLabel(const std::string& code)
        {
            return GenLabel(crc32::calc((const unsigned char*)&code[0], code.size()), 20);
        }
        std::string ChangeLabel(const std::string& orig)
        {
            std::map<std::string, unsigned>::iterator
                j = label_trans.lower_bound(orig);
            if(j != label_trans.end() && j->first == orig)
                return GenLabel(j->second, ShortLabelLength);
            size_t lno = label_trans.size();
            label_trans.insert(j, std::make_pair(orig, lno));
            return GenLabel(unsigned(lno), ShortLabelLength);
        }
    private:
        struct ChainItem
        {
            std::string code;
            std::string nextlabel;
            unsigned n_uses;
        };
        typedef std::map<std::string/*label*/, ChainItem> Chains;
        std::vector<std::string> heads;
        std::map<std::string, unsigned> label_trans;
        Chains code;
    } CodeSeq;

    struct OutCode;
    struct OutLine
    {
        OutLine(OutCode& o) : out(o) { }
        template<typename T>
        OutLine& operator<< (const T& b) { buf << b; return *this; }
        ~OutLine();
    private:
        OutCode& out;
        std::ostringstream buf;
    };
    struct OutCode
    {
        OutCode(std::ostream& o, size_t i)
            : out(o),
              indent(i) { }
        ~OutCode()
        {
            for(size_t a=0; a<seq.size(); )
            {
                if((seq[a][0] == '/' && seq[a][1] == '*')
                || (seq[a][0] == 'n' && seq[a][1] == '_')
                || (seq[a][0] == 'r' && seq[a][1] == 'e' && seq[a][2] == 'p')
                || (seq[a][0] == 'o' && seq[a][1] == 'p' && seq[a][2] == '_')
                  )
                {
                    out << Indent(indent) << seq[a] << "\n";
                    seq.erase(seq.begin()+a);
                }
                else ++a;
            }
            #ifdef USE_CONTINUATIONS
            if(!seq.empty() && seq.back().compare(0,9,"goto Tail")==false)
            {
                size_t skip = 0;
                for(size_t a=0; a<seq.size(); ++a)
                {
                    if(seq[a].compare(0,8,"opcode =") == false)
                    {
                        std::string opcode_assign = seq[a];
                        seq.erase(seq.begin()+a);
                        seq.insert(seq.begin(), opcode_assign);
                        skip = 1;
                        break;
                    }
                }
                if(seq.size() > skip+1)
                {
                    //const std::string tail_label = seq.back().substr(5);
                    //const std::string tail_opcode = tail_label.substr(9, tail_label.size()-10);
                    //CodeSeq.AddContinuation(tail_label, tail_opcode);
                    //std::string continuation_line =
                    //    "/""* Will tailcall " + tail_opcode + " *""/";
                    //seq.insert(seq.begin()+skip, continuation_line);
                    seq.back() = "goto PickContinuation;";
                }
            }
            #endif
            CodeSeq.OutChain(out, seq, indent);
            if(seq.empty())
                out << Indent(indent) << "return;\n";
        }
        void Flush()
        {
            for(size_t a=0; a<seq.size(); ++a)
                out << Indent(indent) << seq[a] << "\n";
            seq.clear();
        }
        bool HasOperations() const { return !seq.empty(); }
    private:
        friend struct OutLine;
        void DidLine(const std::string& line) { seq.push_back(line); }
        std::ostream& out;
        std::vector<std::string> seq;
        size_t indent;
    };
    OutLine::~OutLine()
    {
        out.DidLine(buf.str());
    }

    struct Synther
    {
        Synther(OutCode& o, size_t i)
            : Out(o),
              know_bytecode_offset(true),
              know_immed_offset(true),
              indent(i)
        {
        }

        /* Reset is called before an operation that requires
         * that the vector's size() reflects the state shown
         * in the Ptr variable.
         * It is assumed that what follows is an instruction
         * that may reallocate the vector and invalidate pointers.
         */

        void ResetImmed(int offset = 0)
        {
        #if 0
            if(know_immed_offset)
            {
                OutLine(Out)  << "mData->mImmed.resize( " << (1-offset) << " + ImmedPtr - &mData->mImmed[0] );";
            }
            know_immed_offset = false;
        #else
            know_immed_offset = false;
            PopImmedBy(offset);
        #endif
        }

        void ResetByteCode(int offset = 0)
        {
        #if 0
            if(know_bytecode_offset)
            {
                OutLine(Out)  << "mData->mByteCode.resize( " << (1-offset) << " + ByteCodePtr - &mData->mByteCode[0] );";
            }
            know_bytecode_offset = false;
        #else
            know_bytecode_offset = false;
            PopByteCodeBy(offset);
        #endif
        }

        void ResetBoth(int b_offset, int i_offset)
        {
            ResetImmed(i_offset);
            ResetByteCode(b_offset);
        }

        void PopByteCodeBy(int n)
        {
        #if 0
            if(know_bytecode_offset)
                OutLine(Out)  << "ByteCodePtr -= " << n << ";";
            else
                for(; n > 0; --n)
                    OutLine(Out)  << "mData->mByteCode.pop_back();";
        #else
          #if 0
            for(; n > 0; --n)
                OutLine(Out)  << "mData->mByteCode.pop_back();";
          #else
            if(n == 1)
                OutLine(Out)  << "mData->mByteCode.pop_back();";
            else if(n > 0)
            {
                OutLine(Out) << "for(unsigned tmp=" << n << "; tmp-->0; ) mData->mByteCode.pop_back();";
            }
          #endif
            if(know_bytecode_offset)
            {
                //OutLine(Out)  << "if(mData->mByteCode.empty()) ByteCodePtr = 0; else ByteCodePtr -= " << n << ";";
                OutLine(Out)  << "ByteCodePtr -= " << n << ";";
            }
        #endif
        }

        void PopImmedBy(int n)
        {
        #if 0
            if(know_immed_offset)
                OutLine(Out)  << "ImmedPtr -= " << n << ";";
            else
                for(; n > 0; --n)
                    OutLine(Out)  << "mData->mImmed.pop_back();";
        #else
            if(know_immed_offset)
                OutLine(Out)  << "ImmedPtr -= " << n << ";";
         #if 0
            for(; n > 0; --n)
                OutLine(Out)  << "mData->mImmed.pop_back();";
         #else
            if(n == 1)
                OutLine(Out)  << "mData->mImmed.pop_back();";
            else if(n > 0)
            {
                OutLine(Out) << "for(unsigned tmp=" << n << "; tmp-->0; ) mData->mImmed.pop_back();";
            }
         #endif
        #endif
        }

        bool KnowOffsets() const
        {
            return know_bytecode_offset && know_immed_offset;
        }

    private:
        OutCode& Out;
        bool know_bytecode_offset, know_immed_offset;
        size_t indent;
    };

    bool SynthOperations(
        size_t indent, std::ostream& outstream,
        const std::vector<Match>& so_far,
        const Match& src_node,
        size_t b_used,
        size_t i_used)
    {
        std::vector<Operation> operations ( src_node.operations );

        outstream
            << Indent(indent)
            << "FP_TRACE_BYTECODE_OPTIMIZATION(" << src_node.defined_on_line << ',';

        unsigned n_with_lines = 0;
        std::ostringstream trace_with;
        for(size_t a=0; a<so_far.size(); ++a)
        {
            if(so_far[a].type != Match::Immed
            && so_far[a].type != Match::AnyOpcode) continue;
            if(n_with_lines == 0)
                trace_with << ",\n" << Indent(indent+1) << "\"    with ";
            else
                trace_with << "\n" << Indent(indent+1,4) << "<< \", ";
            trace_with << so_far[a].name << " = \" << ";
            if(so_far[a].type == Match::Immed)
                trace_with << so_far[a].name;
            else
                trace_with << "FP_TRACE_OPCODENAME(" << so_far[a].name << ")";
            ++n_with_lines;
        }
        if(n_with_lines > 0)
            outstream << "\n" << Indent(indent+1);

        outstream << '"';
        for(size_t a=so_far.size(); a-- > 0; )
        {
            if(a+1 != so_far.size()) outstream << " ";
            outstream << so_far[a].name;
            if(!so_far[a].condition.empty())
                outstream << "[" << so_far[a].condition << "]";
        }
        if(n_with_lines > 0)
            outstream << "\",\n" << Indent(indent+1) << "\"";
        else
            outstream << "\", \"";
        for(size_t a=0; a<operations.size(); ++a)
        {
            if(a > 0) outstream << ' ';
            switch(operations[a].type)
            {
                case Operation::ImmedFunc:
                    outstream << '[' << operations[a].result << ']';
                    break;
                case Operation::OpcodeFunc:
                    outstream << '{' << operations[a].result << '}';
                    break;
                case Operation::Opcode:
                    outstream << operations[a].result;
            }
        }
        outstream << "\"";
        if(n_with_lines == 0)
            outstream << ", \"\"";
        else if(n_with_lines == 1)
            trace_with << " << \"\\n\"";
        else
            trace_with << "\n" << Indent(indent+1,4) << "<< \"\\n\"";
        outstream << trace_with.str();
        outstream << ");\n";

        if(!operations.empty() && operations[0].result == "DO_POWI")
        {
            outstream
                << Indent(indent) << "if(TryCompilePowi(" << so_far.back().name << "))\n"
                << Indent(indent) << "  return;\n";
            return false;
        }

        OutCode Out(outstream, indent);

        int n_b_exist  = (int)(b_used-1);
        int n_i_exist  = (int)(i_used  );

        int b_offset = n_b_exist;
        int i_offset = n_i_exist;

        Synther offset_synth(Out, indent);

        bool rep_v_used = false;
        bool op_v_used  = false;

        bool changed = false;
        for(size_t a=0; a<operations.size(); ++a)
        {
            std::string opcode = trim(operations[a].result);
            if(opcode == "DO_STACKPLUS1")
            {
                outstream
                    << Indent(indent) << "incStackPtr();\n"
                    << Indent(indent) << "--mStackPtr;\n";
                continue;
            }

            if(operations[a].type == Operation::ImmedFunc)
            {
                //bool requires_var = false;
                //for(size_t a=0; a<opcode.size(); ++a)
                //    if(isalpha(opcode[a]))
                //        requires_var = true;

                if(i_offset > 0)
                {
                    bool i_redundant = false;
                    const Match& m = so_far[i_offset];
                    if(m.type == Match::Immed && opcode == m.name)
                        i_redundant = true;
                    if(i_redundant)
                    {
                        //requires_var = false;
                        OutLine(Out)
                            << "/* " << Iexpr(i_offset-1) << " = " << opcode << "; */"
                            << " // redundant, matches " << so_far[i_offset].name
                            << " @ " << (i_offset);
                    }
                    else
                    {
                        if(rep_v_used || true)
                            OutLine(Out)  << Iexpr(i_offset-1) << " = " << opcode << ";";
                        else
                        {
                            OutLine(Out)  << "rep_v = " << opcode << ";";
                            OutLine(Out)  << Iexpr(i_offset-1) << " = rep_v;";
                            rep_v_used = true;
                        }
                        changed = true;
                    }
                }
                else
                {
                    offset_synth.ResetImmed();
                    if(rep_v_used || true)
                        OutLine(Out)  << "mData->mImmed.push_back(" << opcode << ");";
                    else
                    {
                        OutLine(Out)  << "rep_v = " << opcode << ";";
                        OutLine(Out)  << "mData->mImmed.push_back(rep_v);";
                        rep_v_used = true;
                    }
                    changed = true;
                }
                //if(requires_var) { Out.Flush(); requires_var = false; }
                --i_offset;
                opcode = "cImmed";
            }

            bool redundant = false;
            if(b_offset > 0 && (!changed || !HasHandlingFor(opcode)))
            {
                const Match& m = so_far[b_offset];
                if(opcode == (m.type == Match::Immed ? "cImmed" : m.name))
                {
                    redundant = true;
                }
            }

            bool requires_var = !(opcode[0] == 'c' && isupper(opcode[1]));

            if(!redundant && HasHandlingFor(opcode))
            {
                if(a+1 == operations.size()
                && opcode[0] == 'c' && isupper(opcode[1])
                && HasTailCallFor(opcode))
                {
                    if(b_offset > 0 && i_offset > 0)
                        offset_synth.ResetBoth(b_offset, i_offset);
                    else
                    {
                        if(b_offset > 0) offset_synth.PopByteCodeBy(b_offset);
                        if(i_offset > 0) offset_synth.PopImmedBy(i_offset);
                    }
                    if(so_far[0].type == Match::FixedOpcode
                    && so_far[0].name == operations.back().result)
                    {
                        OutLine(Out)
                            << "/* opcode = " << operations.back().result << "; */"
                            << " // redundant, matches " << so_far[0].name << " @ 0";
                    }
                    else
                    {
                    #if 0
                        OutLine(Out)
                            << "/* opcode = " << operations.back().result << "; */"
                            << " // redundant, not really needed with tailcalls";
                    #else
                        OutLine(Out)  << "opcode = " << operations.back().result << ";";
                    #endif
                    }

                    if(!offset_synth.KnowOffsets())
                        OutLine(Out)  << "FP_ReDefinePointers();";
                    OutLine(Out) << "FP_TRACE_BYTECODE_ADD(" << opcode << ");";
                    OutLine(Out)  << "goto TailCall_" << opcode << ";";
                    PossiblyUnusedLabelList.erase("TailCall_" + opcode);
                    return true;
                }
                else
                {
                    offset_synth.ResetBoth(b_offset>0 ? b_offset : 0,
                                           i_offset>0 ? i_offset : 0);
                    if(op_v_used || true || opcode[0]!='c' || opcode=="cImmed")
                        OutLine(Out)  << "AddFunctionOpcode(" << opcode << ");";
                    else
                    {
                        OutLine(Out)  << "op_v = " << opcode << ";";
                        OutLine(Out)  << "AddFunctionOpcode(op_v);";
                        op_v_used = true;
                    }
                    //if(requires_var) { Out.Flush(); requires_var = false; }
                    i_offset = b_offset = 0;
                    changed = true;
                }
            }
            else
            {
                if(b_offset > 0)
                {
                    if(redundant)
                    {
                        requires_var = false;
                        OutLine(Out)
                            << "/* " << Bexpr(b_offset) << " = " << opcode << "; */"
                            << " // redundant, matches " << so_far[b_offset].name
                            << " @ " << (b_offset);
                    }
                    else
                    {
                        OutLine(Out)  << Bexpr(b_offset) << " = " << opcode << ";";
                        changed = true;
                    }
                }
                else
                {
                    offset_synth.ResetByteCode();
                    if(op_v_used || true || opcode[0]!='c' || opcode=="cImmed")
                        OutLine(Out)  << "mData->mByteCode.push_back(" << opcode << ");";
                    else
                    {
                        OutLine(Out)  << "op_v = " << opcode << ";";
                        OutLine(Out)  << "mData->mByteCode.push_back(op_v);";
                        op_v_used = true;
                    }
                    changed = true;
                }
                if(requires_var) { Out.Flush(); requires_var = false; }
                --b_offset;
            }
        }
        offset_synth.ResetBoth(b_offset,i_offset);

        if(!Out.HasOperations())
            outstream << "return;\n";
        return true;
    }

    std::set<std::string> declared;

    enum { mode_children = 1, mode_operations = 2 };
    bool Generate(
        const Node& head,
        const std::vector<Match>& so_far,
        size_t indent,
        std::ostream& declarations,
        std::ostream& code,
        size_t b_used,
        size_t i_used,
        int mode = mode_children+mode_operations)
    {
        if(!head.predecessors.empty() && (mode & mode_children))
        {
            std::string last_op_name = BexprName(b_used);
            std::string default_label_name;
            /*
            if(last_op_name != "opcode")
                code << Indent(indent) << "const unsigned " << last_op_name << " = " << Bexpr(b_used) << ";\n";
            code << Indent(indent) << "switch(" << last_op_name << ")\n";
            */
        #ifdef USE_CONTINUATIONS
            if(b_used == 0)
            {
                code << Indent(indent-1,2) << "PickContinuation:\n";
            }
        #endif
            bool needs_default_case = false;
            bool needs_immed_case   = false;
            std::set<std::string> other_cases;
            for(size_t b=0; b<head.predecessors.size(); ++b)
            {
                const Node& n = *head.predecessors[b];
                if(n.opcode.type == Match::AnyOpcode)
                    needs_default_case = true;
                else if(n.opcode.type == Match::Immed)
                    needs_immed_case = true;
                else
                    other_cases.insert(n.opcode.name);
            }
            bool needs_switchcase =
                (needs_immed_case+needs_default_case+other_cases.size()) > 1;

            if(needs_switchcase)
            {
                code << Indent(indent) << "switch(" << Bexpr(b_used) << ")\n";
                code << Indent(indent) << "{\n";
            }

            if(!other_cases.empty())
            {
                unsigned other_indent = needs_switchcase ? 1 : 1;

                for(std::set<std::string>::const_iterator
                    i = other_cases.begin();
                    i != other_cases.end();
                    ++i)
                {
                    const std::string& op = *i;

                #ifndef USE_CONTINUATIONS
                    if(b_used == 0)
                    {
                        code << Indent(indent) << "TailCall_" << op << ":\n";
                        PossiblyUnusedLabelList.insert("TailCall_" + op);

                        /* If ByteCodePtr is null, just add the opcode. */
                        code << Indent(indent) << "if(!ByteCodePtr)\n";
                        code << Indent(indent) << "{\n";
                        { OutCode Out(code, indent+other_indent);
                          Synther(Out, indent+other_indent).ResetBoth(0,0);
                          OutLine(Out) << "mData->mByteCode.push_back(opcode);";
                        }
                        code << Indent(indent) << "}\n";
                    }
                #endif
                    if(needs_switchcase)
                    {
                        code << Indent(indent) << "case " << op << ":\n";
                    }
                    else
                    {
                        code << Indent(indent) << "if(" << Bexpr(b_used) << " == " << op << ")\n"
                             << Indent(indent) << "{\n";
                    }

                    bool returned = false;

                    for(int round=0; round<4; ++round)
                        for(size_t a=0; a<head.predecessors.size(); ++a)
                        {
                            const Node& n = *head.predecessors[a];
                            if(round < 2  && n.opcode.has_operations) continue;
                            if(round >= 2 && !n.opcode.has_operations) continue;
                            if((round & 1) != !!n.opcode.condition.empty()) continue;

                            if(n.opcode.type == Match::FixedOpcode
                            && n.opcode.name == op)
                            {
                                //code << Indent(indent) << "  {\n";
                                std::vector<Match> ref(so_far);
                                ref.push_back(n.opcode);

                                if(n.opcode.condition.empty())
                                {
                                    returned |=
                                        Generate(n, ref, indent+other_indent, declarations,code, b_used+1, i_used, mode_children|(round>=2?mode_operations:0));
                                }
                                else
                                {
                                    code << Indent(indent+other_indent) << "if(" << n.opcode.condition << ")\n";
                                    code << Indent(indent+other_indent) << "{\n";
                                    Generate(n, ref, indent+other_indent+1, declarations,code, b_used+1, i_used, mode_children|(round>=2?mode_operations:0));
                                    code << Indent(indent+other_indent) << "}\n";
                                }
                            }
                        }
                    if(!returned)
                    {
                        if(needs_default_case)
                        {
                            if(default_label_name.empty())
                            {
                                std::ostringstream tmp;
                                tmp << "Default" << DefaultLabelCounter++;
                                default_label_name = tmp.str();
                            }
                            code << Indent(indent+other_indent) << "goto " << default_label_name << ";\n";
                        }
                        else if(needs_switchcase)
                        {
                            code << Indent(indent+other_indent) << "break;\n";
                        }
                    }
                    if(!needs_switchcase)
                        code << Indent(indent) << "}\n";
                }
            }

            std::set<std::string> immed_labels;
            bool immed_returned = false;
            /*
              Round 0:  no-code, has condition,  mode_children
              Round 1:  no-code, no condition,   mode_children
              Round 2:  code, has condition,     mode_children | mode_operations
              Round 3:  code, no condition,      mode_children | mode_operations
            */
            if(needs_immed_case)
            {
                if(needs_switchcase)
                {
                    code << Indent(indent) << "case cImmed:\n";
                }
                else
                {
                    code << Indent(indent) << "if(" << Bexpr(b_used) << " == cImmed)\n"
                         << Indent(indent) << "{\n";
                }
                unsigned imm_indent = needs_switchcase ? 1 : 1;

                for(int round=0; round<4; ++round)
                    for(size_t a=0; a<head.predecessors.size(); ++a)
                    {
                        const Node& n = *head.predecessors[a];
                        if(n.opcode.type == Match::Immed)
                        {
                            if(round < 2  && n.opcode.has_operations) continue;
                            if(round >= 2 && !n.opcode.has_operations) continue;
                            if((round & 1) != !!n.opcode.condition.empty()) continue;

                            //code << Indent(indent) << "  /* round " << round << " a = " << a << " */\n";
                            std::set<std::string>::iterator i = immed_labels.lower_bound(n.opcode.name);
                            if(i == immed_labels.end() || *i != n.opcode.name)
                            {
                                if(declared.find(n.opcode.name) == declared.end())
                                {
                                    declared.insert(n.opcode.name);
                                    declarations << Indent(1) << "Value_t " << n.opcode.name << ";\n";
                                }

                                code << Indent(indent+imm_indent) << n.opcode.name << " = " << Iexpr(i_used) << ";\n";
                                immed_labels.insert(i, n.opcode.name);
                            }
                            std::vector<Match> ref(so_far);
                            ref.push_back(n.opcode);
                            if(n.opcode.condition.empty())
                                immed_returned =
                                    Generate(n, ref, indent+imm_indent, declarations,code, b_used+1, i_used+1, mode_children|(round>=2?mode_operations:0));
                            else
                            {
                                code << Indent(indent+imm_indent) << "if(" << n.opcode.condition << ")\n";
                                code << Indent(indent+imm_indent) << "{\n";
                                Generate(n, ref, indent+imm_indent+1, declarations,code, b_used+1, i_used+1, mode_children|(round>=2?mode_operations:0));
                                code << Indent(indent+imm_indent) << "}\n";
                            }
                        }
                    }

                if(needs_switchcase)
                {
                    if(!immed_returned)
                        code << Indent(indent+imm_indent) << "break;\n";
                }
                else
                {
                    code << Indent(indent) << "}\n";
                }
            }

            std::set<std::string> opcode_labels;

            if(needs_default_case)
            {
                if(needs_switchcase)
                {
                    code << Indent(indent) << "default:";
                    if(!default_label_name.empty())
                        code << " " << default_label_name << ":;";
                    code << "\n";
                }
                unsigned default_indent = needs_switchcase ? 1 : 0;

                for(int round=0; round<4; ++round)
                    for(size_t a=0; a<head.predecessors.size(); ++a)
                    {
                        const Node& n = *head.predecessors[a];
                        if(n.opcode.type == Match::AnyOpcode)
                        {
                            if(round < 2  && n.opcode.has_operations) continue;
                            if(round >= 2 && !n.opcode.has_operations) continue;
                            if((round & 1) != !!n.opcode.condition.empty()) continue;
                            //code << Indent(indent) << "  /* round " << round << " a = " << a << " */\n";
                            std::set<std::string>::iterator i = opcode_labels.lower_bound(n.opcode.name);
                            if(i == opcode_labels.end() || *i != n.opcode.name)
                            {
                                if(declared.find(n.opcode.name) == declared.end())
                                {
                                    declared.insert(n.opcode.name);
                                    declarations << Indent(1) << "unsigned " << n.opcode.name << ";\n";
                                }
                                code << Indent(indent+default_indent) << n.opcode.name << " = " << Bexpr(b_used) << ";\n";
                                opcode_labels.insert(i, n.opcode.name);
                            }
                            std::vector<Match> ref(so_far);
                            ref.push_back(n.opcode);
                            if(n.opcode.condition.empty())
                                Generate(n, ref, indent+default_indent, declarations,code, b_used+1, i_used, mode_children|(round>=2?mode_operations:0));
                            else
                            {
                                code << Indent(indent+default_indent) << "if(" << n.opcode.condition << ")\n";
                                code << Indent(indent+default_indent) << "{\n";
                                Generate(n, ref, indent+default_indent+1, declarations,code, b_used+1, i_used, mode_children|(round>=2?mode_operations:0));
                                code << Indent(indent+default_indent) << "}\n";
                            }
                        }
                    }
            }

            if(needs_switchcase)
            {
                code << Indent(indent) << "}\n";
            }
        }
        if(head.opcode.has_operations && (mode & mode_operations))
        {
            /*if(!head.predecessors.empty())
                std::cout << Indent(indent) << "/""* NOTE: POSSIBLY AMBIGIOUS *""/\n";*/
            if(SynthOperations(indent,code, so_far, head.opcode, b_used, i_used))
            {
                //out << Indent(indent) << "return;\n";
                // ^ now redundant, as it is done by SynthOperations()
                return true;
            }
        }
        return false;
    }

    void Generate(std::ostream& out)
    {
        std::ostringstream code;
        std::ostream& declarations = out;
        Generate(global_head, std::vector<Match>(), 1, declarations,code, 0,0);
        out << code.str();

        { OutCode Out(out, 1);
          Synther(Out, 1).ResetBoth(0,0);
          OutLine(Out) << "mData->mByteCode.push_back(opcode);";
        }
        CodeSeq.Flush(out);

        out << "return;\n";
        out << "// This list of dummy gotos is here to inhibit\n"
               "// compiler warnings on unused labels\n";
        unsigned a=0;
        for(std::set<std::string>::const_iterator
            i = PossiblyUnusedLabelList.begin();
            i != PossiblyUnusedLabelList.end();
            ++i)
        {
            out << "goto " << *i << ";";
            if(a+1 == PossiblyUnusedLabelList.size()
            || a%3 == 2) out << "\n";
            ++a;
        }
    }

    bool isnamech(char c)
    {
        return c=='_' || (c>='a' && c<='z') || (c>='A' && c<='Z');
    }

    void VerifyExpressions(
        const std::set<std::string>& expressions,
        const std::set<std::string>& identifiers,
        unsigned lineno)
    {
        for(std::set<std::string>::const_iterator
            i = expressions.begin();
            i != expressions.end();
            ++i)
        {
            const std::string& e = *i;
            for(size_t a=0; a<e.size(); ++a)
            {
                if( isnamech(e[a])
                && !isnamech(e[a+1])
                && (a==0 || (!isnamech(e[a-1])
                          && !std::isdigit(e[a-1]))
                   ))
                {
                    std::string id(e, a, 1);
                    if(identifiers.find(id) == identifiers.end())
                    {
                        std::cerr
                            << "WARNING: Identifier '" << id
                            << "' used on line " << lineno
                            << " is undefined\n";
                    }
                }
            }
        }
    }

    template<typename T>
    bool SplitConditionString(const std::string& s, T& set, T& unset)
    {
        size_t a=0, b=s.size();
        bool glue_necessary = false;
        while(a < b)
        {
            if(s[a] == ' ') { ++a; continue; }
            if(s[a] == '(' || s[a] == ')') return false; // Cannot evaluate
            bool positive = true;
            if(glue_necessary)
            {
                if(s[a] != '&' || s[a+1] != '&') return false;
                a += 2;
                if(a >= b) return false;
            }
            while(a < b && s[a] == ' ') ++a;
            if(s[a] == '!') { ++a; positive = false; if(a >= b) return false; }
            while(a < b && s[a] == ' ') ++a;
            if(a >= b) return false;
            if((s[a] >= 'A' && s[a] <= 'Z')
            || (s[a] >= 'a' && s[a] <= 'z')
            || s[a]=='_')
            {
                size_t begin = a;
                while(++a < b && ((s[a] >= 'A' && s[a] <= 'Z')
                               || (s[a] >= 'a' && s[a] <= 'z')
                               || (s[a] >= '0' && s[a] <= '9')
                               || s[a]=='_'))
                    { }
                (positive ? set : unset).insert( s.substr(begin, a-begin) );
                glue_necessary = true;
                continue;
            }
            else
                break;
        }
        return true;
    }

    struct ParsingMode
    {
        std::set<std::string> different_preconditions;
        bool collect_preconditions;
        std::set<std::string> allowed_preconditions;
    };
    void Parse(ParsingMode& mode)
    {
        unsigned lineno = 0;
        while(true)
        {
            char Buf[2048];
            if(!std::fgets(Buf, sizeof(Buf), stdin)) break;
            ++lineno;
            char* bufptr = Buf;
            while(*bufptr == ' ' || *bufptr == '\t') ++bufptr;
            if(*bufptr == '#' || *bufptr == '\r' || *bufptr == '\n') continue;

            std::string Precondition;

            if(*bufptr == 'I' && bufptr[1] == 'F' && bufptr[2] == '(')
            {
                bufptr += 3;
                size_t balance = 0;
                while(*bufptr != ')' || balance != 0)
                {
                    if(*bufptr == '\r' || *bufptr == '\n') break;
                    if(*bufptr == '(') ++balance;
                    if(*bufptr == ')') --balance;
                    Precondition += *bufptr++;
                }
                if(*bufptr == ')') ++bufptr;
                while(*bufptr == ' ' || *bufptr == '\t') ++bufptr;
            }

            if(mode.collect_preconditions)
            {
                std::set<std::string> en, dis;
                if(SplitConditionString(Precondition, en, dis))
                {
                    mode.different_preconditions.insert(en.begin(), en.end());
                    mode.different_preconditions.insert(dis.begin(), dis.end());
                }
                else
                {
                    while(!Precondition.empty()
                        && Precondition[0] == '!') Precondition.erase(Precondition.begin());
                    if(!Precondition.empty())
                        mode.different_preconditions.insert(Precondition);
                }
                continue;
            }

            if(!Precondition.empty())
            {
                std::set<std::string> en, dis;
                /*fprintf(stderr, "Rule %u, Condition %s, allowed:",
                    lineno, Precondition.c_str());
                for(std::set<std::string>::const_iterator
                    i = mode.allowed_preconditions.begin();
                    i != mode.allowed_preconditions.end();
                    ++i)
                    fprintf(stderr, " %s", i->c_str());*/

                if(SplitConditionString(Precondition, en, dis))
                {
                    bool ok = true;
                    for(std::set<std::string>::const_iterator
                        i = en.begin(); i != en.end(); ++i)
                    {
                        /*fprintf(stderr, "[en:%s]", i->c_str());*/
                        if(mode.allowed_preconditions.find(*i)
                        == mode.allowed_preconditions.end())
                            { ok = false; break; }
                    }
                    for(std::set<std::string>::const_iterator
                        i = dis.begin(); i != dis.end(); ++i)
                    {
                        /*fprintf(stderr, "[dis:%s]", i->c_str());*/
                        if(mode.allowed_preconditions.find(*i)
                        != mode.allowed_preconditions.end())
                            { ok = false; break; }
                    }
                    /*fprintf(stderr, " - decision: %s\n", ok?"OK":" Not ok");
                    fflush(stderr);*/
                    if(!ok) continue;
                }
                else
                {
                    std::string cond = Precondition;
                    bool inverse = false;
                    while(cond[0] == '!')
                        { inverse = !inverse; cond.erase(cond.begin()); }
                    if((mode.allowed_preconditions.find(cond)
                     == mode.allowed_preconditions.end()) != inverse)
                    {
                        continue;
                    }
                }
            }

            std::set<std::string> identifiers_defined_here;
            std::set<std::string> expressions_used_here;

            std::vector<Match> sequence;
            while(true)
            {
                Match m;
                m.has_operations = false;
                while(*bufptr == ' ' || *bufptr == '\t') ++bufptr;
                if(*bufptr == '-' && bufptr[1] == '>') break;
                while(isalnum(*bufptr)) m.name += *bufptr++;
                m.name = trim(m.name);
                while(*bufptr == ' ' || *bufptr == '\t') ++bufptr;
                if(*bufptr == '#') break;
                if(*bufptr == '[')
                {
                    size_t balance = 0; ++bufptr;
                    while(*bufptr != ']' || balance != 0)
                    {
                        if(*bufptr == '\r' || *bufptr == '\n') break;
                        if(*bufptr == '[') ++balance;
                        if(*bufptr == ']') --balance;
                        m.condition += *bufptr++;
                    }
                    if(*bufptr == ']') ++bufptr;
                    m.condition = trim(m.condition);
                    expressions_used_here.insert(m.condition);
                }
                if(m.name[0] == 'c' && m.name.size() > 1)
                    m.type = Match::FixedOpcode;
                else if(isupper(m.name[0]))
                {
                    m.type = Match::AnyOpcode;
                    identifiers_defined_here.insert(m.name);
                }
                else
                {
                    m.type = Match::Immed;
                    identifiers_defined_here.insert(m.name);
                }

                sequence.push_back(m);
            }

            Node* head = &global_head;
            for(size_t b=sequence.size(); b-->0; )
            {
                const Match& m = sequence[b];
                bool dup = false;
                for(size_t a=0; a< head->predecessors.size(); ++a)
                {
                    if(m == head->predecessors[a]->opcode)
                    {
                        head = head->predecessors[a];
                        dup = true;
                        break;
                    }
                }
                if(!dup)
                {
                    Node* newhead = new Node;
                    newhead->opcode = m;
                    head->predecessors.push_back(newhead);
                    head = newhead;
                }
            }

            if(*bufptr == '-' && bufptr[1] == '>')
            {
                if(head->opcode.has_operations)
                {
                    std::cerr << "WARNING: Duplicate definition on line " << lineno
                              << ", previously defined on line " << head->opcode.defined_on_line
                              << std::endl;
                }
                head->opcode.has_operations = true;
                head->opcode.defined_on_line = lineno;
                bufptr += 2;

                while(true)
                {
                    while(*bufptr == ' ' || *bufptr == '\t') ++bufptr;
                    if(*bufptr == '#' || *bufptr == '\r' || *bufptr == '\n') break;
                    Operation op;
                    if(*bufptr == '[')
                    {
                        size_t balance = 0; ++bufptr;
                        while(*bufptr != ']' || balance != 0)
                        {
                            if(*bufptr == '\r' || *bufptr == '\n') break;
                            if(*bufptr == '[') ++balance;
                            if(*bufptr == ']') --balance;
                            op.result += *bufptr++;
                        }
                        if(*bufptr == ']') ++bufptr;
                        op.type = Operation::ImmedFunc;
                        expressions_used_here.insert(op.result);
                    }
                    else if(*bufptr == '{')
                    {
                        size_t balance = 0; ++bufptr;
                        while(*bufptr != '}' || balance != 0)
                        {
                            if(*bufptr == '\r' || *bufptr == '\n') break;
                            if(*bufptr == '{') ++balance;
                            if(*bufptr == '}') --balance;
                            op.result += *bufptr++;
                        }
                        if(*bufptr == '}') ++bufptr;
                        op.type = Operation::OpcodeFunc;
                        expressions_used_here.insert(op.result);
                    }
                    else
                    {
                        while(isalnum(*bufptr))
                            op.result += *bufptr++;
                        op.type = Operation::Opcode;
                    }
                    head->opcode.operations.push_back(op);
                }
            }

            VerifyExpressions(expressions_used_here,
                              identifiers_defined_here,
                              lineno);
        }
    }

    std::string GetPreconditionMaskName(
        const std::vector<std::string>& cond,
        size_t bitmask)
    {
        std::set<std::string> enabled;
        std::set<std::string> disabled;

        for(size_t b=0; b<cond.size(); ++b)
        {
            if(bitmask & (1 << b))
            {
                std::string c = cond[b];
                /*std::set<std::string> e, d;
                if(!SplitConditionString(c, e, d))*/
                    enabled.insert(c);
                /*else
                {
                    enabled.insert(e.begin(), e.end());
                    disabled.insert(d.begin(), d.end());
                }*/
            }
            else
            {
                disabled.insert(cond[b]);
                // No splitting of conditions here,
                // because !(a && b) yields !a || !b, not !a && !b
            }
        }

        std::string result;
        bool first=true;
        for(std::set<std::string>::const_iterator
            i = enabled.begin(); i != enabled.end(); ++i)
        {
            if(disabled.find(*i) != disabled.end())
                return std::string(); // conflicting conditions
            if(first) first=false; else result += " && ";
            result += "(" + *i + ")";
        }
        for(std::set<std::string>::const_iterator
            i = disabled.begin(); i != disabled.end(); ++i)
        {
            if(first) first=false; else result += " && ";
            result += "!(" + *i + ")";
        }
        return result;
    }
}

int main()
{
    ParsingMode mode;
    mode.collect_preconditions = true;
    Parse(mode);

    std::vector<std::string> different_preconditions(
        mode.different_preconditions.begin(),
        mode.different_preconditions.end() );

    size_t n_different_parsers = 1 << different_preconditions.size();

    std::ostream& out = std::cout;

    out <<
    //  "template<typename Value_t>\n"
    //  "inline void FunctionParserBase<Value_t>::AddFunctionOpcode(unsigned opcode)\n"
    //  "{\n"
        "  unsigned* ByteCodePtr;\n"
        "  Value_t*   ImmedPtr;\n"
        //"  unsigned n_b, n_i, op_v;\n"
        //"  Value_t  rep_v;\n"
        "\n"
        "  #define FP_ReDefinePointers() \\\n"
        "    ByteCodePtr = !mData->mByteCode.empty() ? &mData->mByteCode[0] + mData->mByteCode.size() - 1 : 0; \\\n"
        "    ImmedPtr    = !mData->mImmed.empty()    ? &mData->mImmed[0]    + mData->mImmed.size()    - 1 : 0;\n"
        "  FP_ReDefinePointers();\n";

    out << "  FP_TRACE_BYTECODE_ADD(opcode);\n";

    for(size_t n=0; n<n_different_parsers; ++n)
    {
        mode.collect_preconditions = false;
        mode.allowed_preconditions.clear();
        for(size_t b=0; b<different_preconditions.size(); ++b)
        {
            if(n & (1 << b))
                mode.allowed_preconditions.insert(different_preconditions[b]);
        }
        if(std::fseek(stdin, 0, SEEK_SET) < 0)
        {
            std::perror("fseek");
            return -1;
        }

        global_head = Node();
        PossiblyUnusedLabelList.clear();
        DefaultLabelCounter=0;
        CodeSeq.clear();
        declared.clear();

        Parse(mode);

        out << "\n";

        if(n_different_parsers == 1)
            Generate(out);
        else
        {
            std::string mask = GetPreconditionMaskName(different_preconditions, n);
            if(!mask.empty())
            {
                out << "#if(" << mask << ")\n" << std::flush;
                Generate(out);
                out << "#endif\n";
            }
        }
    }

    out <<
        "\n"
        "#undef FP_ReDefinePointers\n"
    //  "}\n"
        "#undef FP_TRACE_BYTECODE_OPTIMIZATION\n"
        "#undef FP_TRACE_OPCODENAME\n";

    return 0;
}
