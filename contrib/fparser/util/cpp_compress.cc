#include "cpp_compress.hh"

#include <vector>
#include <map>
#include <set>
#include <iostream>

#include <stdio.h>

namespace
{
    std::string macro_prefix_chars = "qghdowam";
    bool verbose = false;

    /* ^List of characters that are _least_ likely to form
     *  a legitimate identifier name in the input code when
     *  a [0-9A-Z] (case sensitive) is appended to it.
     */

    std::set<std::string> parametric_macro_list;

    void reset_parametric_macro_list()
    {
        static const char* const list[] =
        {
            /* IMPORTANT NOTE: HARDCODED LIST OF ALL PREPROCESSOR
             * MACROS THAT TAKE PARAMETERS. ANY EXTERNAL MACROS
             * NOT LISTED HERE WILL BE BROKEN AFTER COMPRESSION.
             * MACROS DEFINED WITHIN THE PROCESSED INPUT ARE FINE.
             */
            "assert", "getc", "putc",
            "FP_TRACE_OPCODENAME",
            "FP_TRACE_BYTECODE_OPTIMIZATION",
            "FP_TRACE_BYTECODE_ADD",
            "N","P", "DBL_ONLY", "LNG_ONLY",
            "incStackPtr", "TryCompilePowi","findName"
        };
        parametric_macro_list.clear();
        for(unsigned p=0; p<sizeof(list)/sizeof(*list); ++p)
            parametric_macro_list.insert(list[p]);
    }

    struct length_rec
    {
        unsigned begin_index;
        unsigned num_tokens;
        unsigned num_occurrences;
    };

    inline bool isnamechar(char c)
    {
        switch(c)
        {
            case '0':case 'a':case 'b':case 'c':case 'd':case 'e':case 'f':
            case '1':case 'g':case 'h':case 'i':case 'j':case 'k':case 'l':
            case '2':case 'm':case 'n':case 'o':case 'p':case 'q':case 'r':
            case '3':case 's':case 't':case 'u':case 'v':case 'w':case 'x':
            case '4':case 'y':case 'z':case '_':
            case '5':case 'A':case 'B':case 'C':case 'D':case 'E':case 'F':
            case '6':case 'G':case 'H':case 'I':case 'J':case 'K':case 'L':
            case '7':case 'M':case 'N':case 'O':case 'P':case 'Q':case 'R':
            case '8':case 'S':case 'T':case 'U':case 'V':case 'W':case 'X':
            case '9':case 'Y':case 'Z': return true;
        }
        return false;
    }

    class StringSeq: private std::string
    {
    private:
        std::vector<bool> sticky;
    public:
        StringSeq()
            : std::string(), sticky() { }

        StringSeq(const std::string& s)
            : std::string(s), sticky() { }

        const std::string& GetString() const { return *this; }

        template<typename T>
        void operator += (const T& k)
        {
            std::string::operator+=(k);
        }
        char operator[] (size_t pos) const
        {
            return std::string::operator[] (pos);
        }

        void SetSticky(size_t n_last_chars)
        {
            if(sticky.size() < size()) sticky.resize(size());
            for(size_t p = size() - n_last_chars; p < size(); ++p)
                sticky[p] = true;
        }
        bool IsSticky(size_t index) const
        {
            return index < sticky.size() && sticky[index];
        }
        bool HasSticky() const { return !sticky.empty(); }

        StringSeq substr(size_t begin, size_t count = std::string::npos) const
        {
            StringSeq result;
            result.std::string::assign(*this, begin,count);
            result.sticky.assign(
                begin < sticky.size()
                    ? sticky.begin() + begin
                    : sticky.end(),
                (count != std::string::npos && begin+count < sticky.size())
                    ? sticky.begin() + begin + count
                    : sticky.end());
            while(!result.sticky.empty() && !result.sticky.back())
                result.sticky.pop_back();
            return result;
        }

        using std::string::empty;
        using std::string::size;
        using std::string::append;
    };

    struct token
    {
        std::string value;
        struct token_meta
        {
            unsigned    hash;
            int         balance;
            unsigned    comma;
            bool        preproc;
        } meta;

        token(const std::string& v) : value(v)
        {
            Rehash();
        }

        token(const StringSeq& v) : value(v.GetString())
        {
            Rehash();
            if(v.HasSticky())
                meta.preproc = true;
        }

        void operator=(const std::string& v) { value=v; Rehash(); }

        bool operator==(const token& b) const
        {
            return meta.hash == b.meta.hash
                && value == b.value
                && meta.preproc == b.meta.preproc;
        }

        bool operator!=(const token& b) const
            { return !operator==(b); }

        void Rehash()
        {
            const char* v = value.data();
            meta.preproc = (v[0] == '#' || v[0] == '\n');
            meta.balance = 0;
            meta.comma   = 0;
            meta.hash    = 0;
            for(size_t a=0; a<value.size(); ++a)
            {
                //meta.hash = meta.hash*0x8088405u + v[a];
                meta.hash = ((meta.hash << (1)) | (meta.hash >> (32-1))) - v[a];
                switch(v[a])
                {
                    case '(': ++meta.balance; break;
                    case ')': --meta.balance; break;
                    case ',': ++meta.comma; break;
                }
            }
        }
        void swap(token& b)
        {
            value.swap(b.value);
            std::swap(meta, b.meta);
        }
    };
    bool Debug = false;
    StringSeq GetSeq(std::vector<token>::const_iterator begin,
                     size_t n, bool NewLines)
    {
        bool InDefineMode = false;

        /* Resequence the input */
        StringSeq result;
        int quotemode = 0;
        while(n-- > 0)
        {
            const token&       tok = *begin;
            const std::string& value = tok.value; ++begin;

            if(Debug)
                result += tok.meta.preproc
                           ? (InDefineMode ? "\xBF" : "\xA1")  // upside-down ? and !
                           : (InDefineMode ? "\x3F" : "\x21"); // ?, !
            else if (!result.empty() && result[result.size()-1] != '\n')
            {
                if (!InDefineMode && NewLines && value[0] == '#') result += '\n';
                else if (isnamechar(value[0])
                     && (isnamechar(result[result.size()-1])
                       || result[result.size()-1] == '"'
                       || result[result.size()-1] == '\''
                            //     || result[result.size()-1]==')'
                            /* An identifier preceded by an identifier character
                             * requires a separator. Also, an identifier
                             * preceded by a macro that expands to something
                             * that ends in an identifier character requires a
                             * separator when using Microsoft C++, thus you may
                             * need to include the ")" check above sometimes.
                             * Also, in C++11 you can't tack an identifier right
                             * after a string or character constant or there will
                             * be an error. 
                             */
                        ))
                {
                    if (!NewLines || InDefineMode)
                        result += ' ';
                    else
                        result += '\n';
                }
            }
            if (value[0] == '#')
            {
                if (value.substr(0,7) == "#define") InDefineMode = true;
                if (value.substr(0,8) == "#include") InDefineMode = true;
            }
            if (value == "\n")
            {
                if (InDefineMode) { result += "\n"; InDefineMode = false; }
                continue;
            }

            switch(quotemode)
            {
                case 0: // prev wasn't a quote, or this is not a quote
                    if (value[0] == '"'
                    && (n>0 && begin->value[0] == '"')) // this and next are quotes
                        { quotemode = 1;
                          result.append(value, 0, value.size()-1);
                          //result += value.substr(0, value.size()-1);
                          continue;
                        }
                    else
                    {
                        result += value;
                        if(tok.meta.preproc || (value[0] == '(' && value.size() > 1))
                        {
                            /* A macro parameter list is an undivisible entity */
                            result.SetSticky(value.size());
                        }
                    }
                    break;
                case 1: // prev was a quote, skip this quote
                    if (n>0 && begin->value[0] == '"')
                        { //result += value.substr(1, value.size()-2);
                          result.append(value, 1, value.size()-2);
                          continue;
                        }
                    else
                        { quotemode = 0;
                          //result += value.substr(1);
                          result.append(value, 1, value.size()-1);
                        }
                    break;
            }
            if (!Debug)
            {
                if (NewLines && !InDefineMode)
                {
                    if (value[0] == '#'
                     || value[0] == '}'
                     || value[0] == '"'
                      )
                    {
                        result += '\n';
                    }
                }
                if (n > 0 && (value.size() == 1 &&
                               (value[0] == '<'  // < < is not <<
                             || value[0] == '>'  // > > is not >>
                             || value[0] == '+'  // + + is not ++
                             || value[0] == '-'  // - - is not --
                             || value[0] == '&') // & & is not &&
                               ) && value == begin->value)
                {
                    result += ' ';
                }
            }
        }
        return result;
    }

    struct DefineParsingMode
    {
        bool Active;
        std::set<std::string> ParamList;

        DefineParsingMode() : Active(false), ParamList() { }

        void Activate() { Active=true; }
        void Deactivate() { Active=false; ParamList.clear(); }
    };

    std::vector<token>
        Tokenize(const StringSeq& input,
                 DefineParsingMode& defmode,
                 bool SplitMacros,
                 bool SplitStrings)
    {
        std::vector<token> result;
        size_t a=0, b=input.size();
        while(a < b)
        {
            /*if(input.IsSticky(a))
            {
                size_t eat = 1;
                while(input.IsSticky(a+eat)) ++eat;
                result.push_back( input.substr(a, eat) );
                a += eat;
                continue;
            }*/
            if (defmode.Active && input[a]=='\\' && input[a+1] == '\n')
            {
                a += 2;
                continue;
            }
            if (defmode.Active && input[a]=='\n')
            {
                defmode.Deactivate();
                result.push_back( token("\n") );
                ++a;
                continue;
            }
            if (input[a]==' ' || input[a]=='\t'
            || input[a]=='\n' || input[a]=='\r') { ++a; continue; }

            if (input[a]=='/' && input[a+1]=='*')
            {
                a += 2;
                while(a < b && (input[a-2]!='*' || input[a-1]!='/')) ++a;
                continue;
            }
            if (input[a]=='/' && input[a+1]=='/')
            {
                while(a < b && input[a]!='\n') ++a;
                continue;
            }
            if (input[a]=='_' || (input[a]>='a' && input[a]<='z')
                              || (input[a]>='A' && input[a]<='Z'))
            {
                size_t name_begin = a;
                while(++a < b)
                {
                    if (isnamechar(input[a])) continue;
                    break;
                }

                std::string name(input.GetString(), name_begin, a-name_begin);
                result.push_back(name);

                if (defmode.Active && defmode.ParamList.find(name)
                                   != defmode.ParamList.end())
                {
                    // Define params are immutable.
                    result.back().meta.preproc = true;
                }

                if (input[a] == '('
                && parametric_macro_list.find(name)
                != parametric_macro_list.end())
                {
                    std::vector<token> remains = Tokenize(input.substr(a), defmode, SplitMacros, SplitStrings);
                    int balance = 1;
                    size_t eat = 1;
                    for(; eat < remains.size() && balance != 0; ++eat)
                        balance += remains[eat].meta.balance;
                    if(SplitMacros)
                    {
                        for(size_t c=0; c<eat; ++c)
                            if(remains[c].meta.balance != 0
                            || remains[c].meta.comma != 0)
                                remains[c].meta.preproc = true;
                        result.insert(result.end(), remains.begin(), remains.end());
                    }
                    else
                    {
                        //result.push_back( GetSeq(remains.begin(), eat, false) );
                        StringSeq tmp = GetSeq(remains.begin(), eat, false);
                        result.back() = result.back().value + tmp.GetString();
                        result.insert(result.end(), remains.begin()+eat, remains.end());
                    }
                    a = b; // done
                }
                continue;
            }
            if (std::isdigit(input[a]) ||
               (input[a] == '.' && std::isdigit(input[a+1])))
            {
                size_t value_begin = a;
                while(++a < b)
                {
                    if ((input[a]>='0' && input[a]<='9')
                    || input[a]=='.' || input[a]=='+' || input[a]=='-'
                    || input[a]=='x' || (input[a]>='a' && input[a]<='f')
                    || input[a]=='p' || (input[a]>='A' && input[a]<='F')
                    || input[a]=='u' || input[a]=='U'
                    || input[a]=='l' || input[a]=='f'
                    || input[a]=='L' || input[a]=='F') continue;
                    break;
                }
                StringSeq s = input.substr(value_begin, a-value_begin );
                /* TODO: Convert hex to decimal */
                result.push_back( s );
                continue;
            }
            if (a+1 < b && input[a] == '>' && input[a+1] == '>')
                { int n = (a+2 < b && input[a+2] == '=') ? 3 : 2;
                  result.push_back(input.substr(a, n)); a += n; continue; }
            if (a+1 < b && input[a] == '<' && input[a+1] == '<')
                { int n = (a+2 < b && input[a+2] == '=') ? 3 : 2;
                  result.push_back(input.substr(a, n)); a += n; continue; }
            if (a+1 < b && input[a] == '+' && input[a+1] == '+')
                { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (a+1 < b && input[a] == '-' && input[a+1] == '-')
                { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (a+1 < b && input[a] == '-' && input[a+1] == '>')
                { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (a+1 < b && input[a] == '#' && input[a+1] == '#')
                { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (a+1 < b && input[a] == '&' && input[a+1] == '&')
                { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (a+1 < b && input[a] == '|' && input[a+1] == '|')
                { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (a+1 < b && (input[a] == '>' || input[a] == '<'
                        || input[a] == '!' || input[a] == '='
                        || input[a] == '+' || input[a] == '-'
                        || input[a] == '*' || input[a] == '/'
                        || input[a] == '&' || input[a] == '|'))
                if (input[a+1] == '=')
                    { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (a+1 < b && (input[a] == ':' && input[a+1] == ':'))
                    { result.push_back(input.substr(a, 2)); a += 2; continue; }
            if (!defmode.Active && input[a] == '#')
            {
                if (input.substr(a,8).GetString() == "#include")
                {
                    size_t p = a;
                    while(p < b && input[p] != '\n') ++p;
                    result.push_back( input.substr(a, p-a) );
                    result.back().meta.preproc = true;
                    result.push_back( token("\n") );
                    a = p;
                    continue;
                }
                if (input.substr(a,7).GetString() == "#define")
                {
                    size_t p = a+7;
                    while(p < b && std::isspace(input[p])) ++p; // skip blank
                    size_t term_begin = p;
                    while(p < b && isnamechar(input[p])) ++p; // skip term
                    if (input[p] != '(')
                    {
                        defmode.Activate();
                        std::string def = input.substr(a, p-a).GetString();
                        if(input[p] != '\n') def += ' ';
                        result.push_back(def); /* #define, term name and a space */
                        a = p;
                        continue;
                    }
                    else
                    {
                        std::string macro(input.GetString(), term_begin, p-term_begin);
                        if(parametric_macro_list.find(macro)
                        == parametric_macro_list.end())
                        {
                            if(verbose) std::cerr << "Detected parametric macro: " << macro << " (ok)\n";
                            parametric_macro_list.insert(macro);
                        }

                        size_t /*param_list_begin = p,*/ param_begin = p+1;
                        int balance = 1;
                        for (++p; true; ++p)
                        {
                            if(input[p]=='(') ++balance;
                            if(input[p]==',' || input[p]==')')
                            {
                                std::string paramname = input.substr(param_begin,p-param_begin).GetString();
                                //std::cerr << "Param name<" << paramname << ">\n";
                                defmode.ParamList.insert(paramname);
                                param_begin = p+1;
                            }
                            if(input[p]==')') { if(--balance == 0) { ++p; break; } }
                        }
                        size_t param_list_end = p;
                        while(input[p] != '\n' && std::isspace(input[p])) ++p;

                        defmode.Activate();
                        std::string def = input.substr(a, param_list_end-a).GetString();
                        if(input[p] != '\n') def += ' ';
                        result.push_back(def); /* #define, term name, params and a space */
                        a = p;
                        continue;
                    }
                }
                size_t preproc_begin = a;
                bool in_quotes = false;
                while(++a < b)
                {
                    if (!in_quotes && input[a]=='"')
                        { in_quotes=true; continue; }
                    if (in_quotes && input[a]=='"' && input[a-1]!='\\')
                        { in_quotes=false; continue; }
                    if (input[a]=='\\' && input[a+1]=='\n') { ++a; continue; }
                    if (input[a]=='\n') break;
                }
                std::string stmt(input.GetString(), preproc_begin, a-preproc_begin);
                if (stmt.substr(0,5) != "#line")
                    result.push_back(stmt);
                result.push_back( token("\n") );
                continue;
            }
            if (input[a] == '"')
            {
                size_t string_begin = a;
                while(++a < b)
                    if (input[a]=='"' &&
                      (input[a-1] != '\\'
                     || input[a-2]=='\\')) { ++a; break; }
                if(input.GetString() == "\"\"") continue; // Don't add an empty string token

                std::string stringconst(input.GetString(), string_begin+1, a-string_begin-2);
                if (SplitMacros)
                    while (true)
                    {
                        size_t p = stringconst.find_first_of(" ,+-", 0, 1);
                        if(p == stringconst.npos) break;
                        if(p > 0)
                            result.push_back( "\""+std::string(stringconst,0,p)+"\"" );
                        result.push_back( "\""+std::string(stringconst,p,1)+"\"" );
                        stringconst.erase(0, p+1);
                    }
                if(!stringconst.empty()) result.push_back("\""+stringconst+"\"");
                continue;
            }
            if (input[a] == '\'')
            {
                size_t char_begin = a; a += 3;
                if (input[a-2] == '\\') ++a;
                result.push_back( std::string(input.GetString(), char_begin, a-char_begin) );
                continue;
            }

            result.push_back( input.substr(a++, 1) );
        }
        return result;
    }

    static const char cbuf[] =
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_";

    unsigned macro_counter = 0;

    std::string GenerateMacroName()
    {
        std::string result;
        unsigned p = macro_counter++;
        result += macro_prefix_chars[(p/36)%macro_prefix_chars.size()];
        result += cbuf[p%36];
        p /= unsigned(macro_prefix_chars.size()*36); // 0-9A-Z
        for(; p != 0; p /= 63)
            result += cbuf[p%63];
        return result;
    }

    bool CompressWithNonparametricMacros
        (std::vector<token>& tokens,
         const std::string& seq_name_buf)
    {
        size_t seq_name_length = seq_name_buf.size();

        /* Find a sub-sequence of tokens for which
         * the occurrence-count times total length is
         * largest and the balance of parentheses is even.
         */
        std::map<unsigned, length_rec> hash_results;
        long best_score=0;
        //size_t best_score_length=0;
        unsigned best_hash=0;

        if(verbose) std::cerr << tokens.size() << " tokens\n";

        std::vector<bool> donttest(tokens.size(), false);

        const size_t lookahead_depth = (140*100000) / tokens.size();
        // ^Lookahead limit. The value chosen is an arbitrary value
        //  to shield against exponential runtime. A larger value
        //  yields better compression, but is slower.

        for(size_t a=0; a<tokens.size(); ++a)
        {
            if (donttest[a]) continue;

            //std::cerr << a << '\t' << best_score << '\t' << best_score_length << '\r' << std::flush;
            size_t cap = a+lookahead_depth;
            for(size_t b=a+1; b<tokens.size() && b<cap; ++b)
            {
                size_t max_match_len = std::min(tokens.size()-b, b-a);
                size_t match_len = 0;
                unsigned hash = 0;
                //int balance = 0;

                while(match_len < max_match_len
                   && tokens[a+match_len] == tokens[b+match_len]
                   && tokens[a+match_len].meta.preproc == false)
                {
                    const token& word = tokens[a+match_len];

                    //balance += word.meta.balance;
                    //if (preserve_params
                    // && (balance < 0 || (word.meta.comma && balance==0))) break;

                    ++match_len;
                    hash = ((hash << (1+0*(match_len&31)))
                          ^ (hash >> (31-0*(match_len&31))))
                          ^ word.meta.hash;
                    //hash = ~hash*0x8088405u + word.meta.hash;

                    donttest[b] = true;

                    //if  (!balance == 0)
                    {
                        std::map<unsigned, length_rec>::iterator i
                            = hash_results.lower_bound(hash);
                        if (i == hash_results.end() || i->first != hash)
                        {
                            length_rec rec;
                            rec.begin_index = (unsigned)a;
                            rec.num_tokens  = (unsigned)match_len;
                            rec.num_occurrences = 1;
                            hash_results.insert(i, std::make_pair(hash,rec));
                            cap = std::max(cap, b+match_len+lookahead_depth);
                        }
                        else if (i->second.begin_index == (unsigned)a)
                        {
                            if (std::equal(
                                tokens.begin()+a, tokens.begin()+a+match_len,
                                tokens.begin() + i->second.begin_index))
                            {
                                long string_len = GetSeq(tokens.begin()+a, match_len, false).size();
                                long n = (i->second.num_occurrences += 1);
                                long define_length = seq_name_length + 9 - long(string_len);
                                long replace_length = long(string_len) - (long(seq_name_length)+1);
                                long score = replace_length * n - define_length;
                                if (score > best_score)
                                {
                                    best_score        = score;
                                    //best_score_length = string_len;
                                    best_hash         = hash;
                                }
                            }
                            cap = std::max(cap, b+match_len+lookahead_depth);
                        }
                    }
                }
            }
        }
        if (best_score > 0)
        {
            const length_rec& rec = hash_results[best_hash];
            if (rec.num_occurrences > 0)
            {
                /* Found a practical saving */
                std::vector<token> sequence
                    (tokens.begin()+rec.begin_index,
                     tokens.begin()+rec.begin_index+rec.num_tokens);
                if(verbose) std::cerr << "#define " << seq_name_buf << " " <<
                    GetSeq(sequence.begin(), sequence.size(), false).GetString()
                        //<< " /* " << rec.num_occurrences
                        //<< " occurrences */"
                          << std::endl;

                /* Replace all occurrences of the sequence with the sequence name */
                std::vector<bool> deletemap(tokens.size(), false);
                for(size_t a=0;//rec.begin_index+rec.num_tokens;
                           a+rec.num_tokens<=tokens.size();
                           ++a)
                {
                    if (std::equal(sequence.begin(),
                                  sequence.end(),
                                  tokens.begin()+a))
                    {
                        tokens[a] = seq_name_buf;
                        for(size_t b=1; b<rec.num_tokens; ++b)
                            deletemap[++a] = true;
                    }
                }
                size_t tgt=0, src=0;
                for(; src < tokens.size(); ++src)
                    if (!deletemap[src])
                        tokens[tgt++].swap(tokens[src]);
                tokens.erase(tokens.begin()+tgt, tokens.end());

                sequence.insert(sequence.begin(),
                    token("#define " + seq_name_buf + " "));
                sequence.push_back( token("\n") );

                tokens.insert(tokens.begin(),
                    sequence.begin(), sequence.end());

                /* Find more repetitions */
                return true;
            }
        }
        return false;
    }

    bool CompressWithParametricMacros(
        std::vector<token>& /*tokens*/,
        const std::string& /*seq_name_buf*/)
    {
        /* TODO: Someone invent an algorithm for this?
         *
         * Find e.g.
         *
         * Cmpeq_i qI x == y ; } q1
         * Cmpge_d qI x >= y ; } q1
         * Cmpge_i qI x >= y ; } q1
         * Cmpgt_d qI x >  y ; } q1
         * Cmpgt_i qI x >  y ; } q1
         * Cmple_d qI x <= y ; } q1
         * Cmple_i qI x <= y ; } q1
         * Cmplt_d qI x <  y ; } q1
         * Cmplt_i qI x <  y ; } q1
         * Cmpne_d qI x != y ; } q1
         * Cmpne_i qI x != y ; } q1
         *
         * Substitute with:
         *
         * #define Zzz(A,B) A qI x B y ; } q1
         * Zzz(Cmpeq_i,==)
         * Zzz(Cmpge_d,>=)
         * Zzz(Cmpge_i,>=)
         * Zzz(Cmpgt_d,>)
         * Zzz(Cmpgt_i,>)
         * Zzz(Cmple_d,<=)
         * Zzz(Cmple_i,<=)
         * Zzz(Cmplt_d,<)
         * Zzz(Cmplt_i,<)
         * Zzz(Cmpne_d,!=)
         * Zzz(Cmpne_i,!=)
         *
         * Which can be further turned into (well, theoretically):
         *
         * #define Zzz(A,B) A qI x B y ; } q1
         * #define Zxy(A,B) Zzz(Cmp##A##_d,B) Zzz(Cmp##A##_i,B)
         * Zzz(Cmpeq_i,==)
         * Zxy(ge,>=)
         * Zxy(gt,>)
         * Zxy(le,<=)
         * Zxy(lt,<)
         * Zxy(ne,!=)
         */
        return false;
    }
}

std::string CPPcompressor::Compress(const std::string& input)
{
    FILE* fp = fopen("cpp_compress_disable", "r");
    if(fp) { fclose(fp); return input; }

    reset_parametric_macro_list();
    macro_counter = 0;

    DefineParsingMode defmode;
    std::vector<token> tokens = Tokenize(input, defmode, false, false);

    int tried_retoken_rounds = 0;
    std::string seq_name_buf = GenerateMacroName();

    //bool preserve_parens = false;

    StringSeq result;
    while (true)
    {
        if (CompressWithNonparametricMacros(tokens, seq_name_buf))
        {
            tried_retoken_rounds = 0;
            seq_name_buf = GenerateMacroName();
        }
        else if (CompressWithParametricMacros(tokens, seq_name_buf))
        {
            tried_retoken_rounds = 0;
            seq_name_buf = GenerateMacroName();
        }
        else
        {
            if (tried_retoken_rounds >= 4) break;
            //preserve_parens = true;

            if(verbose) std::cerr << "Retokenizing\n";
            //static int counter=0; ++counter;
            //if(counter>=1) {Debug=true;}
            result = GetSeq(tokens.begin(), tokens.size(), true);
            //if(counter>=1) break;

            DefineParsingMode defmode;
            tokens = Tokenize(result, defmode, tried_retoken_rounds&1, tried_retoken_rounds&2);
            ++tried_retoken_rounds;
        }
    }
    return result.GetString();
}

std::string CPPcompressor::Compress
    (const std::string& input,
     const std::string& m)
{
    if(!m.empty()) macro_prefix_chars = m;
    return Compress(input);
}
