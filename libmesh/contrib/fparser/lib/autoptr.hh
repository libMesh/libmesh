#ifndef FPOptimizerAutoPtrHH
#define FPOptimizerAutoPtrHH

template<typename Ref>
class FPOPT_autoptr
{
public:
    FPOPT_autoptr()                   : p(0)   { }
    FPOPT_autoptr(Ref*        b) : p(b)   { Birth(); }
    FPOPT_autoptr(const FPOPT_autoptr& b) : p(b.p) { Birth(); }

    inline Ref& operator* () const { return *p; }
    inline Ref* operator->() const { return p; }

    FPOPT_autoptr& operator= (Ref*        b) { Set(b); return *this; }
    FPOPT_autoptr& operator= (const FPOPT_autoptr& b) { Set(b.p); return *this; }
#ifdef __GXX_EXPERIMENTAL_CXX0X__
    FPOPT_autoptr(FPOPT_autoptr&& b)      : p(b.p) { b.p = 0; }
    FPOPT_autoptr& operator= (FPOPT_autoptr&& b) { if(p != b.p) { Forget(); p=b.p; b.p=0; }
                                                   return *this; }
#endif

    ~FPOPT_autoptr() { Forget(); }

    void UnsafeSetP(Ref* newp) { p = newp; }
    void swap(FPOPT_autoptr<Ref>& b) { Ref* tmp=p; p=b.p; b.p=tmp; }

private:
    inline static void Have(Ref* p2);
    inline void Forget();
    inline void Birth();
    inline void Set(Ref* p2);
private:
    Ref* p;
};

//
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Forget()
{
    if(!p) return;
    p->RefCount -= 1;
    if(!p->RefCount) delete p;
    //assert(p->RefCount >= 0);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Have(Ref* p2)
{
    if(p2) ++(p2->RefCount);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Birth()
{
    Have(p);
}
template<typename Ref>
inline void FPOPT_autoptr<Ref>::Set(Ref* p2)
{
    Have(p2);
    Forget();
    p = p2;
}

#endif
