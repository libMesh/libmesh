
template <class T>
inline const T& max (const T& a, const T& b)
{
  return b>a?b:a;
}

template <class T>
inline T& max (T &a, T &b)
{
  return b>a?b:a;
}

template <class T>
inline const T& min (const T& a, const T& b)
{
  return b>a?a:b;
}

template <class T>
inline T& min (T &a, T &b)
{
  return b>a?a:b;
}

