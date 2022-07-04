#define BOOST_PYTHON_MAX_ARITY 7
#define BOOST_SIMD_NO_STRICT_ALIASING 1
#include <pythonic/core.hpp>
#include <pythonic/python/core.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pythonic/types/ndarray.hpp>
#include <pythonic/types/float32.hpp>
#include <pythonic/types/float.hpp>
#include <pythonic/types/numpy_texpr.hpp>
#include <pythonic/types/int.hpp>
#include <pythonic/__builtin__/None.hpp>
#include <pythonic/__builtin__/bool_.hpp>
#include <pythonic/__builtin__/getattr.hpp>
#include <pythonic/__builtin__/list.hpp>
#include <pythonic/__builtin__/print.hpp>
#include <pythonic/__builtin__/str.hpp>
#include <pythonic/__builtin__/sum.hpp>
#include <pythonic/__builtin__/tuple.hpp>
#include <pythonic/__builtin__/xrange.hpp>
#include <pythonic/math/sqrt.hpp>
#include <pythonic/numpy/arccos.hpp>
#include <pythonic/numpy/asarray.hpp>
#include <pythonic/numpy/cos.hpp>
#include <pythonic/numpy/exp.hpp>
#include <pythonic/numpy/float_.hpp>
#include <pythonic/numpy/floor.hpp>
#include <pythonic/numpy/int32.hpp>
#include <pythonic/numpy/isnan.hpp>
#include <pythonic/numpy/sin.hpp>
#include <pythonic/numpy/square.hpp>
#include <pythonic/numpy/zeros.hpp>
#include <pythonic/operator_/mod.hpp>
namespace __pythran_putilse30
{
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  struct tcost
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef double __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type1;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef typename pythonic::assignable<long>::type __type5;
      typedef typename pythonic::assignable<double>::type __type6;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type7;
      typedef indexable<__type4> __type8;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>(), std::declval<__type2>())) __type16;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type16>::type::iterator>::value_type>::type __type17;
      typedef indexable<__type17> __type18;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type19;
      typedef decltype(std::declval<__type19>()[std::declval<__type4>()]) __type20;
      typedef decltype(std::declval<__type19>()[std::declval<__type17>()]) __type21;
      typedef decltype((std::declval<__type20>() - std::declval<__type21>())) __type22;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type23;
      typedef decltype(std::declval<__type23>()[std::declval<__type5>()]) __type25;
      typedef typename pythonic::assignable<decltype((std::declval<__type22>() * std::declval<__type25>()))>::type __type26;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type26>())) __type27;
      typedef container<typename std::remove_reference<__type27>::type> __type28;
      typedef container<typename std::remove_reference<__type0>::type> __type29;
      typedef typename __combined<__type7,__type18,__type8,__type28,__type29>::type __type32;
      typedef decltype(std::declval<__type32>()[std::declval<__type4>()]) __type33;
      typedef decltype((std::declval<__type6>() + std::declval<__type33>())) __type34;
      typedef long __type35;
      typedef decltype((std::declval<__type5>() + std::declval<__type35>())) __type36;
      typedef typename __combined<__type34,__type6>::type __type37;
      typedef typename __combined<__type35,__type36,__type5>::type __type38;
      typedef decltype((std::declval<__type37>() / std::declval<__type38>())) __type39;
      typedef __type0 __ptype0;
      typedef __type4 __ptype1;
      typedef typename pythonic::assignable<typename __combined<__type33,__type34,__type35,__type36,__type39,__type5,__type6>::type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0 const & dists, argument_type1 const & slopes, argument_type2&& cost) const
    ;
  }  ;
  struct expandandadd
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type1;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type5;
      typedef typename std::tuple_element<2,typename std::remove_reference<__type1>::type>::type __type7;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type7>())) __type8;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type10;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type>())) __type11;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type11>::type>::type __type12;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type12>())) __type13;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type13>::type::iterator>::value_type>::type __type14;
      typedef decltype(std::declval<__type10>()[std::declval<__type14>()]) __type15;
      typedef typename std::tuple_element<4,typename std::remove_reference<__type1>::type>::type __type17;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type17>())) __type18;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type18>::type::iterator>::value_type>::type __type19;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type5>(), std::declval<__type9>(), std::declval<__type15>(), std::declval<__type19>())) __type20;
      typedef decltype(std::declval<__type0>()[std::declval<__type20>()]) __type21;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type9>(), std::declval<__type14>(), std::declval<__type19>())) __type22;
      typedef __type21 __ptype16;
      typedef __type22 __ptype17;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type operator()(argument_type0 const & b, argument_type1 const & ref, argument_type2&& index, argument_type3 const & pindex, argument_type4 const & iens, argument_type5&& a, argument_type6 const & sign) const
    ;
  }  ;
  struct had_rebin_pythran_3672_to_1836
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type1;
      typedef long __type2;
      typedef decltype((std::declval<__type1>() - std::declval<__type2>())) __type3;
      typedef decltype((std::declval<__type3>() * std::declval<__type2>())) __type5;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type6;
      typedef decltype((std::declval<__type6>() - std::declval<__type2>())) __type8;
      typedef decltype((std::declval<__type8>() + std::declval<__type2>())) __type10;
      typedef decltype((std::declval<__type10>() * std::declval<__type2>())) __type12;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type5>(), std::declval<__type12>())) __type13;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type13>::type::iterator>::value_type>::type __type14;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type15;
      typedef typename std::tuple_element<2,typename std::remove_reference<__type15>::type>::type __type16;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type16>())) __type17;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type17>::type::iterator>::value_type>::type __type18;
      typedef decltype((std::declval<__type18>() * std::declval<__type2>())) __type20;
      typedef typename pythonic::assignable<decltype((std::declval<__type20>() + std::declval<__type2>()))>::type __type22;
      typedef typename std::tuple_element<3,typename std::remove_reference<__type15>::type>::type __type24;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type24>())) __type25;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type25>::type::iterator>::value_type>::type __type26;
      typedef typename pythonic::lazy<long>::type __type27;
      typedef decltype((std::declval<__type26>() + std::declval<__type27>())) __type28;
      typedef typename pythonic::assignable<decltype((std::declval<__type28>() * std::declval<__type2>()))>::type __type30;
      typedef decltype((std::declval<__type30>() + std::declval<__type2>())) __type32;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type14>(), std::declval<__type22>(), std::declval<__type32>())) __type33;
      typedef typename pythonic::assignable<decltype(std::declval<__type0>()[std::declval<__type33>()])>::type __type34;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type47;
      typedef typename pythonic::assignable<decltype((std::declval<__type47>() - std::declval<__type47>()))>::type __type48;
      typedef decltype((std::declval<__type48>() + std::declval<__type34>())) __type53;
      typedef typename __combined<__type48,__type53>::type __type54;
      typedef decltype((std::declval<__type54>() + std::declval<__type34>())) __type55;
      typedef typename __combined<__type48,__type53,__type55>::type __type56;
      typedef decltype((std::declval<__type56>() + std::declval<__type34>())) __type57;
      typedef typename __combined<__type48,__type53,__type55,__type57>::type __type58;
      typedef decltype((std::declval<__type58>() + std::declval<__type34>())) __type59;
      typedef typename __combined<__type34,__type48,__type53,__type55,__type57,__type59>::type __type60;
      typedef double __type62;
      typedef decltype((std::declval<__type48>() + std::declval<__type62>())) __type64;
      typedef typename __combined<__type48,__type64>::type __type65;
      typedef decltype((std::declval<__type65>() + std::declval<__type62>())) __type67;
      typedef typename __combined<__type48,__type64,__type67>::type __type68;
      typedef decltype((std::declval<__type68>() + std::declval<__type62>())) __type69;
      typedef typename __combined<__type48,__type64,__type67,__type69>::type __type70;
      typedef decltype((std::declval<__type70>() + std::declval<__type62>())) __type72;
      typedef typename __combined<__type48,__type62,__type64,__type67,__type69,__type72>::type __type73;
      typedef decltype((std::declval<__type60>() / std::declval<__type73>())) __type74;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type75;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type15>::type>::type __type77;
      typedef decltype((std::declval<__type77>() - std::declval<__type77>())) __type80;
      typedef typename pythonic::assignable<decltype((std::declval<__type80>() - std::declval<__type2>()))>::type __type82;
      typedef decltype((std::declval<__type82>() + std::declval<__type2>())) __type84;
      typedef typename __combined<__type2,__type82,__type84>::type __type85;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type75>(), std::declval<__type85>(), std::declval<__type18>(), std::declval<__type26>())) __type86;
      typedef __type74 __ptype52;
      typedef __type86 __ptype53;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type operator()(argument_type0 const & mem, argument_type1&& hadens, argument_type2 const & ens, argument_type3 const & startyear, argument_type4 const & endyear, argument_type5 const & nc_miss_val, argument_type6 const & miss_val) const
    ;
  }  ;
  struct __init__
  {
    typedef void callable;
    typedef void pure;
    struct type
    {
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }  ;
    typename type::result_type operator()() const;
    ;
  }  ;
  struct copystride
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type1;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type5;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type>())) __type6;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type7>())) __type8;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
      typedef decltype(std::declval<__type5>()[std::declval<__type9>()]) __type10;
      typedef typename std::tuple_element<2,typename std::remove_reference<__type1>::type>::type __type12;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type12>())) __type13;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type13>::type::iterator>::value_type>::type __type14;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type10>(), std::declval<__type14>())) __type15;
      typedef decltype(std::declval<__type0>()[std::declval<__type15>()]) __type16;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type17;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type18;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type19;
      typedef decltype(std::declval<__type19>()[std::declval<__type14>()]) __type20;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type17>(), std::declval<__type18>(), std::declval<__type4>(), std::declval<__type9>(), std::declval<__type20>())) __type21;
      typedef __type16 __ptype60;
      typedef __type21 __ptype61;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type operator()(argument_type0&& a, argument_type1 const & b, argument_type2 const & index, argument_type3 const & n, argument_type4 const & m, argument_type5 const & pindex, argument_type6 const & fmiss_val) const
    ;
  }  ;
  struct getindex
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef long __type0;
      typedef __type0 __ptype64;
      typedef __type0 __ptype65;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0 const & alldates, argument_type1 const & somedates, argument_type2&& index) const
    ;
  }  ;
  struct rgb
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef pythonic::types::list<__type0> __type1;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type2;
      typedef pythonic::types::list<__type2> __type3;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type4;
      typedef pythonic::types::list<__type4> __type5;
      typedef typename __combined<__type1,__type3,__type5>::type __type6;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::float_())>::type>::type __type7;
      typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::asarray())>::type>::type>()(std::declval<__type6>(), std::declval<__type7>()))>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0 const & r, argument_type1 const & g, argument_type2 const & b) const
    ;
  }  ;
  struct snhtmov
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    struct type
    {
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type0;
      typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type __type1;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type3;
      typedef typename pythonic::assignable<long>::type __type4;
      typedef __type3 __ptype72;
      typedef __type4 __ptype73;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type operator()(argument_type0 const & t, argument_type1&& tsa, argument_type2 const & snhtparas, argument_type3&& index, argument_type4&& count, argument_type5&& tmean, argument_type6&& tsquare) const
    ;
  }  ;
  struct stationaverage
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 >
    struct type
    {
      typedef typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type __type0;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type0>::type>::type __type1;
      typedef typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type __type2;
      typedef typename std::tuple_element<3,typename std::remove_reference<__type0>::type>::type __type3;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type1>(), std::declval<__type2>(), std::declval<__type3>())) __type4;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::float_())>::type>::type __type5;
      typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type4>(), std::declval<__type5>()))>::type __type6;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type8;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type11;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type11>::type::iterator>::value_type>::type __type12;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>())) __type14;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type14>::type::iterator>::value_type>::type __type15;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type9>(), std::declval<__type12>(), std::declval<__type15>())) __type16;
      typedef indexable<__type16> __type17;
      typedef decltype(std::declval<__type6>()[std::declval<__type16>()]) __type23;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type24;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type25;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type25>())) __type26;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type26>::type::iterator>::value_type>::type __type27;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type27>(), std::declval<__type9>(), std::declval<__type12>(), std::declval<__type15>())) __type28;
      typedef decltype(std::declval<__type24>()[std::declval<__type28>()]) __type29;
      typedef decltype((std::declval<__type23>() + std::declval<__type29>())) __type30;
      typedef container<typename std::remove_reference<__type30>::type> __type31;
      typedef typename __combined<__type6,__type17,__type31>::type __type32;
      typedef decltype(std::declval<__type32>()[std::declval<__type16>()]) __type34;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::int32())>::type>::type __type39;
      typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type4>(), std::declval<__type39>()))>::type __type40;
      typedef decltype(std::declval<__type40>()[std::declval<__type16>()]) __type44;
      typedef long __type45;
      typedef decltype((std::declval<__type44>() + std::declval<__type45>())) __type46;
      typedef container<typename std::remove_reference<__type46>::type> __type47;
      typedef typename __combined<__type40,__type17,__type47>::type __type48;
      typedef decltype(std::declval<__type48>()[std::declval<__type16>()]) __type50;
      typedef decltype((std::declval<__type34>() / std::declval<__type50>())) __type51;
      typedef container<typename std::remove_reference<__type51>::type> __type52;
      typedef double __type53;
      typedef container<typename std::remove_reference<__type53>::type> __type54;
      typedef typename __combined<__type6,__type17,__type52,__type54>::type __type55;
      typedef typename pythonic::assignable<decltype(pythonic::types::make_tuple(std::declval<__type55>(), std::declval<__type48>()))>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 >
    typename type<argument_type0, argument_type1>::result_type operator()(argument_type0 const & currentdata, argument_type1 const & thresh) const
    ;
  }  ;
  struct sdist
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type1;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef decltype(std::declval<__type0>()[std::declval<__type4>()]) __type5;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>(), std::declval<__type2>())) __type8;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
      typedef decltype(std::declval<__type0>()[std::declval<__type9>()]) __type10;
      typedef decltype((std::declval<__type5>() * std::declval<__type10>())) __type11;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type12;
      typedef decltype(std::declval<__type12>()[std::declval<__type4>()]) __type13;
      typedef decltype(std::declval<__type12>()[std::declval<__type9>()]) __type14;
      typedef decltype((std::declval<__type13>() * std::declval<__type14>())) __type15;
      typedef decltype((std::declval<__type11>() + std::declval<__type15>())) __type16;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type17;
      typedef decltype(std::declval<__type17>()[std::declval<__type4>()]) __type18;
      typedef decltype(std::declval<__type17>()[std::declval<__type9>()]) __type19;
      typedef decltype((std::declval<__type18>() * std::declval<__type19>())) __type20;
      typedef decltype((std::declval<__type16>() + std::declval<__type20>())) __type21;
      typedef typename pythonic::assignable<long>::type __type22;
      typedef __type21 __ptype116;
      typedef __type22 __ptype117;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0&& dists, argument_type1 const & x, argument_type2 const & y, argument_type3 const & z) const
    ;
  }  ;
  struct snhtmov2
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    struct type
    {
      typedef double __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type1;
      typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef __type0 __ptype120;
      typedef __type4 __ptype121;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type operator()(argument_type0 const & t, argument_type1&& tsa, argument_type2 const & snhtparas, argument_type3&& index, argument_type4&& count, argument_type5&& tmean, argument_type6&& tsquare) const
    ;
  }  ;
  struct zonaltrends
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 >
    struct type
    {
      typedef typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type __type0;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
      typedef pythonic::types::list<__type1> __type2;
      typedef typename std::tuple_element<3,typename std::remove_reference<__type0>::type>::type __type3;
      typedef pythonic::types::list<__type3> __type4;
      typedef typename __combined<__type2,__type4>::type __type5;
      typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type5>()))>::type __type6;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type8;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>())) __type11;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type11>::type::iterator>::value_type>::type __type12;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type9>(), std::declval<__type12>())) __type13;
      typedef indexable<__type13> __type14;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type19;
      typedef pythonic::types::contiguous_slice __type20;
      typedef typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type __type21;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type21>())) __type22;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type22>::type::iterator>::value_type>::type __type23;
      typedef decltype(std::declval<__type19>()(std::declval<__type9>(), std::declval<__type20>(), std::declval<__type23>(), std::declval<__type12>())) __type24;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::isnan())>::type>::type>()(std::declval<__type24>())) __type25;
      typedef typename pythonic::assignable<decltype((~std::declval<__type25>()))>::type __type26;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type9>(), std::declval<__type26>(), std::declval<__type23>(), std::declval<__type12>())) __type27;
      typedef decltype(std::declval<__type19>()[std::declval<__type27>()]) __type28;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type28>())) __type29;
      typedef container<typename std::remove_reference<__type29>::type> __type30;
      typedef double __type31;
      typedef container<typename std::remove_reference<__type31>::type> __type32;
      typedef typename pythonic::assignable<long>::type __type33;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type26>())) __type34;
      typedef decltype((std::declval<__type33>() + std::declval<__type34>())) __type35;
      typedef typename __combined<__type33,__type34,__type35>::type __type36;
      typedef container<typename std::remove_reference<__type36>::type> __type37;
      typedef typename pythonic::assignable<typename __combined<__type6,__type14,__type32,__type37>::type>::type result_type;
    }
    ;
    template <typename argument_type0 >
    typename type<argument_type0>::result_type operator()(argument_type0 const & gslopes) const
    ;
  }  ;
  struct find_gstatindex
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type>())) __type1;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef decltype(std::declval<__type0>()[std::declval<__type4>()]) __type5;
      typedef double __type6;
      typedef decltype((std::declval<__type5>() + std::declval<__type6>())) __type7;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type9;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type9>::type>::type __type10;
      typedef decltype((std::declval<__type6>() / std::declval<__type10>())) __type11;
      typedef decltype((std::declval<__type7>() / std::declval<__type11>())) __type12;
      typedef typename pythonic::lazy<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::floor())>::type>::type>()(std::declval<__type12>()))>::type __type13;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type14;
      typedef decltype(std::declval<__type14>()[std::declval<__type4>()]) __type15;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type17;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type17>::type>::type __type18;
      typedef decltype((std::declval<__type6>() / std::declval<__type18>())) __type19;
      typedef decltype((std::declval<__type15>() / std::declval<__type19>())) __type20;
      typedef typename pythonic::lazy<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::floor())>::type>::type>()(std::declval<__type20>()))>::type __type21;
      typedef long __type22;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type13>(), std::declval<__type21>(), std::declval<__type22>())) __type23;
      typedef __type23 __ptype176;
      typedef __type23 __ptype177;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4>::result_type operator()(argument_type0 const & glons, argument_type1 const & glats, argument_type2 const & lons, argument_type3 const & lats, argument_type4&& gstatindex) const
    ;
  }  ;
  struct statcore
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename pythonic::assignable<double>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type1;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type2;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type2>::type>::type __type3;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>())) __type4;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type __type5;
      typedef decltype(std::declval<__type1>()[std::declval<__type5>()]) __type6;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type7;
      typedef typename pythonic::assignable<long>::type __type8;
      typedef decltype(std::declval<__type7>()[std::declval<__type8>()]) __type9;
      typedef decltype((std::declval<__type6>() * std::declval<__type9>())) __type10;
      typedef decltype((std::declval<__type0>() + std::declval<__type10>())) __type11;
      typedef typename __combined<__type0,__type11>::type __type12;
      typedef decltype((std::declval<__type12>() + std::declval<__type6>())) __type18;
      typedef typename __combined<__type0,__type11,__type18>::type __type19;
      typedef double __type21;
      typedef decltype((std::declval<__type0>() + std::declval<__type9>())) __type23;
      typedef typename __combined<__type0,__type23>::type __type24;
      typedef decltype((std::declval<__type24>() + std::declval<__type21>())) __type25;
      typedef typename __combined<__type0,__type21,__type23,__type25,__type9>::type __type26;
      typedef decltype((std::declval<__type19>() / std::declval<__type26>())) __type27;
      typedef typename __combined<__type0,__type10,__type11,__type18,__type21,__type23,__type25,__type27,__type6,__type9>::type __type28;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type6>())) __type31;
      typedef decltype((std::declval<__type31>() * std::declval<__type9>())) __type33;
      typedef decltype((std::declval<__type0>() + std::declval<__type33>())) __type34;
      typedef typename __combined<__type0,__type34>::type __type35;
      typedef decltype((std::declval<__type35>() + std::declval<__type31>())) __type38;
      typedef typename __combined<__type0,__type34,__type38>::type __type39;
      typedef decltype((std::declval<__type39>() / std::declval<__type26>())) __type41;
      typedef typename __combined<__type0,__type21,__type23,__type25,__type31,__type33,__type34,__type38,__type41,__type9>::type __type42;
      typedef typename pythonic::assignable<decltype(pythonic::types::make_tuple(std::declval<__type28>(), std::declval<__type42>()))>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0 const & feld, argument_type1 const & arr, argument_type2 const & weights, argument_type3 const & dim) const
    ;
  }  ;
  struct expand
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type1;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef long __type5;
      typedef typename std::tuple_element<2,typename std::remove_reference<__type1>::type>::type __type7;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type7>())) __type8;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type10;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type>())) __type11;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type11>::type>::type __type12;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type12>())) __type13;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type13>::type::iterator>::value_type>::type __type14;
      typedef decltype(std::declval<__type10>()[std::declval<__type14>()]) __type15;
      typedef typename std::tuple_element<4,typename std::remove_reference<__type1>::type>::type __type17;
      typedef decltype((std::declval<__type17>() - std::declval<__type5>())) __type19;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type19>())) __type20;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type20>::type::iterator>::value_type>::type __type21;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type5>(), std::declval<__type9>(), std::declval<__type15>(), std::declval<__type21>())) __type22;
      typedef decltype(std::declval<__type0>()[std::declval<__type22>()]) __type23;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type25;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type21>())) __type26;
      typedef decltype(std::declval<__type25>()[std::declval<__type26>()]) __type27;
      typedef decltype((std::declval<__type21>() + std::declval<__type5>())) __type29;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type29>())) __type30;
      typedef decltype(std::declval<__type25>()[std::declval<__type30>()]) __type31;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type5>(), std::declval<__type9>(), std::declval<__type14>(), std::declval<__type27>(), std::declval<__type31>())) __type32;
      typedef __type23 __ptype184;
      typedef __type32 __ptype185;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0 const & b, argument_type1 const & index, argument_type2 const & pindex, argument_type3&& a) const
    ;
  }  ;
  struct getindex2
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef long __type0;
      typedef __type0 __ptype188;
      typedef __type0 __ptype189;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0 const & alldates, argument_type1 const & somedates, argument_type2&& index) const
    ;
  }  ;
  struct thin
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    struct type
    {
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type0;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type3;
      typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type __type5;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type8;
      typedef typename pythonic::lazy<decltype((std::declval<__type1>() / std::declval<__type8>()))>::type __type9;
      typedef __type3 __ptype204;
      typedef __type3 __ptype205;
      typedef typename pythonic::assignable<typename __combined<__type5,__type9>::type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
    typename type<argument_type0, argument_type1, argument_type2>::result_type operator()(argument_type0 const & t, argument_type1&& index, argument_type2 const & n) const
    ;
  }  ;
  struct belttrends
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 >
    struct type
    {
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type0;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
      typedef pythonic::types::list<__type1> __type2;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type3;
      typedef typename std::tuple_element<1,typename std::remove_reference<__type3>::type>::type __type4;
      typedef pythonic::types::list<__type4> __type5;
      typedef typename __combined<__type2,__type5>::type __type6;
      typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type6>()))>::type __type7;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type10;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type10>::type::iterator>::value_type>::type __type11;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>())) __type14;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type14>::type::iterator>::value_type>::type __type15;
      typedef decltype(pythonic::types::make_tuple(std::declval<__type11>(), std::declval<__type15>())) __type16;
      typedef indexable<__type16> __type17;
      typedef double __type20;
      typedef container<typename std::remove_reference<__type20>::type> __type21;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type22;
      typedef pythonic::types::contiguous_slice __type23;
      typedef typename pythonic::assignable<decltype(std::declval<__type22>()(std::declval<__type23>(), std::declval<__type15>()))>::type __type24;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::isnan())>::type>::type>()(std::declval<__type24>())) __type25;
      typedef typename pythonic::assignable<decltype((~std::declval<__type25>()))>::type __type26;
      typedef decltype(std::declval<__type24>()[std::declval<__type26>()]) __type27;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type27>())) __type28;
      typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type26>()))>::type __type29;
      typedef decltype((std::declval<__type28>() / std::declval<__type29>())) __type30;
      typedef container<typename std::remove_reference<__type30>::type> __type31;
      typedef typename pythonic::assignable<typename __combined<__type7,__type17,__type21,__type31>::type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 >
    typename type<argument_type0, argument_type1>::result_type operator()(argument_type0 const & zslopes, argument_type1 const & belts) const
    ;
  }  ;
  struct rmean
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type1;
      typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type>::type __type2;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
      typedef decltype(std::declval<__type0>()[std::declval<__type4>()]) __type5;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type6;
      typedef indexable<__type4> __type9;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type10;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type11;
      typedef typename pythonic::assignable<decltype((std::declval<__type11>() - std::declval<__type11>()))>::type __type12;
      typedef indexable<__type12> __type13;
      typedef container<typename std::remove_reference<__type4>::type> __type14;
      typedef typename __combined<__type10,__type13,__type14>::type __type16;
      typedef long __type17;
      typedef decltype((std::declval<__type12>() + std::declval<__type17>())) __type18;
      typedef typename __combined<__type12,__type17,__type18>::type __type19;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type19>())) __type20;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type20>::type::iterator>::value_type>::type __type21;
      typedef decltype(std::declval<__type16>()[std::declval<__type21>()]) __type22;
      typedef indexable<__type22> __type23;
      typedef container<typename std::remove_reference<__type5>::type> __type32;
      typedef decltype(std::declval<__type0>()[std::declval<__type22>()]) __type35;
      typedef container<typename std::remove_reference<__type35>::type> __type36;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type40;
      typedef typename pythonic::assignable<decltype((std::declval<__type11>() / std::declval<__type17>()))>::type __type42;
      typedef indexable<__type42> __type43;
      typedef decltype((std::declval<__type11>() / std::declval<__type17>())) __type48;
      typedef decltype((std::declval<__type48>() + std::declval<__type17>())) __type50;
      typedef decltype((std::declval<__type19>() - std::declval<__type48>())) __type54;
      typedef decltype((std::declval<__type54>() - std::declval<__type17>())) __type56;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type50>(), std::declval<__type56>())) __type57;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type57>::type::iterator>::value_type>::type __type58;
      typedef indexable<__type58> __type59;
      typedef container<typename std::remove_reference<__type11>::type> __type75;
      typedef double __type76;
      typedef container<typename std::remove_reference<__type76>::type> __type77;
      typedef decltype((-std::declval<__type11>())) __type79;
      typedef decltype((std::declval<__type79>() / std::declval<__type17>())) __type81;
      typedef decltype((std::declval<__type81>() + std::declval<__type17>())) __type83;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type83>(), std::declval<__type48>())) __type86;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type86>::type::iterator>::value_type>::type __type87;
      typedef decltype((std::declval<__type42>() + std::declval<__type87>())) __type88;
      typedef decltype(std::declval<__type16>()[std::declval<__type88>()]) __type89;
      typedef decltype(std::declval<__type0>()[std::declval<__type89>()]) __type90;
      typedef container<typename std::remove_reference<__type90>::type> __type91;
      typedef typename __combined<__type40,__type43,__type77,__type91>::type __type95;
      typedef decltype((std::declval<__type58>() - std::declval<__type17>())) __type97;
      typedef decltype(std::declval<__type95>()[std::declval<__type97>()]) __type98;
      typedef decltype((std::declval<__type98>() * std::declval<__type11>())) __type99;
      typedef decltype((std::declval<__type58>() + std::declval<__type48>())) __type103;
      typedef decltype(std::declval<__type16>()[std::declval<__type103>()]) __type104;
      typedef decltype(std::declval<__type0>()[std::declval<__type104>()]) __type105;
      typedef decltype((std::declval<__type99>() + std::declval<__type105>())) __type106;
      typedef decltype((std::declval<__type58>() - std::declval<__type48>())) __type110;
      typedef decltype(std::declval<__type16>()[std::declval<__type110>()]) __type111;
      typedef decltype(std::declval<__type0>()[std::declval<__type111>()]) __type112;
      typedef decltype((std::declval<__type106>() - std::declval<__type112>())) __type113;
      typedef decltype((std::declval<__type113>() / std::declval<__type11>())) __type114;
      typedef container<typename std::remove_reference<__type114>::type> __type115;
      typedef decltype((std::declval<__type48>() - std::declval<__type17>())) __type124;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type81>(), std::declval<__type124>())) __type125;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type125>::type::iterator>::value_type>::type __type126;
      typedef decltype((std::declval<__type42>() + std::declval<__type126>())) __type127;
      typedef decltype(std::declval<__type16>()[std::declval<__type127>()]) __type128;
      typedef decltype(std::declval<__type0>()[std::declval<__type128>()]) __type129;
      typedef container<typename std::remove_reference<__type129>::type> __type130;
      typedef typename __combined<__type40,__type43,__type59,__type77,__type91>::type __type136;
      typedef decltype(std::declval<__type136>()[std::declval<__type97>()]) __type139;
      typedef decltype((std::declval<__type139>() * std::declval<__type11>())) __type140;
      typedef decltype((std::declval<__type103>() - std::declval<__type17>())) __type146;
      typedef decltype(std::declval<__type16>()[std::declval<__type146>()]) __type147;
      typedef decltype(std::declval<__type0>()[std::declval<__type147>()]) __type148;
      typedef decltype((std::declval<__type140>() + std::declval<__type148>())) __type149;
      typedef decltype((std::declval<__type149>() - std::declval<__type112>())) __type156;
      typedef decltype((std::declval<__type156>() / std::declval<__type11>())) __type157;
      typedef container<typename std::remove_reference<__type157>::type> __type158;
      typedef decltype(std::declval<__type136>()[std::declval<__type21>()]) __type161;
      typedef container<typename std::remove_reference<__type161>::type> __type162;
      typedef __type5 __ptype216;
      typedef __type4 __ptype217;
      typedef typename pythonic::assignable<typename __combined<__type6,__type23,__type9,__type32,__type36>::type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4>::result_type operator()(argument_type0 const & t, argument_type1&& tret, argument_type2&& tmean, argument_type3&& index, argument_type4 const & runmean) const
    ;
  }  ;
  struct plus
  {
    typedef void callable;
    typedef void pure;
    template <typename argument_type0 , typename argument_type1 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type1;
      typedef typename pythonic::assignable<decltype((std::declval<__type0>() + std::declval<__type1>()))>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 >
    typename type<argument_type0, argument_type1>::result_type operator()(argument_type0 const & a, argument_type1 const & b) const
    ;
  }  ;
  struct tdist
  {
    typedef void callable;
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    struct type
    {
      typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type0;
      typedef double __type1;
      typedef decltype((std::declval<__type0>() * std::declval<__type1>())) __type2;
      typedef decltype((std::declval<__type2>() / std::declval<__type1>())) __type4;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::cos())>::type>::type>()(std::declval<__type4>())) __type5;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type6;
      typedef decltype((std::declval<__type6>() * std::declval<__type1>())) __type8;
      typedef decltype((std::declval<__type8>() / std::declval<__type1>())) __type10;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::cos())>::type>::type>()(std::declval<__type10>())) __type11;
      typedef typename pythonic::assignable<decltype((std::declval<__type5>() * std::declval<__type11>()))>::type __type12;
      typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type13;
      typedef typename std::tuple_element<0,typename std::remove_reference<__type13>::type>::type __type14;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type14>())) __type15;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type15>::type::iterator>::value_type>::type __type16;
      typedef decltype(std::declval<__type12>()[std::declval<__type16>()]) __type17;
      typedef long __type18;
      typedef decltype((std::declval<__type16>() + std::declval<__type18>())) __type19;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type19>(), std::declval<__type14>())) __type22;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type22>::type::iterator>::value_type>::type __type23;
      typedef decltype(std::declval<__type12>()[std::declval<__type23>()]) __type24;
      typedef decltype((std::declval<__type17>() * std::declval<__type24>())) __type25;
      typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::sin())>::type>::type>()(std::declval<__type10>())) __type35;
      typedef typename pythonic::assignable<decltype((std::declval<__type5>() * std::declval<__type35>()))>::type __type36;
      typedef decltype(std::declval<__type36>()[std::declval<__type16>()]) __type37;
      typedef decltype(std::declval<__type36>()[std::declval<__type23>()]) __type38;
      typedef decltype((std::declval<__type37>() * std::declval<__type38>())) __type39;
      typedef decltype((std::declval<__type25>() + std::declval<__type39>())) __type40;
      typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::sin())>::type>::type>()(std::declval<__type4>()))>::type __type45;
      typedef decltype(std::declval<__type45>()[std::declval<__type16>()]) __type46;
      typedef decltype(std::declval<__type45>()[std::declval<__type23>()]) __type47;
      typedef decltype((std::declval<__type46>() * std::declval<__type47>())) __type48;
      typedef decltype((std::declval<__type40>() + std::declval<__type48>())) __type49;
      typedef typename pythonic::assignable<long>::type __type50;
      typedef __type49 __ptype268;
      typedef __type50 __ptype269;
      typedef typename pythonic::assignable<pythonic::types::none_type>::type result_type;
    }
    ;
    template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
    typename type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type operator()(argument_type0&& dists, argument_type1 const & lats, argument_type2 const & lons, argument_type3 const & weight) const
    ;
  }  ;
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename tcost::type<argument_type0, argument_type1, argument_type2>::result_type tcost::operator()(argument_type0 const & dists, argument_type1 const & slopes, argument_type2&& cost) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type0;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef typename pythonic::assignable<long>::type __type6;
    typedef long __type7;
    typedef decltype((std::declval<__type6>() + std::declval<__type7>())) __type8;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type9;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type13;
    typedef decltype(std::declval<__type9>()[std::declval<__type13>()]) __type14;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type13>(), std::declval<__type1>())) __type17;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type17>::type::iterator>::value_type>::type __type18;
    typedef decltype(std::declval<__type9>()[std::declval<__type18>()]) __type19;
    typedef decltype((std::declval<__type14>() - std::declval<__type19>())) __type20;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type21;
    typedef decltype(std::declval<__type21>()[std::declval<__type6>()]) __type23;
    typedef typename pythonic::assignable<double>::type __type24;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type25;
    typedef indexable<__type13> __type27;
    typedef indexable<__type18> __type29;
    typedef typename pythonic::assignable<decltype((std::declval<__type20>() * std::declval<__type23>()))>::type __type30;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type30>())) __type31;
    typedef container<typename std::remove_reference<__type31>::type> __type32;
    typedef double __type33;
    typedef container<typename std::remove_reference<__type33>::type> __type34;
    typedef typename __combined<__type25,__type27,__type29,__type32,__type34>::type __type37;
    typedef decltype(std::declval<__type37>()[std::declval<__type13>()]) __type38;
    typedef decltype((std::declval<__type24>() + std::declval<__type38>())) __type39;
    typedef typename __combined<__type24,__type39>::type __type40;
    typedef typename __combined<__type6,__type7,__type8>::type __type41;
    typedef decltype((std::declval<__type40>() / std::declval<__type41>())) __type42;
    typedef typename __combined<__type6,__type8>::type __type45;
    typedef decltype((std::declval<__type1>() - std::declval<__type13>())) __type48;
    typedef decltype((std::declval<__type45>() + std::declval<__type48>())) __type49;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type>::type l;
    typename pythonic::assignable<typename __combined<__type48,__type49,__type6,__type7,__type8>::type>::type id = 0L;
    typename pythonic::assignable<typename __combined<__type6,__type7,__type8>::type>::type goodstats = 0L;
    typename pythonic::assignable<typename __combined<__type24,__type38,__type39,__type42,__type6,__type7,__type8>::type>::type tcost = 0.0;
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(slopes));
      for (long  l_ = 0L; l_ < __target1; l_ += 1L)
      {
        cost[l_] = 0.0;
      }
    }
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(slopes));
      for ( l = 0L; l < __target1; l += 1L)
      {
        if (slopes[l] == slopes[l])
        {
          goodstats += 1L;
          {
            long  __target2 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(slopes));
            for (long  k = l; k < __target2; k += 1L)
            {
              {
                typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type20>() * std::declval<__type23>()))>::type>::type s;
                if ((slopes[l] == slopes[l]) and (slopes[k] == slopes[k]))
                {
                  s = ((slopes[l] - slopes[k]) * dists[id]);
                  cost[l] += pythonic::numpy::proxy::square{}(s);
                  cost[k] += pythonic::numpy::proxy::square{}(s);
                  if ((l == 24L) and (k == 29L))
                  {
                    pythonic::__builtin__::print(dists[id], slopes[l], slopes[k]);
                  }
                }
              }
              id += 1L;
            }
          }
          if (goodstats > 0L)
          {
            tcost += cost[l];
          }
        }
        else
        {
          id += (std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(slopes)) - l);
        }
      }
      if (l == __target1)
      l -= 1L;
    }
    if (goodstats > 0L)
    {
      tcost /= goodstats;
      {
        long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(slopes));
        for ( l = 0L; l < __target1; l += 1L)
        {
          cost[l] /= goodstats;
        }
        if (l == __target1)
        l -= 1L;
      }
    }
    return tcost;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
  typename expandandadd::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type expandandadd::operator()(argument_type0 const & b, argument_type1 const & ref, argument_type2&& index, argument_type3 const & pindex, argument_type4 const & iens, argument_type5&& a, argument_type6 const & sign) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type0;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type3;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type>())) __type4;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type4>::type>::type __type5;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type5>())) __type6;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type6>::type::iterator>::value_type>::type __type7;
    typedef decltype(std::declval<__type3>()[std::declval<__type7>()]) __type8;
    typedef decltype((std::declval<__type8>() + std::declval<__type7>())) __type9;
    typedef typename pythonic::assignable<decltype((std::declval<__type9>() - std::declval<__type8>()))>::type __type11;
    typedef typename pythonic::assignable<decltype((std::declval<__type9>() - std::declval<__type7>()))>::type __type14;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type15;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type16;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type17;
    typedef typename std::tuple_element<4,typename std::remove_reference<__type0>::type>::type __type19;
    typedef long __type20;
    typedef decltype((std::declval<__type19>() - std::declval<__type20>())) __type21;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type21>())) __type22;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type22>::type::iterator>::value_type>::type __type23;
    typedef decltype((std::declval<__type23>() + std::declval<__type20>())) __type25;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type16>(), std::declval<__type17>(), std::declval<__type25>())) __type26;
    typedef indexable<__type26> __type27;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type28;
    typedef typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type __type30;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type30>())) __type31;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type31>::type::iterator>::value_type>::type __type32;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type19>())) __type35;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type35>::type::iterator>::value_type>::type __type36;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type16>(), std::declval<__type32>(), std::declval<__type7>(), std::declval<__type36>())) __type37;
    typedef indexable<__type37> __type38;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type39;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type16>(), std::declval<__type17>(), std::declval<__type32>(), std::declval<__type8>(), std::declval<__type36>())) __type41;
    typedef decltype(std::declval<__type39>()[std::declval<__type41>()]) __type42;
    typedef container<typename std::remove_reference<__type42>::type> __type43;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename __combined<__type28,__type38,__type43>::type>())) __type44;
    typedef typename std::tuple_element<3,typename std::remove_reference<__type44>::type>::type __type45;
    typedef container<typename std::remove_reference<__type45>::type> __type46;
    typedef typename __combined<__type15,__type27,__type46>::type __type47;
    typedef decltype(std::declval<__type47>()[std::declval<__type26>()]) __type51;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type51>())) __type52;
    typedef container<typename std::remove_reference<__type20>::type> __type62;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type52>::type::iterator>::value_type>::type __type65;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type16>(), std::declval<__type32>(), std::declval<__type7>(), std::declval<__type65>())) __type66;
    typedef indexable<__type66> __type67;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type69;
    typedef typename __combined<__type11,__type14>::type __type70;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type16>(), std::declval<__type32>(), std::declval<__type70>(), std::declval<__type65>())) __type71;
    typedef decltype(std::declval<__type69>()[std::declval<__type71>()]) __type72;
    typedef container<typename std::remove_reference<__type72>::type> __type73;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename __combined<__type28,__type38,__type67,__type43,__type73>::type>())) __type75;
    typedef typename std::tuple_element<3,typename std::remove_reference<__type75>::type>::type __type76;
    typedef container<typename std::remove_reference<__type76>::type> __type77;
    typedef typename __combined<__type15,__type27,__type62,__type77>::type __type79;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type16>(), std::declval<__type17>(), std::declval<__type23>())) __type80;
    typedef decltype(std::declval<__type79>()[std::declval<__type80>()]) __type81;
    typedef decltype(std::declval<__type79>()[std::declval<__type26>()]) __type86;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type81>(), std::declval<__type86>())) __type87;
    ;
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b));
      for (long  l = 0L; l < __target1; l += 1L)
      {
        {
          long  __target2 = std::get<2>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b));
          for (long  k = 0L; k < __target2; k += 1L)
          {
            {
              long  __target3 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(pindex));
              for (long  m = 0L; m < __target3; m += 1L)
              {
                {
                  typename pythonic::assignable<typename __combined<__type11,__type14>::type>::type pm;
                  if (sign == 0.0)
                  {
                    {
                      long  __target4 = std::get<4>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b));
                      for (long  j_ = 0L; j_ < __target4; j_ += 1L)
                      {
                        a[pythonic::types::make_tuple(l, k, m, j_)] = b[pythonic::types::make_tuple(l, iens, k, pindex[m], j_)];
                      }
                    }
                  }
                  else
                  {
                    if (std::get<2>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(ref)) != std::get<2>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(a)))
                    {
                      pm = ((pindex[m] + m) - m);
                    }
                    else
                    {
                      pm = ((pindex[m] + m) - pindex[m]);
                    }
                    {
                      long  __target4 = (std::get<4>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b)) - 1L);
                      for (long  j = 0L; j < __target4; j += 1L)
                      {
                        if (index[pythonic::types::make_tuple(l, iens, (j + 1L))] == 0L)
                        {
                          if (j == 0L)
                          {
                            index[pythonic::types::make_tuple(l, iens, (j + 1L))] = std::get<3>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(a));
                            {
                              long  __target5 = index[pythonic::types::make_tuple(l, iens, (j + 1L))];
                              for (long  i__ = 0L; i__ < __target5; i__ += 1L)
                              {
                                a[pythonic::types::make_tuple(l, k, m, i__)] = ref[pythonic::types::make_tuple(l, k, pm, i__)];
                              }
                            }
                            index[pythonic::types::make_tuple(l, iens, (j + 1L))] = 0L;
                          }
                          else
                          {
                            index[pythonic::types::make_tuple(l, iens, (j + 1L))] = std::get<3>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(a));
                            {
                              long  __target5 = index[pythonic::types::make_tuple(l, iens, (j + 1L))];
                              for (long  i_ = index[pythonic::types::make_tuple(l, iens, j)]; i_ < __target5; i_ += 1L)
                              {
                                a[pythonic::types::make_tuple(l, k, m, i_)] = ref[pythonic::types::make_tuple(l, k, pm, i_)];
                              }
                            }
                            index[pythonic::types::make_tuple(l, iens, (j + 1L))] = 0L;
                          }
                          break;
                        }
                        else
                        {
                          {
                            long  __target5 = index[pythonic::types::make_tuple(l, iens, (j + 1L))];
                            for (long  i = index[pythonic::types::make_tuple(l, iens, j)]; i < __target5; i += 1L)
                            {
                              if (ref[pythonic::types::make_tuple(l, k, pm, i)] == ref[pythonic::types::make_tuple(l, k, pm, i)])
                              {
                                a[pythonic::types::make_tuple(l, k, m, i)] = (ref[pythonic::types::make_tuple(l, k, pm, i)] + (sign * b[pythonic::types::make_tuple(l, iens, k, pindex[m], j)]));
                              }
                              else
                              {
                                a[pythonic::types::make_tuple(l, k, m, i)] = -1e+30;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
  typename had_rebin_pythran_3672_to_1836::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type had_rebin_pythran_3672_to_1836::operator()(argument_type0 const & mem, argument_type1&& hadens, argument_type2 const & ens, argument_type3 const & startyear, argument_type4 const & endyear, argument_type5 const & nc_miss_val, argument_type6 const & miss_val) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type0;
    typedef long __type1;
    typedef decltype((std::declval<__type0>() - std::declval<__type1>())) __type2;
    typedef decltype((std::declval<__type2>() * std::declval<__type1>())) __type4;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type5;
    typedef decltype((std::declval<__type5>() - std::declval<__type1>())) __type7;
    typedef decltype((std::declval<__type7>() + std::declval<__type1>())) __type9;
    typedef decltype((std::declval<__type9>() * std::declval<__type1>())) __type11;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>(), std::declval<__type11>())) __type12;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type13;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type13>::type>::type __type14;
    typedef decltype((std::declval<__type14>() - std::declval<__type14>())) __type17;
    typedef typename pythonic::assignable<decltype((std::declval<__type17>() - std::declval<__type1>()))>::type __type19;
    typedef decltype((std::declval<__type19>() + std::declval<__type1>())) __type21;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type22;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type12>::type::iterator>::value_type>::type __type23;
    typedef typename std::tuple_element<2,typename std::remove_reference<__type13>::type>::type __type25;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type25>())) __type26;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type26>::type::iterator>::value_type>::type __type27;
    typedef decltype((std::declval<__type27>() * std::declval<__type1>())) __type29;
    typedef typename pythonic::assignable<decltype((std::declval<__type29>() + std::declval<__type1>()))>::type __type31;
    typedef typename std::tuple_element<3,typename std::remove_reference<__type13>::type>::type __type33;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type33>())) __type34;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type34>::type::iterator>::value_type>::type __type35;
    typedef typename pythonic::lazy<long>::type __type36;
    typedef decltype((std::declval<__type35>() + std::declval<__type36>())) __type37;
    typedef typename pythonic::assignable<decltype((std::declval<__type37>() * std::declval<__type1>()))>::type __type39;
    typedef decltype((std::declval<__type39>() + std::declval<__type1>())) __type41;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type23>(), std::declval<__type31>(), std::declval<__type41>())) __type42;
    typedef typename pythonic::assignable<decltype(std::declval<__type22>()[std::declval<__type42>()])>::type __type43;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type56;
    typedef typename pythonic::assignable<decltype((std::declval<__type56>() - std::declval<__type56>()))>::type __type57;
    typedef decltype((std::declval<__type57>() + std::declval<__type43>())) __type62;
    typedef typename __combined<__type57,__type62>::type __type63;
    typedef decltype((std::declval<__type63>() + std::declval<__type43>())) __type64;
    typedef typename __combined<__type57,__type62,__type64>::type __type65;
    typedef decltype((std::declval<__type65>() + std::declval<__type43>())) __type66;
    typedef typename __combined<__type57,__type62,__type64,__type66>::type __type67;
    typedef decltype((std::declval<__type67>() + std::declval<__type43>())) __type68;
    typedef double __type70;
    typedef decltype((std::declval<__type57>() + std::declval<__type70>())) __type72;
    typedef typename __combined<__type57,__type72>::type __type73;
    typedef decltype((std::declval<__type73>() + std::declval<__type70>())) __type75;
    typedef typename __combined<__type57,__type72,__type75>::type __type76;
    typedef decltype((std::declval<__type76>() + std::declval<__type70>())) __type77;
    typedef typename __combined<__type57,__type72,__type75,__type77>::type __type78;
    typedef decltype((std::declval<__type78>() + std::declval<__type70>())) __type80;
    typename pythonic::assignable<typename __combined<__type1,__type19,__type21>::type>::type index = ((std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(hadens)) - std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(hadens))) - 1L);
    {
      long  __target1 = (((endyear - 1850L) + 1L) * 12L);
      for (long  itime = ((startyear - 1850L) * 12L); itime < __target1; itime += 1L)
      {
        index += 1L;
        {
          long  __target2 = std::get<2>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(hadens));
          for (long  ilat = 0L; ilat < __target2; ilat += 1L)
          {
            {
              long  __target3 = std::get<3>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(hadens));
              for (long  ilon = 0L; ilon < __target3; ilon += 1L)
              {
                typename pythonic::lazy<typename pythonic::lazy<long>::type>::type ishift;
                {
                  if (ilon < 18L)
                  {
                    ishift = 18L;
                  }
                  else
                  {
                    ishift = -18L;
                  }
                  typename pythonic::assignable<typename __combined<__type43,__type57,__type62,__type64,__type66,__type68>::type>::type sum = (nc_miss_val - nc_miss_val);
                  typename pythonic::assignable<typename __combined<__type57,__type70,__type72,__type75,__type77,__type80>::type>::type n = (nc_miss_val - nc_miss_val);
                  typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type37>() * std::declval<__type1>()))>::type>::type jj = ((ilon + ishift) * 2L);
                  ;
                  typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type29>() + std::declval<__type1>()))>::type>::type ii_ = ((ilat * 2L) + 0L);
                  ;
                  typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<__type22>()[std::declval<__type42>()])>::type>::type x__ = mem[pythonic::types::make_tuple(itime, ii_, (jj + 0L))];
                  if (x__ != nc_miss_val)
                  {
                    sum += x__;
                    n += 1.0;
                  }
                  ;
                  typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<__type22>()[std::declval<__type42>()])>::type>::type x_ = mem[pythonic::types::make_tuple(itime, ii_, (jj + 1L))];
                  if (x_ != nc_miss_val)
                  {
                    sum += x_;
                    n += 1.0;
                  }
                  ;
                  typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type29>() + std::declval<__type1>()))>::type>::type ii = ((ilat * 2L) + 1L);
                  ;
                  typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<__type22>()[std::declval<__type42>()])>::type>::type x = mem[pythonic::types::make_tuple(itime, ii, (jj + 0L))];
                  if (x != nc_miss_val)
                  {
                    sum += x;
                    n += 1.0;
                  }
                  ;
                  typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<__type22>()[std::declval<__type42>()])>::type>::type x___ = mem[pythonic::types::make_tuple(itime, ii, (jj + 1L))];
                  if (x___ != nc_miss_val)
                  {
                    sum += x___;
                    n += 1.0;
                  }
                  if (n > 0L)
                  {
                    hadens[pythonic::types::make_tuple(ens, index, ilat, ilon)] = (sum / n);
                  }
                  else
                  {
                    hadens[pythonic::types::make_tuple(ens, index, ilat, ilon)] = miss_val;
                  }
                }
              }
            }
          }
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  typename __init__::type::result_type __init__::operator()() const
  {
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
  typename copystride::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type copystride::operator()(argument_type0&& a, argument_type1 const & b, argument_type2 const & index, argument_type3 const & n, argument_type4 const & m, argument_type5 const & pindex, argument_type6 const & fmiss_val) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type0;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type __type4;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>())) __type5;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type>())) __type6;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type7>())) __type8;
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b));
      for (long  l = 0L; l < __target1; l += 1L)
      {
        {
          long  __target2 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(pindex));
          for (long  k = 0L; k < __target2; k += 1L)
          {
            {
              long  __target3 = std::get<2>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b));
              for (long  j = 0L; j < __target3; j += 1L)
              {
                if (b[pythonic::types::make_tuple(l, pindex[k], j)] != fmiss_val)
                {
                  a[pythonic::types::make_tuple(n, m, l, k, index[j])] = b[pythonic::types::make_tuple(l, pindex[k], j)];
                }
              }
            }
          }
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename getindex::type<argument_type0, argument_type1, argument_type2>::result_type getindex::operator()(argument_type0 const & alldates, argument_type1 const & somedates, argument_type2&& index) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type0;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type1;
    typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type>::type __type2;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
    typedef typename pythonic::lazy<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type>::type __type4;
    typedef typename pythonic::lazy<long>::type __type5;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type6;
    typedef long __type7;
    typedef indexable<__type7> __type8;
    typedef container<typename std::remove_reference<__type7>::type> __type10;
    typedef typename __combined<__type6,__type8,__type10>::type __type11;
    typedef decltype(std::declval<__type11>()[std::declval<__type5>()]) __type12;
    typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type __type13;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type12>(), std::declval<__type13>())) __type14;
    typename pythonic::lazy<typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type>::type>::type n = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(somedates));
    typename pythonic::lazy<typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type>::type m = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(alldates));
    std::get<0>(index) = 0L;
    typename pythonic::lazy<typename __combined<__type4,__type5>::type>::type jsave = 0L;
    {
      long  __target1 = n;
      for (long  j = 0L; j < __target1; j += 1L)
      {
        {
          long  __target2 = m;
          for (long  k = index[jsave]; k < __target2; k += 1L)
          {
            if (alldates[k] == somedates[j])
            {
              index[j] = k;
              jsave = j;
              break;
            }
          }
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename rgb::type<argument_type0, argument_type1, argument_type2>::result_type rgb::operator()(argument_type0 const & r, argument_type1 const & g, argument_type2 const & b) const
  {
    return pythonic::numpy::proxy::asarray{}(typename pythonic::assignable<typename __combined<pythonic::types::list<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>,pythonic::types::list<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>,pythonic::types::list<typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type>>::type>::type({ r, g, b }), pythonic::numpy::proxy::float_{});
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
  typename snhtmov::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type snhtmov::operator()(argument_type0 const & t, argument_type1&& tsa, argument_type2 const & snhtparas, argument_type3&& index, argument_type4&& count, argument_type5&& tmean, argument_type6&& tsquare) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type0;
    typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type __type1;
    typedef long __type2;
    typedef decltype((std::declval<__type1>() / std::declval<__type2>())) __type3;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type4;
    typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type4>::type>::type>::type __type5;
    typedef typename pythonic::assignable<typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type>::type __type6;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>(), std::declval<__type5>(), std::declval<__type6>())) __type7;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type6>::type>::type __type8;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type7>::type::iterator>::value_type>::type __type9;
    typedef indexable<__type9> __type10;
    typedef double __type12;
    typedef container<typename std::remove_reference<__type12>::type> __type13;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type14;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type5>())) __type15;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type15>::type::iterator>::value_type>::type __type16;
    typedef indexable<__type16> __type17;
    typedef typename pythonic::assignable<long>::type __type18;
    typedef decltype((std::declval<__type18>() + std::declval<__type2>())) __type20;
    typedef typename __combined<__type18,__type2,__type20>::type __type21;
    typedef container<typename std::remove_reference<__type21>::type> __type22;
    typedef typename __combined<__type14,__type17,__type22>::type __type24;
    typedef decltype(std::declval<__type24>()[std::declval<__type9>()]) __type25;
    typedef decltype((std::declval<__type9>() - std::declval<__type3>())) __type29;
    typedef decltype(std::declval<__type24>()[std::declval<__type29>()]) __type30;
    typedef decltype((std::declval<__type25>() - std::declval<__type30>())) __type31;
    typedef container<typename std::remove_reference<__type31>::type> __type32;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type35;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type36;
    typedef indexable<__type18> __type37;
    typedef container<typename std::remove_reference<__type16>::type> __type38;
    typedef typename __combined<__type36,__type37,__type38>::type __type39;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type30>(), std::declval<__type25>())) __type47;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type47>::type::iterator>::value_type>::type __type48;
    typedef decltype(std::declval<__type39>()[std::declval<__type48>()]) __type49;
    typedef typename pythonic::assignable<decltype(std::declval<__type35>()[std::declval<__type49>()])>::type __type50;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type50>())) __type51;
    typedef container<typename std::remove_reference<__type51>::type> __type52;
    typedef typename __combined<__type8,__type10,__type32,__type52>::type __type53;
    typedef decltype(std::declval<__type53>()[std::declval<__type9>()]) __type54;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type55;
    typedef container<typename std::remove_reference<__type50>::type> __type58;
    typedef typename __combined<__type55,__type10,__type32,__type58>::type __type72;
    typedef decltype(std::declval<__type72>()[std::declval<__type9>()]) __type73;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type73>())) __type74;
    typedef typename pythonic::assignable<decltype((std::declval<__type3>() / std::declval<__type6>()))>::type __type80;
    typedef decltype((std::declval<__type80>() * std::declval<__type6>())) __type81;
    typedef decltype((std::declval<__type9>() + std::declval<__type81>())) __type82;
    typedef decltype(std::declval<__type72>()[std::declval<__type82>()]) __type83;
    typedef decltype((std::declval<__type73>() + std::declval<__type83>())) __type84;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type6>())) __type86;
    typename pythonic::assignable<typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type>::type n = std::get<0>(snhtparas);
    typename pythonic::assignable<typename pythonic::assignable<typename std::tuple_element<1,typename std::remove_reference<__type0>::type>::type>::type>::type max_miss = std::get<1>(snhtparas);
    typename pythonic::assignable<typename pythonic::assignable<typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type>::type>::type ninc = std::get<2>(snhtparas);
    typename pythonic::lazy<typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type4>::type>::type>::type>::type ni = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(t));
    typename pythonic::assignable<typename __combined<__type18,__type2,__type20>::type>::type good = 0L;
    {
      long  __target1 = ni;
      for (long  j = 0L; j < __target1; j += 1L)
      {
        if (t[j] == t[j])
        {
          index.fast(good) = j;
          good += 1L;
        }
        count[j] = good;
      }
    }
    {
      typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type3>() / std::declval<__type6>()))>::type>::type fak;
      typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type7>::type::iterator>::value_type>::type>::type k;
      if (good > (n - (2L * max_miss)))
      {
        ;
        ;
        {
          long  __target1 = ni;
          std::function<bool(long, long)>  __cmp1 = std::less<long>();
          if (ninc < 0L)
          __cmp1 = std::greater<long>();
          for ( k = (n / 2L); __cmp1(k, __target1); k += ninc)
          {
            if ((count[k] - count[(k - (n / 2L))]) > ((n / 2L) - max_miss))
            {
              tmean[k] = 0.0;
              tsquare[k] = 0.0;
              {
                long  __target2 = count[k];
                for (long  i_ = count[(k - (n / 2L))]; i_ < __target2; i_ += 1L)
                {
                  typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<__type35>()[std::declval<__type49>()])>::type>::type x = t[index[i_]];
                  tmean[k] += x;
                  tsquare[k] += pythonic::numpy::proxy::square{}(x);
                }
              }
              tmean[k] /= (count[k] - count[(k - (n / 2L))]);
              tsquare[k] /= (count[k] - count[(k - (n / 2L))]);
            }
          }
          if (k == __target1)
          k -= ninc;
        }
        fak = ((n / 2L) / ninc);
        {
          long  __target1 = ni;
          std::function<bool(long, long)>  __cmp1 = std::less<long>();
          if (ninc < 0L)
          __cmp1 = std::greater<long>();
          for ( k = (n / 2L); __cmp1(k, __target1); k += ninc)
          {
            {
              typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type54>() - std::declval<__type74>()))>::type>::type tdiv;
              typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type84>() / std::declval<__type12>()))>::type>::type m;
              if (((count[k] - count[(k - (n / 2L))]) > ((n / 2L) - max_miss)) and ((count[(k + (fak * ninc))] - count[k]) > ((n / 2L) - max_miss)))
              {
                m = ((tmean[k] + tmean[(k + (fak * ninc))]) / 2.0);
                tdiv = (tsquare[k] - pythonic::numpy::proxy::square{}(tmean[k]));
                if (tdiv > 0L)
                {
                  tsa[k] = (((n / 2L) * (pythonic::numpy::proxy::square{}((tmean[k] - m)) + pythonic::numpy::proxy::square{}((tmean[(k + (fak * ninc))] - m)))) / pythonic::math::proxy::sqrt{}(tdiv));
                }
                else
                {
                  tsa[k] = 0.0;
                }
                {
                  long  __target2 = ninc;
                  for (long  i = 0L; i < __target2; i += 1L)
                  {
                    tsa[(k + i)] = tsa[k];
                  }
                }
              }
            }
          }
          if (k == __target1)
          k -= ninc;
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 >
  typename stationaverage::type<argument_type0, argument_type1>::result_type stationaverage::operator()(argument_type0 const & currentdata, argument_type1 const & thresh) const
  {
    typedef typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type __type0;
    typedef typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type3;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>())) __type4;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type0>::type>::type __type5;
    typedef typename std::tuple_element<3,typename std::remove_reference<__type0>::type>::type __type7;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type5>(), std::declval<__type1>(), std::declval<__type7>())) __type8;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::int32())>::type>::type __type9;
    typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type8>(), std::declval<__type9>()))>::type __type10;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type5>())) __type12;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type12>::type::iterator>::value_type>::type __type13;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type14;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type7>())) __type16;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type16>::type::iterator>::value_type>::type __type17;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type13>(), std::declval<__type14>(), std::declval<__type17>())) __type18;
    typedef indexable<__type18> __type19;
    typedef decltype(std::declval<__type10>()[std::declval<__type18>()]) __type21;
    typedef long __type22;
    typedef decltype((std::declval<__type21>() + std::declval<__type22>())) __type23;
    typedef container<typename std::remove_reference<__type23>::type> __type24;
    typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::float_())>::type>::type __type29;
    typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type8>(), std::declval<__type29>()))>::type __type30;
    typedef decltype(std::declval<__type30>()[std::declval<__type18>()]) __type38;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type39;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type __type40;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type40>(), std::declval<__type13>(), std::declval<__type14>(), std::declval<__type17>())) __type41;
    typedef decltype(std::declval<__type39>()[std::declval<__type41>()]) __type42;
    typedef decltype((std::declval<__type38>() + std::declval<__type42>())) __type43;
    typedef container<typename std::remove_reference<__type43>::type> __type44;
    typedef typename __combined<__type30,__type19,__type44>::type __type45;
    typedef decltype(std::declval<__type45>()[std::declval<__type18>()]) __type47;
    typedef typename __combined<__type10,__type19,__type24>::type __type48;
    typedef decltype(std::declval<__type48>()[std::declval<__type18>()]) __type50;
    typedef decltype((std::declval<__type47>() / std::declval<__type50>())) __type51;
    typedef container<typename std::remove_reference<__type51>::type> __type52;
    typedef double __type53;
    typedef container<typename std::remove_reference<__type53>::type> __type54;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type12>::type::iterator>::value_type>::type>::type ipar;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type>::type ip;
    typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type16>::type::iterator>::value_type>::type>::type it;
    typename pythonic::assignable<typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type>::type s = pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(currentdata);
    typename pythonic::assignable<typename __combined<__type30,__type19,__type52,__type54>::type>::type currentdataav = pythonic::numpy::proxy::zeros{}(pythonic::types::make_tuple(std::get<1>(s), std::get<2>(s), std::get<3>(s)), pythonic::numpy::proxy::float_{});
    typename pythonic::assignable<typename __combined<__type10,__type19,__type24>::type>::type avcount = pythonic::numpy::proxy::zeros{}(pythonic::types::make_tuple(std::get<1>(s), std::get<2>(s), std::get<3>(s)), pythonic::numpy::proxy::int32{});
    {
      long  __target1 = std::get<0>(s);
      for (long  istat = 0L; istat < __target1; istat += 1L)
      {
        {
          long  __target2 = std::get<1>(s);
          for ( ipar = 0L; ipar < __target2; ipar += 1L)
          {
            {
              long  __target3 = std::get<2>(s);
              for ( ip = 0L; ip < __target3; ip += 1L)
              {
                {
                  long  __target4 = std::get<3>(s);
                  for ( it = 0L; it < __target4; it += 1L)
                  {
                    if (currentdata[pythonic::types::make_tuple(istat, ipar, ip, it)] == currentdata[pythonic::types::make_tuple(istat, ipar, ip, it)])
                    {
                      currentdataav[pythonic::types::make_tuple(ipar, ip, it)] = (currentdataav[pythonic::types::make_tuple(ipar, ip, it)] + currentdata[pythonic::types::make_tuple(istat, ipar, ip, it)]);
                      avcount[pythonic::types::make_tuple(ipar, ip, it)] = (avcount[pythonic::types::make_tuple(ipar, ip, it)] + 1L);
                    }
                  }
                  if (it == __target4)
                  it -= 1L;
                }
              }
              if (ip == __target3)
              ip -= 1L;
            }
          }
          if (ipar == __target2)
          ipar -= 1L;
        }
      }
    }
    {
      long  __target1 = std::get<1>(s);
      for ( ipar = 0L; ipar < __target1; ipar += 1L)
      {
        {
          long  __target2 = std::get<2>(s);
          for ( ip = 0L; ip < __target2; ip += 1L)
          {
            {
              long  __target3 = std::get<3>(s);
              for ( it = 0L; it < __target3; it += 1L)
              {
                if (avcount[pythonic::types::make_tuple(ipar, ip, it)] >= thresh)
                {
                  currentdataav[pythonic::types::make_tuple(ipar, ip, it)] = (currentdataav[pythonic::types::make_tuple(ipar, ip, it)] / avcount[pythonic::types::make_tuple(ipar, ip, it)]);
                }
                else
                {
                  currentdataav[pythonic::types::make_tuple(ipar, ip, it)] = -1e+30;
                }
              }
              if (it == __target3)
              it -= 1L;
            }
          }
          if (ip == __target2)
          ip -= 1L;
        }
      }
      if (ipar == __target1)
      ipar -= 1L;
    }
    return pythonic::types::make_tuple(currentdataav, avcount);
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  typename sdist::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type sdist::operator()(argument_type0&& dists, argument_type1 const & x, argument_type2 const & y, argument_type3 const & z) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type0;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef long __type3;
    typedef typename pythonic::assignable<long>::type __type4;
    typedef decltype((std::declval<__type4>() + std::declval<__type3>())) __type5;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type6;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type6>(), std::declval<__type1>())) __type9;
    typename pythonic::assignable<typename __combined<__type3,__type4,__type5>::type>::type id = 0L;
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(x));
      for (long  l = 0L; l < __target1; l += 1L)
      {
        {
          long  __target2 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(x));
          for (long  k = l; k < __target2; k += 1L)
          {
            dists.fast(id) = (((x[l] * x[k]) + (y[l] * y[k])) + (z[l] * z[k]));
            id += 1L;
          }
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 , typename argument_type5 , typename argument_type6 >
  typename snhtmov2::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4, argument_type5, argument_type6>::result_type snhtmov2::operator()(argument_type0 const & t, argument_type1&& tsa, argument_type2 const & snhtparas, argument_type3&& index, argument_type4&& count, argument_type5&& tmean, argument_type6&& tsquare) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type5>::type>::type __type0;
    typedef typename pythonic::assignable<long>::type __type1;
    typedef indexable<__type1> __type2;
    typedef long __type3;
    typedef indexable<__type3> __type4;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type6;
    typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type>::type __type7;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type7>())) __type8;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type8>::type::iterator>::value_type>::type __type9;
    typedef indexable<__type9> __type10;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type11;
    typedef decltype(std::declval<__type11>()[std::declval<__type9>()]) __type14;
    typedef container<typename std::remove_reference<__type14>::type> __type15;
    typedef double __type16;
    typedef container<typename std::remove_reference<__type16>::type> __type17;
    typedef typename __combined<__type0,__type10,__type4,__type17>::type __type20;
    typedef decltype((std::declval<__type1>() - std::declval<__type3>())) __type22;
    typedef decltype(std::declval<__type20>()[std::declval<__type22>()]) __type23;
    typedef decltype((std::declval<__type23>() + std::declval<__type14>())) __type25;
    typedef container<typename std::remove_reference<__type25>::type> __type26;
    typedef typename __combined<__type0,__type2,__type4,__type17,__type26>::type __type30;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type31;
    typedef decltype((std::declval<__type1>() + std::declval<__type3>())) __type34;
    typedef typename __combined<__type1,__type3,__type34>::type __type35;
    typedef decltype((std::declval<__type35>() - std::declval<__type3>())) __type37;
    typedef container<typename std::remove_reference<__type37>::type> __type38;
    typedef container<typename std::remove_reference<__type3>::type> __type40;
    typedef typename __combined<__type31,__type10,__type38,__type40>::type __type42;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type43;
    typedef typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type43>::type>::type>::type __type44;
    typedef typename pythonic::assignable<decltype((std::declval<__type44>() / std::declval<__type3>()))>::type __type46;
    typedef typename pythonic::assignable<typename std::tuple_element<1,typename std::remove_reference<__type43>::type>::type>::type __type47;
    typedef decltype((std::declval<__type46>() - std::declval<__type47>())) __type48;
    typedef decltype((std::declval<__type7>() - std::declval<__type48>())) __type50;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type48>(), std::declval<__type50>())) __type51;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type51>::type::iterator>::value_type>::type __type52;
    typedef typename pythonic::assignable<decltype((std::declval<__type52>() + std::declval<__type46>()))>::type __type53;
    typedef typename pythonic::assignable<decltype((std::declval<__type7>() - std::declval<__type3>()))>::type __type55;
    typedef typename __combined<__type53,__type55>::type __type56;
    typedef decltype(std::declval<__type42>()[std::declval<__type56>()]) __type57;
    typedef decltype(std::declval<__type30>()[std::declval<__type57>()]) __type58;
    typedef decltype(std::declval<__type42>()[std::declval<__type52>()]) __type61;
    typedef decltype(std::declval<__type30>()[std::declval<__type61>()]) __type62;
    typedef decltype((std::declval<__type58>() - std::declval<__type62>())) __type63;
    typedef decltype((std::declval<__type57>() - std::declval<__type61>())) __type69;
    typedef typename pythonic::assignable<decltype((std::declval<__type52>() - std::declval<__type46>()))>::type __type76;
    typedef typename __combined<__type1,__type76>::type __type78;
    typedef decltype(std::declval<__type42>()[std::declval<__type78>()]) __type79;
    typedef decltype(std::declval<__type30>()[std::declval<__type79>()]) __type80;
    typedef decltype((std::declval<__type62>() - std::declval<__type80>())) __type81;
    typedef decltype((std::declval<__type61>() - std::declval<__type79>())) __type87;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type6>::type>::type __type88;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type14>())) __type97;
    typedef container<typename std::remove_reference<__type97>::type> __type98;
    typedef typename __combined<__type88,__type10,__type4,__type17>::type __type101;
    typedef decltype(std::declval<__type101>()[std::declval<__type22>()]) __type104;
    typedef decltype((std::declval<__type104>() + std::declval<__type97>())) __type107;
    typedef container<typename std::remove_reference<__type107>::type> __type108;
    typedef typename __combined<__type88,__type2,__type4,__type17,__type98>::type __type109;
    typedef decltype(std::declval<__type109>()[std::declval<__type57>()]) __type113;
    typedef decltype(std::declval<__type109>()[std::declval<__type79>()]) __type118;
    typedef decltype((std::declval<__type113>() - std::declval<__type118>())) __type119;
    typedef decltype((std::declval<__type57>() - std::declval<__type79>())) __type126;
    typedef typename pythonic::assignable<decltype((std::declval<__type119>() / std::declval<__type126>()))>::type __type127;
    typedef decltype((std::declval<__type58>() - std::declval<__type80>())) __type138;
    typedef typename pythonic::assignable<decltype((std::declval<__type138>() / std::declval<__type126>()))>::type __type146;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type146>())) __type147;
    typedef decltype((std::declval<__type127>() - std::declval<__type147>())) __type148;
    typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::math::proxy::sqrt())>::type>::type>()(std::declval<__type148>()))>::type __type149;
    typename pythonic::assignable<typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type43>::type>::type>::type>::type n = std::get<0>(snhtparas);
    typename pythonic::assignable<typename pythonic::assignable<typename std::tuple_element<1,typename std::remove_reference<__type43>::type>::type>::type>::type max_miss = std::get<1>(snhtparas);
    ;
    typename pythonic::assignable<typename pythonic::assignable<typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type>::type>::type ni = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(t));
    typename pythonic::assignable<typename __combined<__type1,__type3,__type34>::type>::type good = 0L;
    {
      long  __target1 = ni;
      for (long  j_ = 0L; j_ < __target1; j_ += 1L)
      {
        tmean[j_] = -1e+30;
        tsquare[j_] = -1e+30;
        tsa[j_] = -1e+30;
      }
    }
    std::get<0>(tmean) = 0.0;
    std::get<0>(tsquare) = 0.0;
    {
      long  __target1 = ni;
      for (long  j = 0L; j < __target1; j += 1L)
      {
        count[j] = 0L;
        if (t[j] == t[j])
        {
          index.fast(good) = j;
          if (good > 0L)
          {
            tmean.fast(good) = (tmean[(good - 1L)] + t[j]);
            tsquare.fast(good) = (tsquare[(good - 1L)] + pythonic::numpy::proxy::square{}(t[j]));
          }
          else
          {
            tmean.fast(good) = t[j];
            tsquare.fast(good) = pythonic::numpy::proxy::square{}(t[j]);
          }
          good += 1L;
        }
        if (good > 0L)
        {
          count[j] = (good - 1L);
        }
      }
    }
    {
      typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type44>() / std::declval<__type3>()))>::type>::type rm;
      if (good > (n - (2L * max_miss)))
      {
        rm = (n / 2L);
        {
          long  __target1 = (ni - (rm - max_miss));
          for (long  k = (rm - max_miss); k < __target1; k += 1L)
          {
            typename pythonic::assignable<typename __combined<__type1,__type76>::type>::type xm = (k - rm);
            if (xm < 0L)
            {
              xm = 0L;
            }
            typename pythonic::assignable<typename __combined<__type53,__type55>::type>::type xp = (k + rm);
            if (xp > (ni - 1L))
            {
              xp = (ni - 1L);
            }
            {
              typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type63>() / std::declval<__type69>()))>::type>::type y;
              typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type81>() / std::declval<__type87>()))>::type>::type x;
              typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type138>() / std::declval<__type126>()))>::type>::type xy;
              typename pythonic::assignable<typename __combined<__type127,__type149>::type>::type sig;
              if (((count[k] - count[xm]) > (rm - max_miss)) and ((count[xp] - count[k]) > (rm - max_miss)))
              {
                x = ((tmean[count[k]] - tmean[count[xm]]) / (count[k] - count[xm]));
                y = ((tmean[count[xp]] - tmean[count[k]]) / (count[xp] - count[k]));
                xy = ((tmean[count[xp]] - tmean[count[xm]]) / (count[xp] - count[xm]));
                sig = ((tsquare[count[xp]] - tsquare[count[xm]]) / (count[xp] - count[xm]));
                if (sig > pythonic::numpy::proxy::square{}(xy))
                {
                  sig = pythonic::math::proxy::sqrt{}((sig - pythonic::numpy::proxy::square{}(xy)));
                  tsa[k] = (((((count[k] - count[xm]) * (x - xy)) * (x - xy)) + (((count[xp] - count[k]) * (y - xy)) * (y - xy))) / sig);
                }
                else
                {
                  tsa[k] = 0.0;
                }
                pythonic::__builtin__::print(xm, k, xp, tsa[k], sig, count[k], (count[xp] - count[k]), (count[k] - count[xm]));
              }
            }
          }
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 >
  typename zonaltrends::type<argument_type0>::result_type zonaltrends::operator()(argument_type0 const & gslopes) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
    typedef typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type __type1;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
    typedef pythonic::types::contiguous_slice __type5;
    typedef typename std::tuple_element<2,typename std::remove_reference<__type1>::type>::type __type6;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type6>())) __type7;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type7>::type::iterator>::value_type>::type __type8;
    typedef typename std::tuple_element<3,typename std::remove_reference<__type1>::type>::type __type9;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type9>())) __type10;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type10>::type::iterator>::value_type>::type __type11;
    typedef decltype(std::declval<__type0>()(std::declval<__type4>(), std::declval<__type5>(), std::declval<__type8>(), std::declval<__type11>())) __type12;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::isnan())>::type>::type>()(std::declval<__type12>())) __type13;
    typedef pythonic::types::list<__type2> __type15;
    typedef pythonic::types::list<__type9> __type17;
    typedef typename __combined<__type15,__type17>::type __type18;
    typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type18>()))>::type __type19;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type11>())) __type20;
    typedef indexable<__type20> __type21;
    typedef typename pythonic::assignable<decltype((~std::declval<__type13>()))>::type __type26;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type4>(), std::declval<__type26>(), std::declval<__type8>(), std::declval<__type11>())) __type27;
    typedef decltype(std::declval<__type0>()[std::declval<__type27>()]) __type28;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type28>())) __type29;
    typedef container<typename std::remove_reference<__type29>::type> __type30;
    typedef double __type31;
    typedef container<typename std::remove_reference<__type31>::type> __type32;
    typedef typename pythonic::assignable<long>::type __type33;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type26>())) __type34;
    typedef decltype((std::declval<__type33>() + std::declval<__type34>())) __type35;
    typedef typename __combined<__type33,__type34,__type35>::type __type36;
    typedef container<typename std::remove_reference<__type36>::type> __type37;
    typename pythonic::assignable<typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type>::type s = pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(gslopes);
    typename pythonic::assignable<typename __combined<__type19,__type21,__type32,__type37>::type>::type zslopes = pythonic::numpy::proxy::zeros{}(typename pythonic::assignable<typename __combined<pythonic::types::list<typename std::tuple_element<0,typename std::remove_reference<typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type>::type>::type>,pythonic::types::list<typename std::tuple_element<3,typename std::remove_reference<typename pythonic::assignable<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type>::type>::type>>::type>::type({ std::get<0>(s), std::get<3>(s) }));
    {
      long  __target1 = std::get<0>(s);
      for (long  k = 0L; k < __target1; k += 1L)
      {
        {
          long  __target2 = std::get<3>(s);
          for (long  ip = 0L; ip < __target2; ip += 1L)
          {
            zslopes[pythonic::types::make_tuple(k, ip)] = 0.0;
            typename pythonic::assignable<typename __combined<__type33,__type34,__type35>::type>::type b = 0L;
            {
              long  __target3 = std::get<2>(s);
              for (long  ipar = 0L; ipar < __target3; ipar += 1L)
              {
                typename pythonic::assignable<typename pythonic::assignable<decltype((~std::declval<__type13>()))>::type>::type mask = (~pythonic::numpy::proxy::isnan{}(gslopes(k,pythonic::types::contiguous_slice(pythonic::__builtin__::None,pythonic::__builtin__::None),ipar,ip)));
                b += pythonic::__builtin__::proxy::sum{}(mask);
                zslopes[pythonic::types::make_tuple(k, ip)] += pythonic::__builtin__::proxy::sum{}(gslopes[pythonic::types::make_tuple(k, mask, ipar, ip)]);
              }
            }
            if (b > 0L)
            {
              zslopes[pythonic::types::make_tuple(k, ip)] /= b;
            }
            else
            {
              zslopes[pythonic::types::make_tuple(k, ip)] = -1e+30;
            }
          }
        }
      }
    }
    return zslopes;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 >
  typename find_gstatindex::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4>::result_type find_gstatindex::operator()(argument_type0 const & glons, argument_type1 const & glats, argument_type2 const & lons, argument_type3 const & lats, argument_type4&& gstatindex) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type0;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type>())) __type1;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type1>::type>::type __type2;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type2>())) __type3;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type3>::type::iterator>::value_type>::type __type4;
    typedef decltype(std::declval<__type0>()[std::declval<__type4>()]) __type5;
    typedef double __type6;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type7;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type7>::type>::type __type8;
    typedef decltype((std::declval<__type6>() / std::declval<__type8>())) __type9;
    typedef decltype((std::declval<__type5>() / std::declval<__type9>())) __type10;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type11;
    typedef decltype(std::declval<__type11>()[std::declval<__type4>()]) __type12;
    typedef decltype((std::declval<__type12>() + std::declval<__type6>())) __type14;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type16;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type16>::type>::type __type17;
    typedef decltype((std::declval<__type6>() / std::declval<__type17>())) __type18;
    typedef decltype((std::declval<__type14>() / std::declval<__type18>())) __type19;
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(lats));
      for (long  l = 0L; l < __target1; l += 1L)
      {
        typename pythonic::lazy<typename pythonic::lazy<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::floor())>::type>::type>()(std::declval<__type10>()))>::type>::type ilon = pythonic::numpy::proxy::floor{}((lons[l] / (360.0 / std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(glons)))));
        typename pythonic::lazy<typename pythonic::lazy<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::floor())>::type>::type>()(std::declval<__type19>()))>::type>::type ilat = pythonic::numpy::proxy::floor{}(((lats[l] + 90.0) / (180.0 / std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(glats)))));
        gstatindex[pythonic::types::make_tuple(ilat, ilon, 0L)] += 1L;
        gstatindex[pythonic::types::make_tuple(ilat, ilon, gstatindex[pythonic::types::make_tuple(ilat, ilon, 0L)])] = l;
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  typename statcore::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type statcore::operator()(argument_type0 const & feld, argument_type1 const & arr, argument_type2 const & weights, argument_type3 const & dim) const
  {
    typedef typename pythonic::assignable<double>::type __type0;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type1;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type2;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type2>::type>::type __type3;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>())) __type4;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type __type5;
    typedef decltype(std::declval<__type1>()[std::declval<__type5>()]) __type6;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type7;
    typedef typename pythonic::assignable<long>::type __type8;
    typedef decltype(std::declval<__type7>()[std::declval<__type8>()]) __type9;
    typedef decltype((std::declval<__type6>() * std::declval<__type9>())) __type10;
    typedef decltype((std::declval<__type0>() + std::declval<__type10>())) __type11;
    typedef typename __combined<__type0,__type11>::type __type12;
    typedef decltype((std::declval<__type12>() + std::declval<__type6>())) __type18;
    typedef typename __combined<__type0,__type11,__type18>::type __type19;
    typedef double __type21;
    typedef decltype((std::declval<__type0>() + std::declval<__type9>())) __type23;
    typedef typename __combined<__type0,__type23>::type __type24;
    typedef decltype((std::declval<__type24>() + std::declval<__type21>())) __type25;
    typedef typename __combined<__type0,__type21,__type23,__type25,__type9>::type __type26;
    typedef decltype((std::declval<__type19>() / std::declval<__type26>())) __type27;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::square())>::type>::type>()(std::declval<__type6>())) __type30;
    typedef decltype((std::declval<__type30>() * std::declval<__type9>())) __type32;
    typedef decltype((std::declval<__type0>() + std::declval<__type32>())) __type33;
    typedef typename __combined<__type0,__type33>::type __type34;
    typedef decltype((std::declval<__type34>() + std::declval<__type30>())) __type37;
    typedef typename __combined<__type0,__type33,__type37>::type __type38;
    typedef decltype((std::declval<__type38>() / std::declval<__type26>())) __type40;
    typedef long __type42;
    typedef decltype((std::declval<__type8>() + std::declval<__type42>())) __type43;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type45;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::NDIM>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type46;
    typedef decltype((std::declval<__type46>() - std::declval<__type42>())) __type48;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type3>::type>::type __type49;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type48>(), std::declval<__type49>(), std::declval<__type42>())) __type51;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type51>::type::iterator>::value_type>::type __type52;
    typedef decltype(std::declval<__type45>()[std::declval<__type52>()]) __type53;
    typedef decltype((std::declval<__type8>() * std::declval<__type53>())) __type54;
    typename pythonic::assignable<typename __combined<__type0,__type10,__type11,__type18,__type21,__type23,__type25,__type27,__type6,__type9>::type>::type s = 0.0;
    typename pythonic::assignable<typename __combined<__type0,__type21,__type23,__type25,__type30,__type32,__type33,__type37,__type40,__type9>::type>::type sq = 0.0;
    typename pythonic::assignable<typename __combined<__type0,__type21,__type23,__type25,__type9>::type>::type count = 0.0;
    {
      typename pythonic::assignable<typename __combined<__type53,__type54,__type8>::type>::type stride;
      typename pythonic::assignable<typename __combined<__type42,__type43,__type8>::type>::type ki;
      typename pythonic::assignable<typename __combined<__type42,__type43,__type8>::type>::type k;
      if (std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(weights)) != 0L)
      {
        stride = 1L;
        {
          long  __target1 = dim;
          for (long  i__ = (pythonic::__builtin__::getattr<pythonic::types::attr::NDIM>(feld) - 1L); i__ > __target1; i__ += -1L)
          {
            stride *= pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(feld)[i__];
          }
        }
        k = 0L;
        ki = 0L;
        {
          long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(arr));
          for (long  i_ = 0L; i_ < __target1; i_ += 1L)
          {
            if (arr[i_] == arr[i_])
            {
              s += (arr[i_] * weights.fast(k));
              sq += (pythonic::numpy::proxy::square{}(arr[i_]) * weights.fast(k));
              count += weights.fast(k);
            }
            ki += 1L;
            if (ki == stride)
            {
              k += 1L;
              ki = 0L;
              if (k == std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(weights)))
              {
                k = 0L;
              }
            }
          }
        }
      }
      else
      {
        {
          long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(arr));
          for (long  i = 0L; i < __target1; i += 1L)
          {
            if (arr[i] == arr[i])
            {
              s += arr[i];
              sq += pythonic::numpy::proxy::square{}(arr[i]);
              count += 1.0;
            }
          }
        }
      }
    }
    if (count > 0L)
    {
      s /= count;
      sq /= count;
    }
    return pythonic::types::make_tuple(s, sq);
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  typename expand::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type expand::operator()(argument_type0 const & b, argument_type1 const & index, argument_type2 const & pindex, argument_type3&& a) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type0;
    typedef typename std::tuple_element<2,typename std::remove_reference<__type0>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type>())) __type3;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type3>::type>::type __type4;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>())) __type5;
    typedef typename std::tuple_element<4,typename std::remove_reference<__type0>::type>::type __type7;
    typedef long __type8;
    typedef decltype((std::declval<__type7>() - std::declval<__type8>())) __type9;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type9>())) __type10;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type12;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type12>())) __type13;
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b));
      for (long  l = 0L; l < __target1; l += 1L)
      {
        {
          long  __target2 = std::get<2>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b));
          for (long  k = 0L; k < __target2; k += 1L)
          {
            {
              long  __target3 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(pindex));
              for (long  m = 0L; m < __target3; m += 1L)
              {
                {
                  long  __target4 = (std::get<4>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(b)) - 1L);
                  for (long  j = 0L; j < __target4; j += 1L)
                  {
                    if (index[pythonic::types::make_tuple(l, (j + 1L))] == 0L)
                    {
                      break;
                    }
                    a[pythonic::types::make_tuple(l, 0L, k, m, index[pythonic::types::make_tuple(l, j)], index[pythonic::types::make_tuple(l, (j + 1L))])] = b[pythonic::types::make_tuple(l, 0L, k, pindex[m], j)];
                  }
                }
              }
            }
          }
        }
      }
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename getindex2::type<argument_type0, argument_type1, argument_type2>::result_type getindex2::operator()(argument_type0 const & alldates, argument_type1 const & somedates, argument_type2&& index) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type0;
    typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef typename pythonic::assignable<typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type>::type __type3;
    typedef typename pythonic::assignable<long>::type __type4;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type5;
    typedef long __type6;
    typedef typename pythonic::assignable<decltype((std::declval<__type5>() + std::declval<__type6>()))>::type __type7;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type8;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type9;
    typedef indexable<__type6> __type11;
    typedef container<typename std::remove_reference<__type6>::type> __type13;
    typedef typename __combined<__type9,__type11,__type13>::type __type14;
    typedef decltype(std::declval<__type14>()[std::declval<__type4>()]) __type15;
    typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type8>::type>::type>::type __type16;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type15>(), std::declval<__type16>())) __type17;
    typename pythonic::lazy<typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type>::type n = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(somedates));
    typename pythonic::lazy<typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type8>::type>::type>::type>::type m = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(alldates));
    std::get<0>(index) = 0L;
    typename pythonic::assignable<typename __combined<__type3,__type4,__type7>::type>::type jsave = 0L;
    {
      long  __target1 = n;
      for (long  j = 0L; j < __target1; j += 1L)
      {
        {
          long  __target2 = m;
          for (long  k = index[jsave]; k < __target2; k += 1L)
          {
            if (alldates[k] == somedates[j])
            {
              index[j] = k;
              jsave = j;
              break;
            }
            if (alldates[k] > somedates[j])
            {
              index[j] = k;
              jsave = (j + 1L);
              break;
            }
          }
        }
      }
    }
    if (index[jsave] == 0L)
    {
      index[jsave] = (m - 1L);
    }
    return pythonic::__builtin__::None;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 >
  typename thin::type<argument_type0, argument_type1, argument_type2>::result_type thin::operator()(argument_type0 const & t, argument_type1&& index, argument_type2 const & n) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type0;
    typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type>::type __type1;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type3;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type4;
    typedef typename pythonic::lazy<decltype((std::declval<__type3>() / std::declval<__type4>()))>::type __type5;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>())) __type6;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>())) __type9;
    typedef typename __combined<__type1,__type5>::type __type10;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type10>())) __type11;
    typename pythonic::lazy<typename __combined<__type1,__type5>::type>::type ni = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(t));
    if (n < 2L)
    {
      {
        long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(t));
        for (long  i_ = 0L; i_ < __target1; i_ += 1L)
        {
          index[i_] = i_;
        }
      }
    }
    else
    {
      ni = (std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(t)) / n);
      {
        long  __target1 = ni;
        for (long  i = 0L; i < __target1; i += 1L)
        {
          {
            long  __target2 = n;
            for (long  j = 0L; j < __target2; j += 1L)
            {
              index[i] = (i * n);
              if (t[((i * n) + j)] == t[((i * n) + j)])
              {
                index[i] = ((i * n) + j);
                break;
              }
            }
          }
        }
      }
    }
    return ni;
  }
  template <typename argument_type0 , typename argument_type1 >
  typename belttrends::type<argument_type0, argument_type1>::result_type belttrends::operator()(argument_type0 const & zslopes, argument_type1 const & belts) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type0;
    typedef pythonic::types::contiguous_slice __type1;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type2;
    typedef typename std::tuple_element<1,typename std::remove_reference<__type2>::type>::type __type3;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type3>())) __type4;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type4>::type::iterator>::value_type>::type __type5;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type6;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type6>::type>::type __type7;
    typedef pythonic::types::list<__type7> __type8;
    typedef pythonic::types::list<__type3> __type11;
    typedef typename __combined<__type11,__type8>::type __type12;
    typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::zeros())>::type>::type>()(std::declval<__type12>()))>::type __type13;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type7>())) __type16;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type16>::type::iterator>::value_type>::type __type17;
    typedef decltype(pythonic::types::make_tuple(std::declval<__type17>(), std::declval<__type5>())) __type18;
    typedef indexable<__type18> __type19;
    typedef double __type22;
    typedef container<typename std::remove_reference<__type22>::type> __type23;
    typedef typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type1>(), std::declval<__type5>()))>::type __type24;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::isnan())>::type>::type>()(std::declval<__type24>())) __type25;
    typedef typename pythonic::assignable<decltype((~std::declval<__type25>()))>::type __type26;
    typedef decltype(std::declval<__type24>()[std::declval<__type26>()]) __type27;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type27>())) __type28;
    typedef typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type26>()))>::type __type29;
    typedef decltype((std::declval<__type28>() / std::declval<__type29>())) __type30;
    typedef container<typename std::remove_reference<__type30>::type> __type31;
    typename pythonic::assignable<typename __combined<__type13,__type19,__type23,__type31>::type>::type beltslopes = pythonic::numpy::proxy::zeros{}(typename pythonic::assignable<typename __combined<pythonic::types::list<typename std::tuple_element<0,typename std::remove_reference<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>()))>::type>::type>,pythonic::types::list<typename std::tuple_element<1,typename std::remove_reference<decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>()))>::type>::type>>::type>::type({ std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(belts)), std::get<1>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(zslopes)) }));
    {
      long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(belts));
      for (long  k = 0L; k < __target1; k += 1L)
      {
        {
          long  __target2 = std::get<1>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(zslopes));
          for (long  ip = 0L; ip < __target2; ip += 1L)
          {
            typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<__type0>()(std::declval<__type1>(), std::declval<__type5>()))>::type>::type hilf = zslopes(pythonic::types::contiguous_slice(belts[pythonic::types::make_tuple(k, 0L)],belts[pythonic::types::make_tuple(k, 1L)]),ip);
            typename pythonic::assignable<typename pythonic::assignable<decltype((~std::declval<__type25>()))>::type>::type mask = (~pythonic::numpy::proxy::isnan{}(hilf));
            typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::sum())>::type>::type>()(std::declval<__type26>()))>::type>::type b = pythonic::__builtin__::proxy::sum{}(mask);
            if (b > 0L)
            {
              beltslopes[pythonic::types::make_tuple(k, ip)] = (pythonic::__builtin__::proxy::sum{}(hilf[mask]) / b);
            }
            else
            {
              beltslopes[pythonic::types::make_tuple(k, ip)] = -1e+30;
            }
          }
        }
      }
    }
    return beltslopes;
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 , typename argument_type4 >
  typename rmean::type<argument_type0, argument_type1, argument_type2, argument_type3, argument_type4>::result_type rmean::operator()(argument_type0 const & t, argument_type1&& tret, argument_type2&& tmean, argument_type3&& index, argument_type4 const & runmean) const
  {
    typedef typename std::remove_cv<typename std::remove_reference<argument_type4>::type>::type __type0;
    typedef long __type1;
    typedef decltype((std::declval<__type0>() / std::declval<__type1>())) __type2;
    typedef decltype((std::declval<__type2>() + std::declval<__type1>())) __type4;
    typedef typename pythonic::assignable<decltype((std::declval<__type0>() - std::declval<__type0>()))>::type __type5;
    typedef decltype((std::declval<__type5>() + std::declval<__type1>())) __type7;
    typedef typename __combined<__type1,__type5,__type7>::type __type8;
    typedef decltype((std::declval<__type8>() - std::declval<__type2>())) __type11;
    typedef decltype((std::declval<__type11>() - std::declval<__type1>())) __type13;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type4>(), std::declval<__type13>())) __type14;
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type>())) __type15;
    typedef typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type15>::type>::type>::type __type16;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type16>())) __type17;
    typedef decltype((-std::declval<__type0>())) __type31;
    typedef decltype((std::declval<__type31>() / std::declval<__type1>())) __type33;
    typedef decltype((std::declval<__type33>() + std::declval<__type1>())) __type35;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type35>(), std::declval<__type2>())) __type38;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type8>())) __type40;
    typedef decltype((std::declval<__type2>() - std::declval<__type1>())) __type47;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type33>(), std::declval<__type47>())) __type48;
    typename pythonic::lazy<typename pythonic::lazy<typename std::tuple_element<0,typename std::remove_reference<__type15>::type>::type>::type>::type ni = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(t));
    typename pythonic::assignable<typename __combined<__type1,__type5,__type7>::type>::type good = (runmean - runmean);
    if (runmean < 2L)
    {
      {
        long  __target1 = ni;
        for (long  i___ = 0L; i___ < __target1; i___ += 1L)
        {
          tret[i___] = t[i___];
        }
      }
    }
    else
    {
      {
        long  __target1 = ni;
        for (long  j = 0L; j < __target1; j += 1L)
        {
          if (t[j] == t[j])
          {
            index[good] = j;
            good += 1L;
          }
          else
          {
            tret[j] = t[j];
          }
        }
      }
      {
        typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type0>() / std::declval<__type1>()))>::type>::type i__;
        if (good > (runmean + 2L))
        {
          i__ = (runmean / 2L);
          {
            typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type0>() / std::declval<__type1>()))>::type>::type i_;
            if ((pythonic::operator_::mod(runmean, 2L)) == 1L)
            {
              tmean[i__] = 0.0;
              {
                long  __target1 = (runmean / 2L);
                for (long  k_ = (((-runmean) / 2L) + 1L); k_ < __target1; k_ += 1L)
                {
                  tmean[i__] += t[index[(i__ + k_)]];
                }
              }
              tmean[i__] /= runmean;
              {
                long  __target1 = ((good - (runmean / 2L)) - 1L);
                for (long  i = ((runmean / 2L) + 1L); i < __target1; i += 1L)
                {
                  tmean[i] = ((((tmean[(i - 1L)] * runmean) + t[index[(i + (runmean / 2L))]]) - t[index[(i - (runmean / 2L))]]) / runmean);
                }
              }
            }
            else
            {
              i_ = (runmean / 2L);
              tmean[i_] = 0.0;
              {
                long  __target1 = ((runmean / 2L) - 1L);
                for (long  k = ((-runmean) / 2L); k < __target1; k += 1L)
                {
                  tmean[i_] += t[index[(i_ + k)]];
                }
              }
              tmean[i_] /= runmean;
              {
                long  __target1 = ((good - (runmean / 2L)) - 1L);
                for (long  i____ = ((runmean / 2L) + 1L); i____ < __target1; i____ += 1L)
                {
                  tmean[i____] = ((((tmean[(i____ - 1L)] * runmean) + t[index[((i____ + (runmean / 2L)) - 1L)]]) - t[index[(i____ - (runmean / 2L))]]) / runmean);
                }
              }
            }
          }
          {
            long  __target1 = good;
            for (long  i_____ = 0L; i_____ < __target1; i_____ += 1L)
            {
              tret[index[i_____]] = tmean[i_____];
            }
          }
        }
        else
        {
          {
            long  __target1 = good;
            for (long  i______ = 0L; i______ < __target1; i______ += 1L)
            {
              tret[index[i______]] = t[index[i______]];
            }
          }
        }
      }
    }
    return tret;
  }
  template <typename argument_type0 , typename argument_type1 >
  typename plus::type<argument_type0, argument_type1>::result_type plus::operator()(argument_type0 const & a, argument_type1 const & b) const
  {
    ;
    return (a + b);
  }
  template <typename argument_type0 , typename argument_type1 , typename argument_type2 , typename argument_type3 >
  typename tdist::type<argument_type0, argument_type1, argument_type2, argument_type3>::result_type tdist::operator()(argument_type0&& dists, argument_type1 const & lats, argument_type2 const & lons, argument_type3 const & weight) const
  {
    typedef decltype(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(std::declval<typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type>())) __type0;
    typedef typename std::tuple_element<0,typename std::remove_reference<__type0>::type>::type __type1;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type1>())) __type2;
    typedef typename pythonic::assignable<long>::type __type3;
    typedef long __type4;
    typedef decltype((std::declval<__type3>() + std::declval<__type4>())) __type5;
    typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type2>::type::iterator>::value_type>::type __type6;
    typedef decltype((std::declval<__type6>() + std::declval<__type4>())) __type8;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::__builtin__::proxy::xrange())>::type>::type>()(std::declval<__type8>(), std::declval<__type1>())) __type11;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type1>::type>::type __type12;
    typedef double __type13;
    typedef decltype((std::declval<__type12>() * std::declval<__type13>())) __type14;
    typedef decltype((std::declval<__type14>() / std::declval<__type13>())) __type16;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::cos())>::type>::type>()(std::declval<__type16>())) __type21;
    typedef typename std::remove_cv<typename std::remove_reference<argument_type2>::type>::type __type22;
    typedef decltype((std::declval<__type22>() * std::declval<__type13>())) __type24;
    typedef decltype((std::declval<__type24>() / std::declval<__type13>())) __type26;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::cos())>::type>::type>()(std::declval<__type26>())) __type27;
    typedef decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::sin())>::type>::type>()(std::declval<__type26>())) __type37;
    {
      typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type21>() * std::declval<__type37>()))>::type>::type y;
      typename pythonic::assignable<typename pythonic::assignable<decltype((std::declval<__type21>() * std::declval<__type27>()))>::type>::type x;
      typename pythonic::assignable<typename pythonic::assignable<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::proxy::sin())>::type>::type>()(std::declval<__type16>()))>::type>::type z;
      typename pythonic::assignable<typename __combined<__type3,__type4,__type5>::type>::type id;
      if (std::get<0>(dists) != std::get<0>(dists))
      {
        x = (pythonic::numpy::proxy::cos{}(((lats * 3.141592653589793) / 180.0)) * pythonic::numpy::proxy::cos{}(((lons * 3.141592653589793) / 180.0)));
        y = (pythonic::numpy::proxy::cos{}(((lats * 3.141592653589793) / 180.0)) * pythonic::numpy::proxy::sin{}(((lons * 3.141592653589793) / 180.0)));
        z = pythonic::numpy::proxy::sin{}(((lats * 3.141592653589793) / 180.0));
        id = 0L;
        {
          long  __target1 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(lats));
          for (long  l = 0L; l < __target1; l += 1L)
          {
            {
              long  __target2 = std::get<0>(pythonic::__builtin__::getattr<pythonic::types::attr::SHAPE>(lats));
              for (long  k = (l + 1L); k < __target2; k += 1L)
              {
                dists.fast(id) = (((x[l] * x[k]) + (y[l] * y[k])) + (z[l] * z[k]));
                id += 1L;
              }
            }
          }
        }
        sdist()(dists, x, y, z);
        dists(pythonic::types::contiguous_slice(pythonic::__builtin__::None,pythonic::__builtin__::None)) = pythonic::numpy::proxy::arccos{}((dists * 0.999999));
        if (weight != 0L)
        {
          dists(pythonic::types::contiguous_slice(pythonic::__builtin__::None,pythonic::__builtin__::None)) = (pythonic::numpy::proxy::exp{}(((((-dists) * 40.0) / 2L) / 3.141592653589793)) * 28.284271247461902);
        }
      }
    }
    return pythonic::__builtin__::None;
  }
}

typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::tdist::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, double>::result_type>::type>::type tdist0(pythonic::types::ndarray<double,1> a0, pythonic::types::ndarray<double,1> a1, pythonic::types::ndarray<double,1> a2, double a3)
{
  return __pythran_putilse30::tdist()(a0, a1, a2, a3);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::stationaverage::type<pythonic::types::ndarray<double,4>, double>::result_type>::type>::type stationaverage0(pythonic::types::ndarray<double,4> a0, double a1)
{
  return __pythran_putilse30::stationaverage()(a0, a1);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::copystride::type<pythonic::types::ndarray<double,5>, pythonic::types::ndarray<double,3>, pythonic::types::ndarray<long,1>, long, long, pythonic::types::ndarray<long,1>, double>::result_type>::type>::type copystride0(pythonic::types::ndarray<double,5> a0, pythonic::types::ndarray<double,3> a1, pythonic::types::ndarray<long,1> a2, long a3, long a4, pythonic::types::ndarray<long,1> a5, double a6)
{
  return __pythran_putilse30::copystride()(a0, a1, a2, a3, a4, a5, a6);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::expandandadd::type<pythonic::types::ndarray<double,5>, pythonic::types::ndarray<double,4>, pythonic::types::ndarray<long,3>, pythonic::types::ndarray<long,1>, long, pythonic::types::ndarray<double,5>, double>::result_type>::type>::type expandandadd0(pythonic::types::ndarray<double,5> a0, pythonic::types::ndarray<double,4> a1, pythonic::types::ndarray<long,3> a2, pythonic::types::ndarray<long,1> a3, long a4, pythonic::types::ndarray<double,5> a5, double a6)
{
  return __pythran_putilse30::expandandadd()(a0, a1, a2, a3, a4, a5, a6);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::plus::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>>::result_type>::type>::type plus0(pythonic::types::ndarray<double,1> a0, pythonic::types::ndarray<double,1> a1)
{
  return __pythran_putilse30::plus()(a0, a1);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::had_rebin_pythran_3672_to_1836::type<pythonic::types::ndarray<float,3>, pythonic::types::ndarray<float,4>, long, long, long, float, float>::result_type>::type>::type had_rebin_pythran_3672_to_18360(pythonic::types::ndarray<float,3> a0, pythonic::types::ndarray<float,4> a1, long a2, long a3, long a4, float a5, float a6)
{
  return __pythran_putilse30::had_rebin_pythran_3672_to_1836()(a0, a1, a2, a3, a4, a5, a6);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::sdist::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>>::result_type>::type>::type sdist0(pythonic::types::ndarray<double,1> a0, pythonic::types::ndarray<double,1> a1, pythonic::types::ndarray<double,1> a2, pythonic::types::ndarray<double,1> a3)
{
  return __pythran_putilse30::sdist()(a0, a1, a2, a3);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::expand::type<pythonic::types::ndarray<double,5>, pythonic::types::ndarray<long,2>, pythonic::types::ndarray<long,1>, pythonic::types::ndarray<double,5>>::result_type>::type>::type expand0(pythonic::types::ndarray<double,5> a0, pythonic::types::ndarray<long,2> a1, pythonic::types::ndarray<long,1> a2, pythonic::types::ndarray<double,5> a3)
{
  return __pythran_putilse30::expand()(a0, a1, a2, a3);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::expand::type<pythonic::types::ndarray<double,5>, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,2>>, pythonic::types::ndarray<long,1>, pythonic::types::ndarray<double,5>>::result_type>::type>::type expand1(pythonic::types::ndarray<double,5> a0, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,2>> a1, pythonic::types::ndarray<long,1> a2, pythonic::types::ndarray<double,5> a3)
{
  return __pythran_putilse30::expand()(a0, a1, a2, a3);
}
typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::statcore::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, long>::result_type>::type>::type statcore0(pythonic::types::ndarray<double,1> a0, pythonic::types::ndarray<double,1> a1, pythonic::types::ndarray<double,1> a2, long a3)
{
  return __pythran_putilse30::statcore()(a0, a1, a2, a3);
}

BOOST_PYTHON_MODULE(putilse30)
{
  #ifdef PYTHONIC_TYPES_NDARRAY_HPP
import_array()
#endif
  #ifdef PYTHONIC_BUILTIN_BASEEXCEPTION_HPP
boost::python::register_exception_translator<pythonic::types::BaseException>(&pythonic::translate_BaseException);
#endif
  #ifdef PYTHONIC_BUILTIN_KEYBOARDINTERRUPT_HPP
boost::python::register_exception_translator<pythonic::types::KeyboardInterrupt>(&pythonic::translate_KeyboardInterrupt);
#endif
  #ifdef PYTHONIC_BUILTIN_SYSTEMEXIT_HPP
boost::python::register_exception_translator<pythonic::types::SystemExit>(&pythonic::translate_SystemExit);
#endif
  #ifdef PYTHONIC_BUILTIN_EXCEPTION_HPP
boost::python::register_exception_translator<pythonic::types::Exception>(&pythonic::translate_Exception);
#endif
  #ifdef PYTHONIC_BUILTIN_STANDARDERROR_HPP
boost::python::register_exception_translator<pythonic::types::StandardError>(&pythonic::translate_StandardError);
#endif
  #ifdef PYTHONIC_BUILTIN_LOOKUPERROR_HPP
boost::python::register_exception_translator<pythonic::types::LookupError>(&pythonic::translate_LookupError);
#endif
  #ifdef PYTHONIC_BUILTIN_INDEXERROR_HPP
boost::python::register_exception_translator<pythonic::types::IndexError>(&pythonic::translate_IndexError);
#endif
  #ifdef PYTHONIC_BUILTIN_KEYERROR_HPP
boost::python::register_exception_translator<pythonic::types::KeyError>(&pythonic::translate_KeyError);
#endif
  #ifdef PYTHONIC_BUILTIN_BUFFERERROR_HPP
boost::python::register_exception_translator<pythonic::types::BufferError>(&pythonic::translate_BufferError);
#endif
  #ifdef PYTHONIC_BUILTIN_ASSERTIONERROR_HPP
boost::python::register_exception_translator<pythonic::types::AssertionError>(&pythonic::translate_AssertionError);
#endif
  #ifdef PYTHONIC_BUILTIN_MEMORYERROR_HPP
boost::python::register_exception_translator<pythonic::types::MemoryError>(&pythonic::translate_MemoryError);
#endif
  #ifdef PYTHONIC_BUILTIN_ENVIRONMENTERROR_HPP
boost::python::register_exception_translator<pythonic::types::EnvironmentError>(&pythonic::translate_EnvironmentError);
#endif
  #ifdef PYTHONIC_BUILTIN_OSERROR_HPP
boost::python::register_exception_translator<pythonic::types::OSError>(&pythonic::translate_OSError);
#endif
  #ifdef PYTHONIC_BUILTIN_TYPEERROR_HPP
boost::python::register_exception_translator<pythonic::types::TypeError>(&pythonic::translate_TypeError);
#endif
  #ifdef PYTHONIC_BUILTIN_SYSTEMERROR_HPP
boost::python::register_exception_translator<pythonic::types::SystemError>(&pythonic::translate_SystemError);
#endif
  #ifdef PYTHONIC_BUILTIN_IMPORTERROR_HPP
boost::python::register_exception_translator<pythonic::types::ImportError>(&pythonic::translate_ImportError);
#endif
  #ifdef PYTHONIC_BUILTIN_RUNTIMEERROR_HPP
boost::python::register_exception_translator<pythonic::types::RuntimeError>(&pythonic::translate_RuntimeError);
#endif
  #ifdef PYTHONIC_BUILTIN_NOTIMPLEMENTEDERROR_HPP
boost::python::register_exception_translator<pythonic::types::NotImplementedError>(&pythonic::translate_NotImplementedError);
#endif
  #ifdef PYTHONIC_BUILTIN_SYNTAXERROR_HPP
boost::python::register_exception_translator<pythonic::types::SyntaxError>(&pythonic::translate_SyntaxError);
#endif
  #ifdef PYTHONIC_BUILTIN_INDENTATIONERROR_HPP
boost::python::register_exception_translator<pythonic::types::IndentationError>(&pythonic::translate_IndentationError);
#endif
  #ifdef PYTHONIC_BUILTIN_VALUEERROR_HPP
boost::python::register_exception_translator<pythonic::types::ValueError>(&pythonic::translate_ValueError);
#endif
  #ifdef PYTHONIC_BUILTIN_EOFERROR_HPP
boost::python::register_exception_translator<pythonic::types::EOFError>(&pythonic::translate_EOFError);
#endif
  #ifdef PYTHONIC_BUILTIN_ATTRIBUTEERROR_HPP
boost::python::register_exception_translator<pythonic::types::AttributeError>(&pythonic::translate_AttributeError);
#endif
  #ifdef PYTHONIC_BUILTIN_WARNING_HPP
boost::python::register_exception_translator<pythonic::types::Warning>(&pythonic::translate_Warning);
#endif
  #ifdef PYTHONIC_BUILTIN_SYNTAXWARNING_HPP
boost::python::register_exception_translator<pythonic::types::SyntaxWarning>(&pythonic::translate_SyntaxWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_BYTESWARNING_HPP
boost::python::register_exception_translator<pythonic::types::BytesWarning>(&pythonic::translate_BytesWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_PENDINGDEPRECATIONWARNING_HPP
boost::python::register_exception_translator<pythonic::types::PendingDeprecationWarning>(&pythonic::translate_PendingDeprecationWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_UNICODEWARNING_HPP
boost::python::register_exception_translator<pythonic::types::UnicodeWarning>(&pythonic::translate_UnicodeWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_USERWARNING_HPP
boost::python::register_exception_translator<pythonic::types::UserWarning>(&pythonic::translate_UserWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_FUTUREWARNING_HPP
boost::python::register_exception_translator<pythonic::types::FutureWarning>(&pythonic::translate_FutureWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_RUNTIMEWARNING_HPP
boost::python::register_exception_translator<pythonic::types::RuntimeWarning>(&pythonic::translate_RuntimeWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_STOPITERATION_HPP
boost::python::register_exception_translator<pythonic::types::StopIteration>(&pythonic::translate_StopIteration);
#endif
  #ifdef PYTHONIC_BUILTIN_IMPORTWARNING_HPP
boost::python::register_exception_translator<pythonic::types::ImportWarning>(&pythonic::translate_ImportWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_DEPRECATIONWARNING_HPP
boost::python::register_exception_translator<pythonic::types::DeprecationWarning>(&pythonic::translate_DeprecationWarning);
#endif
  #ifdef PYTHONIC_BUILTIN_REFERENCEERROR_HPP
boost::python::register_exception_translator<pythonic::types::ReferenceError>(&pythonic::translate_ReferenceError);
#endif
  #ifdef PYTHONIC_BUILTIN_ARITHMETICERROR_HPP
boost::python::register_exception_translator<pythonic::types::ArithmeticError>(&pythonic::translate_ArithmeticError);
#endif
  #ifdef PYTHONIC_BUILTIN_ZERODIVISIONERROR_HPP
boost::python::register_exception_translator<pythonic::types::ZeroDivisionError>(&pythonic::translate_ZeroDivisionError);
#endif
  #ifdef PYTHONIC_BUILTIN_FLOATINGPOINTERROR_HPP
boost::python::register_exception_translator<pythonic::types::FloatingPointError>(&pythonic::translate_FloatingPointError);
#endif
  #ifdef PYTHONIC_BUILTIN_OVERFLOWERROR_HPP
boost::python::register_exception_translator<pythonic::types::OverflowError>(&pythonic::translate_OverflowError);
#endif
  #ifdef PYTHONIC_BUILTIN_UNICODEERROR_HPP
boost::python::register_exception_translator<pythonic::types::UnicodeError>(&pythonic::translate_UnicodeError);
#endif
  #ifdef PYTHONIC_BUILTIN_TABERROR_HPP
boost::python::register_exception_translator<pythonic::types::TabError>(&pythonic::translate_TabError);
#endif
  #ifdef PYTHONIC_BUILTIN_NAMEERROR_HPP
boost::python::register_exception_translator<pythonic::types::NameError>(&pythonic::translate_NameError);
#endif
  #ifdef PYTHONIC_BUILTIN_UNBOUNDLOCALERROR_HPP
boost::python::register_exception_translator<pythonic::types::UnboundLocalError>(&pythonic::translate_UnboundLocalError);
#endif
  #ifdef PYTHONIC_BUILTIN_IOERROR_HPP
boost::python::register_exception_translator<pythonic::types::IOError>(&pythonic::translate_IOError);
#endif
  #ifdef PYTHONIC_BUILTIN_GENERATOREXIT_HPP
boost::python::register_exception_translator<pythonic::types::GeneratorExit>(&pythonic::translate_GeneratorExit);
#endif
  #ifdef _OPENMP
omp_set_max_active_levels(1);
#endif
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::tdist::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, double>::result_type>::type>::type>();
  boost::python::def("tdist", &tdist0);
  pythonic::python_to_pythran<pythonic::types::ndarray<double,3>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,4>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,2>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::stationaverage::type<pythonic::types::ndarray<double,4>, double>::result_type>::type>::type>();
  boost::python::def("stationaverage", &stationaverage0);
  pythonic::python_to_pythran<pythonic::types::ndarray<long,1>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,5>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,4>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,2>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,3>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::copystride::type<pythonic::types::ndarray<double,5>, pythonic::types::ndarray<double,3>, pythonic::types::ndarray<long,1>, long, long, pythonic::types::ndarray<long,1>, double>::result_type>::type>::type>();
  boost::python::def("copystride", &copystride0);
  pythonic::python_to_pythran<pythonic::types::ndarray<long,2>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<long,1>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<long,3>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,5>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,4>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,2>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,3>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::expandandadd::type<pythonic::types::ndarray<double,5>, pythonic::types::ndarray<double,4>, pythonic::types::ndarray<long,3>, pythonic::types::ndarray<long,1>, long, pythonic::types::ndarray<double,5>, double>::result_type>::type>::type>();
  boost::python::def("expandandadd", &expandandadd0);
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::plus::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>>::result_type>::type>::type>();
  boost::python::def("plus", &plus0);
  pythonic::python_to_pythran<pythonic::types::ndarray<float,1>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<float,3>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<float,2>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<float,4>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::had_rebin_pythran_3672_to_1836::type<pythonic::types::ndarray<float,3>, pythonic::types::ndarray<float,4>, long, long, long, float, float>::result_type>::type>::type>();
  boost::python::def("had_rebin_pythran_3672_to_1836", &had_rebin_pythran_3672_to_18360);
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::sdist::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>>::result_type>::type>::type>();
  boost::python::def("sdist", &sdist0);
  pythonic::python_to_pythran<pythonic::types::ndarray<long,2>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<long,1>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,5>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,4>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,2>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,3>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::expand::type<pythonic::types::ndarray<double,5>, pythonic::types::ndarray<long,2>, pythonic::types::ndarray<long,1>, pythonic::types::ndarray<double,5>>::result_type>::type>::type>();
  boost::python::def("expand", &expand0);
  pythonic::python_to_pythran<pythonic::types::ndarray<long,1>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,5>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,4>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,2>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,3>>();
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::python_to_pythran<pythonic::types::numpy_texpr<pythonic::types::ndarray<long,2>>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::expand::type<pythonic::types::ndarray<double,5>, pythonic::types::numpy_texpr<pythonic::types::ndarray<long,2>>, pythonic::types::ndarray<long,1>, pythonic::types::ndarray<double,5>>::result_type>::type>::type>();
  boost::python::def("expand", &expand1);
  pythonic::python_to_pythran<pythonic::types::ndarray<double,1>>();
  pythonic::pythran_to_python<typename std::remove_cv<typename std::remove_reference<typename __pythran_putilse30::statcore::type<pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, pythonic::types::ndarray<double,1>, long>::result_type>::type>::type>();
  boost::python::def("statcore", &statcore0);
  __pythran_putilse30::__init__()();
}