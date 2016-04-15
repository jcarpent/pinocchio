//
// Copyright (c) 2015-2016 CNRS
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#ifndef __se3_force_hpp__
#define __se3_force_hpp__

#include <Eigen/Core>
#include "pinocchio/spatial/spatial.hpp"
#include "pinocchio/assert.hpp"

namespace se3
{
  
  template<class Derived>
  struct traits< ForceBase<Derived> >
  {
    typedef Derived SE3ActionReturnType;
  };

  template<class Derived>
  class ForceBase : public SpatialBase<Derived>
  {
  protected:

    SPATIAL_TYPEDEF_TEMPLATE(Derived);
    
  public:
    
    typedef SpatialBase<Derived> Base;
    using Base::derived;
    
    template<typename MotionDerived>
    inline typename traits<Derived>::DotReturnType dot(const MotionBase<MotionDerived> & m) const
    {
      return derived().dot(m.derived());
    }
    
    template<typename OtherDerived>
    inline bool operator== (const ForceBase<OtherDerived> & other) const
    {
      return derived().isEqual(other.derived());
    }
    
    template<typename OtherDerived>
    inline Derived & operator= (const ForceBase<OtherDerived> & other)
    {
      other.derived().setTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived & operator+= (const ForceBase<OtherDerived> & other)
    {
      other.derived().addTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived & operator-= (const ForceBase<OtherDerived> & other)
    {
      other.derived().subTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived operator+ (const ForceBase<OtherDerived> & other) const
    {
      return derived().add(other.derived());
    }
    
    template<typename OtherDerived>
    inline Derived operator- (const ForceBase<OtherDerived> & other) const
    {
      return derived().sub(other.derived());
    }
    
    inline Derived operator- () const
    {
      return derived().opposite();
    }
    

//    ConstAngular_t angular() const { return derived().angular_impl(); }
//    ConstLinear_t linear() const { return derived().linear_impl(); }
//    Angular_t angular() { return derived().angular_impl(); }
//    Linear_t linear() { return derived().linear_impl(); }
//    
//    template<typename D>
//    void angular(const Eigen::MatrixBase<D> & n) { derived().angular_impl(n); }
//    
//    template<typename D>
//    void linear(const Eigen::MatrixBase<D> & f) { derived().linear_impl(f); }

//    const Vector6 & toVector() const { return derived().toVector_impl(); }
//    Vector6 & toVector() { return derived().toVector_impl(); }
//    operator Vector6 () const { return toVector(); }

//    bool operator== (const Derived & other) const {return derived().isEqual(other);}
//    Derived & operator= (const Derived & other) { return derived().__equl__(other); }
//    Derived & operator+= (const Derived & phi) { return derived().__pequ__(phi); }
//    Derived & operator-= (const Derived & phi) { return derived().__mequ__(phi); }
//    Derived operator+(const Derived & phi) const { return derived().__plus__(phi); }
//    Derived operator*(double a) const    { return derived().__mult__(a); }
//    Derived operator-() const { return derived().__minus__(); }
//    Derived operator-(const Derived & phi) const { return derived().__minus__(phi); }
    

//    Derived se3Action(const SE3 & m) const { return derived().se3Action_impl(m); }
//    Derived se3ActionInverse(const SE3 & m) const { return derived().se3ActionInverse_impl(m); }

    

  }; // class ForceBase


  template<typename T, int U>
  struct traits< ForceTpl<T, U> > : traits < ForceBase< ForceTpl<T,U> > >
  {
    typedef ForceTpl<T,U> Type;
    typedef traits< ForceBase<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
    typedef T DotReturnType;
    typedef T Scalar_t;
    typedef Eigen::Matrix<T,3,1,U> Vector3;
    typedef Eigen::Matrix<T,4,1,U> Vector4;
    typedef Eigen::Matrix<T,6,1,U> Vector6;
    typedef Eigen::Matrix<T,3,3,U> Matrix3;
    typedef Eigen::Matrix<T,4,4,U> Matrix4;
    typedef Eigen::Matrix<T,6,6,U> Matrix6;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type Linear_t;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type Angular_t;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstLinear_t;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstAngular_t;
    typedef Matrix6 ActionMatrix_t;
    typedef Eigen::Quaternion<T,U> Quaternion_t;
    typedef SE3Tpl<T,U> SE3;
    typedef ForceTpl<T,U> Force;
    typedef MotionTpl<T,U> Motion;
    typedef Symmetric3Tpl<T,U> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
  }; // traits ForceTpl

  template<typename _Scalar, int _Options>
  class ForceTpl : public ForceBase< ForceTpl<_Scalar,_Options> >
  {
  public:
    
    typedef class ForceBase<ForceTpl> Base;
    friend class ForceBase< ForceTpl<_Scalar,_Options> >;
    typedef typename traits<ForceTpl>::SE3ActionReturnType SE3ActionReturnType;
    
    SPATIAL_TYPEDEF_TEMPLATE(ForceTpl);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    ForceTpl() : data() {}
    
    using Base::operator=;

    template<typename F3,typename N3>
    inline ForceTpl(const Eigen::MatrixBase<F3> & f,const Eigen::MatrixBase<N3> & n)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(F3,3);
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(N3,3);
      data << f, n;
    }

    template<typename F6>
    explicit ForceTpl(const Eigen::MatrixBase<F6> & f)
    : data(f)
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector6,F6);;
    }

    template<typename OtherScalar, int OtherOptions>
    explicit ForceTpl(const ForceTpl<OtherScalar,OtherOptions> & clone)
    : data(clone.coeffs())
    {}

    static inline ForceTpl Zero() { return ForceTpl(Linear_t::Zero(), Angular_t::Zero()); }
    static inline ForceTpl Random() { return ForceTpl(Linear_t::Random(), Angular_t::Random()); }

    inline ForceTpl & setZero () { data.setZero (); return *this; }
    inline ForceTpl & setRandom () { data.setRandom (); return *this; }

//    const Vector6 & toVector_impl() const { return data; }
//    Vector6 & toVector_impl() { return data; }

    void disp(std::ostream & os) const
    {
      using namespace std;
      os
      << "  f = " << linear().transpose() << endl
      << "  tau = " << angular().transpose() << endl;
    }

    /// af = aXb.act(bf)
    template<typename SE3Scalar, int SE3Options>
    inline ForceTpl SE3ActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    {
      const Vector3 & Rf = M.rotation() * linear();
      return ForceTpl(Rf,M.translation().cross(Rf)+M.rotation()*angular());
    }


    template<typename SE3Scalar, int SE3Options>
    inline ForceTpl SE3InvActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    {
      return ForceTpl(M.rotation().transpose()*linear(),
        M.rotation().transpose()*(angular() - M.translation().cross(linear())) );
    }
    
    inline bool isEqual (const ForceTpl & other) const { return data == other.data; }

    // Arithmetic operators
//    template<typename S2, int O2>
//    ForceTpl & operator= (const ForceTpl<S2,O2> & other)
//    {
//      data = other.toVector();
//      return *this;
//    }

    template<typename F6>
    inline ForceTpl & operator= (const Eigen::MatrixBase<F6> & phi)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(F6,6);
      data = phi;
      return *this;
    }
    
//    template<typename OtherScalar, int OtherOptions>
//    ForceTpl & operator= (const ForceTpl<OtherScalar,OtherOptions> & clone)
//    {
//      data = clone.coeffs();
//      return *this;
//    }
    
    template<typename OtherDerived>
    inline bool isEqual (const ForceBase<OtherDerived> & other) const
    { return other.isEqual(*this); }
    
    template<typename OtherDerived>
    inline void setTo (ForceBase<OtherDerived> &) const
    {
      PINOCCHIO_STATIC_ASSERT(true,YOU_CALLED_A_METHOD_WHICH_CANNOT_BE_APPLIED);
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline void setTo (ForceTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.coeffs() = data;
    }
    
    
    template<typename OtherScalar, int OtherOptions>
    inline void addTo (ForceTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.coeffs() += data;
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline void subTo (ForceTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.coeffs() -= data;
    }
    
    template<typename OtherDerived>
    inline ForceTpl add (const ForceBase<OtherDerived> & other) const
    { return other.derived().add(*this); }
    
    template<typename OtherScalar, int OtherOptions>
    inline ForceTpl<OtherScalar, OtherOptions> add (const ForceTpl<OtherScalar, OtherOptions> & other) const
    {
      typedef ForceTpl<OtherScalar, OtherOptions> ReturnType;
      return ReturnType(data + other.coeffs());
    }
    
    template<typename OtherDerived>
    inline ForceTpl sub (const ForceBase<OtherDerived> & other) const
    { return other.derived().sub(*this); }
    
    template<typename OtherScalar, int OtherOptions>
    inline ForceTpl<OtherScalar, OtherOptions> sub (const ForceTpl<OtherScalar, OtherOptions> & other) const
    {
      typedef ForceTpl<OtherScalar, OtherOptions> ReturnType;
      return ReturnType(data - other.coeffs());
    }
    
    inline ForceTpl opposite () const
    {
      return ForceTpl(-data);
    }
    
//    template<typename SE3Scalar, int SE3Options>
//    SE3ActionReturnType SE3ActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
//    {
//      Vector3 Rf (M.rotation() * linear());
//      return ForceTpl(Rf, M.translation().cross(Rf) + M.rotation()*angular());
//    }
//    
//    template<typename SE3Scalar, int SE3Options>
//    SE3ActionReturnType SE3InvActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
//    {
//      return ForceTpl(M.rotation().transpose()*linear(),
//                      M.rotation().transpose()*(angular() - M.translation().cross(linear())));
//    }

//    ForceTpl & __equl__(const ForceTpl & other) { data = other.data; return *this; }
//    ForceTpl & __pequ__ (const ForceTpl & phi) { data += phi.data; return *this; }
//    ForceTpl & __mequ__ (const ForceTpl & phi) { data -= phi.data; return *this; }
//    ForceTpl __plus__(const ForceTpl & phi) const { return ForceTpl(data + phi.data); }
//    ForceTpl __mult__(const double a) const { return ForceTpl(a*data); }
//    ForceTpl __minus__() const { return ForceTpl(-data); }
//    ForceTpl __minus__(const ForceTpl & phi) const { return ForceTpl(data - phi.data); }


//    ConstAngular_t angular_impl() const { return data.template segment<3> (ANGULAR); }
//    Angular_t angular_impl() { return data.template segment<3> (ANGULAR); }
    Angular_t angular() { return data.template segment<3> (ANGULAR); }
    ConstAngular_t angular() const { return data.template segment<3> (ANGULAR); }
    
    template<typename D>
    inline void angular(const Eigen::MatrixBase<D> & n)
    { EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(D,3); data.template segment<3> (ANGULAR) = n; }
    
//    ConstLinear_t linear_impl() const { return data.template segment<3> (LINEAR);}
//    Linear_t linear_impl() { return data.template segment<3> (LINEAR);}
    Linear_t linear() { return data.template segment<3> (LINEAR);}
    ConstLinear_t linear() const { return data.template segment<3> (LINEAR);}
    template<typename D>
    inline void linear(const Eigen::MatrixBase<D> & f)
    { EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(D,3); data.template segment<3> (LINEAR) = f; }
    
    
    Vector6 & coeffs() { return data; }
    const Vector6 & coeffs() const { return data; }
    operator Vector6 () const { return coeffs(); }

  protected:
    Vector6 data;

  }; // class ForceTpl

  typedef ForceTpl<double,0> Force;

} // namespace se3

#endif // ifndef __se3_force_hpp__

