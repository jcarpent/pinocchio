//
// Copyright (c) 2015-2016 CNRS
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
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

#ifndef __se3_motion_hpp__
#define __se3_motion_hpp__

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "pinocchio/spatial/spatial.hpp"

namespace se3
{
  
  template<class Derived>
  struct traits< MotionBase<Derived> >
  {
    typedef Derived SE3ActionReturnType;
  };

  template<class Derived>
  class MotionBase : public SpatialBase<Derived>
  {
  protected:
    SPATIAL_TYPEDEF_TEMPLATE(Derived);
    
  public:
    
    typedef SpatialBase<Derived> Base;
    using Base::derived;

//    ConstAngular_t angular() const  { return derived().angular_impl(); }
//    ConstLinear_t linear() const  { return derived().linear_impl(); }
//    Angular_t angular()  { return derived().angular_impl(); }
//    Linear_t linear()   { return derived().linear_impl(); }
    
//    template<typename D>
//    void angular(const Eigen::MatrixBase<D> & w) { derived().angular_impl(w); }
//    template<typename D>
//    void linear(const Eigen::MatrixBase<D> & v) { derived().linear_impl(v); }

//    const Vector6 & toVector() const { return derived().toVector_impl(); }
//    Vector6 & toVector() { return derived().toVector_impl(); }
//    operator Vector6 () const { return toVector(); }

    inline ActionMatrix_t toActionMatrix() const { return derived().toActionMatrix_impl(); }
    inline operator Matrix6 () const { return toActionMatrix(); }
    
    template<typename ForceDerived>
    inline typename traits<Derived>::DotReturnType dot(const ForceBase<ForceDerived> & m) const
    {
      return derived().dot(m.derived());
    }
    
    template<typename OtherDerived>
    inline bool operator== (const MotionBase<OtherDerived> & other) const
    {
      return derived().isEqual(other.derived());
    }
    
    template<typename OtherDerived>
    inline Derived & operator= (const MotionBase<OtherDerived> & other)
    {
      other.derived().setTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived & operator+= (const MotionBase<OtherDerived> & other)
    {
      other.derived().addTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived & operator-= (const MotionBase<OtherDerived> & other)
    {
      other.derived().subTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived operator+ (const MotionBase<OtherDerived> & other) const
    {
      return derived().add(other.derived());
    }
    
    template<typename OtherDerived>
    inline Derived add (const MotionBase<OtherDerived> & other) const
    {
//      Derived dest = derived();
//      other.derived().addTo(dest);
//      return dest;
      return derived().add(other.derived());
    }

    
    template<typename OtherDerived>
    inline Derived operator- (const MotionBase<OtherDerived> & other) const
    {
      return derived().sub(other.derived());
    }
    
    template<typename OtherDerived>
    inline Derived sub (const MotionBase<OtherDerived> & other) const
    {
      return derived().sub(other.derived());
    }
    
    inline Derived operator- () const
    {
      return derived().opposite();
    }

//    template <typename OtherDerived>
//    bool operator== (const MotionBase<OtherDerived> & other) const {return derived().isEqual(other);}
//
//    Derived operator-() const { return derived().__minus__(); }
//    template <typename OtherDerived>
//    Derived operator+(const MotionBase<OtherDerived> & v2) const { return derived().__plus__(v2); }
//    template <typename OtherDerived>
//    Derived operator-(const MotionBase<OtherDerived> & v2) const { return derived().__minus__(v2); }
//    template <typename OtherDerived>
//    Derived & operator+=(const MotionBase<OtherDerived> & v2) { return derived().__pequ__(v2); }
    
//    template <typename OtherDerived>
//    void applyPlusEqual(MotionBase<OtherDerived> & other) const { derived().applyPlusEqual(other); }

    template<typename ForceDerived>
    ForceDerived cross(const ForceBase<ForceDerived> & f) const { return derived().cross(f.derived()); }
    
    template<typename OtherDerived>
    OtherDerived cross(const MotionBase<OtherDerived> & m) const { return derived().cross(m.derived()); }
    
    static Derived Zero() { return Derived::Zero(); }
    static Derived Random() { return Derived::Random(); }


  }; // class MotionBase


  template<typename T, int U>
  struct traits< MotionTpl<T,U> > : traits < MotionBase< MotionTpl<T,U> > >
  {
    typedef T Scalar_t;
    typedef MotionTpl<T,U> Type;
    typedef traits< MotionBase<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
    typedef T DotReturnType;
    typedef Eigen::Matrix<T,3,1,U> Vector3;
    typedef Eigen::Matrix<T,4,1,U> Vector4;
    typedef Eigen::Matrix<T,6,1,U> Vector6;
    typedef Eigen::Matrix<T,3,3,U> Matrix3;
    typedef Eigen::Matrix<T,4,4,U> Matrix4;
    typedef Eigen::Matrix<T,6,6,U> Matrix6;
    typedef Matrix6 ActionMatrix_t;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type Linear_t;
    typedef typename Vector6::template FixedSegmentReturnType<3>::Type Angular_t;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstLinear_t;
    typedef typename Vector6::template ConstFixedSegmentReturnType<3>::Type ConstAngular_t;
    typedef Eigen::Quaternion<T,U> Quaternion_t;
    typedef SE3Tpl<T,U> SE3;
    typedef ForceTpl<T,U> Force;
    typedef MotionTpl<T,U> Motion;
    typedef Symmetric3Tpl<T,U> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
  }; // traits MotionTpl


  template<typename _Scalar, int _Options>
  class MotionTpl : public MotionBase< MotionTpl<_Scalar,_Options> >
  {
  public:
    typedef MotionTpl Derived;
    
    SPATIAL_TYPEDEF_TEMPLATE(MotionTpl);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
  public:
    typedef class MotionBase< MotionTpl > Base;
//    friend class MotionBase< MotionTpl<_Scalar,_Options> >;
    typedef typename traits<MotionTpl>::SE3ActionReturnType SE3ActionReturnType;
    
    // Constructors
    MotionTpl() : data() {}

    template<typename V3,typename W3>
    MotionTpl(const Eigen::MatrixBase<V3> & v, const Eigen::MatrixBase<W3> & w)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(V3,3);
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(W3,3);
      data << v, w;
    }

    template<typename V6>
    MotionTpl(const Eigen::MatrixBase<V6> & v)
    : data(v)
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector6,V6);
    }


    template<typename OtherScalar, int OtherOptions>
    MotionTpl(const MotionTpl<OtherScalar, OtherOptions> & clone)
    : data(clone.coeffs())
    {}

    // initializers
    static MotionTpl Zero()   { return MotionTpl(Vector6::Zero());   }
    static MotionTpl Random() { return MotionTpl(Vector6::Random()); }

    MotionTpl & setZero () { data.setZero (); return *this; }
    MotionTpl & setRandom () { data.setRandom (); return *this; }

//    const Vector6 & toVector_impl() const { return data; }
//    Vector6 & toVector_impl() { return data; }

    inline ActionMatrix_t toActionMatrix_impl () const
    {
      ActionMatrix_t X;
      X.block <3,3> (ANGULAR, ANGULAR) = X.block <3,3> (LINEAR, LINEAR) = skew (angular());
      X.block <3,3> (LINEAR, ANGULAR) = skew (linear());
      X.block <3,3> (ANGULAR, LINEAR).setZero ();

      return X;
    }

    // Getters
    ConstAngular_t angular() const { return data.template segment<3> (ANGULAR); }
    ConstLinear_t linear()  const { return data.template segment<3> (LINEAR); }
    
    Angular_t angular() { return data.template segment<3> (ANGULAR); }
    Linear_t linear()  { return data.template segment<3> (LINEAR); }
    
    template<typename W3>
    void angular(const Eigen::MatrixBase<W3> & w)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(W3,3);
      data.template segment<3> (ANGULAR) = w;
    }
    
    template<typename V3>
    void linear(const Eigen::MatrixBase<V3> & v)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(V3,3);
      data.template segment<3> (LINEAR) = v;
    }
    
    template<typename V6>
    inline MotionTpl & operator= (const Eigen::MatrixBase<V6> & v)
    {
      EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(Vector6,V6);
      data = v;
      return *this;
    }
    
    template<typename OtherDerived>
    inline bool isEqual (const MotionBase<OtherDerived> & other) const
    {
      return other.derived().isEqual(*this);
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline bool isEqual (const MotionTpl<OtherScalar,OtherOptions> & other) const
    {
      return coeffs() == other.coeffs();
    }
    
//    template<typename OtherDerived>
//    void setTo (MotionBase<OtherDerived> &) const
//    {
//      PINOCCHIO_STATIC_ASSERT(true,YOU_CALLED_A_METHOD_WHICH_CANNOT_BE_APPLIED);
//    }
    
    template<typename OtherScalar, int OtherOptions>
    inline void setTo (MotionTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.coeffs() = data;
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline void addTo (MotionTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.coeffs() += data;
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline void subTo (MotionTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.coeffs() -= data;
    }
    
    template<typename OtherDerived>
    inline MotionTpl add (const MotionBase<OtherDerived> & other) const
    {
      MotionTpl dest (*this);
      other.derived().addTo(dest);
      return dest;
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline MotionTpl add (const MotionTpl<OtherScalar, OtherOptions> & other) const
    {
      return MotionTpl(data + other.coeffs());
    }
    
    template<typename OtherDerived>
    inline MotionTpl sub (const MotionBase<OtherDerived> & other) const
    {
      MotionTpl dest (*this);
      other.derived().subTo(dest);
      return dest;
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline MotionTpl sub (const MotionTpl<OtherScalar, OtherOptions> & other) const
    {
      return MotionTpl(data - other.coeffs());
    }
    
    inline MotionTpl opposite () const
    {
      return MotionTpl(-data);
    }
    
    template<typename SE3Scalar, int SE3Options>
    inline SE3ActionReturnType SE3ActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    {
      Vector3 Rw (M.rotation() * angular());
      return MotionTpl(M.rotation()*linear() + M.translation().cross(Rw),
                       Rw);
    }
    
    template<typename SE3Scalar, int SE3Options>
    inline SE3ActionReturnType SE3InvActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    {
      return MotionTpl(M.rotation().transpose()*(linear()-M.translation().cross(angular())),
                       M.rotation().transpose()*angular());
    }

//    Derived __minus__() const { return MotionTpl(-data); }
//    template <typename OtherDerived>
//    Derived __plus__(const MotionBase<OtherDerived> & v2) const { return Derived(data + v2.toVector()); }
//    template <typename OtherDerived>
//    Derived __minus__(const MotionBase<OtherDerived> & v2) const { return Derived(data - v2.toVector()); }
//    template <typename OtherDerived>
//    Derived & __pequ__(const MotionBase<OtherDerived> & v2) { v2.applyPlusEqual(*this); return *this; }
    
    template<typename ForceDerived>
    inline Scalar_t dot(const ForceBase<ForceDerived> & f) const
    {
      return f.derived().dot(*this);
    }
    
    template<typename ForceScalar, int ForceOptions>
    inline Scalar_t dot(const ForceTpl<ForceScalar, ForceOptions> & f) const
    {
      return data.dot(f.coeffs());
    }

    template<typename OtherScalar, int OtherOptions>
    inline MotionTpl cross(const MotionTpl<OtherScalar,OtherOptions> & other) const
    {
      return MotionTpl(linear().cross(other.angular())+angular().cross(other.linear()),
                       angular().cross(other.angular()) );
    }

    template<typename ForceScalar, int ForceOptions>
    inline ForceTpl<ForceScalar,ForceOptions> cross(const ForceTpl<ForceScalar,ForceOptions> & f) const
    {
      typedef ForceTpl<ForceScalar,ForceOptions> ReturnType;
      return ReturnType(angular().cross(f.linear()),
                        angular().cross(f.angular())+linear().cross(f.linear()) );
    }


    void disp_impl(std::ostream & os) const
    {
      using namespace std;
      os
      << "  v = " << linear().transpose () << endl
      << "  w = " << angular().transpose () << endl;
    }
    
    Vector6 & coeffs() { return data; }
    const Vector6 & coeffs() const { return data; }
    operator Vector6 () const { return coeffs(); }

//    /** \brief Compute the classical acceleration of point according to the spatial velocity and spatial acceleration of the frame centered on this point
//     */
//    static inline Vector3 computeLinearClassicalAcceleration (const MotionTpl & spatial_velocity, const MotionTpl & spatial_acceleration)
//    {
//      return spatial_acceleration.linear () + spatial_velocity.angular ().cross (spatial_velocity.linear ());
//    }
//
//    /**
//      \brief Compute the spatial motion quantity of the parallel frame translated by translation_vector
//     */
//    MotionTpl translate (const Vector3 & translation_vector) const
//    {
//      return MotionTpl (linear() + angular().cross (translation_vector), angular());
//    }

  protected:
    Vector6 data;

  }; // class MotionTpl

  template<typename S,int O>
  MotionTpl<S,O> operator^( const MotionTpl<S,O> &m1, const MotionTpl<S,O> &m2 ) { return m1.cross(m2); }
  template<typename S,int O>
  ForceTpl<S,O> operator^( const MotionTpl<S,O> &m, const ForceTpl<S,O> &f ) { return m.cross(f); }

  typedef MotionTpl<double,0> Motion;


  ///////////////   BiasZero  ///////////////
  struct BiasZero;

  template<>
  struct traits< BiasZero > : traits < MotionBase<BiasZero> >
  {
    typedef double Scalar_t;
    typedef BiasZero Type;
    typedef traits< MotionBase<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
    typedef Scalar_t DotReturnType;
    typedef Eigen::Matrix<double,3,1,0> Vector3;
    typedef Eigen::Matrix<double,4,1,0> Vector4;
    typedef Eigen::Matrix<double,6,1,0> Vector6;
    typedef Eigen::Matrix<double,3,3,0> Matrix3;
    typedef Eigen::Matrix<double,4,4,0> Matrix4;
    typedef Eigen::Matrix<double,6,6,0> Matrix6;
    typedef Matrix6 ActionMatrix_t;
    typedef Vector3 Angular_t;
    typedef const Vector3 ConstAngular_t;
    typedef Vector3 Linear_t;
    typedef const Vector3 ConstLinear_t;
    typedef Eigen::Quaternion<double,0> Quaternion_t;
    typedef SE3Tpl<double,0> SE3;
    typedef ForceTpl<double,0> Force;
    typedef MotionTpl<double,0> Motion;
    typedef Symmetric3Tpl<double,0> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
  }; // traits BiasZero

  struct BiasZero : public MotionBase<BiasZero>
  {
    SPATIAL_TYPEDEF_NO_TEMPLATE(BiasZero);
    
    template <typename OtherDerived>
    operator typename MotionBase<OtherDerived>::Derived () const
    {
      return OtherDerived::Zero();
    }
    
    template <typename Scalar, int Options>
    operator MotionTpl<Scalar,Options> () const
    {
      return MotionTpl<Scalar,Options>::Zero();
    }
    
    template<typename Derived>
    const MotionBase<Derived> & add (const MotionBase<Derived> & other) const
    {
      return other;
    }
    
    template<typename OtherDerived>
    void addTo (MotionBase<OtherDerived> &) const
    {
    }
    
    template<typename OtherDerived>
    void subTo (MotionBase<OtherDerived> &) const
    {
    }
  }; // struct BiasZero

  template <typename Derived>
  inline const Derived & operator+ (const MotionBase<Derived> & v, const BiasZero &) { return v.derived(); }
  
  template <typename Derived>
  inline const Derived & operator- (const MotionBase<Derived> & v, const BiasZero &) { return v.derived(); }
  
  template <typename Derived>
  inline const Derived & operator+ (const BiasZero &, const MotionBase<Derived> & v) { return v.derived(); }
  
  template <typename Derived>
  inline const Derived & operator- (const BiasZero &, const MotionBase<Derived> & v) { return v.derived().opposite(); }

} // namespace se3

#endif // ifndef __se3_motion_hpp__
