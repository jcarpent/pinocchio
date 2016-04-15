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

#ifndef __se3_se3_hpp__
#define __se3_se3_hpp__

#include <Eigen/Geometry>
#include "pinocchio/spatial/spatial.hpp"
#include "pinocchio/spatial/skew.hpp"
#include "pinocchio/math/sincos.hpp"


namespace se3
{

  /* Type returned by the "se3Action" and "se3ActionInverse" functions. */
//  namespace internal 
//  {
//    template<typename D>
//    struct ActionReturn    { typedef D Type; };
//  }
  
  template<class Derived>
  struct traits< SE3Base<Derived> >
  {
    typedef Derived SE3ActionReturnType;
    typedef Derived SE3ProductReturnType;
  };
  
  template<class LHSDerived, class RHSDerived> class SE3ProductType;

  /** The rigid transform aMb can be seen in two ways: 
   *
   * - given a point p expressed in frame B by its coordinate vector Bp, aMb
   * computes its coordinates in frame A by Ap = aMb Bp.
   * - aMb displaces a solid S centered at frame A into the solid centered in
   * B. In particular, the origin of A is displaced at the origin of B: $^aM_b
   * ^aA = ^aB$.

   * The rigid displacement is stored as a rotation matrix and translation vector by:
   * aMb (x) =  aRb*x + aAB
   * where aAB is the vector from origin A to origin B expressed in coordinates A.
   */
  template<class Derived>
  class SE3Base : public SpatialBase<Derived>
  {
    
  public:
    
    typedef SpatialBase<Derived> Base;
    using Base::derived;
    
    template<typename OtherDerived>
    inline bool operator== (const SE3Base<OtherDerived> & other) const
    {
      return derived().isEqual(other.derived());
    }
    
    template<typename OtherDerived>
    inline Derived & operator= (const SE3Base<OtherDerived> & other)
    {
      other.derived().setTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline typename traits< SE3ProductType<Derived,OtherDerived> >::SE3ProductReturnType operator* (const SE3Base<OtherDerived> & other) const
    {
      typedef typename traits< SE3ProductType<Derived,OtherDerived> >::SE3ProductReturnType ReturnType;
//      return derived().mult(other.derived());
      return ReturnType(derived(),other.derived());
    }
    
    template<typename OtherDerived>
    inline void applyThisOnTheRight (SE3Base<OtherDerived> & dest) const
    {
      dest = dest * derived();
    }
    
    template<typename OtherDerived>
    inline void applyThisOnTheLeft (SE3Base<OtherDerived> & dest) const
    {
      dest = derived() * dest;
    }
    
    inline Derived inverse () const
    {
      return derived().inverse();
    }

    /// ay = aXb.act(by)
    template<typename SpatialDerived>
    inline typename SpatialDerived::SE3ActionReturnType act(const SpatialBase<SpatialDerived> & spatial_obj) const
    {
      return spatial_obj.derived().SE3ActOn(derived());
    }
    
    template<typename EigenDerived>
    inline typename traits<Derived>::SE3ActionOnVecReturnType
    act(const Eigen::MatrixBase<EigenDerived> & vec) const
    {
      return derived().act(vec);
    }
    
    /// by = aXb.actInv(ay)
    template<typename SpatialDerived>
    inline typename SpatialDerived::SE3ActionReturnType actInv(const SpatialBase<SpatialDerived> & spatial_obj) const
    {
      return derived().actInv(spatial_obj.derived());
    }
    
    template<typename EigenDerived, typename OtherEigenDerived>
    inline typename Eigen::MatrixBase<OtherEigenDerived>
    actInv(const Eigen::MatrixBase<EigenDerived> & vec) const
    {
      return derived().actInv(vec);
    }

    template<typename OtherDerived>
    inline bool isApprox(const SE3Base<OtherDerived> & other, const typename traits<Derived>::Scalar_t & prec = Eigen::NumTraits<typename traits<Derived>::Scalar_t>::dummy_precision()) const
    {
      return derived().isApprox(other,prec);
    }
    
    static inline Derived Identity()
    {
      return Derived::Identity();
    }
    
    static inline Derived Random()
    {
      return Derived::Random();
    }
  }; // class SE3Base
  
  namespace internal {
    template<class LHSDerived, class RHSDerived, class Dest>
    inline void SE3_product(const SE3ProductType<LHSDerived,RHSDerived> & prod, Dest & dest);
    
    
    template<typename Scalar, int Options>
    inline void SE3_product(const SE3ProductType< SE3Tpl<Scalar,Options>,SE3Tpl<Scalar,Options> > & prod, SE3Tpl<Scalar,Options> & dest)
    {
      dest.rotation() = prod.lhs().rotation() * prod.rhs().rotation();
      dest.translation() = prod.lhs().translation() + prod.lhs().rotation() * prod.rhs().translation();
    }
  }
  
  template<int axis> struct SE3Revolute;
  
  namespace internal {
    
    template<typename Scalar, int Options>
    inline void SE3_product(const SE3ProductType< SE3Tpl<Scalar,Options>, SE3Revolute<0> > & prod, SE3Tpl<Scalar,Options> & dest)
    {
      const double & ca = prod.rhs().ca();
      const double & sa = prod.rhs().sa();
//      Eigen::Matrix3d R3;
//      R3 <<
//      1,0,0,
//      0,ca,-sa,
//      0,sa,ca;
//      dest.rotation() = R3*prod.lhs().rotation();
      dest.rotation().col(0) = prod.lhs().rotation().col(0);
      dest.rotation().col(1) = ca * prod.lhs().rotation().col(1) + sa * prod.lhs().rotation().col(2);
      dest.rotation().col(2) = -sa * prod.lhs().rotation().col(1) + ca * prod.lhs().rotation().col(2);
      dest.translation() = prod.lhs().translation();
    }
    
    template<typename Scalar, int Options>
    inline void SE3_product(const SE3ProductType< SE3Tpl<Scalar,Options>, SE3Revolute<1> > & prod, SE3Tpl<Scalar,Options> & dest)
    {
      const double & ca = prod.rhs().ca();
      const double & sa = prod.rhs().sa();
      
      dest.rotation().col(1) = prod.lhs().rotation().col(1);
      dest.rotation().col(0) = ca * prod.lhs().rotation().col(0) - sa * prod.lhs().rotation().col(2);
      dest.rotation().col(2) = sa * prod.lhs().rotation().col(0) + ca * prod.lhs().rotation().col(2);
      dest.translation() = prod.lhs().translation();
    }
    
    template<typename Scalar, int Options>
    inline void SE3_product(const SE3ProductType< SE3Tpl<Scalar,Options>, SE3Revolute<2> > & prod, SE3Tpl<Scalar,Options> & dest)
    {
      const double & ca = prod.rhs().ca();
      const double & sa = prod.rhs().sa();
      
      dest.rotation().col(2) = prod.lhs().rotation().col(2);
      dest.rotation().col(0) = ca * prod.lhs().rotation().col(0) + sa * prod.lhs().rotation().col(1);
      dest.rotation().col(1) = -sa * prod.lhs().rotation().col(0) + ca * prod.lhs().rotation().col(1);
      dest.translation() = prod.lhs().translation();
    }
  }
  
  template<class LHSDerived, class RHSDerived>
  struct traits< SE3ProductType<LHSDerived,RHSDerived> > : traits< SE3Base< SE3ProductType<LHSDerived,RHSDerived> > >
  {
    typedef SE3ProductType<LHSDerived,RHSDerived> Type;
    typedef traits< SE3Base<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
  };
  
  template<typename Scalar, int Options, class OtherDerived>
  struct traits< SE3ProductType<SE3Tpl<Scalar,Options>, OtherDerived> > : traits< SE3Base< SE3ProductType<SE3Tpl<Scalar,Options>, OtherDerived> > >
  {
    typedef SE3ProductType<SE3Tpl<Scalar,Options>, OtherDerived> Type;
    typedef traits< SE3Base<Type> > BaseTraits;
//    using typename BaseTraits::SE3ActionReturnType;
    typedef Type SE3ProductReturnType;
    typedef SE3Tpl<Scalar,Options> LHSType;
    typedef OtherDerived RHSType;
    typedef Scalar Scalar_t;
    typedef typename traits<LHSType>::SE3ActionOnVecReturnType SE3ActionOnVecReturnType;
    typedef LHSType SE3ActionReturnType;
  };
  
  template<class LHSDerived, class RHSDerived>
  class SE3ProductType : public SE3Base< SE3ProductType<LHSDerived,RHSDerived> >
  {
  public:
    typedef typename traits<SE3ProductType>::SE3ActionReturnType SE3ActionReturnType;
    
    SE3ProductType(const SE3Base<LHSDerived> & lhs, const SE3Base<RHSDerived> & rhs) : terms(lhs.derived(), rhs.derived()) {}//lhs_(lhs.derived()), rhs_(rhs.derived()) {}
    
    template<typename DestDerived>
    inline void evalTo(SE3Base<DestDerived> & dest) const
    {
      internal::SE3_product(*this,dest.derived());
    }
    
    template<typename DestDerived>
    inline void setTo(SE3Base<DestDerived> & dest) const
    {
      evalTo(dest);
    }
    
    inline SE3ActionReturnType eval() const
    {
//      std::cout << "bad news" << std::endl;
      SE3ActionReturnType res;
      evalTo(res);
      return res;
    }
    
    operator SE3ActionReturnType () const
    {
      return eval();
    }
    
//    const LHSDerived & lhs () const { return lhs_; }
//    const RHSDerived & rhs () const { return rhs_; }
    
    const LHSDerived & lhs () const { return terms.first; }
    const RHSDerived & rhs () const { return terms.second; }
    
  protected:
    
//    const SE3Base<LHSDerived> & lhs_;
//    const SE3Base<RHSDerived> & rhs_;
    
//    const LHSDerived & lhs_;
//    const RHSDerived & rhs_;
    
    std::pair<const LHSDerived &,const RHSDerived &> terms;
  };
  
  
  
  
//  template<typename Scalar, int Options>
//  template<typename DestScalar, int DestOptions>
//  void SE3ProductType< SE3Tpl<Scalar,Options>, SE3Tpl<Scalar,Options> >::evalTo(SE3Tpl<DestScalar, DestOptions> & dest) const
//  {
//    
//  }

  

  template<typename T, int U>
  struct traits< SE3Tpl<T,U> > : traits< SE3Base< SE3Tpl<T,U> > >
  {
    typedef SE3Tpl<T,U> Type;
    typedef traits< SE3Base<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
    
    typedef T Scalar_t;
    typedef Eigen::Matrix<T,3,1,U> SE3ActionOnVecReturnType;
    typedef Eigen::Matrix<T,3,1,U> Vector3;
    typedef Eigen::Matrix<T,4,1,U> Vector4;
    typedef Eigen::Matrix<T,6,1,U> Vector6;
    typedef Eigen::Matrix<T,3,3,U> Matrix3;
    typedef Eigen::Matrix<T,4,4,U> Matrix4;
    typedef Eigen::Matrix<T,6,6,U> Matrix6;
    typedef Matrix3 Angular_t;
    typedef const Matrix3 ConstAngular_t;
    typedef Vector3 Linear_t;
    typedef const Vector3 ConstLinear_t;
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
  }; // traits SE3Tpl

  template<typename _Scalar, int _Options>
  class SE3Tpl : public SE3Base< SE3Tpl< _Scalar, _Options > >
  {

  public:
    typedef _Scalar Scalar;
//    friend class SE3Base< SE3Tpl< _Scalar, _Options > >;
    typedef SE3Base< SE3Tpl<_Scalar,_Options> > Base;
    typedef typename traits<SE3Tpl>::SE3ActionReturnType SE3ActionReturnType;
    
    SPATIAL_TYPEDEF_TEMPLATE(SE3Tpl);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


    SE3Tpl(): rot(), trans() {};
//    using Base::operator*;
    using Base::act;
    using Base::operator=;


    template<typename M3,typename V3>
    SE3Tpl(const Eigen::MatrixBase<M3> & R, const Eigen::MatrixBase<V3> & p)
    : rot(R), trans(p)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(V3,3)
      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(M3,3,3)
    }

    template<typename M4>
    SE3Tpl(const Eigen::MatrixBase<M4> & m) 
    : rot(m.template block<3,3>(LINEAR,LINEAR)), trans(m.template block<3,1>(LINEAR,ANGULAR))
    {
      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(M4,4,4);
    }

    SE3Tpl(int) : rot(Angular_t::Identity()), trans(Linear_t::Zero()) {}

    template<typename OtherScalar, int OtherOptions>
    SE3Tpl(const SE3Tpl<OtherScalar,OtherOptions> & clone)
    : rot(clone.rotation()),trans(clone.translation()) {}


//    template<typename OtherDerived>
//    SE3Tpl & operator= (const SE3Tpl<OtherDerived> & other)
//    {
//      rot = other.rotation (); trans = other.translation (); return *this;
//    }
    
    template<typename OtherScalar, int OtherOptions>
    inline void setTo(SE3Tpl<OtherScalar,OtherOptions> & M) const
    {
      M.translation() = translation();
      M.rotation() = rotation();
    }

    inline static SE3Tpl Identity()
    {
      return SE3Tpl(1);
    }

    inline SE3Tpl & setIdentity () { rot.setIdentity (); trans.setZero (); return *this;}

    /// aXb = bXa.inverse()
    inline SE3Tpl inverse() const
    {
      return SE3Tpl(rot.transpose(), -rot.transpose()*trans);
    }

    inline static SE3Tpl Random()
    {
      Quaternion_t q(Vector4::Random());
      q.normalize();
      return SE3Tpl(q.matrix(),Vector3::Random());
    }

    inline SE3Tpl & setRandom ()
    {
      Quaternion_t q(Vector4::Random());
      q.normalize ();
      rot = q.matrix ();
      trans.setRandom ();

      return *this;
    }
    
    inline Matrix4 toHomogeneousMatrix() const
    {
      Matrix4 M;
      M.template block<3,3>(LINEAR,LINEAR) = rot;
      M.template block<3,1>(LINEAR,ANGULAR) = trans;
      M.template block<1,3>(ANGULAR,LINEAR).setZero();
      M(3,3) = 1;
      return M;
    }
    
    inline operator Matrix4 () const
    {
      return toHomogeneousMatrix();
    }

    /// Vb.toVector() = bXa.toMatrix() * Va.toVector()
    inline Matrix6 toActionMatrix() const
    {
      typedef Eigen::Block<Matrix6,3,3> Block3;
      Matrix6 M;
      M.template block<3,3>(ANGULAR,ANGULAR)
      = M.template block<3,3>(LINEAR,LINEAR) = rot;
      M.template block<3,3>(ANGULAR,LINEAR).setZero();
      Block3 B = M.template block<3,3>(LINEAR,ANGULAR);
      
      B.col(0) = trans.cross(rot.col(0));
      B.col(1) = trans.cross(rot.col(1));
      B.col(2) = trans.cross(rot.col(2));
      return M;
    }
    
    inline operator Matrix6 () const
    {
      return toActionMatrix();
    }

    void disp(std::ostream & os) const
    {
      using namespace std;
      os
      << "  R =\n" << rot << endl
      << "  p = " << trans.transpose() << endl;
    }
    
    template<typename OtherDerived>
    inline bool isEqual (const SE3Base<OtherDerived> & other) const
    {
      return other.derived().isEqual(*this);
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline bool isEqual (const SE3Tpl<OtherScalar,OtherOptions> & other) const
    {
      return rotation() == other.rotation() && translation() == other.translation();
    }

    /// --- GROUP ACTIONS ON M6, F6 and I6 ---
    
//    using Base::operator*;
    
    template<typename OtherDerived>
    inline typename OtherDerived::SE3ActionReturnType act(const SE3Base<OtherDerived> & other) const
    {
      return other.derived().SE3ActOn(*this);
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline SE3Tpl SE3ActOn(const SE3Tpl<OtherScalar,OtherOptions> & other) const
    {
      return (other * (*this)).eval();
    }
    
    template<typename EigenDerived>
    inline typename traits<SE3Tpl>::SE3ActionOnVecReturnType
    act(const Eigen::MatrixBase<EigenDerived> & vec) const
    {
      EIGEN_STATIC_ASSERT(EigenDerived::RowsAtCompileTime==3,THIS_METHOD_IS_ONLY_FOR_OBJECTS_OF_A_SPECIFIC_SIZE);
      return translation() + rotation()*vec;
    }
    
    template<typename EigenDerived>
    inline typename traits<SE3Tpl>::SE3ActionOnVecReturnType
    actInv(const Eigen::MatrixBase<EigenDerived> & vec) const
    {
      EIGEN_STATIC_ASSERT(EigenDerived::RowsAtCompileTime==3,THIS_METHOD_IS_ONLY_FOR_OBJECTS_OF_A_SPECIFIC_SIZE);
      return rotation().transpose()*(vec-translation());
    }
    
    /// by = aXb.actInv(ay)
    template<typename SpatialDerived>
    inline typename SpatialDerived::SE3ActionReturnType actInv(const SpatialBase<SpatialDerived> & spatial_obj) const
    {
      return spatial_obj.SE3InvActOn(*this);
    }

    template<typename OtherDerived>
    inline SE3Tpl mult(const SE3Base<OtherDerived> & other) const
    {
      return act(other);
    }

    template<typename OtherScalar, int OtherOptions>
    inline bool isApprox (const SE3Tpl<OtherScalar,OtherOptions> & M2, const Scalar_t & prec = Eigen::NumTraits<Scalar_t>::dummy_precision()) const
    {
      return rot.isApprox(M2.rot,prec) && trans.isApprox(M2.trans,prec);
    }

    ConstAngular_t & rotation() const { return rot; }
    Angular_t & rotation() { return rot; }
    template<typename D>
    inline void rotation(const Eigen::MatrixBase<D> & R)
    {
      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(D,3,3);
      rot = R;
    }
    
    ConstLinear_t & translation() const { return trans;}
    Linear_t & translation() { return trans;}
    template<typename D>
    inline void translation(const Eigen::MatrixBase<D> & t)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(D,3); trans=t;
    }

  protected:
    Angular_t rot;
    Linear_t trans;
    
  }; // class SE3Tpl

  typedef SE3Tpl<double,0> SE3;

} // namespace se3



#endif // ifndef __se3_se3_hpp__

