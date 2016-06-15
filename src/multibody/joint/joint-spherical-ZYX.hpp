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

#ifndef __se3_joint_spherical_ZYX_hpp__
#define __se3_joint_spherical_ZYX_hpp__
#include <iostream>
#include "pinocchio/multibody/joint/joint-base.hpp"
#include "pinocchio/multibody/joint/joint-dense.hpp"
#include "pinocchio/multibody/constraint.hpp"
#include "pinocchio/math/sincos.hpp"
#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/spatial/skew.hpp"

#include <stdexcept>

namespace se3
{

  struct JointDataSphericalZYX;
  struct JointModelSphericalZYX;
  
  template <typename _Scalar, int _Options=0> struct JointSphericalZYXTpl;

  
  template <typename _Scalar, int _Options>
  struct JointSphericalZYXTpl
  {
    typedef _Scalar Scalar;
    enum { Options = _Options };
    typedef Eigen::Matrix<Scalar,3,1,Options> Vector3;
    typedef Eigen::Matrix<Scalar,3,3,Options> Matrix3;
    typedef Eigen::Matrix<Scalar,6,1,Options> Vector6;
    typedef Eigen::Matrix<Scalar,6,6,Options> Matrix6;
    typedef MotionTpl<Scalar,Options> Motion;
    typedef ForceTpl<Scalar,Options> Force;
    typedef SE3Tpl<Scalar,Options> SE3;

//    struct BiasSpherical
//    {
//      typename MotionTpl<Scalar,Options>::Vector3 c_J;
//
//      BiasSpherical ()  {c_J.fill (NAN);}
//      //BiasSpherical (const Motion::Vector3 & c_J) c_J (c_J) {}
//
//      template<typename OtherScalar, int OtherOptions>
//      operator MotionTpl<OtherScalar,OtherOptions> () const
//      {
//        typedef MotionTpl<OtherScalar,OtherOptions> ReturnType;
//        return  ReturnType(ReturnType::Linear_t::Zero (), c_J);
//      }
//
//      typename MotionTpl<Scalar,Options>::Vector3 & operator() () { return c_J; }
//      const typename MotionTpl<Scalar,Options>::Vector3 & operator() () const { return c_J; }
//      
//      template<typename OtherScalar, int OtherOptions>
//      MotionTpl<OtherScalar, OtherOptions> add (const MotionTpl<OtherScalar, OtherOptions> & other) const
//      {
//        
//        return MotionTpl(data + other.coeffs());
//      }
//
//
//    }; // struct BiasSpherical

//    template<typename OtherScalar, int OtherOptions>
//    inline friend const Motion operator+ (const MotionTpl<OtherScalar, OtherOptions> & v, const BiasSpherical & c)
//    {
//      typedef
//      return Motion (v.linear (), v.angular () + c ());
//    }
//    inline friend const Motion operator+ (const BiasSpherical & c, const MotionTpl<OtherScalar, OtherOptions> & v) { return Motion (v.linear (), v.angular () + c ()); }

    struct MotionSpherical : MotionBase<MotionSpherical>
    {
      MotionSpherical ()                       {w.fill(NAN);}
      MotionSpherical (const typename MotionTpl<Scalar,Options>::Vector3 & w) : w (w)  {}
      typename MotionTpl<Scalar,Options>::Vector3 w;

      typename MotionTpl<Scalar,Options>::Vector3 & operator() () { return w; }
      const typename MotionTpl<Scalar,Options>::Vector3 & operator() () const { return w; }

      template<typename OtherScalar, int OtherOptions>
      inline operator MotionTpl<OtherScalar, OtherOptions> () const
      {
        typedef MotionTpl<OtherScalar, OtherOptions> ReturnType;
        typedef typename ReturnType::Vector3 Vector3;
        //return Motion (typename MotionTpl<_Scalar,_Options>::Vector3::Zero (), w);
        return ReturnType(Vector3::Zero(),w);
      }
      
      template<typename OtherScalar, int OtherOptions>
      MotionTpl<OtherScalar, OtherOptions> add (const MotionTpl<OtherScalar, OtherOptions> & other) const
      {
        typedef MotionTpl<OtherScalar, OtherOptions> ReturnType;
        ReturnType res (other);
        res.angular() += w;
        return res;
      }
      
      template<typename OtherScalar, int OtherOptions>
      MotionTpl<OtherScalar, OtherOptions> sub (const MotionTpl<OtherScalar, OtherOptions> & other) const
      {
        typedef MotionTpl<OtherScalar, OtherOptions> ReturnType;
        ReturnType res (other);
        res.angular() -= w;
        return w;
      }
      
      template<typename OtherScalar, int OtherOptions>
      inline void addTo (MotionTpl<OtherScalar, OtherOptions> & dest) const
      {
        dest.angular() += w;
      }
      
      template<typename OtherScalar, int OtherOptions>
      inline void subTo (MotionTpl<OtherScalar, OtherOptions> & dest) const
      {
        dest.angular() -= w;
      }
    }; // struct MotionSpherical

//    inline friend const MotionSpherical operator+ (const MotionSpherical & m, const BiasSpherical & c)
//    { return MotionSpherical (m.w + c.c_J); }

    template<typename MotionScalar,int MotionOptions>
    friend MotionTpl<MotionScalar,MotionOptions>
    operator+ (const MotionSpherical & m1, const MotionTpl<MotionScalar,MotionOptions> & m2)
    {
      typedef MotionTpl<MotionScalar,MotionOptions> ReturnType;
      return ReturnType (m2.linear(),m2.angular()+m1.w);
    }
    
    struct ConstraintZYXRotationalSubspace;
    

    template<typename EigenDerived>
    friend MotionTpl<Scalar,Options>
    operator* (const ConstraintZYXRotationalSubspace & S, const Eigen::MatrixBase<EigenDerived> & v)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(EigenDerived,3);
      typedef MotionTpl<Scalar,Options> ReturnType;
      typedef typename ReturnType::Vector3 Vector3;
      return ReturnType(Vector3::Zero(),S()*v);
    }

  }; // struct JointSphericalZYXTpl

  template <>
  struct traits< typename JointSphericalZYXTpl<double,0>::MotionSpherical >   {
    typedef double Scalar_t;
    typedef MotionTpl<double,0> SE3ActionReturnType;
    typedef Scalar_t DotReturnType;
    typedef Eigen::Matrix<double,3,1,0> Vector3;
    typedef Eigen::Matrix<double,4,1,0> Vector4;
    typedef Eigen::Matrix<double,6,1,0> Vector6;
    typedef Eigen::Matrix<double,3,3,0> Matrix3;
    typedef Eigen::Matrix<double,4,4,0> Matrix4;
    typedef Eigen::Matrix<double,6,6,0> Matrix6;
    typedef Matrix3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Matrix3 ConstAngular_t;
    typedef const Vector3 ConstLinear_t;
    typedef Matrix6 ActionMatrix_t;
    typedef Eigen::Quaternion<double,0> Quaternion_t;
    typedef SE3Tpl<double,0> SE3;
    typedef ForceTpl<double,0> Force;
    typedef MotionTpl<double,0> Motion;
    typedef Symmetric3Tpl<double,0> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointMotion;
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointForce;
    typedef Eigen::Matrix<Scalar_t,6,3> DenseBase;
  };
  
  template <>
  struct traits< typename JointSphericalZYXTpl<double,2>::MotionSpherical >   {
    typedef double Scalar_t;
    typedef MotionTpl<double,2> SE3ActionReturnType;
    typedef Scalar_t DotReturnType;
    typedef Eigen::Matrix<double,3,1,0> Vector3;
    typedef Eigen::Matrix<double,4,1,0> Vector4;
    typedef Eigen::Matrix<double,6,1,0> Vector6;
    typedef Eigen::Matrix<double,3,3,0> Matrix3;
    typedef Eigen::Matrix<double,4,4,0> Matrix4;
    typedef Eigen::Matrix<double,6,6,0> Matrix6;
    typedef Matrix3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Matrix3 ConstAngular_t;
    typedef const Vector3 ConstLinear_t;
    typedef Matrix6 ActionMatrix_t;
    typedef Eigen::Quaternion<double,0> Quaternion_t;
    typedef SE3Tpl<double,0> SE3;
    typedef ForceTpl<double,0> Force;
    typedef MotionTpl<double,0> Motion;
    typedef Symmetric3Tpl<double,0> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointMotion;
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointForce;
    typedef Eigen::Matrix<Scalar_t,6,3> DenseBase;
  };
  
  template <>
  struct traits< JointSphericalZYXTpl<double,0>::ConstraintZYXRotationalSubspace > : traits < ConstraintBase< JointSphericalZYXTpl<double,0>::ConstraintZYXRotationalSubspace > >
  {
    typedef Eigen::Matrix <double,6,3,0> SE3ActionReturnType;
    
    typedef double Scalar_t;
    typedef Eigen::Matrix<double,3,1,0> Vector3;
    typedef Eigen::Matrix<double,4,1,0> Vector4;
    typedef Eigen::Matrix<double,6,1,0> Vector6;
    typedef Eigen::Matrix<double,3,3,0> Matrix3;
    typedef Eigen::Matrix<double,4,4,0> Matrix4;
    typedef Eigen::Matrix<double,6,6,0> Matrix6;
    typedef Matrix3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Matrix3 ConstAngular_t;
    typedef const Vector3 ConstLinear_t;
    typedef Matrix6 ActionMatrix_t;
    typedef Eigen::Quaternion<double,0> Quaternion_t;
    typedef SE3Tpl<double,0> SE3;
    typedef ForceTpl<double,0> Force;
    typedef MotionTpl<double,0> Motion;
    typedef Symmetric3Tpl<double,0> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointMotion;
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointForce;
    typedef Eigen::Matrix<Scalar_t,6,3> DenseBase;
  };
  
  template <>
  struct traits< JointSphericalZYXTpl<double,2>::ConstraintZYXRotationalSubspace > : traits < ConstraintBase< JointSphericalZYXTpl<double,0>::ConstraintZYXRotationalSubspace > >
  {
    typedef Eigen::Matrix <double,6,3,0> SE3ActionReturnType;
    
    typedef double Scalar_t;
    typedef Eigen::Matrix<double,3,1,0> Vector3;
    typedef Eigen::Matrix<double,4,1,0> Vector4;
    typedef Eigen::Matrix<double,6,1,0> Vector6;
    typedef Eigen::Matrix<double,3,3,0> Matrix3;
    typedef Eigen::Matrix<double,4,4,0> Matrix4;
    typedef Eigen::Matrix<double,6,6,0> Matrix6;
    typedef Matrix3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Matrix3 ConstAngular_t;
    typedef const Vector3 ConstLinear_t;
    typedef Matrix6 ActionMatrix_t;
    typedef Eigen::Quaternion<double,0> Quaternion_t;
    typedef SE3Tpl<double,0> SE3;
    typedef ForceTpl<double,0> Force;
    typedef MotionTpl<double,0> Motion;
    typedef Symmetric3Tpl<double,0> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointMotion;
    typedef Eigen::Matrix<Scalar_t,3,1,0> JointForce;
    typedef Eigen::Matrix<Scalar_t,6,3> DenseBase;
  };

  
  template<typename _Scalar, int _Options>
  struct JointSphericalZYXTpl<_Scalar,_Options>::ConstraintZYXRotationalSubspace : ConstraintBase <ConstraintZYXRotationalSubspace>
  {
    enum { NV = 3, Options = 0 };
    typedef _Scalar Scalar;
    typedef Eigen::Matrix <_Scalar,3,3,_Options> Matrix3;
    typedef Eigen::Matrix <_Scalar,3,1,_Options> Vector3;
    typedef Eigen::Matrix <_Scalar,6,3,_Options> ConstraintDense;
    typedef Eigen::Matrix <_Scalar,6,3,_Options> SE3ActionReturnType;
    
    Matrix3 S_minimal;
    
    Motion operator* (const MotionSpherical & vj) const
    { return Motion (Motion::Vector3::Zero (), S_minimal * vj ()); }
    
    ConstraintZYXRotationalSubspace () : S_minimal () { S_minimal.fill (NAN); }
    ConstraintZYXRotationalSubspace (const Matrix3 & subspace) : S_minimal (subspace) {}
    
    Matrix3 & operator() () { return S_minimal; }
    const Matrix3 & operator() () const { return S_minimal; }
    
    Matrix3 & matrix () { return S_minimal; }
    const Matrix3 & matrix () const { return S_minimal; }
    
    int nv_impl() const { return NV; }
    
    struct ConstraintTranspose
    {
      const ConstraintZYXRotationalSubspace & ref;
      ConstraintTranspose(const ConstraintZYXRotationalSubspace & ref) : ref(ref) {}
      
      
      const typename Eigen::ProductReturnType<
      Eigen::Transpose<const Matrix3>,
      //        typename Motion::ConstAngular_t::Base /* This feature leads currently to a bug */
      Eigen::Block<const typename Force::Vector6,3,1>
      >::Type
      operator* (const Force & phi) const
      {
        return ref.S_minimal.transpose () * phi.angular();
      }
      
      /* [CRBA]  MatrixBase operator* (Constraint::Transpose S, ForceSet::Block) */
      template<typename D>
      const typename Eigen::ProductReturnType<
      typename Eigen::Transpose<const Matrix3>,
      typename Eigen::MatrixBase<const D>::template NRowsBlockXpr<3>::Type
      >::Type
      operator* (const Eigen::MatrixBase<D> & F) const
      {
        EIGEN_STATIC_ASSERT(D::RowsAtCompileTime==6,THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE)
        return ref.S_minimal.transpose () * F.template bottomRows<3> ();
      }
    }; // struct ConstraintTranspose
    
    ConstraintTranspose transpose () const { return ConstraintTranspose(*this); }
    
    operator ConstraintXd () const
    {
      ConstraintDense S;
      (S.template block <3,3> (Inertia::LINEAR, 0)).setZero ();
      S.template block <3,3> (Inertia::ANGULAR, 0) = S_minimal;
      return ConstraintXd(S);
    }
    
    //      const typename Eigen::ProductReturnType<
    //      const ConstraintDense,
    //      const Matrix3
    //      >::Type
    template<typename SE3Scalar, int SE3Options>
    SE3ActionReturnType SE3ActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    {
      //        Eigen::Matrix <Scalar,6,3,Options> X_subspace;
      //        X_subspace.template block <3,3> (Motion::LINEAR, 0) = skew (m.translation ()) * m.rotation ();
      //        X_subspace.template block <3,3> (Motion::ANGULAR, 0) = m.rotation ();
      //
      //        return (X_subspace * S_minimal).eval();
      
      SE3ActionReturnType result;
      result.template block <3,3> (Motion::ANGULAR, 0) = M.rotation () * S_minimal;
      for (int k = 0; k < 3; ++k)
        result.template middleRows<3> (Motion::LINEAR).col(k) =
        M.translation ().cross(result.template middleRows<3> (Motion::ANGULAR).col(k));
      
      return result;
    }
    
  }; // struct ConstraintZYXRotationalSubspace

  typedef JointSphericalZYXTpl<double,0> JointSphericalZYX;

  inline Motion operator^ (const Motion & m1, const JointSphericalZYX::MotionSpherical & m2)
  {
//    const Motion::Matrix3 m2_cross (skew (Motion::Vector3 (-m2.w)));
//    return Motion(m2_cross * m1.linear (), m2_cross * m1.angular ());
    return Motion(m1.linear ().cross (m2.w), m1.angular ().cross (m2.w));
  }

  /* [CRBA] ForceSet operator* (Inertia Y,Constraint S) */

  template <typename _Scalar, int _Options>
//  const typename Eigen::ProductReturnType<
//  Eigen::Matrix < _Scalar, 6, 3, _Options >,
//  Eigen::Matrix < _Scalar, 3, 3, _Options> 
//  >::Type
  Eigen::Matrix <_Scalar,6,3,_Options>
  operator* (const InertiaTpl<_Scalar,_Options> & Y,
             const typename JointSphericalZYXTpl<_Scalar,_Options>::ConstraintZYXRotationalSubspace & S)
  {
    Eigen::Matrix < _Scalar, 6, 3, _Options > M;
    M.template topRows<3>() = alphaSkew ( -Y.mass (),  Y.lever () );
    M.template bottomRows<3> () = 
      (typename InertiaTpl<_Scalar,_Options>::Matrix3)(Y.inertia () - 
       typename Symmetric3::AlphaSkewSquare(Y.mass (), Y.lever ()));

    return (M * S.matrix ()).eval();
  }
  
  /* [ABA] Y*S operator (Inertia Y,Constraint S) */
  //  inline Eigen::Matrix<double,6,3>
  template <typename _Scalar, int _Options>
  typename Eigen::ProductReturnType<
  const Eigen::Block<const Inertia::Matrix6,6,3>,
  const typename JointSphericalZYXTpl<_Scalar,_Options>::ConstraintZYXRotationalSubspace::Matrix3
  >::Type
  operator*(const typename InertiaTpl<_Scalar,_Options>::Matrix6 & Y,
            const typename JointSphericalZYXTpl<_Scalar,_Options>::ConstraintZYXRotationalSubspace & S)
  {
    return Y.template block<6,3> (0,Inertia::ANGULAR,0,3) * S.S_minimal;
  }
  

//  namespace internal
//  {
//    template<>
//    struct ActionReturn<JointSphericalZYX::ConstraintZYXRotationalSubspace >
//    {
////      typedef const typename Eigen::ProductReturnType<
////      Eigen::Matrix <double,6,3,0>,
////      Eigen::Matrix <double,3,3,0>
////      >::Type Type;
//      typedef Eigen::Matrix <double,6,3,0> Type;
//    };
//  }

  template<>
  struct traits<JointSphericalZYX>
  {
    enum {
      NQ = 3,
      NV = 3
    };
    typedef JointDataSphericalZYX JointData;
    typedef JointModelSphericalZYX JointModel;
    typedef JointSphericalZYX::ConstraintZYXRotationalSubspace Constraint_t;
    typedef SE3 Transformation_t;
    typedef JointSphericalZYX::MotionSpherical Motion_t;
    typedef JointSphericalZYX::MotionSpherical Bias_t;
    typedef Eigen::Matrix<double,6,NV> F_t;
    
    // [ABA]
    typedef Eigen::Matrix<double,6,NV> U_t;
    typedef Eigen::Matrix<double,NV,NV> D_t;
    typedef Eigen::Matrix<double,6,NV> UD_t;

    typedef Eigen::Matrix<double,NQ,1> ConfigVector_t;
    typedef Eigen::Matrix<double,NV,1> TangentVector_t;
  };
  template<> struct traits<JointDataSphericalZYX> { typedef JointSphericalZYX Joint; };
  template<> struct traits<JointModelSphericalZYX> { typedef JointSphericalZYX Joint; };

  struct JointDataSphericalZYX : public JointDataBase<JointDataSphericalZYX>
  {
    typedef JointSphericalZYX Joint;
    SE3_JOINT_TYPEDEF;

    typedef Motion::Scalar_t Scalar;

    typedef Eigen::Matrix<Scalar,6,6> Matrix6;
    typedef Eigen::Matrix<Scalar,3,3> Matrix3;
    typedef Eigen::Matrix<Scalar,3,1> Vector3;
    
    Constraint_t S;
    Transformation_t M;
    Motion_t v;
    Bias_t c;

    F_t F;
    
    // [ABA] specific data
    U_t U;
    D_t Dinv;
    UD_t UDinv;

    JointDataSphericalZYX () : M(1), U(), Dinv(), UDinv() {}
    
    JointDataDense<NQ, NV> toDense_impl() const
    {
      return JointDataDense<NQ, NV>(S, M, v, c, F, U, Dinv, UDinv);
    }
  }; // strcut JointDataSphericalZYX

  struct JointModelSphericalZYX : public JointModelBase<JointModelSphericalZYX>
  {
    typedef JointSphericalZYX Joint;
    SE3_JOINT_TYPEDEF;

    using JointModelBase<JointModelSphericalZYX>::id;
    using JointModelBase<JointModelSphericalZYX>::idx_q;
    using JointModelBase<JointModelSphericalZYX>::idx_v;
    using JointModelBase<JointModelSphericalZYX>::setIndexes;
    typedef Motion::Vector3 Vector3;
    typedef double Scalar_t;

    JointData createData() const { return JointData(); }

    void calc (JointData & data,
               const Eigen::VectorXd & qs) const
    {
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q = qs.segment<NQ>(idx_q ());

      double c0,s0; SINCOS (q(0), &s0, &c0);
      double c1,s1; SINCOS (q(1), &s1, &c1);
      double c2,s2; SINCOS (q(2), &s2, &c2);

      data.M.rotation () << c0 * c1,
                c0 * s1 * s2 - s0 * c2,
                c0 * s1 * c2 + s0 * s2,
                s0 * c1,
                s0 * s1 * s2 + c0 * c2,
                s0 * s1 * c2 - c0 * s2,
                -s1,
                c1 * s2,
                c1 * c2;

      data.S.matrix () <<  -s1, 0., 1., c1 * s2, c2, 0, c1 * c2, -s2, 0;
    }

    void calc (JointData & data,
               const Eigen::VectorXd & qs,
               const Eigen::VectorXd & vs ) const
    {
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q = qs.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NV>::Type & q_dot = vs.segment<NQ> (idx_v ());

      double c0,s0; SINCOS (q(0), &s0, &c0);
      double c1,s1; SINCOS (q(1), &s1, &c1);
      double c2,s2; SINCOS (q(2), &s2, &c2);

      data.M.rotation () << c0 * c1,
                c0 * s1 * s2 - s0 * c2,
                c0 * s1 * c2 + s0 * s2,
                s0 * c1,
                s0 * s1 * s2 + c0 * c2,
                s0 * s1 * c2 - c0 * s2,
                -s1,
                c1 * s2,
                c1 * c2;


      data.S.matrix () <<  -s1, 0., 1., c1 * s2, c2, 0, c1 * c2, -s2, 0;

      data.v () = data.S.matrix () * q_dot;

      data.c ()(0) = -c1 * q_dot (0) * q_dot (1);
      data.c ()(1) = -s1 * s2 * q_dot (0) * q_dot (1) + c1 * c2 * q_dot (0) * q_dot (2) - s2 * q_dot (1) * q_dot (2);
      data.c ()(2) = -s1 * c2 * q_dot (0) * q_dot (1) - c1 * s2 * q_dot (0) * q_dot (2) - c2 * q_dot (1) * q_dot (2);
    }
    
    void calc_aba(JointData & data, Inertia::Matrix6 & I, const bool update_I) const
    {
      data.U.noalias() = I.middleCols<3> (Inertia::ANGULAR) * data.S.matrix();
      Inertia::Matrix3 tmp (data.S.matrix().transpose() * data.U.middleRows<3> (Inertia::ANGULAR));
      data.Dinv = tmp.inverse();
      data.UDinv.noalias() = data.U * data.Dinv;
      
      if (update_I)
        I.noalias() -= data.UDinv * data.U.transpose();
    }

    ConfigVector_t integrate_impl(const Eigen::VectorXd & qs,const Eigen::VectorXd & vs) const
    {
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q = qs.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NV>::Type & q_dot = vs.segment<NV> (idx_v ());


      return(q + q_dot);
    }

    ConfigVector_t interpolate_impl(const Eigen::VectorXd & q0,const Eigen::VectorXd & q1, const double u) const
    { 
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_0 = q0.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_1 = q1.segment<NQ> (idx_q ());

      return ((1-u) * q_0 + u * q_1);
    }

    ConfigVector_t random_impl() const
    { 
      ConfigVector_t result(ConfigVector_t::Random());
      return result;
    } 

    ConfigVector_t randomConfiguration_impl(const ConfigVector_t & lower_pos_limit, const ConfigVector_t & upper_pos_limit ) const throw (std::runtime_error)
    { 
      ConfigVector_t result;
      for (int i = 0; i < result.size(); ++i)
      {
        if(lower_pos_limit[i] == -std::numeric_limits<double>::infinity() || 
            upper_pos_limit[i] == std::numeric_limits<double>::infinity() )
        {
          std::ostringstream error;
          error << "non bounded limit. Cannot uniformly sample joint nb " << id() ;
          assert(false && "non bounded limit. Cannot uniformly sample joint spherical ZYX");
          throw std::runtime_error(error.str());
        }
        result[i] = lower_pos_limit[i] + ( upper_pos_limit[i] - lower_pos_limit[i]) * rand()/RAND_MAX;
      }
      return result;
    } 

    TangentVector_t difference_impl(const Eigen::VectorXd & q0,const Eigen::VectorXd & q1) const
    { 
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_0 = q0.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_1 = q1.segment<NQ> (idx_q ());

      return ( q_1 - q_0);

    } 

    double distance_impl(const Eigen::VectorXd & q0,const Eigen::VectorXd & q1) const
    { 
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_0 = q0.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_1 = q1.segment<NQ> (idx_q ());

      return (q_1 - q_0).norm();
    }
    
    JointModelDense<NQ, NV> toDense_impl() const
    {
      return JointModelDense<NQ, NV>( id(),
                                      idx_q(),
                                      idx_v()
                                    );
    }

    static const std::string shortname()
    {
      return std::string("JointModelSphericalZYX");
    }

    template <class D>
    bool operator == (const JointModelBase<D> &) const
    {
      return false;
    }
    
    bool operator == (const JointModelBase<JointModelSphericalZYX> & jmodel) const
    {
      return jmodel.id() == id()
              && jmodel.idx_q() == idx_q()
              && jmodel.idx_v() == idx_v();
    }

  }; // struct JointModelSphericalZYX

} // namespace se3

#endif // ifndef __se3_joint_spherical_ZYX_hpp__
