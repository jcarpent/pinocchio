//
// Copyright (c) 2015 CNRS
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

#ifndef __se3_joint_free_flyer_hpp__
#define __se3_joint_free_flyer_hpp__

#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/multibody/joint/joint-base.hpp"
#include "pinocchio/multibody/constraint.hpp"
#include "pinocchio/spatial/explog.hpp"
#include "pinocchio/math/fwd.hpp"
#include "pinocchio/math/quaternion.hpp"

#include <stdexcept>

namespace se3
{

  struct JointDataFreeFlyer;
  struct JointModelFreeFlyer;

  struct ConstraintIdentity;

  template <>
  struct traits<ConstraintIdentity>
  {
    typedef double Scalar_t;
    
    typedef Eigen::Matrix<double,3,1,0> Vector3;
    typedef Eigen::Matrix<double,4,1,0> Vector4;
    typedef Eigen::Matrix<double,6,1,0> Vector6;
    typedef Eigen::Matrix<double,3,3,0> Matrix3;
    typedef Eigen::Matrix<double,4,4,0> Matrix4;
    typedef Eigen::Matrix<double,6,6,0> Matrix6;
    typedef Matrix6 SE3ActionReturnType;
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
    typedef Eigen::Matrix<Scalar_t,6,1,0> JointMotion;
    typedef Eigen::Matrix<Scalar_t,6,1,0> JointForce;
    typedef Eigen::Matrix<Scalar_t,6,6> DenseBase;
  }; // traits ConstraintRevolute

//  namespace internal
//  {
//    template<>
//    struct ActionReturn<ConstraintIdentity>
//    { typedef SE3::Matrix6 Type; };
//  }

    struct ConstraintIdentity : ConstraintBase <ConstraintIdentity>
    {
      SPATIAL_TYPEDEF_NO_TEMPLATE(ConstraintIdentity);
      enum { NV = 6, Options = 0 };
      typedef traits<ConstraintIdentity>::JointMotion JointMotion;
      typedef traits<ConstraintIdentity>::JointForce JointForce;
      typedef traits<ConstraintIdentity>::DenseBase DenseBase;
      typedef traits<ConstraintIdentity>::SE3ActionReturnType SE3ActionReturnType;

      template<typename SE3Scalar, int SE3Options>
      SE3ActionReturnType SE3ActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
      {
        return M.toActionMatrix();
      }

      int nv_impl() const { return NV; }

      struct TransposeConst 
      {
        Force::Vector6 operator* (const Force & phi)
        {  return phi.coeffs();  }
      };
      
      TransposeConst transpose() const { return TransposeConst(); }
      operator ConstraintXd () const { return ConstraintXd(SE3::Matrix6::Identity()); }
    }; // struct ConstraintIdentity

    template<typename D>
    Motion operator* (const ConstraintIdentity&, const Eigen::MatrixBase<D>& v)
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(D,6);
      return Motion(v);
    }


  /* [CRBA] ForceSet operator* (Inertia Y,Constraint S) */
  inline Inertia::Matrix6 operator* (const Inertia & Y, const ConstraintIdentity &)
  {
    return Y.matrix();
  }
  
  /* [ABA] Y*S operator*/
  inline const Inertia::Matrix6 & operator* (const Inertia::Matrix6 & Y, const ConstraintIdentity &)
  {
    return Y;
  }
  
  inline Inertia::Matrix6 & operator* (Inertia::Matrix6 & Y, const ConstraintIdentity &)
  {
    return Y;
  }
  
  inline Inertia::Matrix6 operator* (const ConstraintIdentity::TransposeConst &, const Inertia & Y)
  {
    return Y.matrix();
  }

  /* [CRBA]  MatrixBase operator* (Constraint::Transpose S, ForceSet::Block) */
  template<typename EigenDerived>
  inline const Eigen::MatrixBase<EigenDerived> &
  operator*(const ConstraintIdentity::TransposeConst &, const Eigen::MatrixBase<EigenDerived> & F )
  {
    EIGEN_STATIC_ASSERT(EigenDerived::RowsAtCompileTime==6,THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE)
    return F;
  }

  
  

  struct JointFreeFlyer;

  template<>
  struct traits<JointFreeFlyer>
  {
    enum {
      NQ = 7,
      NV = 6
    };
    typedef JointDataFreeFlyer JointData;
    typedef JointModelFreeFlyer JointModel;
    typedef ConstraintIdentity Constraint_t;
    typedef SE3 Transformation_t;
    typedef Motion Motion_t;
    typedef BiasZero Bias_t;
    typedef Eigen::Matrix<double,6,NV> F_t;
    
    // [ABA]
    typedef Eigen::Matrix<double,6,NV> U_t;
    typedef Eigen::Matrix<double,NV,NV> D_t;
    typedef Eigen::Matrix<double,6,NV> UD_t;

    typedef Eigen::Matrix<double,NQ,1> ConfigVector_t;
    typedef Eigen::Matrix<double,NV,1> TangentVector_t;
  };
  template<> struct traits<JointDataFreeFlyer> { typedef JointFreeFlyer Joint; };
  template<> struct traits<JointModelFreeFlyer> { typedef JointFreeFlyer Joint; };

  struct JointDataFreeFlyer : public JointDataBase<JointDataFreeFlyer>
  {
    typedef JointFreeFlyer Joint;
    SE3_JOINT_TYPEDEF;

    typedef Eigen::Matrix<double,6,6> Matrix6;
    typedef Eigen::Quaternion<double> Quaternion;
    typedef Eigen::Matrix<double,3,1> Vector3;
    
    Constraint_t S;
    Transformation_t M;
    Motion_t v;
    Bias_t c;

    F_t F; // TODO if not used anymore, clean F_t
    
    // [ABA] specific data
    U_t U;
    D_t Dinv;
    UD_t UDinv;
    
    JointDataFreeFlyer() : M(1), U(), Dinv(), UDinv(UD_t::Identity()) {}

    JointDataDense<NQ, NV> toDense_impl() const
    {
      return JointDataDense<NQ, NV>(S, M, v, c, F, U, Dinv, UDinv);
    }
  }; // struct JointDataFreeFlyer

  struct JointModelFreeFlyer : public JointModelBase<JointModelFreeFlyer>
  {
    typedef JointFreeFlyer Joint;
    SE3_JOINT_TYPEDEF;

    using JointModelBase<JointModelFreeFlyer>::id;
    using JointModelBase<JointModelFreeFlyer>::idx_q;
    using JointModelBase<JointModelFreeFlyer>::idx_v;
    using JointModelBase<JointModelFreeFlyer>::setIndexes;
    typedef Motion::Vector3 Vector3;
    typedef double Scalar_t;

    JointData createData() const { return JointData(); }
    void calc(JointData & data,
              const Eigen::VectorXd & qs) const
    {
      typedef Eigen::Map<const Motion_t::Quaternion_t> ConstQuaternionMap_t;
      
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type q = qs.segment<NQ>(idx_q());
      ConstQuaternionMap_t quat(q.tail<4> ().data());

      data.M.rotation (quat.matrix());
      data.M.translation (q.head<3>());
    }
    
    void calc(JointData & data,
              const Eigen::VectorXd & qs,
              const Eigen::VectorXd & vs ) const
    {
      typedef Eigen::Map<const Motion_t::Quaternion_t> ConstQuaternionMap_t;
      
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type q = qs.segment<NQ>(idx_q());
      data.v = vs.segment<NV>(idx_v());

      ConstQuaternionMap_t quat(q.tail<4> ().data());
      
      data.M.rotation (quat.matrix());
      data.M.translation (q.head<3>());
    }
    
    void calc_aba(JointData & data, Inertia::Matrix6 & I, const bool update_I) const
    {
      data.U = I;
      data.Dinv = I.inverse();
      
      if (update_I)
        I.setZero();
    }

    ConfigVector_t integrate_impl(const Eigen::VectorXd & qs, const Eigen::VectorXd & vs) const
    {
      typedef Eigen::Map<Motion_t::Quaternion_t> QuaternionMap_t;
      typedef Eigen::Map<const Motion_t::Quaternion_t> ConstQuaternionMap_t;
      
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q = qs.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NV>::Type & q_dot = vs.segment<NV> (idx_v ());

      ConfigVector_t result;
      // Translational part
      result.head<3>() =  (q.head<3>() + q_dot.head<3>());

      // Quaternion part
      ConstQuaternionMap_t quat(q.tail<4>().data());
      Motion_t::Quaternion_t pOmega(se3::exp3(q_dot.tail<3>()));
      
      QuaternionMap_t quat_result (result.tail<4>().data());
      
      quat_result = pOmega*quat;

      return result; 
    } 

    ConfigVector_t interpolate_impl(const Eigen::VectorXd & q0, const Eigen::VectorXd & q1, const double u) const
    {
      typedef Eigen::Map<Motion_t::Quaternion_t> QuaternionMap_t;
      typedef Eigen::Map<const Motion_t::Quaternion_t> ConstQuaternionMap_t;
      
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_0 = q0.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_1 = q1.segment<NQ> (idx_q ());

      ConfigVector_t result;
      // Translational part
      result.head<3> () << ((1-u)*q_0.head<3>() + u * q_1.head<3>());

      //Quaternion part
      ConstQuaternionMap_t p0 (q_0.segment<4>(3).data());
      ConstQuaternionMap_t p1 (q_1.segment<4>(3).data());
      QuaternionMap_t quat_result (result.tail<4>().data());
      
      quat_result = p0.slerp(u, p1);

      return result; 
    }

    ConfigVector_t random_impl() const
    { 
      ConfigVector_t q(ConfigVector_t::Random());
      q.segment<4>(3).normalize();// /= q.segment<4>(3).norm();
      return q;
    } 

    ConfigVector_t randomConfiguration_impl(const ConfigVector_t & lower_pos_limit, const ConfigVector_t & upper_pos_limit ) const throw (std::runtime_error)
    {
      ConfigVector_t result;
      // Translational part
      for (Index i = 0; i < 3; ++i)
      {
        if(lower_pos_limit[i] == -std::numeric_limits<double>::infinity() ||
           upper_pos_limit[i] == std::numeric_limits<double>::infinity() )
        {
          std::ostringstream error;
          error << "non bounded limit. Cannot uniformly sample joint nb " << id();
          throw std::runtime_error(error.str());
        }
          
        result[i] = lower_pos_limit[i] + (upper_pos_limit[i] - lower_pos_limit[i]) * (Scalar_t)(rand())/RAND_MAX;
      }
          
      // Rotational part
      const Scalar_t u1 = (Scalar_t)rand() / RAND_MAX;
      const Scalar_t u2 = (Scalar_t)rand() / RAND_MAX;
      const Scalar_t u3 = (Scalar_t)rand() / RAND_MAX;
      
      const Scalar_t mult1 = sqrt (1-u1);
      const Scalar_t mult2 = sqrt (u1);
      
      Scalar_t s2,c2; SINCOS(2.*PI*u2,&s2,&c2);
      Scalar_t s3,c3; SINCOS(2.*PI*u3,&s3,&c3);
      
      
      result.segment<4>(3) << mult1 * s2,
                              mult1 * c2,
                              mult2 * s3,
                              mult2 * c3;
      return result;
    }

    TangentVector_t difference_impl(const Eigen::VectorXd & q0, const Eigen::VectorXd & q1) const
    {
      typedef Eigen::Map<const Motion_t::Quaternion_t> ConstQuaternionMap_t;
      using std::acos;
      
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_0 = q0.segment<NQ> (idx_q ());
      Eigen::VectorXd::ConstFixedSegmentReturnType<NQ>::Type & q_1 = q1.segment<NQ> (idx_q ());

      TangentVector_t result;
      // Translational part
      result.head<3>() << q_1.head<3> () - q_0.head<3> ();

      // Quaternion part
      // Compute relative rotation between q0 and q1.
      ConstQuaternionMap_t quat0 (q_0.segment<4>(3).data());
      ConstQuaternionMap_t quat1 (q_1.segment<4>(3).data());
      
      const Motion_t::Quaternion_t quat_relatif (quat1*quat0.conjugate());
      const Scalar_t theta = 2.*acos(quat_relatif.w());
      
      if (quat0.dot(quat1) >= 0.)
        result.tail<3>() << theta * quat_relatif.vec().normalized();
      else
        result.tail<3>() << -(2*PI-theta) * quat_relatif.vec().normalized();

      return result;
    } 

    double distance_impl(const Eigen::VectorXd & q0,const Eigen::VectorXd & q1) const
    { 
      return difference_impl(q0,q1).norm();
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
      return std::string("JointModelFreeFlyer");
    }

    template <class D>
    bool operator == (const JointModelBase<D> &) const
    {
      return false;
    }
    
    bool operator == (const JointModelBase<JointModelFreeFlyer> & jmodel) const
    {
      return jmodel.id() == id()
              && jmodel.idx_q() == idx_q()
              && jmodel.idx_v() == idx_v();
    }
  }; // struct JointModelFreeFlyer

} // namespace se3

#endif // ifndef __se3_joint_free_flyer_hpp__
