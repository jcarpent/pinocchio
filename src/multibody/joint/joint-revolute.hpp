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

#ifndef __se3_joint_revolute_hpp__
#define __se3_joint_revolute_hpp__

#include "pinocchio/math/sincos.hpp"
#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/multibody/joint/joint-base.hpp"
#include "pinocchio/multibody/joint/joint-dense.hpp"

#include <stdexcept>

namespace se3
{

  template<int axis> struct JointDataRevolute;
  template<int axis> struct JointModelRevolute;
  
  template<int axis> struct SE3Revolute;
  template<int axis> struct MotionRevolute;

  namespace revolute
  {
    template<int axis>
    struct CartesianVector3
    {
      double w; 
      CartesianVector3(const double w) : w(w) {}
      CartesianVector3() : w(NAN) {}
      operator Eigen::Vector3d ();
    };
    template<> inline CartesianVector3<0>::operator Eigen::Vector3d () { return Eigen::Vector3d(w,0,0); }
    template<> inline CartesianVector3<1>::operator Eigen::Vector3d () { return Eigen::Vector3d(0,w,0); }
    template<> inline CartesianVector3<2>::operator Eigen::Vector3d () { return Eigen::Vector3d(0,0,w); }
    inline Eigen::Vector3d operator+ (const Eigen::Vector3d & w1,const CartesianVector3<0> & wx)
    { return Eigen::Vector3d(w1[0]+wx.w,w1[1],w1[2]); }
    inline Eigen::Vector3d operator+ (const Eigen::Vector3d & w1,const CartesianVector3<1> & wy)
    { return Eigen::Vector3d(w1[0],w1[1]+wy.w,w1[2]); }
    inline Eigen::Vector3d operator+ (const Eigen::Vector3d & w1,const CartesianVector3<2> & wz)
    { return Eigen::Vector3d(w1[0],w1[1],w1[2]+wz.w); }
  } // namespace revolute


  template<int axis>
  struct traits< MotionRevolute<axis> > : traits< MotionBase< MotionRevolute<axis> > >
  {
    typedef double Scalar_t;
    typedef MotionRevolute<axis> Type;
    typedef traits< MotionBase<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
    typedef Scalar_t DotReturnType;
    typedef Eigen::Matrix<double,3,1,0> Vector3;
    typedef Eigen::Matrix<double,4,1,0> Vector4;
    typedef Eigen::Matrix<double,6,1,0> Vector6;
    typedef Eigen::Matrix<double,3,3,0> Matrix3;
    typedef Eigen::Matrix<double,4,4,0> Matrix4;
    typedef Eigen::Matrix<double,6,6,0> Matrix6;
    typedef Vector3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Vector3 ConstAngular_t;
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
  }; // traits MotionRevolute
  
  template<int axis> struct SE3Revolute;
  
  template<int axis>
  struct traits< SE3Revolute<axis> > : traits< SE3Base < SE3Revolute<axis> > >
  {
    typedef SE3Revolute<axis> Type;
    typedef traits< SE3Base<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
    
    typedef double Scalar_t;
    typedef Eigen::Matrix<Scalar_t,3,1,0> SE3ActionOnVecReturnType;
  };
  
  
  /// \brief SE3Revolute is a sparse implementation of SE3 transform dedicated
  ///        to the transformation allowed by a Revolute Joint.
  ///        Question: is this class must also heritate from Eigen::AngleAxis?
  template<int _axis>
  struct SE3Revolute : SE3Base < SE3Revolute<_axis> >
  {
    typedef double Scalar;
    enum { axis = _axis };
    
    SE3Revolute(const Scalar & angle) : angle_(angle) {}
    SE3Revolute() : sa_(NAN), ca_(NAN) {}
    
//    Scalar & angle() { return angle_; }
//    const Scalar & angle() const { return angle_; }
    
    Scalar & sa() { return sa_; }
    Scalar & ca() { return ca_; }
    
    const Scalar & sa() const { return sa_; }
    const Scalar & ca() const { return ca_; }
    
   
//    template<typename OtherScalar, int OtherOptions>
//    inline void setTo(SE3Tpl<OtherScalar,OtherOptions> & M) const;
//    {
//      M.translation().setZero();
//      M.rotation() = rotation();
//    }
    
    template<typename OtherScalar, int OtherOptions>
    operator SE3Tpl<OtherScalar,OtherOptions> () const
    {
      typedef SE3Tpl<OtherScalar,OtherOptions> ReturnType;
      ReturnType res;
      
//      Scalar ca,sa;
//      SINCOS(angle_,&sa,&ca);
      
      switch (axis) {
        case 0:
          res.rotation() <<
          1,0,0,
          0,ca_,-sa_,
          0,sa_,ca_;
          break;
        case 1:
          res.rotation() <<
          ca_, 0,  sa_,
          0, 1,   0,
          -sa_, 0,  ca_;
          break;
        case 2:
          res.rotation() <<
          ca_,-sa_,0,
          sa_,ca_,0,
          0,0,1;
          break;
      }
      
      return res;
    }
    
    SE3Tpl<Scalar> dense() const
    {
      typedef SE3Tpl<Scalar> ReturnType;
//      ReturnType res;
//      
//      switch (axis) {
//        case 0:
//          res.rotation() <<
//          1,0,0,
//          0,ca_,-sa_,
//          0,sa_,ca_;
//          break;
//        case 1:
//          res.rotation() <<
//          ca_, 0,  sa_,
//          0, 1,   0,
//          -sa_, 0,  ca_;
//          break;
//        case 2:
//          res.rotation() <<
//          ca_,-sa_,0,
//          sa_,ca_,0,
//          0,0,1;
//          break;
//      }
      
      return (ReturnType) (*this);
    }
    
  protected:
    
    /// \brief Fixed axis of the angle axis transform.
    ///        Value must be 0 (X axis), 1 (Y axis) or 2 (Z axis)
    enum { AXIS = _axis };
    Scalar angle_;
    Scalar sa_, ca_; // TODO remove this quantity
  };

  template<int axis>
  struct MotionRevolute : MotionBase < MotionRevolute<axis> >
  {
    SPATIAL_TYPEDEF_TEMPLATE(MotionRevolute);

    MotionRevolute()                   : w(NAN) {}
    MotionRevolute( const double & w ) : w(w)  {}
    double w;

    operator Motion() const
    {
      return Motion(  Motion::Vector3::Zero(),
        (Vector3)typename revolute::CartesianVector3<axis>(w)
        );
    }
    
    template<typename OtherScalar, int OtherOptions>
    void addTo (MotionTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.angular()[axis] += w;
    }
    
    template<typename OtherScalar, int OtherOptions>
    void subTo (MotionTpl<OtherScalar, OtherOptions> & dest) const
    {
      dest.angular()[axis] -= w;
    }
    
    template<typename OtherScalar, int OtherOptions>
    MotionTpl<OtherScalar, OtherOptions> sub (const MotionTpl<OtherScalar, OtherOptions> & other) const
    {
      MotionTpl<OtherScalar, OtherOptions> res (-other.coeffs());
      res.angular()[axis] += w;
      
      return res;
    }
  }; // struct MotionRevolute

  template<int axis, typename MotionScalar, int MotionOptions>
  MotionTpl<MotionScalar,MotionOptions> operator+ (const MotionRevolute<axis>& m1, const MotionTpl<MotionScalar,MotionOptions> & m2)
  {
    typedef MotionTpl<MotionScalar,MotionOptions> ReturnType;
    ReturnType res(m2);
    res += m1;
    
    return res;
    //    return ReturnType(m2.linear(),m2.angular()+typename revolute::CartesianVector3<axis>(m1.w));
  }
  
  template<int axis, typename MotionScalar, int MotionOptions>
  MotionTpl<MotionScalar,MotionOptions> operator- (const MotionRevolute<axis>& m1, const MotionTpl<MotionScalar,MotionOptions> & m2)
  {
    typedef MotionTpl<MotionScalar,MotionOptions> ReturnType;
    ReturnType res(-m2);
    res += m1;
    
    return res;
    //    return ReturnType(-m2.linear(),-m2.angular()+typename revolute::CartesianVector3<axis>(m1.w));
  }

  template<int axis> struct ConstraintRevolute;

  template<int axis>
  struct traits< ConstraintRevolute<axis> >
  {
    typedef double Scalar_t;
    typedef Eigen::Matrix<double,6,1,0> SE3ActionReturnType;
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
    typedef Eigen::Matrix<Scalar_t,1,1,0> JointMotion;
    typedef Eigen::Matrix<Scalar_t,1,1,0> JointForce;
    typedef Eigen::Matrix<Scalar_t,6,1> DenseBase;
  }; // traits ConstraintRevolute

  template<int axis>
  struct ConstraintRevolute : ConstraintBase < ConstraintRevolute<axis> >
  { 
    SPATIAL_TYPEDEF_TEMPLATE(ConstraintRevolute);
    enum { NV = 1, Options = 0 };
    typedef typename traits<ConstraintRevolute>::JointMotion JointMotion;
    typedef typename traits<ConstraintRevolute>::JointForce JointForce;
    typedef typename traits<ConstraintRevolute>::DenseBase DenseBase;
    typedef typename traits<ConstraintRevolute>::SE3ActionReturnType SE3ActionReturnType;

    template<typename EigenDerived>
    MotionRevolute<axis> operator*(const Eigen::MatrixBase<EigenDerived> & v) const
    {
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(EigenDerived,NV);
      return MotionRevolute<axis>(v[0]);
    }

    template<typename SE3Scalar, int SE3Options>
    SE3ActionReturnType SE3ActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    { 
      SE3ActionReturnType res;
      res.template head<3>() = M.translation().cross(M.rotation().col(axis));
      res.template tail<3>() = M.rotation().col(axis);
      return res;
    }

    int nv_impl() const { return NV; }
    
    struct TransposeConst
    {
      const ConstraintRevolute<axis> & ref; 
      TransposeConst(const ConstraintRevolute<axis> & ref) : ref(ref) {} 

      typename Force::ConstAngular_t::template ConstFixedSegmentReturnType<1>::Type
      operator* (const Force & f) const
      { return f.angular().template segment<1>(axis); }

        /// [CRBA]  MatrixBase operator* (Constraint::Transpose S, ForceSet::Block)
      template<typename EigenDerived>
      friend typename Eigen::MatrixBase<EigenDerived>::ConstRowXpr
      operator*(const TransposeConst &, const Eigen::MatrixBase<EigenDerived> & F )
      {
        EIGEN_STATIC_ASSERT(EigenDerived::RowsAtCompileTime==6,THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);
        return F.row(Inertia::ANGULAR + axis);
      }
    }; // struct TransposeConst

    TransposeConst transpose() const { return TransposeConst(*this); }


    /* CRBA joint operators
     *   - ForceSet::Block = ForceSet
     *   - ForceSet operator* (Inertia Y,Constraint S)
     *   - MatrixBase operator* (Constraint::Transpose S, ForceSet::Block)
     *   - SE3::act(ForceSet::Block)
     */
     operator ConstraintXd () const
     {
      Eigen::Matrix<double,6,1> S;
      S << Eigen::Vector3d::Zero(), (Eigen::Vector3d)revolute::CartesianVector3<axis>(1);
      return ConstraintXd(S);
    }
  }; // struct ConstraintRevolute

  template<int axis> 
  struct JointRevolute {
      static Eigen::Matrix3d cartesianRotation(const double & angle); 
  };

  inline Motion operator^( const Motion& m1, const MotionRevolute<0>& m2)
  {
    /* nu1^nu2    = ( v1^w2+w1^v2, w1^w2 )
     * nu1^(0,w2) = ( v1^w2      , w1^w2 )
     * (x,y,z)^(w,0,0) = ( 0,zw,-yw )
     * nu1^(0,wx) = ( 0,vz1 wx,-vy1 wx,    0,wz1 wx,-wy1 wx)
     */
    const Motion::Vector3& v = m1.linear();
    const Motion::Vector3& w = m1.angular();
    const double & wx = m2.w;
    return Motion( Motion::Vector3(0,v[2]*wx,-v[1]*wx),
                   Motion::Vector3(0,w[2]*wx,-w[1]*wx)
                 );
  }

  inline Motion operator^( const Motion& m1, const MotionRevolute<1>& m2)
  {
    /* nu1^nu2    = ( v1^w2+w1^v2, w1^w2 )
     * nu1^(0,w2) = ( v1^w2      , w1^w2 )
     * (x,y,z)^(0,w,0) = ( -z,0,x )
     * nu1^(0,wx) = ( -vz1 wx,0,vx1 wx,    -wz1 wx,0,wx1 wx)
     */
    const Motion::Vector3& v = m1.linear();
    const Motion::Vector3& w = m1.angular();
    const double & wx = m2.w;
    return Motion(  Motion::Vector3(-v[2]*wx,0, v[0]*wx),
                    Motion::Vector3(-w[2]*wx,0, w[0]*wx)
                 );
  }

  inline Motion operator^( const Motion& m1, const MotionRevolute<2>& m2)
  {
    /* nu1^nu2    = ( v1^w2+w1^v2, w1^w2 )
     * nu1^(0,w2) = ( v1^w2      , w1^w2 )
     * (x,y,z)^(0,0,w) = ( y,-x,0 )
     * nu1^(0,wx) = ( vy1 wx,-vx1 wx,0,    wy1 wx,-wx1 wx,0 )
     */
    const Motion::Vector3& v = m1.linear();
    const Motion::Vector3& w = m1.angular();
    const double & wx = m2.w;
    return Motion( Motion::Vector3(v[1]*wx,-v[0]*wx,0),
                   Motion::Vector3(w[1]*wx,-w[0]*wx,0)
                 );
    }

  template<> inline
  Eigen::Matrix3d JointRevolute<0>::cartesianRotation(const double & angle) 
  {
    Eigen::Matrix3d R3; 
    double ca,sa; SINCOS (angle,&sa,&ca);
    R3 << 
    1,0,0,
    0,ca,-sa,
    0,sa,ca;
    return R3;
  }

  template<> inline
  Eigen::Matrix3d JointRevolute<1>::cartesianRotation(const double & angle)
  {
    Eigen::Matrix3d R3; 
    double ca,sa; SINCOS (angle,&sa,&ca);
    R3 << 
    ca, 0,  sa,
    0, 1,   0,
    -sa, 0,  ca;
    return R3;
  }

  template<> inline
  Eigen::Matrix3d JointRevolute<2>::cartesianRotation(const double & angle) 
  {
    Eigen::Matrix3d R3; 
    double ca,sa; SINCOS (angle,&sa,&ca);
    R3 << 
    ca,-sa,0,
    sa,ca,0,
    0,0,1;
    return R3;
  }

  /* [CRBA] ForceSet operator* (Inertia Y,Constraint S) */
  Eigen::Matrix<double,6,1> inline
  operator*( const Inertia& Y,const ConstraintRevolute<0> & )
  { 
    /* Y(:,3) = ( 0,-z, y,  I00+yy+zz,  I01-xy   ,  I02-xz   ) */
    const double 
    &m = Y.mass(),
    &x = Y.lever()[0],
    &y = Y.lever()[1],
    &z = Y.lever()[2];
    const Inertia::Symmetric3 & I = Y.inertia();
    
    Eigen::Matrix<double,6,1> res; res << 0.0,-m*z,m*y,
    I(0,0)+m*(y*y+z*z),
    I(0,1)-m*x*y,
    I(0,2)-m*x*z ;
    return res;
  }
  /* [CRBA] ForceSet operator* (Inertia Y,Constraint S) */
  Eigen::Matrix<double,6,1> inline
  operator*( const Inertia& Y,const ConstraintRevolute<1> & )
  { 
    /* Y(:,4) = ( z, 0,-x,  I10-xy   ,  I11+xx+zz,  I12-yz   ) */
    const double 
    &m = Y.mass(),
    &x = Y.lever()[0],
    &y = Y.lever()[1],
    &z = Y.lever()[2];
    const Inertia::Symmetric3 & I = Y.inertia();
    Eigen::Matrix<double,6,1> res; res << m*z,0,-m*x,
    I(1,0)-m*x*y,
    I(1,1)+m*(x*x+z*z),
    I(1,2)-m*y*z ;
    return res;
  }
  /* [CRBA] ForceSet operator* (Inertia Y,Constraint S) */
  Eigen::Matrix<double,6,1> inline
  operator*( const Inertia& Y,const ConstraintRevolute<2> & )
  { 
    /* Y(:,5) = (-y, x, 0,  I20-xz   ,  I21-yz   ,  I22+xx+yy) */
    const double 
    &m = Y.mass(),
    &x = Y.lever()[0],
    &y = Y.lever()[1],
    &z = Y.lever()[2];
    const Inertia::Symmetric3 & I = Y.inertia();
    Eigen::Matrix<double,6,1> res; res << -m*y,m*x,0,
    I(2,0)-m*x*z,
    I(2,1)-m*y*z,
    I(2,2)+m*(x*x+y*y) ;
    return res;
  }
  
  /* [ABA] I*S operator (Inertia Y,Constraint S) */
  template<int axis>
  inline const Eigen::MatrixBase<const Inertia::Matrix6>::ColXpr
  operator*(const Inertia::Matrix6 & Y,const ConstraintRevolute<axis> & )
  {
    return Y.col(Inertia::ANGULAR + axis);
  }
  
  template<int axis>
  inline Eigen::MatrixBase<Inertia::Matrix6>::ColXpr
  operator*(Inertia::Matrix6 & Y,const ConstraintRevolute<axis> & )
  {
    return Y.col(Inertia::ANGULAR + axis);
  }

//  namespace internal 
//  {
//    template<int axis>
//    struct ActionReturn<ConstraintRevolute<axis> >
//    { typedef Eigen::Matrix<double,6,1> Type; };
//  }

  template<int axis>
  struct traits< JointRevolute<axis> >
  {
    enum {
      NQ = 1,
      NV = 1
    };
    
    typedef JointDataRevolute<axis> JointData;
    typedef JointModelRevolute<axis> JointModel;
    typedef ConstraintRevolute<axis> Constraint_t;
//    typedef SE3 Transformation_t;
    typedef SE3Revolute<axis> Transformation_t;
    typedef MotionRevolute<axis> Motion_t;
    typedef BiasZero Bias_t;
    typedef Eigen::Matrix<double,6,NV> F_t;
    
    // [ABA]
    typedef Eigen::Matrix<double,6,NV> U_t;
    typedef Eigen::Matrix<double,NV,NV> D_t;
    typedef Eigen::Matrix<double,6,NV> UD_t;

    typedef Eigen::Matrix<double,NQ,1> ConfigVector_t;
    typedef Eigen::Matrix<double,NV,1> TangentVector_t;
  };

  template<int axis> struct traits< JointDataRevolute<axis> > { typedef JointRevolute<axis> Joint; };
  template<int axis> struct traits< JointModelRevolute<axis> > { typedef JointRevolute<axis> Joint; };

  template<int axis>
  struct JointDataRevolute : public JointDataBase< JointDataRevolute<axis> >
  {
    typedef JointRevolute<axis> Joint;
    SE3_JOINT_TYPEDEF_TEMPLATE;

    Constraint_t S;
    Transformation_t M;
    Motion_t v;
    Bias_t c;
    F_t F;

    // [ABA] specific data
    U_t U;
    D_t Dinv;
    UD_t UDinv;

    JointDataRevolute() : M(1), U(), Dinv(), UDinv()
    {}

    JointDataDense<NQ, NV> toDense_impl() const
    {
      return JointDataDense<NQ, NV>(S, M, v, c, F, U, Dinv, UDinv);
    }
  }; // struct JointDataRevolute

  template<int axis>
  struct JointModelRevolute : public JointModelBase< JointModelRevolute<axis> >
  {
    typedef JointRevolute<axis> Joint;
    SE3_JOINT_TYPEDEF_TEMPLATE;

    using JointModelBase<JointModelRevolute>::id;
    using JointModelBase<JointModelRevolute>::idx_q;
    using JointModelBase<JointModelRevolute>::idx_v;
    using JointModelBase<JointModelRevolute>::setIndexes;
    typedef Motion::Vector3 Vector3;
    typedef double Scalar_t;
    
    JointData createData() const { return JointData(); }
    inline void calc(JointData & data,
                     const Eigen::VectorXd & qs ) const
    {
//      const double & q = qs[idx_q()];
//      data.M.rotation(JointRevolute<axis>::cartesianRotation(q));
//      data.M.angle() = qs[idx_q()];
      SINCOS (qs[idx_q()],&data.M.sa(),&data.M.ca());
    }

    inline void calc( JointData & data,
                     const Eigen::VectorXd & qs,
                     const Eigen::VectorXd & vs ) const
    {
//      const double & q = qs[idx_q()];
      const double & v = vs[idx_v()];

      
      SINCOS (qs[idx_q()],&data.M.sa(),&data.M.ca());
//      data.M.rotation(JointRevolute<axis>::cartesianRotation(q));
      data.v.w = v;
    }
    
    inline void calc_aba(JointData & data, Inertia::Matrix6 & I, const bool update_I) const
    {
      data.U = I.col(Inertia::ANGULAR + axis);
      data.Dinv[0] = 1./I(Inertia::ANGULAR + axis,Inertia::ANGULAR + axis);
      data.UDinv.noalias() = data.U * data.Dinv[0];
      
      if (update_I)
        I.noalias() -= data.UDinv * data.U.transpose();
    }

    ConfigVector_t integrate_impl(const Eigen::VectorXd & qs,const Eigen::VectorXd & vs) const
    {
      const Scalar_t & q = qs[idx_q()];
      const Scalar_t & v = vs[idx_v()];

      ConfigVector_t result;
      result << (q + v);
      return result; 
    }

    ConfigVector_t interpolate_impl(const Eigen::VectorXd & q0,const Eigen::VectorXd & q1, const double u) const
    { 
      const Scalar_t & q_0 = q0[idx_q()];
      const Scalar_t & q_1 = q1[idx_q()];

      ConfigVector_t result;
      result << ((1-u) * q_0 + u * q_1);
      return result; 
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
          assert(false && "non bounded limit. Cannot uniformly sample joint revolute");
          throw std::runtime_error(error.str());
        }
        result[i] = lower_pos_limit[i] + ( upper_pos_limit[i] - lower_pos_limit[i]) * rand()/RAND_MAX;
      }
      return result;
    } 

    TangentVector_t difference_impl(const Eigen::VectorXd & q0,const Eigen::VectorXd & q1) const
    { 
      const Scalar_t & q_0 = q0[idx_q()];
      const Scalar_t & q_1 = q1[idx_q()];

      ConfigVector_t result;
      result << (q_1 - q_0);
      return result; 
    } 

    double distance_impl(const Eigen::VectorXd & q0,const Eigen::VectorXd & q1) const
    { 
      const Scalar_t & q_0 = q0[idx_q()];
      const Scalar_t & q_1 = q1[idx_q()];

      return (q_1-q_0);
    }

    JointModelDense<NQ, NV> toDense_impl() const
    {
      return JointModelDense<NQ, NV>( id(),
                                      idx_q(),
                                      idx_v()
                                    );
    }

    static const std::string shortname();

    template <class D>
    bool operator == (const JointModelBase<D> &) const
    {
      return false;
    }
    
    bool operator == (const JointModelBase<JointModelRevolute <axis> > & jmodel) const
    {
      return jmodel.id() == id()
              && jmodel.idx_q() == idx_q()
              && jmodel.idx_v() == idx_v();
    }

  }; // struct JointModelRevolute

  typedef JointRevolute<0> JointRX;
  typedef JointDataRevolute<0> JointDataRX;
  typedef JointModelRevolute<0> JointModelRX;

  template<> inline
  const std::string JointModelRevolute<0>::shortname()
  {
    return std::string("JointModelRX") ;
  }

  typedef JointRevolute<1> JointRY;
  typedef JointDataRevolute<1> JointDataRY;
  typedef JointModelRevolute<1> JointModelRY;

  template<> inline
  const std::string JointModelRevolute<1>::shortname()
  {
    return std::string("JointModelRY") ;
  }

  typedef JointRevolute<2> JointRZ;
  typedef JointDataRevolute<2> JointDataRZ;
  typedef JointModelRevolute<2> JointModelRZ;

  template<> inline
  const std::string JointModelRevolute<2>::shortname()
  {
    return std::string("JointModelRZ") ;
  }

} //namespace se3

#endif // ifndef __se3_joint_revolute_hpp__
