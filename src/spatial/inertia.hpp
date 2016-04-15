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

#ifndef __se3_inertia_hpp__
#define __se3_inertia_hpp__

#include <Eigen/Core>
#include <iostream>

#include "pinocchio/spatial/symmetric3.hpp"
#include "pinocchio/spatial/spatial.hpp"
#include "pinocchio/spatial/force.hpp"
#include "pinocchio/spatial/motion.hpp"
#include "pinocchio/spatial/skew.hpp"

namespace se3
{
  template<class Derived>
  struct traits<InertiaBase <Derived> >
  {
    typedef Derived SE3ActionReturnType;
  };

  template<class Derived>
  class InertiaBase : public SpatialBase<Derived>
  {
  protected:

    SPATIAL_TYPEDEF_TEMPLATE(Derived);

  public:
    typedef SpatialBase<Derived> Base;
    using Base::derived;

//    Scalar_t           mass()    const { return static_cast<const Derived*>(this)->mass(); }
//    Scalar_t &         mass() { return static_cast<const Derived*>(this)->mass(); }
//    const Vector3 &    lever()   const { return static_cast<const Derived*>(this)->lever(); }
//    Vector3 &          lever() { return static_cast<const Derived*>(this)->lever(); }
//    const Symmetric3 & inertia() const { return static_cast<const Derived*>(this)->inertia(); }
//    Symmetric3 &       inertia() { return static_cast<const Derived*>(this)->inertia(); }

//    Matrix6 matrix() const { return derived().matrix_impl(); }
//    operator Matrix6 () const { return matrix(); }

    template<typename OtherDerived>
    inline bool operator== (const InertiaBase<OtherDerived> & other) const
    {
      return derived().isEqual(other.derived());
    }
    
    template<typename OtherDerived>
    inline Derived & operator= (const InertiaBase<OtherDerived> & other)
    {
      other.derived().setTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived & operator+= (const InertiaBase<OtherDerived> & other)
    {
      other.derived().addTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived & operator-= (const InertiaBase<OtherDerived> & other)
    {
      other.derived().subTo(derived());
      return derived();
    }
    
    template<typename OtherDerived>
    inline Derived operator+ (const InertiaBase<OtherDerived> & other) const
    {
      return derived().add(other.derived());
    }
    
    template<typename MotionDerived, typename ForceDerived>
    ForceBase<ForceDerived> operator* (const MotionBase<MotionDerived> & v) const
    {
      return derived().__mult__(v);
    }
    
//    Derived& operator= (const Derived& clone){return derived().__equl__(clone);}
//    bool operator== (const Derived& other) const {return derived().isEqual(other);}
//    Derived& operator+= (const Derived & Yb) { return derived().__pequ__(Yb); }
//    Derived operator+(const Derived & Yb) const { return derived().__plus__(Yb); }
    

    template<typename MotionDerived>
    Scalar_t vtiv(const MotionBase<MotionDerived> & v) const
    {
      return derived().vtiv(v.derived());
    }

    void setZero() { derived().setZero(); }
    void setIdentity() { derived().setIdentity(); }
    void setRandom() { derived().setRandom(); }

    /// aI = aXb.act(bI)
//    Derived se3Action(const SE3 & M) const { return derived().se3Action_impl(M); }

    /// bI = aXb.actInv(aI)
//    Derived se3ActionInverse(const SE3 & M) const { return derived().se3ActionInverse_impl(M); }

  }; // class InertiaBase


  template<typename T, int U>
  struct traits< InertiaTpl<T,U> > : traits< InertiaBase< InertiaTpl<T,U> > >
  {
    typedef T Scalar_t;
    typedef InertiaTpl<T,U> Type;
    typedef traits< InertiaBase<Type> > BaseTraits;
    using typename BaseTraits::SE3ActionReturnType;
    typedef Eigen::Matrix<T,3,1,U> Vector3;
    typedef Eigen::Matrix<T,4,1,U> Vector4;
    typedef Eigen::Matrix<T,6,1,U> Vector6;
    typedef Eigen::Matrix<T,3,3,U> Matrix3;
    typedef Eigen::Matrix<T,4,4,U> Matrix4;
    typedef Eigen::Matrix<T,6,6,U> Matrix6;
    typedef Matrix6 ActionMatrix_t;
    typedef Vector3 Angular_t;
    typedef Vector3 Linear_t;
    typedef const Vector3 ConstAngular_t;
    typedef const Vector3 ConstLinear_t;
    typedef Eigen::Quaternion<T,U> Quaternion_t;
    typedef SE3Tpl<T,U> SE3;
    typedef ForceTpl<T,U> Force;
    typedef MotionTpl<T,U> Motion;
    typedef Symmetric3Tpl<T,U> Symmetric3;
    enum {
      LINEAR = 0,
      ANGULAR = 3
    };
  }; // traits InertiaTpl

  template<typename _Scalar, int _Options>
  class InertiaTpl : public InertiaBase< InertiaTpl<_Scalar,_Options> >
  {
  public:
    friend class InertiaBase< InertiaTpl<_Scalar,_Options> >;
    typedef InertiaBase< InertiaTpl<_Scalar,_Options> > Base;
    typedef typename traits<InertiaTpl>::SE3ActionReturnType SE3ActionReturnType;
    SPATIAL_TYPEDEF_TEMPLATE(InertiaTpl);
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
  public:
    // Constructors
    InertiaTpl() : m(), c(), I() {}

    InertiaTpl(const Scalar_t m_, const Vector3 &c_, const Matrix3 &I_)
    : m(m_), c(c_), I(I_)
    {}
    
    template<typename EigenDerived>
    InertiaTpl(const Eigen::MatrixBase<EigenDerived> & I6)
    {
      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(EigenDerived,6,6);
      assert((I6 - I6.transpose()).isMuchSmallerThan(I6));
      m = I6(LINEAR, LINEAR);
      const Matrix3 & mc_cross = I6.template block <3,3> (ANGULAR,LINEAR);
      c = unSkew(mc_cross);
      c /= m;
      
      Matrix3 I3 (mc_cross * mc_cross);
      I3 /= m;
      I3 += I6.template block<3,3>(ANGULAR,ANGULAR);
      I = Symmetric3(I3);
    }

    InertiaTpl(Scalar_t _m, 
     const Vector3 &_c, 
     const Symmetric3 &_I)
    : m(_m),
    c(_c),
    I(_I)
    {

    }
    InertiaTpl(const InertiaTpl & clone)  // Clone constructor for std::vector 
    : m(clone.m),
    c(clone.c),
    I(clone.I)    
    {

    }

    template<typename S2,int O2>
    InertiaTpl( const InertiaTpl<S2,O2> & clone )
    : m(clone.mass()),
    c(clone.lever()),
    I(clone.inertia().matrix())
    {

    }

    // Initializers
    static InertiaTpl Zero() 
    {
      return InertiaTpl(0., 
                        Vector3::Zero(), 
                        Symmetric3::Zero());
    }
    
    void setZero() { m = 0.; c.setZero(); I.setZero(); }

    static InertiaTpl Identity() 
    {
      return InertiaTpl(1., 
                        Vector3::Zero(), 
                        Symmetric3::Identity());
    }
    
    void setIdentity () { m = 1.; c.setZero(); I.setIdentity(); }

    static InertiaTpl Random()
    {
        // We have to shoot "I" definite positive and not only symmetric.
      return InertiaTpl(Eigen::internal::random<Scalar_t>()+1,
                        Vector3::Random(),
                        Symmetric3::RandomPositive());
    }

    static InertiaTpl FromEllipsoid(
        const Scalar_t m, const Scalar_t x, const Scalar_t y, const Scalar_t z)
    {
      Scalar_t a = m * (y*y + z*z) / 5;
      Scalar_t b = m * (x*x + z*z) / 5;
      Scalar_t c = m * (y*y + x*x) / 5;
      return InertiaTpl(m, Vector3::Zero(), Symmetric3(a, 0, b, 0, 0, c));
    }

    static InertiaTpl FromCylinder(
        const Scalar_t m, const Scalar_t r, const Scalar_t l)
    {
      Scalar_t a = m * (r*r / 4 + l*l / 12);
      Scalar_t c = m * (r*r / 2);
      return InertiaTpl(m, Vector3::Zero(), Symmetric3(a, 0, a, 0, 0, c));
    }

    static InertiaTpl FromBox(
        const Scalar_t m, const Scalar_t x, const Scalar_t y, const Scalar_t z)
    {
      Scalar_t a = m * (y*y + z*z) / 12;
      Scalar_t b = m * (x*x + z*z) / 12;
      Scalar_t c = m * (y*y + x*x) / 12;
      return InertiaTpl(m, Vector3::Zero(), Symmetric3(a, 0, b, 0, 0, c));
    }

    
    void setRandom()
    {
      m = static_cast<Scalar_t> (std::rand()) / RAND_MAX;
      c.setRandom(); I.setRandom();
    }

    inline Matrix6 matrix() const
    {
      Matrix6 M;
      M.template block<3,3>(LINEAR, LINEAR ).setZero();
      M.template block<3,3>(LINEAR, LINEAR ).diagonal().fill(m);
      M.template block<3,3>(ANGULAR,LINEAR ) = alphaSkew(m,c);
      M.template block<3,3>(LINEAR, ANGULAR) = -M.template block<3,3> (ANGULAR, LINEAR);
      M.template block<3,3>(ANGULAR,ANGULAR) = I.matrix();
      for(int k=0; k<3; ++k)
        M.template block<3,3>(ANGULAR,ANGULAR).col(k) -= c.cross(M.template block<3,3>(ANGULAR, LINEAR).col(k));

      return M;
    }
    
    operator Matrix6 () const { return matrix(); }

    // Arithmetic operators
    template<typename OtherScalar, int OtherOptions>
    inline void setTo (InertiaTpl<OtherScalar,OtherOptions> & dest) const
    {
      dest.mass() = m;
      dest.lever() = c;
      dest.inertia() = I;
    }

    // Requiered by std::vector boost::python bindings. 
    inline bool isEqual( const InertiaTpl& Y2 ) const
    { 
      return (m==Y2.m) && (c==Y2.c) && (I==Y2.I);
    }

    template<typename OtherScalar, int OtherOptions>
    inline  InertiaTpl add(const InertiaTpl<OtherScalar,OtherOptions> & Yb) const
    {
      /* Y_{a+b} = ( m_a+m_b,
       *             (m_a*c_a + m_b*c_b ) / (m_a + m_b),
       *             I_a + I_b - (m_a*m_b)/(m_a+m_b) * AB_x * AB_x )
       */

      const double & mab = m+Yb.mass();
      const Vector3 & AB = (c-Yb.lever());
      return InertiaTpl(mab,
                        (m*c+Yb.mass()*Yb.lever())/mab,
                        I+Yb.inertia() - (m*Yb.mass()/mab) * typename Symmetric3::SkewSquare(AB));
    }

    
    template<typename OtherScalar, int OtherOptions>
    inline InertiaTpl & operator+= (const InertiaTpl<OtherScalar,OtherOptions> & Yb)
    {
      const InertiaTpl& Ya = *this;
      const double & mab = Ya.m+Yb.m;
      const Vector3 & AB = (Ya.c-Yb.c);
      c *= m; c += Yb.mass()*Yb.lever(); c /= mab;
      I += Yb.inertia(); I -= (Ya.m*Yb.mass()/mab) * typename Symmetric3::SkewSquare(AB);
      m  = mab;
      return *this;
    }
    
    template<typename OtherScalar, int OtherOptions>
    inline void addTo (InertiaTpl<OtherScalar,OtherOptions> & dest) const
    {
      dest += *this;
    }

//    Force __mult__(const Motion &v) const 
//    {
//      Force f;
//      f.linear() = m*(v.linear() - c.cross(v.angular()));
//      f.angular() = c.cross(f.linear()) + I*v.angular();
//      return f;
//    }
    
    inline ForceTpl<_Scalar,_Options> operator* (const MotionTpl<_Scalar,_Options> & v) const
    {
      typedef ForceTpl<_Scalar,_Options> ForceType;
      ForceType f;
      f.linear() = m*(v.linear() - c.cross(v.angular()));
      f.angular() = c.cross(f.linear()) + I*v.angular();
      return f;
    }
    
    template<typename MotionScalar, int MotionOptions>
    inline Scalar_t vtiv(const MotionTpl<MotionScalar,MotionOptions> & v) const
    {
      const Vector3 cxw (c.cross(v.angular()));
      Scalar_t res = m * (v.linear().squaredNorm() - 2.*v.linear().dot(cxw));
//      const Vector3 mcxcxw (-m*c.cross(cxw));
      res -= m * v.angular().dot(c.cross(cxw));
      res += I.vtiv(v.angular());
      
      return res;
    }

    // Getters
    Scalar_t           mass()    const { return m; }
    const Vector3 &    lever()   const { return c; }
    const Symmetric3 & inertia() const { return I; }
    
    Scalar_t &   mass()    { return m; }
    Vector3 &    lever()   { return c; }
    Symmetric3 & inertia() { return I; }

    /// aI = aXb.act(bI)
    template<typename SE3Scalar, int SE3Options>
    inline  InertiaTpl SE3ActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    {
      /* The multiplication RIR' has a particular form that could be used, however it
       * does not seems to be more efficient, see http://stackoverflow.com/questions/
       * 13215467/eigen-best-way-to-evaluate-asa-transpose-and-store-the-result-in-a-symmetric .*/
       return InertiaTpl(m,
                         M.translation()+M.rotation()*c,
                         I.rotate(M.rotation()));
     }

    ///bI = aXb.actInv(aI)
    template<typename SE3Scalar, int SE3Options>
    inline InertiaTpl SE3InvActOn(const SE3Tpl<SE3Scalar,SE3Options> & M) const
    {
      return InertiaTpl(m,
                        M.rotation().transpose()*(c-M.translation()),
                        I.rotate(M.rotation().transpose()) );
    }

    template<typename MotionScalar, int MotionOptions>
    inline ForceTpl<MotionScalar,MotionOptions> vxiv(const MotionTpl<MotionScalar, MotionOptions> & v) const
    {
      typedef ForceTpl<MotionScalar,MotionOptions> ReturnType;
      const Vector3 & mcxw = m*c.cross(v.angular());
      const Vector3 & mv_mcxw = m*v.linear()-mcxw;
      return ReturnType(v.angular().cross(mv_mcxw),
                        v.angular().cross(c.cross(mv_mcxw)+I*v.angular())-v.linear().cross(mcxw) );
    }

    void disp(std::ostream & os) const
    {
      os
      << "  m = " << m << "\n"
      << "  c = " << c.transpose() << "\n"
      << "  I = \n" << (Matrix3)I << "\n";
    }

  protected:
    Scalar_t m;
    Vector3 c;
    Symmetric3 I;
    
  }; // class InertiaTpl
  

  typedef InertiaTpl<double,0> Inertia;
    
} // namespace se3

#endif // ifndef __se3_inertia_hpp__
