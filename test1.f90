
double precision function runge_kuttaf(f1, g1, i)
use constants
  integer                                          ::i
  double precision                                 ::f1
  complex*16                                       ::g1,p
!  double precision, parameter                      ::gamma = 0.2d12
  double precision, parameter                      ::omega_r = 3.2*1.6d-19/6.63d-34
  double precision                                 ::t, k1, k2, k3, k4
  double precision                                 ::omegat1
  double precision                                 ::funcf
  funcf(p) = aimag(omegat1*p)*2
 runge_kuttaf = funcf(g1)
  k1 = dt * funcf(g1)
  k2 = dt * funcf(g1+k1/2)
  k3 = dt * funcf(g1+k2/2)
  k4 = dt * funcf(g1+k3)
!common omegat(1:10)
write(*,*) gamma
end function runge_kuttaf



program main
  use constants 
  implicit none
  complex*16                                     ::p0 = 0d0
!  double precision                               ::gamma = 0.2d12
  double precision, parameter                    ::f0 = 0d0, t_end = 1d-12, &
                                                   T = 0.1d-12 , tgrid = t_end/dt
  double precision                               ::f(tgrid+1,1)
  complex*16                                     ::p(tgrid+1,1)
  double precision                               ::runge_kuttaf

  gamma = 1
  p0 = runge_kuttaf(1d0, cmplx(1d0,1d0), 1d0, 10)

!write(*,*) omegat
end program main
