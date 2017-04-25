
double precision function runge_kuttaf(f1, g1, i)
use constants
  integer                                          ::i
  double precision                                 ::f1
  complex*16                                       ::g1,p
!  double precision, parameter                      ::gamma = 0.2d12
!  double precision, parameter                      ::omega_r = 3.2*1.6d-19/6.63d-34
  double precision                                 ::t, k1, k2, k3, k4
  double precision                                 ::omegat1
  double precision                                 ::funcf
  funcf(p) = aimag(omegat*p)*2d0
  k1 = dt * funcf(g1)
  k2 = dt * funcf(g1+k1/2d0)
  k3 = dt * funcf(g1+k2/2d0)
  k4 = dt * funcf(g1+k3)
  runge_kuttaf = f1 + k1/6d0 +k2/3d0 + k3/3d0 + k4/6d0
!write(*,*) 'funcf(g1)='
!write(*,*) funcf(g1)
!common omegat(1:10)
!write(*,*) gamma
end function runge_kuttaf

complex*16 function runge_kuttap(g1, i)
use constants
  integer                                          ::i
  complex*16                                       ::g1,p
!  double precision, parameter                      ::gamma = 0.2d12
!  double precision, parameter                      ::omega_r = 3.2*1.6d-19/6.63d-34
  double precision                                 ::t
  complex*16                                       ::funcp, k1, k2, k3, k4
  funcp(p) = -(0d0,1d0)*(-((0d0,1d0)*gamma) * p - omegat)
  k1 = dt * funcp(g1)
  k2 = dt * funcp(g1+k1/2d0)
  k3 = dt * funcp(g1+k2/2d0)
  k4 = dt * funcp(g1+k3)
  runge_kuttap = g1 + k1/6d0 +k2/3d0 + k3/3d0 + k4/6d0
!write(*,*) runge_kuttap
!common omegat(1:10)s
!write(*,*) gamma
end function runge_kuttap

program main
  use constants 
  implicit none
  complex*16                                     ::p0 = (0d0, 0d0)
!  double precision                               ::gamma = 0.2d12
  double precision, parameter                    ::f0 = 0d0, t_end = 1d-12, &
                                                   T = 0.1d-12 , tgrid = nint(t_end/dt)
  double precision                               ::f
  complex*16                                     ::p, f6
  double precision                               ::runge_kuttaf

  double precision                               ::ptest, diff
  integer                                        ::i
  complex*16                                     ::runge_kuttap
  character(80)                                  :: list_file

  p = p0
  f = f0
  do i = 1,tgrid
  p = runge_kuttap(p, i)
  f = runge_kuttaf(f, p, i)
  end do

!  write(*,*) 'tgrid=' , tgrid
  ptest = abs(p)**2;
  diff = abs(abs(p) - sqrt(0.1d0*((1d0-exp(-gamma*t_end)) / (gamma*1d-12))**2)) / sqrt(abs(0.1d0*((1d0-exp(-gamma*t_end)) / (gamma*1d-12))**2))
  write (*,*) 'p(1001) =' , p
  write(*,*) 'theory peak p', sqrt(0.1d0*((1d0-exp(-0.2d12*1d-12)) / (0.2d12*1d-12))**2)
  write(*,*) 'diff=', diff
!  call system('matlab -r test04202017')
!write(*,*) omegat
end program main
