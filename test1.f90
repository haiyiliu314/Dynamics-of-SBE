
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
  k1 = dt * funcf(g1)
  k2 = dt * funcf(g1+k1/2)
  k3 = dt * funcf(g1+k2/2)
  k4 = dt * funcf(g1+k3)
  runge_kuttaf = f1 + k1/6 +k2/3 + k3/3 + k4/6
!common omegat(1:10)
!write(*,*) gamma
end function runge_kuttaf

complex*16 function runge_kuttap(g1, i)
use constants
  integer                                          ::i
  complex*16                                       ::g1,p
!  double precision, parameter                      ::gamma = 0.2d12
  double precision, parameter                      ::omega_r = 3.2*1.6d-19/6.63d-34
  double precision                                 ::t
  complex*16                                       ::funcp, k1, k2, k3, k4
  funcp(p) = -(0,1)*(-((0,1)*gamma) * p - omegat)
  k1 = dt * funcp(g1)
  k2 = dt * funcp(g1+k1/2)
  k3 = dt * funcp(g1+k2/2)
  k4 = dt * funcp(g1+k3)
  runge_kuttap = g1 + k1/6 +k2/3 + k3/3 + k4/6
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
  double precision                               ::f(tgrid+1,1)
  complex*16                                     ::p(tgrid+1,1), f6
  double precision                               ::runge_kuttaf

  double precision                               ::ptest(tgrid+1, 1), diff
  integer                                        ::i
  complex*16                                     ::runge_kuttap
  character(80)                                  :: list_file

  p(1,1) = p0
  f(1,1) = f0
  do i = 1,tgrid
  p(i+1,1) = runge_kuttap(p(i,1), i)
  f(i+1,1) = runge_kuttaf(f(i,1), p(i,1), i)
  end do
  write(*,*) 'tgrid=' , tgrid
  ptest = abs(p)**2;
  diff = abs(abs(p(tgrid+1,1)) - sqrt(0.1d0*((1d0-exp(-gamma*t_end)) / (gamma*1d-12))**2)) / sqrt(abs(0.1d0*((1d0-exp(-gamma*t_end)) / (gamma*1d-12))**2))
  write (*,*) 'p(1001) =' , p(1001,1)
  write(*,*) 'theory peak p', sqrt(0.1d0*((1d0-exp(-0.2d12*1d-12)) / (0.2d12*1d-12))**2)
  write(*,*) 'diff=', diff
  write(list_file, '(A)') 'psquared.dat'
  open(unit=700,file=list_file)

  DO i = 1, tgrid+1

    write(700,*)   ptest(i, 1)
  END DO
  close(700)
!  call system('matlab -r test04202017')
!write(*,*) omegat
end program main
