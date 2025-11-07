!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This code compute the background cosmological quantities as well
!!!! as first order and second order growth factor of matter density
!!!! fluctuations stored in tables for input of RAMSES and MPGRAFIC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! -----------------------------------------------------------------
!!!! Developper: Pier-Stefano Corasaniti
!!!! Version 2: f90 (overrun previous f77 version)
!!!! Version 3: includes D_plus_second_order calculation
!!!! Version 4: increases accuracy of the calculation of the different
!!!!            time definitions as function of scale factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE COSMO

  ! FUNDAMENTAL CONSTANTS & CONVERSION UNITS

  REAL(8), PARAMETER :: c = 2.99792458d+8, pi = 3.14159265358979, G = 6.6738d-11
  REAL(8), PARAMETER :: kappa = 8.d0*pi*G, sigma_B = 5.6704d-8
  REAL(8), PARAMETER :: Mpc = 3.085678d+22
  REAL(8), PARAMETER :: H0_INV = 9.77814484161389D+9

  ! RADIATIVE DOF PARAMETER VALUES

  REAL(8), PARAMETER :: Num_Nu_massless=3.046, Tcmb = 2.7255d0

  ! COSMO PARAMS

  REAL(8) :: h, Omega_m, w0, wa, Omega_r

  ! SAMPLED REDSHIFTS

  INTEGER, PARAMETER :: NZ = 100000

  ! TIME INTEGRATION VARIABLES

  INTEGER, PARAMETER :: NZ_INT = 1000
  REAL(8), PARAMETER :: EPS = 1.0d-7
  INTEGER, PARAMETER :: NMAX = 100, KMAXX = 100
  REAL(8) :: h1, hmin
  INTEGER :: nok, nbad, nvar
  

CONTAINS

  FUNCTION w_DE(a)
    IMPLICIT NONE
    REAL(8) :: w_DE, a
    w_DE = w0 + wa * (1.d0-a)
  END FUNCTION w_DE
  

  FUNCTION f_DE(z)
    IMPLICIT NONE
    REAL(8) :: f_DE, z, a
    a = 1.d0/(1.d0+z)
    f_DE = exp(3.*wa*(a-1.d0))*(1.d0+z)**(3.d0*(1.d0+w0+wa))
  END FUNCTION f_DE


  FUNCTION inv_hubblez_on_oneplusz(z)
    IMPLICIT NONE
    REAL(8) :: inv_hubblez_on_oneplusz, z, h2z
    h2z = Omega_m*(1.d0+z)**3 + &
         & Omega_r*(1.d0+z)**4 + &
         & (1.-Omega_m-Omega_r)*f_DE(z)
    inv_hubblez_on_oneplusz = 1.d0/sqrt(h2z)/(1.d0+z)
  END FUNCTION inv_hubblez_on_oneplusz

  FUNCTION inv_hubblez(z)
    IMPLICIT NONE
    REAL(8) :: inv_hubblez, z, h2z
    h2z = Omega_m*(1.d0+z)**3+ &
         & Omega_r*(1.d0+z)**4 + &
         & (1.-Omega_m-Omega_r)*f_DE(z)
    inv_hubblez = 1.d0/sqrt(h2z)
  END FUNCTION inv_hubblez

  FUNCTION Lookback_time(z)
    IMPLICIT NONE
    REAL(8) :: Lookback_time, z
    REAL(8) :: clever
    Lookback_time = clever(0.d0,z,inv_hubblez_on_oneplusz)
  END FUNCTION Lookback_time

  FUNCTION Conformal_time(z)
    IMPLICIT NONE
    REAL(8) :: Conformal_time, z
    REAL(8) :: clever
    Conformal_time = clever(z,1.d+20,inv_hubblez)
  END FUNCTION Conformal_time

  FUNCTION Cosmic_time(z)
    IMPLICIT NONE
    REAL(8) :: Cosmic_time, z
    REAL(8) :: clever
    Cosmic_time = clever(z,1.d+20,inv_hubblez_on_oneplusz)
  END FUNCTION Cosmic_time

  FUNCTION dsuperconf_time(z)
    IMPLICIT NONE
    REAL(8) :: dsuperconf_time, z, h2z, clever
    h2z = Omega_m*(1.d0+z)**3 + &
         & Omega_r*(1.d0+z)**4 + &
         & (1.d0-Omega_m-Omega_r)*f_DE(z)
    dsuperconf_time=(1.d0+z)/sqrt(h2z)
  END FUNCTION dsuperconf_time

  FUNCTION Superconf_time(z)
    IMPLICIT NONE
    REAL(8) :: Superconf_time, z, clever
    Superconf_time=clever(z,1.d+20,dsuperconf_time)
  END FUNCTION Superconf_time

  SUBROUTINE FDERIVS(n,x,y,dydx)
    IMPLICIT NONE
    INTEGER :: n
    REAL(8) :: z, a, Omega_DE, f_OM, f_OR, f_Q
    REAL(8) :: x, y(4), dydx(4)
    REAL(8) :: D2, f, dplus, dplus2
    
    z = exp(-x)-1.d0
    a = 1.d0/(1.d0+z)
    Omega_DE = 1.d0-Omega_r-Omega_m
    f_OM = Omega_m*exp(-3.d0*x)*inv_hubblez(z)**2
    f_OR = Omega_r*exp(-4.d0*x)*inv_hubblez(z)**2
    f_Q = Omega_DE*f_DE(z)*inv_hubblez(z)**2

    f = y(2)
    D2 = y(3)
    dydx(1) = f
    dydx(2) = -f**2 - &
         & f/2.d0*(1.d0-f_OR-3.d0*w_DE(a)*f_Q)+ &
         & 3.D0/2.d0*f_OM
    dydx(3) = y(4)
    dydx(4) = 3.D0/2.d0*f_OM*D2- &
         & 3.D0/2.d0*f_OM*(exp(y(1)))**2- &
         & 1.D0/2.D0*(1.d0-f_OR-3.d0*w_DE(a)*f_Q)*y(4)

  END SUBROUTINE FDERIVS
    

END MODULE COSMO

PROGRAM NEW_DARK_COSMOS
  USE COSMO
  IMPLICIT NONE
  CHARACTER(LEN=80) :: file_ramses, file_mpgrafic
  CHARACTER(LEN=6) :: SELECT_RAD
  REAL(8) :: grho, grhog, grhor, grhornomass
  INTEGER :: I
  REAL(8) :: almin, almax, a, z
  REAL(8) :: z_ini, a_ini, xini, xfin
  REAL(8) :: alfastart, alfafin, afin, zfin
  REAL(8) :: y(4), dydx(4)
  REAL(8) :: f, dplus, dplus2
  EXTERNAL rkqs

  WRITE(*,*) 'Insert h'
  READ(*,*) h

  WRITE(*,*) 'Insert Omega_m'
  READ(*,*) Omega_m

  WRITE(*,*) 'Insert w'
  READ(*,*) w0

  write(*,*) 'Insert w_a'
  read(*,*) wa

  write(*,*) 'Include radiation (RAD/NO_RAD)'
  read(*,*) SELECT_RAD

  grho = 3*(100.d0*h)**2/c**2*1000**2 
  grhog = kappa/c**2*4*sigma_B/c**3*Tcmb**4*Mpc**2
  grhor = 7.d0/8*(4.d0/11)**(4.d0/3)*grhog
  grhornomass = grhor*Num_Nu_massless

  SELECT CASE (SELECT_RAD)
  CASE("RAD")
     Omega_r = (grhog+grhornomass)/grho
  CASE("NO_RAD")
     Omega_r = 0.d0
  END SELECT


  print*
  WRITE(*,'(A18,T20,E12.6)') 'Omega_photons =  ', grhog/grho
  WRITE(*,'(A18,T20,E12.6)') 'Omega_radiation =', Omega_r
  WRITE(*,'(A18,T20,E12.6)') 'Omega_de =       ', 1.d0-Omega_m-Omega_r
  print*

  WRITE(*,*) 'Insert RAMSES-Input file'
  READ(*,*) file_ramses
  WRITE(*,*) file_ramses

  WRITE(*,*) 'Insert MPGRAFIC-Input file'
  READ(*,*) file_mpgrafic
  WRITE(*,*) file_mpgrafic

  ! COMPUTE RAMSES INPUT FILE - BACKGROUND EVOLUTION

  OPEN(UNIT=2,FILE=TRIM(file_ramses),STATUS='unknown',FORM='formatted')

  almin = -4.d0   ! log10(a_ini)
  almax = 0.d0    ! log10(a_today=1)

  RAMSES: DO I = 1, NZ

     a = 10.d0**(almin+(almax-almin)*DBLE(I-1)/DBLE(NZ-1))
     z = 1.d0/a-1.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT FOR RAMSES INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     WRITE(2,'(5E20.10)') a, 1./inv_hubblez(z)*a**2, &
          & Superconf_time(z)-Superconf_time(0.d0), &
          & Lookback_time(0.d0)-Lookback_time(z), &
          & Conformal_time(z)-Conformal_time(0.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END DO RAMSES

  CLOSE(2) 

  print*,'Universe Age (Gyr):',Cosmic_time(0.d0)*H0_INV*1.d-9/h
  print*

  ! COMPUTE MPGRAFIC INPUT FILE

  OPEN(UNIT=3,FILE=TRIM(file_mpgrafic),STATUS='unknown',FORM='formatted')    

  h1 = 1.d-7
  hmin = 0.d0
  nvar = 4
  nok = 0  
  nbad = 0

  a_ini = 1.d-8              ! initial value of the scale factor
  z_ini = 1.d0/a_ini-1.d0
  xini = log(1.d0/(1.d0+z_ini))
  xfin = 0.d0

  y(1) = log((3.d0/5.d0*a_ini*1d+4))
  y(2) = 0.d0
  y(3) = -3./7.*(exp(y(1)))**2
  y(4) = 0.d0

  MPGRAFIC: DO I = 1, NZ_INT

     alfafin  = xini + DBLE(I)*(xfin-xini)/DBLE(NZ_INT)
     alfastart = xini + DBLE(I-1)*(xfin-xini)/DBLE(NZ_INT)

     call odeint(y,nvar,alfastart,alfafin,EPS,h1,hmin,nok,nbad,FDERIVS,rkqs)

     dplus = exp(y(1))
     f = y(2)
     dplus2 = y(3)

     afin = exp(alfafin)
     zfin = 1.d0/afin-1.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! OUTPUT FOR MPGRAFIC INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!

     WRITE(3,'(5E20.10)')afin,afin**2/inv_hubblez(zfin),dplus,f,dplus2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END DO MPGRAFIC
  CLOSE(3)

END PROGRAM NEW_DARK_COSMOS


FUNCTION CLEVER(A,B,F)
  IMPLICIT NONE
  REAL(8) :: CLEVER
  integer, parameter :: N = 800
  integer :: LORR(N), IFLAG, LVL
  real(8) :: FV(5),FIT(N),F2T(N),F3T(N),DAT(N)
  real(8) ::  ARESTT(N),ESTT(N),EPST(N),PSUM(N), F
  real(8) :: U, ACC, FOURU, EPS, ERROR, ALPHA, DA, AREA, AREST, KOUNT
  real(8) :: WT, EST, DX, ESTL, ESTR, SUM
  real(8) :: ARESTL, ARESTR, DIFF, A, B, S, T,G

  U=1.0E-7
  ACC=1.0E-7
  FOURU=4.0*U
  IFLAG=1
  EPS=ACC
  ERROR=0.0
  LVL=1
  LORR(LVL)=1
  PSUM(LVL)=0.0
  ALPHA=A
  DA=B-A
  AREA=0.0
  AREST=0.0
  FV(1)=F(ALPHA)
  FV(3)=F(ALPHA+0.5*DA)
  FV(5)=F(ALPHA+DA)
  KOUNT=3
  WT=DA/6.
  EST=WT*(FV(1)+4.*FV(3)+FV(5))

1 DX=.5*DA
  FV(2)=F(ALPHA+0.5*DX)
  FV(4)=F(ALPHA+1.5*DX)
  KOUNT=KOUNT+2
  WT=DX/6.
  ESTL=WT*(FV(1)+4.*FV(2)+FV(3))
  ESTR=WT*(FV(3)+4.*FV(4)+FV(5))
  SUM=ESTL+ESTR
  ! Added by Max:
  if (.not.((SUM.ge.0).or.(SUM.le.0))) return ! Sum = NaN, so might as well bail out right now to save time
  ARESTL=WT*(ABS(FV(1))+ABS(4.*FV(2))+ABS(FV(3)))
  ARESTR=WT*(ABS(FV(3))+ABS(4.*FV(4))+ABS(FV(5)))
  AREA=AREA+((ARESTL+ARESTR)-AREST)
  DIFF=EST-SUM

  IF(ABS(DIFF).LE.EPS*ABS(AREA)) GO TO 2
  IF(ABS(DX).LE.FOURU*ABS(ALPHA)) GO TO 51
  IF(LVL.GE.N) GO TO 5
  IF(KOUNT.GE.2000) GO TO 6

  LVL=LVL+1
  LORR(LVL)=0
  FIT(LVL)=FV(3)
  F2T(LVL)=FV(4)
  F3T(LVL)=FV(5)
  DA=DX
  DAT(LVL)=DX
  AREST=ARESTL
  ARESTT(LVL)=ARESTR
  EST=ESTL
  ESTT(LVL)=ESTR
  EPS=EPS/1.4
  EPST(LVL)=EPS
  FV(5)=FV(3)
  FV(3)=FV(2)
  GO TO 1

2 ERROR=ERROR+DIFF/15.
3 IF(LORR(LVL).EQ.0) GO TO 4
  SUM=PSUM(LVL)+SUM
  LVL=LVL-1
  IF(LVL.GT.1) GO TO 3
  CLEVER=SUM
  RETURN

4 PSUM(LVL)=SUM
  LORR(LVL)=1
  ALPHA=ALPHA+DA
  DA=DAT(LVL)
  FV(1)=FIT(LVL)
  FV(3)=F2T(LVL)
  FV(5)=F3T(LVL)
  AREST=ARESTT(LVL)
  EST=ESTT(LVL)
  EPS=EPST(LVL)
  GO TO 1

5 IFLAG=2
  GO TO 2
6 IFLAG=3
  GO TO 2
51 IFLAG=4
  GO TO 2
END FUNCTION CLEVER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine odeint(ystart,nvar,alfa1,alfa2,EPS,hh1,hmin,nok,nbad, &
              & derivs,rkqs)

           implicit double precision (a-h,o-z)

           integer nbad,nok,nvar,KMAXX,MAXSTP,NMAX,kmax,kount
           real*8 EPS,hh1,hmin,ystart(nvar),TINY,alfa1,alfa2
           external derivs,rkqs
           parameter (MAXSTP=1000000,NMAX=100,KMAXX=200,TINY=1.d-30)
           integer i,nstp
           real*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX), &
                & yp(NMAX,KMAXX),yscal(NMAX)
           common /path/ kmax,kount,dxsav,xp,yp


           x=alfa1
           h=sign(hh1,alfa2-alfa1)
           nok=0
           nbad=0
           kount=0

           do i=1,nvar
              y(i)=ystart(i)
           enddo
           if (kmax.gt.0) xsav=x-2.*dxsav
           do nstp=1,MAXSTP
              call derivs(nvar,x,y,dydx)
              do i=1,nvar
                 yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
              enddo
              if(kmax.gt.0)then
                 if(abs(x-xsav).gt.abs(dxsav)) then
                    if(kount.lt.kmax-1)then
                       kount=kount+1
                       xp(kount)=x
                       do i=1,nvar
                          yp(i,kount)=y(i)
                       enddo
                       xsav=x
                    endif
                 endif
              endif
              if((x+h-alfa2)*(x+h-alfa1).gt.0.d0) h=alfa2-x
              call rkqs(y,dydx,nvar,x,h,EPS,yscal,hdid,hnext,derivs)
              if(hdid.eq.h)then
                 nok=nok+1
              else
                 nbad=nbad+1
              endif
              if((x-alfa2)*(alfa2-alfa1).ge.0.)then
                 do i=1,nvar
                    ystart(i)=y(i)
                 enddo

                 if(kmax.ne.0)then
                    kount=kount+1
                    xp(kount)=x
                    do i=1,nvar
                       yp(i,kount)=y(i)
                    enddo
                 endif
                 return
              endif
              h=hnext
              if(abs(hnext).lt.hmin) then
                 h=hmin
              end if
           enddo
           write(*,*) 'too many steps in odeint; x,hmin',x,alfa1,alfa2
           return
         end subroutine odeint

         subroutine rkqs(y,dydx,n,x,htry,EPS,yscal,hdid,hnext,derivs)

           implicit double precision (a-h,o-z)

           integer n,NMAX
           real*8 EPS,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
           external derivs
           parameter (NMAX=100)
           integer i
           real*8 errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW
           real*8   ERRCON,PSHRNK
           parameter (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89d-4)     
           h=htry
1          call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
           errmax=0.
           do i=1,n
              errmax=max(errmax,abs(yerr(i)/yscal(i)))
           enddo
           errmax=errmax/EPS
           if(errmax.gt.1.)then
              h=SAFETY*h*(errmax**PSHRNK)
              if(h.lt.0.1*h)then
                 h=.1*h
              endif
              xnew=x+h
              goto 1
           else
              if(errmax.gt.ERRCON)then
                 hnext=SAFETY*h*(errmax**PGROW)
              else
                 hnext=5.*h
              endif
              hdid=h
              x=x+h

              do i=1,n
                 y(i)=ytemp(i)
              enddo
              return
           endif
         END subroutine rkqs

         subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)

           implicit double precision (a-h,o-z)

           integer n,NMAX
           real*8 h,x,dydx(n),y(n),yerr(n),yout(n)
           external derivs
           parameter (NMAX=100)
           integer i
           real*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX), &
                & ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53, &
                & B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
           parameter (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., &
                & B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, &
                & B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512., &
                & B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378., &
                & C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648., &
                & DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336., &
                & DC6=C6-.25)

           do i=1,n
              ytemp(i)=y(i)+B21*h*dydx(i)
           enddo
           call derivs(n,x+A2*h,ytemp,ak2)
           do i=1,n
              ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
           enddo
           call derivs(n,x+A3*h,ytemp,ak3)
           do i=1,n
              ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
           enddo
           call derivs(n,x+A4*h,ytemp,ak4)
           do i=1,n
              ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53* &
                   & ak3(i)+B54*ak4(i))
           enddo
           call derivs(n,x+A5*h,ytemp,ak5)
           do i=1,n
              ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+ &
                   & B64*ak4(i)+B65*ak5(i))
           enddo
           call derivs(n,x+A6*h,ytemp,ak6)
           do i=1,n
              yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
           enddo
           do i=1,n
              yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6* &
                   & ak6(i))
           enddo
           return
         END subroutine rkck
