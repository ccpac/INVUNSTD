module global
  integer tri,trj,tck
  real(8) dt,a,b,c,drx,dry,dcx,dcy
  parameter(dt=0.02d0,tri=24,trj=24,a=0.12d0,b=0.12d0,c=0.003d0)
  parameter(drx=a/dble(tri),dry=b/dble(trj))
  parameter(dcx=drx,dcy=dry)
end module

program SSKF
  use global
  implicit none
  integer i,j,t,pos,nt,posi,posj
  real(8) T0,sigmay,cmh,cmh2,tbeg,tend,ktr,aux0,aux1,rhocr,tt,kt,rhoc,Tref,x,y
  parameter(nt=100)
  real(8),dimension(2*tri*trj)::vX
  real(8),dimension(tri*trj)::vY
  real(8),dimension(2*tri*trj,2*tri*trj)::mF1,mP,maux,maux2
  real(8),dimension(tri*trj,2*tri*trj)::mH
  real(8),dimension(2*tri*trj,tri*trj)::mK
  !!!!!!!!!!
  !Data I/O!
  !!!!!!!!!!
    open(unit=10,file="test_0.dat",status="replace")
    open(unit=11,file="qest.dat",status="replace")
    open(unit=12,file="test_c.dat",status="replace")
    open(unit=16,file="residual.dat",status="replace")
    open(unit=20,file="time_t.dat",status="replace")
    open(unit=21,file="time_q.dat",status="replace")
    posi=2
    posj=12
  !!!!!!!!!!!!!!!!!!!!!!
  ! Thermal Properties !
  !!!!!!!!!!!!!!!!!!!!!!
    Tref=600.d0                         ! Reference temperature in K
    ktr=kt(Tref)                        ! Constant thermal conductivity, in W/mK
    rhocr=rhoc(Tref)                    ! Constant heat capacity, in J/m3K
  !!!!!!!!!!!!!!!!!!!!
  !Initialization and!
  ! basic parameters !
  !!!!!!!!!!!!!!!!!!!!
    T0=300.d0                           ! Initial condition in K
    sigmay=1.25d0                       ! Std. deviation of temperature measurements
    vX(1:tri*trj)=T0                    ! Avg. temperatura state variables
    vX(tri*trj+1:2*tri*trj)=0.d0        ! Heat flux state variables
  !!!!!!!!!!!!!!!!!!!!!!
  ! Observation Matrix !
  !!!!!!!!!!!!!!!!!!!!!!
    cmh=-c/(6.d0*ktr)
    cmh2=c/(3.d0*ktr)
    mH=0.d0
    do i=1,tri*trj
      mH(i,i)=1.d0
      mH(i,tri*trj+i)=cmh
    enddo
  !!!!!!!!!!!!!!!
  !Steady-State !
  !Kalman Filter!
  !!!!!!!!!!!!!!!
    open(unit=15,file='steady_mPmK.dat',status='unknown',form='unformatted',access='stream')
    read(unit=15)mF1,mP,mK
    close(unit=15)
    maux2=matmul(mK,mH)
    maux2=-maux2
    do i=1,2*tri*trj
      maux2(i,i)=maux2(i,i)+1.d0
    enddo
    maux=matmul(maux2,mF1)
    open(unit=15,file='vy_raw.dat',status="old")
    call cpu_time(tbeg)
    do t=1,nt
      tt=real(t,8)*dt
      write(unit=*,fmt=*)tt
      !!!!!!!!!!!!!!!!!!!!!!
      !Reading Measurements!
      !!!!!!!!!!!!!!!!!!!!!!
        read(unit=15,fmt=*)vy
      !!!!!!!!!!!!
      !Prediction!
      ! & Update !
      !!!!!!!!!!!!
        vX=matmul(maux,vX)+matmul(mK,vY)
      !!!!!!!!!!!!!
      !Data Output!
      !!!!!!!!!!!!!
        write(unit=10,fmt=99)'vT','"Temperature [K], "IC99%-", "IC99%+"',tri,trj,t,tt
        write(unit=11,fmt=99)'vq','"Heat Flux [W/m2]", "IC99%-", "IC99%+"',tri,trj,t,tt
        write(unit=12,fmt=99)'vT','"Temperature [K], "IC99%-", "IC99%+"',tri,trj,t,tt
        write(unit=16,fmt=99)'vY-vT','"Residuals [K]"',tri,trj,t,tt
        do i=1,tri
          x=(real(i,8)-0.5d0)*drx
          do j=1,trj
            pos=(i-1)*trj+j
            y=(real(j,8)-0.5d0)*dry
            write(unit=10,fmt='(5(es14.6,1x))')x,y,vX(pos)+cmh*vX(tri*trj+pos),&
              vX(pos)+cmh*vX(tri*trj+pos)-&
                2.576d0*dsqrt(mP(pos,pos)+cmh*(mP(pos,tri*trj+pos)+&
                mP(tri*trj+pos,pos))+(cmh**2.d0)*mP(tri*trj+pos,tri*trj+pos)),&
              vX(pos)+cmh*vX(tri*trj+pos)+&
                2.576d0*dsqrt(mP(pos,pos)+cmh*(mP(pos,tri*trj+pos)+&
                mP(tri*trj+pos,pos))+(cmh**2.d0)*mP(tri*trj+pos,tri*trj+pos))
            write(unit=12,fmt='(5(es14.6,1x))')x,y,vX(pos)+cmh2*vX(tri*trj+pos),&
              vX(pos)+cmh2*vX(tri*trj+pos)-&
                2.576d0*dsqrt(mP(pos,pos)+cmh2*(mP(pos,tri*trj+pos)+&
                mP(tri*trj+pos,pos))+(cmh2**2.d0)*mP(tri*trj+pos,tri*trj+pos)),&
              vX(pos)+cmh2*vX(tri*trj+pos)+&
                2.576d0*dsqrt(mP(pos,pos)+cmh2*(mP(pos,tri*trj+pos)+&
                mP(tri*trj+pos,pos))+(cmh2**2.d0)*mP(tri*trj+pos,tri*trj+pos))
            pos=tri*trj+pos
            aux0=2.576d0*dsqrt(mP(pos,pos))
            write(unit=11,fmt='(5(es14.6,1x))')x,y,vX(pos),vX(pos)-aux0,vX(pos)+aux0
            write(unit=16,fmt='(3(es14.6,1x))')x,y,vY(trj*(i-1)+j)-(vX(trj*(i-1)+j)+cmh*vX(tri*trj+trj*(i-1)+j))
          enddo
        enddo
        pos=(posi-1)*trj+posj
        aux0=vX(pos)+cmh*vX(pos+tri*trj)
        aux1=2.576d0*dsqrt(mP(pos,pos)+cmh*(mP(pos,tri*trj+pos)+&
             mP(tri*trj+pos,pos))+(cmh**2.d0)*mP(tri*trj+pos,tri*trj+pos))
        write(unit=20,fmt='(5(es14.6,1x))')tt,aux0,aux0-aux1,aux0+aux1,aux0-vy(pos)
        pos=pos+tri*trj
        aux0=vX(pos)
        aux1=2.576d0*dsqrt(mP(pos,pos))
        write(unit=21,fmt='(4(es14.6,1x))')tt,aux0,aux0-aux1,aux0+aux1
    enddo
  close(unit=10)
  close(unit=11)
  close(unit=12)
  close(unit=15)
  close(unit=16)
  close(unit=20)
  close(unit=21)
  call cpu_time(tend)
  open(unit=12,file="time.dat",status='replace')
  write(unit=12,fmt=*)tend-tbeg
  close(unit=12)
  stop
  !!!!!!!!!!!!!!!!!!!!
  ! Formatting Rules !
  ! and cmd template !
  !!!!!!!!!!!!!!!!!!!!
    !write(unit=10,fmt=99)'vT','"Temperature [K]"',nx,ny,t,tt
    99 format ('TITLE = ',A,/,&
      'Variables = "x [m]", "y [m]",',A,/,&
      'ZONE, i=',i3,', j=',i3,', f=point, STRANDID=',i3,',SOLUTIONTIME=',es14.6)
end program SSKF


!!!!!!!!!!!!!!!!
!Thermophysical!
!  Properties  !
!!!!!!!!!!!!!!!!
real(8) function kt(T)
    real(8) T
    kt=12.45d0+.014d0*T+2.517d-6*(T**2.d0)
end function

real(8) function rhoc(T)
    real(8) T
    rhoc=1324.75d0*T+3557900.d0
end function