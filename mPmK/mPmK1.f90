module global
  integer tri,trj,tck
  real(8) dt,a,b,c,drx,dry,dcx,dcy
  parameter(dt=0.02d0,tri=24,trj=24,a=0.12d0,b=0.12d0,c=0.003d0)
  parameter(drx=a/dble(tri),dry=b/dble(trj))
  parameter(dcx=drx,dcy=dry)
end module

program mPmK
  use global
  implicit none
  external SB02MD
  integer i,j,N,M
  integer LDALL,LDR
  parameter(N=2*tri*trj,M=tri*trj)
  parameter(LDALL=N,LDR=N)
  integer,parameter::lwork=100*N
  integer,dimension(M)::ipiv
  real(8),dimension(lwork)::work
  real(8) rhocr,ktr,sigmat,sigmay,sigmaq,&
    r11,r12,r21,r22,rx,ry,rt,dtt,cmh,kt,rhoc,Tref
  real(8),dimension(N,N)::mF1
  real(8),dimension(LDALL,N)::mP,mQ,maux,mG
  real(8),dimension(M,N)::mH
  real(8),dimension(LDALL,M)::mK,maux2
  real(8),dimension(LDR,M)::mR
  real(8),dimension(M,M)::mbInv,mInv
  !!!!!!!!!!!!!!!!!!
  !SLICOT Variables!
  !!!!!!!!!!!!!!!!!!
    real(8) rcond
    integer LDS,LDT,LDU,liwork,ldwork,info
    parameter(LDS=2*N+M,LDT=2*N+M,&
      LDU=2*N,liwork=max(M,2*N),&
      ldwork=max(7*(2*N+1)+16,16*N,2*N+M,3*M))
    integer,dimension(liwork)::iwork
    real(8),dimension(ldwork)::dwork
    real(8),dimension(2*N)::alfar,alfai,beta
    real(8),dimension(LDS,2*N)::S
    real(8),dimension(LDT,2*N)::T
    real(8),dimension(LDU,2*N)::U
    logical,dimension(2*N)::bwork
  !!!!!!!!!!!!!!!!
  !Setting values!
  !!!!!!!!!!!!!!!!
    Tref=600.d0
    rhocr=rhoc(Tref)
    ktr=kt(Tref)
    sigmat=1.d-1
    sigmaq=1.d5
    sigmay=1.25d0
  !!!!!!!!!!!!!!!!!!!!!
  !Assembling Matrices!
  !!!!!!!!!!!!!!!!!!!!!
    dtt=dt
    r11=1.d0-ktr*dtt*(drx**(-2.d0)+dry**(-2.d0))/rhocr
    r12=1.d0-ktr*dtt*(drx**(-2.d0)+2.d0*dry**(-2.d0))/rhocr
    r21=1.d0-ktr*dtt*(2.d0*drx**(-2.d0)+dry**(-2.d0))/rhocr
    r22=1.d0-2.d0*ktr*dtt*(drx**(-2.d0)+dry**(-2.d0))/rhocr
    rx=ktr*dtt/rhocr/drx**2.d0
    ry=ktr*dtt/rhocr/dry**2.d0
    rt=dtt/c/rhocr
    mF1=0.d0
    mF1(1,1)=r11
    do i=2,tri-1
      mF1(i,i)=r21
    enddo
    mF1(tri,tri)=r11
    do i=tri+1,tri*(trj-1),tri
      mF1(i,i)=r12
      do j=2,tri-1
        mF1(i+j-1,i+j-1)=r22
      enddo
      mF1(tri+i-1,tri+i-1)=r12
    enddo
    mF1(tri*(trj-1)+1,tri*(trj-1)+1)=r11
    do i=tri*(trj-1)+2,tri*trj-1
      mF1(i,i)=r21
    enddo
    mF1(tri*trj,tri*trj)=r11
    do i=1,tri
      mF1((i-1)*trj+1,(i-1)*trj+2)=rx
      do j=2,trj-1
        mF1((i-1)*trj+j,(i-1)*trj+j-1)=rx
        mF1((i-1)*trj+j,(i-1)*trj+j+1)=rx
      enddo
      mF1(i*trj,i*trj-1)=rx
    enddo
    do i=1,tri
      mF1(i,i+trj)=ry
    enddo
    do i=tri+1,(tri-1)*trj
      mF1(i,i-trj)=ry
      mF1(i,i+trj)=ry
    enddo
    do i=(tri-1)*trj+1,tri*trj
      mF1(i,i-trj)=ry
    enddo
    do i=1,tri*trj
      mF1(i,tri*trj+i)=rt
      mF1(tri*trj+i,tri*trj+i)=1.d0
    enddo
    cmh=-c/(6.d0*ktr)
    mH=0.d0
    do i=1,tri*trj
      mH(i,i)=1.d0
      mH(i,tri*trj+i)=cmh
    enddo
    mQ=0.d0
    mR=0.d0
    do i=1,M
      mQ(i,i)=sigmat**2
      mQ(tri*trj+i,tri*trj+i)=sigmaq**2
      mR(i,i)=sigmay**2
    enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Calculating Steady-State !
  !Covariance and Kalman Gain!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mP=0.d0
    mK=0.d0
    maux=0.d0
    do i=1,N
      do j=1,N
        maux(i,j)=mF1(j,i)
      enddo
    enddo
    do i=1,M
      do j=1,N
        mK(j,i)=mH(i,j)
      enddo
    enddo
    mG=matmul(transpose(mH),mH)/sigmay**2
    maux2=0.d0
    write(unit=*,fmt=*)"Solving Riccati Equation..."
    CALL SB02OD('D','G','N','L','Z','S', N, M, 1,&
      maux,LDALL,mG,LDALL,mQ,LDALL,mR,LDR,maux2,LDALL,&
      RCOND,mP,LDALL,ALFAR,ALFAI,BETA,S,LDS,&
      T,LDT,U,LDU,0,IWORK,DWORK,LDWORK,BWORK,INFO)
    write(unit=*,fmt=*)"Calculating K matrix..."
    maux2=matmul(mP,transpose(mH))
    mbInv=matmul(mH,maux2)
    do i=1,tri*trj
      mbInv(i,i)=mbInv(i,i)+sigmay**2.d0
    enddo
    mInv=mbInv
    call dgetrf(M,M,mInv,M,ipiv,info)
    call dgetri(M,mInv,M,ipiv,work,lwork,info)
    mK=matmul(maux2,mInv)
  !!!!!!!!!!!!!
  !Data Output!
  !!!!!!!!!!!!!
    write(unit=*,fmt=*)"Writing output files..."
    open(unit=10,file='steady_mPmK.dat',status='unknown',form='unformatted')
    write(unit=10)mF1,mP,mK
    close(unit=10)
end program mPmK

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