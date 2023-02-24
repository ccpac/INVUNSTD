module global
    integer,parameter::nx=24,ny=24,nt=100,nxy=nx*ny
    real(8)::drx,dry,dt
    real(8)::dcx,dcy,dcz
    real(8)::p1,p2,p3,Tref,kref
end module global

program CKF
    use mkl_service
    use global
    implicit none
    include 'mkl_lapack.fi'
    integer::i,j,k,t,pi,pj,pk
    integer::info,nthreads
    integer,parameter::lwork=100*nxy
    integer,dimension(nxy)::ipiv
    real(8),dimension(lwork)::work
    real(8)::T2Th,Th2T
    real(8)::tt,x,y
    real(8)::r11,r12,r21,r22,rx,ry,rt,t1,t2
    real(8)::ktr,rhocr,a,b,c,kt,rhoc,time,cmh,cmh2,T0,sT,sq,sy,auxt,auxq
    real(8)::tf,ti
    real(8),dimension(nxy)::vye,vy,vyp,vyk,vqe
    real(8),dimension(2*nxy)::vxp,vxm
    real(8),dimension(nxy,nxy)::mR,mInv
    real(8),dimension(nxy,2*nxy)::mH
    real(8),dimension(2*nxy,nxy)::mK,mPHT
    real(8),dimension(2*nxy,2*nxy)::mF1,mQ,mPp,mPm,maux

!!!!!!!!!!!!!!!!
! Problem Data !
!!!!!!!!!!!!!!!!
    
    a=0.12d0
    b=0.12d0
    c=0.003d0
    time=2.d0
    T0=300.d0

!!!!!!!!!!!!!!!!!!!!!
! Measurement Noise !
!!!!!!!!!!!!!!!!!!!!!

    sy=1.d-2

!!!!!!!!!!!!!
! Grid size !
!!!!!!!!!!!!!

    drx=a/real(nx,8)
    dry=b/real(ny,8)
    dt=time/real(nt,8)

    dcx=drx
    dcy=dry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Position for Profiling Results !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pi=9
    pj=9
    pk=pj+(pi-1)*ny

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! General Material Properties !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    p1=12.45d0
    p2=0.014d0
    p3=2.517d-6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constant Thermal Properties !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Tref=600.d0
    ktr=kt(Tref)
    kref=ktr
    rhocr=rhoc(Tref)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Testing Kirchhoff Transform !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open(unit=10,file="kirch.dat",status="replace")
    t1=1.d0
    do i=1,100
        t1=t1+25.d0
        t2=T2Th(t1)
        write(unit=10,fmt='(3(es14.6,1x))')t1,t2,Th2T(t2)
    enddo
    close(unit=10)

!!!!!!!!!!!!!!!!!!!!
! Evolution Matrix !
!!!!!!!!!!!!!!!!!!!!

    r11=1.d0-ktr*dt*(drx**(-2.d0)+dry**(-2.d0))/rhocr
    r12=1.d0-ktr*dt*(drx**(-2.d0)+2.d0*dry**(-2.d0))/rhocr
    r21=1.d0-ktr*dt*(2.d0*drx**(-2.d0)+dry**(-2.d0))/rhocr
    r22=1.d0-2.d0*ktr*dt*(drx**(-2.d0)+dry**(-2.d0))/rhocr
    rx=ktr*dt/rhocr/drx**2.d0
    ry=ktr*dt/rhocr/dry**2.d0
    rt=dt/c/rhocr
    mF1=0.d0
    mF1(1,1)=r11
    do i=2,nx-1
        mF1(i,i)=r21
    enddo
    mF1(nx,nx)=r11
    do i=nx+1,nx*(ny-1),nx
        mF1(i,i)=r12
        do j=2,nx-1
            mF1(i+j-1,i+j-1)=r22
        enddo
        mF1(nx+i-1,nx+i-1)=r12
    enddo
    mF1(nx*(ny-1)+1,nx*(ny-1)+1)=r11
    do i=nx*(ny-1)+2,nx*ny-1
        mF1(i,i)=r21
    enddo
    mF1(nx*ny,nx*ny)=r11
    do i=1,nx
        mF1((i-1)*ny+1,(i-1)*ny+2)=rx
        do j=2,ny-1
            mF1((i-1)*ny+j,(i-1)*ny+j-1)=rx
            mF1((i-1)*ny+j,(i-1)*ny+j+1)=rx
        enddo
        mF1(i*ny,i*ny-1)=rx
    enddo
    do i=1,nx
        mF1(i,i+ny)=ry
    enddo
    do i=nx+1,(nx-1)*ny
        mF1(i,i-ny)=ry
        mF1(i,i+ny)=ry
    enddo
    do i=(nx-1)*ny+1,nx*ny
        mF1(i,i-ny)=ry
    enddo
    do i=1,nx*ny
        mF1(i,nx*ny+i)=rt
        mF1(nx*ny+i,nx*ny+i)=1.d0
    enddo

!!!!!!!!!!!!!!!!!!!!!
! Observation Model !
!!!!!!!!!!!!!!!!!!!!!
    
    mH=0.d0
    cmh=-c/(6.d0*ktr)
    cmh2=c/(3.d0*ktr)
    do i=1,nxy
        mH(i,i)=1.d0
        mH(i,nxy+i)=cmh
    enddo

!!!!!!!!!!!!!!!!!!!!
! Noise Covariance !
!!!!!!!!!!!!!!!!!!!!

    sT=1.d-1
    sq=1.d2
    mQ=0.d0
    mR=0.d0
    do i=1,nxy
        mQ(i,i)=sT**2.d0
        mQ(nxy+i,nxy+i)=sq**2.d0
        mR(i,i)=sy**2.d0
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initializing Remaining !
!  Vectors and Matrices  !
!!!!!!!!!!!!!!!!!!!!!!!!!!

    mK=0.d0
    mPp=mQ
    mPm=mQ
    mPHT=0.d0
    mInv=0.d0
    vxp=0.d0
    vxm=0.d0
    vye=0.d0
    vy=0.d0
    vyk=0.d0
    vyp=0.d0
    vqe=0.d0
    do i=1,nxy
        vxp(i)=T2Th(T0)
        vxp(nxy+i)=0.d0
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! State Estimation Problem!
! Classical Kalman Filter !
!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nthreads=8
    call MKL_SET_NUM_THREADS(nthreads)
    ti=DSECND()
    
    open(unit=10,file="vTest.dat",status="replace")
    open(unit=11,file="vqest.dat",status="replace")
    open(unit=12,file="vy_raw.dat",status="old")
    open(unit=13,file="vTest_c.dat",status="replace")
    open(unit=14,file="vR.dat",status="replace")
    open(unit=15,file="profiles.dat",status="replace")
    do t=1,nt
        tt=real(t,8)*dt
        write(unit=*,fmt=*)tt
    !!!!!!!!!!!!!
    ! Read Data !
    !!!!!!!!!!!!!

        read(unit=12,fmt=*)vye,vy,vqe

    !!!!!!!!!!!!!!!!!!!!!!!
    ! Kirchhoff Transform !
    !!!!!!!!!!!!!!!!!!!!!!!

        do i=1,nxy
            vyk(i)=T2Th(vy(i))
        enddo

    !!!!!!!!!!
    ! x = Fx !
    !!!!!!!!!!

        call dgemv('N',2*nxy,2*nxy,1.d0,mF1,2*nxy,vxp,1,0.d0,vxm,1)
        vxp=vxm

    !!!!!!!!!!!!!!!!!
    ! P = FPF^T + Q !
    !!!!!!!!!!!!!!!!!

        mPm=mPp
        call dgemm('N','N',2*nxy,2*nxy,2*nxy,1.d0,mF1,2*nxy,mPm,2*nxy,0.d0,mPp,2*nxy)
        mPm=mQ
        call dgemm('N','T',2*nxy,2*nxy,2*nxy,1.d0,mPp,2*nxy,mF1,2*nxy,1.d0,mPm,2*nxy)

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! K = PH^T(HPH^T + R)^-1 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!

        mPHT=0.d0
        call dgemm('N','T',2*nxy,nxy,2*nxy,1.d0,mPm,2*nxy,mH,nxy,0.d0,mPHT,2*nxy)
        mInv=mR
        call dgemm('N','N',nxy,2*nxy,nxy,1.d0,mH,nxy,mPHT,2*nxy,1.d0,mInv,nxy)
        call dgetrf(nxy,nxy,mInv,nxy,ipiv,info)
        call dgetri(nxy,mInv,nxy,ipiv,work,lwork,info)
        mK=0.d0
        call dgemm('N','N',2*nxy,2*nxy,nxy,1.d0,mPHT,2*nxy,mInv,nxy,0.d0,mK,2*nxy)

    !!!!!!!!!!!!!!!!!!!
    ! x = x + K(y-Hx) !
    !!!!!!!!!!!!!!!!!!!

        call dgemv('N',nxy,2*nxy,1.d0,mH,nxy,vxm,1,0.d0,vyp,1)
        vxp=vxm
        call dgemv('N',2*nxy,nxy,1.d0,mK,2*nxy,vyk-vyp,1,1.d0,vxp,1)

    !!!!!!!!!!!!!!!
    ! P = (I-KH)P !
    !!!!!!!!!!!!!!!
        
        call dgemm('N','N',2*nxy,2*nxy,nxy,1.d0,mK,2*nxy,mH,nxy,0.d0,maux,2*nxy)
        do i=1,2*nxy
            maux(i,i)=maux(i,i)-1.d0
        enddo
        call dgemm('N','N',2*nxy,2*nxy,2*nxy,1.d0,-maux,2*nxy,mPm,2*nxy,0.d0,mPp,2*nxy)

    !!!!!!!!!!
    ! Output !
    !!!!!!!!!!
        write(*,*)vqe(pk),pk
        write(unit=10,fmt=99)'vT','"T [K]", "y [K]"',nx,ny,t,tt
        write(unit=11,fmt=99)'vq','"q [W/m2]"',nx,ny,t,tt
        write(unit=13,fmt=99)'vT','"T [K]"',nx,ny,t,tt
        write(unit=14,fmt=99)'vR','"Y-T [K]"',nx,ny,t,tt
        do i=1,nx
            x=(real(i,8)-0.5d0)*drx
            do j=1,ny
                y=(real(j,8)-0.5d0)*dry
                k=j+(i-1)*ny
                write(unit=10,fmt=*)x,y,Th2T(vxp(k)+cmh*vxp(nxy+k)),vy(k)
                write(unit=11,fmt=*)x,y,vxp(nxy+k)
                write(unit=13,fmt=*)x,y,Th2T(vxp(k)+cmh*vxp(nxy+k))
                write(unit=14,fmt=*)x,y,vy(k)-Th2T(vxp(k)+cmh*vxp(nxy+k))
            enddo
        enddo
        auxt=2.576d0*dsqrt(mPp(pk,pk))
        auxq=2.576d0*dsqrt(mPp(nxy+pk,nxy+pk))
        write(unit=15,fmt='(11(es15.6, 1x))')tt,vye(pk),vy(pk),Th2T(vxp(pk)+cmh*vxp(nxy+pk)),&
            Th2T(vxp(pk)+cmh*vxp(nxy+pk)-auxt),&
            Th2T(vxp(pk)+cmh*vxp(nxy+pk)+auxt),&
            vy(pk)-Th2T(vxp(pk)+cmh*vxp(nxy+pk)),&
            vqe(pk),&
            vxp(nxy+pk),&
            vxp(nxy+pk)-auxq,&
            vxp(nxy+pk)+auxq
    enddo
    close(unit=10)
    close(unit=11)
    close(unit=12)
    close(unit=13)
    close(unit=14)
    close(unit=15)
    tf=DSECND()
    open(unit=20,file="time.dat",status="replace")
    write(unit=20,fmt=*)nthreads,tf-ti
    close(unit=20)

!!!!!!!!!!!!!!!!!!
!Formatting Rules!
!!!!!!!!!!!!!!!!!!

    99 format ('TITLE = ',A,/,&
        'Variables = "x [m]", "y [m]",',A,/,&
        'ZONE, i=',i3,', j=',i3,', f=point, STRANDID=',i3,',SOLUTIONTIME=',es14.6)
end program CKF

subroutine vqf(vq,it)
    use global
    implicit none
    integer::i,j,k,it,i0,i1,j0,j1
    real(8)::tt,x,y
    real(8),dimension(nx*ny)::vq
    tt=real(it,8)*dt
    i0=8
    i1=10
    j0=8
    j1=10
    do i=1,nx
        x=(real(i,8)-0.5d0)*drx
        do j=1,ny
            k=i+(j-1)*nx
            y=(real(j,8)-0.5d0)*dry
            if((i.ge.i0).and.(i.le.i1).and.(j.ge.j0).and.(j.le.j1))then
                vq(k)=1.d5
            else
                vq(k)=0.d0
            endif
        enddo
    enddo
end subroutine

!!!!!!!!!!!!!!!!
!Thermophysical!
!  Properties  !
!!!!!!!!!!!!!!!!
real(8) function kt(T)
    use global
    implicit none
    real(8)::T
    kt=p1+p2*T+p3*(T**2.d0)
end function

real(8) function rhoc(T)
    implicit none
    real(8)::T
    rhoc=1324.75d0*T+3557900.d0
end function

real(8) function T2Th(T)
    use global
    implicit none
    real(8)::T
    T2Th=((p1*T + p2/2.d0*T**2.d0 + p3/3.d0*T**3.d0)&
        -(p1*Tref + p2/2.d0*Tref**2.d0 + p3/3.d0*Tref**3.d0))/kref
end function

real(8) function Th2T(th)
    use global
    implicit none
    real(8)::th,Told,Tnew,T,tol,eps,f,fold,T2Th
    eps=1.d0
    tol=1.d-6
    T=100.d0
    f=T2Th(T)-th
    Told=0.95*T
    fold=T2Th(Told)-th
    do while(eps.gt.tol)
        Tnew=T-f*(T-Told)/(f-fold)
        eps=dabs(Tnew-T)
        Told=T
        fold=f
        T=Tnew
        f=T2Th(T)-th
    enddo
    Th2T=T
end function