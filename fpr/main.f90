module global
    integer,parameter::nx=24,ny=24,nz=6,nt=250,nxy=nx*ny
    real(8)::drx,dry,dt
    real(8)::dcx,dcy,dcz
    real(8)::a,b,c
end module global

program fpr
    use global
    implicit none
    integer::i,j,k,t,posi,posj,pos
    integer,parameter::lwork=100*nxy
    integer,dimension(4)::iseed
    real(8)::tt,x,y,r
    real(8)::time,T0,sy
    real(8),dimension(nxy)::vy,vye,vq
    real(8),allocatable,dimension(:,:,:)::hT

!!!!!!!!!!!!!!!!
! Problem Data !
!!!!!!!!!!!!!!!!

    a=0.12d0
    b=0.12d0
    c=0.003d0
    time=5.d0
    T0=300.d0

!!!!!!!!!!!!!!!!!!!!!
! Measurement Noise !
!!!!!!!!!!!!!!!!!!!!!

    sy=1.d0

!!!!!!!!!!!!!!!!!
! Random Number !
!   Generator   !
!!!!!!!!!!!!!!!!!

    iseed(1)=0
    iseed(2)=1500
    iseed(3)=3000
    iseed(4)=4095

!!!!!!!!!!!!!
! Grid size !
!!!!!!!!!!!!!

    drx=a/real(nx,8)
    dry=b/real(ny,8)
    dt=time/real(nt,8)

    dcx=drx
    dcy=dry
    dcz=c/real(nz,8)

!!!!!!!!!!!!!!!!!!
! Time Evolution !
!!!!!!!!!!!!!!!!!!

    posi=2
    posj=12
    pos=(posi-1)*ny+posj

!!!!!!!!!!!!!!!!!!!!!!!!!!
! Synthetic Measurements !
! & Reference Heat Flux  !
!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(hT(nx,ny,nz))
    hT=T0
    open(unit=10,file="vq.dat",status="replace")
    open(unit=11,file="vy_raw.dat",status="replace")
    open(unit=12,file="vy.dat",status="replace")
    open(unit=13,file="vTc.dat",status="replace")
    open(unit=14,file="time_exa.dat",status="replace")
    do t=1,nt
        tt=real(t,8)*dt
        vq=0.d0
        call vqf(vq,t)
        call modelocompleto(hT,vq)
        write(unit=10,fmt=99)'vq','"Heat Flux [W/m2]"',nx,ny,t,tt
        write(unit=12,fmt=99)'vy','"Temperature [K]"',nx,ny,t,tt
        write(unit=13,fmt=99)'vTc','"Temperature [K]"',nx,ny,t,tt
        do i=1,nx
            x=(real(i,8)-0.5d0)*drx
            do j=1,ny
                k=ny*(i-1)+j
                y=(real(j,8)-0.5d0)*dry
                call dlarnv(3,iseed,1,r)
                vye(k)=hT(i,j,1)
                vy(k)=vye(k)+r*sy
                write(unit=10,fmt=*)x,y,vq(k)
                write(unit=12,fmt=*)x,y,vy(k)
                write(unit=13,fmt=*)x,y,hT(i,j,nz)
            enddo
        enddo
        write(unit=11,fmt=*)vye,vy,vq
        write(unit=14,fmt='(4(es14.6,1x))')tt,hT(posi,posj,1),vy(pos),vq(pos)
    enddo
    close(unit=10)
    close(unit=11)
    close(unit=12)
    close(unit=13)
    close(unit=14)
    deallocate(hT)

!!!!!!!!!!!!!!!!!!
!Formatting Rules!
!!!!!!!!!!!!!!!!!!

    99 format ('TITLE = ',A,/,&
        'Variables = "x [m]", "y [m]",',A,/,&
        'ZONE, i=',i3,', j=',i3,', f=point, STRANDID=',i3,',SOLUTIONTIME=',es14.6)
end program fpr

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
            k=ny*(i-1)+j
            y=(real(j,8)-0.5d0)*dry
            if((i.ge.i0).and.(i.le.i1).and.(j.ge.j0).and.(j.le.j1))then
                vq(k)=5.d6
            else
                vq(k)=0.d0
            endif
        enddo
    enddo
end subroutine

!subroutine vqf(vq,it)
!    use global
!    implicit none
!    integer::i,j,k,it,i0,j0
!    real(8)::tt,x,y,x0,y0,q0,sig,arg
!    real(8),dimension(nx*ny)::vq
!    tt=real(it,8)*dt
!    i0=1
!    x0=(real(i0,8)-0.5d0)*drx
!    j0=1
!    y0=(real(j0,8)-0.5d0)*dry
!    q0=7.5d6
!    sig=5.d-2
!    do i=1,nx
!        x=(real(i,8)-0.5d0)*drx
!        do j=1,ny
!            k=ny*(i-1)+j
!            y=(real(j,8)-0.5d0)*dry
!            arg=(x-x0)**2.d0+(y-y0)**2.d0
!            vq(k)=q0*dexp(-0.5d0*arg/sig**2.d0)
!        enddo
!    enddo
!end subroutine
 
!subroutine vqf(vq,it)
!    use global
!    implicit none
!    integer::i,j,k,it,i0,j0
!    real(8)::tt,x,y,x0,y0,q0,q1,q2,m
!    real(8),dimension(nx*ny)::vq
!    tt=real(it,8)*dt
!    i0=1
!    x0=(real(i0,8)-0.5d0)*drx
!    j0=1
!    y0=(real(j0,8)-0.5d0)*dry
!    q0=1.d6
!    q1=1.d5
!    q2=4.d5
!    do i=1,nx
!        x=(real(i,8)-0.5d0)*drx
!        do j=1,ny
!            k=ny*(i-1)+j
!            y=(real(j,8)-0.5d0)*dry
!            if((x.le.0.06))then
!                m=(q1-q0)/(0.06d0)
!                vq(k)=q0+m*x
!            endif
!            if((x.ge.0.06d0).and.(x.le.0.09))then
!                vq(k)=q1
!            endif
!            if((x.ge.0.09))then
!                m=(q2-q1)/(0.12d0-0.09d0)
!                vq(k)=q1+m*(x-0.09d0)
!            endif
!        enddo
!    enddo
!end subroutine

!subroutine vqf(vq,it)
!    use global
!    implicit none
!    integer::i,j,k,it,i0,j0
!    real(8)::tt,x,y,x0,y0,q0,q1,q2,m
!    real(8),dimension(nx*ny)::vq
!    tt=real(it,8)*dt
!    i0=1
!    x0=(real(i0,8)-0.5d0)*drx
!    j0=1
!    y0=(real(j0,8)-0.5d0)*dry
!    q0=1.d6
!    q1=1.d5
!    q2=4.d5
!    do i=1,nx
!        x=(real(i,8)-0.5d0)*drx
!        do j=1,ny
!            k=ny*(i-1)+j
!            y=(real(j,8)-0.5d0)*dry
!            if((x.le.0.06))then
!                m=(q1-q0)/(0.06d0)
!                vq(k)=q0+m*x
!            endif
!            if((x.ge.0.06d0).and.(x.le.0.09))then
!                vq(k)=q1
!            endif
!            if((x.ge.0.09))then
!                m=(q2-q1)/(0.12d0-0.09d0)
!                vq(k)=q1+m*(x-0.09d0)
!            endif
!            vq(k)=vq(k)*dexp(-4.d0*(y-b/2.d0)**2.d0/b**2.d0)
!        enddo
!    enddo
!end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Solução do Modelo Completo!
!Método dos Volumes Finitos!
!	     Implicito         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine modelocompleto(hTold,qest)
    use global
    implicit none
    integer::i,j,k
    integer,parameter::tci=nx,tcj=ny,tck=nz
    real(8)::tol,eps,kt,rhoc
    real(8),dimension(tci*tcj)::qest
    real(8),dimension(tci,tcj,tck)::hT,hTold,hTnew,hkt,hrhoc,ap0,ap
    real(8),dimension(tci-1,tcj,tck)::ae
    real(8),dimension(tci,tcj-1,tck)::an
    real(8),dimension(tci,tcj,tck-1)::at


!!!!!!!!!!!!
!Tolerância!
!!!!!!!!!!!!
    tol=1.d-8
!!!!!!!!!!!!!!
! Método de  !
!Gauss-Seidel!
!!!!!!!!!!!!!!
    hTnew=hTold
    eps=1.d0
    do while(eps.gt.tol)
        hT=hTnew
    !!!!!!!!!!!!!!
    !Propriedades!
    !Termofísicas!
    !!!!!!!!!!!!!!
        do i=1,tci
            do j=1,tcj
                do k=1,tck
                    hkt(i,j,k)=kt(hTnew(i,j,k))
                    hrhoc(i,j,k)=rhoc(hTnew(i,j,k))
                enddo
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Coeficientes!
    !!!!!!!!!!!!!!
        do i=1,tci-1
            do j=1,tcj
                do k=1,tck
                    ae(i,j,k)=2.d0*hkt(i,j,k)*hkt(i+1,j,k)/(hkt(i,j,k)+hkt(i+1,j,k))*dcy*dcz/dcx
                enddo
            enddo
        enddo
        do i=1,tci
            do j=1,tcj-1
                do k=1,tck
                    an(i,j,k)=2.d0*hkt(i,j,k)*hkt(i,j+1,k)/(hkt(i,j,k)+hkt(i,j+1,k))*dcx*dcz/dcy
                enddo
            enddo
        enddo
        do i=1,tci
            do j=1,tcj
                do k=1,tck-1
                    at(i,j,k)=2.d0*hkt(i,j,k)*hkt(i,j,k+1)/(hkt(i,j,k)+hkt(i,j,k+1))*dcx*dcy/dcz
                enddo
            enddo
        enddo
        do i=1,tci
            do j=1,tcj
                do k=1,tck
                  ap0(i,j,k)=hrhoc(i,j,k)*dcx*dcy*dcz/dt
                enddo
            enddo
        enddo
    !!!!!!!!!!!!!!!!!!!
    !Pontos Interiores!
    !!!!!!!!!!!!!!!!!!!
        do i=2,tci-1
            do j=2,tcj-1
                do k=2,tck-1
                    ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                    hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
                enddo
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !!!!!!!!!!!!!!
        i=1
        do j=2,tcj-1
            do k=2,tck-1
                ap(i,j,k)=ae(i,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                            an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                            at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !!!!!!!!!!!!!!
        i=tci
        do j=2,tcj-1
            do k=2,tck-1
                ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                            an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                            at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=0!
    !!!!!!!!!!!!!!
        j=1
        do i=2,tci-1
            do k=2,tck-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j,k)*hTnew(i,j+1,k)+&
                    at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                    ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=b!
    !!!!!!!!!!!!!!
        j=tcj
        do i=2,tci-1
            do k=2,tck-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j-1,k)*hTnew(i,j-1,k)+&
                    at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                    ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno z=0!
    !!!!!!!!!!!!!!
        k=1
        do i=2,tci-1
            do j=2,tcj-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                    at(i,j,k)*hTnew(i,j,k+1)+&
                    ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno z=c!
    !!!!!!!!!!!!!!
        k=tck
        do i=2,tci-1
            do j=2,tcj-1
                ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
                hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                    an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                    at(i,j,k-1)*hTnew(i,j,k-1)+&
                    ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
            enddo
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !	      y=0!
    !!!!!!!!!!!!!!
        i=1
        j=1
        do k=2,tck-1
            ap(i,j,k)=ae(i,j,k)+an(i,j,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=b!
    !!!!!!!!!!!!!!	
        i=1
        j=tcj
        do k=2,tck-1
            ap(i,j,k)=ae(i,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=1
        k=1
        do j=2,tcj-1
            ap(i,j,k)=ae(i,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=1
        k=tck
        do j=2,tcj-1
            ap(i,j,k)=ae(i,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=0!
    !!!!!!!!!!!!!!
        i=tci
        j=1
        do k=2,tck-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=b!
    !!!!!!!!!!!!!!
        i=tci
        j=tcj
        do k=2,tck-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=tci
        k=1
        do j=2,tcj-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=tci
        k=tck
        do j=2,tcj-1
            ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        j=1
        k=1
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        j=1
        k=tck
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j,k)*hTnew(i,j+1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=b!
    !		  z=0!
    !!!!!!!!!!!!!!
        j=tcj
        k=1
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k)*hTnew(i,j,k+1)+&
                ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno y=b!
    !		  z=c!
    !!!!!!!!!!!!!!
        j=tcj
        k=tck
        do i=2,tci-1
            ap(i,j,k)=ae(i,j,k)+ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
            hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+ae(i-1,j,k)*hTnew(i-1,j,k)+&
                an(i,j-1,k)*hTnew(i,j-1,k)+&
                at(i,j,k-1)*hTnew(i,j,k-1)+&
                ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
        enddo
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=1
        j=1
        k=1
        ap(i,j,k)=ae(i,j,k)+an(i,j,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=1
        j=1
        k=tck
        ap(i,j,k)=ae(i,j,k)+an(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=b!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=1
        j=tcj
        k=1
        ap(i,j,k)=ae(i,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=0!
    !		  y=b!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=1
        j=tcj
        k=tck
        ap(i,j,k)=ae(i,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i,j,k)*hTnew(i+1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=0!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=tci
        j=1
        k=1
        ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=0!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=tci
        j=1
        k=tck
        ap(i,j,k)=ae(i-1,j,k)+an(i,j,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j,k)*hTnew(i,j+1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=b!
    !		  z=0!
    !!!!!!!!!!!!!!
        i=tci
        j=tcj
        k=1
        ap(i,j,k)=ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k)*hTnew(i,j,k+1)+&
            ap0(i,j,k)*hTold(i,j,k))/ap(i,j,k)
    !!!!!!!!!!!!!!
    !Contorno x=a!
    !		  y=b!
    !		  z=c!
    !!!!!!!!!!!!!!
        i=tci
        j=tcj
        k=tck
        ap(i,j,k)=ae(i-1,j,k)+an(i,j-1,k)+at(i,j,k-1)+ap0(i,j,k)
        hTnew(i,j,k)=(ae(i-1,j,k)*hTnew(i-1,j,k)+&
            an(i,j-1,k)*hTnew(i,j-1,k)+&
            at(i,j,k-1)*hTnew(i,j,k-1)+&
            ap0(i,j,k)*hTold(i,j,k)+qest(ny*(i-1)+j)*dcx*dcy)/ap(i,j,k)
    !!!!!!!!!!!!!!!!!!!!!!!
    !Teste da Converg�ncia!
    !!!!!!!!!!!!!!!!!!!!!!!
        eps=sum((hTnew-hT)**2.d0)
    enddo
    hTold=hTnew
end subroutine

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