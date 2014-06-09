c       Arranges particles within the system.
c       
c       Parameters:
c       i = loop index
c       lattice = whether to start in a lattice (0), a displaced lattice
c           (1), or a random configuration (2)
c       
c       Functions:
c       ran2 = a uniformly-distributed random number generator
        
        subroutine readcon()
        include 'param'
        integer i
        real*8 sepx,sepy,alt,variance
        real*8 r,phi
        real*8 ran2
        
        variance=0.05d0
        
        if (lattice .eq. 2) then ! Random configuration
c            Generate the initial positions of x and y.
            do i=1,npart1
                if( circ .eq. 1) then ! Circular bounds
                  r=dsqrt(mod(ran2(seedcon),1.0))*radius
                  phi=ran2(seedcon)*2*pi
                  
                  x1(i)=r*cos(phi)
                  y1(I)=r*sin(phi)
                else
                  x1(i)=ran2(seedcon)
                  y1(i)=ran2(seedcon)
                  
                  x1(i)=x1(i)*boxx
                  y1(i)=y1(i)*boxy
                endif
            enddo

            do i=1,npart2
                if( circ .eq. 1) then ! Circular bonuds
                  r=dsqrt(mod(ran2(seedcon),1.0))*radius
                  phi=ran2(seedcon)*2*pi
                  
                  x2(i)=r*cos(phi)
                  y2(I)=r*sin(phi)
                else 
                  x2(i)=ran2(seedcon)
c                seedcon=seedcon+31
                  y2(i)=ran2(seedcon)
c                seedcon=seedcon+31
                  
                  x2(i)=x2(i)*boxx
                  y2(i)=y2(i)*boxy
                endif
            enddo
        else if (lattice .eq. 1) then ! Displaced lattice configuration.
            call formlattice()
            
            do i=1,npart1
                x1(i)=x1(i)+variance*boxx*(2*ran2(seedcon)-1.0d0)
                y1(i)=y1(i)+variance*boxy*(2*ran2(seedcon)-1.0d0)
                
                x1(i)=x1(i)-boxx*floor(x1(i)*invboxx)
                y1(i)=y1(i)-boxy*floor(y1(i)*invboxy)
            enddo
            
            do i=1,npart2
                x2(i)=x2(i)+variance*boxx*(2*ran2(seedcon)-1.0d0)
                y2(i)=y2(i)+variance*boxy*(2*ran2(seedcon)-1.0d0)
                
                x2(i)=x2(i)-boxx*floor(x2(i)*invboxx)
                y2(i)=y2(i)-boxy*floor(y2(i)*invboxy)
            enddo
        else ! Lattice configuration.
            if(circ .eq. 1) then
              call circlattice()
            else
              call formlattice()
            endif
       endif
        
        do i=1,npart1
            vx1(i)=0.0d0
            vy1(i)=0.0d0
            theta1(i)=0.0d0
        enddo

        do i=1,npart2
            vx2(i)=0.0d0
            vy2(i)=0.0d0
            theta2(i)=0.0d0
        enddo
        
        return
        end
        
c-----------------------------------------------------------------------

c       Generates a random number

        REAL*8 FUNCTION ran2(idum)
        INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
        REAL AM,EPS,RNMX
        PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     &      NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv(NTAB),iy
        SAVE iv,iy
        DATA iv /NTAB*0/, iy /0/
        if (idum.le.0.or.iy.eq.0) then
            idum=max(-idum,1)
            do 11 j=NTAB+8,1,-1
                k=idum/IQ
                idum=IA*(idum-k*IQ)-IR*k
                if (idum.lt.0) idum=idum+IM
                if (j.le.NTAB) iv(j)=idum
11          continue
            iy=iv(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ran2=min(AM*iy,RNMX)
        return
        END
