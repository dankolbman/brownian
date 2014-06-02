c       Forms a lattice structure for the given number of particles.
c       
c       Parameters:
c       ni = number of particle positions in x and y
c       lattx,latty = the lattice separation distance in x and y
c       disti = the distance of the first particle to be placed from the
c           corner in each direction
c       even = 1 if the row is even, 0 if odd
        
        subroutine formlattice()
        include 'param'
        integer i,j,nx,ny,npart,rem1,rem2,formseed
        real*8 lattx,latty,distx,disty,even,prob2
        real*8 ran3
        
        npart=npart1+npart2
        if (floor((dsqrt(dfloat(npart)))**2) .lt. npart) then
            write(*,*) "Warning: npart is not a square integer."
        endif
        rem1=npart1
        rem2=npart2
        prob2=dfloat(rem2)/dfloat(rem1+rem2)
        formseed=89
        
        latty=boxy/dsqrt(dfloat(npart))
        lattx=latty*dsqrt(3.0d0)/2.0d0
c        write(*,*) lattx,latty
        
        nx=int(dsqrt(dfloat(npart)))
        ny=nx
c        write(*,*) nx,ny
        
        even=0.0d0
        
        do i=1, nx
            do j=1,ny
                if (ran3(formseed) .gt. prob2) then
                    y1(npart1-rem1+1)=even*half*latty+dfloat(j-1)*latty
                    x1(npart1-rem1+1)=dfloat(i-1)*lattx
                    rem1=rem1-1
                else
                    y2(npart2-rem2+1)=even*half*latty+dfloat(j-1)*latty
                    x2(npart2-rem2+1)=dfloat(i-1)*lattx
                    rem2=rem2-1
                endif
                
                if (rem1+rem2 .gt. 0) then
                    prob2=dfloat(rem2)/dfloat(rem1+rem2)
                    formseed=formseed+15
                else
                    prob2=0.0d0
                endif
            enddo
            even=1.0d0-even
        enddo
        
 10     return
        end
        
c-----------------------------------------------------------------------

c       Generates a random number

        REAL*8 FUNCTION ran3(idum)
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
        ran3=min(AM*iy,RNMX)
        return
        END
