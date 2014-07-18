c       Moves particles according to overdamped Langevin dynamics,
c       accounting for active forces.
c       
c       Parameters:
c       i = loop index
c       normabc = Gaussian-distributed numbers
c       equil = 1 when equilibriating, 0 otherwise
c       
c       Functions:
c       ran1 = a uniformly-distributed random number generator
c       gasdev = a normally-distributed random number generator
        
        subroutine movepart(equil)
        include 'param'
        integer equil
        integer i
        real*8 norma,normb,normc,phi,ranr,angle
        real*8 ran1,gasdev
        
        call forcecalc()
c        write(*,*) vprop1,vprop2
c         write(*,*) pretrad,prerotd
        
        do i=1,npart1
c           Get random numbers.
            norma=gasdev(seed1a)
            normb=gasdev(seed1b)
            normc=gasdev(seed1c)
            
c           Calculate particle velocity.
            theta1(i)=theta1(i)+prerotd*norma
            
c           Here, during equilibriation, equil=1, and during motion-
c           tracking, equil=0. Therefore, self-propulsion and thermal 
c           fluctuations do not contribute during equilibriation.
c            vx1(i)=dfloat(1-equil)*vprop1*dcos(theta1(i))+fx1(i)+
c     &          dfloat(1-equil)*pretrad*normb
c            vy1(i)=dfloat(1-equil)*vprop1*dsin(theta1(i))+fy1(i)+
c     &          dfloat(1-equil)*pretrad*normc
            vx1(i)=vprop1*dcos(theta1(i))+fx1(i)+
     &          pretrad*normb
            vy1(i)=vprop1*dsin(theta1(i))+fy1(i)+
     &          pretrad*normc

c           Circular conditions
            if(circ .eq. 1) then
c             Propose new coords
              newx=x1(i)+dt*vx1(i)
              newy=y1(i)+dt*vy1(i)
c             Accept move if within the bounds
              if( (newx**2 + newy**2) .lt. (radius-dia*0.5d0)**2 ) then
                x1(i)=newx
                y1(i)=newy
              else
                angle = atan2(newy,newx)
                x1(i)=(radius-dia/2.0)*cos(angle)
                y1(i)=(radius-dia/2.0)*sin(angle)
              endif
            else
c             Shift each particle.
              x1(1)=x1(i)+dt*vx1(i)
              y1(1)=y1(i)+dt*vy1(i)
              counter1x(i)=counter1x(i) + floor(x1(i)*invboxx)
              counter1y(i)=counter1y(i) + floor(y1(i)*invboxy)
c             Adjust particle position to within the box dimensions.c
              x1(i)=x1(i)-boxx*floor(x1(i)*invboxx)
              y1(i)=y1(i)-boxy*floor(y1(i)*invboxy)
            endif
        enddo
        
        do i=1,npart2
c           Get random numbers.
            norma=gasdev(seed2a)
            normb=gasdev(seed2b)
            normc=gasdev(seed2c)
            
c           Calculate particle velocity.
            theta2(i)=theta2(i)+prerotd*norma
            
c           Here, during equilibriation, equil=1, and during motion-
c           tracking, equil=0. Therefore, self-propulsion and thermal 
c           fluctuations do not contribute during equilibriation.
c            vx2(i)=dfloat(1-equil)*vprop2*dcos(theta2(i))+fx2(i)+
c     &          dfloat(1-equil)*pretrad*normb
c            vy2(i)=dfloat(1-equil)*vprop2*dsin(theta2(i))+fy2(i)+
c     &          dfloat(1-equil)*pretrad*normc
            vx2(i)=vprop2*dcos(theta2(i))+fx2(i)+
     &          pretrad*normb
            vy2(i)=vprop2*dsin(theta2(i))+fy2(i)+
     &          pretrad*normc

c           Circular conditions
            if(circ .eq. 1) then
c             Propose new coords
              newx=x2(i)+dt*vx2(i)
              newy=y2(i)+dt*vy2(i)
c             Accept move if within the bounds
              if( (newx**2 + newy**2) .lt. (radius-dia*0.5d0)**2 ) then
                x2(i)=newx
                y2(i)=newy
              else
                angle = atan2(newy,newx)
                x2(i)=(radius-dia/2.0)*cos(angle)
                y2(i)=(radius-dia/2.0)*sin(angle)
              endif
            else
c             Shift each particle.
              x2(1)=x2(i)+dt*vx2(i)
              y2(1)=y2(i)+dt*vy2(i)
              counter2x(i)=counter2x(i) + floor(x2(i)*invboxx)
              counter2y(i)=counter2y(i) + floor(y2(i)*invboxy)
c             Adjust particle position to within the box dimensions.c
              x2(i)=x2(i)-boxx*floor(x2(i)*invboxx)
              y2(i)=y2(i)-boxy*floor(y2(i)*invboxy)
            endif
        enddo
        
        return
        end
        
c----------------------------------------------------------------------
        
        FUNCTION gasdev(idum)
        INTEGER idum
        REAL*8 gasdev, ran1
        INTEGER iset
        REAL fac,gset,rsq,v1,v2
        SAVE iset,gset
        DATA iset/0/
        if (iset.eq.0) then
 1          v1=2.*ran1(idum)-1.
            v2=2.*ran1(idum)-1.
            rsq=v1**2+v2**2
            if(rsq.ge.1..or.rsq.eq.0.)goto 1
            fac=sqrt(-2.*log(rsq)/rsq)
            gset=v1*fac
            gasdev=v2*fac
            iset=1
        else
            gasdev=gset
            iset=0
        endif
        return
        END
c       (C) Copr. 1986-92 Numerical Recipes Software .

c-----------------------------------------------------------------------

c       Generates a random number

        REAL*8 FUNCTION ran1(idum)
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
        ran1=min(AM*iy,RNMX)
        return
        END
