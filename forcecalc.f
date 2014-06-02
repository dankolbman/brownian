c       Calculates all the pairwise forces.
c       
c       Parameters:
c       i = loop indices
c       da = the separation between particles
c       frange = the interaction cutoff
c       f*** = the force/separation distance between two particles
c       netforces = the sum of interaction forces between two particles
c       
c       Functions:
c       forcescp = calculates the position-dependence of the SCP
c       forceljp = calculates the position-dependence of the LJP
        
        subroutine forcecalc()
        include 'param'
        integer i,j
        real*8 dx,dy,dr
        real*8 frange
        real*8 fscp,fljp,frep,fadh,netforces
        real*8 forcescp,forceljp,forcerep,forceadh
        
        frange=6.0d0
        
c       Reset all forces.
        do i=1,npart1
            fx1(i)=0.0d0
            fy1(i)=0.0d0
        enddo
        do i=1,npart2
            fx2(i)=0.0d0
            fy2(i)=0.0d0
        enddo
        
c       For 1-1 particle interactions...
c       ...for each particle i, from each particle j...
        do i=1,npart1-1
            do j=i+1,npart1
c               ...find the distance between the particles.
                dx=x1(i)-x1(j)
                dx=dx-boxx*anint(dx*invboxx)
                dy=y1(i)-y1(j)
                dy=dy-boxy*anint(dy*invboxy)
                dr=dsqrt(dx*dx+dy*dy)
c               For all particles within a given proximity...
                if (dr .le. frange) then
c                   ...find the forces.
                    fscp=prescp1*forcescp(kappa,dr)
                    fljp=preljp1*forceljp(dia,dr)
                    frep=prerep1*forcerep(dia,dr)
                    fadh=preadh1*forceadh(dia,contact,dr)
                    netforces=fscp+fljp+frep+fadh
                    
                    fx1(i)=fx1(i)+dx*netforces
                    fx1(j)=fx1(j)-dx*netforces
                    fy1(i)=fy1(i)+dy*netforces
                    fy1(j)=fy1(j)-dy*netforces
                endif
            enddo
        enddo
        
c       For 2-2 particle interactions...
c       ...for each particle i, from each particle j...
        do i=1,npart2-1
            do j=i+1,npart2
c               ...find the distance between the particles.
                dx=x2(i)-x2(j)
                dx=dx-boxx*anint(dx*invboxx)
                dy=y2(i)-y2(j)
                dy=dy-boxy*anint(dy*invboxy)
                dr=dsqrt(dx*dx+dy*dy)
c               For all particles within a given proximity...
                if (dr .le. frange) then
c                   ...find the forces.
                    fscp=prescp2*forcescp(kappa,dr)
                    fljp=preljp2*forceljp(dia,dr)
                    frep=prerep2*forcerep(dia,dr)
                    fadh=preadh2*forceadh(dia,contact,dr)
                    netforces=fscp+fljp+frep+fadh
                    
                    fx2(i)=fx2(i)+dx*netforces
                    fx2(j)=fx2(j)-dx*netforces
                    fy2(i)=fy2(i)+dy*netforces
                    fy2(j)=fy2(j)-dy*netforces
                endif
            enddo
        enddo
        
c       For 1-2 particle interactions...
c       ...for each particle i, from each particle j...
        do i=1,npart1
            do j=1,npart2
c               ...find the distance between the particles.
                dx=x1(i)-x2(j)
                dx=dx-boxx*anint(dx*invboxx)
                dy=y1(i)-y2(j)
                dy=dy-boxy*anint(dy*invboxy)
                dr=dsqrt(dx*dx+dy*dy)
c               For all particles within a given proximity...
                if (dr .le. frange) then
c                   ...find the forces.
                    fscp=prescp12*forcescp(kappa,dr)
                    fljp=preljp12*forceljp(dia,dr)
                    frep=prerep12*forcerep(dia,dr)
                    fadh=preadh12*forceadh(dia,contact,dr)
                    netforces=fscp+fljp+frep+fadh
                    
                    fx1(i)=fx1(i)+dx*netforces
                    fx2(j)=fx2(j)-dx*netforces
                    fy1(i)=fy1(i)+dy*netforces
                    fy2(j)=fy2(j)-dy*netforces
                endif
            enddo
        enddo
        
        return
        end
        
c-----------------------------------------------------------------------

        function forcescp(kappa,dr)
        real*8 forcescp
        real*8 kappa,dr
        
        forcescp=(1.0d0/dr**3)*dexp(-kappa*dr)*(1.0d0+kappa*dr)
        
        return
        end
        
c-----------------------------------------------------------------------
        
        function forceljp(dia,dr)
        real*8 forceljp
        real*8 dia,dr
        
        forceljp=(1.0d0/dr)*(2.0d0*(dia/dr)**13-(dia/dr)**7)
        
        return
        end
        
c-----------------------------------------------------------------------
        
        function forcerep(dia,dr)
        real*8 forcerep
        real*8 dia,dr
        real*8 x
        
        x=1.0d0-(dia/dr)
        forcerep=(dabs(x)-x)/2.0d0
        
        return
        end
c-----------------------------------------------------------------------
        
        function forceadh(dia,contact,dr)
        real*8 forceadh
        real*8 dia,contact,dr
        real*8 x
        
        !* False! This is the potential energy, not the force.
        x=dabs(dr-dia)-contact
        forceadh=(x-dabs(x))/2.0d0/dr
        
        return
        end
