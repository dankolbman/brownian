c       Generates pair distribution function histograms.
c       
c       Parameters:
c       -x1orig,y1orig = origin coordinates for the pair distribution
c           function
c       -dx,dy,dr = x, y, and radial distances between particles
c       -nideal = a normalization constant from the 1-1/2-2 lattice
c       -i,j = loop indices
c       -indexr,indexx,indexy = distribution histogram indices
c       -histrAA,histxAA,histyAA = local histograms between species 1
c       -histxyAA = local histograms
        
        subroutine grcalc()
        include 'param'
        real*8 x1orig,y1orig,x2orig,y2orig,dx,dy,dr,nideal
        integer i,j,indexr,indexx,indexy
        real*8 histr11(histbinmax),histr22(histbinmax),
     &      histr12(histbinmax)
        
c       Clear the temporary histograms.
        do i=1,histbinmax
            histr11(i)=0.0d0
            histr22(i)=0.0d0
            histr12(i)=0.0d0
        enddo
        
c       For 1-1 particle distributions:
        do i=1,npart1-1
c           Set the origin.
            x1orig=x1(i)
            y1orig=y1(i)
            if(dsqrt(x1orig*x1orig+y1orig*y1orig) .lt. radius/2.0) then
              do j=i+1,npart1
c               Get the distance from the origin of each particle.
                  dx=x1orig-x1(j)
                  dy=y1orig-y1(j)
!                write(*,*) i,j,x1orig,x1(j)
c                 Do periodic boundaries
                  if(circ .eq. 0) then
                    dx=dx-boxx*anint(dx*invboxx)
                    dy=dy-boxy*anint(dy*invboxy)
                  endif
                  dr=dsqrt(dx**2+dy**2)
c                 Convert the radial distances into an array index.
                  indexr=int(dr*invdelta)+1
c                 Add to the histograms.
                  if (indexr .le. histbinmax .and. dr .lt. radius/2.0)
     &              histr11(indexr)=histr11(indexr)+2.0
              enddo
            endif
        enddo
        
c       For 2-2 particle distributions:
        do i=1,npart2-1
c           Set the origin.
            x2orig=x2(i)
            y2orig=y2(i)
            if(dsqrt(x2orig*x2orig+y2orig*y2orig) .lt. radius/2.0) then
              do j=i+1,npart2
c               Get the distance from the origin of each particle.
                  dx=x2orig-x2(j)
                  dy=y2orig-y2(j)
                  if(circ .eq. 0) then
                    dx=dx-boxx*anint(dx*invboxx)
                    dy=dy-boxy*anint(dy*invboxy)
                  endif
                  dr=dsqrt(dx**2+dy**2)
c               Convert the radial distances into an array index.
                  indexr=int(dr*invdelta)+1
c               Add to the histograms.
                  if (indexr .le. histbinmax .and. dr .lt. radius/2.0)
     &              histr22(indexr)=histr22(indexr)+2.0
              enddo
            endif
        enddo
        
c       For 1-2 particle distributions:
        do i=1,npart1
c           Set the origin.
            x1orig=x1(i)
            y1orig=y1(i)
            
            do j=1,npart2
c               Get the distance from the origin of each particle.
                dx=x1orig-x2(j)
                dy=y1orig-y2(j)
                if(circ .eq. 0) then
                  dx=dx-boxx*anint(dx*invboxx)
                  dy=dy-boxy*anint(dy*invboxy)
                endif
                dr=dsqrt(dx**2+dy**2)
c               Convert the radial distances into an array index.
                indexr=int(dr*invdelta)+1
c               Add to the histograms.
                if (indexr .le. histbinmax)
     &              histr12(indexr)=histr12(indexr)+2.0
            enddo
        enddo
        
c       Normalize the distributions to the ideal lattice displacement.
c        nideal=1.0/(2.0*sqrt(3.0))
        nideal=3.0/2.0
        do i=1,histbinmax
            histr11(i)=histr11(i)/(nideal*pi*delta*delta*float(2*i-1))
            gr11(i)=gr11(i)+histr11(i)/float(npart1)
            histr22(i)=histr22(i)/(nideal*pi*delta*delta*float(2*i-1))
            gr22(i)=gr22(i)+histr22(i)/float(npart2)
            histr12(i)=histr12(i)/(nideal*pi*delta*delta*float(2*i-1))
            gr12(i)=gr12(i)+histr12(i)/float(npart2)
        enddo
        
        return
        end
