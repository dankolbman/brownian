c       Clears old data from calculations of the pair distribution.
c       
c       Parameters:
c       -i,j = loop indices

        subroutine grinit()
        include 'param'
        integer i,j
        
c       Set the bin width.
        delta=1.0/20.0
        if(circ .eq. 1) then
          delta=1.0/(radius*2.0)
        endif
        invdelta=1.0/delta
        
c       Clear the radial distribution histograms.
        do i=1,histbinmax
            gr11(i)=0.0d0
            gr22(i)=0.0d0
            gr12(i)=0.0d0
        enddo
        
        return
        end
