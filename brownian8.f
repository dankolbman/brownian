c       Simulates 2D Brownian motion in a binary system of active,
c       interactive spheres susceptible to thermal noise.
c       
c       Motion evolves according to v(t + dt) = vp + f(t) + nu(t)
c           vp :: self-propulsion speed
c           f(t) = D * F(t) / kT :: drift velocity
c           nu(t) = sqrt(2 * D * dt) * R(t) :: Langevin noise
c       dt is the time step, vp is the particle's self-propulsion speed,
c       D is the diffusivity, F is the sum of the external forces on a 
c       particle, kT is the thermal energy unit, and R(t) is a normally-
c       distributed random generator.
c       
c       Parameters:
c       nequil = the number of movements to equilibriate the system
c       ncor = the number of movements between data-taking
c       nitn = the number of trials to simulate
c       nrun = the number of runs of data-taking per trial
c       lattice = whether to start in a lattice (0), a displaced lattice
c           (1), or a random configuration (2)
c       time = a time stamp
c       loadingi = a progress indicator
c       i = loop indices
c       raoldi = the previous particle positions
c       meansqrdisi = the mean square displacement as it evolves for
c           each trial
c       sumsqrdisi = the sum square displacement for each run
c       vaavei = the average velocity for each run
c       vaveruni = the average speed as it evolves for each trial
c       deltaai = the displacement between runs
c       outfile,istring = strings for file naming
        
        program brownian8
        include 'param'
        integer time(3)
        integer loading1,loading2
        integer i,j,k
c        real*8 rxold1(npart1),ryold1(npart1)
        real*8 meansqrdis1,sumsqrdis1,meansqrdis2,sumsqrdis2
        real*8 msdave1(nrun),msdave2(nrun)
        real*8 vxave1,vyave1,vaverun1,vxave2,vyave2,vaverun2
        real*8 deltax1,deltay1,deltax2,deltay2
        character outfile*15,istring*3
        
        call itime(time)
        write (*,"('Time: ',i2,':',i2,':',i2)") time
        
c       Open data files.
        open(unit=10,file='msd1.dat',status='unknown')
        open(unit=11,file='msd2.dat',status='unknown')
        open(unit=12,file='msdave1.dat',status='unknown')
        open(unit=13,file='msdave2.dat',status='unknown')
        open(unit=14,file='gr11.dat',status='unknown')
        open(unit=15,file='gr22.dat',status='unknown')
        open(unit=16,file='gr12.dat',status='unknown')
        open(unit=17,file='gr11ave.dat',status='unknown')
        open(unit=18,file='gr22ave.dat',status='unknown')
        open(unit=19,file='gr12ave.dat',status='unknown')
        open(unit=30,file='vave1.dat',status='unknown')
        open(unit=31,file='vave2.dat',status='unknown')
        
        do i=1,nitn
            write(istring,'(i3)') i
            outfile='itn1pos'//trim(adjustl(istring))//'.dat'
            open(unit=500+i,file=outfile,status='unknown')
            outfile='itn2pos'//trim(adjustl(istring))//'.dat'
            open(unit=700+i,file=outfile,status='unknown')
            outfile='itn1vel'//trim(adjustl(istring))//'.dat'
            open(unit=900+i,file=outfile,status='unknown')
            outfile='itn2vel'//trim(adjustl(istring))//'.dat'
            open(unit=1100+i,file=outfile,status='unknown')
            outfile='config1'//trim(adjustl(istring))//'.dat'
            open(unit=1300+i,file=outfile,status='unknown')
            outfile='config2'//trim(adjustl(istring))//'.dat'
            open(unit=1500+i,file=outfile,status='unknown')
        enddo
        
        do i=1,nrun
            msdave1(i)=0.0d0
            msdave2(i)=0.0d0
        enddo
        
c       Establish system parameters.
        call scale()
        
        do i=1,nitn
c           Arrange particles.
            call readcon()
            do j=1,nequil
                call movepart(1)
c               Run loading bar.
                loading1=20*j/nequil
                loading2=20*(j-1)/nequil
                if (loading1 .gt. loading2) then
                    write(*,'(a4,t6,i3,t10,a1)') 'Equ: ',5*loading1,'%'
                endif
            enddo
            
c           Reset data for every trial.
            meansqrdis1=0.0d0
            meansqrdis2=0.0d0
            call grinit()
            
 10         format(a5,t8,a3,t15,a4,t30,a24)
 20         format(a5,t8,a3,t15,a4,t30,a5,t50,a5,t70,a7)
            write(10,10)'% itn','run','time','mean square displacement'
            write(11,10)'% itn','run','time','mean square displacement'
            write(30,20)'% itn','run','time','vxave','vyave','vaverun'
            write(31,20)'% itn','run','time','vxave','vyave','vaverun'
            
            do j=1,npart1
                x01(j)=x1(j)
                y01(j)=y1(j)
                counter1x(j)=0
                counter1y(j)=0
                write(1300+i,*) x01(j),y01(j)
            enddo
            do j=1,npart2
                x02(j)=x2(j)
                y02(j)=y2(j)
                counter2x(j)=0
                counter2y(j)=0
                write(1500+i,*) x02(j),y02(j)
            enddo
            write(1300+i,*)
            write(1300+i,*)
            write(1500+i,*)
            write(1500+i,*)
            
            do j=1,nrun
                write(500+i,20)'% run','par','time','x','y'
                write(700+i,20)'% run','par','time','x','y'
                write(900+i,20)'% run','par','time','vx','vy'
                write(1100+i,20)'% run','par','time','vx','vy'
                
c               Run loading bar
                loading1=(20*((i-1)*nrun+j))/(nrun*nitn)
                loading2=(20*((i-1)*nrun+j-1))/(nrun*nitn)
                if (loading1 .gt. loading2) then
                    write(*,'(a4,t6,i3,t10,a1)') 'Run: ',5*loading1,'%'
                endif
                
c               Reset data.
                sumsqrdis1=0.0d0
                sumsqrdis2=0.0d0
                vxave1=0.0d0
                vyave1=0.0d0
                vaverun1=0.0d0
                vxave2=0.0d0
                vyave2=0.0d0
                vaverun2=0.0d0
                
c               Move particles.
                do k=1,ncor
                    call movepart(0)
                enddo
                
c               Produce radial distribution histogram data.
                call grcalc()
                
                do k=1,npart1
c                   Calculate statistics of motion.
                    if( circ .eq. 1) then
c                     Don't need to worry about wrap in circular bounds
                      deltax1=x1(k)-x01(k)
                      deltay1=y1(k)-y01(k)
                    else
                      deltax1=counter1x(k)*boxx+x1(k)-x01(k)
                      deltay1=counter1y(k)*boxy+y1(k)-y01(k)
                    endif
                    sumsqrdis1=sumsqrdis1+deltax1*deltax1
                    sumsqrdis1=sumsqrdis1+deltay1*deltay1
                    vxave1=vxave1+vx1(k)
                    vyave1=vyave1+vy1(k)
c                   Write position and velocity data per run.
                    write(500+i,30) j,k,dfloat(j)*ncor*dt,x1(k),y1(k)
                    write(900+i,30) j,k,dfloat(j)*ncor*dt,vx1(k),vy1(k)
                enddo
                vxave1=vxave1/dfloat(npart1)
                vyave1=vyave1/dfloat(npart1)
                vaverun1=dsqrt(vxave1*vxave1+vyave1*vyave1)
                
                do k=1,npart2
c                   Calculate statistics of motion.
                    if( circ .eq. 1) then
c                     Don't need to worry about wrap in circular bounds
                      deltax2=x2(k)-x02(k)
                      deltay2=y2(k)-y02(k)
                    else
                      deltax2=counter2x(k)*boxx+x2(k)-x02(k)
                      deltay2=counter2y(k)*boxy+y2(k)-y02(k)
                    endif
                    sumsqrdis2=sumsqrdis2+deltax2*deltax2
                    sumsqrdis2=sumsqrdis2+deltay2*deltay2
                    vxave2=vxave2+vx2(k)
                    vyave2=vyave2+vy2(k)
c                   Write position and velocity data per run.
                    write(700+i,30) j,k,dfloat(j)*ncor*dt,x2(k),y2(k)
                    write(1100+i,30) j,k,dfloat(j)*ncor*dt,vx2(k),vy2(k)
                enddo
                vxave2=vxave2/dfloat(npart2)
                vyave2=vyave2/dfloat(npart2)
                vaverun2=dsqrt(vxave2*vxave2+vyave2*vyave2)
                
 30             format(i3,t8,i3,t14,e12.5,t29,e17.10,t49,e17.10,t69,
     &              e17.10)
                write(500+i,*)
                write(500+i,*)
                write(700+i,*)
                write(700+i,*)
                write(900+i,*)
                write(900+i,*)
                write(1100+i,*)
                write(1100+i,*)
                
c               Calculate mean square displacement.
                meansqrdis1=sumsqrdis1/dfloat(npart1*ncor)
                msdave1(j)=msdave1(j)+meansqrdis1
                write(10,50) i,j,dfloat(j)*ncor*dt,meansqrdis1
                meansqrdis2=sumsqrdis2/dfloat(npart2*ncor)
                msdave2(j)=msdave2(j)+meansqrdis2
                write(11,50) i,j,dfloat(j)*ncor*dt,meansqrdis2
 50             format(i3,t8,i3,t14,e12.5,t29,e12.5)
                
c               Write average velocity data.
                write(30,30) i,j,dfloat(j)*ncor*dt,vxave1,vyave1,
     &              vaverun1
                write(31,30) i,j,dfloat(j)*ncor*dt,vxave2,vyave2,
     &              vaverun2
                
c               Write radial distribution data.
                write(14,20)'% itn','run','time','r','g'
                write(15,20)'% itn','run','time','r','g'
                write(16,20)'% itn','run','time','r','g'
                do k=1,histbinmax
                    write(14,30) i,j,(dfloat(j)*ncor*dt),(dfloat(k)*
     &                  delta),gr11(k)
                    write(15,30) i,j,(dfloat(j)*ncor*dt),(dfloat(k)*
     &                  delta),gr22(k)
                    write(16,30) i,j,(dfloat(j)*ncor*dt),(dfloat(k)*
     &                  delta),gr12(k)
                enddo
                
                write(14,*)
                write(14,*)
                write(15,*)
                write(15,*)
                write(16,*)
                write(16,*)
            enddo
            
            write(10,*)
            write(10,*)
            write(11,*)
            write(11,*)
            write(30,*)
            write(30,*)
            write(31,*)
            write(31,*)
        enddo
        
        do i=1,nrun
            msdave1(i)=msdave1(i)/nitn
            msdave2(i)=msdave2(i)/nitn
            write(12,50) 0,i,dfloat(i)*ncor*dt,msdave1(i)
            write(13,50) 0,i,dfloat(i)*ncor*dt,msdave2(i)
        enddo
        
c       Close files.
        close(10)
        close(11)
        close(14)
        close(15)
        close(16)
        close(30)
        close(31)
        do i=1,nitn
            close(500+i)
            close(700+i)
            close(900+i)
            close(1100+i)
            close(1300+i)
            close(1500+i)
        enddo
        
        call itime(time)
        write (*,"('Time: ',i2,':',i2,':',i2)") time
        
        stop
        end
