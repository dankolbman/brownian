c       Establishes system parameters and de-dimensionalizes them.
c       
c       Parameters:
c       npart = the total number of particles in the system
c       dielec = the relative permittivity of the system


c       boltz = Boltzman's constant
c       elcharge = the elementary charge of an electron
c       temp = the temperature of the system
c       eta = the viscosity of the system
c       packing = the packing fraction of the particles in the system
c       kapfact = a factor in the screened Coulomb potential
c       charge1 = the charge of a particle (in elementary charges)
c       epsilon1 = the strength of the Lennard-Jones potential
c       repul1 = the strength of the Fily-Marchetti linear repulsion
c       unittime = the simulation time scale
c       unitlength = the simulation length scale
c       unitenergy = the simulation energy scale
        
        subroutine scale()
        include 'param'
        integer npart
        real*8 dielec,boltz,elcharge
        real*8 temp,eta,packing,kapfact
        real*8 charge1,epsilon1,repul1,adhesion1
        real*8 charge2,epsilon2,repul2,adhesion2
        real*8 charge12,epsilon12,repul12,adhesion12
        real*8 unittime,unitlength,unitenergy
        
c       Set constants (in CGS units).
        boltz=1.38d-16 ! erg / K
        elcharge=4.803d-10 ! statC
        
c       Set experimental parameters (in CGS units).
        dt=1.0d-6 ! s
        temp=298.0d0 ! K
        dia=1.07d-4 ! cm
        eta=1.0d-2 ! g/(cm*s)
        packing=0.40d0 ! %
        npart=npart1+npart2
        seed1a=-74300007
        seed1b=-25255273
        seed1c=-83647528
        seed2a=-62534449
        seed2b=-111827837
        seed2c=-92733893 
        seedcon=-46766782
        lattice=0 ! 0 = lattice, 1 = displaced lattice, 2 = random
        circ=1  ! 0 = box bounds, 1 = circular bounds

        diffus=boltz*temp/(3.0d0*pi*eta*dia) ! cm**2/s
        rotdiffus=500*boltz*temp/(pi*eta*dia*dia*dia) ! 1/s !*
c       Set interaction parameters (dimensionless).
c       Screened Coulomb Potential
        kappa=0.5d0 ! 1/cm
        dielec=78.0d0 ! (unitless)
        charge1=0.0d8 ! (elementary charges)
        charge2=0.0d8
        charge12=0.0d8
c       Lennard-Jones Potential
        epsilon1=0.0d0 ! (thermal energy unit)
        epsilon2=0.0d0
        epsilon12=0.0d0
c       Hookean Contact Repulsion (Soft <= ~0.01; Hard >= ~0.1)
        repul1=0.0 !1.0d-2 !0.31623d-2 ! (thermal energy unit / particle length unit)
        repul2=0.0 !repul1
        repul12=0.0 !2*repul1*repul2/(repul1+repul2) !*
c       Contact Adhesion Force
        contact=0.05*dia
        adhesion1=0.0d0 ! (thermal energy unit / particle length unit)
        adhesion2=0.0d0
        adhesion12=0.0d0
c       Self-Propulsion Speed
        vprop1=0.0 !1.0d2 ! (particle length unit / diffusion time unit)
        vprop2=0.0 !vprop1*0.50
        
c       Set dimensionless units.
        unittime=dia*dia/diffus
        unitlength=dia
        unitenergy=boltz*temp
        
c       Nondimensionalize parameters.
        diffus=diffus*unittime/(unitlength*unitlength)
        dia=dia/unitlength
        dt=dt/unittime
        kappa=kappa*unitlength
        contact=contact/unitlength
        
c       Calculate force coefficients (dimensionless).
        pretrad=dsqrt(2.0d0*diffus/dt)
        prerotd=dsqrt(2.0d0*rotdiffus*dt)
        
        kapfact=kappa*half*dia
        prescp1=dexp(kapfact)/(1.0d0+kapfact)
        prescp1=(charge1*prescp1*elcharge)**2/dielec/boltz/temp
        prescp1=prescp1*unitlength
        prescp2=dexp(kapfact)/(1.0d0+kapfact)
        prescp2=(charge2*prescp2*elcharge)**2/dielec/boltz/temp
        prescp2=prescp2*unitlength
        prescp12=dexp(kapfact)/(1.0d0+kapfact)
        prescp12=(charge12*prescp12*elcharge)**2/dielec/boltz/temp
        prescp12=prescp12*unitlength
        
        preljp1=24.0d0*epsilon1
        preljp2=24.0d0*epsilon2
        preljp12=24.0d0*epsilon12
        
        prerep1=repul1/unitlength
        prerep2=repul2/unitlength
        prerep12=repul12/unitlength
        
        preadh1=adhesion1/contact
        preadh2=adhesion2/contact
        preadh12=adhesion12/contact
        
c       Set box dimensions (dimensionless).
        boxy=dsqrt((pi*dia**2*npart)/(2.0d0*dsqrt(3.0d0)*packing))
        boxx=boxy*half*dsqrt(3.0d0)
c       Determine circle radius for circular bounds
        radius=dsqrt( (dia/2)**2*npart/packing )
        
        invboxy=1.0d0/boxy
        invboxx=1.0d0/boxx
        
c       Record parameters.
        open(unit=1,file='sysparam.dat',status='unknown')
        write(1,*) 'Technical Parameters:'
        write(1,'(t5,a10,t34,i6)') 'nequil = ',nequil
        write(1,'(t5,a10,t34,i6)') 'ncor = ',ncor
        write(1,'(t5,a10,t34,i6)') 'nitn = ',nitn
        write(1,'(t5,a10,t34,i6)') 'nrun = ',nrun
        write(1,*)
        write(1,*) 'Experimental Parameters (de-dimensionalized):'
        write(1,'(t5,a14,t34,e17.10)') 'Temperature = ',temp
        write(1,'(t5,a12,t34,e17.10)') 'Diffusion = ',diffus
        write(1,'(t5,a20,t34,e17.10)') 'Angular Diffusion = ',rotdiffus
        write(1,'(t5,a12,t34,e17.10)') 'Time Step = ',dt
        write(1,'(t5,a15,t34,e17.10)') 'Packing Frac = ',packing
        write(1,'(t5,a12,t34,e17.10)') 'Box Width = ',boxx
        write(1,'(t5,a13,t34,e17.10)') 'Box Height = ',boxy
        write(1,'(t5,a13,t34,e17.10)') 'Circle Rad = ',radius
        write(1,'(t5,a10,t34,i6)') 'Circ = ',circ
        write(1,'(t5,a15,t34,e17.10)') 'Packing Frac = ',packing
        write(1,'(t5,a12,t34,e17.10)') 'Time Unit = ',unittime
        write(1,'(t5,a14,t34,e17.10)') 'Length Unit = ',unitlength
        write(1,'(t5,a14,t34,e17.10)') 'Energy Unit = ',unitenergy
        write(1,*)
        write(1,*) 'Particle Parameters (de-dimensionalized):'
        write(1,'(t5,a11,t34,e17.10)') 'Diameter = ',dia
        write(1,'(t5,a9,t34,i4)') 'Npart1 = ',npart1
        write(1,'(t5,a9,t34,i4)') 'Npart2 = ',npart2
        write(1,'(t5,a10,t34,e17.10)') 'Charge1 = ',charge1
        write(1,'(t5,a10,t34,e17.10)') 'Charge2 = ',charge2
        write(1,'(t5,a11,t34,e17.10)') 'Epsilon1 = ',epsilon1
        write(1,'(t5,a11,t34,e17.10)') 'Epsilon2 = ',epsilon2
        write(1,'(t5,a9,t34,e17.10)') 'Repul1 = ',repul1
        write(1,'(t5,a9,t34,e17.10)') 'Repul2 = ',repul2
        write(1,'(t5,a9,t34,e17.10)') 'Vprop1 = ',vprop1
        write(1,'(t5,a9,t34,e17.10)') 'Vprop2 = ',vprop2
        write(1,*)
        write(1,*) 'Calculated Force Coefficients (de-dimensionalized):'
        write(1,'(t5,a10,t34,e17.10)') 'Prescp1 = ',prescp1
        write(1,'(t5,a10,t34,e17.10)') 'Prescp2 = ',prescp2
        write(1,'(t5,a11,t34,e17.10)') 'Prescp12 = ',prescp12
        write(1,'(t5,a10,t34,e17.10)') 'Preljp1 = ',preljp1
        write(1,'(t5,a10,t34,e17.10)') 'Preljp2 = ',preljp2
        write(1,'(t5,a11,t34,e17.10)') 'Preljp12 = ',preljp12
        write(1,'(t5,a10,t34,e17.10)') 'Prerep1 = ',prerep1
        write(1,'(t5,a10,t34,e17.10)') 'Prerep2 = ',prerep2
        write(1,'(t5,a11,t34,e17.10)') 'Prerep12 = ',prerep12
        write(1,'(t5,a9,t34,e17.10)') 'Pretrad = ',pretrad
        write(1,'(t5,a9,t34,e17.10)') 'Prerotd = ',prerotd
        
        close(1)
        
        return
        end
