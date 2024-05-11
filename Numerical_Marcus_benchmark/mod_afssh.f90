Module mod_afssh
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.d0
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Potential
integer nquant,n_el
real*8 g_coup,epsilon
real*8 band_width,gama_coup,omega
real*8 V_exothermicity,omg_B,gamma_B,temperature
real*8,allocatable :: Vc(:),e_metal(:)
real*8 beta,gamma_D,lambda_B,V_reorg,V_barrier
real*8 omg_c,omg_scaled
real*8 s01,s02,x_cr
real*8,allocatable :: mass(:),omg(:),knots(:),weights(:)

!! Input/Output
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:),pop4(:,:)
complex*16,allocatable :: rho(:,:,:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: x(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:),acc_old(:)
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop,flag_ortho
complex*16,allocatable :: delr(:,:,:),delp(:,:,:),delacc(:,:,:)
complex*16,allocatable :: delr_old(:,:,:),delp_old(:,:,:)

!! Quantum
integer nbasis
integer,allocatable :: state(:),state_tentative(:),state_old(:),binary(:)
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_site(:,:),Hamil_diab(:,:),delH_dels(:,:,:),delH_dels_ad(:,:,:)
real*8,allocatable :: pot(:,:),force(:,:,:),force_old(:,:,:),delF(:,:,:)
complex*16,allocatable :: ci(:,:),ci_old(:,:),sigma(:,:,:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:,:),W_overlap(:,:),hop_prob_net(:)
integer ielec_hop

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,dtq,total_time,curr_time,traj_num,tim_eq,E_corr
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun,current,above_barrier,barrier
integer mol_LUMO
!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2

integer,allocatable:: seed(:)

!! spin boson parameters

real*8 exo,gh

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)
  real*8 rnd,c_0,c_e,kt

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) dtq
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) flag_ortho
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) n_el
  read(10,*) temperature
  read(10,*) band_width
  read(10,*) gama_coup
  read(10,*) iforward
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) seed2
  read(10,*) gh
  read(10,*) exo
  read(10,*) barrier
  read(10,*) mol_LUMO
  read(10,*) E_corr
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  !nquant=nquant*nb_vib
  nbasis=nquant

  energy_cutoff=energy_cutoff*wave_to_J
  !temperature=temperature*au2J/kb
  !band_width=band_width*au2J
  !gama_coup=gama_coup*au2J
  !kt=kb*temperature
  !beta=1.d0/(kb*temperature)
  nsteps=nint(total_time/dtc)+1
  !lambda_D=V_reorg/4.d0

  !-----------------------------------------------------------------  
  i=nsteps/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop4(4,i))
  allocate(rho(nquant,nquant,i))
  allocate(x(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass),acc_old(nclass))
  allocate(state(n_el),state_tentative(n_el),state_old(n_el),binary(nquant))
  allocate(mass(nclass),omg(nclass))
  allocate(delr(nquant,nquant,nclass),delp(nquant,nquant,nclass),delacc(nquant,nquant,nclass))
  allocate(delr_old(nquant,nquant,nclass),delp_old(nquant,nquant,nclass))
  allocate(si_adiab(nbasis,nquant),ci(nquant,n_el),V_k(nquant),V_k_old(nquant),sigma(nquant,nquant,n_el))
  allocate(Hamil_site(nbasis,nbasis),Hamil_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nclass),delH_dels_ad(nquant,nquant,nclass))
  allocate(pot(nquant,nquant),force(nquant,nquant,nclass),force_old(nquant,nquant,nclass),delf(nquant,nquant,nclass))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass),vdotd(nquant,nquant),hop_prob(nquant,n_el),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant,n_el),si_adiab_prev(nbasis,nquant))
  allocate(Vc(nquant),e_metal(nquant))
  allocate(knots(n_el-1),weights(n_el-1))

  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0
  !  barrier=1.d-3

  if (mol_LUMO==2) then
    open(30,file="raw_x.txt")
    do i=1,n_el-1
        read(30,*) knots(i)
    enddo
    close(30)

    open(30,file="raw_w.txt")
    do i=1,n_el-1
        read(30,*) weights(i)
    enddo
    close(30)
  end if

end subroutine setup
!---------------------------------------------------------- 

subroutine simple_fermi
integer :: unfilled(n_el)
integer :: n,irnd1,i,irnd2,k,l,tot
real*8 rnd7,rnd8,pr,rnd9,tot_e
integer :: forb,vorb,bin(n_el)



n=1
do i=2,nquant
    if (.not.(any(i==state))) then
        unfilled(n)=i
        n=n+1
    end if
enddo

!write(78,*) state
!write(78,*) unfilled

!do i=1,nquant
!    write(129,*) i,V_k(i)
!enddo

!write(54,*) n_el

tot=1
do k=1,tot
    bin=0
    do l=1,500000
        call random_number(rnd7)
        irnd1=ceiling(n_el*rnd7)
        forb=state(irnd1)
   !     if (irnd1==1) cycle     ! for the case of filled orbital
        call random_number(rnd8)
        irnd2=ceiling(n_el*rnd8)
   !      bin(irnd2)=bin(irnd2)+1

        vorb=unfilled(irnd2)
        
  !      if (unfilled(irnd2)==1) cycle  !for the case of filled orbital

        pr=exp(-(V_k(vorb)-V_k(forb))/temperature)
        call random_number(rnd9)
        if (pr>rnd9) then
     !       write(125,*) 'v',vorb,V_k(vorb),'f',forb,V_k(forb)
            state(irnd1)=vorb
            unfilled(irnd2)=forb
        end if
    enddo

    binary=binary+binstate(state)
  

    !do i=1,n_el
    !    write(217,*) i,bin(i)
    !enddo
enddo
write(133,*) binary

end subroutine
!..............................................................................

subroutine fermi2

integer :: bstate(nquant),binary(nquant)
integer :: i,j,k,l,irnd1,irnd2
real*8 :: pr,rndf,rnd1,rnd2,rnd9
integer :: vorb,forb

bstate=binstate(state)
do k=1,1
  do l=1,200000
    call random_number(rnd1)
    call random_number(rnd2)

    irnd1=ceiling(rnd1*nquant)
    irnd2=ceiling(rnd2*nquant)
        
    if ((irnd1>mol_LUMO).and.(irnd2>mol_LUMO)) then
        if ((bstate(irnd1)==1).and.(bstate(irnd2)==0)) then 
            forb=irnd1
            vorb=irnd2
        else
            cycle
        end if
    else
        cycle
    end if
    
    pr=dexp(-(V_k(vorb)-V_k(forb))/temperature)
    call random_number(rnd9)
    if (pr>rnd9) then
        bstate(forb)=0
        bstate(vorb)=1
    end if
  enddo
  binary=binary+bstate
  write(159,*) bstate
enddo

state=numstate(bstate,1)

write(155,*) binary
!do i=3,nquant
!    write(623,*) V_k(i),binary(i)/real(k),fdist(V_k(i))
!enddo

end subroutine 

!.............................................................................

!subroutine fermi4
!integer i, metal(n_el-1),irnd1,irnd2,j,unmetal(n_el-1)
!real*8 rnd1,rnd2,pr,rnd
!integer swap1,swap2

!metal=state(mol_LUMO:n_el+1)
!do i=1,n_el-1
!    unmetal(i)=n_el+1+i
!enddo


!write(411,*)metal
!write(411,*)unmetal


!do i=1,500000
!    call random_number(rnd1)
!    irnd1=ceiling(rnd1*(n_el-1))

!    call random_number(rnd2)
!    irnd2=ceiling(rnd2*(n_el-1))

!    pr=dexp(-(V_k(unmetal(irnd2))-V_k(metal(irnd1)))/temperature)

   ! write(411,*) metal
   ! write(411,*) unmetal, "before swap"
!    call random_number(rnd)
!    if (rnd<pr) then
!        swap1=metal(irnd1)
!        swap2=unmetal(irnd2)
!        metal(irnd1)=swap2
!        unmetal(irnd2)=swap1
!    end if
    !write(411,*) metal 
   !write(411,*) unmetal, "after swap"

!enddo

!state(1)=2
!state(mol_LUMO:n_el)=metal
!write(411,*) state
!binary=binary+binstate(state)
!state=0

!write(411,*) state
!write(411,*)metal
!write(411,*)unmetal

!write(411,*) state
!end subroutine 
!-----------------------------------------------------------------------------

subroutine fermi4
integer i,irnd1,irnd2,j
real*8 rnd1,rnd2,pr,rnd
integer swap1,swap2
integer, allocatable :: metal(:),unmetal(:)

if (mod(mol_LUMO,2)==0) then
    allocate(metal(n_el-1))
    allocate(unmetal(n_el-1))
else
    allocate(metal(n_el))
    allocate(unmetal(n_el-1))
end if


metal=state(mol_LUMO:n_el)
do i=1,n_el-1
    unmetal(i)=n_el+1+i
enddo


!write(411,*)metal
!write(411,*)unmetal
!write(411,*)

do i=1,10000
    call random_number(rnd1)
    if (mod(mol_LUMO,2)==0) then
        irnd1=ceiling(rnd1*(n_el-1))
    else
        irnd1=ceiling(rnd1*n_el)
    end if

    call random_number(rnd2)
    irnd2=ceiling(rnd2*(n_el-1))

    pr=dexp(-(V_k(unmetal(irnd2))-V_k(metal(irnd1)))/temperature)
 !   write(155,*) curr_time,pr
   ! write(411,*) metal
   ! write(411,*) unmetal, "before swap"
    call random_number(rnd)
    if (rnd<pr) then
        swap1=metal(irnd1)
        swap2=unmetal(irnd2)
        metal(irnd1)=swap2
        unmetal(irnd2)=swap1
    end if
    !write(411,*) metal 

enddo


if (mol_LUMO==2) state(1)=2

!!write(166,*) state(mol_LUMO:n_el)
!write(166,*) metal
!write(166,*)
state(mol_LUMO:n_el)=metal
!write(411,*) state
binary=binary+binstate(state)


end subroutine


!..............................................................................    

function fdist(x)
real*8 x,fdist

fdist=1/(1+exp(x/temperature))


end function


!..............................................................................
subroutine main
  implicit none
  integer i,j,k,n,pstate(nquant)
  real*8 t1,t2
    
  call files(0)

  call cpu_time(t1)

  call setup_parameters
  call initialize_averages
  
  above_barrier=0.d0

  binary=0.d0
  do i=1,N_traj
    traj_num=i
    call init_cond    
    pstate=binstate(state)

    call evolve(nsteps)

    do j=mol_LUMO+1,nquant
     if (V_k(j)>barrier) then
        above_barrier=above_barrier+pstate(j)
    end if
    enddo
  enddo
  above_barrier=above_barrier/N_traj
  call write_average
  write(156,*) above_barrier,n_el,above_barrier/n_el 
  do i=3,nquant
       write(623,*)V_k(i),binary(i)/real(N_traj),fdist(V_k(i))
  enddo
  call cpu_time(t2);tim_tot=tim_tot+t2-t1
  call files(1)

  open(256, file="ended")
    write(256,*) 'ended'
  close(256)

end subroutine main
!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(100,file="pop.out")
    open(1001,file="current.out")
  else
    close(100)
    close(1001)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  pop=0.d0
  pop4=0.d0
  current=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i,j
  real*8 sig_x,sig_p,rnd,ak,su
  real*8 energy_0,pstate(nquant)

  if (mod(mol_LUMO,2)==0) then
    do i=1,n_el
        state(i)=i+mol_LUMO-1
    enddo
  else
      do i=1,n_el
        state(i)=i+1
      enddo
  end if

 ! state(1)=1
 ! state(n_el)=state(n_el)+5
 ! state(n_el-1)=state(n_el-1)+3
 ! state(n_el-2)=state(n_el-2)+10  
  curr_time=0.d0
!  call simple_fermi
  call fermi4
  write(123,*) state
  
!  pstate=binstate(state)

!  above_barrier=0.d0
!  do i=1,nquant
!    if (V_k(i)>barrier) then
!        above_barrier=above_barrier+pstate(i)
!    end if
!  enddo
  

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,j,nstep_sm,iflag_coll,i_do_something
  real*8 t1,t2
  integer iterm

  !call cpu_time(t1)

  iterm=0
  do i=1,nsteps
    call average(i)
  !  call hop_state
    call two_state_hop
    curr_time=curr_time+dtc
  enddo

  !call cpu_time(t2)
  !tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!----------------------------------------------------------------- 

function binstate(nstate)
integer :: nstate(n_el),binstate(nquant)
integer :: i


binstate=0
do i=1,n_el
    binstate(nstate(i))=1
enddo

end function

!....................................................................

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  integer occupied,j,binf(nquant)

  if((mod(i,nstep_avg)==1.or.nstep_avg==1)) then
    if(nstep_avg==1) then
      j=i
    else
      j=i/nstep_avg+1
    endif

    binf=binstate(state)
    pop(:,j)=pop(:,j)+binf

    if  (binf(1)==0) then
        if (binf(2)==0) then
            pop4(1,j)=pop4(1,j)+1
        else
            pop4(2,j)=pop4(2,j)+1
        end if
    else
        if (binf(2)==0) then
           pop4(3,j)=pop4(3,j)+1
        else
           pop4(4,j)=pop4(4,j)+1
        end if
    end if
  
    !write(188,*) curr_time
!    write(188,*) pop4(:,j)
    !write(188,*) binf(1),binf(2)
    !write(188,*)
  endif

end subroutine average
!-----------------------------------------------------------------  

subroutine write_average
  !! Writes the final useful output
  implicit none
  integer i,j,i1,k
  real*8 nf

  nf=dfloat(n_traj)

  pop=pop/nf
  pop4=pop4/nf

  do i=1,nsteps/nstep_avg
   write(100,'(21f15.7)')(i-1)*nstep_avg*dtc,pop(1,i),pop(2,i),pop4(1,i),pop4(2,i),pop4(3,i),pop4(4,i)
   call calchemcurr(pop(:,i),barrier,current)
   write(1001,'(21f15.7)')(i-1)*nstep_avg*dtc,current
  enddo

  

end subroutine write_average
!-----------------------------------------------------------------  

subroutine calchemcurr(pops,sH,pchem)
real*8, intent(in) :: pops(nquant)
real*8, intent(in) :: sH
real*8, intent(out) :: pchem
integer :: i,k 
real*8 :: nex,ex_el

!if (curr_time<10) then
!    do i=2,nquant
!        if (V_k(i)>sH) then
!            ex_el=ex_el+pops(i)
!        end if
!    enddo
!end if
write(65,*) curr_time
pchem=0.0
do i=mol_LUMO+1,nquant
    if (V_k(i)>sH) then
        pchem=pchem+pops(i)!*gfac(V_k(i),sH)
        !write(30,*) V_k(i),gfac(V_k(i),sH)        
    end if
enddo
write(761,*) above_barrier,pchem

!pchem=(above_barrier-pchem)
pchem=pchem-above_barrier
end subroutine 

!.................................................................

function gfac(epi,eps)
real*8 :: epi,eps, gfac
real*8 :: theta, thetac, Dx, lamb,dth
integer i

thetac=acos((eps/epi)**(0.5))
gfac=0.d0
theta=0.d0
dth=thetac/1000
Dx=44.32
lamb=52.91


do i=0,1000
    theta=i*dth
    gfac=gfac+dth*sin(theta)*exp(-Dx/(lamb*cos(theta)))
enddo


end function
!...............................................................

    




!.................................................................
subroutine setup_parameters
  implicit none
  integer i,sp
  real*8 si_diab(nbasis,2),Vb
  real*8 c_0,c_e
  real*8 omg_max,delw
  real*8 rho
  real*8 para_data(4)
  CHARACTER :: fileplace

lambda_B=0.5*2000*0.0002**2*gh**2
V_exothermicity = exo

write(36,*) 1/(1+exp(-V_exothermicity/temperature))

write(56,*) gh,V_exothermicity,gama_coup,temperature
sp=1
if (sp.eq.1) then

   rho=(nquant-1)/(band_width)
   Vc=dsqrt(gama_coup/(2*pi*rho))

    
  
  V_k(1)=exo
  V_k(2)=exo
  
  do i=1,nquant/2-mol_LUMO+1
    V_k(nquant/2-i+mol_LUMO)=-band_width/2.d0*(1/2.d0+knots(i)/2.d0)
    V_k(nquant/2+i+mol_LUMO-1)=band_width/2.d0*(1/2.d0+knots(i)/2.d0)
    Vc(nquant/2-i+2)=sqrt(gama_coup/(2*pi))*sqrt(band_width*weights(i))/2.d0
    Vc(nquant/2+i-1+2)=sqrt(gama_coup/(2*pi))*sqrt(band_width*weights(i))/2.d0
  enddo



else
  rho=(nquant-1)/(band_width)
  Vc=dsqrt(gama_coup/(2*pi*rho))
  do i=mol_LUMO+1,nquant
    e_metal(i)=-band_width/2.d0+(i-mol_LUMO-1)*band_width/real(nquant-mol_LUMO-1)
  enddo


  V_k(1)=V_exothermicity
  V_k(2)=V_exothermicity
  do i=mol_LUMO+1,nquant
    V_k(i)=e_metal(i)
  enddo

end if

!do i=1,nquant
!    write(677,*) i,V_k(i)
!enddo
!write(667,*)
!-----------------------------------------------------------------------------
!V_reorg =0.00027*au2J
!write(6,*) (V_reorg+V_exothermicity)**2/(4*V_reorg)/au2J
!write(6,*) 2*2*pi*gama_coup*kb*temperature/(hbar*omg*sqrt(2*V_reorg*kb*temperature))

!write(6,*) e_metal(nquant/2-1:nquant/2+2)/wave_to_J
!write(6,*) Vc(nquant/2-1:nquant/2+2)/wave_to_J

!write(6,*) e_metal(nquant/2-1:nquant/2+2)/au2J
!write(6,*) Vc(nquant/2-1:nquant/2+2)/au2J

!write(6,*) temperature
!write(6,*) omg_B
!write(6,*) gamma_B

end subroutine setup_parameters
!-----------------------------------------------------------------  

function numstate(bstate,jflag)
integer :: numstate(n_el),bstate(nquant),jflag,i,k

k=1
do i=1,nquant
    if (bstate(i)==jflag) then
        numstate(k)=i
        k=k+1
    end if
enddo 


end function    
!.........................................................................
!---------------------------------------------------------- 

subroutine draw_pes
  implicit none
  integer i

  call init_cond

  do i=1,n_el
    state(i)=i
  enddo
  !state(n_el-1)=state(n_el-1)+2
!  state(n_el)=state(n_el)+2

  do i=1,1000
    x(1)=-5d-10+20.d-10*i/999.d0
    write(20,*) x(1)*1.d10,V_k(1:6)/au2J
    write(21,*) x(1)*1.d10,pot_en/au2J
  enddo
  stop

end subroutine draw_pes
!-----------------------------------------------------------------  

function TSE(istate)
real*8 TSE,E0,E1
integer istate(n_el),d,i

TSE=0.d0


TSE=sum(V_k(istate))


if ((any(1==istate).and.(any(2==istate)))) then
    TSE=TSE+E_corr
!    write(143,*) curr_time,E_corr
end if


end function



!.................................................................
subroutine hop_state
  implicit none
  integer i,state_new(n_el),occupied(mol_LUMO),d,d_comp
  real*8 coup,V_exo,k_Marc,rnd,su,prob(nquant)

!write(1300,*)
!write(1300,*) state
!write(1300,*) curr_time
!write(1300,*)'a','   state(i) ', '   occupied(d) ','   occupied(d_comp) ','   d ','    d_comp '
call check_occupied(occupied)
!write(129,*) occupied
dloop :do d=1,2
    if (d==1) d_comp=2 
    if (d==2) d_comp=1
  
  if(occupied(d)==0) then !! (0,1), (1,0) or (0,0)
    su=0.d0
    call random_number(rnd)
    outerloop: do i=1,n_el
        if (state(i)>mol_LUMO) then
          state_new=state
          state_new(i)=d
          if (occupied(d_comp)==0) then  
                V_exo=exo+lambda_B-V_k(state(i))   !(0,0) to (0,1) or (1,0)
!                write(1300,*)'a',state(i), occupied(d),occupied(d_comp),d,d_comp
          else
                V_exo=exo-lambda_B-V_k(state(i))+E_corr ! (0,1) or (1,0) to (1,1)
!                write(1300,*) 'b',state(i),occupied(d),occupied(d_comp),d,d_comp
           end if

           coup=Vc(state(i))
           k_Marc=2*pi*coup**2/hbar * 1.d0/sqrt(4*pi*lambda_B*temperature) * exp(-(V_exo+lambda_B)**2/(4*lambda_B*temperature))
           su=su+k_Marc*dtc
           if(rnd<su) then
             state=state_new
             exit dloop
           endif
        end if
    enddo outerloop
  else    !! Anionic ! (1,0) or (0,1) or (1,1)
    su=0.d0
    call random_number(rnd)
    outerloop2: do i=1,nquant
        if (i>mol_LUMO) then 
            if(.not.(any(i==state))) then
                state_new=state
                state_new(occupied(d))=i
                coup=Vc(i)

                if (occupied(d_comp)==0) then  
                    V_exo=-exo-lambda_B+V_k(i) ! (1,0) or (0,1) to (0,0)
 !                   write(1300,*) 'c',i,occupied(d),occupied(d_comp),d,d_comp
                else
                    V_exo=-exo+lambda_B+V_k(i)-E_corr !(1,1) to (0,1) or (1,0)
 !                   write(1300,*) 'd',i,occupied(d),occupied(d_comp),d,d_comp
                end if

                k_Marc=2*pi*coup**2/hbar * 1.d0/sqrt(4*pi*lambda_B*temperature) * exp(-(V_exo+lambda_B)**2/(4*lambda_B*temperature))
                su=su+k_Marc*dtc
                if(rnd<su) then
                    state=state_new
                    exit dloop
                endif
            endif
        end if
    enddo outerloop2
  endif
end do dloop 




end subroutine hop_state


!................................................................
subroutine two_state_hop
integer :: AO(mol_LUMO),new_AO(mol_LUMO),i,j,swap
integer :: metal(nquant-mol_LUMO),new_metal(nquant-mol_LUMO)
integer :: binst(nquant),outbin(nquant)
real*8 :: su,V_curr,V_new,V_exo,k_Marc,rnd

binst=binstate(state)

AO=binst(1:mol_LUMO)
metal=binst(mol_LUMO+1:nquant)

su=0.d0
call random_number(rnd)
dloop : do i=1,mol_LUMO
    do j=1,nquant-mol_LUMO
        new_AO=AO
        new_metal=metal
        if (AO(i).ne.metal(j)) then
            swap=new_AO(i)
            new_AO(i)=new_metal(j)
            new_metal(j)=swap
            V_curr=AO_energy(AO)+sum(metal*V_k(mol_LUMO+1:nquant))
            V_new=AO_energy(new_AO)+sum(new_metal*V_k(mol_LUMO+1:nquant))
            V_exo=V_new-V_curr
            k_Marc=2*pi*Vc(j)**2/hbar * 1.d0/sqrt(4*pi*lambda_B*temperature) * exp(-(V_exo+lambda_B)**2/(4*lambda_B*temperature))
            su=su+k_Marc*dtc
            if (su>rnd) then
                outbin(1:mol_LUMO)=new_AO
                outbin(mol_LUMO+1:nquant)=new_metal
                state=numstate(outbin,1)
                exit dloop
            end if
        end if
    end do
enddo dloop
                


end subroutine

!----------------------------------------------------------------
function AO_energy(AOs)
real*8 :: AO_energy
integer AOs(mol_LUMO)


if (mol_LUMO==2) then
    if (AOs(1)==0) then
        if (AOs(2)==0) then
            AO_energy=0.d0
        else
            AO_energy=lambda_B+exo
        end if
    else
        if (AOs(2)==0) then
            AO_energy=lambda_B+exo
        else
            AO_energy=2*exo+E_corr
        end if
    end if
end if


if (mol_LUMO==1) then
    if (AOs(1)==0) then
        AO_energy=0.d0
    else
        AO_energy=exo
    end if
end if



end function
!.................................................................
subroutine check_occupied(occupied)
  implicit none
  integer,intent(out)::occupied(mol_LUMO)
  integer i,d

  occupied=0 !! Neutral state
do d=1,mol_LUMO
  do i=1,n_el
    if(state(i)==d) occupied(d)=i !! Anionic
  enddo
enddo


end subroutine check_occupied

!.................................................................
End Module mod_afssh
