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
real*8 g_coup,epsilon,omega
real*8 band_width,gama_coup
real*8 V_exothermicity,omg_B,gamma_B,temperature
real*8,allocatable :: Vc(:),e_metal(:)
real*8 beta,gamma_D,lambda_B,V_reorg,V_barrier,Ei,ipos
real*8 omg_c,omg_scaled
real*8 s01,s02,x_cr
real*8,allocatable :: mass(:),omg(:)

!! Input/Output
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:),pop4(:,:)
complex*16,allocatable :: rho(:,:,:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: pos(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:),acc_old(:)
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop,flag_ortho
!complex*16,allocatable :: delr(:,:,:),delp(:,:,:),delacc(:,:,:)
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
real*8,allocatable ::hop_prob(:,:),W_overlap(:,:),hop_prob_net(:),knots(:),weights(:)
integer ielec_hop

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,dtq,total_time,curr_time,traj_num,tim_eq,E_corr
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg,old_energy
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun,current,above_barrier,barrier

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate,mol_LUMO
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2

integer,allocatable:: seed(:)

!! spin boson parameters

real*8 exo,gh
integer :: LZi,LZt
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
  read(10,*) lambda_B
  read(10,*) exo
  read(10,*) barrier
  read(10,*) mol_LUMO
  read(10,*) E_corr
  read(10,*) ipos
  read(10,*) Ei
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
  omega=0.0002
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
  allocate(pos(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass),acc_old(nclass))
  allocate(state(n_el),state_tentative(n_el),state_old(n_el),binary(nquant))
  allocate(mass(nclass),omg(nclass))
  allocate(V_k(nquant),V_k_old(nquant))
  allocate(Vc(nquant),e_metal(nquant))
  allocate(knots(nquant/2-mol_LUMO),weights(nquant/2-mol_LUMO))
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
  mass(1)=2000

if (mol_LUMO==2) then
  open(30,file="raw_x.txt")
  do i=1,int(nquant/2)-mol_LUMO+1
    read(30,*) knots(i)
  enddo
  close(30)

  open(30,file="raw_w.txt")
  do i=1,nquant/2-mol_LUMO+1
    read(30,*) weights(i)
  enddo
  close(30)
end if

gh=sqrt(2*lambda_B/(mass(1)*omega**2))

end subroutine setup
!---------------------------------------------------------- 
subroutine gaussian_random_number(rnd0)
!USE IFPORT
   !! generates gaussian distribution with center 0, sigma 1
   !! q0+sig*rnd gives center=q0, sigma=sig
   implicit none
   integer(kind=4) :: n,j,M,O,k
   real*8,intent(out)::rnd0
   real*8 rnd1,rnd2,pi
   pi=dacos(-1.d0)
   call random_number(rnd1)
   call random_number(rnd2)
   rnd0 = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!-----------------------------------------------------------
subroutine setup_initial_values2
integer :: n,m,p,y,i,j
real*8 :: rnd2,rnd1,rnd3,su,rnd4
real*8 :: sop, Ei,momentum(nclass)
integer :: instate(nquant)


call gaussian_random_number(rnd2)
call gaussian_random_number(rnd1)
v(1)=sqrt(2*Ei/mass(1))+sqrt((temperature/mass(1)))*rnd2  


if (mol_LUMO==2) ipos=gh
pos(1)=ipos+rnd1/(sqrt(mass(1)*omega**2/temperature))


end subroutine
!--------------------------------------------------------------------
subroutine test_setup_initial_values2(bins)
integer i,k,bins,ind,n,l,lnd
real*8 quant,phy
integer,allocatable :: bin(:),tin(:)



allocate(bin(bins))
allocate(tin(bins))

n=0
bin=0
tin=0
do i=1,10000
    call setup_initial_values2
    quant=0.d0
    kloop :do k=-int(bins/2),int(bins/2)
        quant=ipos+float(k)/(bins/2)*(1/(sqrt(mass(1)*omega**2/temperature)))*6
        ind=k+int(bins/2)
        if (quant>pos(1)) then
            bin(ind)=bin(ind)+1
            n=n+1
            exit kloop
        end if
    end do kloop

    lloop : do l=-int(bins/2),int(bins/2)
        phy=sqrt(2*Ei/mass(1))+float(l)/(bins/2)*sqrt(temperature/mass(1))*6
        lnd=l+int(bins/2)
        if(phy>v(1)) then
            tin(lnd)=tin(lnd)+1
            exit lloop
        end if
    end do lloop



end do



do i=1,bins
    write(166,*)ipos+float(i-bins/2)/(bins/2)*(1/(sqrt(mass(1)*omega**2/temperature))),bin(i)
    write(168,*)sqrt(2*Ei/mass(1))+float(i-bins/2)/(bins/2)*sqrt(temperature/mass(1)),tin(i)
enddo

end subroutine
!.........................................................
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



do i=1,100000
    call random_number(rnd1)
    if (mod(mol_LUMO,2)==0) then
        irnd1=ceiling(rnd1*(n_el-1))
    else
        irnd1=ceiling(rnd1*n_el)
    end if

    call random_number(rnd2)
    irnd2=ceiling(rnd2*(n_el-1))

    pr=dexp(-(V_k(unmetal(irnd2))-V_k(metal(irnd1)))/temperature)
    call random_number(rnd)
    if (rnd<pr) then
        swap1=metal(irnd1)
        swap2=unmetal(irnd2)
        metal(irnd1)=swap2
        unmetal(irnd2)=swap1
    end if

enddo

if (mol_LUMO==2) state(1)=2

state(mol_LUMO:n_el)=metal
binary=binary+binstate(state)

end subroutine 
!...........................................................
function fdist(x)
real*8 x,fdist

fdist=1/(1+exp(x/temperature))


end function
!-----------------------------------------------------------------------------
subroutine two_state_marcus
real*8 k1,k2
real*8 V_exo1,V_exo2,pl
integer tim

V_exo1=exo-V_k(2)
V_exo2=V_k(2)-exo
k1=2*pi*Vc(2)**2/hbar * 1.d0/sqrt(4*pi*lambda_B*temperature) * exp(-(lambda_B+V_exo1)**2/(4*lambda_B*temperature))
k2=2*pi*Vc(2)**2/hbar * 1.d0/sqrt(4*pi*lambda_B*temperature) * exp(-(lambda_B+V_exo2)**2/(4*lambda_B*temperature))

open(167,file="2state_marcus")
do tim=1,total_time,500
    pl=k1*dexp(-(k1+k2)*tim)/(k1+k2)+k2/(k1+k2)
    write(167,*) tim,pl
enddo
close(167)

end subroutine
!..............................................................................
subroutine main
  implicit none
  integer i,j,k,n,pstate(nquant)
  real*8 t1,t2
    
  call files(0)

  call cpu_time(t1)

  call initialize_averages

  !call test_setup_initial_values2(300)
  above_barrier=0.d0

  LZi=0
  LZt=0
  binary=0.d0
  do i=1,N_traj
    write(145,*) i,'traj_no'
    traj_num=i
    call init_cond
    acc(1)=-st_grad(state)/mass(1)  

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
  do i=mol_LUMO+1,nquant
       write(623,*)V_k(i),binary(i)/real(N_traj),fdist(V_k(i))
  enddo

  call cpu_time(t2);tim_tot=tim_tot+t2-t1
  call files(1)

  write(566,*) tim_tot
  open(256, file="ended")
    write(256,*) 'ended'
  close(256)

  if (nquant==2) call two_state_marcus
  write(1400,*) LZi,LZt
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

  call setup_initial_values2

  
  if (mol_LUMO==2) then
      do i=1,n_el
        state(i)=i+1
      end do
  else
    do i=1,n_el
        state(i)=i+1
    enddo
  end if


      
  curr_time=0.d0
  call setup_parameters(pos,state)

  if (nquant>2) call fermi4
  

end subroutine init_cond
!-----------------------------------------------------------------
subroutine velocity_verlet
real*8 :: delr(nclass),delv(nclass)
real*8 :: gama_dt,gamma_B,c0,c1,c2


x_old=pos
V_k_old=V_k

gamma_B=2*omega
gama_dt=gamma_B*dtc
c0=dexp(-gama_dt)
c1=1.d0/gama_dt*(1.d0-c0)
c2=1.d0/gama_dt*(1.d0-c1)
call stochastic_force(delr,delv)
pos=pos+c1*dtc*v+c2*dtc*dtc*acc+delr
acc_old=acc


call setup_parameters(pos,state)
acc(1)=-st_grad(state)/mass(1)
v=c0*v+(c1-c2)*dtc*acc_old+c2*dtc*acc+delv




end subroutine
!................................................................................

subroutine velocity_verlet2


x_old=pos
V_k_old=V_k
pos=pos+v*dtc+0.5*acc*dtc*dtc
acc_old=acc
call setup_parameters(pos,state)
acc(1)=-st_grad(state)/mass(1)
v=v+0.5*(acc+acc_old)*dtc
end subroutine

!...................................................................................

subroutine stochastic_force(delr,delv)
real*8, intent(out) :: delr(nclass),delv(nclass)
integer :: i
real*8 :: rnd1,rnd2,sig_r,sig_v,sig_rv,gdt,gamma_B
gamma_B=2*omega
gdt=gamma_B*dtc

do i=1,nclass

sig_r=dtc*dsqrt(temperature/mass(i)*1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
sig_v=dsqrt(temperature/mass(i)*(1-dexp(-2*gdt)))
sig_rv=(dtc*temperature/mass(i)*1.d0/gdt*(1-dexp(-gdt))**2)/(sig_r*sig_v) !correlation coefficient

call gaussian_random_number(rnd1)
call gaussian_random_number(rnd2)
delr(i)=sig_r*rnd1
delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
enddo

end subroutine stochastic_force
!........................................................................................



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
    if (curr_time>0.d0) call hop_state
    call velocity_verlet
    write(122,*)curr_time,TSE(state,pos),0.5*sum(mass*v*v)
    write(141,*)curr_time,TSE(state,pos)+0.5*sum(mass*v*v)
    
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
    pop(:,j)=pop(:,j)+binstate(state)
    
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

   write(1002,'(21f15.7)')(i-1)*nstep_avg*dtc,sum(pop(int(nquant/2+1):nquant,i))

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

pchem=0.0
do i=mol_LUMO+1,nquant
    if (V_k(i)>sH) then
        pchem=pchem+pops(i)!*gfac(V_k(i),sH)
        !write(30,*) V_k(i),gfac(V_k(i),sH)        
    end if
enddo
write(1335,*) pchem,above_barrier
!pchem=above_barrier-pchem
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
!.................................................................
subroutine setup_parameters(x_para,istate)
  implicit none
  integer i,sp,d
  real*8 si_diab(nbasis,2),Vb
  real*8 c_0,c_e,pi
  real*8 omg_max,delw
  real*8 rho,V_fil,V_empty
  CHARACTER :: fileplace
  real*8, intent(in) :: x_para(nclass)
  integer, intent(in) :: istate(n_el)

V_exothermicity = exo
pi=acos(-1.d0)
write(36,*) 1/(1+exp(-V_exothermicity/temperature))

if (mol_LUMO==2) sp=1
if (mol_LUMO==1) sp=2
if (sp.eq.1) then


    
  
  do i=1,nquant/2-mol_LUMO+1
    V_k(nquant/2-i+mol_LUMO)=-band_width/2.d0*(1/2.d0+knots(i)/2.d0)
    V_k(nquant/2+i+mol_LUMO-1)=band_width/2.d0*(1/2.d0+knots(i)/2.d0)
    Vc(nquant/2-i+2)=sqrt(gama_coup/(2*pi))*sqrt(band_width*weights(i))/2.d0
    Vc(nquant/2+i-1+2)=sqrt(gama_coup/(2*pi))*sqrt(band_width*weights(i))/2.d0
  enddo


V_empty=0.5*omega**2*sum(mass*x_para*x_para)
V_fil=0.5*omega**2*sum(mass*(x_para-gh)*(x_para-gh))+V_exothermicity

do d=1,mol_LUMO  
    do i=1,n_el
      if (state(i)==d) then    
        V_k(d)=V_fil
       else
        V_k(d)=V_empty
      end if
    enddo
enddo




else
  rho=(nquant-1)/(band_width)
  Vc=dsqrt(gama_coup/(2*pi*rho))
  V_empty=0.5*omega**2*sum(mass*x_para*x_para)
  V_fil=0.5*omega**2*sum(mass*(x_para-gh)*(x_para-gh))+V_exothermicity

  do i=mol_LUMO+1,nquant
    if (nquant>2) then
        e_metal(i)=-band_width/2.d0+(i-mol_LUMO-1)*band_width/real(nquant-mol_LUMO-1)
    else
        e_metal(i)=-band_width/2.d0+(i-mol_LUMO)*band_width/real(nquant-mol_LUMO)
    end if
  enddo

  do d=1,mol_LUMO  
    do i=1,n_el
      if (state(i)==d) then    
        V_k(d)=V_fil
       else
        V_k(d)=V_empty
      end if
    enddo
  enddo




  do i=mol_LUMO+1,nquant
    V_k(i)=e_metal(i)
  enddo

end if 

!write(455,*) Vc
!write(156,*)
!do i=1,nquant/2-mol_LUMO+1
!    write(156,*) nquant/2-i+mol_LUMO,V_k(nquant/2-i+mol_LUMO), 'up'
!    write(156,*) nquant/2+i+mol_LUMO-1, V_k(nquant/2+i+mol_LUMO-1),'down'
!enddo


!do i=1,nquant
!    write(159,*) i,V_k(i)
!enddo
!write(159,*)

!  do i=1,int(nquant/2)-mol_LUMO+1
!    write(33,*) i,weights(i)
!  enddo
!  close(33)

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

!stop

end subroutine setup_parameters
!-----------------------------------------------------------------  

function st_grad(istate)

integer :: istate(n_el),i,j,d
real*8 :: st_grad,grad_U1,grad_E_corr,E_max,grad_U0,grad_U2

grad_U0=omega**2*sum(mass*pos)
grad_U1=omega**2*sum(mass*(pos-gh))
grad_U2=omega**2*sum(mass*(pos-2*gh))


st_grad=0.d0

if (mol_LUMO==1) then 
    do d=1,mol_LUMO
        if (any(d==istate)) then
            st_grad=st_grad+grad_U1
        else
            st_grad=st_grad+grad_U0
        end if
    enddo
end if


if (mol_LUMO==2) then
    if (any(1==istate).and.any(2==istate)) then
        st_grad=st_grad+grad_U2
    else if (any(1==istate).or.any(2==istate)) then
        st_grad=st_grad+grad_U1
    else
        st_grad=st_grad+grad_U0
    end if    
end if



end function  
!---------------------------------------------------------- 
!..............................................................
subroutine check_energy_conservation(cstate,nstate,if_hop)
integer, intent(in) :: cstate(n_el),nstate(n_el)
integer, intent(out) :: if_hop
real*8 :: TE1,Eh,TE2


Eh=0.5*omega**2*sum(mass*pos*pos)

TE1=Eh+sum(V_k(cstate))+0.5*sum(mass*v*v)

TE2=Eh+sum(V_k(nstate))

if (TE1>TE2) then
    if_hop=1
    v(1)=sqrt(v(1)**2+2*(sum(V_k(cstate))-sum(V_k(nstate)))/mass(1)) ! why reverse velocity to get energy conservation
end if

end subroutine

!-----------------------------------------------------------
subroutine draw_pes(cstate,nstate)
  implicit none
  integer i
  integer, intent(in) :: cstate(n_el),nstate(n_el)
  real*8 x(nclass),dx,pot_U1,pot_U0



dx=60/1000.d0

do i=1,1000
   x(1)=-20+i*dx
   write(int(curr_time),*) x(1),TSE(cstate,x),TSE(nstate,x)
!   write(756,*)pos(1),ist_energy(gstate,pos(1))
!   write(757,*)pos(1),st_energy(estate),st_energy(nstate)
enddo

end subroutine draw_pes
!----------------------------------------------------------------- 
function TSE(istate,px)
real*8 TSE,px(nclass),E0,E1,E2
integer istate(n_el),d,i

TSE=0.d0


E0=0.5*omega**2*sum(mass*px*px)
E1=0.5*omega**2*sum(mass*(px-gh)*(px-gh))+exo
E2=0.5*omega**2*sum(mass*(px-2*gh)*(px-2*gh))+2*exo+E_corr

if (mol_LUMO==1) then
    do d=1,mol_LUMO  
        if (any(d==istate)) then
            TSE=TSE+E1
        else
            TSE=TSE+E0
        end if
    end do
end if

if (mol_LUMO==2) then
    if (any(1==istate).and.any(2==istate)) then
        TSE=TSE+E2
    else if (any(1==istate).or.any(2==istate)) then
        TSE=TSE+0.5*omega**2*mass(1)*(px(1)**2-2*gh*px(1)+2*gh**2)+exo
    else
        TSE=TSE+E0
    end if    
end if



do i=1,n_el
    if (istate(i)>mol_LUMO) TSE=TSE+V_k(istate(i))
enddo



end function
!-----------------------------------------------------------------
function crossing_prob(cstate,nstate,V_coup)
integer :: cstate(n_el),nstate(n_el)
real*8 :: Gam,den,crossing_prob
real*8 :: pi
real*8 :: V_coup

!pi=acos(-1.d0)
pi=3.1416

Gam=2*pi*V_coup**2

den=abs(v(1)*(st_grad(cstate)-st_grad(nstate)))

crossing_prob=1-dexp(-Gam/den)



end function

!-----------------------------------------------------------------
subroutine hop_state
  implicit none
  integer i,state_new(n_el),occupied(mol_LUMO),d,ihop,temp
  real*8 coup,exo_old,LZ,rnd,su,prob(nquant),exo_new,Gam,pi
  integer rep

  
  pi=acos(-1.d0)
call check_occupied(occupied)

call random_number(rnd)
ihop=0
LZ=0.d0
rep=0

dloop : do  d=1,mol_LUMO 
 if(occupied(d)==0) then !! Neutral
      outerloop: do i=1,n_el
            if (state(i)>mol_LUMO) then
                 state_new=state
                 temp=state_new(i)
                 state_new(i)=d
                 exo_new=TSE(state,pos)-TSE(state_new,pos)
                 old_energy=TSE(state,x_old)-TSE(state_new,x_old)

                if (exo_new*old_energy<0) then  
                    coup=Vc(temp)
                    Gam=abs(-mass(1)*omega**2*gh*v(1))
                    LZ=2*pi*coup**2/Gam
                    if (LZ<0.d0) LZ=0.d0
                    if(rnd<LZ) then
                         state=state_new
                         if (nquant==2) then
                            LZt=LZt+1
                            if (LZ>0.2) LZi=LZi+1
                         end if

                         exit dloop
                    endif
                end if
            end if
    enddo outerloop
 else    !! Anionic
    outerloop2: do i=mol_LUMO+1,nquant
        if(.not.(any(i==state))) then
            state_new=state
            temp=state_new(occupied(d))
            state_new(occupied(d))=i
            exo_new=TSE(state,pos)-TSE(state_new,pos)
            old_energy=TSE(state,x_old)-TSE(state_new,x_old)

            if (exo_new*old_energy<0) then
                coup=Vc(i)
                Gam=abs(-mass(1)*omega**2*gh*v(1))
                LZ=2*pi*coup**2/Gam
                if (LZ<0.d0) LZ=0.d0
                if(rnd<LZ) then
                    state=state_new
                    if (nquant==2) then
                        LZt=LZt+1
                        if (LZ>0.2) LZi=LZi+1
                    end if
                    exit dloop
                endif
            end if
        endif
    enddo outerloop2
 endif
enddo dloop

end subroutine hop_state
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
!-----------------------------------------------------------------  

End Module mod_afssh
