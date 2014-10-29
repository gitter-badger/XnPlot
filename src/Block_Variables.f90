!> @brief This module contains the block allocatable variables
module Block_Variables
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision          ! Integers and reals precision definition.
USE Data_Type_PostProcess ! Definition of Type_PostProcess.
USE Data_Type_Vector      ! Definition of Type_Vector.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
public
save
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type(Type_Vector), allocatable:: node(:,:,:)      ! Nodes coordinates.
type(Type_Vector), allocatable:: NFi(:,:,:)       ! Face i normal versor.
type(Type_Vector), allocatable:: NFj(:,:,:)       ! Face j normal versor.
type(Type_Vector), allocatable:: NFk(:,:,:)       ! Face k normal versor.
type(Type_Vector), allocatable:: NFiS(:,:,:)      ! Face i normal versor with surface area module.
type(Type_Vector), allocatable:: NFjS(:,:,:)      ! Face j normal versor with surface area module.
type(Type_Vector), allocatable:: NFkS(:,:,:)      ! Face k normal versor with surface area module.
real(R8P),         allocatable:: Si(:,:,:)        ! Face i area.
real(R8P),         allocatable:: Sj(:,:,:)        ! Face j area.
real(R8P),         allocatable:: Sk(:,:,:)        ! Face k area.
real(R8P),         allocatable:: volume(:,:,:)    ! Volumes of cells.
integer(I4P),      allocatable:: icc(:,:,:)       ! Cell centered icc values.
integer(I4P),      allocatable:: ricc(:,:,:)      ! Cell centered rcc values.
integer(I4P),      allocatable:: vicc(:,:,:)      ! Node centered rcc values.
type(Type_Vector), allocatable:: momentum(:,:,:)  ! Momentum.
real(R8P),         allocatable:: pressure(:,:,:)  ! Pressure.
real(R8P),         allocatable:: f(:,:,:)         ! Level set function.
real(R8P),         allocatable:: f0(:,:,:)        ! Level 0 (level set).
real(R8P),         allocatable:: visc(:,:,:)      ! Viscosity.
real(R8P),         allocatable:: vitl(:,:,:)      ! Turbulent viscosity.
real(R8P),         allocatable:: ken(:,:,:)       ! Turbulent kinetic energy.
real(R8P),         allocatable:: eps(:,:,:)       ! Turbulent kinetic energy dissipation.
real(R8P),         allocatable:: vord(:,:,:)      ! Variable to identify vortices (lambda 2).
real(R8P),         allocatable:: qfac(:,:,:)      ! Variable to identify vortices (q factor).
real(R8P),         allocatable:: heli(:,:,:)      ! Helicity.
type(Type_Vector), allocatable:: vorticity(:,:,:) ! Vorticity.
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Procedure for allocating block variables.
  subroutine block_allocate(pp,gc,Ni,Nj,Nk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(IN):: pp       !< Post-processor data.
  integer(I4P),           intent(IN):: gc(1:6)  !< Number of ghost cells.
  integer(I4P),           intent(IN):: Ni,Nj,Nk !< Number of cells.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call block_deallocate()
  allocate(node(0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; node = 0._R8P
  if (pp%fcc) then
                      allocate(icc (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; icc  = 0_I4P
                      allocate(ricc(0      :Ni+1    ,0      :Nj+1    ,0      :Nk+1    )) ; ricc = 0_I4P
    if (.not.pp%cell) allocate(vicc(0      :Ni+1    ,0      :Nj+1    ,0      :Nk+1    )) ; vicc = 0_I4P
  endif
  if (pp%sol) then
    allocate(momentum(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; momentum = 0._R8P
    allocate(pressure(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; pressure = 0._R8P
    if (pp%level_set) then
      allocate(f (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; f  = 0._R8P
      allocate(f0(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; f0 = 0._R8P
    endif
    if (pp%zeroeq.or.pp%oneeq.or.pp%twoeq) then
      allocate(visc(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; visc = 0._R8P
    endif
    if (pp%oneeq) then
      allocate(vitl(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; vitl = 0._R8P
    endif
    if (pp%twoeq) then
      allocate(ken(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; ken = 0._R8P
      allocate(eps(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; eps = 0._R8P
    endif
    if (pp%vordet) then ! computing vord variable
      allocate(NFi      (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; NFi       = 0._R8P
      allocate(NFj      (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; NFj       = 0._R8P
      allocate(NFk      (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; NFk       = 0._R8P
      allocate(NFiS     (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; NFiS      = 0._R8P
      allocate(NFjS     (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; NFjS      = 0._R8P
      allocate(NFkS     (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; NFkS      = 0._R8P
      allocate(Si       (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; Si        = 0._R8P
      allocate(Sj       (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; Sj        = 0._R8P
      allocate(Sk       (0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6))) ; Sk        = 0._R8P
      allocate(volume   (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; volume    = 0._R8P
      allocate(vord     (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; vord      = 0._R8P
      allocate(qfac     (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; qfac      = 0._R8P
      allocate(heli     (1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; heli      = 0._R8P
      allocate(vorticity(1-gc(1):Ni+gc(2),1-gc(3):Nj+gc(4),1-gc(5):Nk+gc(6))) ; vorticity = 0._R8P
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine block_allocate

  !> @brief Procedure for deallocating block variables.
  subroutine block_deallocate()
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(node))      deallocate(node)
  if (allocated(icc))       deallocate(icc)
  if (allocated(ricc))      deallocate(ricc)
  if (allocated(vicc))      deallocate(vicc)
  if (allocated(momentum))  deallocate(momentum)
  if (allocated(pressure))  deallocate(pressure)
  if (allocated(f))         deallocate(f)
  if (allocated(f0))        deallocate(f0)
  if (allocated(visc))      deallocate(visc)
  if (allocated(vitl))      deallocate(vitl)
  if (allocated(ken))       deallocate(ken)
  if (allocated(eps))       deallocate(eps)
  if (allocated(NFi))       deallocate(NFi)
  if (allocated(NFj))       deallocate(NFj)
  if (allocated(NFk))       deallocate(NFk)
  if (allocated(NFiS))      deallocate(NFiS)
  if (allocated(NFjS))      deallocate(NFjS)
  if (allocated(NFkS))      deallocate(NFkS)
  if (allocated(Si))        deallocate(Si)
  if (allocated(Sj))        deallocate(Sj)
  if (allocated(Sk))        deallocate(Sk)
  if (allocated(volume))    deallocate(volume)
  if (allocated(vord))      deallocate(vord)
  if (allocated(qfac))      deallocate(qfac)
  if (allocated(heli))      deallocate(heli)
  if (allocated(vorticity)) deallocate(vorticity)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine block_deallocate

  !> @brief Procedure for interpolating cell centered variable into node centered one.
  subroutine varinterpolation(Ni,Nj,Nk,var,vari)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  Ni,Nj,Nk                   !< Block dimensions.
  real(R_P),    intent(IN)::  var (0:Ni+1,0:Nj+1,0:Nk+1) !< Cell centered variable.
  real(R_P),    intent(OUT):: vari(0:Ni  ,0:Nj  ,0:Nk  ) !< Node centered interpolated variable.
  integer(I_P)::              i,j,k                      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k)         &
  !$OMP SHARED(Ni,Nj,Nk,var,vari)
  !$OMP DO
  do k=0,Nk
    do j=0,Nj
      do i=0,Ni
        vari(i,j,k) = var(i+1,j+1,k+1)  + var(i,j+1,k+1) &
                    + var(i+1,j  ,k+1)  + var(i,j,  k+1) &
                    + var(i+1,j+1,k  )  + var(i,j+1,k  ) &
                    + var(i+1,j  ,k  )  + var(i,j  ,k  )
        vari(i,j,k) = 0.125_R_P*vari(i,j,k)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine varinterpolation

  !> @brief Procedure for computing block metrics.
  subroutine compute_metrics(gc,Ni,Nj,Nk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: gc(1:6)           !< Number of ghost cells.
  integer(I4P), intent(IN):: Ni,Nj,Nk          !< Number of cells.
  type(Type_Vector)::        NFS,s1,s2,nd,db   !< Dummy vector variables.
  real(R8P)::                signi,signj,signk !< Dummy variables for checking the directions of normals.
  real(R8P)::                Vx,Vy,Vz          !< Dummy variables for computing volume.
  real(R8P)::                xp,yp,zp          !< Dummy variables for computing face coordinates.
  real(R8P)::                xm,ym,zm          !< Dummy variables for computing face coordinates.
  integer(I4P)::             i,j,k             !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing faces normals
  ! positioning at the middle of the block
  i = max(1,Ni/2)
  j = max(1,Nj/2)
  k = max(1,Nk/2)
  ! checking the direction of i normals
  s1 = node(i,j  ,k) - node(i,j-1,k-1)
  s2 = node(i,j-1,k) - node(i,j,  k-1)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(node(i,  j,k)+node(i,  j-1,k)+node(i,  j,k-1)+node(i,  j-1,k-1))
  s2 = 0.25_R_P*(node(i-1,j,k)+node(i-1,j-1,k)+node(i-1,j,k-1)+node(i-1,j-1,k-1))
  db = s1 - s2
  signi = sign(1._R_P,(nd.dot.db))
  ! checking the direction of j normals
  s1 = node(i,j,k  ) - node(i-1,j,k-1)
  s2 = node(i,j,k-1) - node(i-1,j,k  )
  nd = s1.cross.s2
  s1 = 0.25_R_P*(node(i,j,  k)+node(i-1,j,  k)+node(i,j,  k-1)+node(i-1,j,  k-1))
  s2 = 0.25_R_P*(node(i,j-1,k)+node(i-1,j-1,k)+node(i,j-1,k-1)+node(i-1,j-1,k-1))
  db = s1 - s2
  signj = sign(1._R_P,(nd.dot.db))
  ! checking the direction of k normals
  s1 = node(i,  j,k) - node(i-1,j-1,k)
  s2 = node(i-1,j,k) - node(i,  j-1,k)
  nd = s1.cross.s2
  s1 = 0.25_R_P*(node(i,j,k  )+node(i-1,j,k  )+node(i,j-1,k  )+node(i-1,j-1,k  ))
  s2 = 0.25_R_P*(node(i,j,k-1)+node(i-1,j,k-1)+node(i,j-1,k-1)+node(i-1,j-1,k-1))
  db = s1 - s2
  signk = sign(1._R_P,(nd.dot.db))
  !$OMP PARALLEL DEFAULT(NONE)                        &
  !$OMP PRIVATE(i,j,k,NFS,Vx,Vy,Vz,xp,yp,zp,xm,ym,zm) &
  !$OMP SHARED(Ni,Nj,Nk,gc,signi,signj,signk,node,Si,Sj,Sk,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=0-gc(1),Ni+gc(2)
        NFS = face_normal4(pt1 = node(i,j-1,k-1), &
                           pt2 = node(i,j  ,k-1), &
                           pt3 = node(i,j  ,k  ), &
                           pt4 = node(i,j-1,k  ))
        NFS = NFS*signi
        NFiS(i,j,k) = NFS
        NFi (i,j,k) = normalize(NFS)
        Si  (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=0-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        NFS = face_normal4(pt1 = node(i-1,j,k-1), &
                           pt2 = node(i-1,j,k  ), &
                           pt3 = node(i  ,j,k  ), &
                           pt4 = node(i  ,j,k-1))
        NFS = NFS*signj
        NFjS(i,j,k) = NFS
        NFj (i,j,k) = normalize(NFS)
        Sj  (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  !$OMP DO
  do k=0-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        NFS = face_normal4(pt1 = node(i-1,j-1,k), &
                           pt2 = node(i  ,j-1,k), &
                           pt3 = node(i  ,j  ,k), &
                           pt4 = node(i-1,j  ,k))
        NFS = NFS*signk
        NFkS(i,j,k) = NFS
        NFk (i,j,k) = normalize(NFS)
        Sk  (i,j,k) =    normL2(NFS)
      enddo
    enddo
  enddo
  ! computing cells volumes
  !$OMP DO
  do k=1-gc(5),Nk+gc(6)
    do j=1-gc(3),Nj+gc(4)
      do i=1-gc(1),Ni+gc(2)
        Vx = 0._R_P
        Vy = 0._R_P
        Vz = 0._R_P

        xp = 0.25_R_P*(node(i  ,j  ,k  )%x + node(i  ,j  ,k-1)%x + &
                       node(i  ,j-1,k  )%x + node(i  ,j-1,k-1)%x)
        yp = 0.25_R_P*(node(i  ,j  ,k  )%y + node(i  ,j  ,k-1)%y + &
                       node(i  ,j-1,k  )%y + node(i  ,j-1,k-1)%y)
        zp = 0.25_R_P*(node(i  ,j  ,k  )%z + node(i  ,j  ,k-1)%z + &
                       node(i  ,j-1,k  )%z + node(i  ,j-1,k-1)%z)
        xm = 0.25_R_P*(node(i-1,j  ,k  )%x + node(i-1,j  ,k-1)%x + &
                       node(i-1,j-1,k  )%x + node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(node(i-1,j  ,k  )%y + node(i-1,j  ,k-1)%y + &
                       node(i-1,j-1,k  )%y + node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(node(i-1,j  ,k  )%z + node(i-1,j  ,k-1)%z + &
                       node(i-1,j-1,k  )%z + node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*NFiS(i,j,k)%x - xm*NFiS(i-1,j,k)%x
        Vy = Vy + yp*NFiS(i,j,k)%y - ym*NFiS(i-1,j,k)%y
        Vz = Vz + zp*NFiS(i,j,k)%z - zm*NFiS(i-1,j,k)%z

        xp = 0.25_R_P*(node(i  ,j  ,k  )%x + node(i  ,j  ,k-1)%x + &
                       node(i-1,j  ,k  )%x + node(i-1,j  ,k-1)%x)
        yp = 0.25_R_P*(node(i  ,j  ,k  )%y + node(i  ,j  ,k-1)%y + &
                       node(i-1,j  ,k  )%y + node(i-1,j  ,k-1)%y)
        zp = 0.25_R_P*(node(i  ,j  ,k  )%z + node(i  ,j  ,k-1)%z + &
                       node(i-1,j  ,k  )%z + node(i-1,j  ,k-1)%z)
        xm = 0.25_R_P*(node(i  ,j-1,k  )%x + node(i  ,j-1,k-1)%x + &
                       node(i-1,j-1,k  )%x + node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(node(i  ,j-1,k  )%y + node(i  ,j-1,k-1)%y + &
                       node(i-1,j-1,k  )%y + node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(node(i  ,j-1,k  )%z + node(i  ,j-1,k-1)%z + &
                       node(i-1,j-1,k  )%z + node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*NFjS(i,j,k)%x - xm*NFjS(i,j-1,k)%x
        Vy = Vy + yp*NFjS(i,j,k)%y - ym*NFjS(i,j-1,k)%y
        Vz = Vz + zp*NFjS(i,j,k)%z - zm*NFjS(i,j-1,k)%z

        xp = 0.25_R_P*(node(i  ,j  ,k  )%x + node(i  ,j-1,k  )%x + &
                       node(i-1,j  ,k  )%x + node(i-1,j-1,k  )%x)
        yp = 0.25_R_P*(node(i  ,j  ,k  )%y + node(i  ,j-1,k  )%y + &
                       node(i-1,j  ,k  )%y + node(i-1,j-1,k  )%y)
        zp = 0.25_R_P*(node(i  ,j  ,k  )%z + node(i  ,j-1,k  )%z + &
                       node(i-1,j  ,k  )%z + node(i-1,j-1,k  )%z)
        xm = 0.25_R_P*(node(i  ,j  ,k-1)%x + node(i  ,j-1,k-1)%x + &
                       node(i-1,j  ,k-1)%x + node(i-1,j-1,k-1)%x)
        ym = 0.25_R_P*(node(i  ,j  ,k-1)%y + node(i  ,j-1,k-1)%y + &
                       node(i-1,j  ,k-1)%y + node(i-1,j-1,k-1)%y)
        zm = 0.25_R_P*(node(i  ,j  ,k-1)%z + node(i  ,j-1,k-1)%z + &
                       node(i-1,j  ,k-1)%z + node(i-1,j-1,k-1)%z)

        Vx = Vx + xp*NFkS(i,j,k)%x - xm*NFkS(i,j,k-1)%x
        Vy = Vy + yp*NFkS(i,j,k)%y - ym*NFkS(i,j,k-1)%y
        Vz = Vz + zp*NFkS(i,j,k)%z - zm*NFkS(i,j,k-1)%z

        volume(i,j,k) = max(Vx,Vy,Vz)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_metrics

  !> @brief Subroutine for correcting the metrics of natural (and negative volume) boundary conditions cells.
  subroutine bc_metrics_correction(Ni,Nj,Nk,rcc)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN):: Ni,Nj,Nk           !< Number of cells.
  real(R4P),    intent(IN):: rcc(1:)            !< rcc array.
  logical::                  bc_correct         !< Flag for inquiring if the bc metrics must be corrected.
  logical::                  bc_wall            !< Flag for inquiring if the bc is "wall-type": different corrections must be used.
  real(R8P)::                tm                 !< Tangential metrics parameter (-1 for wall-type bc).
  real(R8P)::                sn                 !< Normal     metrics coefficient correction.
  integer(I4P)::             i,j,k              !< counters.
  integer(I4P), parameter::  wall         = -1  !< Wall boundary condition.
  integer(I4P), parameter::  simmetry     = -2  !< Simmetry boundary condition.
  integer(I4P), parameter::  movingwall   = -10 !< Moving wall boundary condition.
  integer(I4P), parameter::  passivewall  = -11 !< Passive wall boundary condition.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE)                  &
  !$OMP PRIVATE(i,j,k,bc_correct,bc_wall,tm,sn) &
  !$OMP SHARED(Ni,Nj,Nk,rcc,icc,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  ! left i
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      i = nint(rcc(icc(0,j,k)))
      bc_correct = ((i<0).OR.(volume(0,j,k)<(0.2_R_P*volume(1,j,k))))
      bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFiS(1,j,k).dot.NFi(0,j,k))
         NFiS( -1,j,  k  ) = -NFiS(1,j,k)+sn*NFi(0,j,k)
         ! tangential metrics
         NFjS(  0,j  ,k  ) = tm*NFjS(1,j  ,k  )
         NFjS(  0,j-1,k  ) = tm*NFjS(1,j-1,k  )
         NFjS(  0,j  ,k-1) = tm*NFjS(1,j  ,k-1)
         NFjS(  0,j-1,k-1) = tm*NFjS(1,j-1,k-1)

         NFkS(  0,j  ,k  ) = tm*NFkS(1,j  ,k  )
         NFkS(  0,j-1,k  ) = tm*NFkS(1,j-1,k  )
         NFkS(  0,j  ,k-1) = tm*NFkS(1,j  ,k-1)
         NFkS(  0,j-1,k-1) = tm*NFkS(1,j-1,k-1)
         ! volume
         volume(0,j,  k  ) = volume(1,j,k)
      end if
    enddo
  enddo
  ! right i
  !$OMP DO
  do k=1,Nk
    do j=1,Nj
      i = nint(rcc(icc(Ni+1,j,k)))
      bc_correct = ((i<0).OR.(volume(Ni+1,j,k)<(0.2_R_P*volume(Ni,j,k))))
      bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFiS(Ni-1,j,k).dot.NFi(Ni,j,k))
         NFiS(  Ni+1,j,  k  ) = -NFiS(Ni-1,j,k)+sn*NFi(Ni,j,k)
         ! tangential metrics
         NFjS(  Ni+1,j  ,k  ) = tm*NFjS(Ni,j  ,k  )
         NFjS(  Ni+1,j-1,k  ) = tm*NFjS(Ni,j-1,k  )
         NFjS(  Ni+1,j  ,k-1) = tm*NFjS(Ni,j  ,k-1)
         NFjS(  Ni+1,j-1,k-1) = tm*NFjS(Ni,j-1,k-1)

         NFkS(  Ni+1,j  ,k  ) = tm*NFkS(Ni,j  ,k  )
         NFkS(  Ni+1,j-1,k  ) = tm*NFkS(Ni,j-1,k  )
         NFkS(  Ni+1,j  ,k-1) = tm*NFkS(Ni,j  ,k-1)
         NFkS(  Ni+1,j-1,k-1) = tm*NFkS(Ni,j-1,k-1)
         ! volume
         volume(Ni+1,j,  k  ) = volume(Ni,j,k)
      end if
    enddo
  enddo
  ! left j
  !$OMP DO
  do k=1,Nk
    do i=1,Ni
      j = nint(rcc(icc(i,0,k)))
      bc_correct = ((j<0).OR.(volume(i,0,k)<(0.2_R_P*volume(i,1,k))))
      bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFjS(i,1,k).dot.NFj(i,0,k))
         NFjS(  i, -1,k  ) = -NFjS(i,1,k)+sn*NFj(i,0,k)
         ! tangential metrics
         NFiS(  i  ,0,k  ) = tm*NFiS(i  ,1,k  )
         NFiS(  i-1,0,k  ) = tm*NFiS(i-1,1,k  )
         NFiS(  i  ,0,k-1) = tm*NFiS(i  ,1,k-1)
         NFiS(  i-1,0,k-1) = tm*NFiS(i-1,1,k-1)

         NFkS(  i  ,0,k  ) = tm*NFkS(i  ,1,k  )
         NFkS(  i-1,0,k  ) = tm*NFkS(i-1,1,k  )
         NFkS(  i  ,0,k-1) = tm*NFkS(i  ,1,k-1)
         NFkS(  i-1,0,k-1) = tm*NFkS(i-1,1,k-1)
         ! volume
         volume(i,  0,k  ) = volume(i,1,k)
      end if
    enddo
  enddo
  ! right j
  !$OMP DO
  do k=1,Nk
    do i=1,Ni
      j = nint(rcc(icc(i,Nj+1,k)))
      bc_correct = ((j<0).OR.(volume(i,Nj+1,k)<(0.2_R_P*volume(i,Nj,k))))
      bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFjS(i,Nj-1,k).dot.NFj(i,Nj,k))
         NFjS(  i,Nj+1,  k  ) = -NFjS(i,Nj-1,k)+sn*NFj(i,Nj,k)
         ! tangential metrics
         NFiS(  i  ,Nj+1,k  ) = tm*NFiS(i  ,Nj,k  )
         NFiS(  i-1,Nj+1,k  ) = tm*NFiS(i-1,Nj,k  )
         NFiS(  i  ,Nj+1,k-1) = tm*NFiS(i  ,Nj,k-1)
         NFiS(  i-1,Nj+1,k-1) = tm*NFiS(i-1,Nj,k-1)

         NFkS(  i  ,Nj+1,k  ) = tm*NFkS(i  ,Nj,k  )
         NFkS(  i-1,Nj+1,k  ) = tm*NFkS(i-1,Nj,k  )
         NFkS(  i  ,Nj+1,k-1) = tm*NFkS(i  ,Nj,k-1)
         NFkS(  i-1,Nj+1,k-1) = tm*NFkS(i-1,Nj,k-1)
         ! volume
         volume(i,Nj+1,  k  ) = volume(i,Nj,k)
      end if
    enddo
  enddo
  ! left k
  !$OMP DO
  do j=1,Nj
    do i=1,Ni
      k = nint(rcc(icc(i,j,0)))
      bc_correct = ((k<0).OR.(volume(i,j,0)<(0.2_R_P*volume(i,j,1))))
      bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFkS(i,j,1).dot.NFk(i,j,0))
         NFkS(  i,  j, -1) = -NFkS(i,j,1)+sn*NFk(i,j,0)
         ! tangential metrics
         NFiS(  i  ,j  ,0) = tm*NFiS(i  ,j  ,1)
         NFiS(  i-1,j  ,0) = tm*NFiS(i-1,j  ,1)
         NFiS(  i  ,j-1,0) = tm*NFiS(i  ,j-1,1)
         NFiS(  i-1,j-1,0) = tm*NFiS(i-1,j-1,1)

         NFjS(  i  ,j  ,0) = tm*NFjS(i  ,j  ,1)
         NFjS(  i-1,j  ,0) = tm*NFjS(i-1,j  ,1)
         NFjS(  i  ,j-1,0) = tm*NFjS(i  ,j-1,1)
         NFjS(  i-1,j-1,0) = tm*NFjS(i-1,j-1,1)
         ! volume
         volume(i,  j,  0) = volume(i,j,1)
      end if
    enddo
  enddo
  ! right k
  !$OMP DO
  do j=1,Nj
    do i=1,Ni
      k = nint(rcc(icc(i,j,Nk+1)))
      bc_correct = ((k<0).OR.(volume(i,j,Nk+1)<(0.2_R_P*volume(i,j,Nk))))
      bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFkS(i,j,Nk-1).dot.NFk(i,j,Nk))
         NFkS(  i,  j,  Nk+1) = -NFkS(i,j,Nk-1)+sn*NFk(i,j,Nk)
         ! tangential metrics
         NFiS(  i  ,j  ,Nk+1) = tm*NFiS(i  ,j  ,Nk)
         NFiS(  i-1,j  ,Nk+1) = tm*NFiS(i-1,j  ,Nk)
         NFiS(  i  ,j-1,Nk+1) = tm*NFiS(i  ,j-1,Nk)
         NFiS(  i-1,j-1,Nk+1) = tm*NFiS(i-1,j-1,Nk)

         NFjS(  i  ,j  ,Nk+1) = tm*NFjS(i  ,j  ,Nk)
         NFjS(  i-1,j  ,Nk+1) = tm*NFjS(i-1,j  ,Nk)
         NFjS(  i  ,j-1,Nk+1) = tm*NFjS(i  ,j-1,Nk)
         NFjS(  i-1,j-1,Nk+1) = tm*NFjS(i-1,j-1,Nk)
         ! volume
         volume(i,  j,  Nk+1) = volume(i,j,Nk)
      end if
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine bc_metrics_correction

  !> @brief Procedure for computing vorticity variables.
  subroutine compute_vorticity(gc,Ni,Nj,Nk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN)::  gc(1:6)                                                        !< Number of ghost cells.
  integer(I4P), intent(IN)::  Ni,Nj,Nk                                                       !< Number of cells.
  real(R8P)::                 Fi(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)) !< Fluxes i direction.
  real(R8P)::                 Fj(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)) !< Fluxes j direction.
  real(R8P)::                 Fk(1:3,1:3,0-gc(1):Ni+gc(2),0-gc(3):Nj+gc(4),0-gc(5):Nk+gc(6)) !< Fluxes k direction.
  type(Type_Vector)::         um                                                             !< Dummy vector variables.
  real(R8P)::                 c(0:2),emin,emax,eval,fval,mu                                  !< Dummy reals.
  real(R8P), dimension(3,3):: IDEN,G,S,O                                                     !< Matrices.
  integer(I4P)::              i,j,k,ii,jj,kk,iter                                            !< Counters.
  real(R8P), parameter ::     eps6=1d-6, eps9=1d-9                                           !< Tolerances.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Fi   = 0._R8P
  Fj   = 0._R8P
  Fk   = 0._R8P
  vord = 0._R8P
  qfac = 0._R8P
  heli = 0._R8P
  IDEN = 0._R8P
  do i=1,3
     IDEN(i,i) = 1._R8P
  end do
  ! extrapolating momentum from inner cells to ghost cells
  do k=1,Nk
    do j=1,Nj
      momentum(1- gc(1),j,k) = 2._R_P*momentum(0,   j,k)-momentum(1, j,k)
      momentum(Ni+gc(2),j,k) = 2._R_P*momentum(Ni+1,j,k)-momentum(Ni,j,k)
    enddo
  enddo
  do k=1,Nk
    do i=0,Ni+1
      momentum(i,1- gc(3),k) = 2._R_P*momentum(i,0,   k)-momentum(i,1, k)
      momentum(i,Nj+gc(4),k) = 2._R_P*momentum(i,Nj+1,k)-momentum(i,Nj,k)
    enddo
  enddo
  do j=0,Nj+1
    do i=0,Ni+1
      momentum(i,j,1- gc(5)) = 2._R_P*momentum(i,j,0   )-momentum(i,j,1 )
      momentum(i,j,Nk+gc(6)) = 2._R_P*momentum(i,j,Nk+1)-momentum(i,j,Nk)
    enddo
  enddo
  ! computing fluxes
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k,um)      &
  !$OMP SHARED(Ni,Nj,Nk,Fi,Fj,Fk,NFiS,NFjS,NFkS,momentum)
  !$OMP DO
  do k=0,Nk+1
    do j=0,Nj+1
      do i=-1,Ni+1
        um = 0.5_R_P*(momentum(i,j,k)+momentum(i+1,j,k))
        Fi(1,1,i,j,k) = um%x*NFiS(i,j,k)%x
        Fi(1,2,i,j,k) = um%x*NFiS(i,j,k)%y
        Fi(1,3,i,j,k) = um%x*NFiS(i,j,k)%z
        Fi(2,1,i,j,k) = um%y*NFiS(i,j,k)%x
        Fi(2,2,i,j,k) = um%y*NFiS(i,j,k)%y
        Fi(2,3,i,j,k) = um%y*NFiS(i,j,k)%z
        Fi(3,1,i,j,k) = um%z*NFiS(i,j,k)%x
        Fi(3,2,i,j,k) = um%z*NFiS(i,j,k)%y
        Fi(3,3,i,j,k) = um%z*NFiS(i,j,k)%z
      enddo
    enddo
  enddo
  !$OMP DO
  do k=0,Nk+1
    do j=-1,Nj+1
      do i=0,Ni+1
        um = 0.5_R_P*(momentum(i,j,k)+momentum(i,j+1,k))
        Fj(1,1,i,j,k) = um%x*NFjS(i,j,k)%x
        Fj(1,2,i,j,k) = um%x*NFjS(i,j,k)%y
        Fj(1,3,i,j,k) = um%x*NFjS(i,j,k)%z
        Fj(2,1,i,j,k) = um%y*NFjS(i,j,k)%x
        Fj(2,2,i,j,k) = um%y*NFjS(i,j,k)%y
        Fj(2,3,i,j,k) = um%y*NFjS(i,j,k)%z
        Fj(3,1,i,j,k) = um%z*NFjS(i,j,k)%x
        Fj(3,2,i,j,k) = um%z*NFjS(i,j,k)%y
        Fj(3,3,i,j,k) = um%z*NFjS(i,j,k)%z
      enddo
    enddo
  enddo
  !$OMP DO
  do k=-1,Nk+1
    do j=0,Nj+1
      do i=0,Ni+1
        um = 0.5_R_P*(momentum(i,j,k)+momentum(i,j,k+1))
        Fk(1,1,i,j,k) = um%x*NFkS(i,j,k)%x
        Fk(1,2,i,j,k) = um%x*NFkS(i,j,k)%y
        Fk(1,3,i,j,k) = um%x*NFkS(i,j,k)%z
        Fk(2,1,i,j,k) = um%y*NFkS(i,j,k)%x
        Fk(2,2,i,j,k) = um%y*NFkS(i,j,k)%y
        Fk(2,3,i,j,k) = um%y*NFkS(i,j,k)%z
        Fk(3,1,i,j,k) = um%z*NFkS(i,j,k)%x
        Fk(3,2,i,j,k) = um%z*NFkS(i,j,k)%y
        Fk(3,3,i,j,k) = um%z*NFkS(i,j,k)%z
      enddo
    enddo
  enddo
  !$OMP END PARALLEL
  do k=0,Nk+1
    do j=0,Nj+1
      do i=0,Ni+1
         ! computing velocity gradient
         do jj=1,3
            do ii=1,3
               G(ii,jj) = FI(ii,jj,i,j,k)-FI(ii,jj,i-1,j,k)&
                        + FJ(ii,jj,i,j,k)-FJ(ii,jj,i,j-1,k)&
                        + FK(ii,jj,i,j,k)-FK(ii,jj,i,j,k-1)
               G(ii,jj) = G(ii,jj)/max(eps6*eps6,volume(i,j,k))
            enddo
         enddo

         ! computing vorticity vector
         um%x = G(3,2) - G(2,3)
         um%y = G(1,3) - G(3,1) ! control signs!!!!
         um%z = G(2,1) - G(1,2)

         ! tensor S^2 + O^2 (saved in G)
         S = 0._R_P
         O = 0._R_P
         do kk=1,3
            do jj=1,3
               S(jj,kk) = 0.5_R_P*(G(jj,kk)+G(kk,jj))
               O(jj,kk) = 0.5_R_P*(G(jj,kk)-G(kk,jj))
            enddo
         enddo
         G = matmul(S,S) + matmul(O,O)

         ! coefficients of characterist polynomial: lamda^3 + c(2)*lamda^2 + c(1)*lamda + c(0) = 0
         c(2) = -(G(1,1) + G(2,2) + G(3,3))
         c(1) = G(1,1)*G(2,2) + G(1,1)*G(3,3) + G(2,2)*G(3,3) - G(2,3)**2 - G(1,3)**2 - G(1,2)**2
         c(0) = G(1,1)*G(2,3)**2 + G(2,2)*G(1,3)**2 + G(3,3)*G(1,2)**2 - 2._R_P*G(2,3)*G(1,3)*G(1,2) - G(1,1)*G(2,2)*G(3,3)

         ! computing second eigenvalue of characteristic polynomial
         mu = sqrt(c(2)**2 - 3._R_P*c(1))
         emin = (-c(2)-mu)/3._R_P
         emax = (-c(2)+mu)/3._R_P
         do iter=1,100
            eval = 0.5_R_P*(emin+emax)
            fval = eval**3 + c(2)*eval**2 + c(1)*eval + c(0)
            if (fval<0._R_P) then
               emax = eval
            else
               emin = eval
            end if
            if (abs(fval)<eps9 .and.((emax-emin)/eval)<eps6) exit
         end do

         ! saving vordet, qfactor and helicity
         vord(     i,j,k) = eval
         qfac(     i,j,k) = 0.5_R_P*(dot_product(O(1,:),O(1,:)) + dot_product(O(2,:),O(2,:)) + dot_product(O(3,:),O(3,:)) - &
                                    (dot_product(S(1,:),S(1,:)) + dot_product(S(2,:),S(2,:)) + dot_product(S(3,:),S(3,:))))
         heli(     i,j,k) = momentum(i,j,k).dot.um / (max(eps6*eps6,normL2(momentum(i,j,k))*normL2(um)))
         vorticity(i,j,k) = um
      enddo
    enddo
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute_vorticity
endmodule Block_Variables
