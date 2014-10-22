program XnPlot
!-----------------------------------------------------------------------------------------------------------------------------------
! Code Name XnPlot
! Author    Stefano Zaghi
! Version   0.0.1
! Date      14-02-2012
!
! XnPlot loads Xnavis mesh and solution files and produces post-processed plot files.
!
! TODO: Graphical User Interface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                    ! Integers and reals precision definition.
USE Data_Type_Tensor                                                                ! Definition of Type_Tensor.
USE Data_Type_Vector                                                                ! Definition of Type_Vector.
USE Data_Type_OS                                                                    ! Definition of Type_OS.
USE Lib_IO_Misc                                                                     ! Library for miscellanea IO procedures.
USE Lib_VTK_IO                                                                      ! Library for I/O VTK files.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
character(100)::            File_grd ='unset'          ! GRD file name.
character(100)::            File_icc ='unset'          ! ICC file name.
character(100)::            File_sol ='unset'          ! Solution file name.
character(100)::            File_out ='unset'          ! Output file name.
integer(I_P)::              unit_grd                   ! Unit of GRD file.
integer(I_P)::              unit_icc                   ! Unit of ICC file.
integer(I_P)::              unit_sol                   ! Unit of Solution file.
integer(I_P)::              unit_out                   ! Unit of Output file.
logical::                   is_file                    ! Flag for inquiring the presence of file.
integer(I_P)::              i,j,k,b,v                  ! Counters.
integer(I_P)::              err                        ! Error traping flag: 0 no errors, >0 error occours.
integer(I_P)::              nprocs                     ! number of processes
integer(I_P)::              myrank = 0_I_P             ! Rank of current processed file.
integer(I_P), allocatable:: blockmap(:,:),procmap(:,:) ! Blocks maps.
logical::                   global_blk_num = .false.   ! Flag for inquiring the activation of global blocks numeration.
integer(I_P)::              Nbtot                      ! Number of total blocks.
! mesh variables
integer(I_P)::                   Nb                ! Number of blocks.
integer(I_P),      allocatable:: gc(:,:)           ! Number of ghost cells.
integer(I_P),      allocatable:: Ni(:),Nj(:),Nk(:) ! Number of cells of each block.
type(Type_Vector), allocatable:: node(:,:,:)       ! Nodes coordinates.
! metrics variables
type(Type_Vector), allocatable:: NFi(:,:,:)        ! Face i normal versor.
type(Type_Vector), allocatable:: NFj(:,:,:)        ! Face j normal versor.
type(Type_Vector), allocatable:: NFk(:,:,:)        ! Face k normal versor.
type(Type_Vector), allocatable:: NFiS(:,:,:)       ! Face i normal versor with surface area module.
type(Type_Vector), allocatable:: NFjS(:,:,:)       ! Face j normal versor with surface area module.
type(Type_Vector), allocatable:: NFkS(:,:,:)       ! Face k normal versor with surface area module.
real(R_P),         allocatable:: Si(:,:,:)         ! Face i area.
real(R_P),         allocatable:: Sj(:,:,:)         ! Face j area.
real(R_P),         allocatable:: Sk(:,:,:)         ! Face k area.
real(R_P),         allocatable:: volume(:,:,:)     ! Volumes of cells.
! icc variables
integer(I_P),      allocatable:: icc(:,:,:)  ! Cell centered icc values.
integer(I_P),      allocatable:: ricc(:,:,:) ! Cell centered rcc values.
integer(I_P),      allocatable:: vicc(:,:,:) ! Node centered rcc values.
real(R4P),         allocatable:: rcc(:)      ! rcc array.
integer(I_P)::                   Nrcc        ! rcc array dimension.
integer(I_P),      allocatable:: Ncc(:)      ! Number internal chimera cells for each block [1:Nb].
! solution variables
integer(I_P)::                   Nvar = 3        ! Number of variables.
type(Type_Vector), allocatable:: momentum(:,:,:) ! Momentum.
real(R_P),         allocatable:: pressure(:,:,:) ! Pressure.
real(R_P),         allocatable:: f(:,:,:)        ! Level set function.
real(R_P),         allocatable:: f0(:,:,:)       ! Level 0 (level set).
real(R_P),         allocatable:: visc(:,:,:)     ! Viscosity.
real(R_P),         allocatable:: vitl(:,:,:)     ! Turbulent viscosity.
real(R_P),         allocatable:: ken(:,:,:)      ! Turbulent kinetic energy.
real(R_P),         allocatable:: eps(:,:,:)      ! Turbulent kinetic energy dissipation.
real(R_P),         allocatable:: vord(:,:,:)     ! Variable to identify vortices (lambda 2).
real(R_P),         allocatable:: qfac(:,:,:)     ! Variable to identify vortices (q factor).
real(R_P),         allocatable:: heli(:,:,:)     ! Helicity.
type(Type_Vector), allocatable:: vorticity(:,:,:)! Vorticity.
real(R_P),         allocatable:: Fi(:,:,:,:,:)   ! Fluxes i direction.
real(R_P),         allocatable:: Fj(:,:,:,:,:)   ! Fluxes j direction.
real(R_P),         allocatable:: Fk(:,:,:,:,:)   ! Fluxes k direction.
real(R_P)::                      t               ! Solution time.
! misc variables
logical::                 fcc       = .false. ! Inquiring flag for icc file.
logical::                 ngc       = .false. ! Inquiring flag for grd without ghosts cells.
logical::                 sol       = .false. ! Inquiring flag for solution file.
logical::                 cell      = .false. ! Inquiring flag for interpolation (or not) variables at nodes.
logical::                 level_set = .false. ! Inquiring flag for level set variable.
logical::                 laminar   = .false. ! Inquiring flag for no turbulent model.
logical::                 zeroeq    = .false. ! Inquiring flag for zero equations turbulent variables.
logical::                 oneeq     = .true.  ! Inquiring flag for one  equations turbulent variables.
logical::                 twoeq     = .false. ! Inquiring flag for two  equations turbulent variables.
logical::                 vordet    = .false. ! Inquiring flag for vordet variable computing.
logical::                 binary    = .true.  ! Inquiring flag for binary output file.
logical::                 tec       = .true.  ! Inquiring flag for tecplot file format.
logical::                 vtk       = .false. ! Inquiring flag for vtk file format.
integer(I_P), external::  tecini112,    &     ! |
                          tecauxstr112, &     ! |
                          teczne112,    &     ! | Tecplot external binary functions.
                          tecdat112,    &     ! |
                          tecend112           ! |
character(1), parameter:: tecendrec = char(0) ! End-character for binary-record finalize.
character(800)::          tecvarname          ! Variables name for tecplot header file.
character(100)::          vtkbdir             ! VTK base name of output directory.
character(100)::          vtkbfile            ! VTK base name of output files.
type(Type_OS)::           OS                  ! Running architecture.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
! parsing command line arguments
call parse_command_arguments

! opening input files
inquire(file=adjustl(trim(File_grd)),exist=is_file,iostat=err)
if (.NOT.is_file) then
  err = File_Not_Found(filename=File_grd,cpn='XnPlot')
else
  open(unit=Get_Unit(unit_grd),file=adjustl(trim(File_grd)),form='UNFORMATTED',action='READ')
endif
if (fcc) then
  inquire(file=adjustl(trim(File_icc)),exist=is_file,iostat=err)
  if (.NOT.is_file) then
    err = File_Not_Found(filename=File_icc,cpn='XnPlot')
  else
    open(unit=Get_Unit(unit_icc),file=adjustl(trim(File_icc)),form='UNFORMATTED',action='READ')
  endif
endif
if (sol) then
  inquire(file=adjustl(trim(File_sol)),exist=is_file,iostat=err)
  if (.NOT.is_file) then
    err = File_Not_Found(filename=File_sol,cpn='XnPlot')
  else
    open(unit=Get_Unit(unit_sol),file=adjustl(trim(File_sol)),form='UNFORMATTED',action='READ')
  endif
endif

! opening output file
if (tec) then
  if (binary) then
    err = tecini112(tecendrec,                   &
                    trim(tecvarname)//tecendrec, &
                    adjustl(trim(File_out))//".plt"//tecendrec,'.'//tecendrec,0,0,1)
  else
    open(unit=Get_Unit(unit_out),file=adjustl(trim(File_out))//".dat")
    write(unit_out,'(A)',iostat=err)trim(tecvarname)
  endif
endif
!if (vtk) then
!  ! making VTS_files directory
!  err = make_dir(directory=vtkbdir)
!endif

! loading mesh dimensions
read(unit_grd)Nb
if (.not.global_blk_num) then
  ! there is no global blocks numerations: the blocks map is identity
  if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:2,1:Nb)) ; blockmap = 0_I_P
  do b=1,Nb
    blockmap(2,b) = b
  enddo
endif
if (fcc) then
  read(unit_icc)b
  if (b/=Nb) then
    write(stderr,'(A)')' Inconsistent ICC file with GRD one:'
    write(stderr,'(A)')' Nb GRD file: '//trim(str(.true.,Nb))
    write(stderr,'(A)')' Nb ICC file: '//trim(str(.true.,b))
    stop
  endif
endif
if (sol) then
  read(unit_sol)t
endif
write(stdout,'(A)')' Number of block of input files: '//trim(str(.true.,Nb))
if (allocated(gc)) deallocate(gc) ; allocate(gc(1:6,1:Nb))
if (allocated(Ni)) deallocate(Ni) ; allocate(Ni(    1:Nb))
if (allocated(Nj)) deallocate(Nj) ; allocate(Nj(    1:Nb))
if (allocated(Nk)) deallocate(Nk) ; allocate(Nk(    1:Nb))
if (ngc) then
  ! mesh without ghosts cells as the output of geogrd
  gc = 0
else
  ! mesh with ghosts cells as the output of overset
  gc = 2
endif
write(stdout,'(A)')' Blocks dimensions'
do b=1,Nb
  read(unit_grd,err=1)Ni(b),Nj(b),Nk(b),gc(1,b),gc(2,b),gc(3,b),gc(4,b),gc(5,b),gc(6,b)
  1 continue
  write(stdout,'(A)')'   b,Ni,Nk,Nk,gc(6): '//trim(str(.true.,b))//', '//       &
                                              trim(str(.true.,Ni(b)))//', '//   &
                                              trim(str(.true.,Nj(b)))//', '//   &
                                              trim(str(.true.,Nk(b)))//', '//   &
                                              trim(str(.true.,gc(1,b)))//', '// &
                                              trim(str(.true.,gc(2,b)))//', '// &
                                              trim(str(.true.,gc(3,b)))//', '// &
                                              trim(str(.true.,gc(4,b)))//', '// &
                                              trim(str(.true.,gc(5,b)))//', '// &
                                              trim(str(.true.,gc(6,b)))
  if (fcc) then
    read(unit_icc,err=2)i,j,k,v,v,v,v,v,v
    2 continue
    if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
      write(stderr,'(A)')' Inconsistent ICC file with GRD one:'
      write(stderr,'(A)')' Block: '//trim(str(.true.,b))
      write(stderr,'(A)')' Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//trim(str(.true.,Nj(b)))//','//trim(str(.true.,Nk(b)))
      write(stderr,'(A)')' Ni,Nj,Nk ICC file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
      stop
    endif
  endif
  if (sol) then
    read(unit_sol)i,j,k
    if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
      write(stderr,'(A)')' Inconsistent SOL file with GRD one:'
      write(stderr,'(A)')' Block: '//trim(str(.true.,b))
      write(stderr,'(A)')' Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//trim(str(.true.,Nj(b)))//','//trim(str(.true.,Nk(b)))
      write(stderr,'(A)')' Ni,Nj,Nk SOL file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
      stop
    endif
  endif
enddo
if (fcc) then ! loading rcc
  if (allocated(Ncc)) deallocate(Ncc) ; allocate(Ncc(1:Nb)) ; Ncc = 0
  do b=1,Nb
    read(unit_icc)(((v,i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
  enddo
  read(unit_icc)Nrcc
  if (allocated(rcc)) deallocate(rcc) ; allocate(rcc(1:Nrcc))
  read(unit_icc)(rcc(v),v=1,Nrcc)
  ! rewind icc file at icc pointer record
  rewind(unit_icc)
  read(unit_icc)b
  do b=1,Nb
    read(unit_icc,err=3)i,j,k,v,v,v,v,v,v
    3 continue
  enddo
endif

do b=1,Nb ! post processing blocks
  ! allocating mesh variables
  if (allocated(node)) deallocate(node)
  allocate(node(0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
  ! loading block nodes coordinates
  read(unit_grd)(((node(i,j,k)%x,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(unit_grd)(((node(i,j,k)%y,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(unit_grd)(((node(i,j,k)%z,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  if (fcc) then ! allocating and loading icc variables
    if (allocated(icc))  deallocate(icc)
    if (allocated(ricc)) deallocate(ricc)
    if (allocated(vicc)) deallocate(vicc)
                   allocate(icc (1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
                   allocate(ricc(0        :Ni(b)+1      ,0        :Nj(b)+1      ,0        :Nk(b)+1      ))
    if (.not.cell) allocate(vicc(0        :Ni(b)+1      ,0        :Nj(b)+1      ,0        :Nk(b)+1      ))
    read(unit_icc)(((icc(i,j,k),i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
    Ncc(b) = count(icc(1:Ni(b),1:Nj(b),1:Nk(b))>0)
    write(stdout,'(A)')' Block:'//trim(str(.true.,b))//' Number of internal chimera cells: '//trim(str(.true.,Ncc(b)))
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(i,j,k)         &
    !$OMP SHARED(b,Ni,Nj,Nk,icc,rcc,ricc)
    ! computing ricc for blanking
    !$OMP DO
    do k=1,Nk(b)
      do j=1,Nj(b)
        do i=1,Ni(b)
          if (icc(i,j,k)>0) then
            ricc(i,j,k) = nint(rcc(icc(i,j,k)))
          else
            ricc(i,j,k) = 0
          end if
          if (ricc(i,j,k)>28) ricc(i,j,k) = 0
        enddo
      enddo
    enddo
    !$OMP END PARALLEL
    if (.not.cell) then ! ricc must be interpolated at nodes
      ricc(0      ,:      ,:      ) = ricc(1    ,:    ,:    )
      ricc(Ni(b)+1,:      ,:      ) = ricc(Ni(b),:    ,:    )
      ricc(:      ,0      ,:      ) = ricc(:    ,1    ,:    )
      ricc(:      ,Nj(b)+1,:      ) = ricc(:    ,Nj(b),:    )
      ricc(:      ,:      ,0      ) = ricc(:    ,:    ,1    )
      ricc(:      ,:      ,Nk(b)+1) = ricc(:    ,:    ,Nk(b))
      do k=0,Nk(b)
        do j=0,Nj(b)
          do i=0,Ni(b)
            vicc(i,j,k) = max(0,maxval(ricc(i:i+1,j:j+1,k:k+1)))
          enddo
        enddo
      enddo
    endif
  endif
  if (sol) then ! allocating and loading solution variables
    if (allocated(momentum)) deallocate(momentum)
    if (allocated(pressure)) deallocate(pressure)
    allocate(momentum(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
    allocate(pressure(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
    read(unit_sol)(((momentum(i,j,k)%x,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(unit_sol)(((momentum(i,j,k)%y,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(unit_sol)(((momentum(i,j,k)%z,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(unit_sol)(((pressure(i,j,k)  ,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    if (level_set) then
      if (allocated(f))  deallocate(f)  ; allocate(f (1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      if (allocated(f0)) deallocate(f0) ; allocate(f0(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((f (i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(unit_sol)(((f0(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (zeroeq.or.oneeq.or.twoeq) then
      if (allocated(visc)) deallocate(visc)
      allocate(visc(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((visc(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (oneeq) then
      if (allocated(vitl)) deallocate(vitl)
      allocate(vitl(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((vitl(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (twoeq) then
      if (allocated(ken)) deallocate(ken) ; allocate(ken(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      if (allocated(eps)) deallocate(eps) ; allocate(eps(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      read(unit_sol)(((ken(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(unit_sol)(((eps(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (vordet) then ! computing vord variable
      ! allocating metrics variables
      if (allocated(NFi))    deallocate(NFi)
      if (allocated(NFj))    deallocate(NFj)
      if (allocated(NFk))    deallocate(NFk)
      if (allocated(NFiS))   deallocate(NFiS)
      if (allocated(NFjS))   deallocate(NFjS)
      if (allocated(NFkS))   deallocate(NFkS)
      if (allocated(Si))     deallocate(Si)
      if (allocated(Sj))     deallocate(Sj)
      if (allocated(Sk))     deallocate(Sk)
      if (allocated(volume)) deallocate(volume)
      allocate(NFi   (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(NFj   (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(NFk   (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(NFiS  (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(NFjS  (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(NFkS  (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(Si    (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(Sj    (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(Sk    (0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(volume(1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      ! computing the metrics of cells
      call compute_metrics(b=b)
      ! correcting the metrics of boundary conditions cells
      call bc_metrics_correction(b=b)
      ! allocating vordet variables
      if (allocated(vord     )) deallocate(vord     )
      if (allocated(qfac     )) deallocate(qfac     )
      if (allocated(heli     )) deallocate(heli     )
      if (allocated(vorticity)) deallocate(vorticity)
      if (allocated(Fi       )) deallocate(Fi       )
      if (allocated(Fj       )) deallocate(Fj       )
      if (allocated(Fk       )) deallocate(Fk       )
      allocate(vord     (        1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      allocate(qfac     (        1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      allocate(heli     (        1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      allocate(vorticity(        1-gc(1,b):Ni(b)+gc(2,b),1-gc(3,b):Nj(b)+gc(4,b),1-gc(5,b):Nk(b)+gc(6,b)))
      allocate(Fi       (1:3,1:3,0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(Fj       (1:3,1:3,0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      allocate(Fk       (1:3,1:3,0-gc(1,b):Ni(b)+gc(2,b),0-gc(3,b):Nj(b)+gc(4,b),0-gc(5,b):Nk(b)+gc(6,b)))
      vord      = 0._R_P
      qfac      = 0._R_P
      heli      = 0._R_P
      Fi        = 0._R_P
      Fj        = 0._R_P
      Fk        = 0._R_P
      ! computing vordet variable
      call compute_vorticity(b=b)
    endif
  endif
  ! saving the block
  err = block_save(b)
enddo

close(unit_grd)
if (fcc) then
  close(unit_icc)
  if (allocated(rcc)) deallocate(rcc)
endif
if (sol) then
  close(unit_sol)
endif
if (tec) then
  if (binary) then
    err = tecend112()
  else
    close(unit_out)
  endif
endif
if (vtk) then
  ! saving the multiblock VTM wrapper
  err = VTM_INI_XML(filename=adjustl(trim(File_out))//'.vtm')
  err = VTM_BLK_XML(block_action='open')
  err = VTM_WRF_XML(flist=(/(adjustl(trim(OS%basename(File_out)))//'-vts_files'//OS%sep//        &
                             adjustl(trim(OS%basename(File_out)))//"."//trim(strz(4,b))//'.vts', &
                             b=1,Nb)/))
  err = VTM_BLK_XML(block_action='close')
  err = VTM_END_XML()
endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine print_usage()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Printing XnPlot usage.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stdout,*)
  write(stdout,'(A)') ' XnPlot'
  write(stdout,'(A)') ' Post processing code for Xnavis code'
  write(stdout,'(A)') ' Usage:'
  write(stdout,'(A)') '   XnPlot -g file_grd'
  write(stdout,'(A)') '         [-o file_output'
  write(stdout,'(A)') '          -i file_icc'
  write(stdout,'(A)') '          -s file_solution'
  write(stdout,'(A)') '          -ngc'
  write(stdout,'(A)') '          -cell'
  write(stdout,'(A)') '          -ls'
  write(stdout,'(A)') '          -eq #turbulent_eq_model (0,1,2)'
  write(stdout,'(A)') '          -vordet'
  write(stdout,'(A)') '          -ascii'
  write(stdout,'(A)') '          -tec yes/no'
  write(stdout,'(A)') '          -vtk yes/no'
  write(stdout,'(A)') '          -proc #proc'
  write(stdout,'(A)') '          -os UIX/WIN]'
  write(stdout,*)
  write(stdout,'(A)') ' Optional arguments and default values:'
  write(stdout,'(A)') '  -o file_output   => output file name; default is basename of grd file with the proper extension'
  write(stdout,'(A)') '  -i file_icc      => icc file name; if passed the icc variable is saved at cell centers'
  write(stdout,'(A)') '  -s file_solution => solution file name; if passed the solution variables are saved'
  write(stdout,'(A)') '  -ngc             => mesh without ghosts cells, as geogrd output (default no, grd with ghosts cells)'
  write(stdout,'(A)') '  -cell            => all variables other than nodes coord. are cell centered (default no, node centered)'
  write(stdout,'(A)') '  -ls              => solution with level set'
  write(stdout,'(A)') '  -nt              => no turbulent model, laminar flow'
  write(stdout,'(A)') '  -eq #0/1/2       => # equations turbulent model (default 1)'
  write(stdout,'(A)') '  -vordet          => computing "vordet" variable for vortices identification (default no)'
  write(stdout,'(A)') '  -ascii           => write ascii output file (default no, write binary one)'
  write(stdout,'(A)') '  -tec yes/no      => write (or not) Tecplot file format (default yes)'
  write(stdout,'(A)') '  -vtk yes/no      => write (or not) VTK file format (default no)'
  write(stdout,'(A)') '  -proc #proc      => if file "proc.input" if found global blocks numeration is used; #proc is the process'
  write(stdout,'(A)') '                      number of the current processed file'
  write(stdout,'(A)') '  -os UIX/WIN      => type of Operating System write (default *UIX OS type)'
  write(stdout,*)
  write(stdout,'(A)') ' Examples:'
  write(stdout,'(A)') '  XnPlot -g xship.grd                       -o mesh.01 (process only grd file)'
  write(stdout,'(A)') '  XnPlot -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01  (solution variables are saved)'
  write(stdout,*)
  write(stdout,'(A)') ' Notes:'
  write(stdout,'(A)') '   1) the output file name extension is not necessary because it assigned according to the type of output:'
  write(stdout,'(A)') '      binary       Tecplot => .plt'
  write(stdout,'(A)') '      ascii        Tecplot => .dat'
  write(stdout,'(A)') '      binary/ascii VTK     => .vtm'
  write(stdout,'(A)') '   2) if a file name "mb.par" is present into the directory where XnPlot is executed the values of ls'
  write(stdout,'(A)') '      turbulence model are loaded from this file thus they can be omitted from the command'
  write(stdout,'(A)') '      line arguments list.'
  write(stdout,'(A)') '   3) all the variables other than the nodes coordinates are saved at cell center if "-cell" option is used;'
  write(stdout,'(A)') '      if blanking is used the blanking mode must be "any corners" or "primary values".'
  write(stdout,*)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine print_usage

  subroutine parse_command_arguments()
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for parsing command line arguments.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P)::   Nca = 0         ! Number of command line arguments.
  character(8)::   switch          ! Switch string.
  character(100):: string          ! Dummy string.
  character(3)::   os_type = 'UIX' ! OS type 'UIX' or 'WIN'.
  integer(I_P)::   b,c,e           ! Counter.
  integer(I_P)::   unitfree        ! Free logical unit.
  ! turbulence models inquiring flags
  logical:: balom = .false.
  logical:: sgs   = .false.
  logical:: spall = .false.
  logical:: des   = .false.
  logical:: ddes  = .false.
  logical:: lambr = .false.
  logical:: chang = .false.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Nca = command_argument_count()
  if (Nca==0) then
    call print_usage
    stop
  endif
  c = 0
  do while (c<Nca)
    c = c + 1
    call get_command_argument(c,switch)
    select case(adjustl(trim(switch)))
    case('-g') ! grd file
      call get_command_argument(c+1,File_grd) ; c = c + 1
    case('-i') ! icc file
      call get_command_argument(c+1,File_icc) ; c = c + 1
      fcc = .true.
    case('-s') ! solution file
      call get_command_argument(c+1,File_sol) ; c = c + 1
      sol = .true.
    case('-o') ! output file
      call get_command_argument(c+1,File_out) ; c = c + 1
    case('-ngc') ! mesh without ghots cells (as geogrd output)
      ngc = .true.
    case('-cell') ! dependent variables are saved cell centered instead of interpolated at nodes
      cell = .true.
    case('-ls') ! solution with level set
      level_set = .true.
    case('-nt') ! no turbulent model
      laminar = .true.
       zeroeq = .false.
        oneeq = .false.
        twoeq = .false.
    case('-eq') ! turbulent equations model
      call get_command_argument(c+1,string) ; c = c + 1
      read(string,*)e
      select case(e)
      case(0)
        zeroeq = .true.
      case(1)
        oneeq = .true.
      case(2)
        twoeq = .true.
      endselect
    case('-vordet') ! vordet computing
      vordet = .true.
    case('-ascii') ! ascii output
      binary = .false.
    case('-tec') ! Tecplot file format
      call get_command_argument(c+1,string) ; c = c + 1
      string = Upper_Case(adjustl(trim(string)))
      tec = (adjustl(trim(string))=='YES')
    case('-vtk') ! VKT file format
      call get_command_argument(c+1,string) ; c = c + 1
      string = Upper_Case(adjustl(trim(string)))
      vtk = (adjustl(trim(string))=='YES')
    case('-proc') ! global blocks numeration actived
      call get_command_argument(c+1,string) ; c = c + 1
      read(string,*)myrank
      global_blk_num = .true.
    case('-os') ! OS type
      call get_command_argument(c+1,os_type) ; c = c + 1
      os_type = Upper_Case(os_type)
    case default
      write(stderr,'(A)') ' Unknown switch '
      write(stderr,'(A)') adjustl(trim(switch))
      write(stderr,*)
      call print_usage
      stop
    endselect
  enddo
  if ( laminar .and. (zeroeq .or. oneeq .or. twoeq) ) then
     write(stderr,'(A)') ' Incompatible switches, laminar disables turbulent model'
     stop
  end if
  if (vordet.and.(.not.fcc)) then
    write(stderr,'(A)') ' In order to compute "vordet" variable the icc file must be provided because the metrics must be computed.'
    stop
  endif
  ! set OS type
  call OS%init(c_id=os_type)
  if (adjustl(trim(File_grd))=='unset') then
    write(stderr,'(A)') ' File name of grd file must be provided!'
    write(stderr,*)
    call print_usage
    stop
  endif
  if (adjustl(trim(File_out))=='unset') then
    ! the name of grd file is used as output file base name
    File_out = adjustl(trim(File_grd))
  endif
  ! converting the directory separators of files names according to the OS using
  File_grd = adjustl(trim(File_grd)) ; File_grd = OS%string_separator_fix(string=File_grd) ! GRD file name.
  File_icc = adjustl(trim(File_icc)) ; File_icc = OS%string_separator_fix(string=File_icc) ! ICC file name.
  File_sol = adjustl(trim(File_sol)) ; File_sol = OS%string_separator_fix(string=File_sol) ! Solution file name.
  File_out = adjustl(trim(File_out)) ; File_out = OS%string_separator_fix(string=File_out) ! Output file name.
  if (vtk) then
    vtkbdir  = adjustl(trim(OS%basedir(File_out)))//adjustl(trim(OS%basename(File_out)))//'-vts_files'//OS%sep
    vtkbfile = adjustl(trim(vtkbdir))//adjustl(trim(OS%basename(File_out)))
  endif
           write(stdout,'(A)')' GRD File '//trim(File_grd)
  if (fcc) write(stdout,'(A)')' ICC File '//trim(File_icc)
  if (sol) write(stdout,'(A)')' SOL File '//trim(File_sol)
           write(stdout,'(A)')' OUT File '//trim(File_out)
  ! updating number of variables saved and tecplot variables name
  if (binary) then
    tecvarname = 'x y z'
  else
    tecvarname = ' VARIABLES ="x" "y" "z"'
  endif
  if (fcc) then
    Nvar = Nvar + 1
    if (binary) then
      tecvarname = trim(tecvarname)//' icc'
    else
      tecvarname = trim(tecvarname)//' "icc"'
    endif
  endif
  if (sol) then
    Nvar = Nvar + 4
    if (binary) then
      tecvarname = trim(tecvarname)//' u v w p'
    else
      tecvarname = trim(tecvarname)//' "u" "v" "w" "p"'
    endif
    if (zeroeq) then
      Nvar = Nvar + 1 ! visc must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' visc'
      else
        tecvarname = trim(tecvarname)//' "visc"'
      endif
    elseif (oneeq) then
      Nvar = Nvar + 2 ! visc and vitl must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' visc vitl'
      else
        tecvarname = trim(tecvarname)//' "visc" "vitl"'
      endif
    elseif (twoeq) then
      Nvar = Nvar + 3 ! visc, ken and eps must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' visc ken eps'
      else
        tecvarname = trim(tecvarname)//' "visc" "ken" "eps"'
      endif
    endif
    if (level_set) then
      Nvar = Nvar + 2 ! f and f0 must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' f f0'
      else
        tecvarname = trim(tecvarname)//' "f" "f0"'
      endif
    endif
    if (vordet) then
      Nvar = Nvar + 3 ! vorticity variables must be saved
      if (binary) then
        tecvarname = trim(tecvarname)//' vordet qfactor helicity vorticity-x vorticity-y vorticity-z'
      else
        tecvarname = trim(tecvarname)//' "vordet" "qfactor" "helicity" "vorticity-x" "vorticity-y" "vorticity-z"'
      endif
    endif
  endif
  ! inquiring the presence of "mb.par"
  inquire(file="mb.par",exist=is_file,iostat=err)
  if (is_file) then
    open(unit=Get_Unit(unitfree),file="mb.par",action='READ')
    do e=1,8 ! skip the first 8 records
      read(unitfree,*)
    enddo
    read(unitfree,*) e      ! number grid levels
    do e=1,e+9 ! skip the e+9 records
      read(unitfree,*)
    enddo
    read(unitfree,*)        ! Reynolds number
    read(unitfree,*)        ! Foude number
    read(unitfree,*)
    read(unitfree,*) string ! turbulence model
    close(unitfree)
    ! checking turbulence model
    if (string(1:3)=='bal'.or.string(1:3)=='BAL') balom = .true.
    if (string(1:3)=='sgs'.or.string(1:3)=='SGS') sgs   = .true.
    if (string(1:3)=='spa'.or.string(1:3)=='SPA') spall = .true.
    if (string(1:3)=='lam'.or.string(1:3)=='LAM') lambr = .true.
    if (string(1:3)=='cha'.or.string(1:3)=='CHA') chang = .true.
    if (string(1:3)=='des'.or.string(1:3)=='DES') then
       spall = .true.
       des   = .true.
    end if
    if (string(1:3)=='dde'.or.string(1:3)=='DDE') then
       spall = .true.
       des   = .true.
       ddes  = .true.
    end if
    zeroeq = (sgs.OR.balom)
    oneeq  = spall
    twoeq  = (lambr.OR.chang)
    write(stdout,'(A)') ' Found mb.par'
    if (level_set) then
      write(stdout,'(A)') ' Free sruface present'
    else
      write(stdout,'(A)') ' Free sruface not present'
    endif
    if (zeroeq) then
      write(stdout,'(A)') ' Zero equation turbulence model'
    elseif (oneeq) then
      write(stdout,'(A)') ' One equation turbulence model'
    elseif(twoeq) then
      write(stdout,'(A)') ' Two equations turbulence model'
    endif
  endif
  ! inquiring the presence of "proc.input" for global blocks numeration
  if (global_blk_num) then
    inquire(file='proc.input',exist=is_file)
    if (.not.is_file) then
      err = File_Not_Found(filename='proc.input',cpn='parse_command_arguments')
    else
      open(unit=Get_Unit(unitfree),file="proc.input",action='READ')
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) Nbtot ! number of total blocks
      if (allocated(procmap)) deallocate(procmap) ; allocate(procmap(1:2,1:Nbtot)) ; procmap = 0_I_P
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      read(unitfree,*) ! record skipped
      do b=1,Nbtot
        read(unitfree,*) c,procmap(1,b),c,procmap(2,b)
      end do
      close(unitfree)
      nprocs = maxval(procmap(2,:))
      ! computing the local (of myrank) number of blocks
      Nb = 0
      do b=1,Nbtot
        if (procmap(2,b)==myrank) Nb = Nb + 1
      enddo
      if (allocated(blockmap)) deallocate(blockmap) ; allocate(blockmap(1:2,1:Nb)) ; blockmap = 0_I_P
      c = 0
      do b=1,Nbtot
        if (procmap(2,b)==myrank) then
          c = c + 1
          blockmap(1,c) = procmap(1,b)
          blockmap(2,c) = b
        endif
      enddo
    endif
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parse_command_arguments

  subroutine compute_metrics(b)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing block metrics.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: b                 ! Actual block number.
  type(Type_Vector)::        NFS,s1,s2,nd,db   ! Dummy vector variables.
  real(R_P)::                signi,signj,signk ! Dummy variables for checking the directions of normals.
  real(R_P)::                Vx,Vy,Vz          ! Dummy variables for computing volume.
  real(R_P)::                xp,yp,zp          ! Dummy variables for computing face coordinates.
  real(R_P)::                xm,ym,zm          ! Dummy variables for computing face coordinates.
  integer(I_P)::             i,j,k             ! Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  ! computing faces normals
  ! positioning at the middle of the block
  i = max(1,Ni(b)/2)
  j = max(1,Nj(b)/2)
  k = max(1,Nk(b)/2)
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
  !$OMP SHARED(b,Ni,Nj,Nk,gc,signi,signj,signk,node,Si,Sj,Sk,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  !$OMP DO
  do k=1-gc(5,b),Nk(b)+gc(6,b)
    do j=1-gc(3,b),Nj(b)+gc(4,b)
      do i=0-gc(1,b),Ni(b)+gc(2,b)
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
  do k=1-gc(5,b),Nk(b)+gc(6,b)
    do j=0-gc(3,b),Nj(b)+gc(4,b)
      do i=1-gc(1,b),Ni(b)+gc(2,b)
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
  do k=0-gc(5,b),Nk(b)+gc(6,b)
    do j=1-gc(3,b),Nj(b)+gc(4,b)
      do i=1-gc(1,b),Ni(b)+gc(2,b)
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
  do k=1-gc(5,b),Nk(b)+gc(6,b)
    do j=1-gc(3,b),Nj(b)+gc(4,b)
      do i=1-gc(1,b),Ni(b)+gc(2,b)
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

  subroutine bc_metrics_correction(b)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for correcting the metrics of natural (and negative volume) boundary conditions cells.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: b          ! Actual block number.
  logical::                  bc_correct ! Flag for inquiring if the bc metrics must be corrected.
  logical::                  bc_wall    ! Flag for inquiring if the bc is "wall-type": different corrections must be used.
  real(R_P)::                tm         ! Tangential metrics parameter (-1 for wall-type bc).
  real(R_P)::                sn         ! Normal     metrics coefficient correction.
  integer(I_P)::             i,j,k      ! counters.
  integer(I_P), parameter::  wall         = -1
  integer(I_P), parameter::  simmetry     = -2
  integer(I_P), parameter::  movingwall   = -10
  integer(I_P), parameter::  passivewall  = -11
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  !$OMP PARALLEL DEFAULT(NONE)                  &
  !$OMP PRIVATE(i,j,k,bc_correct,bc_wall,tm,sn) &
  !$OMP SHARED(b,Ni,Nj,Nk,rcc,icc,NFi,NFj,NFk,NFiS,NFjS,NFkS,volume)
  ! left i
  !$OMP DO
  do k=1,Nk(b)
    do j=1,Nj(b)
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
  do k=1,Nk(b)
    do j=1,Nj(b)
      i = nint(rcc(icc(Ni(b)+1,j,k)))
      bc_correct = ((i<0).OR.(volume(Ni(b)+1,j,k)<(0.2_R_P*volume(Ni(b),j,k))))
      bc_wall    = ((i==wall).OR.(i==simmetry).OR.(i==movingwall).OR.(i==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFiS(Ni(b)-1,j,k).dot.NFi(Ni(b),j,k))
         NFiS(  Ni(b)+1,j,  k  ) = -NFiS(Ni(b)-1,j,k)+sn*NFi(Ni(b),j,k)
         ! tangential metrics
         NFjS(  Ni(b)+1,j  ,k  ) = tm*NFjS(Ni(b),j  ,k  )
         NFjS(  Ni(b)+1,j-1,k  ) = tm*NFjS(Ni(b),j-1,k  )
         NFjS(  Ni(b)+1,j  ,k-1) = tm*NFjS(Ni(b),j  ,k-1)
         NFjS(  Ni(b)+1,j-1,k-1) = tm*NFjS(Ni(b),j-1,k-1)

         NFkS(  Ni(b)+1,j  ,k  ) = tm*NFkS(Ni(b),j  ,k  )
         NFkS(  Ni(b)+1,j-1,k  ) = tm*NFkS(Ni(b),j-1,k  )
         NFkS(  Ni(b)+1,j  ,k-1) = tm*NFkS(Ni(b),j  ,k-1)
         NFkS(  Ni(b)+1,j-1,k-1) = tm*NFkS(Ni(b),j-1,k-1)
         ! volume
         volume(Ni(b)+1,j,  k  ) = volume(Ni(b),j,k)
      end if
    enddo
  enddo
  ! left j
  !$OMP DO
  do k=1,Nk(b)
    do i=1,Ni(b)
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
  do k=1,Nk(b)
    do i=1,Ni(b)
      j = nint(rcc(icc(i,Nj(b)+1,k)))
      bc_correct = ((j<0).OR.(volume(i,Nj(b)+1,k)<(0.2_R_P*volume(i,Nj(b),k))))
      bc_wall    = ((j==wall).OR.(j==simmetry).OR.(j==movingwall).OR.(j==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFjS(i,Nj(b)-1,k).dot.NFj(i,Nj(b),k))
         NFjS(  i,Nj(b)+1,  k  ) = -NFjS(i,Nj(b)-1,k)+sn*NFj(i,Nj(b),k)
         ! tangential metrics
         NFiS(  i  ,Nj(b)+1,k  ) = tm*NFiS(i  ,Nj(b),k  )
         NFiS(  i-1,Nj(b)+1,k  ) = tm*NFiS(i-1,Nj(b),k  )
         NFiS(  i  ,Nj(b)+1,k-1) = tm*NFiS(i  ,Nj(b),k-1)
         NFiS(  i-1,Nj(b)+1,k-1) = tm*NFiS(i-1,Nj(b),k-1)

         NFkS(  i  ,Nj(b)+1,k  ) = tm*NFkS(i  ,Nj(b),k  )
         NFkS(  i-1,Nj(b)+1,k  ) = tm*NFkS(i-1,Nj(b),k  )
         NFkS(  i  ,Nj(b)+1,k-1) = tm*NFkS(i  ,Nj(b),k-1)
         NFkS(  i-1,Nj(b)+1,k-1) = tm*NFkS(i-1,Nj(b),k-1)
         ! volume
         volume(i,Nj(b)+1,  k  ) = volume(i,Nj(b),k)
      end if
    enddo
  enddo
  ! left k
  !$OMP DO
  do j=1,Nj(b)
    do i=1,Ni(b)
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
  do j=1,Nj(b)
    do i=1,Ni(b)
      k = nint(rcc(icc(i,j,Nk(b)+1)))
      bc_correct = ((k<0).OR.(volume(i,j,Nk(b)+1)<(0.2_R_P*volume(i,j,Nk(b)))))
      bc_wall    = ((k==wall).OR.(k==simmetry).OR.(k==movingwall).OR.(k==passivewall))
      tm = 1._R_P
      if (bc_wall) tm = -1._R_P
      if (bc_correct) then
         ! normal metrics
         sn = 2._R_P*(NFkS(i,j,Nk(b)-1).dot.NFk(i,j,Nk(b)))
         NFkS(  i,  j,  Nk(b)+1) = -NFkS(i,j,Nk(b)-1)+sn*NFk(i,j,Nk(b))
         ! tangential metrics
         NFiS(  i  ,j  ,Nk(b)+1) = tm*NFiS(i  ,j  ,Nk(b))
         NFiS(  i-1,j  ,Nk(b)+1) = tm*NFiS(i-1,j  ,Nk(b))
         NFiS(  i  ,j-1,Nk(b)+1) = tm*NFiS(i  ,j-1,Nk(b))
         NFiS(  i-1,j-1,Nk(b)+1) = tm*NFiS(i-1,j-1,Nk(b))

         NFjS(  i  ,j  ,Nk(b)+1) = tm*NFjS(i  ,j  ,Nk(b))
         NFjS(  i-1,j  ,Nk(b)+1) = tm*NFjS(i-1,j  ,Nk(b))
         NFjS(  i  ,j-1,Nk(b)+1) = tm*NFjS(i  ,j-1,Nk(b))
         NFjS(  i-1,j-1,Nk(b)+1) = tm*NFjS(i-1,j-1,Nk(b))
         ! volume
         volume(i,  j,  Nk(b)+1) = volume(i,j,Nk(b))
      end if
    enddo
  enddo
  !$OMP END PARALLEL
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine bc_metrics_correction

  subroutine compute_vorticity(b)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for computing vorticity variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  b                 ! Actual block number.
  type(Type_Vector)::         um                ! Dummy vector variables.
  real(R_P)::                 c(0:2),emin,emax,eval,fval,mu
  real(R_P), dimension(3,3):: IDEN,G,S,O
  integer(I_P)::              i,j,k,ii,jj,kk,iter    ! Counters.
  real(R_P), parameter :: eps6=1d-6, eps9=1d-9
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  IDEN = 0._R_P
  do i=1,3
     IDEN(i,i) = 1._R_P
  end do
  ! extrapolating momentum from inner cells to ghost cells
  do k=1,Nk(b)
    do j=1,Nj(b)
      momentum(1-    gc(1,b),j,k) = 2._R_P*momentum(0,      j,k)-momentum(1,    j,k)
      momentum(Ni(b)+gc(2,b),j,k) = 2._R_P*momentum(Ni(b)+1,j,k)-momentum(Ni(b),j,k)
    enddo
  enddo
  do k=1,Nk(b)
    do i=0,Ni(b)+1
      momentum(i,1-    gc(3,b),k) = 2._R_P*momentum(i,0,      k)-momentum(i,1,    k)
      momentum(i,Nj(b)+gc(4,b),k) = 2._R_P*momentum(i,Nj(b)+1,k)-momentum(i,Nj(b),k)
    enddo
  enddo
  do j=0,Nj(b)+1
    do i=0,Ni(b)+1
      momentum(i,j,1-    gc(5,b)) = 2._R_P*momentum(i,j,0      )-momentum(i,j,1    )
      momentum(i,j,Nk(b)+gc(6,b)) = 2._R_P*momentum(i,j,Nk(b)+1)-momentum(i,j,Nk(b))
    enddo
  enddo
  ! computing fluxes
  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(i,j,k,um)      &
  !$OMP SHARED(b,Ni,Nj,Nk,Fi,Fj,Fk,NFiS,NFjS,NFkS,momentum)
  !$OMP DO
  do k=0,Nk(b)+1
    do j=0,Nj(b)+1
      do i=-1,Ni(b)+1
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
  do k=0,Nk(b)+1
    do j=-1,Nj(b)+1
      do i=0,Ni(b)+1
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
  do k=-1,Nk(b)+1
    do j=0,Nj(b)+1
      do i=0,Ni(b)+1
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
  do k=0,Nk(b)+1
    do j=0,Nj(b)+1
      do i=0,Ni(b)+1
         ! computing velocity gradient
         do jj=1,3
            do ii=1,3
               G(ii,jj) = FI(ii,jj,i,j,k)-FI(ii,jj,i-1,j,k)&
                        + FJ(ii,jj,i,j,k)-FJ(ii,jj,i,j-1,k)&
                        + FK(ii,jj,i,j,k)-FK(ii,jj,i,j,k-1)
               G(ii,jj) = G(ii,jj)/volume(i,j,k)
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

  function block_save(b) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Function for saving the patch.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: b                 ! Actual block number.
  integer(I_P)::             err               ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::             tecnull(1:Nvar)   ! Tecplot null array.
  integer(I_P)::             tecvarloc(1:Nvar) ! Tecplot array of variables location.
  integer(I_P)::             nnode,ncell       ! Number of nodes and cells.
  real(R_P), allocatable::   vari(:,:,:)       ! Interpolated generic variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.cell) then ! variables must be interpolated at nodes, allocating dummy variable
    if (allocated(vari)) deallocate(vari) ; allocate(vari(0:Ni(b),0:Nj(b),0:Nk(b)))
  endif
  nnode = (Ni(b)+1)*(Nj(b)+1)*(Nk(b)+1) ! computing number of nodes
  ncell =  Ni(b)   * Nj(b)   * Nk(b)    ! computing number of cells
  if (tec) then
    if (binary) then
      tecnull = 0
      tecvarloc = 1
      if (Nvar>3.and.cell) tecvarloc(4:Nvar)= 0
      err = teczne112('blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//tecendrec, &
                      0,                                                                                    &
                      Ni(b)+1,                                                                              &
                      Nj(b)+1,                                                                              &
                      Nk(b)+1,                                                                              &
                      0,                                                                                    &
                      0,                                                                                    &
                      0,                                                                                    &
                      0.0,                                                                                  &
                      0,                                                                                    &
                      0,                                                                                    &
                      1,                                                                                    & !1=>block,0=>point
                      0,                                                                                    &
                      0,                                                                                    &
                      0,                                                                                    &
                      0,                                                                                    &
                      0,                                                                                    &
                      tecnull,                                                                              &
                      tecvarloc,                                                                            &
                      tecnull,                                                                              &
                      0)
    else
      if (Nvar>3.and.cell) then
        write(unit_out,'(A)',iostat=err)' ZONE  T="blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//'"'// &
                                        ', I='//trim(str(no_sign=.true.,n=Ni(b)+1))//                                              &
                                        ', J='//trim(str(no_sign=.true.,n=Nj(b)+1))//                                              &
                                        ', K='//trim(str(no_sign=.true.,n=Nk(b)+1))//                                              &
                                        ', DATAPACKING=BLOCK'//                                                                    &
                                        ', VARLOCATION=([1-3]=NODAL,[4-'//trim(str(.true.,Nvar))//']=CELLCENTERED)'
      else
        write(unit_out,'(A)',iostat=err)' ZONE  T="blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//'"'// &
                                        ', I='//trim(str(no_sign=.true.,n=Ni(b)+1))//                                              &
                                        ', J='//trim(str(no_sign=.true.,n=Nj(b)+1))//                                              &
                                        ', K='//trim(str(no_sign=.true.,n=Nk(b)+1))//                                              &
                                        ', DATAPACKING=BLOCK'//                                                                    &
                                        ', VARLOCATION=([1-'//trim(str(.true.,Nvar))//']=NODAL)'
      endif
    endif
    err = tec_var(n=nnode,var=node(0:Ni(b),0:Nj(b),0:Nk(b))%x,d=1)
    err = tec_var(n=nnode,var=node(0:Ni(b),0:Nj(b),0:Nk(b))%y,d=1)
    err = tec_var(n=nnode,var=node(0:Ni(b),0:Nj(b),0:Nk(b))%z,d=1)
    if (fcc) then
      if (cell) then
        err = tec_var(n=ncell,var=real(ricc(1:Ni(b),1:Nj(b),1:Nk(b)),R_P),d=1)
      else
        err = tec_var(n=nnode,var=real(vicc(0:Ni(b),0:Nj(b),0:Nk(b)),R_P),d=1)
      endif
    endif
    if (sol) then
      if (cell) then
        err = tec_var(n=ncell,var=momentum(1:Ni(b),1:Nj(b),1:Nk(b))%x,d=1)
        err = tec_var(n=ncell,var=momentum(1:Ni(b),1:Nj(b),1:Nk(b))%y,d=1)
        err = tec_var(n=ncell,var=momentum(1:Ni(b),1:Nj(b),1:Nk(b))%z,d=1)
        err = tec_var(n=ncell,var=pressure(1:Ni(b),1:Nj(b),1:Nk(b))  ,d=1)
        if (zeroeq) then
          err = tec_var(n=ncell,var=visc(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
        elseif (oneeq) then
          err = tec_var(n=ncell,var=visc(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
          err = tec_var(n=ncell,var=vitl(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
        elseif (twoeq) then
          err = tec_var(n=ncell,var=visc(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
          err = tec_var(n=ncell,var=ken (1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
          err = tec_var(n=ncell,var=eps (1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
        endif
        if (level_set) then
          err = tec_var(n=ncell,var=f (1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
          err = tec_var(n=ncell,var=f0(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
        endif
        if (vordet) then
          err = tec_var(n=ncell,var=vord(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
          err = tec_var(n=ncell,var=qfac(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
          err = tec_var(n=ncell,var=heli(1:Ni(b),1:Nj(b),1:Nk(b)),d=1)
          err = tec_var(n=ncell,var=vorticity(1:Ni(b),1:Nj(b),1:Nk(b))%x,d=1)
          err = tec_var(n=ncell,var=vorticity(1:Ni(b),1:Nj(b),1:Nk(b))%y,d=1)
          err = tec_var(n=ncell,var=vorticity(1:Ni(b),1:Nj(b),1:Nk(b))%z,d=1)
        endif
      else
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=momentum(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%x,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=momentum(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%y,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=momentum(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%z,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=pressure(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)  ,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        if (zeroeq) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=visc(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        elseif (oneeq) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=visc(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vitl(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        elseif (twoeq) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=visc(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=ken (0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=eps (0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        endif
        if (level_set) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=f (0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=f0(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        endif
        if (vordet) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vord(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=qfac(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=heli(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vorticity(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%x,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vorticity(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%y,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vorticity(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%z,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = tec_var(n=nnode,var=vari(0:Ni(b),0:Nj(b),0:Nk(b)),d=1)
        endif
      endif
    endif
  endif
  if (vtk) then
    if (binary) then
      err = VTK_INI_XML(output_format = 'binary',                                                             &
                        filename      = adjustl(trim(vtkbfile))//".blk"//trim(strz(4,blockmap(2,b)))//        &
                                                                 '.grp'//trim(strz(3,blockmap(1,b)))//'.vts', &
                        mesh_topology = 'StructuredGrid',                                                     &
                        nx1 = 0, nx2 = Ni(b),                                                                 &
                        ny1 = 0, ny2 = Nj(b),                                                                 &
                        nz1 = 0, nz2 = Nk(b))
    else
      err = VTK_INI_XML(output_format = 'ascii',                                                              &
                        filename      = adjustl(trim(vtkbfile))//".blk"//trim(strz(4,blockmap(2,b)))//        &
                                                                 '.grp'//trim(strz(3,blockmap(1,b)))//'.vts', &
                        mesh_topology = 'StructuredGrid',                                                     &
                        nx1 = 0, nx2 = Ni(b),                                                                 &
                        ny1 = 0, ny2 = Nj(b),                                                                 &
                        nz1 = 0, nz2 = Nk(b))
    endif
    err = VTK_GEO_XML(nx1 = 0, nx2 = Ni(b),                                 &
                      ny1 = 0, ny2 = Nj(b),                                 &
                      nz1 = 0, nz2 = Nk(b),                                 &
                      NN = nnode,                                           &
                      X=reshape(node(0:Ni(b),0:Nj(b),0:Nk(b))%x,(/nnode/)), &
                      Y=reshape(node(0:Ni(b),0:Nj(b),0:Nk(b))%y,(/nnode/)), &
                      Z=reshape(node(0:Ni(b),0:Nj(b),0:Nk(b))%z,(/nnode/)))
    if (cell) then
      err = VTK_DAT_XML(var_location = 'cell', var_block_action = 'open')
    else
      err = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
    endif
    if (fcc) then
      if (cell) then
        err = VTK_VAR_XML(NC_NN=ncell,varname='icc',var=reshape(ricc(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
      else
        err = VTK_VAR_XML(NC_NN=nnode,varname='icc',var=reshape(vicc(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
      endif
    endif
    if (sol) then
      if (cell) then
        err = VTK_VAR_XML(NC_NN=ncell,varname='u',var=reshape(momentum(1:Ni(b),1:Nj(b),1:Nk(b))%X,(/ncell/)))
        err = VTK_VAR_XML(NC_NN=ncell,varname='v',var=reshape(momentum(1:Ni(b),1:Nj(b),1:Nk(b))%y,(/ncell/)))
        err = VTK_VAR_XML(NC_NN=ncell,varname='w',var=reshape(momentum(1:Ni(b),1:Nj(b),1:Nk(b))%z,(/ncell/)))
        err = VTK_VAR_XML(NC_NN=ncell,varname='p',var=reshape(pressure(1:Ni(b),1:Nj(b),1:Nk(b))  ,(/ncell/)))
        if (zeroeq) then
          err = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=reshape(visc(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
        elseif (oneeq) then
          err = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=reshape(visc(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='vitl',var=reshape(visc(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
        elseif (twoeq) then
          err = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=reshape(visc(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='ken', var=reshape(ken (1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='eps', var=reshape(eps (1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
        endif
        if (level_set) then
          err = VTK_VAR_XML(NC_NN=ncell,varname='f', var=reshape(f (1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='f0',var=reshape(f0(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
        endif
        if (vordet) then
          err = VTK_VAR_XML(NC_NN=ncell,varname='vordet',     var=reshape(vord(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='qfactor',    var=reshape(qfac(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='helicity',   var=reshape(heli(1:Ni(b),1:Nj(b),1:Nk(b)),(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='vorticity-x',var=reshape(vorticity(1:Ni(b),1:Nj(b),1:Nk(b))%X,(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='vorticity-y',var=reshape(vorticity(1:Ni(b),1:Nj(b),1:Nk(b))%Y,(/ncell/)))
          err = VTK_VAR_XML(NC_NN=ncell,varname='vorticity-z',var=reshape(vorticity(1:Ni(b),1:Nj(b),1:Nk(b))%Z,(/ncell/)))
        endif
      else
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=momentum(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%x,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = VTK_VAR_XML(NC_NN=nnode,varname='u',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=momentum(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%y,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = VTK_VAR_XML(NC_NN=nnode,varname='v',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=momentum(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%z,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = VTK_VAR_XML(NC_NN=nnode,varname='w',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                              var=pressure(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)  ,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
        err = VTK_VAR_XML(NC_NN=nnode,varname='p',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        if (zeroeq) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=visc(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='visc',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        elseif (oneeq) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=visc(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='visc',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vitl(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='vitl',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        elseif (twoeq) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=visc(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='visc',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=ken (0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='ken', var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=eps (0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='eps', var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        endif
        if (level_set) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=f (0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='f', var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=f0(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='f0',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        endif
        if (vordet) then
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vord(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='vordet',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=qfac(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='qfactor',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=heli(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1),vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='helicity',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vorticity(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%X,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='vorticity-x',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vorticity(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%Y,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='vorticity-y',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
          call varinterpolation(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b), &
                                var=vorticity(0:Ni(b)+1,0:Nj(b)+1,0:Nk(b)+1)%Z,vari=vari(0:Ni(b),0:Nj(b),0:Nk(b)))
          err = VTK_VAR_XML(NC_NN=nnode,varname='vorticity-z',var=reshape(vari(0:Ni(b),0:Nj(b),0:Nk(b)),(/nnode/)))
        endif
      endif
    endif
    if (cell) then
      err = VTK_DAT_XML(var_location = 'cell', var_block_action = 'close')
    else
      err = VTK_DAT_XML(var_location = 'node', var_block_action = 'close')
    endif
    err = VTK_GEO_XML()
    err = VTK_END_XML()
  endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction block_save

  subroutine varinterpolation(Ni,Nj,Nk,var,vari)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Subroutine for interpolating celle centered variable into node centered one.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN)::  Ni,Nj,Nk                   ! Block dimensions.
  real(R_P),    intent(IN)::  var (0:Ni+1,0:Nj+1,0:Nk+1) ! Cell centered variable.
  real(R_P),    intent(OUT):: vari(0:Ni  ,0:Nj  ,0:Nk  ) ! Node centered interpolated variable.
  integer(I_P)::              i,j,k                      ! Counters.
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

  function tec_var(n,var,d) result(err)
  !---------------------------------------------------------------------------------------------------------------------------------
  ! Interface function for saving variables into Tecplot file.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I_P), intent(IN):: n        ! Number of var elements.
  real(R_P),    intent(IN):: var(1:n) ! Variable to be saved.
  integer(I_P), intent(IN):: d        ! Double precision output (1 yes, 0 no).
  integer(I_P)::             err      ! Error traping flag: 0 no errors, >0 error occours.
  integer(I_P)::             e        ! Element counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (binary) then
    err = tecdat112(n,var,d)
  else
    write(unit_out,FR_P,iostat=err)(var(e),e=1,n)
  endif
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction tec_var
endprogram XnPlot
