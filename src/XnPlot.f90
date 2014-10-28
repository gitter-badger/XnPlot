!> @addtogroup Program Programs
!> List of excutable programs.
!> @addtogroup GlobalVarPar Global Variables and Parameters
!> List of global variables and parameters.
!> @addtogroup PrivateVarPar Private Variables and Parameters
!> List of private variables and parameters.
!> @addtogroup PrivateProcedure Private Procedures
!> List of private procedures.

!> @ingroup Program
!> @{
!> @defgroup XnPlotProgram XnPlot
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup XnPlotPrivateVarPar XnPlot
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup XnPlotPrivateProcedure XnPlot
!> @}

!> @brief XnPlot loads Xnavis mesh and solution files and produces post-processed plot files.
!> @author    Stefano Zaghi
!> @version   0.0.3
!> @date      2014-10-28
!> @copyright GNU Public License version 3.
!> @ingroup XnPlotProgram
program XnPlot
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                    ! Integers and reals precision definition.
USE Block_Variables                                                                 ! Block variables definition.
USE Data_Type_Command_Line_Interface                                                ! Definition of Type_Command_Line_Interface.
USE Data_Type_PostProcess                                                           ! Definition of Type_PostProcess.
USE Lib_Metrics                                                                     ! Library for metrics computations.
USE Lib_Vorticity                                                                   ! Library for vorticity computations.
USE Lib_TEC                                                                         ! Library for I/O Tecplot files.
!USE Lib_VTK_IO                                                                      ! Library for I/O VTK files.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout => OUTPUT_UNIT, stderr => ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
type(Type_PostProcess)::    pp                !< Post-processor data.
integer(I4P)::              Nb                !< Number of blocks.
integer(I4P), allocatable:: gc(:,:)           !< Number of ghost cells [1:6,1:Nb].
integer(I4P), allocatable:: Ni(:),Nj(:),Nk(:) !< Number of cells of each block [1:Nb].
real(R4P),    allocatable:: rcc(:)            !< rcc array.
integer(I4P)::              Nrcc              !< rcc array dimension.
integer(I4P), allocatable:: Ncc(:)            !< Number internal chimera cells for each block [1:Nb].
integer(I4P)::              i,j,k,b,v         !< Counters.
real(R8P)::                 t                 !< Solution time.
integer(I4P)::              error             !< Error trapping flag.
! misc variables
!character(100)::          vtkbdir             ! VTK base name of output directory.
!character(100)::          vtkbfile            ! VTK base name of output files.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
call parse_command_arguments(pp=pp)

!if (vtk) then
!  ! making VTS_files directory
!  err = make_dir(directory=vtkbdir)
!endif

! loading headers
read(pp%unit_grd)Nb
if (.not.pp%global_blk_num) then
  call pp%set_blockmap(Nb=Nb)
endif
if (pp%fcc) then
  read(pp%unit_icc)b
  if (b/=Nb) then
    write(stderr,'(A)')'+--> Inconsistent ICC file with GRD one:'
    write(stderr,'(A)')'|--> Nb GRD file: '//trim(str(.true.,Nb))
    write(stderr,'(A)')'|--> Nb ICC file: '//trim(str(.true.,b))
    stop
  endif
endif
if (pp%sol) then
  read(pp%unit_sol)t
endif
write(stdout,'(A)')'+--> Number of block of input files: '//trim(str(.true.,Nb))
if (allocated(gc)) deallocate(gc) ; allocate(gc(1:6,1:Nb))
if (allocated(Ni)) deallocate(Ni) ; allocate(Ni(    1:Nb))
if (allocated(Nj)) deallocate(Nj) ; allocate(Nj(    1:Nb))
if (allocated(Nk)) deallocate(Nk) ; allocate(Nk(    1:Nb))
if (pp%ngc) then ! mesh without ghosts cells as the output of geogrd
  gc = 0
else ! mesh with ghosts cells as the output of overset
  gc = 2
endif
if (pp%verbose) write(stdout,'(A)')'+--> Blocks dimensions'
do b=1,Nb
  read(pp%unit_grd,err=1)Ni(b),Nj(b),Nk(b),gc(1,b),gc(2,b),gc(3,b),gc(4,b),gc(5,b),gc(6,b)
  1 continue
  if (pp%verbose) write(stdout,'(A)')'|-->   b,Ni,Nk,Nk,gc(6): '//trim(str(.true.,b))//', '//       &
                                                                  trim(str(.true.,Ni(b)))//', '//   &
                                                                  trim(str(.true.,Nj(b)))//', '//   &
                                                                  trim(str(.true.,Nk(b)))//', '//   &
                                                                  trim(str(.true.,gc(1,b)))//', '// &
                                                                  trim(str(.true.,gc(2,b)))//', '// &
                                                                  trim(str(.true.,gc(3,b)))//', '// &
                                                                  trim(str(.true.,gc(4,b)))//', '// &
                                                                  trim(str(.true.,gc(5,b)))//', '// &
                                                                  trim(str(.true.,gc(6,b)))
  if (pp%fcc) then
    read(pp%unit_icc,err=2)i,j,k,v,v,v,v,v,v
    2 continue
    if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
      write(stderr,'(A)')'+--> Inconsistent ICC file with GRD one:'
      write(stderr,'(A)')'|--> Block: '//trim(str(.true.,b))
      write(stderr,'(A)')'|--> Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//&
                                                     trim(str(.true.,Nj(b)))//','//&
                                                     trim(str(.true.,Nk(b)))
      write(stderr,'(A)')'|--> Ni,Nj,Nk ICC file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
      stop
    endif
  endif
  if (pp%sol) then
    read(pp%unit_sol)i,j,k
    if (i/=Ni(b).OR.j/=Nj(b).OR.k/=Nk(b)) then
      write(stderr,'(A)')'+--> Inconsistent SOL file with GRD one:'
      write(stderr,'(A)')'|--> Block: '//trim(str(.true.,b))
      write(stderr,'(A)')'|--> Ni,Nj,Nk GRD file: '//trim(str(.true.,Ni(b)))//','//&
                                                     trim(str(.true.,Nj(b)))//','//&
                                                     trim(str(.true.,Nk(b)))
      write(stderr,'(A)')'|--> Ni,Nj,Nk SOL file: '//trim(str(.true.,i))//','//trim(str(.true.,j))//','//trim(str(.true.,k))
      stop
    endif
  endif
enddo
if (pp%fcc) then ! loading rcc
  if (allocated(Ncc)) deallocate(Ncc) ; allocate(Ncc(1:Nb)) ; Ncc = 0
  do b=1,Nb
    read(pp%unit_icc)(((v,i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
  enddo
  read(pp%unit_icc)Nrcc
  if (allocated(rcc)) deallocate(rcc) ; allocate(rcc(1:Nrcc))
  read(pp%unit_icc)(rcc(v),v=1,Nrcc)
  ! rewind icc file at icc pointer record
  rewind(pp%unit_icc)
  read(pp%unit_icc)b
  do b=1,Nb
    read(pp%unit_icc,err=3)i,j,k,v,v,v,v,v,v
    3 continue
  enddo
endif

do b=1,Nb ! post processing blocks
  call block_allocate(pp=pp,gc=gc(:,b),Ni=Ni(b),Nj=Nj(b),Nk=Nk(b)) ! allocating block variables
  ! loading block nodes coordinates
  read(pp%unit_grd)(((node(i,j,k)%x,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(pp%unit_grd)(((node(i,j,k)%y,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  read(pp%unit_grd)(((node(i,j,k)%z,i=0-gc(1,b),Ni(b)+gc(2,b)),j=0-gc(3,b),Nj(b)+gc(4,b)),k=0-gc(5,b),Nk(b)+gc(6,b))
  if (pp%fcc) then ! loading icc variables
    read(pp%unit_icc)(((icc(i,j,k),i=1-gc(1,b),Ni(b)+gc(2,b)),j=1-gc(3,b),Nj(b)+gc(4,b)),k=1-gc(5,b),Nk(b)+gc(6,b))
    if (pp%verbose) then
      Ncc(b) = count(icc(1:Ni(b),1:Nj(b),1:Nk(b))>0)
      write(stdout,'(A)')'+--> Block:'//trim(str(.true.,b))//' Number of internal chimera cells: '//trim(str(.true.,Ncc(b)))
    endif
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
    if (.not.pp%cell) then ! ricc must be interpolated at nodes
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
  if (pp%sol) then ! allocating and loading solution variables
    read(pp%unit_sol)(((momentum(i,j,k)%x,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(pp%unit_sol)(((momentum(i,j,k)%y,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(pp%unit_sol)(((momentum(i,j,k)%z,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    read(pp%unit_sol)(((pressure(i,j,k)  ,i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    if (pp%level_set) then
      read(pp%unit_sol)(((f (i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(pp%unit_sol)(((f0(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (pp%zeroeq.or.pp%oneeq.or.pp%twoeq) then
      read(pp%unit_sol)(((visc(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (pp%oneeq) then
      read(pp%unit_sol)(((vitl(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (pp%twoeq) then
      read(pp%unit_sol)(((ken(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
      read(pp%unit_sol)(((eps(i,j,k),i=0,Ni(b)+1),j=0,Nj(b)+1),k=0,Nk(b)+1)
    endif
    if (pp%vordet) then ! computing vord variable
      ! computing the metrics of cells
      call compute_metrics(gc=gc(:,b),Ni=Ni(b),Nj=Nj(b),Nk=Nk(b))
      ! correcting the metrics of boundary conditions cells
      call bc_metrics_correction(Ni=Ni(b),Nj=Nj(b),Nk=Nk(b),rcc=rcc)
      ! computing vordet variable
      call compute_vorticity(gc=gc(:,b),Ni=Ni(b),Nj=Nj(b),Nk=Nk(b))
    endif
  endif
  if (pp%tec) error = file_tec%save_block(pp=pp,b=b,Ni=Ni(b),Nj=Nj(b),Nk=Nk(b))
enddo

call pp%finalize()
if (pp%tec) call file_tec%finalize(pp=pp)
if (allocated(gc )) deallocate(gc )
if (allocated(Ni )) deallocate(Ni )
if (allocated(Nj )) deallocate(Nj )
if (allocated(Nk )) deallocate(Nk )
if (allocated(rcc)) deallocate(rcc)
if (allocated(Ncc)) deallocate(Ncc)
!if (vtk) then
!  ! saving the multiblock VTM wrapper
!  err = VTM_INI_XML(filename=adjustl(trim(File_out))//'.vtm')
!  err = VTM_BLK_XML(block_action='open')
!  err = VTM_WRF_XML(flist=(/(adjustl(trim(OS%basename(File_out)))//'-vts_files'//OS%sep//        &
!                             adjustl(trim(OS%basename(File_out)))//"."//trim(strz(4,b))//'.vts', &
!                             b=1,Nb)/))
!  err = VTM_BLK_XML(block_action='close')
!  err = VTM_END_XML()
!endif
stop
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  !> @brief Subroutine for parsing command line arguments.
  subroutine parse_command_arguments(pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  type(Type_PostProcess), intent(INOUT):: pp    !< Post-processor data.
  type(Type_Command_Line_Interface)::     cli   !< Command Line Interface (CLI).
  integer(I4P)::                          error !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  write(stdout,'(A)')'+--> XnPlot, post-processor for Xnavis '
  ! initializing CLI
  call cli%init(progname='XnPlot',                                                            &
                version ='v0.0.2',                                                            &
                help    =' The XnPlot Command Line Interface (CLI) has the following options',&
                examples=["XnPlot -g xship.grd -o grid",                                      &
                          "XnPlot -g cc.01.grd -i cc.01 -o mesh.01",                          &
                          "XnPlot -g cc.01.grd -i cc.01 -s sol.00.01 -o sol.01"])
  ! setting CLAs
  call cli%add(pref='|-->',switch='-g',     help='Grid file (.grd)',required=.true.,act='store',error=error)
  call cli%add(pref='|-->',switch='-o',     help='output file name; default is basename of grd file with the proper extension',&
                                            required=.false.,act='store',def='unset',error=error)
  call cli%add(pref='|-->',switch='-i',     help='ICC file',required=.false.,act='store',def='unset',error=error)
  call cli%add(pref='|-->',switch='-s',     help='solution file name; if passed the solution variables are saved',&
                                            required=.false.,act='store',def='unset',error=error)
  call cli%add(pref='|-->',switch='-ngc',   help='mesh without ghosts cells, as geogrd output',&
                                            required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-cell',  help='variables other than nodes coord. are cell centered',&
                                            required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-ls',    help='solution with level set',&
                                            required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-nt',    help='no turbulent model used',&
                                            required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-eq',    help='equations turbulent model',required=.false.,act='store',&
                                            def='1',choices='0,1,2',error=error)
  call cli%add(pref='|-->',switch='-vordet',help='computing variables for vortices identification',&
                                            required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-ascii', help='write ascii output files',&
                                            required=.false.,act='store_true',def='.false.',error=error)
  call cli%add(pref='|-->',switch='-tec',   help='write output Tecplot files',&
                                            required=.false.,act='store',def='yes',choices='yes,no',error=error)
  call cli%add(pref='|-->',switch='-vtk',   help='write output VTK files',&
                                            required=.false.,act='store',def='no',choices='yes,no',error=error)
  call cli%add(pref='|-->',switch='-proc',  help='process number for global block numeration if proc.input is found',&
                                            required=.false.,act='store',def='-1',error=error)
  call cli%add(pref='|-->',switch='-os',    help='type of Operating System',&
                                            required=.false.,act='store',def='UIX',choices='UIX,WIN',error=error)
  call cli%add(pref='|-->',switch='-vb',    help='Verbose output',required=.false.,act='store_true',def='.false.',error=error)
  ! parsing CLI
  write(stdout,'(A)')'+--> Parsing Command Line Arguments'
  call cli%parse(error=error,pref='|-->'); if (error/=0) stop
  ! setting post-processing options accordingly to CLI
  call pp%set_from_cli(cli=cli)
  call pp%set_from_mbpar()
  call pp%set_from_procinput()
  call pp%input_files_init()
  if (pp%tec) call file_tec%init(pp=pp)
  !if (vtk) then
  !  vtkbdir  = adjustl(trim(OS%basedir(File_out)))//adjustl(trim(OS%basename(File_out)))//'-vts_files'//OS%sep
  !  vtkbfile = adjustl(trim(vtkbdir))//adjustl(trim(OS%basename(File_out)))
  !endif
              write(stdout,'(A)')'+--> Files'
              write(stdout,'(A)')'|--> GRD File '//trim(pp%File_grd)
  if (pp%fcc) write(stdout,'(A)')'|--> ICC File '//trim(pp%File_icc)
  if (pp%sol) write(stdout,'(A)')'|--> SOL File '//trim(pp%File_sol)
              write(stdout,'(A)')'|--> OUT File '//trim(pp%File_out)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine parse_command_arguments
endprogram XnPlot
