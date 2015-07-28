module Lib_VTK
!> @brief This module contains procedures for saving VTK output.
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision          ! Integers and reals precision definition.
USE Block_Variables       ! Block variables definition.
USE Data_Type_PostProcess ! Definition of Type_PostProcess.
USE Data_Type_OS          ! Definition of Type_OS.
USE Data_Type_Vector      ! Definition of Type_Vector.
USE Lib_IO_Misc           ! Library for miscellanea IO procedures.
USE Lib_VTK_IO            ! Library for I/O VTK files.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public:: file_vtk
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type:: Type_File_VTK
  type(Type_OS)::                 OS      !< OS definition.
  character(len=:), allocatable:: rootdir !< Root of directory.
  contains
    procedure:: init       ! Procedure for initializing file writing.
    procedure:: save_block ! Procedure for saving one block.
    procedure:: finalize   ! Procedure for finalizing file writing.
endtype Type_File_VTK
type(Type_File_VTK):: file_vtk
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  subroutine init(file_d,pp)
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @brief Procedure for initializing file writing.
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_VTK),   intent(INOUT):: file_d !< File data.
  type(Type_PostProcess), intent(IN)::    pp     !< Post-processor data.
  integer(I4P)::                          error  !< Error trapping flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  file_d%rootdir = trim(adjustl(pp%File_out))//'-vtk'//file_d%OS%sep
  error = file_d%OS%make_dir(directory=file_d%rootdir)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine init

  function save_block(file_d,pp,b,Ni,Nj,Nk) result(error)
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @brief Procedure for saving block.
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_VTK),   intent(INOUT):: file_d      !< File data.
  type(Type_PostProcess), intent(IN)::    pp          !< Post-processor data.
  integer(I4P),           intent(IN)::    b           !< Actual block number.
  integer(I4P),           intent(IN)::    Ni,Nj,Nk    !< Number of cells.
  integer(I4P)::                          error       !< Error trapping flag.
  integer(I4P)::                          nnode,ncell !< Number of nodes and cells.
  real(R8P), allocatable::                vari(:,:,:) !< Interpolated generic variables.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.pp%cell) then ! variables must be interpolated at nodes, allocating dummy variable
    if (allocated(vari)) deallocate(vari) ; allocate(vari(0:Ni,0:Nj,0:Nk))
  endif
  nnode = (Ni+1)*(Nj+1)*(Nk+1) ! computing number of nodes
  ncell =  Ni   * Nj   * Nk    ! computing number of cells
  associate(binary=>pp%binary,cell=>pp%cell,blockmap=>pp%blockmap,sol=>pp%sol,fcc=>pp%fcc,vordet=>pp%vordet,&
            zeroeq=>pp%zeroeq,oneeq=>pp%oneeq,twoeq=>pp%twoeq,level_set=>pp%level_set)
  if (binary) then
    error = VTK_INI_XML(output_format='binary',                                                                                    &
                        filename=file_d%rootdir//'blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//'.vts',&
                        mesh_topology='StructuredGrid',                                                                            &
                        nx1=0,nx2=Ni,ny1=0,ny2=Nj,nz1=0,nz2=Nk)
  else
    error = VTK_INI_XML(output_format='ascii',                                                                                     &
                        filename=file_d%rootdir//'blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//'.vts',&
                        mesh_topology='StructuredGrid',                                                                            &
                        nx1=0,nx2=Ni,ny1=0,ny2=Nj,nz1=0,nz2=Nk)
  endif
  error = VTK_GEO_XML(nx1=0,nx2=Ni,ny1=0,ny2=Nj,nz1=0,nz2=Nk,NN=nnode,&
                     X=node(0:Ni,0:Nj,0:Nk)%x,Y=node(0:Ni,0:Nj,0:Nk)%y,Z=node(0:Ni,0:Nj,0:Nk)%z)
  if (cell) then
    error = VTK_DAT_XML(var_location='cell',var_block_action='open')
    if (fcc) then
      error = VTK_VAR_XML(NC_NN=ncell,varname='icc',var=real(ricc(1:Ni,1:Nj,1:Nk),R8P))
    endif
    if (sol) then
      error = VTK_VAR_XML(NC_NN=ncell,varname='u',var=momentum(1:Ni,1:Nj,1:Nk)%x)
      error = VTK_VAR_XML(NC_NN=ncell,varname='v',var=momentum(1:Ni,1:Nj,1:Nk)%y)
      error = VTK_VAR_XML(NC_NN=ncell,varname='w',var=momentum(1:Ni,1:Nj,1:Nk)%z)
      error = VTK_VAR_XML(NC_NN=ncell,varname='p',var=pressure(1:Ni,1:Nj,1:Nk)  )
      if (zeroeq) then
        error = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=visc(1:Ni,1:Nj,1:Nk))
      elseif (oneeq) then
        error = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=visc(1:Ni,1:Nj,1:Nk))
        error = VTK_VAR_XML(NC_NN=ncell,varname='vitl',var=vitl(1:Ni,1:Nj,1:Nk))
      elseif (twoeq) then
        error = VTK_VAR_XML(NC_NN=ncell,varname='visc',var=visc(1:Ni,1:Nj,1:Nk))
        error = VTK_VAR_XML(NC_NN=ncell,varname='ken', var=ken (1:Ni,1:Nj,1:Nk))
        error = VTK_VAR_XML(NC_NN=ncell,varname='eps', var=eps (1:Ni,1:Nj,1:Nk))
      endif
      if (level_set) then
        error = VTK_VAR_XML(NC_NN=ncell,varname='f', var=f (1:Ni,1:Nj,1:Nk))
        error = VTK_VAR_XML(NC_NN=ncell,varname='f0',var=f0(1:Ni,1:Nj,1:Nk))
      endif
      if (vordet) then
        error = VTK_VAR_XML(NC_NN=ncell,varname='vord', var=vord(     1:Ni,1:Nj,1:Nk))
        error = VTK_VAR_XML(NC_NN=ncell,varname='qfac', var=qfac(     1:Ni,1:Nj,1:Nk))
        error = VTK_VAR_XML(NC_NN=ncell,varname='heli', var=heli(     1:Ni,1:Nj,1:Nk))
        error = VTK_VAR_XML(NC_NN=ncell,varname='vorX', var=vorticity(1:Ni,1:Nj,1:Nk)%x)
        error = VTK_VAR_XML(NC_NN=ncell,varname='vorY', var=vorticity(1:Ni,1:Nj,1:Nk)%y)
        error = VTK_VAR_XML(NC_NN=ncell,varname='vorZ', var=vorticity(1:Ni,1:Nj,1:Nk)%z)
      endif
    endif
    error = VTK_DAT_XML(var_location='cell',var_block_action='close')
  else
    error = VTK_DAT_XML(var_location='node',var_block_action='open')
    if (fcc) then
      error = VTK_VAR_XML(NC_NN=nnode,varname='icc',var=real(vicc(0:Ni,0:Nj,0:Nk),R8P))
    endif
    if (sol) then
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=momentum(0:Ni+1,0:Nj+1,0:Nk+1)%x,vari=vari)
      error = VTK_VAR_XML(NC_NN=nnode,varname='u', var=vari(0:Ni,0:Nj,0:Nk))
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=momentum(0:Ni+1,0:Nj+1,0:Nk+1)%y,vari=vari)
      error = VTK_VAR_XML(NC_NN=nnode,varname='v', var=vari(0:Ni,0:Nj,0:Nk))
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=momentum(0:Ni+1,0:Nj+1,0:Nk+1)%z,vari=vari)
      error = VTK_VAR_XML(NC_NN=nnode,varname='w', var=vari(0:Ni,0:Nj,0:Nk))
      call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=pressure(0:Ni+1,0:Nj+1,0:Nk+1)  ,vari=vari)
      error = VTK_VAR_XML(NC_NN=nnode,varname='p', var=vari(0:Ni,0:Nj,0:Nk))
      if (zeroeq) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=visc(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='visc', var=vari(0:Ni,0:Nj,0:Nk))
      elseif (oneeq) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=visc(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='visc', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vitl(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='vitl', var=vari(0:Ni,0:Nj,0:Nk))
      elseif (twoeq) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=visc(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='visc', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=ken (0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='ken', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=eps (0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='eps', var=vari(0:Ni,0:Nj,0:Nk))
      endif
      if (level_set) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=f (0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='f', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=f0(0:Ni+1,0:Nj+1,0:Nk+1),vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='f0', var=vari(0:Ni,0:Nj,0:Nk))
      endif
      if (vordet) then
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vord(     0:Ni+1,0:Nj+1,0:Nk+1),  vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='vord', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=qfac(     0:Ni+1,0:Nj+1,0:Nk+1),  vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='qfac', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=heli(     0:Ni+1,0:Nj+1,0:Nk+1),  vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='heli', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vorticity(0:Ni+1,0:Nj+1,0:Nk+1)%x,vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='vorX', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vorticity(0:Ni+1,0:Nj+1,0:Nk+1)%y,vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='vorY', var=vari(0:Ni,0:Nj,0:Nk))
        call varinterpolation(Ni=Ni,Nj=Nj,Nk=Nk,var=vorticity(0:Ni+1,0:Nj+1,0:Nk+1)%z,vari=vari)
        error = VTK_VAR_XML(NC_NN=nnode,varname='vorZ', var=vari(0:Ni,0:Nj,0:Nk))
      endif
    endif
    error = VTK_DAT_XML(var_location='node',var_block_action='close')
  endif
  error = VTK_GEO_XML()
  error = VTK_END_XML()
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction save_block

  subroutine finalize(file_d,pp,Nb)
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @brief Procedure for finalizing file writing.
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_File_VTK),   intent(IN):: file_d  !< File data.
  type(Type_PostProcess), intent(IN):: pp      !< Post-processor data.
  integer(I4P),           intent(IN):: Nb      !< Number of blocks.
  integer(I4P)::                       b       !< Blocks counter.
  integer(I4P)::                       error   !< Error trapping flag.
  character(len=:), allocatable::      rootdir !< Root dir.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rootdir = file_d%rootdir
  rootdir = trim(file_d%OS%basename(file_d%rootdir(1:len(file_d%rootdir)-1)))//file_d%OS%sep
  ! saving the multiblock VTM wrapper
  error = VTM_INI_XML(filename=adjustl(trim(pp%File_out))//'.vtm')
  error = VTM_BLK_XML(block_action='open')
  associate(blockmap=>pp%blockmap)
  error = VTM_WRF_XML(flist=[(rootdir//'blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b)))//'.vts',b=1,Nb)],&
                      nlist=[('blk_'//trim(strz(4,blockmap(2,b)))//'_grp_'//trim(strz(3,blockmap(1,b))),b=1,Nb)])
  endassociate
  error = VTM_BLK_XML(block_action='close')
  error = VTM_END_XML()
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine finalize
endmodule Lib_VTK
