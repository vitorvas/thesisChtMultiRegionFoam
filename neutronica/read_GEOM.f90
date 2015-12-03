! ALVARO BERNAL GARCIA
! PhD Candidate
! ISIRYM, UNIVERSITAT POLITECNICA DE VALENCIA
! 8/10/2013



! This subroutine reads the file GEOM containing the material compositions for each node


 subroutine read_GEOM

  use variables_xs_module

  integer state_file, aux, nx, ny, nz, z, i, j, ini, fin, error, finx, mult
  integer, allocatable :: vector(:)
  character*133 line, comp_line
  double precision, allocatable :: axial_length_nodes(:)


  ! The file "GEOM" is opened
  open (unit=8,file='GEOM',status='OLD',iostat=state_file)

  ! Checking errors at the opening of the file
  if (state_file/=0) then
     print '(////)'
     print *,'ERROR AT THE OPENING OF THE FILE "GEOM"'
     print '(////)'
     stop
  end if

  ! Look for the line starting by "geo_dim", which contains the mesh: nx*ny*nz
  aux=0
  comp_line=' geo_dim'
  do while (aux==0)
     read(8,'(A133)') line
     if (line(1:8)==comp_line) then
        aux=1
        read(line(9:133),*)  nx, ny, nz
        number_zplane=nz ! "number_zplane" is defined in "variables_xs_module"
     end if
  end do


  ! Go to the line starting by "grid_z"
  aux=0
  comp_line=' grid_z'
  do while (aux==0)
     read(8,'(A133)') line
     if (line(1:7)==comp_line) then
        aux=1
     end if
  end do

  ! Get the axial lengths of the nodes
  allocate (axial_length_nodes(number_zplane))
  fin=len_trim(line)
  line=adjustl(line(8:fin))
  fin=0
  i=1
  aux=0
  do while (aux==0)
     comp_line=' '
     ini=fin+1
     fin=scan(line(ini:len_trim(line)),comp_line)+ini-1
     if (fin<ini) then
        fin=len_trim(line)
        aux=1
     end if
     finx=scan(line(ini:fin),'*')+ini-1
     if (finx<ini) then
        read(line(ini:fin),*,iostat=error) axial_length_nodes(i)
     else
        read(line(ini:finx-1),*) mult
        read(line(finx+1:fin),*,iostat=error) axial_length_nodes(i)
        axial_length_nodes(i:i+mult-1)=axial_length_nodes(i)
        i=i+mult-1
     end if
     if (error==0) i=i+1
  end do

  ! Calculate the axial mesh of the nodes
  allocate(axial_mesh_nodes(number_zplane+1))
  axial_mesh_nodes(1)=0
  do i=2,number_zplane+1
     axial_mesh_nodes(i)=axial_mesh_nodes(i-1)+axial_length_nodes(i-1)
  end do
  deallocate(axial_length_nodes)


  ! Go to the first line starting by " planar_reg"
  aux=0
  comp_line=' planar_reg'
  do while (aux==0)
     read(8,'(A133)') line
     if (line(1:11)==comp_line) then
        aux=1
     end if
  end do

  ! Determine the number of nodes per z-plane, without knowing previously the number of nodes in each row
  allocate (vector(nx))
  number_nodes_pz=0 ! "number_nodes_pz" is defined in "variables_xs_module"
  do j=1,ny
     vector(:)=0
     read(8,'(A133)') line
     line=trim(adjustl(line))
     comp_line=' '
     fin=0
     i=1
     aux=0
     do while (aux==0)
        ini=fin+1
        fin=scan(line(ini:len_trim(line)),comp_line)+ini-1
        if (fin<ini) then
           fin=len_trim(line)
           aux=1
        end if
        read(line(ini:fin),*,iostat=error) vector(i)
        if (error==0) i=i+1
     end do
     number_nodes_pz=number_nodes_pz+i-1
  end do
  deallocate(vector)



  ! Instructions to get the material composition for each z-plane
  ! Firstly, the file must be closed and opened another time
  close(8)
  open (unit=8,file='GEOM',status='OLD')
  ! Then, go to the first line starting by " planar_reg"
  aux=0
  comp_line=' planar_reg'
  do while (aux==0)
     read(8,'(A133)') line
     if (line(1:11)==comp_line) then
        aux=1
     end if
  end do
  ! Finally, get the material composition for each z-plane
  allocate (materials_nodes(number_zplane,number_nodes_pz))  ! "materials_nodes" is defined in "variables_xs_module"
  do z=1,nz
     read(8,*) (materials_nodes(z,i),i=1,number_nodes_pz)
     read(8,*)
  end do


  ! Go to the line where control rods data is defined : CR_axinfo
  aux=0
  comp_line=' CR_axinfo'
  do while (aux==0)
     read(8,'(A133)') line
     if (line(1:10)==comp_line) then
        aux=1
     end if
  end do

  ! Get fully inserted position and step size of control rods
  fin=len_trim(line)
  read(line(11:fin),*) fully_inserted_z,step_size, comp_line


  ! Go to the line where control rods are defined :  bank_conf
  aux=0
  comp_line=' bank_conf'
  do while (aux==0)
     read(8,'(A133)') line
     if (line(1:10)==comp_line) then
        aux=1
     end if
  end do

  ! Get the control rod map
  allocate (control_rod_map(number_nodes_pz))
  read(8,*) (control_rod_map(i),i=1,number_nodes_pz)

  ! Close the file "GEOM"
  close(8)


  return
 end subroutine read_GEOM
