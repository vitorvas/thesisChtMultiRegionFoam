! ALVARO BERNAL GARCIA
! Phd Candidate
! ISIRYM, UNIVERSIDAD POLITECNICA DE VALENCIA
! VALENCIA, SPAIN
! 24/09/2013



module variables_xs_module

implicit none

public ini_variables_pbttinput



! Nodes and their material composition that will be read from the file "GEOM"
integer number_nodes_pz, number_zplane
integer, dimension (:,:), allocatable :: materials_nodes ! (number_zplane,number_nodes_pz)
double precision fully_inserted_z, step_size
integer, dimension (:), allocatable :: control_rod_map ! (number_nodes_pz)
double precision, dimension (:), allocatable :: axial_mesh_nodes ! (number_zplane+1)


! Variables defined in file "pbttinput.f90"
integer ntfentry,ndmentry,ntentry,ntdentry,ibor
double precision ppml1,ppml2,ppml3
logical PR_exist


! Variables used to store cross section data   -->  pbtt.h
integer ncmpur13,ncmpr13
double precision, dimension(:,:,:,:), allocatable :: xspbttur,&
& xspbttur2, xspbttur3   ! (ntdentry,6,2,ncmpur13)
                    ! 6: number of cross sections , 2: number of groups
double precision, dimension(:,:,:,:), allocatable :: xspbttr, &
& xspbttr2, xspbttr3    ! (ntdentry,6,2,ncmpr13)
double precision, allocatable, dimension(:,:) :: rvelpbur ! (2,ncmpur13)
double precision, allocatable, dimension(:,:) :: rvelpbr ! (2,ncmpr13)
double precision, allocatable, dimension(:,:,:,:) :: adfpbur, detpbur ! (ntdentry,2,2,ncmpur13)
double precision, allocatable, dimension(:,:,:,:) :: adfpbr, detpbr ! (ntdentry,2,2,ncmpr13)
double precision, allocatable, dimension(:,:) :: betapbur, lambdapbur ! (6,ncmpur13)
                                                                       ! 6: number of groups of precursors
double precision, allocatable, dimension(:,:) :: betapbr,lambdapbr ! (6,ncmpr13)



! Variables used to store cross section data   -->  xsec.h
integer, allocatable, dimension(:) :: icrcomp ! (ncmpur13)



! Variables defined in file "xsecfb13.f90"
double precision :: timeedt, eps
double precision, dimension (:,:), allocatable :: ppml
logical :: rotadf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  abernal, 18/10/14
integer thcoupling ! thermal-hydraulic coupling. If thcoupling=0 no couling, if thcoupling=1 coupling
double precision, allocatable, dimension (:) :: temp_coupling,dens_coupling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! Variables defines in file "pbttxsec.f90"
! THESE VARIABLES ARE THE FINAL CROSS SECTIONS
double precision, dimension (:,:,:), allocatable :: xsd, xst, xsa, xstr, &
& xsnf, xskf, xsf, xsxea, adfn, adfe, adfs, adfw ! (2,number_nodes_pz,number_zplane)
                                                                 ! 2: number of energy groups
double precision, dimension (:,:), allocatable :: xss, rnxe ! (number_nodes_pz,number_zplane)



! Variable used to calculate ADF for the different materials
double precision, dimension (:,:,:), allocatable :: mat_adf




contains




! The following subroutine is used to initizalize the variables used in NEMTAB XS processing
subroutine ini_variables_pbttinput


allocate (xspbttur(ntdentry,6,2,ncmpur13),xspbttur2(ntdentry,6,2,ncmpur13),&
& xspbttur3(ntdentry,6,2,ncmpur13))
allocate (xspbttr(ntdentry,6,2,ncmpr13),xspbttr2(ntdentry,6,2,ncmpr13),&
& xspbttr3(ntdentry,6,2,ncmpr13))
allocate (rvelpbur(2,ncmpur13))
allocate (rvelpbr(2,ncmpr13))
allocate (adfpbur(ntdentry,2,2,ncmpur13), detpbur(ntdentry,2,2,ncmpur13))
allocate (adfpbr(ntdentry,2,2,ncmpr13), detpbr(ntdentry,2,2,ncmpr13))
allocate (betapbur(6,ncmpur13), lambdapbur(6,ncmpur13))
allocate (betapbr(6,ncmpr13),lambdapbr(6,ncmpr13))
allocate (icrcomp(ncmpur13))


end subroutine ini_variables_pbttinput








! The following subroutine is used to initizalize the variables used to store cross sections
subroutine ini_variables_pbttxsec


allocate (xsd(2,number_nodes_pz,number_zplane), xst(2,number_nodes_pz,number_zplane),&
& xsa(2,number_nodes_pz,number_zplane), xstr(2,number_nodes_pz,number_zplane),&
& xsnf(2,number_nodes_pz,number_zplane), xskf(2,number_nodes_pz,number_zplane),&
& xsf(2,number_nodes_pz,number_zplane), xsxea(2,number_nodes_pz,number_zplane))
allocate (xss(number_nodes_pz,number_zplane),rnxe(number_nodes_pz,number_zplane))
allocate (adfn(2,number_nodes_pz,number_zplane), adfe(2,number_nodes_pz,number_zplane),&
& adfs(2,number_nodes_pz,number_zplane), adfw(2,number_nodes_pz,number_zplane))


end subroutine ini_variables_pbttxsec






end module variables_xs_module
