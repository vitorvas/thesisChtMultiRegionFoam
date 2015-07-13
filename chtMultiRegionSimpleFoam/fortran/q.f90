subroutine heatflux(q)
  
! ---------------------------------------------------------------
! 5760: number of cells in the fuel region of the blockMesh mesh
! ---------------------------------------------------------------

      real, dimension(9)  :: q
!Se fosse uma funcao      real, dimension(9) :: heatflux

      do i = 1, 9
! -----------------------------------------
! 1e6 is the energy related to h (kg/m s^3)
! and can be considered Watt/volume
! ----------------------------------------- 
         q(i) = i
      end do
end subroutine

