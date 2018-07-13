!==============================================================================|
subroutine read_mesh(domain)
  ! Read bathymetry and grid coefficients from text file

  use MOD_SIM, only : LAG_SIM
  use ALL_VARS
  use LIMS
  
  implicit none
  
  class(LAG_SIM), intent(inout) :: domain
  integer :: ii, index, ionode=100, ioelem=101

  write(*, "(A)", advance='no') "Opening mesh files... "
  open(ioelem, file = "./"//trim(folderprefix)//"/mesh_elem.dat")
  open(ionode, file = "./"//trim(folderprefix)//"/mesh_node.dat")
  write(*, *) "Finished"
  
  
  write(*, "(A)", advance='no') "Reading headers... "
  read(ioelem, *) N, KB ! get number of nodes, elements, and sigma layers
  read(ionode, *) M
  write(*, *) "Finished"
  
  
  domain%nnodes = M; domain%nelements = N; domain%nlayers = KB
  KBM1 = KB - 1
  KBM2 = KB - 2
  
  write(*, "(A)", advance='no') "Allocating mesh based variables... "
  call ALLOC_VARS
  write(*, *) "Finished"
  
  write(*, "(A)", advance='no') "Getting element data... "
  element_loop: do ii = 1, N 
     read(ioelem, *) index, NV(ii, 1), NV(ii, 3), NV(ii, 2) ! get node indices
  end do element_loop
  NV(:, 4) = NV(:, 1) ! duplicate node for computation
  write(*, *) "Finished"

  write(*, "(A)", advance='no') "Getting node data... "
  node_loop: do ii = 1, M
     read(ionode, *) index, VX(ii), VY(ii), H(ii) ! get node position and depth
  end do node_loop
  write(*, *) "Finished"

  close(ioelem)
  close(ionode)

end subroutine read_mesh
