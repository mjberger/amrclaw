! ============================================================================
!  Program:     AMRClaw
!  File:        resize_nodes.f90
!  Author:      Mandli and Marsha Berger
! ============================================================================
!  Description:  Resize the node space if have run out 
! ============================================================================

subroutine resize_nodes(new_size,status)
    
    use amr_module
    implicit none
    
    integer, intent(out) :: status
    integer, intent(in) :: new_size

    integer :: i
    
    real(kind=8), allocatable, target, dimension(:,:) :: new_rnode
    integer, allocatable, target, dimension(:,:) :: new_node
    
    print *, "Expanding maximum number of grids from ", maxgr," to ", new_size

    ! first for rnode
    allocate(new_rnode(rsize,new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_rnode(1:rsize,1:maxgr) = rnode     ! new way, use allocatable, not pointer       

    call move_alloc(new_rnode,rnode)

    ! next for node
    allocate(new_node(nsize,new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_node(1:nsize,1:maxgr) = node     ! new way, use allocatable, not pointer       

    call move_alloc(new_node,node)

    !! need to rethread new space to be able to use it when new grids requested
    do i = maxgr, new_size
       node(nextfree,i) = i+1
    end do
    ! reset last one to null
    node(nextfree, new_size) = null


    ! reset maxgr and next free node,  to continue
    maxgr = new_size
    ndfree = maxgr + 1     
    
    return
    
end subroutine resize_nodes
