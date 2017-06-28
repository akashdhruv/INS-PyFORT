subroutine solver(ps,ps_rhs,ps_new,dx,dy,Nxb,Nyb)

   implicit none
   integer, intent(in) :: Nxb,Nyb
   real, intent(in) :: dx,dy
   real, intent(in), dimension(Nxb,Nyb) :: ps_rhs
   real, intent(in), dimension(Nxb+2,Nyb+2) :: ps
   real, intent(out), dimension(Nxb+2,Nyb+2) :: ps_new

   integer :: i,j

   ps_new = ps

   do j=2,Nyb+1
      do i=2,Nxb+1

      ps_new(i,j)=((ps(i,j+1)/(dy*dy))+(ps_new(i,j-1)/(dy*dy))&
                  +(ps(i+1,j)/(dx*dx))+(ps_new(i-1,j)/(dx*dx))&
                  +ps_RHS(i-1,j-1))&
                  *(1/((1/(dx*dx))+(1/(dy*dy))+&
                  (1/(dx*dx))+(1/(dy*dy))))

      end do
   end do

end subroutine solver
