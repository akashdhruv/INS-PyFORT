subroutine tempsolver(T_new,T,u,v,dx,dy,dt,inRe,Pr,Nxb,Nyb)

  implicit none

  integer,intent(in) :: Nxb,Nyb
  real, intent(in), dimension (Nxb+2,Nyb+2) :: T,u,v
  real, intent(in) :: dx,dy,dt,inRe,Pr
  real, intent(out), dimension(Nxb+2,Nyb+2) :: T_new

  real :: u_plus, u_mins, v_plus, v_mins, u_conv, v_conv
  real :: Tx_plus, Tx_mins, Ty_plus, Ty_mins
  integer :: i,j

  T_new = T

  do j=2,Nyb+1
     do i=2,Nxb+1

     u_conv = (u(i,j)+u(i-1,j))/2.
     v_conv = (v(i,j)+v(i,j-1))/2.

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     Tx_plus = (T(i+1,j)-T(i,j))/dx
     Tx_mins = (T(i,j)-T(i-1,j))/dx

     Ty_plus = (T(i,j+1)-T(i,j))/dy
     Ty_mins = (T(i,j)-T(i,j-1))/dy

     T_new(i,j) = T(i,j)+((dt*inRe)/(Pr*dx*dx))*(T(i+1,j)+T(i-1,j)-2*T(i,j))&
                        +((dt*inRe)/(Pr*dy*dy))*(T(i,j+1)+T(i,j-1)-2*T(i,j))&
                        -((dt))*(u_plus*Tx_mins + u_mins*Tx_plus)&
                        -((dt))*(v_plus*Ty_mins + v_mins*Ty_plus)

     end do
  end do


end subroutine tempsolver
