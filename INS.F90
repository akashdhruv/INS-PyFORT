subroutine predictor(u,v,G1_old,G2_old,G1_new,G2_new,inRe,dx,dy,dt,u_star,v_star,tstep,Nxb,Nyb)


    implicit none
    integer, intent(in) :: Nxb,Nyb,tstep
    real,dimension(Nxb+2,Nyb+2),intent(in) :: u,v
    real,intent(in) :: dx,dy,inRe,dt
    real,dimension(Nxb+2,Nyb+2),intent(out) :: u_star,v_star
    real,dimension(Nxb,Nyb),intent(in) :: G1_old,G2_old   
    real,dimension(Nxb,Nyb),intent(out):: G1_new,G2_new


    real,dimension(Nxb,Nyb) :: C1,D1,C2,D2,G1,G2

    call Convective_U(u,v,dx,dy,C1,Nxb,Nyb)
    call Diffusive_U(u,dx,dy,inRe,D1,Nxb,Nyb)

    G1 = C1 + D1

    if (tstep == 0) then

         u_star(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dt/1)*(G1)
         G1_new = G1
    else

         u_star(2:Nxb+1,2:Nyb+1)=u(2:Nxb+1,2:Nyb+1)+(dt/2)*(3*G1_old-G1)
         G1_new = G1
    endif

    call Convective_V(u,v,dx,dy,C2,Nxb,Nyb)
    call Diffusive_V(v,dx,dy,inRe,D2,Nxb,Nyb)

    G2 = C2 + D2

    if (tstep == 0) then

         v_star(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dt/1)*(G2)
         G2_new = G2
    else

         v_star(2:Nxb+1,2:Nyb+1)=v(2:Nxb+1,2:Nyb+1)+(dt/2)*(3*G2_old-G2)
         G2_new = G2
    endif

end subroutine predictor

subroutine corrector(ut,vt,p,u,v,dt,dx,dy,Nxb,Nyb)

    implicit none

    integer,intent(in) :: Nxb,Nyb
    real, intent(in) :: dt,dx,dy
    real, intent(in), dimension(Nxb+2,Nyb+2) :: ut,vt,p

    real, intent(out), dimension(Nxb+2,Nyb+2) :: u,v

    u(2:Nxb+1,2:Nyb+1) = ut(2:Nxb+1,2:Nyb+1) - (dt/dx)*(p(3:Nxb+2,2:Nyb+1)-p(2:Nxb+1,2:Nyb+1))
    v(2:Nxb+1,2:Nyb+1) = vt(2:Nxb+1,2:Nyb+1) - (dt/dy)*(p(2:Nxb+1,3:Nyb+2)-p(2:Nxb+1,2:Nyb+1))
   
end subroutine corrector


subroutine Convective_U(ut,vt,dx,dy,C1,Nxb,Nyb)

      implicit none

      integer,intent(in) :: Nxb,Nyb
      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt
      real,intent(in) :: dx
      real,intent(in) :: dy
      real, dimension(Nxb,Nyb), intent(out) :: C1

      real, dimension(Nxb,Nyb) :: ue,uw,us,un,vs,vn

      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(3:Nxb+2,2:Nyb+1))/2
      uw = (ut(2:Nxb+1,2:Nyb+1)+ut(1:Nxb,2:Nyb+1))/2
      us = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,1:Nyb))/2
      un = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      vs = (vt(2:Nxb+1,1:Nyb)+vt(3:Nxb+2,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2

      C1 = -((ue**2)-(uw**2))/dx - ((un*vn)-(us*vs))/dy

end subroutine Convective_U

subroutine Convective_V(ut,vt,dx,dy,C2,Nxb,Nyb)

      implicit none

      integer, intent(in) :: Nxb,Nyb
      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt
      real, intent(in) :: dx
      real, intent(in) :: dy
      real,dimension(Nxb,Nyb) :: vn,vs,ve,vw,ue,uw
      real, dimension(Nxb,Nyb), intent(out) :: C2

      vs = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,1:Nyb))/2
      vn = (vt(2:Nxb+1,2:Nyb+1)+vt(2:Nxb+1,3:Nyb+2))/2
      ve = (vt(2:Nxb+1,2:Nyb+1)+vt(3:Nxb+2,2:Nyb+1))/2
      vw = (vt(2:Nxb+1,2:Nyb+1)+vt(1:Nxb,2:Nyb+1))/2
      ue = (ut(2:Nxb+1,2:Nyb+1)+ut(2:Nxb+1,3:Nyb+2))/2
      uw = (ut(1:Nxb,2:Nyb+1)+ut(1:Nxb,3:Nyb+2))/2

      C2 = -((ue*ve)-(uw*vw))/dx - ((vn**2)-(vs**2))/dy

end subroutine Convective_V

subroutine Diffusive_U(ut,dx,dy,inRe,D1,Nxb,Nyb)

      implicit none

      integer, intent(in) :: Nxb,Nyb
      real,dimension(Nxb+2,Nyb+2), intent(in) :: ut
      real, intent(in) :: dx
      real, intent(in) :: dy
      real, intent(in) :: inRe
      real,dimension(Nxb,Nyb) :: uP,uN,uS,uE,uW
      real, dimension(Nxb,Nyb), intent(out) :: D1

      uP = ut(2:Nxb+1,2:Nyb+1)
      uE = ut(3:Nxb+2,2:Nyb+1)
      uW = ut(1:Nxb,2:Nyb+1)
      uN = ut(2:Nxb+1,3:Nyb+2)
      uS = ut(2:Nxb+1,1:Nyb)

      D1 = ((inRe*(uE-uP)/dx)-(inRe*(uP-uW)/dx))/dx + &
           ((inRe*(uN-uP)/dy)-(inRe*(uP-uS)/dy))/dy

end subroutine Diffusive_U

subroutine Diffusive_V(vt,dx,dy,inRe,D2,Nxb,Nyb)

      implicit none
 
      integer, intent(in) :: Nxb,Nyb
      real,dimension(Nxb+2,Nyb+2), intent(in) :: vt
      real, intent(in) :: dx
      real, intent(in) :: dy
      real, intent(in) :: inRe
      real,dimension(Nxb,Nyb) :: vP,vE,vW,vN,vS
      real, dimension(Nxb,Nyb), intent(out) :: D2

      vP = vt(2:Nxb+1,2:Nyb+1)
      vE = vt(3:Nxb+2,2:Nyb+1)
      vW = vt(1:Nxb,2:Nyb+1)
      vN = vt(2:Nxb+1,3:Nyb+2)
      vS = vt(2:Nxb+1,1:Nyb)

      D2 = ((inRe*(vE-vP)/dx)-(inRe*(vP-vW)/dx))/dx + &
           ((inRe*(vN-vP)/dy)-(inRe*(vP-vS)/dy))/dy

end subroutine Diffusive_V

