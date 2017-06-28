subroutine solver(ps,ps_rhs,maxit,&
                  x_comm,y_comm,solver_comm,x_id,y_id,myid,&
                  x_procs,y_procs,procs,nblockx,nblocky,&
                  ps_new,p_counter,p_res,dx,dy,Nxb,Nyb)
   use mpi

   implicit none
   integer, intent(in) :: Nxb,Nyb
   real, intent(in) :: dx,dy
   real, intent(in), dimension(Nxb,Nyb) :: ps_rhs
   real, intent(in), dimension(Nxb+2,Nyb+2) :: ps
   real, intent(out), dimension(Nxb+2,Nyb+2) :: ps_new
   integer, intent(out) :: p_counter
   real, intent(out) :: p_res
   integer, intent(in) :: maxit
   integer, intent(in) :: solver_comm,x_comm,y_comm
   integer, intent(in) :: x_id,y_id,myid
   integer, intent(in) :: x_procs,y_procs,procs,nblockx,nblocky

   integer :: i,j
   real :: res
   integer :: counter = 0
   integer :: ierr
   real, dimension(Nxb+2,Nyb+2) :: ps_ex

   ps_new = ps
     

  do while(counter<maxit)

   do j=2,Nyb+1
      do i=2,Nxb+1

      ps_new(i,j)=((ps(i,j+1)/(dy*dy))+(ps_new(i,j-1)/(dy*dy))&
                  +(ps(i+1,j)/(dx*dx))+(ps_new(i-1,j)/(dx*dx))&
                  +ps_RHS(i-1,j-1))&
                  *(1/((1/(dx*dx))+(1/(dy*dy))+&
                  (1/(dx*dx))+(1/(dy*dy))))

      end do
   end do

   ps_ex = ps_new

   call apply_BC(ps_ex,ps_new,x_comm,y_comm,x_id,y_id,x_procs,y_procs,nblockx,nblocky,Nxb,Nyb)

   res = res + sum(sum((ps-ps_new)**2,1))

   call MPI_ALLREDUCE(res,p_res,1,MPI_REAL,MPI_SUM,solver_comm,ierr)

   p_res = sqrt(p_res/((Nxb+2)*(Nyb+2)*procs))

   counter = counter +1

   if(p_res .lt. 0.000001 .and. p_res .ne. 0.) exit

  end do

  p_counter = counter

end subroutine solver


subroutine apply_BC(ps_ex,ps_new,x_comm,y_comm,x_id,y_id,x_procs,y_procs,nblockx,nblocky,Nxb,Nyb)


      use mpi

      implicit none
      integer, intent(in) :: Nxb,Nyb
      real, intent(in), dimension(Nxb+2,Nyb+2) :: ps_ex
      real, intent(out),dimension(Nxb+2,Nyb+2) :: ps_new
      integer,intent(in) :: x_comm,y_comm,x_id,y_id,x_procs,y_procs,nblockx,nblocky


      integer :: status(MPI_STATUS_SIZE)
      integer :: ierr

      ps_new = ps_ex

      if (x_procs > 1) then

       if(mod(x_id,2) == 0) then

             if(x_id == 0) then

                 call MPI_SENDRECV(ps_new(Nxb+1,:), Nyb+2, MPI_REAL,mod(x_id+1,x_procs), 1,&
                                   ps_new(Nxb+2,:), Nyb+2, MPI_REAL,mod(x_id+1,x_procs), 2,x_comm, status, ierr)

             else if(x_id == nblockx-1) then

                call MPI_SENDRECV(ps_new(2,:), Nyb+2, MPI_REAL,mod(x_id-1+x_procs,x_procs), 3,&
                                  ps_new(1,:), Nyb+2, MPI_REAL,mod(x_id-1+x_procs,x_procs), 4,x_comm,status, ierr)

             else
                call MPI_SENDRECV(ps_new(Nxb+1,:), Nyb+2, MPI_REAL,mod(x_id+1,x_procs), 1,&
                                  ps_new(Nxb+2,:), Nyb+2, MPI_REAL,mod(x_id+1,x_procs), 2, x_comm, status, ierr)

                call MPI_SENDRECV(ps_new(2,:), Nyb+2, MPI_REAL,mod(x_id-1+x_procs,x_procs), 3,&
                                  ps_new(1,:), Nyb+2, MPI_REAL,mod(x_id-1+x_procs,x_procs), 4,x_comm,status, ierr)


             end if

       else if (mod(x_id,2) == 1) then

             if(x_id == nblockx-1) then

               call MPI_SENDRECV(ps_new(2,:), Nyb+2,MPI_REAL,mod(x_id-1+x_procs,x_procs), 2,&
                                 ps_new(1,:), Nyb+2,MPI_REAL,mod(x_id-1+x_procs,x_procs), 1,x_comm, status, ierr)

             else

               call MPI_SENDRECV(ps_new(2,:), Nyb+2, MPI_REAL,mod(x_id-1+x_procs,x_procs), 2,&
                                 ps_new(1,:), Nyb+2, MPI_REAL,mod(x_id-1+x_procs,x_procs), 1,x_comm,status, ierr)

               call MPI_SENDRECV(ps_new(Nxb+1,:), Nyb+2, MPI_REAL,mod(x_id+1,x_procs), 4,&
                                 ps_new(Nxb+2,:), Nyb+2, MPI_REAL,mod(x_id+1,x_procs), 3,x_comm,status,ierr)


             end if

       end if

     end if

     !! Second dimension !!

     if (y_procs > 1) then

       if(mod(y_id,2) == 0) then

             if(y_id == 0) then

                  call MPI_SENDRECV(ps_new(:,Nyb+1), Nxb+2, MPI_REAL,mod(y_id+1,y_procs), 5,&
                                    ps_new(:,Nyb+2), Nxb+2, MPI_REAL,mod(y_id+1,y_procs), 6,y_comm,status, ierr)

             else if(y_id == nblocky-1) then

                  call MPI_SENDRECV(ps_new(:,2), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 7,&
                                    ps_new(:,1), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 8,y_comm,status,ierr)

             else
                  call MPI_SENDRECV(ps_new(:,Nyb+1), Nxb+2, MPI_REAL,mod(y_id+1,y_procs), 5,&
                                    ps_new(:,Nyb+2), Nxb+2, MPI_REAL,mod(y_id+1,y_procs), 6,y_comm,status,ierr)

                  call MPI_SENDRECV(ps_new(:,2), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 7,&
                                    ps_new(:,1), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 8,y_comm,status,ierr)


             end if

       else if (mod(y_id,2) == 1) then

             if(y_id == nblocky-1) then

                  call MPI_SENDRECV(ps_new(:,2), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 6,&
                                    ps_new(:,1), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 5,y_comm,status, ierr)

             else

                  call MPI_SENDRECV(ps_new(:,2), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 6,&
                                    ps_new(:,1), Nxb+2, MPI_REAL,mod(y_id-1+y_procs,y_procs), 5,y_comm,status, ierr)


                  call MPI_SENDRECV(ps_new(:,Nyb+1), Nxb+2, MPI_REAL,mod(y_id+1,y_procs), 8,&
                                    ps_new(:,Nyb+2), Nxb+2, MPI_REAL,mod(y_id+1,y_procs), 7,y_comm,status, ierr)


             end if

       end if

    end if

    if ( x_id == 0) then

           ps_new(1,:)=ps_new(2,:)

    end if

    if ( x_id == nblockx-1) then

           ps_new(Nxb+2,:)=ps_new(Nxb+1,:)

     end if


    if ( y_id == 0) then

           ps_new(:,1)=ps_new(:,2)

    end if

    if ( y_id == nblocky-1) then

           ps_new(:,Nyb+2)=ps_new(:,Nyb+1)

    end if


end subroutine
