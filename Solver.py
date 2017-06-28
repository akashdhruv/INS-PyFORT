# Importing Libraries
import INS
import POISSON
import HEAT
import numpy as np
from math import *
from mpi4py import MPI
import time
#import matplotlib.pyplot as plt


#__________________Defining MPI communication function______________________#

def MPI_applyBC(u,x_id,y_id,x_procs,y_procs,x_comm,y_comm):

	if(x_procs > 1):

		if(x_id%2 == 0):

			if(x_id == 0):

				x_comm.send(u[-2,:],dest=(x_id+1)%x_procs,tag=1)
				u[-1,:] = x_comm.recv(source=(x_id+1)%x_procs,tag=2)

			elif(x_id == nblockx - 1):

				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=3)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=4)

			else:
				x_comm.send(u[-2,:],dest=(x_id+1)%x_procs,tag=1)
				u[-1,:] = x_comm.recv(source=(x_id+1)%x_procs,tag=2)

				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=3)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=4)

		elif(x_id%2 == 1):

			if(x_id == nblockx - 1):
				
				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=2)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=1)

			else:

				x_comm.send(u[1,:],dest=(x_id-1+x_procs)%x_procs,tag=2)
				u[0,:] = x_comm.recv(source=(x_id-1+x_procs)%x_procs,tag=1)
				
				x_comm.send(u[-2,:],dest=(x_id+1)%x_procs,tag=4)	
				u[-1,:] = x_comm.recv(source=(x_id+1)%x_procs,tag=3)

	if(y_procs > 1):

		if(y_id%2 == 0):

			if(y_id == 0):

				y_comm.send(u[:,-2],dest=(y_id+1)%y_procs,tag=5)
				u[:,-1] = y_comm.recv(source=(y_id+1)%y_procs,tag=6)

			elif(y_id == nblocky - 1):

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=7)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=8)

			else:

				y_comm.send(u[:,-2],dest=(y_id+1)%y_procs,tag=5)
				u[:,-1] = y_comm.recv(source=(y_id+1)%y_procs,tag=6)

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=7)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=8)


		elif(y_id%2 == 1):

			if(y_id == nblocky - 1):

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=6)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=5)
		
			else:

				y_comm.send(u[:,1],dest=(y_id-1+y_procs)%y_procs,tag=6)
				u[:,0] = y_comm.recv(source=(y_id-1+y_procs)%y_procs,tag=5)


				y_comm.send(u[:,-2],dest=(y_id+1)%y_procs,tag=8)
				u[:,-1] = y_comm.recv(source=(y_id+1)%y_procs, tag=7)

	return u

##_______________________________________MAIN_______________________________________________#

#_____________________________Initializing MPI environment__________________________________#

nblockx = 2
nblocky = 2

comm = MPI.COMM_WORLD
myid = comm.Get_rank()
procs = comm.Get_size()

x_comm = comm.Split(myid/nblockx,myid%nblockx)
y_comm = comm.Split(myid%nblockx,myid/nblockx)

x_id = x_comm.Get_rank()
x_procs = x_comm.Get_size()

y_id = y_comm.Get_rank()
y_procs = y_comm.Get_size()

t1 = MPI.Wtime()

#______________________________Domain Length and Limits_____________________________________#

Dx_min = -0.5
Dx_max =  0.5

Dy_min = -0.5
Dy_max =  0.5

Lx = Dx_max - Dx_min
Ly = Dy_max - Dy_min

gr_Lx = Lx/nblockx
gr_Ly = Ly/nblocky

#______________________________________Block size__________________________________________#

Nxb = 20
Nyb = 20

dx = gr_Lx/Nxb
dy = gr_Ly/Nyb

#_______________________________________Constants__________________________________________#

PRES_VAR = 0
TEMP_VAR = 1
PNEW_VAR = 2
TNEW_VAR = 3

CENT_VAR = 4

VELC_VAR = 0
VSTR_VAR = 1
VOLD_VAR = 2

FACE_VAR = 3

GONE_VAR = 0
GTWO_VAR = 1
G1NW_VAR = 2
G2NW_VAR = 3
PRHS_VAR = 4

WORK_VAR = 5

#___________________________________physical variables___________________________________#


x = Dx_min + (myid%nblockx)*gr_Lx + dx*np.linspace(0,Nxb,Nxb+1)

y = Dy_min + (myid/nblockx)*gr_Ly + dy*np.linspace(0,Nyb,Nyb+1)

[X,Y] = np.meshgrid(x,y)

center = np.zeros((CENT_VAR,Nxb+2,Nyb+2),dtype=float)
facex  = np.zeros((FACE_VAR,Nxb+2,Nyb+2),dtype=float)
facey  = np.zeros((FACE_VAR,Nxb+2,Nyb+2),dtype=float)
work   = np.zeros((WORK_VAR,Nxb,Nyb),dtype=float)

center[TEMP_VAR,:,:] = 313.0

#___________________________________ins parameters______________________________________#

ins_inRe = 0.001  
ins_sig  = 0.01
ins_cfl  = 0.15

#__________________________________heat parameters______________________________________#

ht_Pr = 0.7

#_________________________________driver parameters_____________________________________#

dt_sig = ins_sig*(min(dx,dy)**2)/ins_inRe
dt_cfl = ins_cfl*min(dx,dy)
dt_temp = dt_sig*ht_Pr


dt = min(dt_sig,dt_cfl)
dt = min(dt,dt_temp)

t = 60.0

nt = int(t/dt)

Maxit = 1500

p_res  = 0.
u_res  = 0.
v_res  = 0.
T_res  = 0.
maxdiv = 0.
mindiv = 0.

ins_p_res  = 0.
ins_v_res  = 0.
ins_v_res  = 0.
ins_T_res  = 0.
ins_maxdiv = 0.
ins_mindiv = 0.

#________________________________Physics Squence_____________________________________#

tstep = 0

while(tstep<=nt):

	facex[VOLD_VAR,:,:] = facex[VELC_VAR,:,:]
	facey[VOLD_VAR,:,:] = facey[VELC_VAR,:,:]
	
	#_____________________________Predictor_____________________________#

	#G1_new,G2_new,ut,vt = INS.predictor(u,v,G1,G2,ins_inRe,dx,dy,dt,tstep,Nxb,Nyb)
	work[G1NW_VAR,:,:],work[G2NW_VAR,:,:],facex[VSTR_VAR,:,:],facey[VSTR_VAR,:,:] = INS.predictor(facex[VELC_VAR,:,:],facey[VELC_VAR,:,:],work[GONE_VAR,:,:],work[GTWO_VAR,:,:],ins_inRe,dx,dy,dt,tstep,Nxb,Nyb)

        work[GONE_VAR,:,:] = work[G1NW_VAR,:,:]
        work[GTWO_VAR,:,:] = work[G2NW_VAR,:,:]

        #__________________Predictor Boundary Conditions_____________________#

	facex[VSTR_VAR,:,:] = MPI_applyBC(facex[VSTR_VAR,:,:],x_id,y_id,x_procs,y_procs,x_comm,y_comm)
	facey[VSTR_VAR,:,:] = MPI_applyBC(facey[VSTR_VAR,:,:],x_id,y_id,x_procs,y_procs,x_comm,y_comm)

        # LOW X
	if(x_id == 0):
		facex[VSTR_VAR,0,:]  =  0.0
        	facey[VSTR_VAR,0,:]  = -facey[VSTR_VAR,1,:]

        # HIGH X
	if(x_id == nblockx-1):
        	facex[VSTR_VAR,-2,:] =  0.0
        	facex[VSTR_VAR,-1,:] =  0.0	
        	facey[VSTR_VAR,-1,:] = -facey[VSTR_VAR,-2,:]

        # LOW Y
	if(y_id == 0):
        	facey[VSTR_VAR,:,0]  =  0.0
        	facex[VSTR_VAR,:,0]  = -facex[VSTR_VAR,:,1]

        # HIGH Y
	if(y_id == nblocky-1):
        	facey[VSTR_VAR,:,-1] =  0.0
        	facey[VSTR_VAR,:,-2] =  0.0
        	facex[VSTR_VAR,:,-1] = 2.0 -facex[VSTR_VAR,:,-2]

        #_____________________________Poisson Solver________________________#

        work[PRHS_VAR,:,:] = -((1/(dy*dt))*(facey[VSTR_VAR,1:-1,1:-1]-facey[VSTR_VAR,1:-1,:-2]))-((1/(dx*dt))*(facex[VSTR_VAR,1:-1,1:-1]-facex[VSTR_VAR,:-2,1:-1]))

        p_counter = 0

        while(p_counter < Maxit):

        	#p_new = POISSON.solver(p,p_RHS,dx,dy,Nxb,Nyb)
		center[PNEW_VAR,:,:] = POISSON.solver(center[PRES_VAR,:,:],work[PRHS_VAR,:,:],dx,dy,Nxb,Nyb)

		#___________________Pressure Boundary Conditions____________#

		center[PNEW_VAR,:,:] = MPI_applyBC(center[PNEW_VAR,:,:],x_id,y_id,x_procs,y_procs,x_comm,y_comm)

                # LOW X
		if(x_id == 0):	
			center[PNEW_VAR,0,:]  =  center[PNEW_VAR,1,:]
                       
		# HIGH X
		if(x_id == nblockx-1):
			center[PNEW_VAR,-1,:] =  center[PNEW_VAR,-2,:]

		# LOW Y
		if(y_id == 0):
			center[PNEW_VAR,:,0]  =  center[PNEW_VAR,:,1] 

		# HIGH Y
		if(y_id == nblocky-1):
			center[PNEW_VAR,:,-1] =  center[PNEW_VAR,:,-2]  
                 
                #_________________Residuals and Convergence Check__________#

                p_res = np.sum((center[PNEW_VAR,:,:]-center[PRES_VAR,:,:])**2)

                center[PRES_VAR,:,:] = center[PNEW_VAR,:,:]

                p_counter += 1
              
		ins_p_res = comm.allreduce(p_res, op=MPI.SUM)
		ins_p_res = sqrt(ins_p_res/((Nxb+2)*(Nyb+2)*procs))

                if(ins_p_res<10**-6 and ins_p_res != 0.):
			break

	#________________________________Corrector____________________________#

	facex[VELC_VAR,:,:],facey[VELC_VAR,:,:] = INS.corrector(facex[VSTR_VAR,:,:],facey[VSTR_VAR,:,:],center[PRES_VAR,:,:],dt,dx,dy,Nxb,Nyb)

        #__________________Corrector Boundary Conditions_____________________#

	facex[VELC_VAR,:,:] = MPI_applyBC(facex[VELC_VAR,:,:],x_id,y_id,x_procs,y_procs,x_comm,y_comm)
	facey[VELC_VAR,:,:] = MPI_applyBC(facey[VELC_VAR,:,:],x_id,y_id,x_procs,y_procs,x_comm,y_comm)

        # LOW X
	if(x_id == 0):
		facex[VELC_VAR,0,:]  =  0.0
		facey[VELC_VAR,0,:]  = -facey[VELC_VAR,1,:]

        # HIGH X
	if(x_id == nblockx - 1):
		facex[VELC_VAR,-2,:] =  0.0
		facex[VELC_VAR,-1,:] =  0.0
		facey[VELC_VAR,-1,:] = -facey[VELC_VAR,-2,:]

        # LOW Y
	if(y_id == 0):
		facey[VELC_VAR,:,0]  =  0.0
		facex[VELC_VAR,:,0]  = -facex[VELC_VAR,:,1]

        # HIGH Y
	if(y_id == nblocky - 1):
		facey[VELC_VAR,:,-1] =  0.0
		facey[VELC_VAR,:,-2] =  0.0
		facex[VELC_VAR,:,-1] = 2.0 - facex[VELC_VAR,:,-2]

        #___________________________Residuals_______________________________#

	u_res = np.sum((facex[VOLD_VAR,:,:]-facex[VELC_VAR,:,:])**2)
	v_res = np.sum((facey[VOLD_VAR,:,:]-facey[VELC_VAR,:,:])**2)

	ins_u_res = comm.allreduce(u_res, op=MPI.SUM)
	ins_u_res = sqrt(ins_u_res/((Nxb+2)*(Nyb+2)*procs))

	ins_v_res = comm.allreduce(v_res, op=MPI.SUM)
	ins_v_res = sqrt(ins_v_res/((Nxb+2)*(Nyb+2)*procs))

        #____________________________Divergence_____________________________#

	maxdiv = -10.0**(10)
	mindiv =  10.0**(10)

	maxdiv = max(maxdiv,np.max(((1/(dy))*(facey[VELC_VAR,1:-1,1:-1]-facey[VELC_VAR,1:-1,:-2])) + ((1/(dx))*(facex[VELC_VAR,1:-1,1:-1]-facex[VELC_VAR,:-2,1:-1]))))
	mindiv = min(mindiv,np.min(((1/(dy))*(facey[VELC_VAR,1:-1,1:-1]-facey[VELC_VAR,1:-1,:-2])) + ((1/(dx))*(facex[VELC_VAR,1:-1,1:-1]-facex[VELC_VAR,:-2,1:-1]))))

	ins_maxdiv = comm.allreduce(maxdiv, op=MPI.MAX)
	ins_mindiv = comm.allreduce(mindiv, op=MPI.MIN)

        #_______________________Heat Advection Diffusion____________________#

        center[TNEW_VAR,:,:] = HEAT.tempsolver(center[TEMP_VAR,:,:],facex[VELC_VAR,:,:],facey[VELC_VAR,:,:],dx,dy,dt,ins_inRe,ht_Pr,Nxb,Nyb)

        #____________________Temperature Boundary Conditions________________#

	center[TNEW_VAR,:,:] = MPI_applyBC(center[TNEW_VAR,:,:],x_id,y_id,x_procs,y_procs,x_comm,y_comm)

	# LOW X
	if(x_id == 0):
		center[TNEW_VAR,0,:]  =  center[TNEW_VAR,1,:]

	# HIGH X
	if(x_id == nblockx - 1):
		center[TNEW_VAR,-1,:] =  center[TNEW_VAR,-2,:]

	# LOW Y
	if(y_id == 0):
		center[TNEW_VAR,:,0]  =  center[TNEW_VAR,:,1]

	# HIGH Y
	if(y_id == nblocky - 1):
		center[TNEW_VAR,:,-1] = 2*383.15 - center[TNEW_VAR,:,-2]

        #___________________________Residuals_______________________________#

        T_res = np.sum((center[TNEW_VAR,:,:]-center[TEMP_VAR,:,:])**2)

	ins_T_res = comm.allreduce(T_res, op=MPI.SUM)
	ins_T_res = sqrt(ins_T_res/((Nxb+2)*(Nyb+2)*procs))
	
	center[TEMP_VAR,:,:] = center[TNEW_VAR,:,:]

        #____________________________Display_________________________________#

	if(myid == 0 and tstep%5 == 0):

		print "---------------------------PARAMETER DISPLAY-----------------------"
		print "Simulation Time     : ",tstep*dt," s"
		print "U velocity Residual : ",ins_u_res
		print "V velocity Residual : ",ins_v_res
        	print "Temperature Residual: ",ins_T_res
		print "Pressure Residual   : ",ins_p_res
		print "Poisson Counter     : ",p_counter
		print "MAXDIV : ",ins_maxdiv," MINDIV: ",ins_mindiv

        #__________________________Convergence Check_________________________#

	tstep += 1

	if(ins_u_res<10**-7 and ins_u_res != 0. and ins_v_res<10**-7 and ins_v_res != 0.):
		break



#_________________________Post Processing and Writing Data to File_____________#

uu = 0.5*(facex[VELC_VAR,:-1,:-1] + facex[VELC_VAR,:-1,1:])
vv = 0.5*(facey[VELC_VAR,:-1,:-1] + facey[VELC_VAR,1:,:-1])
pp = 0.25*(center[PRES_VAR,:-1,:-1] + center[PRES_VAR,1:,:-1] + center[PRES_VAR,:-1,1:] + center[PRES_VAR,1:,1:])
tt = 0.25*(center[TEMP_VAR,:-1,:-1] + center[TEMP_VAR,1:,:-1] + center[TEMP_VAR,:-1,1:] + center[TEMP_VAR,1:,1:])


uu = uu.T
vv = vv.T
pp = pp.T
tt = tt.T

X = np.reshape(X,np.size(X))
Y = np.reshape(Y,np.size(Y))
uu = np.reshape(uu,np.size(uu))
vv = np.reshape(vv,np.size(vv))
pp = np.reshape(pp,np.size(pp))
tt = np.reshape(tt,np.size(tt))

DataOut = np.column_stack((X.T,Y.T,uu.T,vv.T,pp.T,tt.T))
np.savetxt('LidData0%d.dat' % myid,DataOut)

t2 = MPI.Wtime()
print t2-t1

