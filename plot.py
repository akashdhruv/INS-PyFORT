import matplotlib.pyplot as plt
import numpy as np

M=20+1
N=20+1

k = 2
d = 2

X=np.zeros((N*d,M*k),dtype=float)
Y=np.zeros((N*d,M*k),dtype=float)
U=np.zeros((N*d,M*k),dtype=float)
V=np.zeros((N*d,M*k),dtype=float)
P=np.zeros((N*d,M*k),dtype=float)
T=np.zeros((N*d,M*k),dtype=float)

for i in range(k*d):

	data=np.loadtxt('LidData0%d.dat' % i)

	x = data[:,0]
	y = data[:,1]
	u = data[:,2]
	v = data[:,3]
	p = data[:,4]
	t = data[:,5]

	x = np.reshape(x,[N,M])
	y = np.reshape(y,[N,M])
	u = np.reshape(u,[N,M])
	v = np.reshape(v,[N,M])
	p = np.reshape(p,[N,M])
	t = np.reshape(t,[N,M])

	X[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=x
	Y[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=y
	U[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=u
	V[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=v
	P[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=p
	T[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=t

plt.figure()
plt.title('Resultant Velocity')
plt.contourf(X,Y,np.sqrt(U**2+V**2),density=5)
plt.axis('equal')


plt.figure()
plt.title('Pressure')
plt.contourf(X,Y,P,density=5)
plt.axis('equal')

plt.figure()
plt.title('Temperature')
plt.contourf(X,Y,T,density=5)
plt.axis('equal')

"""
plt.figure()
plt.contourf(X,Y,U,density=5)
plt.axis('equal')

plt.figure()
plt.contourf(X,Y,V,density=5)
plt.axis('equal')
"""

plt.show()
