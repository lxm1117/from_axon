import numpy as np

fw=open('soma_v_burst.txt','w')
v0=-76.9248
vclamp=0
dt=0.02
tstop=200
total_len=int(tstop/dt)
#v=np.zeros(total_len)

t1=50
t2=150
for i in range(0,int(t1/dt)):
	fw.write(str(i*dt)+'\t'+str(v0)+'\n')

for i in range(int(t1/dt), int(t2/dt)):
	fw.write(str(i*dt)+'\t'+str(vclamp)+'\n')

for i in range(int(t2/dt), total_len):
	fw.write(str(i*dt)+'\t'+str(v0)+'\n')

fw.close()
