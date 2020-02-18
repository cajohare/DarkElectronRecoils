from numpy import *
from HaloFuncs import *
from Params import *
from AtomicFuncs import *

# s orbitals
n_s = array([1]+[2]*2+[3]*3+[4]*3+[5]*4)
Z_s = array([54.9179,47.2500,26.0942,68.1771,16.8296,12.0759,31.9030,8.0145,5.8396,14.7123,3.8555,2.6343,1.8124])
c_1s = array([-0.965401,-0.040350,0.001890,-0.003868,-0.000263,0.000547,-0.000791,0.000014,-0.000013,-0.000286,0.000005,-0.000003,0.000001])
c_2s = array([0.313912,0.236118,-0.985333,0.000229,-0.346825,0.345786,-0.120941,-0.005057,0.001528,-0.151508,-0.000281,0.000134,-0.000040])
c_3s = array([-0.140382,-0.125401,0.528161,-0.000435,0.494492,-1.855445,0.128637,-0.017980,0.000792,0.333907,-0.000228,0.000191,-0.000037])
c_4s = array([0.064020,0.059550,-0.251138,0.000152,-0.252274,1.063559,-0.071737,-0.563072,-0.697466,-0.058009,-0.018353,0.002292,-0.000834])
c_5s = array([-0.022510,-0.021077,0.088978,-0.000081,0.095199,-0.398492,0.025623,0.274471,0.291110,0.011171,-0.463123,-0.545266,-0.167779])

# p orbitals
n_p = array([2]*2+[3]*3+[4]*3+[5]*4)
Z_p = array([58.7712,22.6065,48.9702,13.4997,9.8328,40.2591,7.1841,5.1284,21.5330,3.4469,2.2384,1.14588])
c_2p = array([0.051242,0.781070,0.114910,-0.000731,0.000458,0.083993,-0.000265,0.000034,0.009061,-0.000014,0.000006,-0.000002])
c_3p = array([0.000264,0.622357,-0.009861,-0.952677,-0.337900,-0.026340,-0.000384,-0.001665,0.087491,0.000240,-0.000083,0.000026])
c_4p = array([0.013769,-0.426955,0.045088,0.748434,0.132850,0.059406,-0.679569,-0.503653,-0.149635,-0.014193,0.000528,-0.000221])
c_5p = array([-0.005879,0.149040,-0.018716,-0.266839,-0.031096,-0.024100,0.267374,0.161460,0.059721,-0.428353,-0.542284,-0.201667])

# d orbitals
n_d = array([3]*3+[4]*5)
Z_d = array([19.9787,12.2129,8.6994,27.7398,15.9410,6.0580,4.0990,2.5857])
c_4d = array([-0.013758,-0.804573,0.260624,0.00749,0.244109,0.597018,0.395554,0.039786])
c_3d = array([0.220185,0.603140,0.194682,-0.014369,0.049865,-0.000300,0.000418,-0.000133])

n = int(input("Number of values"))

E_r_vals = logspace(-1.0,3.0,n)/1000.0 # keV
q_vals = logspace(0.0,4.0,n)

np = 50

F1 = zeros(shape=(n,n))
F2 = zeros(shape=(n,n))
F3 = zeros(shape=(n,n))
for i in range(0,n):
    F1[i,:] = f_nl_ion_sq(q_vals,E_r_vals[i],1,c_5p,n_p,Z_p,np=np)
    F2[i,:] = f_nl_ion_sq(q_vals,E_r_vals[i],0,c_5s,n_s,Z_s,np=np)
    F3[i,:] = f_nl_ion_sq(q_vals,E_r_vals[i],2,c_4d,n_d,Z_d,np=np)
    print(i,'of',n)

savetxt('../data/fion/fion_Xe.txt',vstack((E_r_vals,q_vals,log10(F1),log10(F2),log10(F3))),delimiter='\t',fmt="%1.16f")