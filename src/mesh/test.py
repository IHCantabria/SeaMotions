
# import matplotlib.pyplot as plt
# import numpy as np

# r0 = np.array( [ -3.0, -2.5 ] )
# r1 = np.array( [  1.0,  1.0 ] )

# r0r1 = r1 - r0

# lm = -r0[1] / ( r1[1] - r0[1] )

# r0r1v = r0 + lm * r0r1
# # r0r1v = r0 + r0r1

# print( "r0:    ", r0 )
# print( "r0r1:  ", r0r1 )
# print( "r0r1v: ", r0r1v )
# print( "lm: ", lm )

# # Use only X (index 0) and Z (index 2)
# x0, z0 = r0[0], r0[1]
# x1, z1 = r1[0], r1[1]

# plt.figure()

# # Points
# plt.scatter([x0, x1], [z0, z1])

# # Vector (arrow)
# plt.arrow(0.0, 0.0, x0, z0, length_includes_head=True, head_width=0.05)
# plt.arrow(0.0, 0.0, x1, z1, length_includes_head=True, head_width=0.05)
# plt.arrow(x0, z0, r0r1[0], r0r1[1], length_includes_head=True, head_width=0.05)
# plt.arrow(0.0, 0.0, r0r1v[0], r0r1v[1], length_includes_head=True, head_width=0.05)

# plt.xlabel("X")
# plt.ylabel("Z")
# plt.axis("equal")
# plt.grid(True)

# plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Parameters
rho = 1000
g = 9.81
R = 1.0
draft = -0.052
N = 180

theta = np.linspace(0, 2*np.pi, N+1)

# 4‑point Gauss–Legendre quadrature
xi = np.array([-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116])
weights = np.array([0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451])

panel_centers = []
panel_pressures = []

Fx_total = 0
Fy_total = 0

fig = plt.figure( )
for i in range(N):
    th0, th1 = theta[i], theta[i+1]
    L = R * (th1 - th0)

    th_mid = 0.5*(th0 + th1)
    th_half = 0.5*(th1 - th0)

    Fx = 0
    Fy = 0
    p_sum = 0
    plt.plot( R*np.cos( th0 ), R*np.sin( th0 ), "o", color="blue" )
    plt.plot( R*np.cos( th1 ), R*np.sin( th1 ), "o", color="blue" )
    for k in range(4):
        th = th_mid + xi[k]*th_half
        x = R*np.cos(th)
        y = R*np.sin(th)

        p = rho*g*(draft - y)
        if ( p > 0 ):
            plt.plot( x, y, "x", color="red" )
        
        p = max(p, 0)

        nx = np.cos(th)
        ny = np.sin(th)

        Fx += p * nx * L * weights[k] / 2.0
        Fy += p * ny * L * weights[k] / 2.0

        p_sum += p

    Fx_total += Fx
    Fy_total += Fy

    panel_centers.append([R*np.cos(th_mid), R*np.sin(th_mid)])
    panel_pressures.append(p_sum/4)
plt.show( )

panel_centers = np.array(panel_centers)
panel_pressures = np.array(panel_pressures)

fig, axs = plt.subplots(1, 2, figsize=(12,6))

axs[0].set_title("Panel Mesh on Circle")
axs[0].plot(np.cos(theta), np.sin(theta))
axs[0].scatter(panel_centers[:,0], panel_centers[:,1])
axs[0].set_aspect('equal')

axs[1].set_title("Pressure Distribution on Panels")
sc = axs[1].scatter(panel_centers[:,0], panel_centers[:,1], c=panel_pressures)
axs[1].set_aspect('equal')
plt.colorbar(sc, ax=axs[1])

plt.tight_layout()
plt.show()

print(Fx_total, Fy_total)
