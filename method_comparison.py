import numpy as np
import matplotlib.pyplot as plt

def solve_2d_ftcs_unstable(Nx, Ny, Nt, alpha, T, Lx=1.0, Ly=1.0):
    dx, dy = Lx/Nx, Ly/Ny
    dt = T/Nt
    rx, ry = alpha*dt/dx**2, alpha*dt/dy**2

    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    X, Y = np.meshgrid(x, y, indexing='ij')
    u = np.sin(np.pi*X)*np.sin(np.pi*Y)
    u_new = np.copy(u)

    for n in range(Nt):
        u_new[1:-1,1:-1] = u[1:-1,1:-1] + rx*(u[2:,1:-1] - 2*u[1:-1,1:-1] + u[:-2,1:-1]) \
                                        + ry*(u[1:-1,2:] - 2*u[1:-1,1:-1] + u[1:-1,:-2])
        u_new[0,:] = u_new[-1,:] = u_new[:,0] = u_new[:,-1] = 0.0
        u[:] = u_new[:]

    return X, Y, u

def solve_2d_cn_stable(Nx, Ny, Nt, alpha, T, Lx=1.0, Ly=1.0, max_iter=50, tol=1e-6):
    dx, dy = Lx/Nx, Ly/Ny
    dt = T/Nt
    rx, ry = alpha*dt/(2*dx**2), alpha*dt/(2*dy**2)

    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    X, Y = np.meshgrid(x, y, indexing='ij')
    u = np.sin(np.pi*X)*np.sin(np.pi*Y)

    for n in range(Nt):
        rhs = u.copy()
        rhs[1:-1,1:-1] += rx*(u[2:,1:-1] - 2*u[1:-1,1:-1] + u[:-2,1:-1]) \
                        + ry*(u[1:-1,2:] - 2*u[1:-1,1:-1] + u[1:-1,:-2])

        u_new = u.copy()
        for it in range(max_iter):
            u_old = u_new.copy()
            u_new[1:-1,1:-1] = (rhs[1:-1,1:-1] + rx*(u_old[2:,1:-1] + u_old[:-2,1:-1]) \
                                             + ry*(u_old[1:-1,2:] + u_old[1:-1,:-2])) / (1 + 2*rx + 2*ry)
            u_new[0,:] = u_new[-1,:] = u_new[:,0] = u_new[:,-1] = 0.0
            if np.max(np.abs(u_new - u_old)) < tol:
                break
        u[:] = u_new[:]

    return X, Y, u

def analytical_2d(X, Y, t, alpha=1.0):
    return np.exp(-2*alpha*np.pi**2*t)*np.sin(np.pi*X)*np.sin(np.pi*Y)

Nx, Ny = 60, 60
Lx, Ly, T, alpha = 1.0, 1.0, 0.1, 1.0
dx, dy = Lx/Nx, Ly/Ny

target_rx_ry = 0.64
dt = target_rx_ry * dx**2 / (2*alpha)
Nt = int(np.ceil(T/dt))

X, Y, u_ftcs = solve_2d_ftcs_unstable(Nx, Ny, Nt, alpha, T)
X, Y, u_cn = solve_2d_cn_stable(Nx, Ny, Nt, alpha, T)
u_exact = analytical_2d(X, Y, T, alpha)

err_ftcs_signed = u_ftcs - u_exact
err_ftcs_signed = np.nan_to_num(err_ftcs_signed, nan=1.0, posinf=1.0, neginf=-1.0)
err_ftcs_signed = np.clip(err_ftcs_signed, -1.0, 1.0)

err_cn = np.abs(u_cn - u_exact)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

im0 = axes[0].imshow(err_ftcs_signed, cmap='RdBu_r', extent=[0,1,0,1], origin='lower', vmin=-1, vmax=1)
axes[0].set_title('FTCS: Unstable')
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')
axes[0].text(0.5, 0.9, 'Checkerboard oscillations',
             ha='center', color='black', fontsize=10, bbox=dict(facecolor='white', alpha=0.7))

im1 = axes[1].imshow(err_cn, cmap='hot', extent=[0,1,0,1], origin='lower')
axes[1].set_title('Crank-Nicolson: Stable')
axes[1].set_xlabel('x')
axes[1].set_ylabel('y')

cbar = fig.colorbar(im1, ax=axes, shrink=0.8)
cbar.set_label('') # removed label

plt.suptitle(f'Stability Comparison at t={T}, rx+ry={target_rx_ry:.2f}', fontsize=12)
plt.subplots_adjust(wspace=0.3, top=0.85)
plt.savefig('method_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"Used: Nx={Nx}, dx={dx:.4f}, dt={dt:.5f}, rx+ry={target_rx_ry:.2f}")