import numpy as np
import matplotlib.pyplot as plt

def solve_2d_ftcs(Nx, Ny, Nt, alpha, T, Lx=1.0, Ly=1.0):
    dx, dy = Lx/Nx, Ly/Ny
    dt = T/Nt
    rx, ry = alpha*dt/dx**2, alpha*dt/dy**2
    
    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    X, Y = np.meshgrid(x, y)
    
    u = np.sin(np.pi*X)*np.sin(np.pi*Y)
    u_new = np.copy(u)
    
    for n in range(Nt):
        u_new[1:-1,1:-1] = u[1:-1,1:-1] + rx*(u[1:-1,2:] - 2*u[1:-1,1:-1] + u[1:-1,:-2]) \
                                        + ry*(u[2:,1:-1] - 2*u[1:-1,1:-1] + u[:-2,1:-1])
        u_new[0,:] = u_new[-1,:] = u_new[:,0] = u_new[:,-1] = 0.0
        u[:] = u_new[:]
    
    return X, Y, u

def solve_2d_cn(Nx, Ny, Nt, alpha, T, Lx=1.0, Ly=1.0):
    dx, dy = Lx/Nx, Ly/Ny
    dt = T/Nt
    rx, ry = alpha*dt/(2*dx**2), alpha*dt/(2*dy**2)
    
    x = np.linspace(0, Lx, Nx+1)
    y = np.linspace(0, Ly, Ny+1)
    X, Y = np.meshgrid(x, y)
    
    u = np.sin(np.pi*X)*np.sin(np.pi*Y)
    
    for n in range(Nt):
        u_temp = np.copy(u)
        for it in range(50):  # bump to 50 for convergence
            u_new = np.copy(u_temp)
            u_new[1:-1,1:-1] = (u[1:-1,1:-1] + rx*(u_temp[1:-1,2:] + u_temp[1:-1,:-2]) 
                                           + ry*(u_temp[2:,1:-1] + u_temp[:-2,1:-1])) / (1+2*rx+2*ry)
            u_temp[:] = u_new[:]
        u[:] = u_new[:]
        u[0,:] = u[-1,:] = u[:,0] = u[:,-1] = 0.0
    
    return X, Y, u

def analytical_2d(X, Y, t, alpha=1.0):
    return np.exp(-2*alpha*np.pi**2*t)*np.sin(np.pi*X)*np.sin(np.pi*Y)

def convergence_study():
    Nx_list = [20, 40, 80, 160]
    T = 0.1
    alpha = 1.0
    
    errors_ftcs = []
    dx_list = []
    
    print(f"{'Nx':<6} {'dx':<8} {'L2 Error':<12} {'Rate'}")
    print("-"*35)
    
    prev_err = None
    for Nx in Nx_list:
        Ny = Nx
        dx = 1.0/Nx
        dx_list.append(dx)
        
        # Use dt = 0.25*dx² for stability and fair comparison
        dt = 0.25*alpha*dx**2
        Nt = int(np.ceil(T/dt))
        
        X, Y, u_num = solve_2d_ftcs(Nx, Ny, Nt, alpha, T)
        u_exact = analytical_2d(X, Y, T, alpha)
        err_l2 = np.sqrt(np.mean((u_num - u_exact)**2))
        errors_ftcs.append(err_l2)
        
        if prev_err is None:
            rate_str = "-"
        else:
            rate = np.log(prev_err/err_l2)/np.log(2)
            rate_str = f"{rate:.3f}"
        
        print(f"{Nx:<6} {dx:<8.4f} {err_l2:<12.3e} {rate_str}")
        prev_err = err_l2
    
    # Plot 1: Convergence
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.loglog(dx_list, errors_ftcs, 'o-', label='FTCS', linewidth=2)
    plt.loglog(dx_list, 1.234e-2*(np.array(dx_list)/0.05)**2, 'k--', label='O(dx²)')
    plt.xlabel('dx'); plt.ylabel('L2 Error'); plt.legend(); plt.grid(True, which='both', ls='--')
    plt.title('Spatial Convergence at t=0.1')
    
    # Plot 2: Heatmap of error - fix orientation
    Nx, Ny = 80, 80
    dt = 0.25*alpha*(1.0/Nx)**2
    Nt = int(np.ceil(T/dt))
    X, Y, u_ftcs = solve_2d_ftcs(Nx, Ny, Nt, alpha, T)
    u_exact = analytical_2d(X, Y, T, alpha)
    error_field = np.abs(u_ftcs - u_exact)
    
    plt.subplot(1,2,2)
    im = plt.imshow(error_field, cmap='hot', extent=[0,1,0,1], origin='lower', interpolation='nearest')
    plt.colorbar(im, label='|Error|')
    plt.xlabel('x'); plt.ylabel('y')
    plt.title(f'FTCS Error, Nx={Nx}, t=0.1')
    
    plt.tight_layout()
    plt.savefig('method_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    convergence_study()