"""
First order IMEX solver for nonlinear FDEs
"""
struct IMEX_I <: FractionalDiffEqAlgorithm end

function solve(Rfun, t0, T, y0, N, alpha, lambda, Mlu, sigma1, Mf, sigma2,Mu,sigma3,Me,sigma4,tol,epsil, ::IMEX_I)

    
    h=(T-t0)/N
    t=collect(t0:h:T)
    d=length(y0)
    N = Int64(N)
    y = zeros(Int64(N+1), d)
    y[1, :] .= y0 # Handle the initial conditions
    
    #  Computing A and B coefficients for history load, which involves 
    # the computation of hypergeometric functions 
    #  a_{k,j}=\int_{t_k}^{t_k+1} (t_k+1 - s)^{alpha-1} (s - t_j)^{1-alpha} ds
    a = zeros(Int64(N-1), 5)
    # Gauss-Jacobi quadrature for hypergeometric functions
    a[1, 1]=1/alpha*HypergeometricFunctions._₂F₁(alpha-1, 1, alpha+1, -1)
    a[1, 2]=gamma(alpha)*gamma(2-alpha)/gamma(2)
    
    # Evaluating coefficients and employing symmetry relationships
    for k=2:3
        a[k, 1]=k^(1-alpha)*1/alpha*HypergeometricFunctions._₂F₁(alpha-1, 1, alpha+1, -1/k)
        a[k, 2]=a[k-1, 1]
        a[k, 3]=a[k-1, 2]
    end
    a[3, 4]=a[1, 2]
    for k=4:N-1
        a[k, 1]=k^(1-alpha)*1/alpha*HypergeometricFunctions._₂F₁(alpha-1, 1, alpha+1, -1/k)
        a[k, 2]=a[k-1, 1]
        a[k, 3]=a[k-1, 2]
        a[k, 4]=a[k-1, 3]
        a[k, 5]=a[k-1, 4]
    end
    
    # Checking for correction terms
    if Mlu>0 | Mf>0 | Mu>0 | Me>0
       # Computing the weights for lambda u
       V1 = zeros(Mlu)
       for i=1:Mlu
           for k=1:Mlu
               V1[i, k] = i^sigma1[k]
           end
       end
       rhs1 = zeros(N, Mlu)
       Flu=zeros(N, Mlu)
       auxgamma = gamma(alpha+1)
       for k=1:Mlu
           Flu[1, k]=gamma(sigma1[k]+1)/gamma(sigma1[k]+alpha+1) - 1/gamma(alpha+1)
           for j=2:N
               Flu[j, k]=(j-1)^sigma1(k)/auxgamma * HypergeometricFunctions._₂F₁(-sigma1[k], 1, alpha+1, -1/(j-1)) - j^sigma1[k]/auxgamma
           end
           rhs1[:, k]=Flu[:, k]
       end
       Wlu = rhs1*inv(V1) # Obtaining correction weights for linear term
       
       # Computing the weights for the force
       V2 = zeros(Mf)
       for i=1:Mf
           for k=1:Mf
               V2[i, k] = i^sigma2[k]
           end
       end
       rhs2 = zeros(N, Mf)
       Force=zeros(N, Mf)
       for k=1:Mf
           Force[1, k]=gamma(sigma2[k]+1)/gamma(sigma2[k]+alpha+1) - 1/gamma(alpha+1)
           for j=2:N
               Force[j, k]=(j-1)^sigma2[k]/auxgamma * HypergeometricFunctions._₂F₁(-sigma2[k], 1, alpha+1, -1/(j-1)) - j^sigma2[k]/auxgamma
           end
           rhs2[:, k]=Force[:, k]
       end
       Wf = rhs2*inv(V2) # Getting correction weights
       
       # Computing the weights for the history load
       V3 = zeros(Mu, Mu)
       for i=1:Mu
           for k=1:Mu
               V3[i, k] = i^sigma3[k]
           end
       end
       rhs3 = zeros(N, Mu)
       History=zeros(N, Mu)
       r=zeros(1, N)
       r[2]=a[1, 1]-a[1, 2]
       for i=3:N
           r[i]=a[i-1, 1]-2*a[i-1, 2]+a[i-1, 3]
       end
       
       auxgamma = gamma(alpha)*gamma(2-alpha)
       
       for k=1:Mu
           auxbeta = beta(sigma3[k]-alpha+1, alpha)
           
           # Performing fast history computation through FFT
           raux = [r[2:N]; 0; zeros(N-2,1)] # Circulant embedding (vector of unique coefficients)
           jaux = [((1:N-1).^sigma3[k]); zeros(N-1,1)]
           res = real.(ifft(fft(raux).*fft(jaux)))
           
           auxgamma1 = gamma(alpha)*gamma(sigma3[k]-alpha+2)
           auxgamma2 = gamma(1+sigma3[k])
           
           for j=2:N
               FF3 = res[j-1, 1]/auxgamma
    
               History[j, k] = -auxgamma2*j^(alpha-1)*((sigma3[k]-alpha+1)*(j)^(sigma3[k]-alpha+1)*beta_inc(sigma3(k)-alpha+1, alpha, (j-1)/j)*auxbeta)/auxgamma1 + (j-1)^sigma3[k] - FF3
    
           end
           rhs3[:, k]=History[:, k]
       end
       Wu = rhs3*inv(V3) # Correction weights for history term
       
       # Computing the weights for the extrapolated f_(k+1) 
       V4 = zeros(Me, Me)
       for i=1:Me
           for k=1:Me
               V4[i, k] = i^sigma4[k]
           end
       end
       rhs4 = zeros(N, Me)
       F4=zeros(N, Me)
       for k=1:Me
           for j=1:N
               F4[j, k]=j^sigma4[k] - (j-1)^sigma4[k]
           end
           rhs4[:, k]=F4[:, k]
       end
       We = rhs4*inv(V4)
        
       #  computing y_1,y_2,...,y_N
       #  (A-lambda*B)*Y=h^alpha*(B + C)*F + D (Eq 4.1 in Zhou et al.)
    
       A=zeros(N, 1)
       A[1]=1
       A[2]=-1+(a[1, 1]-a[1, 2])/auxgamma
       for i=3:N
           A[i]=(a[i-1, 1]-2*a[i-1, 2]+a[i-1, 3])/auxgamma
       end
       
       # B1 is stored as a sparse matrix
       B1col=h^alpha/gamma(alpha+1)*[1; zeros(N-1,1)]
       A[1] = A[1] - lambda * B1col[1]
       
       # Sparse Toeplitz allocation
       B1 = sptoeplitz(B1col, [B1col[1, 1], zeros(1, N-1)])
       
       # Similar procedure for B2
       B2col=h^alpha/gamma(alpha+1)*[-1; 1; zeros(N-2, 1)]
       B2 = sptoeplitz(B2col, [B2col[1, 1], zeros(1, N-1)])
       
       # Vector with initial conditions and force
       C=zeros(N, d)
       C[1, :] = y0 + h^alpha/gamma(alpha+1)*Rfun(t0, y0, alpha, lambda)#FIXME:
       for i=2:N
           C[i, :]=1/auxgamma*y0*(a[i-1, 1]-a[i-1, 2])
       end
       
       # Diagonal matrix for fast inverse approximation
       # D_delt = diag(1,delt,...,delt^(N-1))
       delt=(epsil)^(1/N)
       D=zeros(N, 1)
       for i=1:N
           D[i]=delt^(i-1)
       end
    
       # Computing Lambda = FN*D*K for fast Toeplitz inversion
       lmt=fft(D.*A)
       
       # Picard iteration algorithm
       e=ones(N, d)                 # Residual vector
       F=zeros(N+1, d)              # no need for circulant embedding here
       FC=zeros(N,d)               # Vector with correction terms
       y1=kron(ones(N, 1), y[1, :])
       auxgamma = gamma(alpha+1)
       countpicard = 0
    
       while norm(e)>tol
    
           Y=[y[1, :]; y1]    # initial iteration value set as initial condition
    
           for i=1:N+1
               F[i, :] = Rfun(t[i], Y[i, :], alpha, lambda)
           end
           
           # Vector with correction terms
           FC[1:N, :]=lambda*h^alpha*Wlu[1:N, 1:Mlu]*(Y[2:Mlu+1, :] - y0) + h^alpha/auxgamma*We[1:N, 1:Me]*(F[2:Me+1, :]-F[1, :]) + h^alpha*Wf[1:N, 1:Mf]*(F[2:Mf+1, :]-F[1, :]) - Wu[1:N, 1:Mu]*(Y[2:Mu+1, :] - y0)
           
           # Right-hand-side vector with sparse Toeplitz Matvec, vector of
           # initial conditions, and vector with correction terms
           b = B1*F[2:end, :] + B2*F[2:end, :] + C + FC
    
           # FFT of rhs for fast inversion (Eq.4.8 from Zhou et al.)
           b1=fft(D.*b)
    
           # diag(Lambda^-1)*(FN*D)*R
           b11=b1./lmt
    
           # Solution for current Picard iteration
           y2=real(ifft(b11)./D)
    
           # Residual
           e=y2-y1;
    
           # Update solution vector
           y1=y2;
    
           countpicard = countpicard + 1
       end
    
       # Attributing converged solution to output vector
       for i=2:N+1
           y[i, :]=y1[i-1, :]
       end

    else
       
       # No correction terms
    
       #  computing y_1,y_2,...,y_N
       #%%%%  (A-lambda*B)*Y=h^alpha*(B + C)*F (Eq 4.1 in Zhou et al. without D)
       auxgamma = gamma(alpha)*gamma(2-alpha)
       A=zeros(N, 1)
       A[1]=1
       A[2]=-1+(a[1, 1]-a[1, 2])/auxgamma
       for i=3:N
           A[i] =(a[i-1, 1]-2*a[i-1, 2]+a[i-1, 3])/auxgamma
       end
       
       # B1 is stored as a sparse matrix
       B1col=h^alpha/gamma(alpha+1)*[1; zeros(N-1)]
       A[1] = A[1] - lambda * B1col[1]
       
       # Sparse Toeplitz allocation
       B1 = Toeplitz(B1col, [B1col[1, 1]; zeros(N-1)])
       
       # Similar procedure for B2
       B2col=h^alpha/gamma(alpha+1)*[-1; 1; zeros(N-2)]
       B2 = Toeplitz(B2col, [B2col[1, 1]; zeros(N-1)])
       
       # Vector with initial conditions and force
       C=zeros(N, d)
       C[1, :] .= y0 + h^alpha/gamma(alpha+1)*Rfun(t0, y0, alpha, lambda)
       for i=2:N
           C[i, :] .= 1/auxgamma*y0*(a[i-1, 1]-a[i-1, 2])
       end
       
       # Diagonal matrix for fast inverse approximation
       #%%% D_delt = diag(1,delt,...,delt^(N-1))
       delt=(epsil)^(1/N)
       D=zeros(N, 1)
       for i=1:N
           D[i]=delt^(i-1)
       end
    
       # Computing Lambda = FN*D*K for fast Toeplitz inversion
       lmt=fft(D.*A)
       
       #%%%% Picard iteration algorithm
       e=ones(N, d)       # Residual vector
       F=zeros(N+1, d)    #%%%% no need for circulant embedding here
       y1=kron(ones(N, 1), y[1, :])
       countpicard = 0
       while (norm(e)>tol)
    
           Y=[y[1, :]; y1]   #%%%% initial iteration value set as initial condition
    
           for i=1:N+1
               F[i, :]=Rfun(t[i], Y[i, :], alpha, lambda)
           end
           
           # Right-hand-side vector with sparse Toeplitz Matvec and vector of
           # initial conditions
           b = B1*F[2:end, :] + B2*F[2:end, :] + C
    
           # FFT of rhs for fast inversion (Eq.4.8 from Zhou et al.)
           b1=fft(D.*b)
    
           # diag(Lambda^-1)*(FN*D)*R
           b11=b1./lmt;
    
           # Solution for current Picard iteration
           y2=real(ifft(b11)./D);
    
           # Residual vector
           e=y2-y1;
    
           # Updating solution
           y1=y2
    
           countpicard = countpicard + 1;
       end
    
       # Attributing convergec solution to output vector
       for i=2:N+1
           y[i, :]=y1[i-1, :]
       end

    end
    return t, y
end
#=
t0 = 0;              # Initial time
T  = 1;              # Final time
y0 = 1;              # Initial condition  
alpha = 0.2;         # Fractional order value
epsil = 0.5*10^-8;   # Epsilon value for fast Topelitz inversion
tol = 1e-7;          # Numerical tolerance for the Picard iteration (nonlinear cases)

Ntests = 5;          # Number of runs to be performed at different dt
lambdamrho = -1;     # linear term coefficient (lambda <= 0)
imex_type = 2;       # Type of IMEX solver (1 = first order, 2 = second order)

#% Correction terms %%

# linear term
Mlu =2;
sigma1 = [1-alpha, 1, 3 * alpha, 4 * alpha, 5 * alpha];

# force term
Mf = 2;
sigma2 = [1-alpha, 1, 3 * alpha, 4 * alpha, 5 * alpha];

# history term
Mu = 1;
sigma3 = [1, 2 - alpha, 3 * alpha, 4 * alpha, 5 * alpha];

# force term linearization
Me = 2;
sigma4 = [1-alpha, 1, 3 * alpha, 4 * alpha, 5 * alpha];

# Analytical solution powers
Nterms = 4;
sigma = [1, 2, 3, 4];

# Initializing arrays
N     = zeros(Ntests, 1);    # Number of time steps
dt    = zeros(Float64, Ntests, 1);    # Time-step size
Err   = zeros(Ntests, 1);    # Error array
pErr  = zeros(Ntests, 1);    # Convergence rate array
ctime = zeros(Ntests, 1);    # CPU time array
#% Right-hand-side definition %%
function fun_f1(t, y, alpha, lambda)

    # Fabricated RHS from defined analytical solution
    # Analytical solution powers
    Nterms = 4;
    sigma = [1, 2, 3, 4];  
    # Coefficient from nonlinear term
    rho = 0.1;
                
    # Fractional derivative of fabricated solution
    dut = 0;    
    for i=1:Nterms
        dut = dut + (gamma(sigma[i]+1)/gamma(sigma[i]+1-alpha))*t^(sigma[i]-alpha)
    end
    
    # Linear term of fabricated solution
    ut = 1;
    for i=1:Nterms
        ut = ut + t^sigma[i]
    end
    
    # Analytical RHS
    ft = dut - lambda * ut + rho*ut^2
    
    # Constructed RHS
    f = -rho*y.^2 .+ ft
    return f
end

t=Any
y11=Any
# Loop over the number of tests
for k=1:Ntests
        
    # Computing time-step size and number of time steps
    dt[k] = 2.0^(-(5+k))
    N[k]  = T/dt[k]
    
        
    (t, y11) = IMEX_I(fun_f1, t0, T, y0, N[k],alpha,lambdamrho, Mlu,sigma1,Mf,sigma2,Mu,sigma3,Me, sigma4,tol,epsil);
#=
    # Evaluating known analytical solution
    y1 = y0*ones(1, size(y11,2));
    # Power-law form of analytical solution
    for i=1:Nterms
        y1 = y1 + t.^sigma(i);
    end    
=#
end

# Printing the results
Data_Fast = [N, dt, Err, pErr, ctime]

using Plots

plot(t, y11)
=#
