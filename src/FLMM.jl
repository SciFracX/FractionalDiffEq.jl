function FastConv(x, y)
    Lx = length(x); Ly = size(y, 2) ; problem_size = size(y, 1)

    r = Lx
    z = zeros(problem_size, r)
    X = ourfft(x, r)
    for i = 1:problem_size
        Y = ourfft(y[i, :]', r)
        Z = X.*Y
        z[i, :] = ourifft(Z, r)
    end
    return z
end

function ourfft(x, n)
    s=length(x)
    x=x[:]
    if s > n
        return fft(x[1:n])
    elseif s < n
        return fft([x; zeros(n-s)])
    else
        return fft(x)
    end
end

function ourifft(x, n)
    s=length(x)
    x=x[:]
    if s > n
        return ifft(x[1:n])
    elseif s < n
        return ifft([x; zeros(n-s)])
    else
        return ifft(x)
    end
end

function Weights(alpha, N)
    """
switch method
    % Trapezoidal method with generating function ((1+x)/2/(1-x))^alpha
    case 1 
    """
    """
        omega1 = zeros(1, N+1); omega2 = copy(omega1)
        omega1[1] = 1 ; omega2[1] = 1
        alpha_minus_1 = alpha - 1 ; alpha_plus_1 = alpha + 1
        for n = 1 : N
            omega1[n+1] = (alpha_plus_1/n - 1)*omega1[n]
            omega2[n+1] = (1 + alpha_minus_1/n)*omega2[n]
        end
        x = fft([omega1  zeros(1, length(omega1))])
        y = fft([omega2  zeros(1, length(omega2))])
        omega = ifft(x.*y)
        omega = omega[1:N+1]/2^alpha
        """

        """
    % Newton-Gregory formula with generating function (1-x)^(-alpha)*(1-alpha/2*(1-x))
    case 2 
        omega1 = zeros(1,N+1) ; omega = omega1 ;
        alphameno1 = alpha - 1 ;
        omega1(1) = 1 ;
        for n = 1 : N
            omega1(n+1) = (1 + alphameno1/n)*omega1(n) ;
        end
        omega(1) = 1-alpha/2 ;
        omega(2:N+1) = (1-alpha/2)*omega1(2:N+1) + alpha/2*omega1(1:N) ;
     
	% BDF-2 with generating function (2/3/(1-4x/3+x^2/3))^alpha
    case 3 
    """
        omega = zeros(1, N+1)
        onethird = 1/3; fourthird = 4/3
        twothird_oneminusalpha = 2/3*(1-alpha)
        fourthird_oneminusalpha = 4/3*(1-alpha)
        omega[1] = 1 ; omega[2] = fourthird*alpha*omega[1]
        for n = 2 : N
            omega[n+1] = (fourthird - fourthird_oneminusalpha/n)*omega[n] + (twothird_oneminusalpha/n - onethird)*omega[n-1]
        end
        omega = omega.*((2/3)^(alpha))
        """   
    end
 """ 
    k = floor(1/abs(alpha))
    if abs(k - 1/alpha) < 1.0e-12
        A = collect(0:k).*abs(alpha)
    else
        A = [collect(0:k).*abs(alpha) 1]
    end
    s = length(A) - 1

    nn = collect(0:N)
    V = zeros(s+1, s+1) ; jj_nu = zeros(s+1, N+1) ; nn_nu_alpha = copy(jj_nu)
    for i = 0 : s
        nu = A[i+1]
        V[i+1,:] = collect(0:s).^nu
        jj_nu[i+1, :] = nn.^nu
        if alpha > 0
            nn_nu_alpha[i+1,:] = gamma(nu+1)/gamma(nu+1+alpha)*nn.^(nu+alpha)
        else
            if i == 0
                nn_nu_alpha[i+1,:] = zeros(1,N+1)
            else
                nn_nu_alpha[i+1,:] = gamma(nu+1)/gamma(nu+1+alpha)*nn.^(nu+alpha)
            end
        end
    end
    temp = FastConv([omega zeros(size(omega)[1])], [jj_nu zeros(size(jj_nu)[1])])
    b = nn_nu_alpha - temp[:, 1:N+1]

    w = V\b
    return omega, w, s
end