function fast_conv(x, y)
    Lx = length(x)
    #Ly = size(y, 2)
    problem_size = size(y, 1)

    r = Lx
    z = zeros(Number, problem_size, r)
    X = ourfft(x, r)
    for i in 1:problem_size
        Y = ourfft(y[i, :]', r)
        Z = X .* Y
        z[i, :] = ourifft(Z, r)
    end
    return z
end

function ourfft(x, n)
    s = length(x)
    x = x[:]
    if s > n
        return fft(x[1:n])
    elseif s < n
        return fft([x; zeros(n - s)])
    else
        return fft(x)
    end
end

function ourifft(x, n)
    s = length(x)
    x = x[:]
    if s > n
        return ifft(x[1:n])
    elseif s < n
        return ifft([x; zeros(n - s)])
    else
        return ifft(x)
    end
end

function rowfft(x::AbstractMatrix, n)
    result = zeros(Complex, size(x)[1], n)
    for i in 1:size(x)[1]
        result[i, :] = ourfft(x[i, :], n)
    end
    return result
end
