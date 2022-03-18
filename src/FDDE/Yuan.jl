struct DelayABMYuan <: FractionalDiffEqAlgorithm end

function solve(f, α, τ, T, h, ::DelayABMYuan)
  N = round(Int, T/h)
  Ndelay = round(Int, τ/h)
  q1=q2=q3=α

  x1=zeros(Ndelay+N+1)
  y1=zeros(Ndelay+N+1)
  z1=zeros(Ndelay+N+1)

  x1[Ndelay+N+1]=0
  y1[Ndelay+N+1]=0
  z1[Ndelay+N+1]=0

  x=zeros(Ndelay+N+1)
  y=zeros(Ndelay+N+1)
  z=zeros(Ndelay+N+1)

  x[Ndelay+N+1]=0
  y[Ndelay+N+1]=0
  z[Ndelay+N+1]=0

  for i=1:Ndelay
      x[i]=0.2
      y[i]=0
      z[i]=0.5
  end

  x0=copy(x[Ndelay])
  y0=copy(y[Ndelay])
  z0=copy(z[Ndelay])

  x1[Ndelay+1]=x0+h^q1*(f(x0, y0, z0, x[1], y[1], z[1], 1))/(gamma(q1)*q1)
  y1[Ndelay+1]=y0+h^q2*(f(x0, y0, z0, x[1], y[1], z[1], 2))/(gamma(q2)*q2)
  z1[Ndelay+1]=z0+h^q3*(f(x0, y0, z0, x[1], y[1], z[1], 3))/(gamma(q3)*q3)

  x[Ndelay+1]=x0+h^q1*(f(x1[Ndelay+1], y1[Ndelay+1], z1[Ndelay+1], x[2], y[2], z[2], 1)+q1*f(x0, y0, z0, x[1], y[1], z[1], 1))/gamma(q1+2)
  y[Ndelay+1]=y0+h^q2*(f(x1[Ndelay+1], y1[Ndelay+1], z1[Ndelay+1], x[2], y[2], z[2], 2)+q2*f(x0, y0, z0, x[1], y[1], z[1], 2))/gamma(q2+2)
  z[Ndelay+1]=z0+h^q3*(f(x1[Ndelay+1], y1[Ndelay+1], z1[Ndelay+1], x[2], y[2], z[2], 3)+q3*f(x0, y0, z0, x[1], y[1], z[1], 3))/gamma(q3+2)

  for n=1:N 

    M1=(n^(q1+1)-(n-q1)*(n+1)^q1)*(f(x0, y0, z0, x[1], y[1], z[1], 1))
    M2=(n^(q2+1)-(n-q2)*(n+1)^q2)*(f(x0, y0, z0, x[1], y[1], z[1], 2))
    M3=(n^(q3+1)-(n-q3)*(n+1)^q3)*(f(x0, y0, z0, x[1], y[1], z[1], 3))

    N1=((n+1)^q1-n^q1)*(f(x0, y0, z0, x[1], y[1], z[1], 1))
    N2=((n+1)^q2-n^q2)*(f(x0, y0, z0, x[1], y[1], z[1], 2))
    N3=((n+1)^q3-n^q3)*(f(x0, y0, z0, x[1], y[1], z[1], 3))

    for j=1:n  
      M1=M1+((n-j+2)^(q1+1)+(n-j)^(q1+1)-2*(n-j+1)^(q1+1))*(f(x[Ndelay+j], y[Ndelay+j], z[Ndelay+j], x[j], y[j], z[j], 1))
      M2=M2+((n-j+2)^(q2+1)+(n-j)^(q2+1)-2*(n-j+1)^(q2+1))*(f(x[Ndelay+j], y[Ndelay+j], z[Ndelay+j], x[j], y[j], z[j], 2))
      M3=M3+((n-j+2)^(q3+1)+(n-j)^(q3+1)-2*(n-j+1)^(q3+1))*(f(x[Ndelay+j], y[Ndelay+j], z[Ndelay+j], x[j], y[j], z[j], 3))

      N1=N1+((n-j+1)^q1-(n-j)^q1)*(f(x[Ndelay+j], y[Ndelay+j], z[Ndelay+j], x[j], y[j], z[j], 1))
      N2=N2+((n-j+1)^q2-(n-j)^q2)*(f(x[Ndelay+j], y[Ndelay+j], z[Ndelay+j], x[j], y[j], z[j], 2))
      N3=N3+((n-j+1)^q3-(n-j)^q3)*(f(x[Ndelay+j], y[Ndelay+j], z[Ndelay+j], x[j], y[j], z[j], 3))
    end

  x1[Ndelay+n+1]=x0+h^q1*N1/(gamma(q1)*q1)
  y1[Ndelay+n+1]=y0+h^q2*N2/(gamma(q2)*q2)
  z1[Ndelay+n+1]=z0+h^q3*N3/(gamma(q3)*q3)

  x[Ndelay+n+1]=x0+h^q1*(f(x[Ndelay+n+1], y[Ndelay+n+1], z[Ndelay+n+1], x[n+1], y[n+1], z[n+1], 1)+M1)/gamma(q1+2)
  y[Ndelay+n+1]=y0+h^q2*(f(x[Ndelay+n+1], y[Ndelay+n+1], z[Ndelay+n+1], x[n+1], y[n+1], z[n+1], 2)+M2)/gamma(q2+2)
  z[Ndelay+n+1]=z0+h^q3*(f(x[Ndelay+n+1], y[Ndelay+n+1], z[Ndelay+n+1], x[n+1], y[n+1], z[n+1], 3)+M3)/gamma(q3+2)
  end
  return x, y, z
end