using SpecialFunctions

h=0.01; tfinal=20; t=collect(0:h:tfinal)
N=ceil(Int, (tfinal)/h);
alpha=0.96;

k=40; m=28; l=3;
x1(t, x, y, z) = k*(y-x);
y1(t, x, y, z) = x-x.*z+m*y+sin(x);
z1(t, x, y, z) = x.*y-l*z;

x, y, z = zeros(N+1), zeros(N+1), zeros(N+1)

x[1]=0.1; y[1]=0.1; z[1]=0.1;


x[2]=x[1]+h.*x1(t[1], x[1], y[1], z[1]);
y[2]=y[1]+h.*y1(t[1], x[1], y[1], z[1]);
z[2]=z[1]+h.*z1(t[1], x[1], y[1], z[1]);

x[3]=x[2]+(h/2).*(3*x1(t[2], x[2], y[2], z[2])-x1(t[1], x[1], y[1], z[1])); 
y[3]=y[2]+(h/2).*(3*y1(t[2], x[2], y[2], z[2])-y1(t[1], x[1], y[1], z[1])); 
z[3]=z[2]+(h/2).*(3*z1(t[2], x[2], y[2], z[2])-z1(t[1], x[1], y[1], z[1])); 
for n=3:N
j=collect(3:n)

temp1 = @. x[1]+((h.^alpha)/(gamma(alpha+1)))*sum(((n+1-j)^alpha-(n-j)^alpha)*x1(t[j-2], x[j-2], y[j-2], z[j-2]))+((h.^alpha)./(gamma(alpha+2))).*sum((x1(t[j-1],x[j-1], y[j-1], z[j-1])-x1(t[j-2], x[j-2], y[j-2], z[j-2])).*((n+1-j).^alpha.*(n-j+3+2*alpha)-(n-j).^alpha.*(n-j+3+3*alpha)))+((h.^alpha)./(2*gamma(alpha+3))).*sum((x1(t[j], x[j], y[j], z[j])-2*x1(t[j-1], x[j-1], y[j-1], z[j-1]) +x1(t[j-2], x[j-2], y[j-2], z[j-2])).*((n-j+1).^alpha.*(2*(n-j).^2+(3*alpha+10).*(n-j)+2*(alpha).^2+9*alpha+12)-(n-j).^alpha.*(2*(n-j).^2+(5*alpha+10).*(n-j)+6*alpha.^2+18*alpha+12)));
x[n+1]=temp1[1]

temp2 = @. y[1]+((h.^alpha)./(gamma(alpha+1))).*sum(((n+1-j).^alpha-(n-j).^alpha).*y1(t[j-2], x[j-2], y[j-2], z[j-2]))+((h.^alpha)./(gamma(alpha+2))).*sum((y1(t[j-1],x[j-1], y[j-1], z[j-1])-y1(t[j-2], x[j-2], y[j-2], z[j-2])).*((n+1-j).^alpha.*(n-j+3+2*alpha)-(n-j).^alpha.*(n-j+3+3*alpha)))+((h.^alpha)./(2*gamma(alpha+3))).*sum((y1(t[j], x[j], y[j], z[j])-2*y1(t[j-1], x[j-1], y[j-1], z[j-1]) +y1(t[j-2], x[j-2], y[j-2], z[j-2])).*((n-j+1).^alpha.*(2*(n-j).^2+(3*alpha+10).*(n-j)+2*(alpha).^2+9*alpha+12)-(n-j).^alpha.*(2*(n-j).^2+(5*alpha+10).*(n-j)+6*alpha.^2+18*alpha+12)));
y[n+1]=temp2[1]

temp3 = @. z[1]+((h.^alpha)./(gamma(alpha+1))).*sum(((n+1-j).^alpha-(n-j).^alpha).*z1(t[j-2], x[j-2], y[j-2], z[j-2]))+((h.^alpha)./(gamma(alpha+2))).*sum((z1(t[j-1],x[j-1], y[j-1], z[j-1])-z1(t[j-2], x[j-2], y[j-2], z[j-2])).*((n+1-j).^alpha.*(n-j+3+2*alpha)-(n-j).^alpha.*(n-j+3+3*alpha)))+((h.^alpha)./(2*gamma(alpha+3))).*sum((z1(t[j], x[j], y[j], z[j])-2*z1(t[j-1], x[j-1], y[j-1], z[j-1]) +z1(t[j-2], x[j-2], y[j-2], z[j-2])).*((n-j+1).^alpha.*(2*(n-j).^2+(3*alpha+10).*(n-j)+2*(alpha).^2+9*alpha+12)-(n-j).^alpha.*(2*(n-j).^2+(5*alpha+10).*(n-j)+6*alpha.^2+18*alpha+12)));
z[n+1]=temp3[1]
end

using Plots
plot3d(x, y, z)