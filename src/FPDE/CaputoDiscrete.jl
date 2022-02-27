function [D]=CaputoD(f,fy)
    r=f[1,1];T=f[1,2];n=f[1,3];a=f[1,4];
    p=[n,r,a];
    A=wj(p);
    e=[n,T];
    if  ischar(fy)==1
        B=fj(e,fy);
    else
        B=fy;
    end
    C = ones(1,n);
    hh=((n/T)^a)/(gamma(1-a));
    D=C*A*B*hh;
    return
    function [A]=wj(p)

    n=p(1,1);r=p(1,2);a=p(1,3);
    A=zeros(n,n+1);
    for iw=1:r-2;
    for jw=1:iw+1;
    A(iw,jw)=w([iw+1-jw,iw+1,iw,n,a]);
    end
    end
    for iw=r-1:n;
    for jw=iw-r+2:iw+1;
    A(iw,jw)=w([iw+1-jw,r,iw,n,a]);
    end
    end
    return
    function [s]=w(q)
    i=q[1,1];r=q[1,2];j=q[1,3];n=q[1,4];a=q[1,5]

    ar=ones(1,r-1);
    br=ones(1,r-1);
    for lj=1:r-2;
        jj=r-lj-1;
        kj=r-2;
        tj=i-1;
        aj=0:tj;
        bj=tj+2:kj+1;
        cj=[aj,bj];
        dj=-1:tj-1;
        ej=tj+1:kj;
        fj=[dj,ej];
        yj=nchoosek(cj,jj);
        pj=nchoosek(fj,jj);
        sj=1;
        tj=1;
        for m=1:jj;
            sj=sj.*yj(:,m);
            tj=tj.*pj(:,m);
            ar(1,lj)=sum(sj);
            br(1,lj)=sum(tj);
        end
    end
    s=0;
    for l=1:r-1;
        k=1;
        for m1=1:l;
            k=k/(m1-a);
        end
        k=k*factorial(l); 
        k=k*(ar(1,l)*((n-j)^(l-a))-br(1,l)*((n-j+1)^(l-a)));
        s=s+k;
    end
    s=s*((-1)^(i+1))/(factorial(i)*factorial(r-1-i));
    return s
end
function fj(e,fy)
    n=e(1,1);T=e(1,2);

    B=zeros(n+1,1);
    for i2=1:n+1;
        B[i2,1] = feval(fy,T*(i2-1)/n);
    end
    return B
end

function [f]=fq(x)
    f=x^8
return