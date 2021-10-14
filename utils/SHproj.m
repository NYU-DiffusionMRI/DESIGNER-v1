function [s] = SHproj(coef, dirs, lmax)

    ndirs =  size(dirs,1);
    ncoef = (lmax+1)*(lmax+2)/2;

    for i = 1:ndirs;
        dirs(i, :) = dirs(i, :) / norm(dirs(i, :));
    end

    el = acos(dirs(:,3));
    az = atan2(dirs(:,2), dirs(:,1));

    Y = zeros(ndirs, ncoef);

    for N=0:2:lmax
                   
        q = legendre(N,cos(el))';
        for m = 1:N
            if mod(m, 2)==0
                q(:,m+1) = q(:,m+1)*(-1)^m * sqrt(2*factorial(N-m)/factorial(N+m));
            else
                q(:,m+1) = -q(:,m+1)*(-1)^m * sqrt(2*factorial(N-m)/factorial(N+m));
            end
        end  
        q = q.*sqrt((2*N+1)/(4*pi));
     
       
        loff = N*(N+1)/2 + 1;
        for m=0:N
            if m
                Y(:,loff+m) = q(:,m+1).*cos(m*az);
                Y(:,loff-m) = q(:,m+1).*sin(m*az);

            else
                Y(:,loff) = q(:,1);
            end
        end
    end

    s = Y*coef;
end