%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duwinddx  = ifft2(ddx.*ddz.*Cp,'symmetric');
dwwinddz  = ifft2(-ddz.*ddx.*Cp,'symmetric');
duwinddz  = ifft2(ddz.*ddz.*Cp,'symmetric');
dwwinddx = ifft2(-ddx.*ddx.*Cp,'symmetric');
theta     = ifft2(CT1,'symmetric');
dthetadz  = ifft2(ddz.*CT1,'symmetric');
dthetadx  = ifft2(ddx.*CT1,'symmetric');

%%%%%%%%%%%%%%%%%%%%%%%%

A = duwinddx - dwwinddz;
B = duwinddz + dwwinddx;

N2 = g*( theta0z + dthetadz )./( theta0 + theta );
    
Km = zeros(Nz,Nx);
for k = 1:Nz/2+1
    for j = 1:Nx
        Ri = N2(k,j)/( A(k,j)^2 + (U0z(k,j)+B(k,j))^2 );
        Km(k,j) = kmix^2*dx*dz ...
                * sqrt( A(k,j)^2 + B(k,j)^2 ) ...
                * sqrt( max( 1-3*Ri,0 ) );
    end
end

Km(1,:) = 0;
Km(Nz/2+1,:) = 0;
for k = Nz/2+2:Nz
    Km(k,:) = Km(Nz-k+2,:);
end

Kh = 3*Km;  % 3KM or 1KM

%%%%%%%%%%%%%%%%%%%%%%%%%%% in Fourior space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Du = ddx.*fft2( Km.*A ) + ddz.*fft2( Km.*B );
Dw = ddx.*fft2( Km.*B ) - ddz.*fft2( Km.*A );
Dt = ddx.*fft2( Kh.*dthetadx ) + ddz.*fft2( Kh.*dthetadz );
Dv = ddx.*Du - ddz.*Dw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
