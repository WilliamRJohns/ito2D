clear ddx ddz inv_laplacian
for k = 1:Nz

    for j = 1:Nx/2+1
      ddx(k,j) = i*2*pi*(j-1)/(xF-xS);
    end

    for j = Nx/2+2:Nx
      ddx(k,j) = -ddx(k,Nx-j+2);
    end

end

for j = 1:Nx

    for k = 1:Nz/2+1
      ddz(k,j) = i*2*pi*(k-1)/(zF-zS);
    end

    for k = Nz/2+2:Nz
      ddz(k,j) = -ddz(Nz-k+2,j);
    end
    
end

laplacian = ddx.^2 + ddz.^2;
laplacian4 = ddx.^4 + ddz.^4;

for k = 1:Nz
    for j = 1:Nx
        if ddx(k,j) == 0 && ddz(k,j) == 0
         inv_laplacian(k,j) = 0;
        else
         inv_laplacian(k,j) = 1./laplacian(k,j);
        end
    end
end

