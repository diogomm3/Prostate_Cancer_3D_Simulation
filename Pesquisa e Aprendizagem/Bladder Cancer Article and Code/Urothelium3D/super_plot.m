close all
clear all

%Nx = 400;  Ny = 300;  Nz = 80;
Nx = 400;  Ny = 200;  Nz = 80;
x = zeros(Nx*Ny*Nz,1);  y = x;  z = x;  c = x;

fileID = fopen('tests_ur3D\cells9252.dat','r');

A = fscanf(fileID,'%d %d %d %d');
fclose(fileID);

iv = 0;
for ix=1:4:length(A)
        iv = iv + 1;
        x(iv) = A(ix)+1;  y(iv) = A(ix+1)+1;  z(iv) = A(ix+2)+1;  c(iv) = A(ix+3);
        CC(x(iv), y(iv), z(iv)) = c(iv);
end

map = [0 0 0
       1 1 1
       0 1 1
       0 0 1
       1 0 1
       0.7 0.7 0.7
       0 1 0
       1 1 0
       1 0.5 0
       1 0 0];

xx = 1:Nx;  yy = 1:Ny;  zz = 1:Nz;

for ic=0:8
    %isosurface(yy, xx, zz, CC, ic)
    isosurface(zz, xx, yy, permute(CC, [1,3,2]), ic)
    colormap(map)
    caxis([0 8])
    pause(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(gcf,'super9252', 'pdf')
