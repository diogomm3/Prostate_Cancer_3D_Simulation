close all
%x = (8,1);
Nx = 400;  Ny = 300;  Nz = 80;
x = zeros(Nx*Ny*Nz,1);  y = x;  z = x;  c = x;

%fileID = fopen('test7310.dat','r');
fileID = fopen('tests_ur3D\test9252.dat','r');

%%%fprintf(fileID,'%d %d %d %d\n',x,y,z,c);
%[x y z c] = fscanf(fileID,'%d %d %d %d')
A = fscanf(fileID,'%d %d %d %d');
fclose(fileID);

iv = 0;
for ix=1:4:length(A)
    iv = iv + 1;
    x(iv) = A(ix);    y(iv) = A(ix+1);
    z(iv) = A(ix+2);  c(iv) = A(ix+3);
end
%scatter3(x,y,z,10,c)
A00 = zeros(Nx,Ny); A05 = zeros(Nx,Ny); A10 = zeros(Nx,Ny); A15 = zeros(Nx,Ny); 
A20 = zeros(Nx,Ny); A25 = zeros(Nx,Ny); A30 = zeros(Nx,Ny); A45 = zeros(Nx,Ny);
A40 = zeros(Nx,Ny); A45 = zeros(Nx,Ny); A50 = zeros(Nx,Ny); A55 = zeros(Nx,Ny);
A60 = zeros(Nx,Ny); A65 = zeros(Nx,Ny); A70 = zeros(Nx,Ny); A75 = zeros(Nx,Ny);

for ix=1:iv
    if z(ix)==0
        A00(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==5
        A05(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==10
        A10(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==15
        A15(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==20
        A20(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==25 
        A25(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==30
        A30(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==35
        A35(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==40
        A40(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==45
        A45(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==50
        A50(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==55
        A55(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==60
        A60(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==65
        A65(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==70
        A70(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==75
        A75(x(ix)+1,y(ix)+1) = c(ix);
    end
end

%colormap(prism(10))
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

subplot(3,3,1)
imagesc(A00')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])


subplot(3,3,2)
imagesc(A10')
set(gca,'YDir','normal')
%colormap(lines(10))
colormap(map)
caxis([0 9])
%colorbar

subplot(3,3,3)
imagesc(A20')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])

subplot(3,3,4)
imagesc(A30')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
%colorbar

subplot(3,3,5)
imagesc(A40')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])

subplot(3,3,6)
imagesc(A50')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])


subplot(3,3,7)
imagesc(A60')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])


subplot(3,3,8)
imagesc(A70')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])

%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(gcf,'test701', 'pdf')
