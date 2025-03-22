close all
%x = (8,1);
Nx = 400; Ny = 400;  Nz = 60;

x = zeros(Nx*Ny*Nz,1);  y = x;  z = x;  c = x;

fileID = fopen('Simulation_Data/Test_TumorDynamics9.dat','r');
fileEVO = fopen('Simulation_Data/Cell_Number_EVO_9.txt','r');

A = fscanf(fileID,'%d %d %d %d');
fclose(fileID);

B = fscanf(fileEVO, '%d %d %d %d');
fclose(fileEVO);

iv = 0;
for ix=1:4:length(A)
    iv = iv + 1;
    x(iv) = A(ix);    y(iv) = A(ix+1);
    z(iv) = A(ix+2);  c(iv) = A(ix+3);
end

%scatter3(x,y,z,10,c)
A00 = zeros(Nx,Ny); A01 = zeros(Nx,Ny); A02 = zeros(Nx,Ny); A03 = zeros(Nx,Ny); A04 = zeros(Nx,Ny);
A10 = zeros(Nx,Nz); A11 = zeros(Nx,Nz); A12 = zeros(Nx,Nz); A13 = zeros(Nx,Nz); A14 = zeros(Nx,Nz);
A20 = zeros(Ny,Nz); A21 = zeros(Ny,Nz); A22 = zeros(Ny,Nz); A23 = zeros(Ny,Nz); A24 = zeros(Ny,Nz);

for ix=1:iv
    if y(ix)==50
        A10(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==125
        A11(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==200
        A12(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==275
        A13(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==350
        A14(x(ix)+1,z(ix)+1) = c(ix);
    end
end
for ix=1:iv
    if x(ix)==50
        A20(y(ix)+1,z(ix)+1) = c(ix);
    elseif x(ix)==125
        A21(y(ix)+1,z(ix)+1) = c(ix);
    elseif x(ix)==200
        A22(y(ix)+1,z(ix)+1) = c(ix);
    elseif x(ix)==275
        A23(y(ix)+1,z(ix)+1) = c(ix);
    elseif x(ix)==350
        A24(y(ix)+1,z(ix)+1) = c(ix);
    end
end
for ix=1:iv
    if z(ix)==0
        A00(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==15
        A01(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==30
        A02(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==45
        A03(x(ix)+1,y(ix)+1) = c(ix);
    elseif z(ix)==60
        A04(x(ix)+1,y(ix)+1) = c(ix);
    end
end

B0 = zeros(length(B)/4,1);
B1 = zeros(length(B)/4,1);
B2 = zeros(length(B)/4,1);
B3 = zeros(length(B)/4,1);
j = 0;
for i=1:4:length(B)
    j = j + 1;
    B0(j) = B(i);
    B1(j) = B(i+1);
    B2(j) = B(i+2);
    B3(j) = B(i+3);
end

C0 = zeros(Nx,Ny,Nz);

l=1;
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            C0(i,j,k) = c(l);
            l=l+1;
        end
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

figure(1)
subplot(4,5,1)
imagesc(A00')
set(gca,'YDir','normal')
colormap(gca,map)
caxis([0 9])
% axis([0 400 0 400]);
xlabel('X')
ylabel('Y')
title('Z = 0')


subplot(4,5,2)
imagesc(A01')
set(gca,'YDir','normal')
%colormap(lines(10))
colormap(map)
caxis([0 9])
%colorbar
title('Z = 15')
xlabel('X')
ylabel('Y')

subplot(4,5,3)
imagesc(A02')
set(gca,'YDir','normal')
%colormap(lines(10))
colormap(map)
caxis([0 9])
%colorbar
title('Z = 30')
xlabel('X')
ylabel('Y')

subplot(4,5,4)
imagesc(A03')
set(gca,'YDir','normal')
%colormap(lines(10))
colormap(map)
caxis([0 9])
%colorbar
title('Z = 45')

subplot(4,5,5)
imagesc(A04')
set(gca,'YDir','normal')
%colormap(lines(10))
colormap(map)
caxis([0 9])
%colorbar
title('Z = 60')

xlabel('X')
ylabel('Y')

subplot(4,5,6)
imagesc(A10')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('Y = 50')
xlabel('X')
ylabel('Z')

subplot(4,5,7)
imagesc(A11')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
% colorbar
title('Y = 125')
xlabel('X')
ylabel('Z')

subplot(4,5,8)
imagesc(A12')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('Y = 200')
xlabel('X')
ylabel('Z')

subplot(4,5,9)
imagesc(A13')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('Y = 275')
xlabel('X')
ylabel('Z')

subplot(4,5,10)
imagesc(A14')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('Y = 350')
xlabel('X')
ylabel('Z')

subplot(4,5,11)
imagesc(A20')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('X = 50')
xlabel('Y')
ylabel('Z')

subplot(4,5,12)
imagesc(A21')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('X = 125')
xlabel('Y')
ylabel('Z')

subplot(4,5,13)
imagesc(A22')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('X = 200')
xlabel('Y')
ylabel('Z')

subplot(4,5,14)
imagesc(A23')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('X = 275')
xlabel('Y')
ylabel('Z')

subplot(4,5,15)
imagesc(A24')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
title('X = 350')
xlabel('Y')
ylabel('Z')

subplot(4,5,17)
plot(B0,B1);
title('TUMORAL CELL NUMBER EVOLUTION')
xlabel('MCS')
ylabel('CELLS NUMBER')

subplot(4,5,18)
plot(B0,B2);
title('LUMINAL CELL NUMBER EVOLUTION')
xlabel('MCS')
ylabel('CELLS NUMBER') 

subplot(4,5,19)
plot(B0,B3);
title('BASAL CELL NUMBER EVOLUTION')
xlabel('MCS')
ylabel('CELLS NUMBER') 

% minC0 = min(C0(:));
% maxC0 = max(C0(:));
% fprintf('Minimum value in C0: %f\n', minC0);
% fprintf('Maximum value in C0: %f\n', maxC0);

figure(2)
[xx,yy,zz] = meshgrid(1:Nx,1:Ny,1:Nz);
s = isosurface(xx, yy, zz, C0, 4);
p = patch(s);
isonormals(xx, yy, zz, C0, p)
view(3);
set(p,'FaceColor',[0.5 1 0.5]);  
set(p,'EdgeColor','none');
camlight;
lighting gouraud;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('TUMOR SURFACE')

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 11, 7]);
saveas(gcf,'Simulation Graphs Results\Tumor_Dynamics9_Slides', 'pdf')

figure(2)
saveas(gcf,'Simulation Graphs Results\Tumor_Dynamics9_3D_View', 'pdf')