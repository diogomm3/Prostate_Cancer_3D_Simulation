s1 = 'tests_ur3D\test925';  s1BM = 'testBM_925';
Nx = 400;  Ny = 400;  Nz = 16;

for ifi=2 % 8

    close all

x = zeros(Nx*Ny*Nz,1);  y = x;  z = x;  c = x;

if ifi==1
    s2 = '1.dat';
elseif ifi==2
    s2 = '2.dat';
elseif ifi==3
    s2 = '3.dat';
elseif ifi==4
    s2 = '4.dat';
elseif ifi==5
    s2 = '5.dat';
elseif ifi==6
    s2 = '6.dat';
elseif ifi==7
    s2 = '7.dat';
elseif ifi==8
    s2 = '8.dat';
end

    s = strcat(s1,s2);
    fileID = fopen(s,'r');
    A = fscanf(fileID,'%d %d %d %d');
    fclose(fileID);

iv = 0;
for ix=1:4:length(A)
    iv = iv + 1;
    x(iv) = A(ix);    y(iv) = A(ix+1);
    z(iv) = A(ix+2);  c(iv) = A(ix+3);
end
%scatter3(x,y,z,10,c)
BM46 = zeros(Nx,Nz); BM47 = zeros(Nx,Nz); BM48 = zeros(Nx,Nz); BM49 = zeros(Nx,Nz); 
BM50 = zeros(Nx,Nz); BM51 = zeros(Nx,Nz); BM52 = zeros(Nx,Nz); BM53 = zeros(Nx,Nz); 
BM45 = zeros(Nx,Nz); BM54 = zeros(Nx,Nz); 

for ix=1:iv
    if y(ix)==45
        BM45(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==46
        BM46(x(ix)+1,z(ix)+1) = c(ix);    
    elseif y(ix)==47
        BM47(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==48
        BM48(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==49
        BM49(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==50
        BM50(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==51
        BM51(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==52
        BM52(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==53
        BM53(x(ix)+1,z(ix)+1) = c(ix);
    elseif y(ix)==54
        BM54(x(ix)+1,z(ix)+1) = c(ix);
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

subplot(5,2,1)
imagesc(BM45')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

subplot(5,2,2)
imagesc(BM46')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

subplot(5,2,3)
imagesc(BM47')
set(gca,'YDir','normal')
%colormap(lines(10))
colormap(map)
caxis([0 9])
%colorbar
daspect([1 1 1])

subplot(5,2,4)
imagesc(BM48')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

subplot(5,2,5)
imagesc(BM49')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
%colorbar
daspect([1 1 1])

subplot(5,2,6)
imagesc(BM50')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

subplot(5,2,7)
imagesc(BM51')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

subplot(5,2,8)
imagesc(BM52')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

subplot(5,2,9)
imagesc(BM53')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

subplot(5,2,10)
imagesc(BM54')
set(gca,'YDir','normal')
colormap(map)
caxis([0 9])
daspect([1 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%

if ifi==1
    s2 = '1';
elseif ifi==2
    s2 = '2';
elseif ifi==3
    s2 = '3';
elseif ifi==4
    s2 = '4';
elseif ifi==5
    s2 = '5';
elseif ifi==6
    s2 = '6';
elseif ifi==7
    s2 = '7';
elseif ifi==8
    s2 = '8';
end
    sBM = strcat(s1BM,s2);
    saveas(gcf,'Simulation Results\testBM_9252', 'pdf')
end
