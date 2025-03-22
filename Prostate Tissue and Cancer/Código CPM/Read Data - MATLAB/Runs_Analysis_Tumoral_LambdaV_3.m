close all;
clear;
clc;

change = 'Tumoral_LambdaV_3';

type = 'Tumor_Dynamics_Tumoral_LambdaV_x09';

Run11 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_1.txt",'r');
Run21 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_2.txt",'r');
Run31 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_3.txt",'r');
Run41 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_4.txt",'r');
Run51 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_5.txt",'r');
Run61 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_6.txt",'r');

A11 = fscanf(Run11,'%d %d %d %d %d %d');
A21 = fscanf(Run21,'%d %d %d %d %d %d');
A31 = fscanf(Run31,'%d %d %d %d %d %d');
A41 = fscanf(Run41,'%d %d %d %d %d %d');
A51 = fscanf(Run51,'%d %d %d %d %d %d');
A61 = fscanf(Run61,'%d %d %d %d %d %d');

MCS = A11(1:6:end);
Med11 = (A11(2:6:end) + A21(2:6:end) + A31(2:6:end) + A41(2:6:end) + A51(2:6:end) + A61(2:6:end))/6;
Med21 = (A11(3:6:end) + A21(3:6:end) + A31(3:6:end) + A41(3:6:end) + A51(3:6:end) + A61(3:6:end))/6;
Med31 = (A11(4:6:end) + A21(4:6:end) + A31(4:6:end) + A41(4:6:end) + A51(4:6:end) + A61(4:6:end))/6;
Med41 = (A11(5:6:end) + A21(5:6:end) + A31(5:6:end) + A41(5:6:end) + A51(5:6:end) + A61(5:6:end))/6;
Med51 = (A11(6:6:end) + A21(6:6:end) + A31(6:6:end) + A41(6:6:end) + A51(6:6:end) + A61(6:6:end))/6;
std11 = std([A11(2:6:end), A21(2:6:end), A31(2:6:end), A41(2:6:end), A51(2:6:end), A61(2:6:end)], [], 2)/sqrt(6);
std21 = std([A11(3:6:end), A21(3:6:end), A31(3:6:end), A41(3:6:end), A51(3:6:end), A61(3:6:end)], [], 2)/sqrt(6);
std31 = std([A11(4:6:end), A21(4:6:end), A31(4:6:end), A41(4:6:end), A51(4:6:end), A61(4:6:end)], [], 2)/sqrt(6);
std41 = std([A11(5:6:end), A21(5:6:end), A31(5:6:end), A41(5:6:end), A51(5:6:end), A61(5:6:end)], [], 2)/sqrt(6);
std51 = std([A11(6:6:end), A21(6:6:end), A31(6:6:end), A41(6:6:end), A51(6:6:end), A61(6:6:end)], [], 2)/sqrt(6);

type = 'Tumor_Dynamics_Normal_Luminal';

Run12 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_1.txt",'r');
Run22 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_2.txt",'r');
Run32 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_3.txt",'r');
Run42 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_4.txt",'r');
Run52 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_5.txt",'r');
Run62 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_6.txt",'r');

A12 = fscanf(Run12,'%d %d %d %d %d %d');
A22 = fscanf(Run22,'%d %d %d %d %d %d');
A32 = fscanf(Run32,'%d %d %d %d %d %d');
A42 = fscanf(Run42,'%d %d %d %d %d %d');
A52 = fscanf(Run52,'%d %d %d %d %d %d');
A62 = fscanf(Run62,'%d %d %d %d %d %d');

Med12 = (A12(2:6:end) + A22(2:6:end) + A32(2:6:end) + A42(2:6:end) + A52(2:6:end) + A62(2:6:end))/6;
Med22 = (A12(3:6:end) + A22(3:6:end) + A32(3:6:end) + A42(3:6:end) + A52(3:6:end) + A62(3:6:end))/6;
Med32 = (A12(4:6:end) + A22(4:6:end) + A32(4:6:end) + A42(4:6:end) + A52(4:6:end) + A62(4:6:end))/6;
Med42 = (A12(5:6:end) + A22(5:6:end) + A32(5:6:end) + A42(5:6:end) + A52(5:6:end) + A62(5:6:end))/6;
Med52 = (A12(6:6:end) + A22(6:6:end) + A32(6:6:end) + A42(6:6:end) + A52(6:6:end) + A62(6:6:end))/6;
std12 = std([A12(2:6:end), A22(2:6:end), A32(2:6:end), A42(2:6:end), A52(2:6:end), A62(2:6:end)], [], 2)/sqrt(6);
std22 = std([A12(3:6:end), A22(3:6:end), A32(3:6:end), A42(3:6:end), A52(3:6:end), A62(3:6:end)], [], 2)/sqrt(6);
std32 = std([A12(4:6:end), A22(4:6:end), A32(4:6:end), A42(4:6:end), A52(4:6:end), A62(4:6:end)], [], 2)/sqrt(6);
std42 = std([A12(5:6:end), A22(5:6:end), A32(5:6:end), A42(5:6:end), A52(5:6:end), A62(5:6:end)], [], 2)/sqrt(6);
std52 = std([A12(6:6:end), A22(6:6:end), A32(6:6:end), A42(6:6:end), A52(6:6:end), A62(6:6:end)], [], 2)/sqrt(6);

figure(1)
subplot(3,3,1)
mseb(MCS,[Med11';Med12'],[std11';std12'])
legend('lambdaV = 7920','lambdaV = 8800 (std)')
xlabel('MCS')
ylabel('CELLS')
title('TUMORAL CELL NUMBER')
legend('Location', 'northwest')

subplot(3,3,2)
mseb(MCS,[Med21';Med22'],[std21';std22'])
legend('lambdaV = 7920','lambdaV = 8800 (std)')
xlabel('MCS')
ylabel('CELLS')
title('LUMINAL CELL NUMBER')
legend('Location', 'northwest')

subplot(3,3,3)
mseb(MCS,[Med31';Med32'],[std31';std32'])
legend('lambdaV = 7920','lambdaV = 8800 (std)')
xlabel('MCS')
ylabel('CELLS')
title('BASAL CELL NUMBER')
legend('Location', 'northwest')

subplot(3,3,4)
mseb(MCS,[Med41';Med42'],[std41';std42'])
legend('lambdaV = 7920','lambdaV = 8800 (std)')
xlabel('MCS')
ylabel('VOLUME')
title('TUMOR VOLUME')
legend('Location', 'northwest')

subplot(3,3,5)
mseb(MCS,[Med51';Med52'],[std51';std52'])
legend('lambdaV = 7920','lambdaV = 8800 (std)')
xlabel('MCS')
ylabel('AREA')
title('TUMOR AREA')
legend('Location', 'northwest')

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 11, 7]);
saveas(gcf,"Simulation Graphs Results\Graph_Results_"+change,'pdf')