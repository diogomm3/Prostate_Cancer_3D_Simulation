close all;
clear;
clc;

change = 'Adhesion_Energy';

type = 'Tumor_Dynamics_Tumoral_AdhesionE_2_3';

Run10 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_1.txt",'r');
Run20 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_2.txt",'r');
Run30 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_3.txt",'r');
Run40 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_4.txt",'r');
Run50 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_5.txt",'r');
Run60 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_6.txt",'r');

A10 = fscanf(Run10,'%d %d %d %d %d %d');
A20 = fscanf(Run20,'%d %d %d %d %d %d');
A30 = fscanf(Run30,'%d %d %d %d %d %d');
A40 = fscanf(Run40,'%d %d %d %d %d %d');
A50 = fscanf(Run50,'%d %d %d %d %d %d');
A60 = fscanf(Run60,'%d %d %d %d %d %d');

MCS = A10(1:6:end);
Med10 = (A10(2:6:end) + A20(2:6:end) + A30(2:6:end) + A40(2:6:end) + A50(2:6:end) + A60(2:6:end))/6;
Med20 = (A10(3:6:end) + A20(3:6:end) + A30(3:6:end) + A40(3:6:end) + A50(3:6:end) + A60(3:6:end))/6;
Med30 = (A10(4:6:end) + A20(4:6:end) + A30(4:6:end) + A40(4:6:end) + A50(4:6:end) + A60(4:6:end))/6;
Med40 = (A10(5:6:end) + A20(5:6:end) + A30(5:6:end) + A40(5:6:end) + A50(5:6:end) + A60(5:6:end))/6;
Med50 = (A10(6:6:end) + A20(6:6:end) + A30(6:6:end) + A40(6:6:end) + A50(6:6:end) + A60(6:6:end))/6;
std10 = std([A10(2:6:end), A20(2:6:end), A30(2:6:end), A40(2:6:end), A50(2:6:end), A60(2:6:end)], [], 2)/sqrt(6);
std20 = std([A10(3:6:end), A20(3:6:end), A30(3:6:end), A40(3:6:end), A50(3:6:end), A60(3:6:end)], [], 2)/sqrt(6);
std30 = std([A10(4:6:end), A20(4:6:end), A30(4:6:end), A40(4:6:end), A50(4:6:end), A60(4:6:end)], [], 2)/sqrt(6);
std40 = std([A10(5:6:end), A20(5:6:end), A30(5:6:end), A40(5:6:end), A50(5:6:end), A60(5:6:end)], [], 2)/sqrt(6);
std50 = std([A10(6:6:end), A20(6:6:end), A30(6:6:end), A40(6:6:end), A50(6:6:end), A60(6:6:end)], [], 2)/sqrt(6);


type = 'Tumor_Dynamics_Normal_Luminal';

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


type = 'Tumor_Dynamics_Tumoral_AdhesionE_2_8';

Run12 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_2.txt",'r');
Run22 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_3.txt",'r');
Run32 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_7.txt",'r');
Run42 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_8.txt",'r');
Run52 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_9.txt",'r');
Run62 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_10.txt",'r');

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


type = 'Tumor_Dynamics_Tumoral_AdhesionE_3_06';

Run13 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_1.txt",'r');
Run23 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_2.txt",'r');
Run33 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_3.txt",'r');
Run43 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_4.txt",'r');
Run53 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_5.txt",'r');
Run63 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_6.txt",'r');

A13 = fscanf(Run13,'%d %d %d %d %d %d');
A23 = fscanf(Run23,'%d %d %d %d %d %d');
A33 = fscanf(Run33,'%d %d %d %d %d %d');
A43 = fscanf(Run43,'%d %d %d %d %d %d');
A53 = fscanf(Run53,'%d %d %d %d %d %d');
A63 = fscanf(Run63,'%d %d %d %d %d %d');

Med13 = (A13(2:6:end) + A23(2:6:end) + A33(2:6:end) + A43(2:6:end) + A53(2:6:end) + A63(2:6:end))/6;
Med23 = (A13(3:6:end) + A23(3:6:end) + A33(3:6:end) + A43(3:6:end) + A53(3:6:end) + A63(3:6:end))/6;
Med33 = (A13(4:6:end) + A23(4:6:end) + A33(4:6:end) + A43(4:6:end) + A53(4:6:end) + A63(4:6:end))/6;
Med43 = (A13(5:6:end) + A23(5:6:end) + A33(5:6:end) + A43(5:6:end) + A53(5:6:end) + A63(5:6:end))/6;
Med53 = (A13(6:6:end) + A23(6:6:end) + A33(6:6:end) + A43(6:6:end) + A53(6:6:end) + A63(6:6:end))/6;
std13 = std([A13(2:6:end), A23(2:6:end), A33(2:6:end), A43(2:6:end), A53(2:6:end), A63(2:6:end)], [], 2)/sqrt(6);
std23 = std([A13(3:6:end), A23(3:6:end), A33(3:6:end), A43(3:6:end), A53(3:6:end), A63(3:6:end)], [], 2)/sqrt(6);
std33 = std([A13(4:6:end), A23(4:6:end), A33(4:6:end), A43(4:6:end), A53(4:6:end), A63(4:6:end)], [], 2)/sqrt(6);
std43 = std([A13(5:6:end), A23(5:6:end), A33(5:6:end), A43(5:6:end), A53(5:6:end), A63(5:6:end)], [], 2)/sqrt(6);
std53 = std([A13(6:6:end), A23(6:6:end), A33(6:6:end), A43(6:6:end), A53(6:6:end), A63(6:6:end)], [], 2)/sqrt(6);

figure(1)
% subplot(3,3,1)
% mseb(MCS,[Med10';Med11';Med12';Med13'],[std10';std11';std12';std13'])
% legend('AdhesionE = 2.3','AdhesionE = 2.55 (std)','AdhesionE = 2.8','AdhesionE = 3.06')
% xlabel('MCS')
% ylabel('CELLS')
% title('TUMORAL CELL NUMBER')
% legend('Location', 'northwest')

% subplot(3,3,2)
% mseb(MCS,[Med20';Med21';Med22';Med23'],[std20';std21';std22';std23'])
% legend('AdhesionE = 2.3','AdhesionE = 2.55 (std)','AdhesionE = 2.8','AdhesionE = 3.06')
% xlabel('MCS')
% ylabel('CELLS')
% title('LUMINAL CELL NUMBER')
% legend('Location', 'northeast')

% subplot(3,3,3)
% mseb(MCS,[Med30';Med31';Med32';Med33'],[std30';std31';std32';std33'])
% legend('AdhesionE = 2.3','AdhesionE = 2.55 (std)','AdhesionE = 2.8','AdhesionE = 3.06')
% xlabel('MCS')
% ylabel('CELLS')
% title('BASAL CELL NUMBER')
% legend('Location', 'northwest')
% 
% subplot(3,3,4)
mseb(MCS,[Med40';Med41';Med42';Med43'],[std40';std41';std42';std43'])
legend('AdhesionE = 2.3','AdhesionE = 2.55 (std)','AdhesionE = 2.8','AdhesionE = 3.06')
xlabel('MCS')
ylabel('VOLUME (VOXELS)')
title('TUMOR VOLUME')
legend('Location', 'northwest')

% subplot(3,3,5)
% mseb(MCS,[Med50';Med51';Med52';Med53'],[std50';std51';std52';std53'])
% legend('AdhesionE = 2.3','AdhesionE = 2.55 (std)','AdhesionE = 2.8','AdhesionE = 3.06')
% xlabel('MCS')
% ylabel('AREA')
% title('TUMOR AREA')
% legend('Location', 'northwest')
% 
% set(gcf, 'PaperOrientation', 'landscape');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 11, 7]);
% saveas(gcf,"Simulation Graphs Results\Graph_Results_"+change,'pdf')