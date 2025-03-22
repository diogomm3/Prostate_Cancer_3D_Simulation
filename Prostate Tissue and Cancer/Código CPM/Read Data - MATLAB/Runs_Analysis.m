close all;
clear;
clc;

type = 'Tumor_Dynamics_Normal_Luminal';

Run1 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_1.txt",'r');
Run2 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_2.txt",'r');
Run3 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_3.txt",'r');
Run4 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_4.txt",'r');
Run5 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_5.txt",'r');
Run6 = fopen("Simulation_Data/"+type+"/SecondGeomLarge_Cell_Number_EVO_6.txt",'r');

A1 = fscanf(Run1,'%d %d %d %d %d %d');
A2 = fscanf(Run2,'%d %d %d %d %d %d');
A3 = fscanf(Run3,'%d %d %d %d %d %d');
A4 = fscanf(Run4,'%d %d %d %d %d %d');
A5 = fscanf(Run5,'%d %d %d %d %d %d');
A6 = fscanf(Run6,'%d %d %d %d %d %d');

Med1 = zeros(length(A1)/6,1);
Med2 = zeros(length(A1)/6,1);
Med3 = zeros(length(A1)/6,1);
Med4 = zeros(length(A1)/6,1);
Med5 = zeros(length(A1)/6,1);
MCS = zeros(length(A1)/6,1);
minus_1 = zeros(length(A1)/6,1);
plus_1 = zeros(length(A1)/6,1);
minus_2 = zeros(length(A1)/6,1);
plus_2 = zeros(length(A1)/6,1);
minus_3 = zeros(length(A1)/6,1);
plus_3 = zeros(length(A1)/6,1);
minus_4 = zeros(length(A1)/6,1);
plus_4 = zeros(length(A1)/6,1);
minus_5 = zeros(length(A1)/6,1);
plus_5 = zeros(length(A1)/6,1);
j = 0;
for i=1:6:length(A1)
    j = j + 1;
    MCS(j) = A1(i);
    Med1(j) = (A1(i+1)+A2(i+1)+A3(i+1)+A4(i+1)+A5(i+1)+A6(i+1))/6;
    Med2(j) = (A1(i+2)+A2(i+2)+A3(i+2)+A4(i+2)+A5(i+2)+A6(i+2))/6;
    Med3(j) = (A1(i+3)+A2(i+3)+A3(i+3)+A4(i+3)+A5(i+3)+A6(i+3))/6;
    Med4(j) = (A1(i+4)+A2(i+4)+A3(i+4)+A4(i+4)+A5(i+4)+A6(i+4))/6;
    Med5(j) = (A1(i+5)+A2(i+5)+A3(i+5)+A4(i+5)+A5(i+5)+A6(i+5))/6;
    std1 = std([A1(i+1),A2(i+1),A3(i+1),A4(i+1),A5(i+1),A6(i+1)])/sqrt(6);
    std2 = std([A1(i+2),A2(i+2),A3(i+2),A4(i+2),A5(i+2),A6(i+2)])/sqrt(6);
    std3 = std([A1(i+3),A2(i+3),A3(i+3),A4(i+3),A5(i+3),A6(i+3)])/sqrt(6);
    std4 = std([A1(i+4),A2(i+4),A3(i+4),A4(i+4),A5(i+4),A6(i+4)])/sqrt(6);
    std5 = std([A1(i+5),A2(i+5),A3(i+5),A4(i+5),A5(i+5),A6(i+5)])/sqrt(6);
    minus_1(j) = Med1(j)-std1;
    plus_1(j) = Med1(j)+std1;
    minus_2(j) = Med2(j)-std2;
    plus_2(j) = Med2(j)+std2;
    minus_3(j) = Med3(j)-std3;
    plus_3(j) = Med3(j)+std3;
    minus_4(j) = Med4(j)-std4;
    plus_4(j) = Med4(j)+std4;
    minus_5(j) = Med5(j)-std5;
    plus_5(j) = Med5(j)+std5;
end

figure(1)
% subplot(3,3,1)
% plot(MCS,Med1,'-k','LineWidth',1);
% hold on
% plot(MCS,plus_1,MCS,minus_1,'r')
% legend('RunMed','Std Deviation')
% xlabel('MCS')
% ylabel('CELLS')
% title('TUMORAL CELL NUMBER')
% legend('Location', 'northwest')
% 
% subplot(3,3,2)
% plot(MCS,Med2,'-k','LineWidth',1);
% hold on
% plot(MCS,plus_2,MCS,minus_2,'r')
% legend('RunMed','Std Deviation')
% xlabel('MCS')
% ylabel('CELLS')
% title('LUMINAL CELL NUMBER')
% legend('Location', 'northeast')
% 
% subplot(3,3,3)
% plot(MCS,Med3,'-k','LineWidth',1);
% hold on
% plot(MCS,plus_3,MCS,minus_3,'r')
% legend('RunMed','Std Deviation')
% xlabel('MCS')
% ylabel('CELLS')
% title('BASAL CELL NUMBER')
% legend('Location', 'northeast')
% 
% subplot(3,3,4)
plot(MCS,Med4,'-k','LineWidth',1);
hold on
plot(MCS,plus_4,MCS,minus_4,'r')
legend('RunMed','Std Deviation')
xlabel('MCS')
ylabel('VOLUME (VOXELS)')
title('TUMOR VOLUME')
legend('Location', 'northwest')
% 
% subplot(3,3,5)
% plot(MCS,Med5,'-k','LineWidth',1);
% hold on
% plot(MCS,plus_5,MCS,minus_5,'r')
% legend('RunMed','Std Deviation')
% xlabel('MCS')
% ylabel('AREA')
% title('TUMOR AREA')
% legend('Location', 'northwest')
% 
% set(gcf, 'PaperOrientation', 'landscape');
% set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 11, 7]);
% saveas(gcf,"Simulation Graphs Results\"+type+"\Graph_Results",'pdf')