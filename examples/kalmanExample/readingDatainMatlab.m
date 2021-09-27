%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reading data in matlab

clear all;
close all;
clc;

fileID = fopen('MCKResult.txt','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
systemResult = fscanf(fileID,formatSpec,sizeA);
systemResult = systemResult';

fclose(fileID);

fileID = fopen('MCKKalmanResult.txt','r');
formatSpec = '%f %f %f';
sizeA = [3 Inf];
kalmanResult = fscanf(fileID,formatSpec,sizeA);
kalmanResult = kalmanResult';

fclose(fileID);

%% new creation

x = systemResult(:,1);
x2 = systemResult(:,2);
t = systemResult(:,3);

xKalman = kalmanResult(:,1);
x2Kalman = kalmanResult(:,2);
tKalman = kalmanResult(:,3);

figure
plot(t,x)
hold on
plot(tKalman,xKalman)
legend("Pos","Filtered Pos")
xlabel("Time (sec)")
ylabel("Position Comparison")

figure
plot(t,x2)
hold on
plot(tKalman,x2Kalman)
legend("Vel","Filtered Vel")
xlabel("Time (sec)")
ylabel("Velocity Comparison")


