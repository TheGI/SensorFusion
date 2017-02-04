addpath('libQuaternion');
clc;
clear all;
close all;
pkg load signal;

% Load Sensor Data
filename = "Log_1424138500562.csv";

temp = real(csvread(filename));
trim = 2:length(temp);
w = temp(trim,4:6);
m = temp(trim,7:9);
a = temp(trim,1:3);
t = temp(trim, 13);

dt = zeros(length(trim),1);
v = zeros(length(trim),3);
p = zeros(length(trim),3);
Ea = zeros(length(trim),3);
q = [ones(length(trim),1),zeros(length(trim),3)];
Ie = [0, 0, 0];

a_mag = sqrt(a(:,1).*a(:,1) + a(:,2).*a(:,2) + a(:,3).*a(:,3));
filtCutOff = 0.001;
[bc, ac] = butter(1, (2*filtCutOff)/(256), 'high');
acc_magFilt = filtfilt(bc, ac, a_mag);
acc_magFilt = abs(acc_magFilt);
filtCutOff = 5;
[bc, ac] = butter(1, (2*filtCutOff)/(256), 'low');
acc_magFilt = filtfilt(bc, ac, acc_magFilt);
stationary = acc_magFilt < 0.58;

for i = 1:2000
  Kp = 1;
  Ki = 0;
  [q(1,:) Ie] = CPFilter(q(1,:), mean(a(1:30,:)), [0, 0, 0], m(1,:), 0.01, Kp, Ki, Ie);
end

for i = 1:length(t)
  if(stationary(i))
    Kp = 0.5;
    Ki = 0;
  else
    Kp = 0;
    Ki = 0;
  end
  if i > 1
    dt(i) = (t(i) - t(i-1))/1e9;
    [q(i,:) Ie] = CPFilter(q(i-1,:), a(i,:), w(i,:), m(i,:), dt(i), Kp, Ki, Ie);
  end
  
  Ea(i,:) = quaternRotate(a(i,:), q(i,:));
%  Ea(i,:) = Ea(i,:) *9.81;
%  Ea(i,3) = Ea(i,3) - 9.81;


  if i > 1
    v(i,:) = v(i-1,:) + Ea(i,:)*dt(i);
    if(stationary(i) == 1)
      v(i,:) = [0 0 0];
    end
  end
end

% Compute integral drift during non-stationary periods
velDrift = zeros(size(v));
stationaryStart = find([0; diff(stationary)] == -1);
stationaryEnd = find([0; diff(stationary)] == 1);
for i = 1:numel(stationaryEnd)
    driftRate = v(stationaryEnd(i)-1, :) / (stationaryEnd(i) - stationaryStart(i));
    enum = 1:(stationaryEnd(i) - stationaryStart(i));
    drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
    velDrift(stationaryStart(i):stationaryEnd(i)-1, :) = drift;
end


v = v - velDrift;

for i = 2:length(v)
  v_abs(i) = norm(v(i,:));
  p(i,:) = p(i-1,:) + v(i,:) * dt(i);
end

[b_coef, a_coef] = butter(1, (2*0.1)/100, 'high');
p(:,3) = filtfilt(b_coef, a_coef, p(:,3));

figure();
title("3D Swing Trajectory and Velocity (km/h)");
xlabel("x (m)");
ylabel("y (m)");
zlabel("z (m)");
grid on;

v_abs = v_abs*3.6;
surface('XData',[p(:,1) p(:,1)],'YData',[p(:,2) p(:,2)],
        'ZData',[p(:,3) p(:,3)],'CData',[v_abs v_abs],
        'FaceColor','none','EdgeColor','flat','Marker','none',
        'LineWidth', 2.5);
colorbar;