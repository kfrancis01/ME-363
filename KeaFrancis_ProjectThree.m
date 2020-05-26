%%ME363
clear all
clc

%Given
p = 2;
q = 3;
r = 9;

m1 = p/1000 ; %kg
m2 = q/1000; %kg
m3 = r/1000; %kg

k = q/150; %N/m
c = p/1200; %
w = 3*r; % rad/s

Zmax = 13.3/100; %mm to m
Fmax = 0.41/100; %mN to N
Y = r/100; %mm to m

Z0 = 0;
Z_dot0 = 0;

damp1 = c/(2*sqrt(m1*k)); %damp ratio
damp2 = c/(2*sqrt(m2*k)); %damp ratio
damp3 = c/(2*sqrt(m3*k)); %damp ratio

wn1 = sqrt(k/m1); %natural frequency
wn2 = sqrt(k/m2); %natural frequency
wn3 = sqrt(k/m3); %natural frequency

r1 = w/wn1;
r2 = w/wn2;
r3 = w/wn3;

wd1 = wn1*sqrt(1-damp1^2);
wd2 = wn2*sqrt(1-damp2^2);
wd3 = wn3*sqrt(1-damp3^2);

%Find the total response and forces for each mass
tmax = 2.5; %s
i = 1;
t = 0;
while t < tmax
    phi1 = atand((2*damp1*r1)/(1-r1^2));
    phi2 = atand((2*damp2*r2)/(1-r2^2));
    phi3 = atand((2*damp3*r3)/(1-r3^2));
    Zp1(i) = ((Y*r1^2))*(1/sqrt(((1-r1^2)^2)+((2*damp1*r1)^2))*(cos((w*t)-phi1)));
    Zp2(i) = ((Y*r2^2))*(1/sqrt(((1-r2^2)^2)+((2*damp2*r2)^2))*(cos((w*t)-phi2)));
    Zp3(i) = ((Y*r3^2))*(1/sqrt(((1-r2^2)^2)+((2*damp3*r3)^2))*(cos((w*t)-phi3)));
    
    A1 = Z0 - Zp1(1);
    A2 = Z0 - Zp2(1);
    A3 = Z0 - Zp3(1);
    
    B1 = Z_dot0 - Zp1(1)+(damp1*wn1*A1);
    B2 = Z_dot0 - Zp2(1)+(damp2*wn2*A2);
    B3 = Z_dot0 - Zp3(1)+(damp3*wn3*A3);
    
    Zh1(i) = (exp(-damp1*wn1*t))*((A1*cos(wd1*t))+((B1/wd1)*sin(wd1*t)));
    Zh2(i) = (exp(-damp2*wn2*t))*((A2*cos(wd2*t))+((B2/wd2)*sin(wd2*t)));
    Zh3(i) = (exp(-damp3*wn3*t))*((A3*cos(wd3*t))+((B3/wd3)*sin(wd3*t)));
    
    Z1 = Zp1 + Zh1; %m
    Z2 = Zp2 + Zh2; %m
    Z3 = Zp3 + Zh3; %m
    
    Z_dot1 = diff(Z1);
    Z_dot2 = diff(Z2);
    Z_dot3 = diff(Z3);
    
    F1 = (k*Z1)+(c*[Z_dot1,zeros(1,1)]); %N
    F2 = (k*Z2)+(c*[Z_dot2,zeros(1,1)]); %N
    F3 = (k*Z3)+(c*[Z_dot3,zeros(1,1)]); %N
    
    t = t + 0.005;
    i = i + 1;
end

%If the Z value exceeds the constraint the value will be displayed in the
%command window
ii = 1;
t = 0;
while t < tmax
    if Z1(ii) > Zmax
        disp("Z1 = " + Z1(ii))
        t = tmax;
    elseif Z2(ii) > Zmax
        disp("Z2 = " + Z2(ii))
        t = tmax;
    elseif Z3(ii) > Zmax
        disp("Z3 = " + Z3(ii))
        t = tmax;
    else
        t = t + 0.005;
        ii = ii + 1;
    end
end

%Just in case the first test doesn't work this will allow you to see the
%Zmax compared to the maximum Z values for each mass
if ii < i
    disp('Check Z values')
    Zmax
    disp("Z1max = " + max(Z1))
    disp("Z2max = " + max(Z2))
    disp("Z3max = " + max(Z3))
end

%Display F that exceeds the constraint
iii = 1;
t = 0;
while t < tmax
    if F1(iii) > Fmax
        disp("F1 = " + F1(iii))
        t = tmax;
    elseif F2(iii) > Fmax
        disp("F2 = " + F2(iii))
        t = tmax;
    elseif F3(iii) > Fmax
        disp("F3 = " + F3(iii))
        t = tmax;
    else
        t = t + 0.005;
        iii = iii + 1;
    end
end

%Display the F maximum value compared to the maximum forces for each mass
if iii < i
    disp('Check F values')
    Fmax
    disp("F1max = " + max(F1))
    disp("F2max = " + max(F2))
    disp("F3max = " + max(F3))
end


%a)
%Energy
E1 = pi()*c*w.*Z1.^2; %Joules
E2 = pi()*c*w.*Z2.^2; %Joules
E3 = pi()*c*w.*Z3.^2; %Joules
disp("E1 = " + max(E1))
disp("E2 = " + max(E2))
disp("E3 = " + max(E3))

%b)
plot([0:0.05:25],Z1)
hold on
plot([0:0.05:25],Z2)
% hold on
% plot(Z3)
hold on
plot(Zmax*ones(1,i)) %create line for Zmax 
legend('Z1','Z2','Z3','Zmax')
grid ON
grid MINOR

ylabel('Displacement (m)');
xlabel('Time (s)');
