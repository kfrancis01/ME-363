clear
clc
close
a = 2;
b = 3;
c = 9;

m1 = a; %kg
m2 = b; %kg
k1 = 10*c; %N/m
k2 = 10*b; %N/m
k3 = 10*a; %N/m
x1_dot = -1; %m/s
x2_dot = 2; %m/s

ksum1 = (k1-k2)+k2;
ksum2 = (k2+k3)-k2;

wn1 = sqrt(ksum1/m1);
wn2 = sqrt(ksum2/m2);
wn = [wn1;wn2];

mass = [m1 0;0 m2];
spring = [(k1-k2) k2;-k2 (k2+k3)];

xo = [0;0];
vo = [-1;2];

lamda = (mass*wn.^2)./spring;
u = (-lamda*mass)+spring;
A = (u(1)^-1)*xo;
B = ([wn(1) 0;0 wn(2)]^-1)*(u^-1)*vo;

i = 1;
for t = 0:0.005:6
 x1 =u(:,1).*((A(1).*cos(wn(1).*t))+(B(1).*sin(wn(1).*t)));
 x2 = u(:,2).*((A(2).*cos(wn(2).*t))+(B(2).*sin(wn(2).*t)));
 xtot = x1+x2;
 xm1(i) = xtot(1);
 xm2(i) = xtot(2);
i = i +1;
end

y = [0:0.005:6];
plot(y,xm1)
hold on
plot(y,xm2)
xlabel('Time (s)')
ylabel('Displacement (m)')
legend('Mass One','Mass Two')