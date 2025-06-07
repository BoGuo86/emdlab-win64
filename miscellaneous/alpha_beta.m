function alpha_beta()

x0 = [1e-2,1e-4,1.5,2];

% loss at 50Hz
p50 = [0
    0.03
    0.07
    0.15
    0.25
    0.37
    0.49
    0.63
    0.79
    0.96
    1.14
    1.35
    1.58
    1.83
    2.25
    2.62
    2.95
    3.21
    3.46];

% loss at 100Hz
p100 = [0
    0.04
    0.17
    0.36
    0.60
    0.88
    1.21
    1.57
    1.98
    2.43
    2.93
    3.49
    4.12
    4.86
    5.78
    6.83];

% loss at 200Hz
p200 = [0
    0.10
    0.44
    0.93
    1.55
    2.31
    3.21
    4.24
    5.41
    6.75
    8.25
    9.94
    11.8
    14.0
    16.6
    19.5];

% loss at 400Hz
p400 = [0
    0.30
    1.18
    2.51
    4.25
    6.40
    9.01
    12.1
    15.7
    20.0
    24.9
    30.6
    37.1
    44.5
    53.0
    62.5];

% loss at 1000Hz
p1000 = [0
    1.41
    5.10
    10.4
    17.7
    26.5
    37.9
    52.1
    69.5
    90.8
    116
    147];

% loss at 2500Hz
p2500 = [0
    5.79
    20.3
    42.6
    74.1
    117
    174
    249
    346
    470
    629];

subplot(321)
b50 = 0:0.1:0.1*length(p50)-0.1;
plot(b50,p50,'*');title('f = 50 Hz');xlabel('B [T]');ylabel('Loss [W/kg]')
subplot(322)
b100 = 0:0.1:0.1*length(p100)-0.1;
plot(b100,p100,'*');title('f = 100 Hz');xlabel('B [T]');ylabel('Loss [W/kg]')
subplot(323)
b200 = 0:0.1:0.1*length(p200)-0.1;
plot(b200,p200,'*');title('f = 200 Hz');xlabel('B [T]');ylabel('Loss [W/kg]')
subplot(324)
b400 = 0:0.1:0.1*length(p400)-0.1;
plot(b400,p400,'*');title('f = 400 Hz');xlabel('B [T]');ylabel('Loss [W/kg]')
subplot(325)
b1000 = 0:0.1:0.1*length(p1000)-0.1;
plot(b1000,p1000,'*');title('f = 1000 Hz');xlabel('B [T]');ylabel('Loss [W/kg]')
subplot(326)
b2500 = 0:0.1:0.1*length(p2500)-0.1;
plot(b2500,p2500,'*');title('f = 2500 Hz');xlabel('B [T]');ylabel('Loss [W/kg]')

[k,res] = fmincon(@evalObjective,x0,[],[],[],[],[0 0 0 1.5],[1 1 2 2]);

    function y = evalObjective(x)
        f = 50;
        y = sum(abs(p50-x(1)*f^x(3)*b50'.^x(4)-x(2)*f^2*b50'.^2));
        f = 100;
        y = y + sum(abs(p100-x(1)*f^x(3)*b100'.^x(4)-x(2)*f^2*b100'.^2));
        f = 200;
        y = y + sum(abs(p200-x(1)*f^x(3)*b200'.^x(4)-x(2)*f^2*b200'.^2));
        f = 400;
        y = y + sum(abs(p400-x(1)*f^x(3)*b400'.^x(4)-x(2)*f^2*b400'.^2));
        f = 1000;
        y = y + sum(abs(p1000-x(1)*f^x(3)*b1000'.^x(4)-x(2)*f^2*b1000'.^2));
        f = 2500;
        y = y + sum(abs(p2500-x(1)*f^x(3)*b2500'.^x(4)-x(2)*f^2*b2500'.^2));
    end

hold on
b = 0:0.1:1.8;
subplot(321);hold on
f = 50;
plot(b,k(1)*f^k(3)*b.^k(4)+k(2)*f^2*b.^2)
subplot(322);hold on
f = 100;
plot(b,k(1)*f^k(3)*b.^k(4)+k(2)*f^2*b.^2)
subplot(323);hold on
f = 200;
plot(b,k(1)*f^k(3)*b.^k(4)+k(2)*f^2*b.^2)
subplot(324);hold on
f = 400;
plot(b,k(1)*f^k(3)*b.^k(4)+k(2)*f^2*b.^2)
subplot(325);hold on
f = 1000;
plot(b,k(1)*f^k(3)*b.^k(4)+k(2)*f^2*b.^2)
subplot(326);hold on
f = 2500;
plot(b,k(1)*f^k(3)*b.^k(4)+k(2)*f^2*b.^2)

clc
disp(k)
disp(res)

end