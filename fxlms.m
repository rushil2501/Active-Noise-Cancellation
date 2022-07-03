T = 1;
N = 30;
fs = 5120;
Ts = 1/fs;
t = T*fs +1;
Tpp = 10;
Tsp = 10;
mu = 10^-2.5;
%X = sin(2*pi*500*(0:Ts:T));
X = rand(1, T*fs +1);
PP = IMPULSE1([1,-.3,0.2],[1,0,0,0,0,0,0,0],0,Ts,Tpp);
SP = IMPULSE1([1, 1.5, -1],[1,0,0,0,0],0,Ts,Tpp);
PP = PP/max(PP);
SP = SP/max(SP);

Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t,1);                      %Control Signal
e_fxlms = zeros(t,1);                     %error
Cw1 =  zeros(1, N);
Xw1 =  zeros(1, N);
Cw_sum = zeros(length(SP), 1);

del = 0.9;

for n=1:t
    for i=1:min(n, length(PP))
        Yd(n) = Yd(n) + PP(i)*X(n-i +1);
    end
    Cy = 0;
    for i=1:min(n,N)
        Cy = Cy + Cw1(i)*X(n-i+1);
    end
    Cw_sum=[Cy; Cw_sum(1: end-1)];
    
    Ys(n) = sum(Cw_sum.*SP);
    e_fxlms(n)=Yd(n)+Ys(n);
    
    temp = 0;
    for i=1:min(n, N)
        temp = temp + SP(i)*X(n-i+1);
    end             
    Xw1=[temp Xw1(1:end-1)];
    Cw1 = Cw1 - mu*e_fxlms(n)*Xw1;
    disp(n);
end

figure(1);
plot(e_fxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise residue')

figure(2);
plot(Yd) 
hold on 
plot(Yd-e_fxlms, 'r');
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal')
hold off

figure(5);
plot(Yd) 
hold on 
plot(Yd-e_fxlms, 'r')
hold on
plot(e_fxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal','errror residual')
hold off



function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)
    sys = tf(num, den, Ts);
    sys3 = impulse(sys,Ti:Ts:Tf);
end