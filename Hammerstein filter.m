T = 10;
N1 = 30;
p = 6;
fs = 5120;
Ts = 1/fs;
t = T*fs+1;

Tpp = 10;
Tsp = 10;

mu = 10^-1;
%X = downsample(filling.VarName2,10);
%X = rand(t,1);
X = downsample(exp1.VarName2,10);
X = X/max(X);

PP = IMPULSE1([1,-.3,0.2],[1,0,0,0,0,0,0,0],0,Ts,Tpp);
SP = IMPULSE1([1, 1.5, -1],[1,0,0,0,0],0,Ts,Tpp);
PP = PP(1:40);
SP = SP(1:40);

PP = PP/max(PP);
SP = SP/max(SP);

Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t, 1);                      %Control Signal
e_hxlms = zeros(t,1);                    %error

Cw=zeros(N1, p);                          %weights
Xw = zeros(N1, p);
 
Cw_sum = zeros(length(SP), 1);

tic;
for n=1:t
    
    for i=1:min(n,length(PP))
        Yd(n) = Yd(n) + PP(i)*X(n-i +1);
    end
    
    Cy = 0;
    for i=1:p
        for j=1:min(n, N1)
               Cy = Cy + Cw(j,i)*(X(n-j+1)^i);
        end
    end
    
    
    Cw_sum=[Cy; Cw_sum(1: end-1)];             
    Ys(n) = sum(Cw_sum.*SP);
    
    e_hxlms(n)=Yd(n)+Ys(n);

    for i=1:p
        temp =0;
        for j=1:min(n,length(SP))
            temp = temp + SP(j)*(X(n-j+1)^i);
        end
        Xw(:,i)= [temp; Xw(1:end-1,i)];
    end
    
    for i=1:p
        Cw(:,i) = Cw(:,i) - mu*e_hxlms(n)*Xw(:,i);
    end
end
toc;

figure(4);
plot(e_hxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise residue')

figure(2);
plot(Yd) 
hold on 
plot(Yd-e_hxlms, 'r');
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal')

figure(3);
plot(Cw);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('weights')

figure(5);
plot(Yd) 
hold on 
plot(Yd-e_hxlms, 'r')
hold on
plot(e_hxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal','errror residual')



function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den, Ts);
    
    sys3 = impulse(sys,Ti:Ts:Tf);

end
