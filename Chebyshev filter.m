T = 1;
N1 = 30;
p = 6;
%fs = 51200;
fs = 5120;
Ts = 1/fs;
t = T*fs +1;

Tpp = 10;
Tsp = 10;

mu = 10^-2.7;
X = rand(t,1);
%x = downsample(UniversalMil.VarName2,10);

PP = IMPULSE1([1,-.3,0.2],[1,0,0,0,0,0,0,0],0,Ts,Tpp);
SP = IMPULSE1([1, 1.5, -1],[1,0,0,0,0],0,Ts,Tpp);

PP = PP/max(PP);
SP = SP/max(SP);

Yd = zeros(t,1);                      %Recorded noise
Ys = zeros(t,1);                      %Control Signal

Cw=zeros(N1, p);                          %weights
e_cxlms=zeros(t,1);                    %error

Xhx=zeros(N1, p);

Xw =  zeros(t, 1);

for n=1:t
    
    for i=1:min(n,length(PP))
        Yd(n) = Yd(n) + PP(i)*X(n-i +1);
    end
    
    for i=1:p
        for j=1:min(n,N1)
            Xw(n) = Xw(n) + Cw(j,i)*T_func(X(n-j+1), i);
        end
    end
    
    
    for i=1:min(n,N1+1)
        Ys(n) = Ys(n) + SP(i)*Xw(n-i+1);
    end
    
    e_cxlms(n)=Yd(n)+Ys(n);

    for i=1:p
        temp =0;
        for j=1:min(n,length(SP))
            temp = temp + SP(j)*T_func(X(n-j+1), i);
        end
        Xhx(:,i)= [temp; Xhx(1:end-1,i)];
    end

    for i=1:p
        Cw(:,i) = Cw(:,i) - mu*e_cxlms(n)*Xhx(:,i);
    end
end



figure(4);
plot(e_cxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise residue')

figure(2);
plot(Yd) 
hold on 
plot(Yd-e_cxlms, 'r');
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
plot(Yd-e_cxlms, 'r')
hold on
plot(e_cxlms);
ylabel('Amplitude');
xlabel('Discrete time k');
legend('Noise signal', 'Control signal','error residual')

function sys3 = IMPULSE1(num,den,Ti,Ts,Tf)

    sys = tf(num, den, Ts);
    sys3 = impulse(sys,Ti:Ts:Tf);

end
function res = T_func(x, n)
    if n<=1
        res = x^n;
    else
        res = 2*x*T_func(x, n-1) - T_func(x, n-2);
    end
end