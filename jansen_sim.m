%% 1 node JR model with known parameters

dt=0.001; % sampling time
F=[10,20,40,80,100,130,160]; % stimulation frequency
stim_time=400;   % stimulation train duration
Tbase=5000;      % Samples before stimulation starts.
Ttot=8000;       % Total number of samples simulated

% the first row saves 1 node model simulation, 2nd and 3rd row saves 2 node
% model simulation

sim_lfp=cell(3,length(F));

% Model parameters 
v0=6;
vm=5;
r=0.3;
C=135;
C1=C;
C2=0.8*C;
C3=0.25*C;
C4=0.25*C;

A=3.25;
B=22;
a=100;
b=50;
ka=1;kA=1;

% Simulating 1 node JR model using difference equations
for ff=1:length(F)
    stim_freq=F(ff);
    cycle=[1,-1,zeros(1,floor(1000/stim_freq)-2)];
    stim=[];
    for ss=1:(floor(stim_time*stim_freq/1000))
        stim=[stim,cycle];
    end
    I=[zeros(1,Tbase),stim,zeros(1,Ttot-Tbase-length(stim))];
    Ip=60*I;
    Ii=60*r*I;
    
    L=length(I);
    
    y=zeros(6,L);
    p = 0.1.*randn(1,L);
    
    for ii=2:length(y)
        y(1,ii) = y(1,ii-1)+y(4,ii-1)*dt;
        y(2,ii) = y(2,ii-1)+y(5,ii-1)*dt;
        y(3,ii) = y(3,ii-1)+y(6,ii-1)*dt;
        y(4,ii) = y(4,ii-1)+ dt* (A*a*(Ii(ii-1)        +(vm/(1+exp(r*(v0-y(2,ii-1)+y(3,ii-1))))))   -(2*a*y(4,ii-1))-(a^2*y(1,ii-1)));
        y(5,ii) = y(5,ii-1)+ dt* (kA*ka*A*a*(p(1,ii-1)+Ip(ii-1)+(C2*vm/(1+exp(r*(v0-C1*y(1,ii-1))))))       -(2*ka*a*y(5,ii-1))-(ka*ka*a^2*y(2,ii-1)));
        y(6,ii) = y(6,ii-1)+ dt* (B*b*(Ip(ii-1)         +(C4*vm/(1+exp(r*(v0-C3*y(1,ii-1))))))        -(2*b*y(6,ii-1))-(b^2*y(3,ii-1)));
    end
    
    outp=y(2,Tbase-1000+1:end)-y(3,Tbase-1000+1:end); % We ignore the first 1000 samples to let the model settle down from the initial values
    sim_lfp{1,ff}=outp-mean(outp);
end

%% 2 node JR model with known parameters

% Nodes 1 and 2 parameters
C11=C;
C12=0.8*C11;
C13=0.25*C11;
C14=0.25*C11;

C21=C;
C22=0.8*C21;
C23=0.25*C21;
C24=0.25*C21;


A1=A;
B1=B;
a1=a;
b1=b;
ka1=ka;
kA1=kA;

A2=5;
B2=26;
a2=50;
b2=10;
ka2=0.2;
kA2=0.4;

P1=50;  
P2=100;

% Inter node coupling parameters
ad=10;
K1=1000; K2=800;

% Simulating 2 node JR model using difference equations
for ff=1:length(F)
    stim_freq=F(ff);
    cycle=[1,-1,zeros(1,floor(1000/stim_freq)-2)];
    stim=[];
    for ss=1:(floor(stim_time*stim_freq/1000))
        stim=[stim,cycle];
    end
    I=[zeros(1,Tbase),stim,zeros(1,Ttot-Tbase-length(stim))];
    Ip=[60*I;0*I];
    Ii=r.*Ip;
    Is=Ip(1,:).*0.4;
    L=length(I);
    
    y0=zeros(1,L);
    y1=y0; y2=y0; y3=y0; y4=y0; y5=y0;
    y6=y0;y7=y0;y8=y0;y9=y0;
    y10=y0;y11=y0;y12=y0;y13=y0;y14=y0; y15=y0;
    
    p(1,:) = P1+0.1.*randn(1,L);
    p(2,:) = P2+0.1.*randn(1,L);
    
    for ii=2:length(y0)
        % Node 1
        y0(ii) = y0(ii-1)+y3(ii-1)*dt;
        y3(ii) = y3(ii-1)+ dt* (A1*a1*(Ii(1,ii-1)+(vm/(1+exp(r*(v0-y1(ii-1)-y13(ii-1)+y2(ii-1))))))   -(2*a1*y3(ii-1))-(a1^2*y0(ii-1)));
        
        y1(ii) = y1(ii-1)+ y4(ii-1)*dt;
        y4(ii) = y4(ii-1)+ dt* (kA1*A1*ka1*a1*(p(1,ii-1)+Ip(1,ii-1)+(C12*vm/(1+exp(r*(v0-C11*y0(ii-1))))))       -(2*ka1*a1*y4(ii-1))-(ka1^2*a1^2*y1(ii-1)));
        y13(ii)=y13(ii-1)+dt*y15(ii-1);
        y15(ii)=y15(ii-1)+dt*(A2*ad*K2*(vm/(1+exp(r*(v0-(y7(ii-1)+y12(ii-1)-y8(ii-1))))))-2*ad*y15(ii-1)-ad*ad*y13(ii-1));
        
        y2(ii) = y2(ii-1)+y5(ii-1)*dt;
        y5(ii) = y5(ii-1)+ dt* (B1*b1*(Ip(1,ii-1)+(C14*vm/(1+exp(r*(v0-C13*y0(ii-1))))))        -(2*b1*y5(ii-1))-(b1^2*y2(ii-1)));
        
        
        % Node 2
        y6(ii) = y6(ii-1)+y9(ii-1)*dt;
        y9(ii) = y9(ii-1)+ dt* (A2*a2*(Ii(2,ii-1)    +(vm/(1+exp(r*(v0-y7(ii-1)-y12(ii-1)+y8(ii-1))))))   -(2*a2*y9(ii-1))-(a2^2*y6(ii-1)));
        
        y7(ii) = y7(ii-1)+y10(ii-1)*dt;
        
        y10(ii) = y10(ii-1)+ dt* (kA2*A2*ka2*a2*(p(2,ii-1)+Ip(2,ii-1)+(K1*(y12(ii-1)))+(C22*vm/(1+exp(r*(v0-C21*y6(ii-1)))))) -(2*a2*ka2*y10(ii-1))-((a2*ka2)^2*y7(ii-1)));
        y12(ii)=y12(ii-1)+dt*y14(ii-1);
        y14(ii)=y14(ii-1)+dt*(A1*ad*(vm/(1+exp(r*(v0-(y1(ii-1)-y2(ii-1))))))-2*ad*y14(ii-1)-ad*ad*y12(ii-1));
        
        y8(ii) = y8(ii-1)+y11(ii-1)*dt;
        y11(ii) = y11(ii-1)+ dt* (B2*b2*(Ip(2,ii-1)        +(C24*vm/(1+exp(r*(v0-C23*y6(ii-1))))))        -(2*b2*y11(ii-1))-(b2^2*y8(ii-1)));
        
    end
    outp1=y1(Tbase-1000+1:end)+y13(Tbase-1000+1:end)-y2(Tbase-1000+1:end);
    outp2=y7(Tbase-1000+1:end)+y12(Tbase-1000+1:end)-y8(Tbase-1000+1:end);
    sim_lfp{2,ff}=outp1-mean(outp1);
    sim_lfp{3,ff}=outp2-mean(outp2);
end

