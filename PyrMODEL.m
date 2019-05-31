clear;
C=0.01; %orig 0.01
s = 30;
dt=0.01;
t=0:dt:35;
I = zeros(s,length(t));
k=1:1:length(t);




gKmax=0.015;%0.015
gNamax = .35;%orig .32
gKv3max = 0.055; %orig 0.055
gleak=0.001; %orig 0.001





init = .06;
step = .002;

Ek = -85;
ENa = 50;
Eleak = -70;
Vc = zeros(s,length(t));
Vc(:,:) = -60;
Ekb = zeros(s,length(t));
Ekb(:,:)=-85;
gleak = zeros(s,length(t)-1);
gleak(:,:)=0.001;
alphan = zeros(s,length(t));
alphan(:,:) = 76.4.*exp(0.082.*-70);

alphn = zeros(s,length(t));
alphn(:,:) = 76.4.*exp(0.082.*-70);

alphm = zeros(s,length(t));
alphm(:,:) = 1.2.*exp(0.0512.*-70);

alphh = zeros(s,length(t));
alphh(:,:) = 0.0013.*exp(-0.1016.*-70);

betan = zeros(s,length(t));
betan(:,:) = 0.0508.*exp(-0.0519.*-70);

betam = zeros(s,length(t));
betam(:,:) = 0.0383.*exp(-0.093.*-70);

betah = zeros(s,length(t));
betah(:,:) = (1.9999*1.2).*exp(0.0384.*-70);

alphk = zeros(s,length(t));
betak = zeros(s,length(t));

beta(:,:) = (-4.1476/(1+exp(((-70-6)-57.57235)/16.88916))+4.1476);
alphk(:,:) = (0.12611*exp(-(-70-6)/32.30984));








taun = zeros(s,length(t));
taun(:,:) = 1./(alphn(1,1)+betan(1,1));
taum = zeros(s,length(t));
taum(:,:) = 1./(alphm(1,1)+betam(1,1));
tauh = zeros(s,length(t));
tauh(:,:) = 1./(alphh(1,1)+betah(1,1));
tauk = zeros(s,length(t));
tauk(:,:) = 1./(alphk(1,1)+betak(1,1));




ninf=zeros(s,length(t));
minf=zeros(s,length(t));
hinf=ones(s,length(t));
kinf=zeros(s,length(t));
mtinf=zeros(s,length(t));
htinf=zeros(s,length(t));


n = zeros(s,length(t));
m = zeros(s,length(t));
h = ones(s,length(t));
k = zeros(s,length(t));
mt = zeros(s,length(t));
ht = zeros(s,length(t));
nb = zeros(s,length(t));
mb = zeros(s,length(t));
hb = ones(s,length(t));
h(:,:) = .2;
kb = zeros(s,length(t));
Ip = zeros(s,length(t));
Ids = zeros(s,length(t));
IKv3 = zeros(s,length(t));
IK = zeros(s,length(t));
INa = zeros(s,length(t));
Ileak = zeros(s,length(t));
Ileak(:,:) = zeros(s,length(t));
ICat = zeros(s,length(t));
ns=randn(1,length(t));
[b a]=butter(8,0.01);%[b a]=butter(8,0.024);
nsb=filtfilt(b,a,ns);


dndt = ones(s,length(t));
%----------------------------------------------------------------SOMA

for j = 1:1:s;

% I(j,1:length(k)) = (init + step*j)+0.075*nsb(i);

count = 1;
counts=1;
spike=1;
threshold =1;
trough=1;
peak =1;
diffwave= 1:1:(4./dt);
for i = 2:1:length(t)-1;
%     I(j,i) = (init + step*j)+0.075*sin(i*.1);

% I(j,1:length(k)) = -0.025 +i*0.000008;

      I(j,i) = (init + step*j)+(abs(0.06*nsb(i)))*-0;  



        alphm(j,i+1) = 136*exp(0.082.*(Vc(j,i)+4));
        betam(j,i+1) = 0.0383.*exp(-0.093.*(Vc(j,i)-4));
        
%         
        alphh(j,i+1) = 0.00013.*exp(-0.1016.*Vc(j,i));
        betah(j,i+1) = (1.9999.*1.2)*exp(0.0384.*Vc(j,i));
        

        taum(j,i+1) = 1./(alphm(j,i) + betam(j,i));
        minf(j,i+1) = (alphm(j,i)./(alphm(j,i)+betam(j,i)));
        m(j,i+1) = (minf(j,i)-((minf(j,i)-m(j,i)).*(exp(-(dt)./(taum(j,i))))));
        
        
        tauh(j,i+1) = (1./(alphh(j,i) + betah(j,i)))*1;
        hinf(j,i+1) = (alphh(j,i)./(alphh(j,i)+betah(j,i)));
        h(j,i+1) = (hinf(j,i)-((hinf(j,i)-h(j,i)).*(exp(-(dt)./(tauh(j,i))))));
        
         INa(j,i+1) = gNamax.*(m(j,i).^3*h(j,i)).*(Vc(j,i)-ENa); 
         

        tauk(j,i+1) = (0.4065 + (2*69.88/pi)*(30.49/(4*(Vc(j,i)+40.45)^2 + 30.49^2)));
       % (0.4065 + (2*69.88/pi)*(30.49/(4*(Vc(j,i)+40.45)^2 + 30.49^2)))+.0;
        kinf(j,i+1) = (-0.95411/(1+exp((Vc(j,i)+15.329)/10.4)) + .95217)^.25;
        k(j,i+1) = (kinf(j,i)-((kinf(j,i)-k(j,i)).*(exp(-(dt)./(tauk(j,i))))));
        
        
        
        IKv3(j,i+1) = gKv3max.*(k(j,i).^4).*(Vc(j,i)-Ek); 
        Ileak(j,i+1) = gleak(j,i).*(Vc(j,i)-Eleak); 
 
        

        
        taun(j,i+1) = (.5507 + (2*139.57/pi)*(22.73/(4*(Vc(j,i)+36.86)^2 + 22.73^2)))-.0;
        ninf(j,i+1) = (-0.9611/(1+exp((Vc(j,i)+36.46)/9.14)) + .9849)^.25;
        n(j,i+1) = (ninf(j,i)-((ninf(j,i)-n(j,i)).*(exp(-(dt)./(taun(j,i))))));
        
        IK(j,i+1) = gKmax.*(n(j,i).^4).*(Vc(j,i)-Ekb(j,i)); 
        
           dndt(j,i+1) = (n(j,i+1)-n(j,i))/dt;    
   
%    if dndt(j,i) > -0.1 && dndt(j,i-1) < -0.1;
% %             gleak(j,i:i+.5/dt)=.3; 
%            gleak(j,i:i+1.5/dt)=0.0012*exp(j/6.06)+0.0045;
% %            h(j,i+1)=0.085;
% %            ggg(j)=gleak(j,i);
% %            leakky(j,i)=(0.0012*exp(j/6.06)+0.0045)*(Vc(j,i)-Eleak);
% %             I(j,i+140)=I(j,i)-(0.0012*exp(j/6.06)+0.0045)*(Vc(j,i)-Eleak);
% %            I(j,i+1.5/dt)=I(j,i)-.125;
% % %        nn(j,counts)=n(j,i);
% % %        hh(j,counts)=h(j,i);
% % %        counts=1+counts;
% %         Ip(j,i:i+1/dt)=.05;
% % n(j,i+1)=0.56;
%    end;
% %            n(j,i+1)=0.56;
%            h(j,i+1)=0.085;
%         end;


    
   VA   = ((I(j,i) + Ip(j,i) - gNamax*m(j,i)^3*(h(j,i))*(Vc(j,i) - ENa) - gKv3max*k(j,i).^4*(Vc(j,i) - Ek) - gKmax*n(j,i).^4*(Vc(j,i) - Ek) - gleak(j,i)*(Vc(j,i) - Eleak))/C)*dt;
   VB    = ((I(j,i) + Ip(j,i)- gNamax*m(j,i)^3*(h(j,i))*((Vc(j,i)+ VA/2) - ENa) - gKv3max*k(j,i).^4*((Vc(j,i)+ VA/2) - Ek) - gKmax*n(j,i).^4*((Vc(j,i) + VA/2) - Ek) - gleak(j,i)*((Vc(j,i)+VA/2) - Eleak))/C)*dt;
   VC    = ((I(j,i) + Ip(j,i)- gNamax*m(j,i)^3*(h(j,i))*((Vc(j,i)+ VB/2) - ENa) - gKv3max*k(j,i).^4*((Vc(j,i) + VB/2) - Ek) - gKmax*n(j,i).^4*((Vc(j,i) + VB/2) - Ek) - gleak(j,i)*((Vc(j,i)+VB/2) - Eleak))/C)*dt;
   VD    = ((I(j,i) + Ip(j,i) - gNamax*m(j,i)^3*(h(j,i))*((Vc(j,i)+ VC) - ENa)- gKv3max*k(j,i).^4*((Vc(j,i) + VC) - Ek) - gKmax*n(j,i).^4*((Vc(j,i) + VC) - Ek) - gleak(j,i)*((Vc(j,i)+VC) - Eleak))/C)*dt;
   Vc(j,i+1) = Vc(j,i) + (VA + 2*VB + 2*VC + VD)./6; 
   
%                    if Vc(j,i) > -35 && Vc(j,i-1) < -35; %spike threshold,
%             spike(count) = i;
%            
%                  if (spike(count) + 2.5./dt) <= length(Vc(j));
%                      threshold(count) = Vc(j,i);
%                      start(count)=p;
%                      diffwave = Vc(start(count):start(count)+(2.5./dt));
% %                      trough(count) = min(diffwave);
% %                      peak(count) = max(diffwave);
%                      diffwave = 1:1:(2.5./dt);
%                end;
%              
%                count=count+1;
%             end; 
%  
        end;

        

        
      
 

        
        
        
      
%       for p = 1:1:length(Vc)-1;
%           dVdt = diff(Vc(j,:))./dt;
%           
%           
%           if p > 3 && dVdt(p) > 36 && dVdt(p-1) < 36 ; %spike threshold,
%             spike(count) = p;
%            
%                  if (spike(count) + 2.5./dt) <= length(dVdt);
%                      threshold(count) = Vc(j,p);
%                      start(count)=p;
%                      diffwave = Vc(start(count):start(count)+(2.5./dt));
%                      trough(count) = min(diffwave);
%                      peak(count) = max(diffwave);
%                      diffwave = 1:1:(2.5./dt);
%                end;
%              
%                count=count+1;
%             end; 
%         end;
%                 

           
 
       
    

  
    
    

      

      
 AHP(j) = sum(threshold - trough)/count;
    ISI{j} = (spike((3:end))-spike((2:end-1)))*dt;
    Freq(j) = (1./((sum(ISI{j}))/(length(ISI{j}))))*1000;


end;


figure
plot (t,Vc,'b');

figure
plot (I(:,i),Freq);
