%euler method simulator

deltat=0.1 ; %timestep
Tmax=10;

tlist=linspace(0,Tmax,Tmax/deltat +1) ;% all our time values
Vlist=zeros(1,length(tlist));

%initialize
V0= -1;
Vlist(1)=V0;


for n=1:length(tlist)-1
    t=tlist(n);
    Vlist(n+1)= Vlist(n) + (sin(t))*deltat;
end

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

figure
set(gca,'FontSize',16)
plot(tlist,Vlist,'.-','LineWidth',2,'MarkerSize',26); hold on
xlabel('t','Fontsize',20); ylabel('V(t)','Fontsize',20); 
plot(linspace(0,Tmax,1000), -cos( linspace(0,Tmax,1000) ))
legend('Euler Approx','Exact Solution')


%%

%euler method simulator

deltat=0.1 ; %timestep
Tmax=10;

tlist=linspace(0,Tmax,Tmax/deltat +1) ;% all our time values
Vlist=zeros(1,length(tlist));

%initialize
V0= 1;
Vlist(1)=V0;

R = 1;
C = 1;


for n=1:length(tlist)-1
    t=tlist(n);
    Vlist(n+1)= Vlist(n) + (-1*Vlist(n)/ (R*C) + sin(t) / C)*deltat;
end

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

figure
set(gca,'FontSize',16)
plot(tlist,Vlist,'.-','LineWidth',2,'MarkerSize',26); hold on
xlabel('t','Fontsize',20); ylabel('V(t)','Fontsize',20); 
% plot(linspace(0,Tmax,1000), -cos( linspace(0,Tmax,1000) ))
% legend('Euler Approx','Exact Solution')


%%

%euler method simulator

deltat=0.1 ; %timestep
Tmax=10;

tlist=linspace(0,Tmax,Tmax/deltat +1) ;% all our time values
Vlist=zeros(1,length(tlist));

%initialize
V0= 1;
Vlist(1)=V0;


for n=1:length(tlist)-1
    t=tlist(n);
    Vlist(n+1)= Vlist(n) + (sin(t))*deltat;
end

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

figure
set(gca,'FontSize',16)
plot(tlist,Vlist,'.-','LineWidth',2,'MarkerSize',26); hold on
xlabel('t','Fontsize',20); ylabel('V(t)','Fontsize',20); 
plot(linspace(0,Tmax,1000), -cos( linspace(0,Tmax,1000) ))
legend('Euler Approx','Exact Solution')



