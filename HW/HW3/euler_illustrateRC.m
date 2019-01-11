% euler_illustrateRC.m

% NOT MY CODE
tlist=linspace(0,Tmax,Tmax/deltat +1) ;
Vlist=zeros(1,length(tlist));

%initialize

Vlist(1)=V0;

for n=1:length(tlist)-1
    t=tlist(n);
    Vlist(n+1)=Vlist(n) + (-Vlist(n)/(R*C) + Iapplist(n)/C )*deltat;
end

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

set(gca,'FontSize',16)
subplot(211)
plot(tlist,Vlist,'.-','LineWidth',2,'MarkerSize',26); hold on
xlabel('time (ms)','Fontsize',20); ylabel('V(t) (mV)','Fontsize',20);
ylim([0, thresh + 2]);

subplot(212)
plot(tlist,Iapplist,'o-','LineWidth',2); hold on
xlabel('time (ms)','Fontsize',20); ylabel('Iapp(t)','Fontsize',20); 
ylim([0 20]);
hold on;