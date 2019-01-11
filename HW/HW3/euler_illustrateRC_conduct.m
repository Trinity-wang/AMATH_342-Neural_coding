%euler_illustrateRC_conduct.m

Vlist=zeros(1,length(tlist));

%initialize
Vlist(1)=V0;


for n=1:length(tlist)-1
    t=tlist(n);
    Vlist(n+1)=Vlist(n) + ( -Vlist(n)/(R*C) + gapplist(n)*(E-Vlist(n)) )*deltat;
end

set(gca,'FontSize',16)
subplot(211)
plot(tlist,Vlist,'-','LineWidth',2,'MarkerSize',26), hold on;
xlabel('time(ms)','Fontsize',20); ylabel('V(t)','Fontsize',20);
ylim([0 thresh + 5]);


subplot(212)
plot(tlist,gapplist,'-','LineWidth',2); hold on
xlabel('time(ms)','Fontsize',20); ylabel('gapp(t)','Fontsize',20); 
ylim([0 4*magn]);