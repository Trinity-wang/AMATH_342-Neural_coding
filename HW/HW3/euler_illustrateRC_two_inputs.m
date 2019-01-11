%euler method simulator
tlist=linspace(0,Tmax,Tmax/deltat +1) ;

%to test linearity, we will solve our system 3 times!
Vlist1=zeros(1,length(tlist));
Vlist2=zeros(1,length(tlist));
Vlist3=zeros(1,length(tlist));


%initialize
V0=0;
Vlist1(1)=V0;
Vlist2(1)=V0;
Vlist3(1)=V0;

for n=1:length(tlist)-1
    t=tlist(n);
    Vlist1(n+1)=Vlist1(n) + (-Vlist1(n)/(R*C) + Iapplist1(n)/C )*deltat;
    Vlist2(n+1)=Vlist2(n) + (-Vlist2(n)/(R*C) + Iapplist2(n)/C )*deltat;
end

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

figure
set(gca,'FontSize',16)
subplot(211)
plot(tlist,Vlist1,'.-','LineWidth',2,'MarkerSize',26); hold on
plot(tlist,Vlist2,'.-','LineWidth',2,'MarkerSize',26);
xlabel('time (ms)','Fontsize',20); ylabel('V(t)','Fontsize',20); 
legend('V1','V2')

subplot(212)
plot(tlist,Iapplist1,'.-','LineWidth',2,'MarkerSize',26); hold on
plot(tlist,Iapplist2,'.-','LineWidth',2,'MarkerSize',26);
xlabel('time (ms)','Fontsize',20); ylabel('I(t)','Fontsize',20); 
legend('I1','I2')