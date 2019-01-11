%euler method simulator

deltat=0.2 ; %timestep
Tmax=10;

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


%define input currents

%Way 1:  general test of linearity
%Iapplist1=ones(1,length(tlist));
%Iapplist2=sin(tlist);
%Iapplist3=Iapplist1+Iapplist2;

%Way 2:  successive input impulses:  watch them summate LINEARLY over time
Iapplist1=zeros(1,length(tlist));
    T1=1;
    T2=2;
    Iapplist1(find (tlist>T1 & tlist<T2) )=1;    
Iapplist2=zeros(1,length(tlist));
    T3=3;
    T4=4;
    Iapplist1(find (tlist>T3 & tlist<T4) )=1;    
Iapplist3=Iapplist1+Iapplist2;


%circuit parameters
R=2;
C=1;


for n=1:length(tlist)-1
    t=tlist(n);
    Vlist1(n+1)=Vlist1(n) + (-Vlist1(n)/(R*C) + Iapplist1(n)/C )*deltat;
    Vlist2(n+1)=Vlist2(n) + (-Vlist2(n)/(R*C) + Iapplist2(n)/C )*deltat;
    Vlist3(n+1)=Vlist3(n) + (-Vlist3(n)/(R*C) + Iapplist3(n)/C )*deltat;
end

set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

figure
set(gca,'FontSize',16)
subplot(211)
plot(tlist,Vlist1,'.-','LineWidth',2,'MarkerSize',26); hold on
plot(tlist,Vlist2,'.-','LineWidth',2,'MarkerSize',26);
plot(tlist,Vlist3,'.-','LineWidth',2,'MarkerSize',26);
plot(tlist,Vlist1+Vlist2,'o-','LineWidth',2,'MarkerSize',20);
xlabel('t','Fontsize',20); ylabel('V(t)','Fontsize',20); 
legend('V1','V2','V3','V1+V2')

subplot(212)
plot(tlist,Iapplist1,'.-','LineWidth',2,'MarkerSize',26); hold on
plot(tlist,Iapplist2,'.-','LineWidth',2,'MarkerSize',26);
plot(tlist,Iapplist3,'.-','LineWidth',2,'MarkerSize',26);
xlabel('t','Fontsize',20); ylabel('I(t)','Fontsize',20); 
legend('I1','I2','I3')