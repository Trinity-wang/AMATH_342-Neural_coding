% Hodgkin/Huxley Equations (V_HH = -V-65), with current definition 
% of membrane potential (V=Vin-Vout)


vna=50;  %set the constants
vk=-77;
vl=-54.4;
gna=120;
gk=36;
gl=.3;
c=1;

v_init=-65;  %the initial conditions
m_init=.052;
h_init=.596;
n_init=.317;

m=zeros(npoints,1); %initialize everything to zero
n=zeros(npoints,1);
h=zeros(npoints,1);
v=zeros(npoints,1);
time=zeros(npoints,1);

m(1)=m_init; %set the initial conditions to be the first entry in the vectors
n(1)=n_init;
h(1)=h_init;
v(1)=v_init;
time(1)=0.0;

% changed I to be a vector rather than a scalar.
for step=1:npoints-1
    v(step+1)=v(step)+((I(step) - gna*h(step)*(v(step)-vna)*m(step)^3 ... 
               -gk*(v(step)-vk)*n(step)^4-gl*(v(step)-vl))/c)*dt;
    m(step+1)=m(step)+ (alpha_m(v(step))*(1-m(step))-beta_m(v(step))*m(step))*dt;
    h(step+1)=h(step)+ (alpha_h(v(step))*(1-h(step))-beta_h(v(step))*h(step))*dt;
    n(step+1)=n(step)+ (alpha_n(v(step))*(1-n(step))-beta_n(v(step))*n(step))*dt;
    time(step+1)=time(step)+dt;
end


set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20); 

figure
subplot(211);
plot(time,v, 'r'), hold on;
plot(time, thresh*ones(npoints,1), 'g');
xlabel('t');
ylabel('V');
ylim([-80,80]);

subplot(212);
plot(time, I, 'b'), hold on;
plot(time, ones(1, length(time))* Imagn, 'g');
