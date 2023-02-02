clear all;
thk = [0.500];
vs = [2.600 3.200];
vp=0.9409+2.0947*vs-0.8206*vs.^2+0.2683*vs.^3-0.0251*vs.^4;
dns=1.6612*vp-0.4721*vp.^2+0.0671*vp.^3-0.0043*vp.^4+0.000106*vp.^5;

thk0 = [0.5];
vs0 = [3.200 3.200];
vp0=0.9409+2.0947*vs0-0.8206*vs0.^2+0.2683*vs0.^3-0.0251*vs0.^4;
dns0=1.6612*vp0-0.4721*vp0.^2+0.0671*vp0.^3-0.0043*vp0.^4+0.000106*vp0.^5;

freq = logspace(-1,1,100);
%
u0 = site_resp_surf(thk0,dns0,vp0,vs0,freq,'R');
u = site_resp_surf(thk,dns,vp,vs,freq,'R');
%
rayp=0.00001;
ub0 = site_resp_body(thk0,dns0,vp0,vs0,freq,rayp,'SV');
ub = site_resp_body(thk,dns,vp,vs,freq,rayp,'SV');

%%
figure('position',[100,100,600,300]);
semilogx(freq,abs(u(:,2))./abs(u0(:,2)),'m--','linewidth',2);
hold on;
semilogx(freq,abs(u(:,1))./abs(u0(:,1)),'r:','linewidth',2);
semilogx(freq,abs(ub(:,1))./abs(ub0(:,1)),'c','linewidth',2);
ylim([0.8,2]);
grid on; box on;
legend('Rayleigh Vertical', 'Rayleigh Horizontal', 'Vertical S (Horizontal)');
xlabel('Frequency (Hz)');
ylabel('Amplification')