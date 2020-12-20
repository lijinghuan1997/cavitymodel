%% 
jieshu=4;
v0 = sqrt(params.te*params.e*2/params.me);
v0i = sqrt(params.ti*params.e*2/params.mp);
x=data{jieshu,2};
E=data{jieshu,3};
B=data{jieshu,4};
y=data{jieshu+1,1};
xnew=data{jieshu+1,2};
Enew=data{jieshu+1,3};
Bnew=data{jieshu+1,4};
rr=@(xx,yy)sqrt(xx.^2+yy.^2);
E_e=@(xx)interp1(x,E,xx,'spline');
B_e=@(xx)interp1(x,B,xx,'spline');
% these handles are the same as the ODE function
% to calculate the vD of the guiding center
% to calculate the mu invariant
rc=@(vx,vy,xx,yy)params.me*v0/B_e(rr(xx,yy))*1e9/params.e/params.D*100*sqrt(vx.^2+vy.^2); %guiyihua
rc_x=@(vx,vy,xx,yy)xx+rc(vx,vy,xx,yy).*(-vy).*judge_0(vx,vy);
rc_y=@(vx,vy,xx,yy)yy+rc(vx,vy,xx,yy).*vx.*judge_0(vx,vy);
r_=@(vx,vy,xx,yy)sqrt(rc_x(vx,vy,xx,yy).^2+rc_y(vx,vy,xx,yy).^2);
        % function handles to calcute miu
vE=@(vx,vy,xx,yy)-E_e(r_(vx,vy,xx,yy))./B_e(r_(vx,vy,xx,yy))*1e9;
vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
B_tidu=dif2(B,x)/params.D*100*1e-9;
B_tidu_e=@(r)interp1(x,B_tidu,r,'spline');
vD=@(vx,vy,xx,yy)1/2*params.me*(vx.^2+vy.^2)*v0^2./params.e./B_e(r_(vx,vy,xx,yy)).^2 ...
    .*B_tidu_e(r_(vx,vy,xx,yy))*1e18;
vdrift=@(vx,vy,xx,yy)vE(vx,vy,xx,yy)-vD(vx,vy,xx,yy);
vdrift_x=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(-rc_y(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
vdrift_y=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(rc_x(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
miu_guiyi1=@(vx,vy,xx,yy)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
    *params.be*params.B0./interp1(xnew,Bnew,rr(xx,yy),'spline')/v0^2;
miu_guiyi2=@(vx,vy,xx,yy)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
    *params.be_*params.B0./interp1(xnew,Bnew,rr(xx,yy),'spline')/v0^2;
vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
% the distribution function for the current-carrying population (electron)
f_e_rec_1= @(vx,vy,xx,yy)params.n0e*(1-params.delta)*(1/pi)*exp(...
        -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy) ...
        + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
        (1+params.tau)*(interp1(xnew,y(3,:),rr(xx,yy),'spline')-1/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')));
% the distribution function for the background population (electron)
f_e_rec_2=@(vx,vy,xx,yy)params.n0e*params.delta*(params.lambda/pi)*exp(...
        -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy)) ...
        + params.lambda*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline'));
% total
f_e_rec__=@(vx,vy,xx,yy)f_e_rec_1(vx,vy,xx,yy)+f_e_rec_2(vx,vy,xx,yy);
% the distribution function for the protons
f_i_rec__=@(vx,vy,xx,yy)params.n0i*(1-params.delta)*(1/pi)^1*exp(...
        -(vx.^2+vy.^2)+sign(params.Omegai)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ai*params.epsilon*(1+params.tau)/params.tau)+...
        (1+params.tau)/params.tau*(params.omega/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')-interp1(xnew,y(3,:),rr(xx,yy),'spline')))+...
            params.n0i*params.delta*(params.lambda*params.tau/pi/params.tau_)^1*exp(...
        -params.lambda*params.tau/params.tau_*(vx.^2+vy.^2)-...
        params.lambda/params.tau_*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline'));
% handle for number density
nn_e=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-3,3,-3,3);
params.n0i=nn_e(0,0);
nn_i=@(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,xx,yy),-3,3,-3,3); 
nn_e_1=@(xx,yy)integral2(@(vx,vy)f_e_rec_1(vx,vy,xx,yy),-3,3,-3,3); % current-carrying
nn_e_2=@(xx,yy)integral2(@(vx,vy)f_e_rec_2(vx,vy,xx,yy),-3,3,-3,3); % background
nn_=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-3,3,-3,3);
% handle for bulk velocity
vvreal= @(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,rr(xx,yy),0).*vy*v0,-3,3,-3,3,'AbsTol',1e2,'RelTol',1e-4)/nn_e(rr(xx,yy),0);
vvreali= @(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,rr(xx,yy),0).*vy*v0i,-3,3,-3,3,'AbsTol',1e2,'RelTol',1e-4)/nn_i(rr(xx,yy),0);

%% calculate
xplot=[0:0.02:1,1.1:0.2:5,5.5:1:50,55:5:105];
for i=1:length(xplot)
    Ne1(i)=nn_e_1(xplot(i),0);  % number density of current-carrying population
    Ne2(i)=nn_e_2(xplot(i),0);  % number density of background population
    Ne(i)=Ne1(i)+Ne2(i);        % number density of electrons
    Ni(i)=nn_i(xplot(i),0);     % number density of ions
    int_ve(i)=vvreal(xplot(i),0);  % bulk velocity of electrons
    int_vi(i)=vvreali(xplot(i),0); % bulk velocity of ions
    % pressure tensor (Prr, Pff--Pphiphi)
    Prr(i)=integral2(@(vx,vy)f_e_rec__(vx,vy,xplot(i),0).*vx.^2*v0^2,-3,3,-3,3,'AbsTol',1e10,'RelTol',1e-4)*params.me*10^6; 
    Pff(i)=integral2(@(vx,vy)f_e_rec__(vx,vy,xplot(i),0).*(vy-int_ve(i)/v0).^2*v0^2,-3,3,-3,3,'AbsTol',1e10,'RelTol',1e-4)*params.me*10^6; 
    Pi(i)=integral2(@(vx,vy)f_i_rec__(vx,vy,xplot(i),0).*vx.^2*v0i^2,-3,3,-3,3,'AbsTol',1e10,'RelTol',1e-4)*params.mp*10^6;
    % temperatures
    Terr(i)=Prr(i)/Ne(i)/1e6/params.e;
    Teff(i)=Pff(i)/Ne(i)/1e6/params.e;
    Tepara(i)=(params.te_*params.lambda*Ne1(i)+params.te_*Ne2(i))/Ne(i);
    Tiperp(i)=Pi(i)/Ni(i)/1e6/params.e;
    Tipara(i)=Tiperp(i);
end
% parallel pressures
Pepara=Ne.*1e6*params.e.*Tepara; 
Pipara=Ni.*1e6*params.e.*Tipara;
% gradient of pressure in the radial direction
Pe_tidu=dif2(Prr,xplot)./params.D*1e2+(Prr-Pff)./xplot/params.D*1e2;
% diamagnetic drift of ions
vPi= 1e-3*(dif2(Ni,xplot)./Ni.*Tiperp+dif2(Tiperp,xplot))*1.6e-19 ./interp1(xnew,y(2,:),xplot,'spline')*(1/params.B0*1e9/params.D*100/1.6e-19);
vPi=-vPi*1e3;
% EB drift
vEB=-interp1(xnew,Enew,xplot,'spline')./interp1(xnew,Bnew,xplot,'spline')*1e9;
% diamagnetic drift of electrons
vPe=-Pe_tidu./interp1(xnew,Bnew,xplot,'spline')*1e9./Ne/1e6/params.e;
% plus  (proved to be the same as the integrated int_ve)
vplus=vEB+vPe;

% generalized Ohm's law
% hall term
j=params.e*1e6*Ne.*((vEB+vPi)-(vEB+vPe));
hall=j.*interp1(xnew,Bnew,xplot,'spline')*1e-9./Ne/1e6/params.e;

% convection term
vv=(params.mp*(vEB+vPi)+params.me*(vEB+vPe))/(params.mp+params.me);
convection=-vv.*interp1(xnew,Bnew,xplot,'spline')*10^-9;

% inertial term
inertial=params.me/params.e.*int_ve.^2./xplot/params.D*1e2;

% divergence P term
divp=-Pe_tidu./Ne/1e6/params.e;
% total E        (proved to be the same as the Enew (model results))
totalE=hall+convection+divp;
%%
line_data{1,1}='Ne1 cm^-3'; line_data{1,2}=Ne1;
line_data{2,1}='Ne2 cm^-3'; line_data{2,2}=Ne2;
line_data{3,1}='Ne cm^-3';  line_data{3,2}=Ne;
line_data{4,1}='Ni cm^-3';  line_data{4,2}=Ni;

line_data{5,1}='int_ve m/s';line_data{5,2}=int_ve;
line_data{6,1}='int_vi m/s';line_data{6,2}=int_vi;
line_data{7,1}='vPi m/s'; line_data{7,2}=vPi; 
line_data{8,1}='vPe m/s'; line_data{8,2}=vPe; 
line_data{9,1}='vEB m/s'; line_data{9,2}=vEB;
line_data{10,1}='vplus m/s'; line_data{10,2}=vplus;

line_data{11,1}='Prr Pa';  line_data{11,2}=Prr;
line_data{12,1}='Pff Pa';  line_data{12,2}=Pff;
line_data{13,1}='Pepara Pa';line_data{13,2}=Pepara;
line_data{14,1}='Pi Pa';  line_data{14,2}=Pi;
line_data{15,1}='Pipara Pa';line_data{15,2}=Pipara;
line_data{16,1}='Terr eV'; line_data{16,2}=Terr;
line_data{17,1}='Teff eV';line_data{17,2}=Teff;
line_data{18,1}='Tepara eV';line_data{18,2}=Tepara;
line_data{19,1}='Tiperp eV';line_data{19,2}=Tiperp;
line_data{20,1}='Tipara eV';line_data{20,2}=Tipara;

line_data{21,1}='j A/m^2';     line_data{21,2}=j;
line_data{22,1}='hall V/m';  line_data{22,2}=hall;
line_data{23,1}='conve V/m';  line_data{23,2}=convection;
line_data{24,1}='divp V/m'; line_data{24,2}=divp; 
line_data{25,1}='inertial V/m'; line_data{25,2}=inertial;

line_data{26,1}='xplot'; line_data{26,2}=xplot;
save line_data.mat line_data
%% example plot
xwidth=0.7;
ywidth=0.1;
% magnetic field (nT)
subplot('position',[0.15 0.97-ywidth xwidth ywidth]);
semilogx(xnew*params.D/1e5,Bnew)
xlim([1,1000])
ylim([10,40])
set(gca,'xtick',[])
set(gca,'ytick',[20 30])
ylabel('B/nT')
% electric field (V/m)
subplot('position',[0.15 0.97-2*ywidth xwidth ywidth]);
semilogx(xnew*params.D/1e5,Enew)
xlim([1,1000])
set(gca,'xtick',[])
set(gca,'ytick',[-1e-3 0 1e-3])
ylabel('E/ V/m')
% number density cm^-3
subplot('position',[0.15 0.97-3*ywidth xwidth ywidth]);
semilogx(xplot*params.D/1e5,Ne)
hold on
semilogx(xplot*params.D/1e5,Ni)
xlim([1,1000])
set(gca,'xtick',[])
set(gca,'ytick',[18 22])
ylabel('N/cm^-3')
% electron temperatures eV
subplot('position',[0.15 0.97-4*ywidth xwidth ywidth]);
semilogx(xplot*params.D/1e5,Terr)
hold on
semilogx(xplot*params.D/1e5,Teff)
hold on
semilogx(xplot*params.D/1e5,Tepara)
xlim([1,1000])
set(gca,'xtick',[])
set(gca,'ytick',[45 55])
ylabel('Te/eV')
% electron bulk velcity m/s
subplot('position',[0.15 0.97-5*ywidth xwidth ywidth]);
semilogx(xplot*params.D/1e5,vEB)
hold on
semilogx(xplot*params.D/1e5,vPe)
hold on
semilogx(xplot*params.D/1e5,int_ve)
xlim([1,1000])
ylim([-0.5e5, 2e5])
set(gca,'xtick',[])
set(gca,'ytick',[0 1e5])
ylabel('Ve/ m/s')
% ion bulk velocity m/s
subplot('position',[0.15 0.97-6*ywidth xwidth ywidth]);
semilogx(xplot*params.D/1e5,vEB)
hold on
semilogx(xplot*params.D/1e5,vPi)
hold on
semilogx(xplot*params.D/1e5,int_vi)
xlim([1,1000])
set(gca,'xtick',[])
set(gca,'ytick',[-5e4 0 5e4])
ylabel('Vi/ m/s')
% Ohm's law  electric field  V/m
subplot('position',[0.15 0.97-7*ywidth xwidth ywidth]);
semilogx(xplot*params.D/1e5,hall)
hold on
semilogx(xplot*params.D/1e5,divp)
hold on
semilogx(xplot*params.D/1e5,convection)
hold on
semilogx(xnew*params.D/1e5,Enew)
xlim([1,1000])
set(gca,'xtick',[])
set(gca,'ytick',[-2e-3 2e-3])
ylabel('E/ V/m')
%% 3-D velocity distributions
% PSD
f_e_= @(vx,vy,vz,xx,yy)params.n0e*(1-params.delta)*(1/pi)^(3/2)*exp(...
        -vx.^2 - vy.^2 - vz.^2 + miu_guiyi1(vx,vy,xx,yy) ...
        + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
        (1+params.tau)*(interp1(xnew,y(3,:),rr(xx,yy),'spline')-1/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')))... 
      +params.n0e*params.delta*(params.lambda/pi)^(3/2)*exp(  ...
        -params.lambda*(vx.^2+vy.^2+vz.^2-miu_guiyi2(vx,vy,xx,yy)) ...
        + params.lambda*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline'));
f_e=@(vx,vy,vz,xx)f_e_(vx,vy,vz,xx,0);
f_i_=@(vx,vy,vz,xx,yy)params.n0i*(1-params.delta)*(1/pi)^(3/2)*exp(...
        -(vx.^2+vy.^2+vz.^2)+sign(params.Omegai)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ai*params.epsilon*(1+params.tau)/params.tau)+...
        (1+params.tau)/params.tau*(params.omega/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')-interp1(xnew,y(3,:),rr(xx,yy),'spline')))+...
            params.n0i*params.delta*(params.lambda*params.tau/pi/params.tau_)^(3/2)*exp(...
        -params.lambda*params.tau/params.tau_*(vx.^2+vy.^2+vz.^2)-...
        params.lambda/params.tau_*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline')); 
f_i=@(vx,vy,vz,xx)f_i_(vx,vy,vz,xx,0);
% energy flux (v-speed; t-theta angle; p-phi angle)
g_e_gse = @(v,t,p,xx,phi)1e-4*1e6*v0^1*v.^4.*f_e(v.*sin(t).*cos(p) - v_mh_e(1)*cos(phi),...
    v.*sin(t).*sin(p) - v_mh_e(1)*sin(phi), v.*cos(t), xx)/2;%-v_mh_e(3)
%% functions
function result=judge_0(x,y)
    if x==0 & y==0
        result=0;
    else
        result=1./sqrt(x.^2+y.^2);
    end
end

function y=dif2(x,t)
    x_diff=diff(x)./diff(t);
    y=interp1((t(1:end-1)+t(2:end))/2,x_diff,t,'spline',0);
end
