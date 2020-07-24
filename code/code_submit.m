%% params_input
params.c=2.99792458e8; % light spEed
params.e=1.60218e-19;  % qe
params.mp=1.6726231e-27;  %proton mass
params.me=9.1093897e-31;  %electron mass
params.m0=params.me/params.mp;  % me/mp
params.bz = 17;           % Bz at the cavity center
params.n0e = 23.7;        % electron number density (not same as at the center)
params.k=0.3;             % k=D/di;  D: normalized length  di:ion inertial length
% in fact, the D is related to the angular velocity.
% the k helps to fix the D, then fix the omega.
% details shown in [Shustov et al., 2016]
params.delta = 0.88;      % the proportion of background plasma
params.omega=-0.003;      % ratio of angular velocity: omegai/omegae 
% 0 represent background; 1 represent current-carrying 
params.lambda=0.87;       % ratio of temperature : te1/te0
params.te_=52;            % te0  background electron temperature
params.tau_=265/params.te_;   % ratio: ti0/te0
params.tau=(290-params.delta*params.tau_*params.te_)/(1-params.delta)/params.te_/params.lambda;  % ratio: ti1/te1
params.te = params.te_*params.lambda;  % te1
params.ti_ = params.te_*params.tau_;   % ti0
params.ti = params.te*params.tau;      % ti1
params.z=1;               % qi : proton  z=1
params.m = params.m0/params.z;    % m0/z
params.T0 = (params.te_*params.lambda*(params.z+params.tau))*params.e*1e7;  %ti1+te1
params.B0 = 1e5*sqrt(4*pi*params.n0e*params.T0);  % normalized B
params.be=10.3619/params.B0;  % μ coefficient for current-carrying
params.be_=-5.8874/params.B0;  % μ coefficient for background

params.di = sqrt(2.998e10^2*params.mp*1e3/(4*pi*params.n0e*(2.998e9*1.6e-19)^2)); % ion inertial length
params.D=params.di*params.k; % normalized length
params.Omega0=-1/sqrt((params.D)^4/(3*10^10)^2/(params.T0)*4*pi*(params.e*3*10^9)^2*params.n0e);  % omegai-omegae
params.Omegai = params.Omega0 * params.omega/(params.omega-1); % omegai  ion angular velocity
params.Omegae = params.Omega0/(params.omega-1); % omegae  electron angular velocity
params.ae=params.m/(params.m + params.omega^2); 
params.ai=params.omega^2/(params.m + params.omega^2);
params.epsilon = (params.omega^2+params.m)*params.z/params.k^2/(params.omega-1)^2;   
params.nita=params.c*1e2*params.T0/params.e/params.c/10/params.Omega0;   % normalized from Shustov 2016;
v0 = sqrt(params.te*params.e*2/params.me);  % electron thermal velocity
v0i = sqrt(params.ti*params.e*2/params.mp);  % proton thermal velocity
%% load 初始电磁场
inputs=[0:0.001:0.2,0.21:0.01:5 5.5:0.5:500 550:50:10000];  % inputs=ρ^2/(2*params.D^2)
% ode45 sloving the equation
% first run, with uniform B and zero E, that is the mu=(vr^2+vphi^2)/B
sol = ode45(@(x,y)vpt2(x,y,params),inputs, [0;-params.bz/params.B0;0]);
% y(1,:) A*rou/params.nita
% y(2,:) normalized B
% y(3,:) electrostatic potential   phi*qe/params.T0
% details for these generalized variables, shown in Shustov 2016;

[y, yp] = deval(inputs,sol);  % y:3*x matrix  yp:the derivative
x = sqrt(2*inputs);  % ρ  ;  unit:params.D
E = -yp(3,:).*x*(params.te_*params.lambda+params.te_*params.lambda*params.tau)/params.D*100; % E ; unit V/m
B = -y(2,:)*params.B0; % B ; unit nT
data={y,x,E,B}; % save into cell data;
disp('initial state finished')
%% ode+电磁场存储
jieshumax=4;  % iteration times
for i=1:jieshumax
    a=max(data{i,2});
    a=a-0.9;
    a=a^2/2;     
    inputs=[0:0.001:0.2,0.21:0.01:5 5.5:0.5:500 550:50:a];  % adjusting the right limit of the inputs 
    %( for the next iteration, the guiding center may not be in the inputs before)
    sol = ode45(@(xx,y)vptl3(xx,y,data{i,2},data{i,3},data{i,4},params),inputs, [0;-params.bz/params.B0;0]);
    [y, yp] = deval(inputs,sol);
    x = sqrt(2*inputs); % ρ  /params.D(cm)
    E = -yp(3,:).*x*(params.te_*params.lambda+params.te_*params.lambda*params.tau)/params.D*100; %E/ V/m
    B = -y(2,:)*params.B0; % B/nT
    data{i+1,1}=y;  
    data{i+1,2}=x;
    data{i+1,3}=E;
    data{i+1,4}=B;
    disp(i)
end
disp('ode finished')
%plot(data{2,2},data{2,4});
jieshu=jieshumax;
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
        % calculate the position of guiding center
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
% distubution function for current-carrying population
f_e_rec_1= @(vx,vy,xx,yy)params.n0e*(1-params.delta)*(1/pi)*exp(...
        -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy) ...
        + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
        (1+params.tau)*(interp1(xnew,y(3,:),rr(xx,yy),'spline')-1/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline'))); 
% distribution function for background population
f_e_rec_2=@(vx,vy,xx,yy)params.n0e*params.delta*(params.lambda/pi)*exp(...
        -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy)) ...
        + params.lambda*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline'));
f_e_rec__=@(vx,vy,xx,yy)f_e_rec_1(vx,vy,xx,yy)+f_e_rec_2(vx,vy,xx,yy);
nn_e=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-3,3,-3,3);
% to make ne=ni at r=0;
params.n0i=nn_e(0,0);
% distribution function for proton
f_i_rec__=@(vx,vy,xx,yy)params.n0i*(1-params.delta)*(1/pi)^1*exp(...
        -(vx.^2+vy.^2)+sign(params.Omegai)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ai*params.epsilon*(1+params.tau)/params.tau)+...
        (1+params.tau)/params.tau*(params.omega/(params.omega-1)*interp1(xnew,y(1,:),rr(xx,yy),'spline')-interp1(xnew,y(3,:),rr(xx,yy),'spline')))+...
            params.n0i*params.delta*(params.lambda*params.tau/pi/params.tau_)^1*exp(...
        -params.lambda*params.tau/params.tau_*(vx.^2+vy.^2)-...
        params.lambda/params.tau_*(1+params.tau)*interp1(xnew,y(3,:),rr(xx,yy),'spline')); 
nn_i=@(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,xx,yy),-3,3,-3,3); 
%% save
save data_end data
%% function    
function dydx=vptl3(xx,y,x,E,B,params)
        dydx=zeros(3,1);  
        lambda = params.lambda;
        tau = params.tau;
        tau_ = params.tau_;
        delta = params.delta;
        k = params.k;
        m = params.m/1;
        z = params.z;
        be = params.be;
        be_ = params.be_;
        % variable qujian: integrate limit (to make sure μ conserved)
        qujian=3;
        v0 = sqrt(params.te*1.6e-19*2/params.me);
        rr=@(xx,yy)sqrt(xx.^2+yy.^2);
        E_e=@(xx)interp1(x,E,xx,'spline');
        B_e=@(xx)interp1(x,B,xx,'spline');
        % calculate the position of guiding center
        rc=@(vx,vy,xx,yy)params.me*v0/B_e(rr(xx,yy))*1e9/params.e/params.D*100*sqrt(vx.^2+vy.^2); %guiyihua
        rc_x=@(vx,vy,xx,yy)xx+rc(vx,vy,xx,yy).*(-vy).*judge_0(vx,vy);
        rc_y=@(vx,vy,xx,yy)yy+rc(vx,vy,xx,yy).*vx.*judge_0(vx,vy);
        r_=@(vx,vy,xx,yy)sqrt(rc_x(vx,vy,xx,yy).^2+rc_y(vx,vy,xx,yy).^2);
        % function handles to calculate drift velocity 
        % E*B drift
        vE=@(vx,vy,xx,yy)-E_e(r_(vx,vy,xx,yy))./B_e(r_(vx,vy,xx,yy))*1e9;
        vf=@(vx,vy,vz,xx,yy)(vx*(-yy)+vy*xx).*judge_0(xx,yy);
        B_tidu=dif2(B,x)/params.D*100*1e-9;
        B_tidu_e=@(r)interp1(x,B_tidu,r,'spline');
        % gradient B drift
        vD=@(vx,vy,xx,yy)1/2*params.me.*(vx.^2+vy.^2)*v0^2./params.e./B_e(r_(vx,vy,xx,yy)).^2 ...
            .*B_tidu_e(r_(vx,vy,xx,yy))*1e18; 
        vdrift=@(vx,vy,xx,yy)vE(vx,vy,xx,yy)-vD(vx,vy,xx,yy);
        vdrift_x=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(-rc_y(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        vdrift_y=@(vx,vy,xx,yy)vdrift(vx,vy,xx,yy).*(rc_x(vx,vy,xx,yy)).*judge_0(rc_x(vx,vy,xx,yy),rc_y(vx,vy,xx,yy));
        
        % calculate miu
        miu_guiyi1=@(vx,vy,xx,yy,y2)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
                *params.be./(-y2)/v0^2;
        miu_guiyi2=@(vx,vy,xx,yy,y2)((vx*v0-vdrift_x(vx,vy,xx,yy)).^2+(vy*v0-vdrift_y(vx,vy,xx,yy)).^2) ...
                *params.be_./(-y2)/v0^2;
            
        % handle for electron distribution functions
        % current-carrying 
        f_e_rec_1= @(vx,vy,xx,yy,y2,y1,y3)params.n0e*(1-params.delta)*(1/pi)*exp(...
            -vx.^2 - vy.^2 + miu_guiyi1(vx,vy,xx,yy,y2) ...
            + sign(params.Omegae)*vf(vx,vy,0,xx,yy).*rr(xx,yy)*sqrt(2*params.ae*params.epsilon*(1+params.tau)) +...
            (1+params.tau)*(y3-1/(params.omega-1)*y1)); 
        % background
        f_e_rec_2= @(vx,vy,xx,yy,y2,y1,y3)params.n0e*params.delta*(params.lambda/pi)*exp(...
            -params.lambda*(vx.^2+vy.^2-miu_guiyi2(vx,vy,xx,yy,y2)) ...
            + params.lambda*(1+params.tau)*y3);
        % total
        f_e_rec__=@(vx,vy,xx,yy,y2,y1,y3)f_e_rec_1(vx,vy,xx,yy,y2,y1,y3)+f_e_rec_2(vx,vy,xx,yy,y2,y1,y3);
    % Ee: number density of curreny-carrying electrons
    Ee=integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    params.n0i=integral2(@(vx,vy)f_e_rec__(vx,vy,1e-7,0,-params.bz/params.B0,0,0),-qujian,qujian,-qujian,qujian);
    nie = params.n0i/params.n0e;
    Ee_func=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*t),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    % the handles below helps to calculate the different derivations
    Ee_func_y2=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,t,y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y1=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),t,y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    Ee_func_y3=@(t)integral2(@(vx,vy)f_e_rec_1(vx,vy,sqrt(2*xx),0,y(2),y(1),t),-qujian,qujian,-qujian,qujian)/params.n0e/(1-params.delta);
    % Ee_: number density of background electrons
    Ee_=integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*t),0,y(2),y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func_y3=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,y(2),y(1),t),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ee__func_y2=@(t)integral2(@(vx,vy)f_e_rec_2(vx,vy,sqrt(2*xx),0,t,y(1),y(3)),-qujian,qujian,-qujian,qujian)/params.n0e/(params.delta);
    Ei=exp((z+tau)/tau*(params.ai*params.epsilon*xx+(params.omega/(params.omega-1))*y(1)-y(3)));
    Ei_=exp(-lambda*(1+tau)/tau_*y(3));
    dydx(1)=y(2);
    if xx==0   % to avoid the denominator=0
        xx_=1e-7;
        % first-order momentum, electron current density
        je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx_),0,y(2),y(1),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
        % curl of magnetic field
        % the right side: electron current and proton current 
        % normalized process dividing (params.Omega0*ρ)
        dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx_*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
        
        % quasi-neutral equation
        Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei)-(Ee_func_y1(y(1)+1e-7)-Ee)/1e-7; % ?（ni-ne)/?y(1)
        PPsi=-(Ee_func_y2(y(2)+1e-7)-Ee)/1e-7-(Ee__func_y2(y(2)+1e-7)-Ee_)/1e-7*delta/(1-delta); % ?（ni-ne)/?y(2)
        Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau)-(Ee_func(xx+1e-7)-Ee)/1e-7-(Ee__func(xx+1e-7)-Ee_)/1e-7*delta/(1-delta);  % ?（ni-ne)/?x
        K=(Ee_func_y3(y(3)+1e-7)-Ee)/1e-7+(Ee__func_y3(y(3)+1e-7)-Ee_)/1e-7*delta/(1-delta)+nie/tau*Ei*(1+tau)+(1+tau)*nie*lambda/tau_*delta/(1-delta)*Ei_; % ?（ne-ni)/?y(3) 
        dydx(3) = (Psi*y(2)+PPsi*dydx(2)+Theta)/K;   % ?y(3)/?x=(?（ni-ne)/?y(1)*?y(1)/?x  + ?（ni-ne)/?y(2)*?y(2)/?x
                                                     %   +?（ni-ne)/?x)/?（ne-ni)/?y(3)
                                                     % or  d(ni-ne)/dx=0
    else
        % first-order momentum
        je =integral2(@(vx,vy)f_e_rec__(vx,vy,sqrt(2*xx),0,y(2),y(1),y(3)).*vy*v0,-qujian,qujian,-qujian,qujian);
        % curl of magnetic field
        dydx(2)=(1-delta)*(je/params.Omega0/sqrt(xx*params.D^2*2/1e4)/params.n0e/(1-delta) - nie*Ei*params.omega/(params.omega-1));
        % quasi-neutral
        Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei)-(Ee_func_y1(y(1)+1e-7)-Ee)/1e-7;  % ?（ni-ne)/?y(1)
        PPsi=-(Ee_func_y2(y(2)+1e-7)-Ee)/1e-7-(Ee__func_y2(y(2)+1e-7)-Ee_)/1e-7*delta/(1-delta); % ?（ni-ne)/?y(2)
        Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau)-(Ee_func(xx+1e-7)-Ee)/1e-7-(Ee__func(xx+1e-7)-Ee_)/1e-7*delta/(1-delta);   % ?（ni-ne)/?x
        K=(Ee_func_y3(y(3)+1e-7)-Ee)/1e-7+(Ee__func_y3(y(3)+1e-7)-Ee_)/1e-7*delta/(1-delta)+nie/tau*Ei*(1+tau)+(1+tau)*nie*lambda/tau_*delta/(1-delta)*Ei_;  % ?（ni-ne)/?y(3)
        dydx(3) = (Psi*y(2)+PPsi*dydx(2)+Theta)/K;  % ?y(3)/?x=(?（ni-ne)/?y(1)*?y(1)/?x  + ?（ni-ne)/?y(2)*?y(2)/?x
                                                     %   +?（ni-ne)/?x)/?（ne-ni)/?y(3)
                                                     % or  d(ni-ne)/dx=0
    end    
end
function dydx=vpt2(x,y,params)
    % params: 常数结构体，含有下列参数
    lambda = params.lambda;
    tau = params.tau;
    tau_ = params.tau_;
    delta = params.delta;
    k = params.k;
    m = params.m/1;
    z = params.z;
    be = params.be;
    be_ = params.be_;
    params.n0i = params.n0e*(1-params.delta)/(1-params.be*params.B0/params.bz) +...
    params.n0e*params.delta/(1-params.be_*params.B0/params.bz);
    nie = params.n0i/params.n0e;
    
    dydx=zeros(3,1);
    xie=be/y(2);
    xie_ = be_/y(2);
    % the number density for four populations
    Ee=exp((z+tau)*(params.ae*params.epsilon*x/(1+xie)-(1/(params.omega-1))*y(1)+y(3)));
    Ei=exp((z+tau)/tau*(params.ai*params.epsilon*x+(params.omega/(params.omega-1))*y(1)-y(3)));
    Ee_=exp(lambda*(1+tau)*y(3));
    Ei_=exp(-lambda*(1+tau)/tau_*y(3));
    
    dydx(1)=y(2);
    
    % curl of magnetic field
    dydx(2)=(1-delta)*(Ee/(params.omega-1)/(1+xie)^2 - nie*Ei*params.omega/(params.omega-1));
    
    % quasi-neutral
    % partial y(1) y(2) x repectively
    Psi=(1+tau)/(params.omega-1)*(nie*params.omega/tau*Ei+Ee/(1+xie));  % % ?（ni-ne)/?y(1)
    PPsi=-Ee*(1+(1+tau)*params.ae*params.epsilon*x/(1+xie))*(be/(be+y(2))^2)...  % ?（ni-ne)/?y(2)
        - Ee_*delta/(1-delta)*(be_/(be_+y(2))^2);
    Theta=params.epsilon*(1+tau)*(nie*params.ai*Ei/tau-params.ae*Ee/(1+xie)^2);   % % ?（ni-ne)/?x
    K=(1+tau)*(Ee/(1+xie)+lambda*delta/(1-delta)/(1+xie_)*Ee_+nie/tau*Ei+nie*lambda/tau_*delta/(1-delta)*Ei_); % % ?（ni-ne)/?y(3)
    dydx(3) = (Psi*y(2)+PPsi*dydx(2)+Theta)/K;       % ?y(3)/?x=(?（ni-ne)/?y(1)*?y(1)/?x  + ?（ni-ne)/?y(2)*?y(2)/?x
                                                     %   +?（ni-ne)/?x)/?（ne-ni)/?y(3)
                                                     % or  d(ni-ne)/dx=0
end
% to avoid the error (divide sqrt(x^2+y^2))
function result=judge_0(x,y)
    if x==0 & y==0
        result=0;
    else
        result=1./sqrt(x.^2+y.^2);
    end
end
% derive differential coefficient
function y=dif2(x,t)
    x_diff=diff(x)./diff(t);
    y=interp1((t(1:end-1)+t(2:end))/2,x_diff,t,'spline',0);
end