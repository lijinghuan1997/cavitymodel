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
nn_i=@(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,xx,yy),-3,3,-3,3); 
nn_e_1=@(xx,yy)integral2(@(vx,vy)f_e_rec_1(vx,vy,xx,yy),-3,3,-3,3); % current-carrying
nn_e_2=@(xx,yy)integral2(@(vx,vy)f_e_rec_2(vx,vy,xx,yy),-3,3,-3,3); % background
nn_=@(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy),-3,3,-3,3);
% handle for bulk velocity
vvreal= @(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,rr(xx,yy),0).*vy*v0,-3,3,-3,3)/nn_e(rr(xx,yy),0);
vvrealx= @(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,rr(xx,yy),0).*vx*v0,-3,3,-3,3)/nn_e(rr(xx,yy),0);
vvreali= @(xx,yy)integral2(@(vx,vy)f_i_rec__(vx,vy,rr(xx,yy),0).*vy*v0i,-3,3,-3,3)/nn_i(rr(xx,yy),0);
 % compare f_tidu+E¡ÁB and int_v
    yn=@(t)interp1(xnew,y(2,:),t,'spline');
    vvx=@(xx,yy)vvreal(xx,yy)*(-yy)/rr(xx,yy);
    vvy=@(xx,yy)vvreal(xx,yy)*xx/rr(xx,yy);
    vf_=@(xx,yy)(vvx(xx,yy)*(-yy)+vvy(xx,yy)*xx)/rr(xx,yy);
    % distance from the cavity center
    xplot_f=[0:0.02:1,1.1:0.2:5,5.5:1:50,55:5:105];
    %xplot_f=xplot_f(1:65);
    vvreal= @(xx,yy)integral2(@(vx,vy)f_e_rec__(vx,vy,rr(xx,yy),0).*vy*v0,-3,3,-3,3)/nn_e(rr(xx,yy),0);
    for i=2:length(xplot_f)
    rr_=xplot_f(i);
    yy=1e-5;
    xx=sqrt(rr_^2-yy.^2);
    % calculate the corresponding velocity in xyz coordinate
    v_1=[vvreal(xx,yy+0.00001),vvreal(xx,yy),vvreal(xx+0.00001,yy)];
    v_2=[vvreal(xx,yy-0.00001),vvreal(xx,yy),vvreal(xx-0.00001,yy)];
    vx_1=[v_1(1)*(-yy-0.00001)/rr(xx,yy+0.00001),v_1(2)*(-yy)/rr(xx,yy),v_1(3)*(-yy)/rr(xx+0.00001,yy)];
    vx_2=[v_2(1)*(-yy+0.00001)/rr(xx,yy-0.00001),v_2(2)*(-yy)/rr(xx,yy),v_2(3)*(-yy)/rr(xx-0.00001,yy)];
    vy_1=[v_1(1)*xx/rr(xx,yy+0.00001),v_1(2)*xx/rr(xx,yy),v_1(3)*(xx+0.00001)/rr(xx+0.00001,yy)];
    vy_2=[v_2(1)*xx/rr(xx,yy-0.00001),v_2(2)*xx/rr(xx,yy),v_2(3)*(xx-0.00001)/rr(xx-0.00001,yy)];
    % calculate the 4 points (pressure tensor) around the point (xx,yy), to get the gradient
    % it's not convenient to make a hendle for the calculation of pressure tensor (we need to calculate the bulk velocity at first)
    ppxl1=integral2(@(vx,vy)f_e_rec__(vx,vy,xx+0.00001,yy).*(vx-vx_1(3)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    ppxr1=integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy+0.00001).*(vx-vx_1(1)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    ppyl1=integral2(@(vx,vy)f_e_rec__(vx,vy,xx+0.00001,yy).*(vy-vy_1(3)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    ppyr1=integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy+0.00001).*(vy-vy_1(1)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    
    ppxl2=integral2(@(vx,vy)f_e_rec__(vx,vy,xx-0.00001,yy).*(vx-vx_2(3)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    ppxr2=integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy-0.00001).*(vx-vx_2(1)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    ppyl2=integral2(@(vx,vy)f_e_rec__(vx,vy,xx-0.00001,yy).*(vy-vy_2(3)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    ppyr2=integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy-0.00001).*(vy-vy_2(1)/v0).^2*v0^2,-3,3,-3,3)*params.me*10^6;
    
    ppxyl1=integral2(@(vx,vy)f_e_rec__(vx,vy,xx+0.00001,yy).*(vx-vx_1(3)/v0).*(vy-vy_1(3)/v0).*v0^2,-3,3,-3,3)*params.me*10^6;
    ppxyl2=integral2(@(vx,vy)f_e_rec__(vx,vy,xx-0.00001,yy).*(vx-vx_2(3)/v0).*(vy-vy_2(3)/v0).*v0^2,-3,3,-3,3)*params.me*10^6;
    ppxyr1=integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy).*(vx-vx_1(2)/v0).*(vy-vy_1(2)/v0).*v0^2,-3,3,-3,3)*params.me*10^6;
    ppxyr2=integral2(@(vx,vy)f_e_rec__(vx,vy,xx,yy-0.00001).*(vx-vx_2(1)/v0).*(vy-vy_2(1)/v0).*v0^2,-3,3,-3,3)*params.me*10^6;
    
    % calculate the gradient
    ppx_tidu_1=(ppxl1-ppxl2)/params.D*100/0.00002; %?Pxx/?x
    ppx_tidu_2=(ppxyr1-ppxyr2)/params.D*100/0.00001; %?Pxy/?y
    ppx_tidu_=ppx_tidu_1+ppx_tidu_2; % we mainly care about
    ppy_tidu_1=(ppyr1-ppyr2)/params.D*100/0.00002; %?Pyy/?y
    ppy_tidu_2=(ppxyl1-ppxyl2)/params.D*100/0.00002; %?Pyx/?x
    ppy_tidu_=ppy_tidu_1+ppy_tidu_2;
    
    % calculate diamagnetic drift
    % gradient P to diamagnetic drift
    vvx_tidu_=   ppy_tidu_./nn_(xx,yy)/params.B0/(-yn(rr(xx,yy)))*1e9/1e6/params.e/1e3;
    vvx_tidu_1= ppy_tidu_1./nn_(xx,yy)/params.B0/(-yn(rr(xx,yy)))*1e9/1e6/params.e/1e3;
    vvx_tidu_2= ppy_tidu_2./nn_(xx,yy)/params.B0/(-yn(rr(xx,yy)))*1e9/1e6/params.e/1e3;
    vvy_tidu_=  -ppx_tidu_./nn_(xx,yy)/params.B0/(-yn(rr(xx,yy)))*1e9/1e6/params.e/1e3;
    vvy_tidu_1=-ppx_tidu_1./nn_(xx,yy)/params.B0/(-yn(rr(xx,yy)))*1e9/1e6/params.e/1e3;
    vvy_tidu_2=-ppx_tidu_2./nn_(xx,yy)/params.B0/(-yn(rr(xx,yy)))*1e9/1e6/params.e/1e3;
    f_tidu_ =(vvx_tidu_*(-yy)+vvy_tidu_*xx)*judge_0(xx,yy);
    [xx,yy]
    [ppx_tidu_,ppy_tidu_]
    Pe_tidu_ = (ppx_tidu_*xx+vvy_tidu_*yy)*judge_0(xx,yy);
    % E * B drift
    vEB=-interp1(xnew,Enew,rr_,'spline')/interp1(xnew,Bnew,rr_,'spline')*1e9/1e3;
    f_tidu_+vEB-v_1(2)/1e3;
    [f_tidu_,vEB,v_1(2)/1e3]
    f_tidu(i)=f_tidu_;
    E_B(i)=vEB;
    Pe_tidu(i)=Pe_tidu_;
    plus(i)=f_tidu(i)+E_B(i);
    % momentum int_ve electron bulk velocity
    int_v(i)=v_1(2)/1e3;
    % Prr
    prr(i)=(ppxl1+ppxl2)/2;
    % Pphiphi
    pff(i)=(ppyl1+ppyl2)/2;
    % ?Pxx/?x
    prrtidu(i)=ppx_tidu_1;
    % ?Pyy/?y
    pfftidu(i)=ppy_tidu_1;
    % ?Pxy/?y
    prftidu(i)=ppx_tidu_2;
    % Pxy/Prphi = 0 in fact (calculation error)
    prf(i)=ppxyr2;
    i
 end
prr(1)=prr(2);
pff(1)=pff(2);
prf(1)=prf(2);
prrtidu(1)=prrtidu(2);
pfftidu(1)=pfftidu(2);
prftidu(1)=prftidu(2);
momentum={'f_tidu',f_tidu;
    'E*B',E_B;
    'plus',plus;
    'int_v',int_v;
    'Pe_tidu',Pe_tidu;
    'xplot_f',xplot_f;
    'prr',prr;
    'pff',pff;
    'prf',prf;
    'prrtidu',prrtidu;
    'pfftidu',pfftidu;
    'prftidu',prftidu};  
% the momentum including:   1) diamagnetic drift velocity
%                    2) E ¡Á B drift
%                    3) diamagnetic + E * B
%                    4) the first momentum result for bulk velocity
%                    5) gradient of pressure (total)
%                    6) xplot_f  (distance from the center)
%                    7) pressure tensor Pxx/Prr
%                    8) pressure tensor Pyy/Pphiphi
%                    9) pressure tensor Pxy/Prphi
%                    10) ?Pxx/?x
%                    11) ?Pyy/?y
%                    12) ?Pxy/?y
%                   
%   ?Pxy/?y is equal to the 1/¦Ñ*(Prr-Pphiphi£©in the cylindrical coordinates
save momentum_end momentum

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
