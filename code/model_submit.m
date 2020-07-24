%% data calculate
xplot=[0:0.01:2,2.1:0.1:4.9,5:5:100];
%xplot=[0:0.01:2];
for i =1:length(xplot)
    Ne1(i)=nn_e_1(xplot(i),0);
    Ne2(i)=nn_e_2(xplot(i),0);
    Ne(i)=Ne1(i)+Ne2(i);
    Ni(i)=nn_i(xplot(i),0);
    int_ve(i)=vvreal(xplot(i),0);
    int_vi(i)=vvreali(xplot(i),0);
    Pe(i)=integral2(@(vx,vy)f_e_rec__(vx,vy,xplot(i),0).*vx.^2*v0^2,-3,3,-3,3)*params.me*10^6; %Prr
    Pi(i)=integral2(@(vx,vy)f_i_rec__(vx,vy,xplot(i),0).*vx.^2*v0i^2,-3,3,-3,3)*params.mp*10^6;
    Teperp(i)=Pe(i)/Ne(i)/1e6/params.e;
    Tepara(i)=(params.te_*params.lambda*Ne1(i)+params.te_*Ne2(i))/Ne(i);
    Tiperp(i)=Pi(i)/Ni(i)/1e6/params.e;
    Tipara(i)=Tiperp(i);
end  
% Prr P=nkT (not correct)
Pepara=Ne.*1e6*params.e.*Tepara; 
Pipara=Ni.*1e6*params.e.*Tipara;
% magnetic pressure 
Pb=interp1(xnew,-(-(y(2,:)*params.B0*1e-5).^2/(8*pi))*1e-1*1e9,xplot,'spline');
% drfit velocity for proton and electron (diamagnetic/E*B) 
vPi= 1e-3*(dif2(Ni,xplot)./Ni.*Tiperp+dif2(Tiperp,xplot))*1.6e-19 ./interp1(xnew,y(2,:),xplot,'spline')*(1/params.B0*1e9/params.D*100/1.6e-19);
vPe= interp1(momentum{6,1},momentum{1,1},xplot,'spline');
vE = interp1(momentum{6,1},momentum{2,1},xplot,'spline');
% current
j=params.e*1e6*Ne.*((vE-vPi)-(vE+vPe))*1e3*1e12*1e3;
% hall term
hall=-j.*interp1(xnew,y(2,:),xplot,'spline')./Ne*params.B0*10^-9/params.e/10^6/10^12;

vv=(params.mp*(vE-vPi)+params.me*(vE+vPe))*1e3/(params.mp+params.me);
% convection term
convection=-vv.*interp1(xnew,y(2,:),xplot,'spline').*params.B0*10^-9;

% inertial term
inertial=params.me/params.e.*int_ve.^2./xplot/params.D*1e2;

% divergence P term
divp=interp1(momentum{6,1},momentum{1,1}*1e3.*interp1(xnew,Bnew,A{6,1},'spline')*1e-9,xplot,'spline');
totalE=hall-convection-divp;
% temperature Trr  Tphiphi
Terr=interp1(momentum{6,1},momentum{7,1},xplot,'spline')./Ne/1e6/params.e;
Teff=interp1(momentum{6,1},momentum{8,1},xplot,'spline')./Ne/1e6/params.e;
%% save data
line_data{1,1}='Ne1 cm^-3'; line_data{1,2}=Ne1;
line_data{2,1}='Ne2 cm^-3'; line_data{2,2}=Ne2;
line_data{3,1}='Ne cm^-3';  line_data{3,2}=Ne;
line_data{4,1}='Ni cm^-3';  line_data{4,2}=Ni;
line_data{5,1}='int_ve m/s';line_data{5,2}=int_ve;
line_data{6,1}='int_vi m/s';line_data{6,2}=int_vi;
line_data{7,1}='Pe Pa';  line_data{7,2}=Pe;
line_data{8,1}='Pi Pa';  line_data{8,2}=Pi;
line_data{9,1}='Teperp eV'; line_data{9,2}=Teperp;
line_data{10,1}='Tepara eV';line_data{10,2}=Tepara;
line_data{11,1}='Tiperp eV';line_data{11,2}=Tiperp;
line_data{12,1}='Tipara eV';line_data{12,2}=Tipara;
line_data{13,1}='Pepara Pa';line_data{13,2}=Pepara;
line_data{14,1}='Pipara Pa';line_data{14,2}=Pipara;
line_data{15,1}='Pb Pa ';    line_data{15,2}=Pb/1e9;
line_data{16,1}='j A/m^2';     line_data{16,2}=j/1e15;
line_data{17,1}='hall V/m';  line_data{17,2}=hall/1e3;
line_data{18,1}='conve V/m';  line_data{18,2}=convection;
line_data{19,1}='divp V/m'; line_data{19,2}=divp; 
line_data{20,1}='vPi km/s'; line_data{20,2}=vPi; 
line_data{21,1}='vPe km/s'; line_data{21,2}=vPe; 
line_data{22,1}='vE km/s'; line_data{22,2}=vE;
line_data{23,1}='Terr eV'; line_data{23,2}=Terr;
line_data{24,1}='Teff eV'; line_data{24,2}=Teff;
line_data{25,1}='inertial V/m'; line_data{25,2}=inertial;
line_data{25,1}='xplot'; line_data{25,2}=xplot;

save model_end line_data

function y=dif2(x,t)
    x_diff=diff(x)./diff(t);
    y=interp1((t(1:end-1)+t(2:end))/2,x_diff,t,'spline',0);
end