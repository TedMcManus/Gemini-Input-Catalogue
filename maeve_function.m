direcgrid = 'test_run/';
if ~exist('xg','var')
    xg = gemini3d.read.grid(direcgrid);
end
cfg = gemini3d.read.config(direcgrid);
[X2,X3,IT] = ndgrid(xg.x2(3:end-2),xg.x3(3:end-2),0:1:cfg.tdur);

% Magnetospheric Accelerated Energetic inverted-V Electrons
pars.maeve.driftE = usrq('Eastward drift (0 m s^-1): ',0);
pars.maeve.driftN = usrq('Northward drift (-120 m s^-1): ',-120);
pars.maeve.ctr_spn = usrq('Contour North-South span (1e5 m): ',1e5);
pars.maeve.ctr_slp = usrq('Contour central slope (0.5): ',0.5);
pars.maeve.bar_frc = usrq('Precip. loading bar fraction (1): ',1);
pars.maeve.bar_pos = usrq('Precip. loading bar starting position (4.5e5 m): ',4.5e5);
pars.maeve.bar_vel = usrq('Precip. loading bar velocity (-3000 m s^-1): ',-3000);
pars.maeve.bar_gsl = usrq('Precip. loading bar gradient scale length (1e5 m): ',1e5);
pars.maeve.dim_frc = usrq('Precip. fade in fraction (0): ',0);
pars.maeve.dim_tim = usrq('Precip. fade in time (100 s): ',100);
pars.maeve.dim_del = usrq('Precip. fade in delay (100 s): ',100);
pars.maeve.Q_amp_h = usrq('Total energy flux peak hi (15 erg s^-1 cm^-2): ',15);
pars.maeve.Q_amp_l = usrq('Total energy flux peak lo (3 erg s^-1 cm^-2): ',3);
pars.maeve.Q_wth_h = usrq('Total energy flux width hi (2e4 m): ',2e4);
pars.maeve.Q_wth_l = usrq('Total energy flux width lo (1.8e5 m): ',2e5);
pars.maeve.Q_gsl_h = usrq('Total energy flux gradient scale length hi (0.4): ',0.4);
pars.maeve.Q_gsl_l = usrq('Total energy flux gradient scale length lo (0.1): ',0.1);
pars.maeve.Q_off_h = usrq('Total energy flux offset hi (-0.6): ',-0.6);
pars.maeve.Q_off_l = usrq('Total energy flux offset lo (5e4 m): ',5e4);
pars.maeve.E_amp_h = usrq('Char. energy peak hi (1e4 eV): ',1e4);
pars.maeve.E_amp_l = usrq('Char. energy peak lo (3e3 eV): ',3e3);
pars.maeve.E_wth_h = usrq('Char. energy width hi (2e4 m): ',2e4);
pars.maeve.E_wth_l = usrq('Char. energy width lo (2e5 m): ',2.2e5);
pars.maeve.E_gsl_h = usrq('Char. energy gradient scale length hi (0.4): ',0.4);
pars.maeve.E_gsl_l = usrq('Char. energy gradient scale length lo (0.1): ',0.1);
pars.maeve.K_amp = usrq('Line integrated current peak (0.3 A m^-1): ',0.3);
pars.maeve.J_wth = usrq('FAC sheet width (3e5 m): ',3e5);
pars.maeve.J_gsl = usrq('FAC sheet gradient scale length (0.1): ',0.1);
pars.maeve.F_amp = usrq('Flow peak (2000 m s^-1): ',2000);
pars.maeve.F_wth = usrq('Flow width (5e4 m): ',5e4);
pars.maeve.F_gsl = usrq('Flow gradient scale length (0.3): ',0.3);

[mapQit,mapE0it,mapJ,mapU,mapV] = MAEVE_map(X2,X3,IT,pars);

t=150;

figure(1)
pcolor(X2(:,:,t),X3(:,:,t),mapQit(:,:,t))
shading flat
colorbar

figure(2)
pcolor(X2(:,:,t),X3(:,:,t),mapJ(:,:,t))
shading flat
colorbar

figure(3)
hold on
plot(squeeze(X3(144/2,:,t)),squeeze(mapQit(144/2,:,t))./pars.maeve.Q_amp_h,'k')
plot(squeeze(X3(144/2,:,t)),squeeze(mapE0it(144/2,:,t))./pars.maeve.E_amp_h,'--k')
plot(squeeze(X3(144/2,:,t)),squeeze(mapJ(144/2,:,t))./(pars.maeve.K_amp/pars.maeve.J_wth),'r')
plot(squeeze(X3(144/2,:,t)),squeeze(mapU(144/2,:,t))./(pars.maeve.F_amp),'b')
plot(squeeze(X3(144/2,:,t)),squeeze(mapV(144/2,:,t))./(pars.maeve.F_amp),'--b')
hold off
legend(...
    'Qit / '+string(pars.maeve.Q_amp_h)+' erg/s/cm^2'...
    ,'E0it / '+string(pars.maeve.E_amp_h)+' eV'...
    ,'j_{||} / '+string(pars.maeve.K_amp/pars.maeve.J_wth)+' A/m^2'...
    ,'v_{East} / '+string(pars.maeve.F_amp)+' m/s'...
    ,'v_{North} / '+string(pars.maeve.F_amp)+' m/s'...
)
grid on

function [Qit,E0it,J,U,V] = MAEVE_map(x2,x3,it,pars)
    p = pars.maeve;
    x2 = x2 - p.driftE*it;
    x3 = x3 - p.driftN*it;
    c = (p.ctr_spn/2)*tanh(2*p.ctr_slp*x2/p.ctr_spn);
    dcdx = p.ctr_slp*sech(2*p.ctr_slp*x2/p.ctr_spn).^2;
%     s = sqrt(1+dcdx.^2);
    s = 1;
    b = bar(x2,p.bar_pos+p.bar_vel*it,p.bar_frc,p.bar_gsl); % loading bar
    d = (2-p.dim_frc*(1-tanh(2*(it-p.dim_del)/p.dim_tim)))/2; % dimming
    J_amp = p.K_amp/p.J_wth;
    Qit = (p.Q_amp_h-p.Q_amp_l)*d.*b.*...
         sheet(x3,c+p.Q_wth_l/2+p.Q_off_l+p.Q_off_h*p.Q_wth_l/2,p.Q_wth_h*s,p.Q_gsl_h)...
          +p.Q_amp_l*...
         sheet(x3,c+p.Q_wth_l/2+p.Q_off_l,p.Q_wth_l*s,p.Q_gsl_l);
    E0it = (p.E_amp_h-p.E_amp_l)*d.*b.*...
         sheet(x3,c+p.Q_wth_l/2+p.Q_off_l+p.Q_off_h*p.Q_wth_l/2,p.E_wth_h*s,p.E_gsl_h)...
          +p.E_amp_l*...
         sheet(x3,c+p.Q_wth_l/2+p.Q_off_l,p.E_wth_l*s,p.E_gsl_l); % h and l in same pos as Q
    J = J_amp*(...
         sheet(x3,c+p.J_wth/2,p.J_wth*s,p.J_gsl)...
        -sheet(x3,c-p.J_wth/2,p.J_wth*s,p.J_gsl)...
        );
    F = p.F_amp*(...
         sheet(x3,c-p.J_wth,p.F_wth*s,p.F_gsl)...
        -sheet(x3,c          ,p.F_wth*s,p.F_gsl)...
        +sheet(x3,c+p.J_wth,p.F_wth*s,p.F_gsl)...
        );
    U = F.*(1/sqrt(1+dcdx.^2));
    V = F.*(dcdx./sqrt(1+dcdx.^2));
end

function [v] = sheet(x3,pos,wdth,gsl)
    v = (tanh(2*(x3-pos+wdth/2)./(gsl*wdth))-tanh((x3-pos-wdth/2)./(gsl*wdth)))/2;
end

function [v] = bar(x2,pos,frac,gsl)
    v = (2-frac*(1-tanh(2*(x2-pos)/gsl)))/2;
end

function [v] = usrq(q,defv)
%     vi = input(q);
    vi='';
        if isempty(vi)
            v = defv;
        else
            v = vi;
        end
end