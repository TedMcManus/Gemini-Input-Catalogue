function pot_out = flow_spec(ref,wtarget,displace,mlonsig,index,flagdirich,ExBg,EyBg,lbound,ubound,steps,xg,vtarg)


%REFERENCE GRID TO USE
if exist('ref','var')
    a = char(ref);
    if a(length(a))=='/'
        ind = length(a)-1;
        a = a(1:ind);
        ref = char(a);
    end
    direcconfig = ref;
    direcgrid = ref;
else
    error('reference sim not found')
end


%READ IN THE SIMULATION INFORMATION (MEANS WE NEED TO CREATE THIS FOR THE SIMULATION WE WANT TO DO)
if (~exist('ymd0','var'))
    dat=gemini3d.read.config(direcconfig);
    ymd0=dat.ymd;
    UTsec0=dat.UTsec0;
    tdur=dat.tdur;
    fprintf('Input config.dat file loaded.\n');
end



%CHECK WHETHER WE NEED TO RELOAD THE GRID (SO THIS ALREADY NEEDS TO BE MADE, AS WELL)
if (~exist('xg','var'))
  %WE ALSO NEED TO LOAD THE GRID FILE
  xg=gemini3d.read.grid(direcgrid);
  lx1=xg.lx(1); lx2=xg.lx(2); lx3=xg.lx(3);
  fprintf('Grid loaded.\n');
end


MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
MLON=squeeze(xg.phi(1,:,:))*180/pi;
llon=xg.lx(2);
llat=xg.lx(3);
mlon=MLON(:,1);    
mlat=MLAT(1,:);
mlonmean=mean(mlon);
mlatmean=mean(mlat);

%TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin=0;
tmax=tdur;
time=tmin:10:tmax;
lt=numel(time);


%SET UP TIME VARIABLES
ymd=ymd0;
UTsec=UTsec0+time;     %time given in file is the seconds from beginning of hour
UThrs=UTsec/3600;
expdate=cat(2,repmat(ymd,[lt,1]),UThrs(:),zeros(lt,1),zeros(lt,1));
t=datenum(expdate);


%CREATE DATA FOR BACKGROUND ELECTRIC FIELDS
Exit=zeros(llon,llat,lt);
Eyit=zeros(llon,llat,lt);
if ~exist('ExBg','var')
    ExBg = 0;
end
if ~exist('EyBg','var')
    EyBg = 0;
end
for it=1:lt
  Exit(:,:,it)=ExBg*ones(llon,llat);   %V/m
  Eyit(:,:,it)=EyBg*ones(llon,llat);
end



%CREATE DATA FOR BOUNDARY CONDITIONS FOR POTENTIAL SOLUTION
if strcmp(flagdirich,'potential')||strcmp(flagdirich,'Potential')
    flagdirich = 1;
end
if strcmp(flagdirich,'FAC')
    flagdirich = 0;
end
if ~exist('flagdirich','var') %if 0 data is interpreted as FAC, else we interpret it as potential
    flagdirich = 1; 
end
Vminx1it=zeros(llon,llat,lt);
Vmaxx1it=zeros(llon,llat,lt);
Vminx2ist=zeros(llat,lt);
Vmaxx2ist=zeros(llat,lt);
Vminx3ist=zeros(llon,lt);
Vmaxx3ist=zeros(llon,lt);


if ~exist('Jpk','var')
    Jpk=400;
end
if ~exist('mlonsig','var')
    mlonsig=30;
end
%MLATSIG is how beefy in latitude. It is determined by the boundary spec
%MLONSIG is how beefy in longitude. Just make it like 20
%Displace is how much latitude the arc spans. 1-2 range
mlatctr=mlatmean+displace*tanh((MLON-mlonmean)/(mlonsig));

[mlatsig,~,~,~,~,~] = Field_Boundary_Spec(wtarget,MLAT,MLON,mlonmean,mlonsig,mlatctr,Jpk,index);

pot_out = Vestimate(lbound,ubound,steps,MLAT,MLON,mlonmean,mlonsig,mlatsig,mlatctr,xg,vtarg);