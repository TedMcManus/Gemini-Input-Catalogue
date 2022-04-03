function [E0,Q] = STEVE_particles_new(center,E0max,E0BG,Qmax,QBG,ref,outdir,mlatsig,force)

arguments
    center (1,1) double %Center of steep tanh shapefunction
    E0max (1,1) double %Max E0 value
    E0BG (1,1) double %Min E0 value
    Qmax (1,1) double %Max Q value
    QBG (1,1) double %Min Q value
    %NOTE - the backgrounds are all imposed as floors, not additively
    
    ref (1,1) string %Reference directory with config and grid files
    outdir (1,1) string %Output directory
    %Custom is './path/OSSE_particles'
    
    mlatsig (1,1) double = 0.15/2 %Distance from high to low (deg MLAT)
    force (1,1) int16 = 1 %Delete things without prompting the user (Default is yes)
end

ref=char(ref);
outdir=char(outdir);

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

%CREATE SOME SPACE FOR OUTPUT FILES
if ~exist('outdir','var')
    error('please specify output directory')
else
    a = char(outdir);
    if a(length(a))=='/'
        ind = length(a)-1;
        a = a(1:ind);
        outdir = char(a);
    end
    out = string(outdir);
    if ~exist(out,'dir')
        mkdir(out)
    elseif ~exist('force','var')||force==false
        x = input(strjoin(['Delete directory','"'+out+'"','recursively?','\n[Y/N]\n']),'s');
        if strcmp(x,"Y")
            rmdir(out,'s')
            mkdir(out)
        else
            error('Operation terminated')
        end
    elseif force==true
        rmdir(out,'s')
        mkdir(out)
    end
end


%READ IN THE SIMULATION INFORMATION (MEANS WE NEED TO CREATE THIS FOR THE SIMULATION WE WANT TO DO)
if ~exist('ymd0','var')
    dat=gemini3d.read.config(direcconfig);
    ymd0=dat.ymd;
    UTsec0=dat.UTsec0;
    tdur=dat.tdur;
    fprintf('Input config file loaded.\n');
end


%CHECK WHETHER WE NEED TO RELOAD THE GRID 
if (~exist('xg','var'))
    %WE ALSO NEED TO LOAD THE GRID FILE
    xg=gemini3d.read.grid([direcgrid,'/']);
    fprintf('Grid loaded.\n');
end


%Create a MLAT-MLON grid 
MLAT=90-squeeze(xg.theta(1,:,:))*180/pi;
MLON=squeeze(xg.phi(1,:,:))*180/pi;
llon=xg.lx(2);
llat=xg.lx(3);
mlon=MLON(:,1);    
mlat=MLAT(1,:);

%mlonmean=mean(mlon);
mlatmean=mean(mlat);


%TIME VARIABLE (SECONDS FROM SIMULATION BEGINNING)
tmin=0;
tmax=tdur;
%lt=tdur+1;
%time=linspace(tmin,tmax,lt)';
time=tmin:1:300;
lt=numel(time);

%SET UP TIME VARIABLES
ymd=ymd0;
UTsec=UTsec0+time;     %time given in file is the seconds from beginning of hour
UThrs=UTsec/3600;
expdate=cat(2,repmat(ymd,[lt,1]),UThrs(:),zeros(lt,1),zeros(lt,1));
t=datenum(expdate);


%CREATE THE PRECIPITAITON INPUT DATA
Qit=zeros(llon,llat,lt);
E0it=zeros(llon,llat,lt);

E0scale = E0max/2;
Qscale = Qmax/2;

%Create a steep tanh function that is high for high latitudes and low for
%low latitudes. Then, multiply the max values by the tanh shapefunction,
%and replace all values less than the background with the background value. 
[ind,~] = STEVE_centerpoint(center,MLAT,mlatmean,mlatsig);
shapefn = (tanh((MLAT-mlatmean+ind)/(mlatsig*2))+1);
for it=1:lt
    E0it(:,:,it)=(shapefn.*E0scale);
    E0now = squeeze(E0it(:,:,it));
    inds = find(E0now(:)<E0BG);
    E0now(inds) = E0BG;
    E0it(:,:,it) = E0now;
    %Qit(:,:,it)=((tanh((MLON-mlonmean)/(mlonsig/4))+1).*Qpk/2.25)+QBG;
    Qit(:,:,it)=(shapefn.*Qscale);
    Qnow = squeeze(Qit(:,:,it));
    inds = find(Qnow(:)<QBG);
    Qnow(inds) = QBG;
    Qit(:,:,it) = Qnow;
end

%%find the center
A = shapefn;
abval = abs(A-1);
val = min(abval(:));
inds = find(abval==val);
latctr = MLAT(inds);
fprintf(['Latitudinal center is at ',num2str(latctr(1)),'\n']);

E0=E0it;
Q=Qit;

