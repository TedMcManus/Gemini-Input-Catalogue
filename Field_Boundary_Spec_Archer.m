function [mlatsig_out,Vmaxx1it,width_out,start_out,top_out,bottom_out,wcase] = Field_Boundary_Spec_Archer(wtarget,MLAT,MLON,mlonmean,mlonsig,mlatctr,Jpk,index)
%This function is kinda gross and I don't think anyone will want to use it
%after me. That said, I'lll write some cursory comments

difference = 10; %maximum precision tolerance
for mlatsig = 0:.001:1.2 %this is the fit parameter
    %first, calculate the potential for the given fit parameter
    %Vmaxx1it(:,:)=Jpk.*exp(-(MLON-mlonmean).^2/2/mlonsig^2).*exp(-(MLAT-mlatctr-1.5*mlatsig).^2/2/mlatsig^2);
    Vmaxx1it(:,:)=Jpk.*exp(-(MLON-mlonmean).^2/2/mlonsig^2).*exp(-(MLAT-mlatctr).^2/2/mlatsig^2);
    %then index through the potential. When we detect a slope increase,
    %save this as the "start," save the place where the slope stops
    %increaasing as the "top," and save the place where the slope goes back
    %to zero as the "bottom." This is a quick and dirty numerical method to
    %estimate the size of the "peak" in flow
    for i =1:127
        if Vmaxx1it(index,i)>10&&~exist('start','var') 
            start = i;
        end
        if Vmaxx1it(index,i)>0&&Vmaxx1it(index,i)>Vmaxx1it(index,i+1)&&~exist('top','var')
            top = i;
        end
        if Vmaxx1it(index,i)<0&&Vmaxx1it(index,i)<Vmaxx1it(index,i+1)&&exist('top','var')&&~exist('bottom','var')
            bottom = i;
        end
    end
    %calculate width of the flow feature using the stored values. We have a
    %"start to top" and "top to bottom" distance
    if exist('top','var')&&exist('start','var')
        if start~=1
            width = abs(MLAT(index,top)-MLAT(index,start));
            widthcase = 'S_T';
        end
    end
    if exist('top','var')&&exist('bottom','var')&&~exist('width','var')
        width = abs(MLAT(index,top)-MLAT(index,bottom));
        widthcase = 'T_B';
    end
    if ~exist('width','var')
        width = 11;
        widthcase = 'null';
    end
    dif = abs(wtarget-width); %check if the fit is good for the desired width

    %if the fit is good, save the fit as "mlatsig_out" and change the fit parameter to signify we have a new best fit 
    if dif < difference&&exist('start','var')&&exist('top','var')&&strcmp(widthcase,'S_T')
        difference = dif;
        mlatsig_out = mlatsig;
        width_out = width;
        start_out = MLAT(index,start);
        top_out = MLAT(index,top);
        wcase = 'S_T';
    elseif dif < difference&&exist('top','var')&&exist('bottom','var')&&strcmp(widthcase,'T_B')
        difference = dif;
        mlatsig_out = mlatsig;
        width_out = width;
        top_out = MLAT(index,top);
        bottom_out = MLAT(index,bottom);
        wcase = 'T_B';
    end

    clear start top bottom width widthcase

end

mlatsig = mlatsig_out; %this is the final fit parameter

%now calculate the output flow and show some plots
Vmaxx1it(:,:)=Jpk.*exp(-(MLON-mlonmean).^2/2/mlonsig^2).*exp(-(MLAT-mlatctr).^2/2/mlatsig^2);
%Vmaxx1it(:,:)=Vmaxx1it(:,:)-Jpk.*exp(-(MLON-mlonmean).^2/2/mlonsig^2).*exp(-(MLAT-mlatctr+1.5*mlatsig).^2/2/mlatsig^2);
%Vmaxx1it=-Vmaxx1it;
plot(MLAT(index,:),Vmaxx1it(index,:))
if strcmp(wcase,'S_T')
    bottom_out = 'null';
elseif strcmp(wcase,'T_B')
    start_out = 'null';
end
%surf(MLON,MLAT,Vmaxx1it(:,:))
