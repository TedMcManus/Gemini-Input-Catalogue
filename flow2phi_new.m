% reconstructing 2D ion flow field from distributed flow dataset
% provides potential from input flow field on gemini grid
% uses parameterized, longitudinally inclined, gaussian ridges as pseudo-basis functions
% see reconstructor_readme.pdf for more details
% contact: jules.van.irsel.GR@dartmouth.edu

% dependencies:
% matlab R2020b
% optimization toolbox
% parallel computing toolbox (optional)
% image processing toolbox

function [recon_phi,recon_v2,recon_v3,P] = flow2phi_new(...
    gemini_v2... % bag of eastward flow vectors
    ,gemini_v3... % bag of northward flow vectors
    ,cad2... % cadence over longitude
    ,cad3... % cadence over latitude
    ,xg... % gemini grid
    ,outdir... 
    ,numf... % number of basis functions used in reconstruction (32)
    ,showplt... % show plot (false)
    ,saveplt... % save plot (false)
    ,showboundary... % show boundary plot (false)
    ,usepar... % use parallel computation (false)
    ,isVerbose...
    ,it...
    ,maxiter... % lsqcurvefit option value for "MaxIterations"
    ,maxfuneval... % lsqcurvefit option value for "MaxFunctionEvaluations"
    )
    arguments
        gemini_v2 (:,:)
        gemini_v3 (:,:)
        cad2 (1,1) int16
        cad3 (1,1) int16
        xg (1,1) struct
        outdir (1,1) string
        numf (1,1) int16 = 32
        showplt (1,1) logical = false
        saveplt (1,1) logical = false
        showboundary (1,1) logical = false
        usepar (1,1) logical = false
        isVerbose (1,1) logical = false
        it (1,1) int16 = 1
        maxiter (1,1) int16 = 400
        maxfuneval (1,1) int32 = 1e3
    end
    global boundaryc Bmag Nm

    outdir=char(outdir);
    %% grid metadata
    lx1 = xg.lx(1); lx2 = xg.lx(2); lx3 = xg.lx(3);
    Bmag = abs(mean(xg.Bmag,'all'));
    [gemini_x2,gemini_x3] = ndgrid(xg.x2(3:end-2),xg.x3(3:end-2)); % E-N distance from grid
    
    %% find arc boundary from v2 input
    if isVerbose
        fprintf('Fitting arc boundary...\n')
    end
    edges = edge(gemini_v2,'Sobel',0.1);
    boundary = zeros(1,lx2);
    for ix2 = 1:lx2
        boundary(ix2) = gemini_x3(1,find(edges(ix2,:),1,'first')); % grab southern-most detected edges
    end

    boundaryf = fit(double(gemini_x2(:,1))*1e-5,double(boundary')*1e-5,'a+b*tanh(c+d*x)','Start',[0 1 0 1]);
    boundaryc = coeffvalues(boundaryf).*[1e5,1e5,1,1e-5];
    arc_bound_x3 = boundaryc(1)+boundaryc(2)*tanh(boundaryc(3)+boundaryc(4)*gemini_x2(:,1));

    if showboundary
        figure(1) % plot arc boundary
        pcolor(gemini_x2,gemini_x3,gemini_v2); shading flat; hold on
        plot(gemini_x2(:,1),boundary,'b','LineWidth',3.0)
        plot(gemini_x2(:,1),arc_bound_x3,'r','LineWidth',3.0); hold off
    end

    %% reconstruct
    if isVerbose
        fprintf('Reconstruction setup...\n')
    end
    if isequal(size(gemini_v2),[lx2 lx3]) && isequal(size(gemini_v3),[lx2 lx3])
        bag_x2 = gemini_x2(1:cad2:end,1:cad3:end);
        bag_x3 = gemini_x3(1:cad2:end,1:cad3:end);
        bag_v2 = gemini_v2(1:cad2:end,1:cad3:end);
        bag_v3 = gemini_v3(1:cad2:end,1:cad3:end);
    else
        error('Bag of vectors does not match grid size, [' + string(lx2) + ',' + string(lx3) + '].')
    end
    
    xdata = double(reshape(cat(3,bag_x2,bag_x3),[numel(bag_x2),2])); % reshape to len(bag) x 2 array
    ydata = double(reshape(cat(3,bag_v2,bag_v3),[numel(bag_v2),2]));
    xdata(isnan(ydata(:,1)),:) = []; % remove nans from bag
    ydata(isnan(ydata(:,1)),:) = [];
    if isVerbose
        fprintf('Number of fitting elements: ' + string(size(ydata,1)) + '\n')
    end

    %% optimization
    if isVerbose
        fprintf('Reconstructing flow...\n')
    end
    if usepar 
        %parpool('local'); 
    end
    tic
    Nm = numf; % number of basis functions
    P_0 = [linspace(-5,5,Nm); ones(1,Nm); ones(1,Nm); ones(1,Nm)]; % initial parameter matrix: [x3pos (100 km), x3sig (100 km), x2inc (kV / 1000 km), x2amp (kV)]
    options = optimoptions('lsqcurvefit'...
        ,'Algorithm','levenberg-marquardt'...
        ,'UseParallel',usepar...
        ,'StepTolerance',1e-6...
        ,'MaxIterations',maxiter...
        ,'MaxFunctionEvaluations',2*numel(P_0)*maxfuneval...
        );
    P = lsqcurvefit(@(P, xdata) F(P,xdata),P_0,xdata,ydata,[],[],options); % optimized parameter matrix
    recon_time = toc;
    %delete(gcp('nocreate'))
    disp(options)
    if isVerbose
        fprintf('Reconstruction time: ' + string(recon_time) + ' seconds\n')

    %% creating output + plot arrays
        fprintf('Constructing output arrays...\n')
    end
    plt_xdata = double(reshape(cat(3,gemini_x2,gemini_x3),[numel(gemini_x2),2]));
    recon_v = F(P,xdata); % s/c locations only
    recon_vt = F(P,plt_xdata);
    recon_v2 = single(reshape(recon_vt(:,1),[lx2,lx3]));
    recon_v3 = single(reshape(recon_vt(:,2),[lx2,lx3]));
    recon_phi = single(phi(P,gemini_x2,gemini_x3));
    recon_phi = recon_phi - mean(recon_phi,'all');
    
    %% Calculating error
    if isVerbose
        fprintf('Determining goodness of fit...\n')
    end
    reg_buf = 0; % potential region buffer
    minx2 = min(bag_x2(:))-reg_buf; maxx2 = max(bag_x2(:))+reg_buf;
    minx3 = min(bag_x3(:))-reg_buf; maxx3 = max(bag_x3(:))+reg_buf;
    reg = (gemini_x2>minx2 & gemini_x2<maxx2 & gemini_x3>minx3 & gemini_x3<maxx3); % region of interest around s/c
    reg_d = double(reg); % for plotting purposes
    reg_d(reg_d==0) = nan;
    error_a_v2 = (recon_v2-gemini_v2).^2; % determine square differences in region of interest
    error_a_v3 = (recon_v3-gemini_v3).^2;
    error_p_v2 = ((recon_v2-gemini_v2)./max(gemini_v2(:))).^2; % determine square percent error in region
    error_p_v3 = ((recon_v3-gemini_v3)./max(gemini_v3(:))).^2;
    if isVerbose
        fprintf('Root median square difference in v2 around s/c region is ' + string(sqrt(median(error_a_v2(reg),'all'))) + ' m/s\n')
        fprintf('Root median square difference in v3 around s/c region is ' + string(sqrt(median(error_a_v3(reg),'all'))) + ' m/s\n')
        fprintf('Root median square percent error in v2 around s/c region is ' + string(100*sqrt(median(error_p_v2(reg),'all'))) + ' %%\n')
        fprintf('Root median square percent error in v3 around s/c region is ' + string(100*sqrt(median(error_p_v3(reg),'all'))) + ' %%\n')
    end
    if ~strcmp(outdir,'')
        if ~exist([outdir,'/reconstructor'],'dir')
            mkdir([outdir,'/reconstructor']);
        end
        fid=fopen([outdir,'/reconstructor/reconstructor_error.txt'],'a');
        fprintf(fid,['--------- IT=',num2str(it),' ---------\n\n']);
        fprintf(fid,'Root median square difference in v2 around s/c region is ' + string(sqrt(median(error_a_v2(reg),'all'))) + ' m/s\n');
        fprintf(fid,'Root median square difference in v3 around s/c region is ' + string(sqrt(median(error_a_v3(reg),'all'))) + ' m/s\n');
        fprintf(fid,'Root median square percent error in v2 around s/c region is ' + string(100*sqrt(median(error_p_v2(reg),'all'))) + ' %%\n');
        fprintf(fid,'Root median square percent error in v3 around s/c region is ' + string(100*sqrt(median(error_p_v3(reg),'all'))) + ' %%\n');
    end
    %% plotting results
    spw = 0.24;
    sph = 0.25;
    spws = (1-3*spw)/4;
    sphs = (1-2*sph)/4;
    spho = 0.07;
    crange_phi = 1.1.*[min(recon_phi(:)),max(recon_phi(:))].*1e-3;
    crange_v2 = 1.1.*[min(gemini_v2(:)),max(gemini_v2(:))];
    crange_v3 = 1.1.*[min(gemini_v3(:)),max(gemini_v3(:))];
    crange_dv2 = [0 20];
    fz = 16;
if showplt
        figure(2)
        set(gcf,'units','inches','OuterPosition',[1 1 14/1.5 10/1.5])

        sp1 = subplot(2,3,1);
        pcolor(gemini_x2,gemini_x3,gemini_v2); hold on
        quiver(bag_x2,bag_x3,bag_v2,bag_v3,'r'); hold off
        shading flat
        cb1 = colorbar('southoutside');
        caxis(crange_v2)
        set(sp1,'position',[spws 3*sphs+sph+spho spw sph])
        set(sp1,'fontsize',fz)
        title('model flow east')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb1.Label.String = 'flow east [m/s]';
        
        sp2 = subplot(2,3,2);
        pcolor(gemini_x2,gemini_x3,gemini_v3); hold on
        quiver(bag_x2,bag_x3,bag_v2,bag_v3,'r'); hold off
        shading flat
        if(~all(gemini_v3(:)==0))
            cb2 = colorbar('southoutside');
            caxis(crange_v3);
        end
        set(sp2,'position',[2*spws+spw 3*sphs+sph+spho spw sph])
        set(sp2,'fontsize',fz)
        title('model flow north')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb2.Label.String = 'flow north [m/s]';

        sp3 = subplot(2,3,3);
        pcolor(gemini_x2,gemini_x3,recon_phi.*1e-3); hold on
        quiver(bag_x2,bag_x3,bag_v2,bag_v3,'r'); hold off
        shading flat
        cb3 = colorbar('southoutside');
        caxis(crange_phi)
        set(sp3,'position',[3*spws+2*spw 3*sphs+sph+spho spw sph])
        set(sp3,'fontsize',fz)
        title('recon. potential')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb3.Label.String = 'potential [kV]';

        sp4 = subplot(2,3,4);
        pcolor(gemini_x2,gemini_x3,recon_v2); hold on
        quiver(xdata(:,1),xdata(:,2),recon_v(:,1),recon_v(:,2),'r'); hold off
        shading flat
        cb4 = colorbar('southoutside');
        caxis(crange_v2)
        set(sp4,'position',[spws sphs+spho spw sph])
        set(sp4,'fontsize',fz)
        title('recon. flow east')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb4.Label.String = 'flow east [m/s]';

        sp5 = subplot(2,3,5);
        pcolor(gemini_x2,gemini_x3,recon_v3); hold on
        quiver(xdata(:,1),xdata(:,2),recon_v(:,1),recon_v(:,2),'r'); hold off
        shading flat
        if(~all(crange_v3==0))
            cb5 = colorbar('southoutside');
            caxis(crange_v3);
        end
        set(sp5,'position',[2*spws+spw sphs+spho spw sph])
        set(sp5,'fontsize',fz)
        title('recon. flow north')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb5.Label.String = 'flow north [m/s]';

        sp6 = subplot(2,3,6);
        pcolor(gemini_x2,gemini_x3,100.*sqrt(error_p_v2)); hold on
        alpha 0.2
        pcolor(gemini_x2,gemini_x3,100.*sqrt(error_p_v2).*reg_d); hold off
        shading flat
        cb6 = colorbar('southoutside');
        caxis(crange_dv2)
        set(sp6,'position',[3*spws+2*spw sphs+spho spw sph])
        set(sp6,'fontsize',fz)
        title('flow east error')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb6.Label.String = 'rms error [%]';

        if saveplt
            saveas(gcf,[outdir,'/reconstructor/reconstructor_it=',num2str(it),'.png']);
        end
    else
        figure(2)
        set(gcf,'units','inches','OuterPosition',[1 1 14/1.5 10/1.5])

        sp1 = subplot(2,3,1);
        pcolor(gemini_x2,gemini_x3,gemini_v2); hold on
        quiver(bag_x2,bag_x3,bag_v2,bag_v3,'r'); hold off
        shading flat
        cb1 = colorbar('southoutside');
        caxis(crange_v2)
        set(sp1,'position',[spws 3*sphs+sph+spho spw sph])
        set(sp1,'fontsize',fz)
        title('model flow east')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb1.Label.String = 'flow east [m/s]';
        
        sp2 = subplot(2,3,2);
        pcolor(gemini_x2,gemini_x3,gemini_v3); hold on
        quiver(bag_x2,bag_x3,bag_v2,bag_v3,'r'); hold off
        shading flat
        if(~all(gemini_v3(:)==0))
            cb2 = colorbar('southoutside');
            caxis(crange_v3);
        end
        set(sp2,'position',[2*spws+spw 3*sphs+sph+spho spw sph])
        set(sp2,'fontsize',fz)
        title('model flow north')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb2.Label.String = 'flow north [m/s]';

        sp3 = subplot(2,3,3);
        pcolor(gemini_x2,gemini_x3,recon_phi.*1e-3); hold on
        quiver(bag_x2,bag_x3,bag_v2,bag_v3,'r'); hold off
        shading flat
        cb3 = colorbar('southoutside');
        caxis(crange_phi)
        set(sp3,'position',[3*spws+2*spw 3*sphs+sph+spho spw sph])
        set(sp3,'fontsize',fz)
        title('recon. potential')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb3.Label.String = 'potential [kV]';

        sp4 = subplot(2,3,4);
        pcolor(gemini_x2,gemini_x3,recon_v2); hold on
        quiver(xdata(:,1),xdata(:,2),recon_v(:,1),recon_v(:,2),'r'); hold off
        shading flat
        cb4 = colorbar('southoutside');
        caxis(crange_v2)
        set(sp4,'position',[spws sphs+spho spw sph])
        set(sp4,'fontsize',fz)
        title('recon. flow east')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb4.Label.String = 'flow east [m/s]';

        sp5 = subplot(2,3,5);
        pcolor(gemini_x2,gemini_x3,recon_v3); hold on
        quiver(xdata(:,1),xdata(:,2),recon_v(:,1),recon_v(:,2),'r'); hold off
        shading flat
        if(~all(crange_v3==0))
            cb5 = colorbar('southoutside');
            caxis(crange_v3);
        end
        set(sp5,'position',[2*spws+spw sphs+spho spw sph])
        set(sp5,'fontsize',fz)
        title('recon. flow north')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb5.Label.String = 'flow north [m/s]';

        sp6 = subplot(2,3,6);
        pcolor(gemini_x2,gemini_x3,100.*sqrt(error_p_v2)); hold on
        alpha 0.2
        pcolor(gemini_x2,gemini_x3,100.*sqrt(error_p_v2).*reg_d); hold off
        shading flat
        cb6 = colorbar('southoutside');
        caxis(crange_dv2)
        set(sp6,'position',[3*spws+2*spw sphs+spho spw sph])
        set(sp6,'fontsize',fz)
        title('flow east error')
        xlabel('distance east [m]')
        ylabel('distance north [m]')
        cb6.Label.String = 'rms error [%]';

        if saveplt
            saveas(gcf,[outdir,'/reconstructor/reconstructor_it=',num2str(it),'.png']);
        end
    end
    %% functions
    % phi basis function
    function phi = phi(P,x2,x3)
        x3pos = P(1,:)*1e5;
        x3sig = P(2,:)*1e5;
        x2inc = P(3,:)*1e-4;
        x2int = P(4,:)*1e3;
        phi = 0;
        for m = 1:Nm
%             phi = phi + (x2inc(m).*x2 + x2int(m)).*exp(-((x3-x3pos(m)-polyval(boundaryc,(x2(:,1))))./x3sig(m)).^2);
            b = boundaryc(1) + boundaryc(2)*tanh(boundaryc(3)+boundaryc(4)*x2(:,1));
            phi = phi + (x2inc(m).*x2 + x2int(m)).*exp(-((x3-x3pos(m)-b)./x3sig(m)).^2);
        end
    end

    % lsqcurvefit fitting function
    function v = F(P,xdata)
        Ni = size(xdata,1); % number of vectors in bag
        x3pos = P(1,:)*1e5; % latitudinal positions [m] (P entries are near unity)
        x3sig = P(2,:)*1e5; % latitudinal widths [m]
        x2inc = P(3,:)*1e-4; % longitudenal slope of potential ridge [V/m]
        x2amp = P(4,:)*1e3; % central amplitude of potential ridge [V]
        E = zeros(Ni,2);
        v = zeros(Ni,2);
        for i = 1:Ni % iterate through bag of vectors
            x2 = xdata(i,1);
            x3 = xdata(i,2);
%             b = 0;
%             db = 0;
%             for j = 1:Nj % define boundary position
%                 b = b + boundaryc(j)*x2^(Nj+1-j); % polynomial of order Nj
%                 db = db + (Nj+1-j)*boundaryc(j)*x2^(Nj-j); % derivative of polynomial
%             end
            b = boundaryc(1) + boundaryc(2)*tanh(boundaryc(3)+boundaryc(4)*x2);
            db = boundaryc(2)*boundaryc(4)*sech(boundaryc(3)+boundaryc(4)*x2)^2;
            % caluclate the elements of -grad(phi) = -sum_m grad(phi_m) (see doccumentation for details)
            for m = 1:Nm % iterate through number of basis functions
                expf = exp(-((x3-x3pos(m)-b)/x3sig(m))^2);
                E(i,1) = E(i,1) + (-x2inc(m)-(2/x3sig(m)^2)*(x2inc(m)*x2+x2amp(m))*(x3-x3pos(m)-b)*db)*expf;
                 E(i,2) = E(i,2) + (2/x3sig(m)^2)*(x2inc(m)*x2+x2amp(m))*(x3-x3pos(m)-b)*expf;
            end
            v(:,1) = -E(:,2)./Bmag;
            v(:,2) =  E(:,1)./Bmag;
        end
    end

end
