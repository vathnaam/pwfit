classdef pwfit
    %PWFIT    Fit a single or a piecewise curve/surface to data
    %
    %   pwfitObject = pwfit(x,y,fitType) create a fit to the data x and y with
    %   the model fitType.
    %
    %   pwfitObject = pwfit(x,{y,z},fitType) create a surface fit to the data
    %   x, y, and z (all must be matrices).
    %
    %   pwfitObject = pwfit(x,y,fitType,breaks) create a piecewise fit to the
    %   data x and y at the specified breaks.
    %
    %   pwfitObject = pwfit(x,y,fitType,breaks,mode) create a piecewise fit
    %   where the break locations are either 'fixed' or 'optimized'.
    %
    %   pwfitObject = pwfit(x,y,fitType,breaks,mode,wt) allows the user to
    %   specify the weight of the spacing between each break.
    %
    %   Valid Model Structures:
    %
    %   - Curve Fitting
    %       - poly1:  a*x + b
    %       - poly2:  a*x^2 + b*x + c
    %       - poly3:  a*x^3 + b*x^2 + c*x + d
    %       - poly4:  a*x^4 + b*x^3 + c*x^2 + d*x + e
    %       - poly5:  a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f
    %
    %   - Surface Fitting
    %       - poly11: a*x + b*y + c
    %       - poly12: a*y^2 + b*x*y + c*x + d*y + e
    %       - poly13: a*y^3 + b*y^2 + c*x*y^2 + d*x*y + e*x + f*y + g
    %       - poly21: a*x^2 + b*x*y + c*x + d*y + e
    %       - poly22: a*x^2 + b*y^2 + c*x*y + d*x + e*y + f
    %       - poly23: a*y^3 + b*x^2 + c*y^2 + d*x^2*y + e*x*y^2 + f*x*y + g*x + h*y + i
    %       - poly31: a*x^3 + b*x^2 + c*x^2*y + d*x*y + e*x + f*y + g
    %       - poly32: a*x^3 + b*x^2 + c*y^2 + d*x^2*y + e*x*y^2 + f*x*y + g*x + h*y + i
    %       - poly33: a*x^3 + b*y^3 + c*x^2 + d*y^2 + e*x^2*y + f*x*y^2 + g*x*y + h*x + i*y + j
    
    % Copyright 2019 Vathna Am
    
    properties (SetAccess = private)
        theta
        breaks
        stats
        data
        fun
        grad
        optiObj
        is3D
        isPW
    end
    
    methods
        %Constructor
        function fitObj = pwfit(xData,yData,fitType,breaks,mode,wt)
            % Create a pwfit object and fit a single or piecewise curve/surface to the input data.
            
            %Check input arguments
            if nargin < 6 || isempty(wt)
                wt = 0.01;
            end
            if nargin < 5 || isempty(mode)
                mode = 'fixed';
            end
            
            %Pre-processing data
            [fitObj.data.xd,fitObj.data.yd,fitObj.data.zd] = fitObj.ppData(xData,yData);
            if(~isempty(fitObj.data.zd))
                fitObj.is3D = true;
            else
                fitObj.is3D = false;
            end
            
            %Get the fit model
            [fitObj.fun,fitObj.grad,fitObj.stats.modelStructure] = fitObj.fitmodel(fitObj,fitType);
            
            if nargin < 4 %if breaks are not provided, it is not a piecewise fit
                fitObj.isPW = false;
                fitObj.breaks = [min(min(fitObj.data.xd));max(max(fitObj.data.xd))];
                [~, fitObj.theta, fitObj.stats.sse] = fitObj.sumAllSse(fitObj.data.xd,fitObj.data.yd,fitObj.data.zd,fitObj.fun,fitObj.grad,fitObj.breaks,1);
            else %piecewise fit if breaks are provided
                fitObj.breaks = breaks;
                fitObj.isPW = true;
                switch mode
                    case 'optimized' %optimize the break locations
                        [x,~] = fitObj.optibreak(fitObj,fitObj.data.xd,fitObj.data.yd,fitObj.data.zd,fitObj.fun,fitObj.grad,fitObj.breaks(:),wt); %Optimize the breaks with constrained fit
                        fitObj.breaks = [min(fitObj.breaks);x;max(fitObj.breaks)];
                        [~, fitObj.theta, fitObj.stats.sse] = fitObj.sumAllSse(fitObj.data.xd,fitObj.data.yd,fitObj.data.zd,fitObj.fun,fitObj.grad,fitObj.breaks(:),wt); %Getting coefficients
                    case 'fixed' %fixed break locations
                        [~, fitObj.theta, fitObj.stats.sse] = fitObj.sumAllSse(fitObj.data.xd,fitObj.data.yd,fitObj.data.zd,fitObj.fun,fitObj.grad,fitObj.breaks(:),1); %Only constrained fit
                end
            end
            %By convention coefficients are row vectors
            fitObj.theta = reshape(fitObj.theta,[length(fitObj.fun(1,1)),length(fitObj.breaks)-1]); fitObj.theta = fitObj.theta';
            fitObj.theta = round(fitObj.theta,10,'significant');
        end
        
        function funEval = feval(fitObj,varargin)
            %Evaluate the fitted function at a specified parameter(s) or using the original data
            
            funEval = [];
            if(isempty(varargin))%if parameters are not provided, use the data provided
                xd = fitObj.data.xd(:);
                yd = fitObj.data.yd(:);
                if(~fitObj.isPW)
                    funEval = fitObj.fun(xd,yd)*fitObj.theta(1,:)';
                elseif(fitObj.isPW)
                    for i = 1:length(fitObj.breaks)-1
                        idx = find(xd <= fitObj.breaks(i+1));
                        funEval = [funEval;fitObj.fun(xd(idx),yd(idx))*fitObj.theta(i,:)'];
                        xd(idx) = []; yd(idx) = [];
                    end
                end
            elseif(~isempty(varargin))%evaluate the function at the specified parameters
                xd = varargin{1}; xd = xd(:);
                if(length(varargin)==2), yd = varargin{2}; yd = yd(:); end
                
                if(~fitObj.isPW)
                    if(fitObj.is3D)
                        funEval = fitObj.fun(xd,yd)*fitObj.theta(1,:)';
                    elseif(~fitObj.is3D)
                        funEval = fitObj.fun(xd,zeros(size(xd)))*fitObj.theta(1,:)';
                    end
                elseif(fitObj.isPW)
                    if(round(min(xd),3,'significant') < round(fitObj.breaks(1),3,'significant') || round(max(xd),3,'significant') > round(fitObj.breaks(end),3,'significant'))
                        error('The x variable does not meet the piecewise condition');
                    end
                    for i = 1:length(fitObj.breaks)-1
                        idx = find(xd <= fitObj.breaks(i+1)+1e-3);
                        if(~isempty(idx))
                            if(fitObj.is3D)
                                funEval = [funEval;fitObj.fun(xd(idx),yd(idx))*fitObj.theta(i,:)'];
                            elseif(~fitObj.is3D)
                                funEval = [funEval;fitObj.fun(xd(idx),size(xd(idx)))*fitObj.theta(i,:)'];
                            end
                            xd(idx) = [];
                            if(length(varargin)==2), yd(idx) = []; end
                        end
                    end
                end
            end
        end
        
        function pwSubFun = print(fitObj,varargin)
            %Print Print the fit model (or all the sub-functions) as a string
            
            %Checking if the user specify the function precision
            if ~isempty((strfind(varargin{end},'%')))
                precision = varargin{end};
                varLength = length(varargin)-1;
            else
                precision = '%1.15d'; %default precision
                varLength = length(varargin);
            end
            
            [rowLength,colLength] = size(fitObj.theta);
            pwSubFun = cell(rowLength,1);
            for i = 1:rowLength
                str = fitObj.stats.modelStructure;
                for j = colLength:-1:1
                    t = ['t' int2str(j)];
                    str = strrep(str,t,num2str(fitObj.theta(i,j),precision));
                end
                pwSubFun{i} = str;
            end
            if ~isempty(varargin) %changing the variable name
                for i = 1:rowLength
                    if varLength < 2
                        if ischar(varargin{1})
                            pwSubFun{i} = regexprep(pwSubFun{i},'xd',varargin(1:varLength));
                        end
                    elseif varLength == 2
                        if ischar(varargin{1}) && ischar(varargin{2})
                            pwSubFun{i} = regexprep(pwSubFun{i},{'xd','yd'},varargin(1:varLength));
                        end
                    end
                end
            end
        end
        
        function plot(objFit)
            % Plot Plot the piecewise fit
            if (~objFit.is3D)
                %Plot 2D fit
                plot(objFit.data.xd,objFit.data.yd,'.k','MarkerSize',15);hold on;
                for i = 1:length(objFit.breaks)-1
                    xlin = linspace(objFit.breaks(i),objFit.breaks(i+1),100); xlin = xlin(:);
                    ylin = objFit.fun(xlin)*objFit.theta(i,:)'; ylin = ylin(:);
                    plot(xlin,ylin);
                    if i < length(objFit.breaks)-1
                        %Plot breaks
                        plot(objFit.breaks(i+1),objFit.fun(objFit.breaks(i+1))*objFit.theta(i,:)','.r','MarkerSize', 20);
                    end
                end
                if(length(objFit.breaks) > 2)
                    vline(objFit.breaks(2:end-1),':r');
                end
                
                xlim([min(objFit.data.xd(:)),max(objFit.data.xd(:))]);
                ylim([min(objFit.data.yd(:)),max(objFit.data.yd(:))]);
                hold off
            elseif (objFit.is3D)
                %Plotting 3D fit
                plot3(objFit.data.xd,objFit.data.yd,objFit.data.zd,'.k','MarkerSize',13); hold on;
                
                [numRow,~] = size(objFit.data.xd); % Find the number of rows in the xdata matrix
                
                for i = 1:length(objFit.breaks)-1
                    idx = find((objFit.breaks(i) <= objFit.data.xd) & (objFit.data.xd <= objFit.breaks(i+1)));
                    x = objFit.data.xd(idx);
                    y = objFit.data.yd(idx);
                    % reshape the xdata and ydata into matrices
                    x = reshape(x,numRow,[]);
                    y = reshape(y,numRow,[]);
                    
                    
                    % Generating extra data points at the break so it can
                    % be plotted
                    if  i > 1
                        lambda = (x(:,1)-objFit.breaks(i).*ones(numRow,1))./(x(:,1)-objFit.data.xd(idx(1)-numRow:idx(1)-1)');
                        x = [objFit.breaks(i).*ones(numRow,1),x];
                        prevColY = y(:,1).*(1-lambda) + objFit.data.yd(idx(1)-numRow:idx(1)-1)'.*lambda;
                        y = [prevColY,y];
                        
                        
                    end
                    
                    if i < length(objFit.breaks)-1
                        lambda = (x(:,end)-objFit.breaks(i+1).*ones(numRow,1))./(x(:,end)-objFit.data.xd(idx(end)+1:idx(end)+numRow)');
                        x = [x, objFit.breaks(i+1).*ones(numRow,1)];
                        nextColY = y(:,end).*(1-lambda) + objFit.data.yd(idx(end)+1:idx(end)+numRow)'.*lambda;
                        y = [y,nextColY];
                    end
                    
                    if i~=1
                        % Plot the break(s)
                        plot3(objFit.breaks(i)*ones(numRow,1),y(:,1),objFit.feval(objFit.breaks(i)*ones(numRow,1),y(:,1)),'-r','linewidth',3);
                    end
                    
                    z = feval(objFit,x,y);
                    z = reshape(z,numRow,[]); % Evaluate the piecewise function
                    hs = surf(x,y,z,'LineStyle','none'); % Plot the piecewise function
                    shading(gca,'interp');
                    colormap winter
                    hs.FaceAlpha = 0.8;
                end
                xlim([min(objFit.data.xd(:)),max(objFit.data.xd(:))]);
                ylim([min(objFit.data.yd(:)),max(objFit.data.yd(:))]);
                hold off;
            end
            grid on;
        end
        
        function ploterr(objFit)
            %Plot error
            if (~objFit.is3D)
                % 2D error
                plot(objFit.data.xd,zeros(1,length(objFit.data.xd)),'-.','color',0.8*[1 1 1]); hold on;
                error = objFit.data.yd(:)-feval(objFit);
                plot(objFit.data.xd,error,'color',[0 0.4470 0.7410]); hold off;
                ylabel('Error');
                xlim([min(objFit.data.xd),max(objFit.data.xd)]);
                ylim([min(error),max(error)]);
            elseif (objFit.is3D)
                % 3D error
                error = objFit.data.zd(:) - feval(objFit);
                [C,h] = contourf(objFit.data.xd,objFit.data.yd,reshape(error,size(objFit.data.xd)),5); hold on;
                %Only plot some of the x and y points
                sub2 = plot(objFit.data.xd,objFit.data.yd,'.','MarkerSize',15/2,'color', 0.8*[1,1,1]);
                title('Error');
                view([-0, 90]);
                colorbar;
                clabel(C,h);
                w = warning ('off','all');
                cmap=cbrewer('div','RdBu',100);
                maxErr = max(abs(error));
                caxis([-maxErr maxErr]);
                warning(w);
                colormap(cmap)
                xlim([min(objFit.data.xd(:)),max(objFit.data.xd(:))]);
                ylim([min(objFit.data.yd(:)),max(objFit.data.yd(:))]);
            end
        end
    end
    
    
    methods (Static, Access=private)
        function [xd,yd,zd] = ppData(xdata,ydata)
            %Pre-processing input data
            if iscell(ydata)%Check if 3D fit
                if(length(ydata) ~= 2)
                    error('For 3D fit, ydata must be {ydata,zdata}');
                end
                xd = xdata; yd = ydata{1}; zd = ydata{2};
                %Check if matrices are supplied
                if(isvector(xd)||isvector(yd)||isvector(zd))
                    error('For 3D fit, your xdata, ydata and zdata must be matrices');
                end
                %Check if the size of all matrices are the same
                if(numel(xd) ~= numel(yd) || numel(xd) ~= numel(zd) || numel(yd) ~= numel(zd))
                    error('xdata, ydata and zdata must be the same size');
                end
            elseif ~iscell(ydata)%Check if 2D fit
                xd = xdata; xd = xd(:); yd = ydata; yd = yd(:); zd = [];
                %Check if the size of all vectors are the same
                if(numel(xd) ~= numel(yd))
                    error('xdata, ydata and zdata must be the same size');
                end
            end
        end
        
        function [ffun,gfun,str] = fitmodel(fitObj,fitType)
            if (~fitObj.is3D)
                %2D fit
                % Setup curve fitting function
                switch(fitType)
                    case('poly1')
                        %poly1: t1*xd + t2
                        ffun = @(x,y) [x,ones(size(x))];
                        gfun = @(x,y) [zeros(size(x)),zeros(size(x))];
                        str = 't1*xd + t2';
                        
                    case('poly2')
                        %poly2: t1*xd^2 + t2*xd + t3
                        ffun = @(x,y) [x.^2,x,ones(size(x))];
                        gfun = @(x,y) [2.*x,ones(size(x)),zeros(size(x))];
                        str = 't1*xd^2 + t2*xd + t3';
                        
                    case('poly3')
                        %poly3: t1*xd^3 + t2*xd^2 + t3*xd + t4
                        ffun = @(x,y) [x.^3,x.^2,x,ones(size(x))];
                        gfun = @(x,y) [3.*x.^2,2.*x,ones(size(x)),zeros(size(x))];
                        str = 't1*xd^3 + t2*xd^2 + t3*xd + t4';
                        
                    case('poly4')
                        %poly4: t1*xd^4 + t2*xd^3 + t3*xd^2 + t4*xd + t5
                        ffun = @(x,y) [x.^4,x.^3,x.^2,x,ones(size(x))];
                        gfun = @(x,y) [4.*x.^3,3.*x.^2,2.*x,ones(size(x)),zeros(size(x))];
                        str = 't1*xd^4 + t2*xd^3 + t3*xd^2 + t4*xd + t5';
                        
                    case('poly5')
                        %poly5: t1*xd^5 + t2*xd^4 + t3*xd^3 + t4*xd^2 + t5*xd + t6
                        ffun = @(x,y) [x.^5,x.^4,x.^3,x.^2,x,ones(size(x))];
                        gfun = @(x,y) [5.*x.^4,4.*x.^3,3.*x.^2,2.*x,ones(size(x)),zeros(size(x))];
                        str = 't1*xd^5 + t2*xd^4 + t3*xd^3 + t4*xd^2 + t5*xd + t6';
                        
                    case('power2') %not working
                        ffun = @(x,y) [log(x),ones(size(x))];
                        gfun = @(x,y) [x*(log(x) - 1),zeros(size(x))];
                        str = 'log(x.^t1) + t2';
                        
                    otherwise
                        error('Unknown 2D fit type');
                end
                
            elseif (fitObj.is3D)
                %3D fit
                % Setup surface fitting function
                switch(fitType)
                    case('poly11')
                        %poly11: t1*xd + t2*yd + t3
                        ffun = @(x,y) [x,y,ones(size(x))];
                        gfun = @(x,y) [zeros(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),zeros(size(x)),zeros(size(x))];
                        str = 't1*xd + t2*yd + t3';
                        
                    case('poly12')
                        %poly12 = t1*yd^2 + t2*xd*yd + t3*yd + t4*xd + t5
                        ffun = @(x,y) [y.^2,x.*y,y,x,ones(size(x))];
                        %gfun = @(x,y) [zeros(size(x)),y,zeros(size(x)),ones(size(x)),zeros(size(x));
                        %2.*y,x,ones(size(x)),zeros(size(x)),zeros(size(x))];
                        gfun = @(x,y) [zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x))];
                        str = 't1*yd^2 + t2*xd*yd + t3*yd + t4*xd + t5';
                        
                    case('poly13')
                        %poly13 = t1*yd^3 + t2*yd^2 + t3*xd*yd^2 + t4*xd*yd + t5*xd + t6*yd + t7
                        ffun = @(x,y) [y.^3,y.^2,x.*y.^2,x.*y,x,y,ones(size(x))];
                        %gfun = @(x,y) [zeros(size(x)),zeros(size(x)),y.^2,y,ones(size(x)),zeros(size(x)),zeros(size(x));
                        %3.*y.^2,2.*y,2.*x.*y,x,zeros(size(x)),ones(size(x)),zeros(size(x))];
                        gfun = @(x,y) [zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x)),zeros(size(x))];
                        str = 't1*yd^3 + t2*yd^2 + t3*xd*yd^2 + t4*xd*yd + t5*xd + t6*yd + t7';
                        
                    case('poly21')
                        %poly21: t1*xd^2 + t2*xd*yd + t3*xd + t4*yd + t5
                        ffun = @(x,y) [x.^2,x.*y,x,y,ones(size(x))];
                        gfun = @(x,y) [2.*x,y,ones(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),x,zeros(size(x)),ones(size(x)),zeros(size(x))];
                        str = 't1*xd^2 + t2*xd*yd + t3*xd + t4*yd + t5';
                        
                    case('poly22')
                        %poly22: t1*xd^2 + t2*yd^2 + t3*xd*yd + t4*xd + t5*yd + t6
                        ffun = @(x,y) [x.^2,y.^2,x.*y,x,y,ones(size(x))];
                        gfun = @(x,y) [2.*x,zeros(size(x)),y,ones(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),2.*y,x,zeros(size(x)),ones(size(x)),zeros(size(x))];
                        str = 't1*xd^2 + t2*yd^2 + t3*xd*yd + t4*xd + t5*yd + t6';
                        
                    case('poly23')
                        %poly23: t1*yd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9
                        ffun = @(x,y) [y.^3,x.^2,y.^2,(x.^2).*y,x.*(y.^2),x.*y,x,y,ones(size(x))];
                        gfun = @(x,y) [zeros(size(x)),2.*x,zeros(size(x)),2.*x.*y,y.^2,y,ones(size(x)),zeros(size(x)),zeros(size(x));
                            3.*(y.^2),zeros(size(x)),2.*y,x.^2,2.*x.*y,x,zeros(size(x)),ones(size(x)),zeros(size(x))];
                        str = 't1*yd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9';
                        
                    case('poly31')
                        %poly31: 	t1*xd^3 + t2*xd^2 + t3*xd^2*yd + t4*xd*yd + t5*xd + t6*yd + t7
                        ffun = @(x,y) [x.^3,x.^2,(x.^2).*y,x.*y,x,y,ones(size(x))];
                        gfun = @(x,y) [3.*(x.^2),2.*x,2.*x.*y,y,ones(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),zeros(size(x)),x.^2,x,zeros(size(x)),ones(size(x)),zeros(size(x))];
                        str = '	t1*xd^3 + t2*xd^2 + t3*xd^2*yd + t4*xd*yd + t5*xd + t6*yd + t7';
                        
                    case('poly32')
                        %poly32: t1*xd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9
                        ffun = @(x,y) [x.^3,x.^2,y.^2,(x.^2).*y,x.*(y.^2),x.*y,x,y,ones(size(x))];
                        gfun = @(x,y) [3.*(x.^2),2.*x,zeros(size(x)),2.*x.*y,y.^2,y,ones(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),zeros(size(x)),2.*y,x.^2,2.*y.*x,x,zeros(size(x)),ones(size(x)),zeros(size(x))];
                        str = 't1*xd^3 + t2*xd^2 + t3*yd^2 + t4*xd^2*yd + t5*xd*yd^2 + t6*xd*yd + t7*xd + t8*yd + t9';
                        
                    case('poly33')
                        %poly33: t1*xd^3 + t2*yd^3 + t3*xd^2 + t4*yd^2 + t5*xd^2*yd + t6*xd*yd^2 + t7*xd*yd + t8*xd + t9*yd + t10
                        ffun = @(x,y) [x.^3,y.^3,x.^2,y.^2,(x.^2).*y,x.*(y.^2),x.*y,x,y,ones(size(x))];
                        gfun = @(x,y) [3.*(x.^2),zeros(size(x)),2.*x,zeros(size(x)),2.*x.*y,y.^2,y,ones(size(x)),zeros(size(x)),zeros(size(x));
                            zeros(size(x)),3.*(y.^2),zeros(size(x)),2.*y,x.^2,2.*y.*x,x,zeros(size(x)),ones(size(x)),zeros(size(x))];
                        str = 't1*xd^3 + t2*yd^3 + t3*xd^2 + t4*yd^2 + t5*xd^2*yd + t6*xd*yd^2 + t7*xd*yd + t8*xd + t9*yd + t10';
                        
                    otherwise
                        error('Unknown 3D fit type');
                end
            end
        end
        
        function [x,feval,Opt] = optibreak(fitObj,xd,yd,zd,ffun,gfun,breaks,wt)
            %Obtaining SSE from lsqlin and use it as objective function
            obj = @(x)  fitObj.sumAllSse(xd,yd,zd,ffun,gfun,[min(breaks);x;max(breaks)],wt);
            
            %Bounds
            lb = min(xd(:))*ones((length(breaks)-2),1);
            ub = max(xd(:))*ones((length(breaks)-2),1);
            
            %Inequality constraints
            A = zeros((length(breaks)-2)-1,length(breaks)-2);
            for i = 1:(length(breaks)-2)-1
                A(i,i) = 1;
                A(i,i+1) = -1;
            end
            b = zeros((length(breaks)-2)-1,1);
            
            %Initial guess
            x0 = breaks(2:end-1);
            
            %OPTI setting
            nloptOpt = nloptset('algorithm','AUGLAG','subalgorithm','LN_COBYLA');
            opts = optiset('solver','nlopt','maxiter',1e10,'maxfeval',1e6,'display','off','tolrfun',1e-3,'tolafun',1e-3,'solverOpts',nloptOpt);
            
            %Create OPTI Object
            Opt = opti('fun',obj,'bounds',lb,ub,'x0',x0,'ineq',A,b,'options',opts);
            
            %Solve
            [x,feval,exitflag,info] = solve(Opt);
            
            if(exitflag ~= 1)
                disp(['      Exitflag: ' num2str(exitflag)]);
                disp(info);
                fprintf(2,'The optimal breakpoint solution was not found. Please check your input arguments.');
            end
            
        end
        
        function [obj, coeffs, sse] = sumAllSse(xd,yd,zd,ffun,gfun,breaks,wt)
            breaks = sort(breaks);
            V = diagVandBlk(xd,yd,ffun,breaks);
            [Aeq,beq] = eqCont(yd,ffun,gfun,breaks);
            
            ws = warning('off','all');
            options = optimoptions('lsqlin','Algorithm','interior-point','display','off');
            if isempty(zd)
                [coeffs,resnorm] = lsqlin(V,yd(:),[],[],Aeq,beq,[],[],[],options);
            elseif ~isempty(zd)
                [coeffs,resnorm] = lsqlin(V,zd(:),[],[],Aeq,beq,[],[],[],options);
            end
            warning(ws);
            
            sse = resnorm;
            d = sum(1./diff(breaks));%use to penalize the spacing between two adjacent breaks
            obj = sse + d*wt; %cost function for optimizing breaks (outer optimization)
            
            
            function V = diagVandBlk(xd,yd,ffun,breaks)
                %Construct the diagonal Vandermonde matrix block
                xd = xd(:);
                yd = yd(:);
                vand = cell(1,length(breaks)-1);
                for i = 1:length(breaks)-2
                    idx = find(xd<round(breaks(i+1),10));
                    vand{i} = ffun(xd(idx),yd(idx));
                    if isempty(vand{i}) %if breakpoints are the same
                        vand{i} = ffun(zeros([],1),zeros([],1));
                    end
                    xd(idx) = []; yd(idx) = [];
                end
                vand{end} = ffun(xd,yd);
                
                V = blkdiag(vand{:});
                
            end
            function [Aeq,beq] = eqCont(yd,ffun,gfun,breaks)
                %Equality constraints
                bp = breaks(2:end-1);
                k = 1;
                j = 1;
                [r,~] = size(gfun(1,1));
                if r == 1
                    Aeq = zeros(2*length(bp),length(ffun(1,1))*(length(breaks)-1));
                    for i = 1:length(breaks)-2
                        Aeq(k,j:(j-1)+2*length(ffun(1))) = [ffun(bp(i)) -ffun(bp(i))]; k = k + 1;
                        Aeq(k,j:(j-1)+2*length(ffun(1))) = [gfun(bp(i)) -gfun(bp(i))]; k = k + 1;
                        j = j + length(ffun(1));
                    end
                    beq = zeros(2*length(bp),1);
                elseif r == 2 %(Bad form, need further work)
                    Aeq = zeros(3*length(bp)*length(yd),length(ffun(1,1))*(length(breaks)-1));
                    for i = 1:length(breaks)-2
                        Aeq(k:i*length(yd),j:(j-1)+2*length(ffun(1,1)))  = [ffun(bp(i)*ones(length(yd),1),yd(:,1)) -ffun(bp(i)*ones(length(yd),1),yd(:,1))]; k = k + length(yd);
                        j = j + length(ffun(1,1));
                    end
                    j = 1;
                    for i = 1:length(breaks)-2
                        Aeq(k:(k-1) + 2*length(yd),j:(j-1)+2*length(ffun(1,1))) = [gfun(bp(i)*ones(length(yd),1),yd(:,1)) -gfun(bp(i)*ones(length(yd),1),yd(:,1))]; k = k + 2*length(yd);
                        j = j + length(ffun(1,1));
                    end
                    beq = zeros(3*length(bp)*length(yd),1);
                end
            end
            
        end
    end
end
