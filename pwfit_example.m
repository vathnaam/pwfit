%% pwfit Examples
%% Univariate Function (Single Fit)
clear
clc
x = linspace(1,3,10);
y = exp(x);
pwfitTest = pwfit(x,y,'poly2');
plot(pwfitTest); % Plot the fit
% Making the plot look pretty
legend('Data','\fontname{Courier}pwfit \fontname{Helvetica}curve','Location','northwest');
ylabel('y'); xlabel('x');
set(gcf,'position',[0 0 900 500]);
set(gca,'LineWidth',2.5,'FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2.5,'MarkerSize',30);
%% Univariate Function (Piecewise Fit)
close all
pwfitTest = pwfit(x,y,'poly2',[1,2,3]);
% subplot(5,1,[1,3]);
plot(pwfitTest); % Plot the fit
% Making the plot look pretty
legend('Data','1^{st} \fontname{Courier}pwfit \fontname{Helvetica}curve','\fontname{Courier}pwfit \fontname{Helvetica}breakpoint','2^{nd} \fontname{Courier}pwfit \fontname{Helvetica}curve','Location','northwest');
xlabel('x'); ylabel('y');
set(gcf,'position',[0 0 900 500]);
set(gca,'LineWidth',2.5,'FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2.5,'MarkerSize',30);

%% Univariate Function (Piecewise Fit with optimized breakpoint)
close all
pwfitTest = pwfit(x,y,'poly2',[1,2,3], 'optimized');
plot(pwfitTest); % Plot the fit
% Making the plot look pretty
legend('Data','1^{st} \fontname{Courier}pwfit \fontname{Helvetica}curve','\fontname{Courier}pwfit \fontname{Helvetica}breakpoint','2^{nd} \fontname{Courier}pwfit \fontname{Helvetica}curve','Location','northwest');
xlabel('x'); ylabel('y');
set(gcf,'position',[0 0 900 500]);
set(gca,'LineWidth',2.5,'FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2.5,'MarkerSize',30);

%% Bivariate Function (Piecewise Fit with optimized breakpoint)
clc
clear
close all
x = linspace(2,5,10);
y = linspace(1,3,10);
[X,Y] = meshgrid(x,y);
Z =  -sin(Y)^2 + exp(X);
pwfitTest = pwfit(X,{Y,Z},'poly21',linspace(x(1),x(end),3),'optimized');
plot(pwfitTest); % Plot the fit
% Making the plot look pretty
xlabel('x');ylabel('y');zlabel('z');
set(gcf,'position',[0 0 900 600]);
set(gca,'LineWidth',2.5,'FontSize',15);
set(findall(gca, 'Type', 'Line'),'LineWidth',2.5,'MarkerSize',20);
