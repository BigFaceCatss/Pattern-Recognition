%start code for project #1: linear regression
%pattern recognition, CSE583/EE552
%Weina Ge, Aug 2008
%Christopher Funk, Jan 2018

addpath export_fig/

%load the data points
load data.mat

%plot the groud truth curve
figure(1)
clf
hold on;
xx = linspace(1,4*pi,npts);
yy = sin(.5*xx);
err = ones(size(xx))*0.3;
% plot the x and y color the area around the line by err (here the std)
h = shadedErrorBar(x,y,err,{'b-','color','b','LineWidth',2},0);
%plot the noisy observations
plot(x,t,'ro','MarkerSize',8,'LineWidth',1.5);
hold off; 
% Make it look good
grid on;
set(gca,'FontWeight','bold','LineWidth',2)
xlabel('x')
ylabel('t')

% Save the image into a decent resolution
% export_fig sampleplot -png -transparent -r150

%% Start your curve fitting program here
%Creat matrix X
m=9;
X=[];
for i=0:1:m
    X=[X,x'.^i];
end
%Regression using error minimization of SSE
%Calculate optimal W1
W1=X\(t');
Y1=(X*W1)';
%plot the regression curve and the noisy observation
figure(2)
clf
hold on;
plot(x,y,'b-','color','g','LineWidth',2)
plot(x,Y1,'b-','color','b','LineWidth',2);
plot(x,t,'ro','MarkerSize',8,'LineWidth',1.5);
hold off;
% Make it look good
grid on;
title('M=9         N=50');
% legend('origin','curve regression')
set(gca,'FontWeight','bold','LineWidth',2)
xlabel('x')
ylabel('t')

%Regression using error minimization with the regularization term
%Create coefficient lambda
lambda=exp(-18);
%Calculate optimal W2
W2=(X'*X+lambda*eye(m+1))\(X'*t');
Y2=(X*W2)';
%plot the regression curve and the noisy observation
figure(3)
clf
hold on;
plot(x,y,'b-','color','g','LineWidth',2)
plot(x,Y2,'b-','color','b','LineWidth',2);
plot(x,t,'ro','MarkerSize',8,'LineWidth',1.5);
hold off;
% Make it look good
grid on;
title('M=9         N=10             ln\lambda=20');
% legend('origin','curve regression')
set(gca,'FontWeight','bold','LineWidth',2)
xlabel('x')
ylabel('t')

%Regression using the ML (maximal likelihood) estimator of the Bayesian approach
%Calculate B=\beta_ML^-1
B=((X*W1-t')'*(X*W1-t'))/npts;
beta=1/B
%Calculate the regression function
Y3 = (X*W1)'+sqrt(B).*randn(1,npts);
%Plot the regression curve and the noisy observation
figure(4)
clf
hold on;
err = ones(size(x))*sqrt(B);
% plot the x and y color the area around the line by err (here the std)
h = shadedErrorBar(x,Y3,err,{'b-','color','r','LineWidth',2},0);
plot(x,y,'b-','color','g','LineWidth',2)
plot(x,Y3,'b-','color','b','LineWidth',2);
plot(x,t,'ro','MarkerSize',8,'LineWidth',1.5);
hold off;
% Make it look good
grid on;
title('M=9         N=50');
% legend('origin','curve regression')
set(gca,'FontWeight','bold','LineWidth',2)
xlabel('x')
ylabel('t')

%Regression using the MAP (maximum a posteriori) estimator of the Bayesian approach
%create x_new
xn = linspace(1,4*pi,npts);
%Create \alpha
a=5*10^(36)
%calculate mean(u) and variance(Sigama)
U=[];
Sigama=[];
phix=[];
phixn=[]
phix1=[];
phix2=zeros(m+1);
phixt=zeros(m+1,1);
for k=1:1:npts
    for n=0:1:m
        phixn=[phixn;xn(1,k)^n];
    end
    for j=1:1:npts
        for i=0:1:m
        phix=[phix;x(1,j)^i]
        end
    phix1=phix*phixn'
    phix2=phix2+phix1;
    phixt=phixt+phix*t(1,j);
    phix=[];
    end
S=inv(a*eye(m+1)+beta*phix2);
u=beta*phixn'*S*phixt;
s=B+phixn'*S*phixn;
U=[U,u];
Sigama=[Sigama,s];
phixn=[];
end
%Calculate the regression function
Y4 = U+sqrt(Sigama).*randn(1,npts);
%Plot the regression curve and the noisy observation
figure(5)
clf
hold on;
err = ones(size(x)).*sqrt(Sigama);
% plot the x and y color the area around the line by err (here the std)
h = shadedErrorBar(x,Y4,err,{'b-','color','r','LineWidth',2},0);
plot(x,y,'b-','color','g','LineWidth',2)
plot(x,Y4,'b-','color','b','LineWidth',2);
plot(x,t,'ro','MarkerSize',8,'LineWidth',1.5);
hold off;
% Make it look good
grid on;
title('M=9         N=10');
% legend('origin','curve regression')
set(gca,'FontWeight','bold','LineWidth',2)
xlabel('x')
ylabel('t')