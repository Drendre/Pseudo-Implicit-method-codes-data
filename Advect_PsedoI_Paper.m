clear; fprintf('\n')
Mh=90; %How many different timestep we want
fact=1.1386; % h is divided with this
D=1; a=2; K=1; x0=1; xf=2; ti=0; tf=0.1; TIME=tf-ti; %Diffusion coefficient, lenght of the space and time interval, 
N=101; Dx=(xf-x0)/(N-1); %space discretization
MaxD=zeros(Mh,20);Axhstep=zeros(Mh,1); rfun=zeros(Mh,1);%%Errors and axis for the timestep-sizes
h=TIME/10; % intitial time step size
xaxis=zeros(N,1); for i=1:N xaxis(i)=x0+(i-1)*Dx; end  %x axis
u0=zeros(N,1); uex=zeros(N,1); %%initial values and exact solution
E=D/Dx^2; E2=1/Dx/2; %Coefficients for the discretization of space derivatives
for i=1:N 
 u0(i)=exp(-xaxis(i)); 
 uex(i)=exp((D+a-K)*tf-xaxis(i));   %%Initial and exact final results
 end
for ih=1:Mh  %% Big loop for the timestep, ih is the Serial Number of run, 1 for the largest timestep; 
T = ceil(TIME/h-0.1)+1;
tax=zeros(T,1); %%physical time axis
for t=1:T
      tax(t)=ti+(t-1)*h; 
end
if(h>0.00000001 && h<0.0002) rfun(ih)=2000000*h*h; end
UUP=zeros(N,1);UPI=zeros(N,1);UTEMP=zeros(N,1);   %%Final U and temporary results
Axhstep(ih)=h; %%Time step size axis
r=h*E; A=a*h*E2; A2=2*A;  
%%%%%%%%%%%%%%%%%%%%%%%%%% Orig UPFD
UUP(:)=u0(:); %% Initialization
tic
for t=1:T-1
    for i=2:N-1  %%% UPFD 
    UTEMP(i)=(UUP(i)+2*A*UUP(i-1)+r*(UUP(i-1)+UUP(i+1)))/(1+2*A+2*r+K*h);
    end
    UTEMP(1)=exp((D+a-K)*tax(t+1)-x0);
    UTEMP(N)=exp((D+a-K)*tax(t+1)-xf);
 UUP(:)=UTEMP(:);
end
toc
MaxD(ih,1)=max(abs(UUP-uex));%% errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5
%% Pesudo-Imp, 1-2 symm
p=1;lam=1/2/p;UPI=u0; 
tic
for t=1:T-1
        %% first stage
    for i=2:N-1  %%%UPFD
     UTEMP(i)=((1-2*p*r*(1-lam))*UPI(i)+p*r*(UPI(i-1)+UPI(i+1))-p*A*(UPI(i+1)-UPI(i-1)))/(1+2*p*r*lam+p*h*K);
    end
     UTEMP(1)=exp((D+a-K)*(tax(t+1)-h/2)-x0);
     UTEMP(N)=exp((D+a-K)*(tax(t+1)-h/2)-xf);
      UTEMP=lam*UTEMP+(1-lam)*UPI;
          %%  2. stage
   for i=2:N-1  %%% theta=1/2
         UPI(i)=((1-r)*UPI(i)+r*(UTEMP(i-1)+UTEMP(i+1))-A*(UTEMP(i+1)-UTEMP(i-1))+K*h*(UTEMP(i)-UPI(i)))/(1+r+h*K);
   end
   UPI(1)=exp((D+a-K)*tax(t+1)-x0);UPI(N)=exp((D+a-K)*tax(t+1)-xf);
end
toc
MaxD(ih,2)=max(abs(UPI-uex)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 6
%% Pesudo-Imp, 1 asymm imp, 2 symm
UPI=u0; 
tic
for t=1:T-1
        %% first stage
    for i=2:N-1  %%%UPFD
     UTEMP(i)=((1-r*(1/lam-1))*UPI(i)+p*r*(UPI(i-1)+UPI(i+1))+p*2*A*UPI(i-1))/(1+2*p*r*lam+2*p*A+p*h*K);
    end
     UTEMP(1)=exp((D+a-K)*(tax(t+1)-h/2)-x0);
     UTEMP(N)=exp((D+a-K)*(tax(t+1)-h/2)-xf);
     UTEMP=lam*UTEMP+(1-lam)*UPI;
          %%  2. stage
   for i=2:N-1  %%% theta=1/2
         UPI(i)=((1-r)*UPI(i)+r*(UTEMP(i-1)+UTEMP(i+1))-A*(UTEMP(i+1)-UTEMP(i-1))+K*h*(UTEMP(i)-UPI(i)))/(1+r+h*K);
   end
   UPI(1)=exp((D+a-K)*tax(t+1)-x0);UPI(N)=exp((D+a-K)*tax(t+1)-xf);
end
toc
MaxD(ih,3)=max(abs(UPI-uex));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 7
%% Pesudo-Imp, 1 symm, 2 asym imp
UPI=u0; 
tic
for t=1:T-1
        %% first stage
    for i=2:N-1  
     UTEMP(i)=((1-r*(1/lam-1))*UPI(i)+p*r*(UPI(i-1)+UPI(i+1))-p*A*(UPI(i+1)-UPI(i-1)))/(1+2*p*r*lam+p*h*K);
    end
     UTEMP(1)=exp((D+a-K)*(tax(t+1)-h/2)-x0);
     UTEMP(N)=exp((D+a-K)*(tax(t+1)-h/2)-xf);
     UTEMP=lam*UTEMP+(1-lam)*UPI;
          %%  2. stage
   for i=2:N-1  
         UPI(i)=((1-r)*UPI(i)+r*(UTEMP(i-1)+UTEMP(i+1))+2*A*UTEMP(i-1)+K*h*(UTEMP(i)-UPI(i)))/(1+r+2*A+h*K);
   end
   UPI(1)=exp((D+a-K)*tax(t+1)-x0);UPI(N)=exp((D+a-K)*tax(t+1)-xf);
end
toc
MaxD(ih,4)=max(abs(UPI-uex));    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 8    
%% Pesudo-Imp, 1 and 2 asymm imp
UPI=u0; 
tic
for t=1:T-1
        %% first stage
    for i=2:N-1  %%%UPFD
     UTEMP(i)=((1-r*(1/lam-1))*UPI(i)+p*r*(UPI(i-1)+UPI(i+1))+p*2*A*UPI(i-1))/(1+2*p*r*lam+2*p*A+p*h*K);
    end
     UTEMP(1)=exp((D+a-K)*(tax(t+1)-h/2)-x0);
     UTEMP(N)=exp((D+a-K)*(tax(t+1)-h/2)-xf);
      UTEMP=lam*UTEMP+(1-lam)*UPI;
          %%  2. stage
   for i=2:N-1  %%% theta=1/2
         UPI(i)=((1-r)*UPI(i)+r*(UTEMP(i-1)+UTEMP(i+1))+2*A*UTEMP(i-1)+K*h*(UTEMP(i)-UPI(i)))/(1+r+2*A+h*K);
   end
   UPI(1)=exp((D+a-K)*tax(t+1)-x0);UPI(N)=exp((D+a-K)*tax(t+1)-xf);
end
toc
MaxD(ih,5)=max(abs(UPI-uex));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=h/fact;

end
 green         = [0 0.9 0.3]; 
    dark_grey     = [0.5,0.5,0.5]; 
    orange        = [1 0.5 0];
    dark_red      = [0.75, 0.0780, 0.1840]; 
    dark_yellow   = [0.9290, 0.6940, 0.1250];
    yellow   = [0.990, 0.8940, 0.250];  
    dark_green    = [0, 0.45, 0];
    dark_purple   = [0.55, 0.0840, 0.5560];
    light_blue     = [0.3, 0.65, 1];
    dark_dark_grey = [0.2, 0.2, 0.25]; 
    Magenta	 =    [1, 0, 1];
    MilkChocolate	= [165/255, 42/255, 42/255];
    Black		= [0, 0, 0];
    Red	= [1, 0, 0];

	%% figure 1:
  figure('Name', 'Errors as a fuction of time step h');   

	plot(Axhstep(:),MaxD(:,1), "--o", 'Color',Red, 'LineWidth', 2, 'MarkerSize',6); %% UPFD 
	hold on;
    plot(Axhstep(:),MaxD(:,2), ':',  'Color', orange, 'LineWidth', 5); %% cncn p=1/3
	hold on;
    plot(Axhstep(:),MaxD(:,3), '-x',  'Color', light_blue, 'LineWidth', 2.7); %% cncn p=2/3
	hold on;
    plot(Axhstep(:),MaxD(:,4), '-.O',  'Color', dark_yellow, 'LineWidth', 4); %% cncn p=1/2
	hold on;
	plot(Axhstep(:),MaxD(:,5), '-bd',  'Color', MilkChocolate, 'LineWidth', 1.1,'MarkerSize',6); %% cncn p=3/4
	hold on;
    plot(Axhstep(:),rfun, '--',  'Color', dark_grey, 'LineWidth', 1, 'MarkerSize',8); %% reference
	hold on;
	hold off;
	
XLabel = xlabel('Time step size \it{h}');
set(XLabel, 'FontSize', 16);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('Errors');
set(YLabel, 'FontSize', 16); 
set(YLabel, 'FontWeight', 'bold');

% Get handle to current axes.
Ax = gca;
Ax.XAxis.FontSize = 18;
Ax.YAxis.FontSize = 18;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({'UPFD A1','PI A5', 'PI A6', 'PI A7','PI A8'},'Location','southeast', 'FontSize', 14);
set(Legend, 'FontSize', 14); 
set(Legend, 'FontWeight', 'bold');
set(gca,'xScale','log');
set(gca,'yScale','log');
grid on;

