%% Test against exact solution in 1D, 3 stage CCL
%% loop for the timesteps 
clear; fprintf('\n')
Mh=14; %How many different timestep we want 
D=1; sig=3; %% diffusion and radiation coefficient
x0=-1; xf=1; N=101; Dx=(xf-x0)/(N-1); %%Space-discretization
ti=0.5; tf=1; TIME=tf-ti; % Time discretzation
h=TIME/10;  %%initial time step size 
dp=1/2; p0=1/2; Np=3; pf=p0+dp*Np; %Np=ceil((pf-p0)/dp+0.1); %How many different parameters we want
MaxD=zeros(Mh,Np+1);Axhstep=zeros(Mh,1); rfun=zeros(Mh,1);%% Errors, axis for the timestep-sizes, reference function
xaxis=zeros(N,1); for i=1:N xaxis(i)=x0+(i-1)*Dx; end  % x axis
u0=zeros(N,1); uex=zeros(N,1); %%initial values and exact solution
E=D/(Dx)^2; %%matrix elements
for i=1:N 
u0(i)=ti*exp(xaxis(i)-ti);  %% initial function
uex(i)=tf*exp(xaxis(i)-tf);  %% exact final results
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ih=1:Mh  %% Big loop for the timestep, ih is the Serial Number of run, 1 for the largest timestep; 
T = ceil(TIME/h-0.1)+1;
taxis=zeros(T,1);B1=zeros(T,1);B2=zeros(T,1); %%physical time
for t=1:T
   taxis(t)=ti+(t-1)*h; %% time axis
 B1(t)=taxis(t)*exp(x0-taxis(t)); B2(t)=taxis(t)*exp(xf-taxis(t)); %% Left and right boundary  
end
if(h>0.0000004 && h<0.01) rfun(ih)=50000*h*h; end
UU=zeros(N,1); UPI=zeros(N,1); UTEMP=zeros(N,1);  %%Final U-s and temporary results
Axhstep(ih)=h; %%Timestep axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UPFD starts here
UU(:)=u0(:); r=h*E;
tic
for t=1:T-1
        %% first stage
    for i=2:N-1  %%%UPFD
        kt=sig*(taxis(t)+h)^4*exp(4*xaxis(i)-4*(taxis(t)+h))+exp(xaxis(i)-(taxis(t)+h));
     UTEMP(i)=(UU(i)+r*(UU(i-1)+UU(i+1))+h*kt)/(1+2*r+2*h+h*sig*UU(i)^3);
    end
     UTEMP(1)=B1(t+1);UTEMP(N)=B2(t+1);
     UU=UTEMP;
end
toc
%%Max errors
MaxD(ih,1)=max(abs(UU-uex)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pseudo-Iplicit starts here
for ip=1:Np
UPI(:)=u0(:); r=h*E;
p=p0+(ip-1)*dp; 
lam=1/2/p  
the=1-lam;
tic
for t=1:T-1
        %% first stage, fractional time step p*h
    for i=2:N-1  %%% Theta = the
        kt=sig*(taxis(t)+p*h)^4*exp(4*xaxis(i)-4*(taxis(t)+p*h))+exp(xaxis(i)-(taxis(t)+p*h));
     UTEMP(i)=((1-2*p*r*the)*UPI(i)+p*r*(UPI(i-1)+UPI(i+1))+p*h*kt)/(1+2*p*r*(1-the)+2*p*h+p*h*sig*UPI(i)^3);
    end
     UTEMP=lam*UTEMP+(1-lam)*UPI;
     UTEMP(1)=B1(t+1);UTEMP(N)=B2(t+1);
          %%  2. stage
   for i=2:N-1  %%% Theta= 1/2
       kt=sig*taxis(t+1)^4*exp(4*xaxis(i)-4*taxis(t+1))+exp(xaxis(i)-taxis(t+1));
         UPI(i)=((1-r)*UPI(i)+r*(UTEMP(i-1)+UTEMP(i+1))+2*h*(UTEMP(i)-UPI(i))+h*kt)/(1+r+2*h+h*sig*UPI(i)*UTEMP(i)^2);
    end
   UPI(1)=B1(t+1);UPI(N)=B2(t+1);
end
toc
%%Max errors
MaxD(ih,1+ip)=max(abs(UPI-uex)); 
end
h=h/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Green         = [0 0.9 0.3]; 
dark_grey     = [0.5,0.5,0.5]; 
Magenta	 =    [1, 0, 1];
dark_yellow   = [0.9290, 0.6940, 0.1250]; 
Red	= [1, 0, 0];
	%% figure 1:
  figure('Name', 'Errors as a fuction of time step h');   

	plot(Axhstep(:),MaxD(:,1), "--o", 'Color',Red, 'LineWidth', 2, 'MarkerSize',6); %% UPFD 
	hold on;
    plot(Axhstep(:),MaxD(:,4), ':',  'Color', dark_yellow, 'LineWidth', 7, 'MarkerSize',10); %% Lambda=1/3
 	hold on;
    plot(Axhstep(:),MaxD(:,3), '->',  'Color', Green, 'LineWidth', 1.2, 'MarkerSize',15); %% Lambda=1/2
	hold on;
    plot(Axhstep(:),MaxD(:,2), ':p',  'Color', Magenta, 'LineWidth', 2, 'MarkerSize',14); %% Lambda=1
	hold on;
    plot(Axhstep(:),rfun, '--',  'Color', dark_grey, 'LineWidth', 1, 'MarkerSize',8); %% reference
	hold on;
  	hold off;
	
XLabel = xlabel('Time step size h');
set(XLabel, 'FontSize', 14);
set(XLabel, 'FontWeight', 'bold');

YLabel = ylabel('Errors');
set(YLabel, 'FontSize', 14); 
set(YLabel, 'FontWeight', 'bold');

% Get handle to current axes.
Ax = gca;
Ax.XAxis.FontSize = 16;
Ax.YAxis.FontSize = 16;
Ax.FontWeight = 'bold';
Ax.TickLength = [0.018 0.035];
Ax.LineWidth = 1.2;

Legend=legend({' UPFD orig',' PI, \lambda=1/3',' PI, \lambda=1/2',' PI, \lambda=1',' ~h^2'},'Location','southeast');
set(Legend, 'FontSize', 13); 
set(Legend, 'FontWeight', 'bold');
set(gca,'xScale','log');
set(gca,'yScale','log');
grid on;