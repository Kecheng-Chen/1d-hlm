%%
%2D axissymmetric and 2D Vinsome
clc;clear;
Q=0.0075.*2;
H=20;
time_array=[1,2,3,4,5,10,20];
poro=0.25;
rho_s=2600;
rho_w=1000;
c_s=1040;
c_w=4000;
rho_c=(1-poro)*rho_s*c_s + poro*rho_w*c_w;
lambda=2.2;
figure;hold on;
for i=1:7
    
    time=time_array(i);
    t=time*365*24*60*60;
    r_0=sqrt(Q*rho_w*c_w*t/(rho_c*pi*H));
    
    A=dlmread('2D axissymmetric output.txt', '', 9, 0);
    A_aq=A(find(A(:,2)<=H/2),:);
    x0 = A_aq(:,1);
    y = A_aq(:,2);
    rangex=sort(uniquetol(x0,0.00001));
    rangey=sort(uniquetol(y,0.01));
    [X,Y] = meshgrid(rangex,rangey);
    Z = griddata(A(:,1),A(:,2),A(:,i+2),X,Y);
    avg_T=zeros(length(rangex),1);
    for dp=1:length(rangex)
        B=Z(:,dp);
        avg_T(dp)=trapz(rangey,B)./(H/2);
    end
    plot(rangex./r_0,(avg_T-14.8)./(25.8-14.8))
    
    A2=dlmread('2D Vinsome output.txt', '', 9, 0);
    A2_aq=A2(find(A2(:,2)<=H/2),:);
    x02 = A2_aq(:,1);
    y2 = A2_aq(:,2);
    rangex2=sort(uniquetol(x02,0.00001));
    rangey2=sort(uniquetol(y2,0.01));
    [X2,Y2] = meshgrid(rangex2,rangey2);
    Z2 = griddata(A2(:,1),A2(:,2),A2(:,i+2),X2,Y2);
    avg_T2=zeros(length(rangex2),1);
    for dp2=1:length(rangex2)
        B2=Z2(:,dp2);
        avg_T2(dp2)=trapz(rangey2,B2)./(H/2);
    end
    plot(rangex2./r_0,(avg_T2-14.8)./(25.8-14.8),'--');
    
end
xlim([0,1.5])
ylim([0,1])
hold off;
xlabel('r^*');
ylabel('T^*');


%%
%2D axissymmetric and 1D Vinsome
clc;clear;
Q=0.0075.*2;
H=20;
time_array=[1,2,3,4,5,10,20];
poro=0.25;
rho_s=2600;
rho_w=1000;
c_s=1040;
c_w=4000;
rho_c=(1-poro)*rho_s*c_s + poro*rho_w*c_w;
lambda=2.2;
figure;hold on;
for i=1:7
    
    time=time_array(i);
    t=time*365*24*60*60;
    r_0=sqrt(Q*rho_w*c_w*t/(rho_c*pi*H));
    
    A=dlmread('2D axissymmetric output.txt', '', 9, 0);
    A_aq=A(find(A(:,2)<=H/2),:);
    x0 = A_aq(:,1);
    y = A_aq(:,2);
    rangex=sort(uniquetol(x0,0.00001));
    rangey=sort(uniquetol(y,0.01));
    [X,Y] = meshgrid(rangex,rangey);
    Z = griddata(A(:,1),A(:,2),A(:,i+2),X,Y);
    avg_T=zeros(length(rangex),1);
    for dp=1:length(rangex)
        B=Z(:,dp);
        avg_T(dp)=trapz(rangey,B)./(H/2);
    end
    plot(rangex./r_0,(avg_T-14.8)./(25.8-14.8))
    
    r_star_array=xlsread(strcat('1D Vinsome output r',...
    '.xlsx'));
    T_star_array=xlsread(strcat('1D Vinsome output T',...
    '.xlsx'));
    plot(r_star_array(2:end,i),T_star_array(2:end,i),'--');
    
end
xlim([0,1.5])
ylim([0,1])
hold off;
xlabel('r^*');
ylabel('T^*');


%%
%2D axissymmetric and 1D Newton
clc;clear;
syms h_star x r_star0;
Q=0.0075.*2;
H=20;
time_array=[1,2,3,4,5,10,20];
poro=0.25;
rho_s=2600;
rho_w=1000;
c_s=1040;
c_w=4000;
rho_c=(1-poro)*rho_s*c_s + poro*rho_w*c_w;
lambda=2.2;
figure;hold on;
for i=1:7
    
    time=time_array(i);
    t=time*365*24*60*60;
    r_0=sqrt(Q*rho_w*c_w*t/(rho_c*pi*H));
    
    A=dlmread('2D axissymmetric output.txt', '', 9, 0);
    A_aq=A(find(A(:,2)<=H/2),:);
    x0 = A_aq(:,1);
    y = A_aq(:,2);
    rangex=sort(uniquetol(x0,0.00001));
    rangey=sort(uniquetol(y,0.01));
    [X,Y] = meshgrid(rangex,rangey);
    Z = griddata(A(:,1),A(:,2),A(:,i+2),X,Y);
    avg_T=zeros(length(rangex),1);
    for dp=1:length(rangex)
        B=Z(:,dp);
        avg_T(dp)=trapz(rangey,B)./(H/2);
    end
    plot(rangex./r_0,(avg_T-14.8)./(25.8-14.8))
    
    rhoc_f=rho_w*c_w;
    D=lambda/rho_c;
    t_0=H^2/(4*D);
    alpha=Q*rhoc_f/(4*pi*H*lambda);  
    fun= @(x,h_star,r_star0) exp(-alpha*h_star*t/t_0*r_star0^2/x-x-gammaln(alpha)+(alpha-1)*log(x));
    T = vpaintegral(fun,x,[alpha*r_star0^2,Inf], 'RelTol', 1e-12, 'AbsTol', 0);
    T2 = @(h0,r0) double(subs(subs(T,r_star0,r0),h_star,h0));
    [h_best,resnorm] = lsqcurvefit(T2,2,rangex./r_0,(avg_T-14.8)./(25.8-14.8));
    plot(rangex./r_0,T2(h_best,rangex./r_0),'--');
    
end
xlim([0,1.5])
ylim([0,1])
hold off;
xlabel('r^*');
ylabel('T^*');