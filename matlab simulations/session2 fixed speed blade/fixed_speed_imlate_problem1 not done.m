clear;clc;
vw=[5 8 11 14]; %wind speed values to study (from the exercise)
ss=0;           %neglecting slip=0
rho=1.225;
N=80;
D=100;
poles=2;        %pair of poles
Ugrid=960;
f=50;

we=2*pi*f;      %electrical speed
ws=we/poles;    %synchronous speed
wg=ws*(1-ss);   %generator speed fast axis
wt=wg/N;        %turbine speed slow axis
lam= wt*(D/2) ./ vw;    %tip speed ratio
Cpp= 0.0045 * (100 - (lam-10).^2);
P1=0.5*rho*Cpp*(pi*(D^2)/4).*vw.^3;

vw = [5 8 11 14];
lam = [19.634954084936204 12.271846303085129 8.924979129516457 7.012483601762931];
Cpp = [0.032254469015270 0.426774214688213 0.444799485576112 0.409836355966191 ];
P1 = 1.0e+006 * [0.019395272430808 1.051148578283138 2.847988990930852 5.409916510373950];

%% Develop a Matlab program to locate the operating points in the Cp-lambda curve.
subplot(2,2,1);
lam1=0:.1:20;
Cpp1= 0.0045 * (100 - (lam1-10).^2);
plot(lam1,Cpp1);hold on;grid on;
plot(lam,Cpp,'ko');
for ii=1:1:4
 txt{ii}=['v_w = ' num2str(vw(ii)) ' m/s'];
 text(lam(ii),Cpp(ii),txt{ii},'FontSize',18);
end
xlabel('\lambda','FontSize',18);
ylabel('C_p','FontSize',18);


%% Develop a Matlab program to plot the power generated for different wind
% speeds. Include a comparison between the obtained power and the maximum available power.
Pn = 3.2e6;%nominal power limit of theturbine in W
h=subplot(2,2,2);
vw2=0:.1:15;
lam2= min(20,max(wt*(D/2) ./ vw2,0));
Cpp2= 0.0045 * (100 - (lam2-10).^2);
P2=0.5*rho*Cpp2*(pi*(D^2)/4).*vw2.^3;
P2a=min(0.5*rho*Cpp2*(pi*(D^2)/4).*vw2.^3,Pn);
plot(vw2,P2,':k');hold on;grid on;
plot(vw2,P2a,'k');
plot(vw,P1,'ko');
for ii=1:1:4
 txt{ii}=['v_w = ' num2str(vw(ii)) ' m/s'];
 text(vw(ii),P1(ii),txt{ii},'FontSize',18);
end
xlabel('v_w','FontSize',18);
ylabel('P','FontSize',18);
set(h,'FontSize',18);


h=subplot(2,2,4);
lam2b= 10;      %impose
Cpp2b= 0.45;    %maximum cpp (maximum power extraction), means perfect variable speed has been archivied
P3=0.5*rho*Cpp2b*(pi*(D^2)/4).*vw2.^3;
P3_limited=min(0.5*rho*Cpp2b*(pi*(D^2)/4).*vw2.^3,Pn);
plot(vw2,P2,':k');hold on;grid on;
plot(vw2,P2a,'k','LineWidth',3);
plot(vw,P1,'ko');
plot(vw2,P3,':r');
plot(vw2,P3_limited,'r');
for ii=1:1:4
 txt{ii}=['v_w = ' num2str(vw(ii)) ' m/s'];
 text(vw(ii),P1(ii),txt{ii},'FontSize',18);
end
xlabel('v_w','FontSize',18);
ylabel('P','FontSize',18);
set(h,'FontSize',18);
axis([3 15 0 6e6]);


%% How should we modify the scripts in order to consider appropriately the maximum power limitation (with implied pitch control to limit the power)?

%use the limited power curve and calculate the Cpp back from P

subplot(2,2,3);
Cpp_limited=P3_limited./(0.5*rho*(pi*(D^2)/4).*vw2.^3);
plot(lam1,Cpp1);hold on;grid on;
plot(lam2,Cpp_limited,':r');
for ii=1:1:4
 txt{ii}=['v_w = ' num2str(vw(ii)) ' m/s'];
 text(lam(ii),Cpp_limited(ii),txt{ii},'FontSize',18);
end
xlabel('\lambda','FontSize',18);
ylabel('C_p','FontSize',18);