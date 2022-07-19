set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultLegendFontName','Times New Roman');

%Banda A

s11_s21a = sparameters('s11_s21_abs44_502cLA.s2p');                          %Measured data
s11a = rfparam(s11_s21a,1,1);
s21a = rfparam(s11_s21a,2,1); 
s12a = rfparam(s11_s21a,1,2); 
s22a = rfparam(s11_s21a,2,2);
phase_s11a = rad2deg(angle(s11a));                                              
phase_s21a = rad2deg(angle(s21a));            
phase_s22a = rad2deg(atan2(imag(s22a),real(s22a)));  
phase_s12a = rad2deg(atan2(imag(s12a),real(s12a)));
freqa = transpose(linspace(5.0,5.50,length(s11a)));
s11a(201) = [];
s21a(201) = []; 
s12a(201) = []; 
s22a(201) = [];
phase_s11a(201) = [];                                              
phase_s21a(201) = [];            
phase_s22a(201) = [];  
phase_s12a(201) = []; 
%Banda B

s11_s21b = sparameters('s11_s21_abs44_502cLB.s2p');                          %Measured data
s11b = rfparam(s11_s21b,1,1);
s21b = rfparam(s11_s21b,2,1); 
s12b = rfparam(s11_s21b,1,2); 
s22b = rfparam(s11_s21b,2,2);
phase_s11b = rad2deg(angle(s11b));                                              
phase_s21b = rad2deg(angle(s21b));            
phase_s22b = rad2deg(atan2(imag(s22b),real(s22b)));  
phase_s12b = rad2deg(atan2(imag(s12b),real(s12b))); 
freqb = transpose(linspace(5.50,6.0,length(s11b)));  
s11b(201) = [];
s21b(201) = []; 
s12b(201) = []; 
s22b(201) = [];
phase_s11b(201) = [];                                              
phase_s21b(201) = [];            
phase_s22b(201) = [];  
phase_s12b(201) = []; 

%Banda C

s11_s21c = sparameters('s11_s21_abs44_502cLC.s2p');                          %Measured data
s11c = rfparam(s11_s21c,1,1);
s21c = rfparam(s11_s21c,2,1); 
s12c = rfparam(s11_s21c,1,2); 
s22c = rfparam(s11_s21c,2,2);
phase_s11c = rad2deg(angle(s11c));                                              
phase_s21c = rad2deg(angle(s21c));            
phase_s22c = rad2deg(atan2(imag(s22c),real(s22c)));  
phase_s12c = rad2deg(atan2(imag(s12c),real(s12c)));
freqc = transpose(linspace(6.0,6.50,length(s11c))); 
s11c(201) = [];
s21c(201) = []; 
s12c(201) = []; 
s22c(201) = [];
phase_s11c(201) = [];                                              
phase_s21c(201) = [];            
phase_s22c(201) = [];  
phase_s12c(201) = []; 

%Banda D

s11_s21d = sparameters('s11_s21_abs44_502cLD.s2p');                          %Measured data
s11d = rfparam(s11_s21d,1,1);
s21d = rfparam(s11_s21d,2,1); 
s12d = rfparam(s11_s21d,1,2); 
s22d = rfparam(s11_s21d,2,2);
phase_s11d = rad2deg(angle(s11d));                                              
phase_s21d = rad2deg(angle(s21d));            
phase_s22d = rad2deg(atan2(imag(s22d),real(s22d)));  
phase_s12d = rad2deg(atan2(imag(s12d),real(s12d)));
freqd = transpose(linspace(6.50,7.0,length(s11d))); 

%Banda total

s11 = [s11a;s11b;s11c;s11d];
s21= [s21a;s21b;s21c;s21d]; 
s12 = [s12a;s12b;s12c;s12d]; 
s22 = [s22a;s22b;s22c;s22d];
phase_s11 = [phase_s11a;phase_s11b;phase_s11c;phase_s11d];                                              
phase_s21 = [phase_s21a;phase_s21b;phase_s21c;phase_s21d];            
phase_s22 = [phase_s22a;phase_s22b;phase_s22c;phase_s22d];  
phase_s12 = [phase_s12a;phase_s12b;phase_s12c;phase_s12d]; 

%generacion de frecuencias    

freqd = transpose(linspace(6.50,7.0,length(s11d)));    %Frequency vector
freqa(201)=[];
freqb(201)=[];
freqc(201)=[];
freq=[freqa;freqb;freqc;freqd];

k = 2.*pi.*freq.*1e9.*sqrt(4*pi*10^-7*8.8541878128e-12*1);              % free space wavenumber
a = 40.39e-3;                                                              % Waveguide dimension 
kc = pi/a;                                                                 % Cutoff wavenumber       
beta = sqrt(k.^2-kc^2);                                                  % Phase constant of the TE10 Mode    
z = 0e-3;                                                                  % Distance of the sample to the reference plane of the VNA (Ideally 0)                                            
d_holder = 7e-3;                                                           % Thickness of the sample

mag_ref_s11 = abs(s11); 
phase_ref_s11 = wrapToPi((deg2rad(phase_s11)-(-2.*z.*beta)));    
phase_ref_s11 = rad2deg(phase_ref_s11); 

mag_ref_s21 = abs(s21); 
phase_ref_s21 = wrapToPi(deg2rad(phase_s21)-(-2.*z.*beta));
phase_ref_s21 = rad2deg(phase_ref_s21); 

%genera graficos
figure
yyaxis left
plot(freq,abs(s11),'LineWidth',2); 
ylabel('Magnitud (-)','FontSize',15)
yyaxis right
plot(freq,phase_ref_s11,'LineWidth',2);
ylabel('Fase (°)','FontSize',15)
xlabel('Frecuencia (GHz)','FontSize',15)
hold on 
xlim([5.0 7.0])
figure
yyaxis left
plot(freq,abs(s21),'LineWidth',2); 
ylabel('Magnitud (-)','FontSize',15)
yyaxis right
plot(freq,phase_ref_s21,'LineWidth',2); 
hold on 
ylabel('Fase (°)','FontSize',15)
xlabel('Frecuencia (GHz)','FontSize',15)
hold on 
xlim([5.0 7.0])
%% 
s11_vec = mag_ref_s11.*cosd(phase_ref_s11)+j.*mag_ref_s11.*sind(phase_ref_s11);  
s21_vec =  mag_ref_s21.*cosd(phase_ref_s21)+j.*mag_ref_s21.*sind(phase_ref_s21); 

plot(freq,real(s11_vec),'LineWidth',2);
hold on 
plot(freq,imag(s11_vec),'LineWidth',2); 
hold on 
plot(freq,real(s21_vec),'LineWidth',2); 
hold on
plot(freq,imag(s21_vec),'LineWidth',2); 
xlim([4.95 7.0])
legend('Re S11','Imag S11','Re S21','Imag S21'); 


% constantes
c = 3e8; 
fc = 3e8/(2*a); 
lambda_c = c/fc; 
lambda_0 = c./(freq.*1e9); 
lambda_g = lambda_0./(sqrt(1-(lambda_0./lambda_c).^2));                    % guided wavelength 
d = 7e-3;
n = 0;                                                                     % integrador para muestras eléctricamente pequeñas 
s11 = s11_vec;                                                             % s11 vector complejo 
s21 = s21_vec;                                                             % s21 vector complejo 

X = (s11.^2 - s21.^2 +1)./(2.*s11); 

v1 = s21 + s11; 

for i = 1:length(X)
  Gamma(i,1) = X(i,1) + sqrt(X(i,1)^2 - 1);
  if abs(Gamma(i,1)) > 1
    Gamma(i,1) = X(i,1) - sqrt(X(i,1)^2 - 1);
  end
end

T_perp = (s11 + s21 - Gamma)./(1 - (s11 + s21).*Gamma);

tau = 1i./(2*pi*d)*log(T_perp);
mu_r = ones(length(s11),1);                                                 %lambda_g.*tau_perp.*(1+Gamma_perp)./(1-Gamma_perp);  %mu_r = ones(600,1); 
eps_r = lambda_0.^2.*(tau.^2 + (1/lambda_c.^2)) .* 1./mu_r; 
eps_r_p = polyfit(freq,eps_r,1); 
eps_r_pol = polyval(eps_r_p,freq); 

figure
plot(freq,real(eps_r),'->','LineWidth',2,'MarkerIndices',1:50:length(eps_r)); 
figure
plot(freq,abs(imag(eps_r)),'--<','LineWidth',2,'MarkerIndices',1:50:length(eps_r));

ylabel('\epsilon_r','FontSize',15)
xlabel('Frequency (GHz)','FontSize',15)
ylim([-1 1])
grid on
xlim([5.0 7.0])
legend('Parte imaginaria')


tan_delta = (2.*pi.*freq.*1e9.*abs(imag(eps_r)))./(2.*pi.*freq.*1e9.*real(eps_r));

mean_a = mean(tan_delta);
mean_tan_delta = mean_a.*ones(1,length(tan_delta)).'; 

figure
plot(freq,tan_delta,'LineWidth',2)
hold on 
plot(freq,mean_tan_delta,'--','LineWidth',2)
hold on
%yline(0.033,'--','LineWidth',2); 
hold on 
%yline(0.017,'--s','LineWidth',2); 
ylabel('tan\delta','FontSize',15); 
xlabel('Frequency (GHz)','FontSize',15); 
xlim([5.0 7.0])
%%
plot(freq,abs(imag(eps_r)),'--<','LineWidth',2,'MarkerIndices',1:50:length(eps_r));
hold on
plot(freqan,abs(imag(eps_an)),'--<','LineWidth',2,'MarkerIndices',1:50:length(eps_an));
xlim([5.0 7.0])
ylim([-2 2])
ylabel('\epsilon_r','FontSize',15)
xlabel('Frequency (GHz)','FontSize',15)
grid on
legend('Medicion Banda Estrecha','Medicion Banda Ancha')
title('\epsilon_r imaginaria Banda estrecha vs Banda Ancha')