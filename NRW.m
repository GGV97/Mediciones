set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultTextFontName','Times New Roman');
set(0,'defaultLegendFontName','Times New Roman');

s11_s21 = sparameters('s11_s21_abs10_100.s2p');                          %Measured data
s11 = rfparam(s11_s21,1,1);
s21 = rfparam(s11_s21,2,1); 
s12 = rfparam(s11_s21,1,2); 
s22 = rfparam(s11_s21,2,2);
phase_s11 = rad2deg(angle(s11));                                              
phase_s21 = rad2deg(angle(s21));            
phase_s22 = rad2deg(atan2(imag(s22),real(s22)));  
phase_s12 = rad2deg(atan2(imag(s12),real(s12)));  

freq = transpose(linspace(4.9,7.05,length(s11)));                           %Frequency vector


k = 2.*pi.*freq.*1e9.*sqrt(4*pi*10^-7*8.8541878128e-12*1);                 % free space wavenumber
a = 40.39e-3;                                                              % Waveguide dimension 
kc = pi/a;                                                                 % Cutoff wavenumber       
beta = sqrt(k.^2-kc^2);                                                    % Phase constant of the TE10 Mode    
z = 0e-3;                                                                  % Distance of the sample to the reference plane of the VNA (Ideally 0)                                            
d_holder = 7e-3;                                                           % Thickness of the sample

mag_ref_s11 = abs(s11); 
phase_ref_s11 = wrapToPi((deg2rad(phase_s11)-(-2.*z.*beta)));    
phase_ref_s11 = rad2deg(phase_ref_s11); 

mag_ref_s21 = abs(s21); 
phase_ref_s21 = wrapToPi(deg2rad(phase_s21)-(-2.*z.*beta));
phase_ref_s21 = rad2deg(phase_ref_s21); 

figure
yyaxis left
plot(freq,abs(s11),'LineWidth',2); 
ylabel('Magnitud (-)','FontSize',15)
yyaxis right
plot(freq,phase_ref_s11,'LineWidth',2);
ylabel('Fase (°)','FontSize',15)
xlabel('Frecuencia (GHz)','FontSize',15)
hold on 

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
hold on
plot(freq,abs(imag(eps_r)),'--<','LineWidth',2,'MarkerIndices',1:50:length(eps_r));

ylabel('\epsilon_r','FontSize',15)
xlabel('Frequency (GHz)','FontSize',15)
ylim([0 10])
grid on
xlim([4.9 7.05])



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
xlim([4.9 7.05])

%%
%epse=eps_r
%epse(max(size(epse)):max(size(epse))+200)=eps_r
figure
freq = transpose(linspace(5.0,7.0,length(epse)));   
plot(freq,real(epse),'->','LineWidth',2,'MarkerIndices',1:50:length(epse)); 
hold on
plot(freq,abs(imag(epse)),'--<','LineWidth',2,'MarkerIndices',1:50:length(epse));

ylabel('\epsilon_r','FontSize',15)
xlabel('Frequency (GHz)','FontSize',15)
ylim([0 15])
grid on
xlim([5.0 7.0])
