% Spatial Modulation in MIMO Systems

clc 
close all
clear all

%% Defining constants

pe_max=10^-5;                                % max ser allowed
N0=3.981*10^(-21);                           % Noise PSD in W/Hz. = -174 in dBm/Hz
B=10^6;                                      % bandwidth = 1 MHz.
sig_n2=B*N0;                                 % sigma n square in watt 
Ga=1;                                        % Antenna Gain
eta=0.1;                                     % efficiency of amplifier
Pc=10^(-8);                                  % Power consumption in channel coding
Pb=4.09*10^(-9);                             % Power consumption in baseband processing
Pf=0.1;                                      % fixed circuit power consumption
Pc1=0.04;                                    % circuit power consumption proportional to system parameters
Pc2=1.3*10^(-8);                             % circuit power consumption proportional to system parameters
d_loss=1.58489*10^10;                        % large-scale channel fading loss = 102dB
Nt=8;                                        % Number of transmitters
Nr=4;                                        % Number of receivers
%Ptr_max=0.251;                               % Maximum allowed Txd Power = 24dBm
Rmin=8*10^6;                                 % Minimum required data rate = 8 Mbps


%% Generate Rayleigh fading channel coefficients

H = 1/sqrt(2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));    % Channel matrix
Rr = eye(Nr,Nr);                                 % Receiver correlation matrix
Rt = eye(Nt,Nt);                                 % Transmit correlation matrix
P=rank(Rr);
k=0.5;                                         % fraction of SE
%% 

Nm=[2;4;8;16;32;64;128;256];                     
Ptm=0.01*[1;2;3;4;5;6;7;8];
chi_cm=(10^4)*[3.5; 14; 58; 220; 930; 3700; 15000;  60000];    % Assuming Nam=2 initially
%P_minm=[0.00316 0.00501187 0.00758577 0.0107152 0.01412537 0.02041738 0.0302 0.05012].';
% for i=1:8
% rho_minm(i,1)=(er_max.*Nm(i).*Nm(i).*16.*1.*factorial(P).*factorial(P-1))/(factorial(2*P-1).*chi_cm(i));
% end

Km = [1;1;1.6;1.6;1.57;2.14 ;2.4;2.9];     %PAPR

Nam=[1;2;3;4];

temp_modes=[Nm Ptm chi_cm Km];

for i=1:size(temp_modes,1)
    modes((i-1)*4+1,:)=[temp_modes(i,:) Nam(1)];
    modes((i-1)*4+2,:)=[temp_modes(i,:) Nam(2)];
    modes((i-1)*4+3,:)=[temp_modes(i,:) Nam(3)];
    modes((i-1)*4+4,:)=[temp_modes(i,:) Nam(4)];
end

Nm=modes(:,1);
Ptm=modes(:,2);
chi_cm=modes(:,3);
Km=modes(:,4);
Nam=modes(:,5);

Ptrm=Ptm.*Nam;
rho_m=(Ga*Ptrm)/(sig_n2 * d_loss);

for i=1:size(Nam)
    Nlm(i)=2^(floor(log2(nchoosek(Nt,Nam(i)))));
end
Nlm=Nlm.';
Ncm=(Nm.^Nam).*Nlm;

pe_m=(nchoosek(2*P-1,P))*(((rho_m.^(-P)).*chi_cm)./Ncm);



% %Eme=mean(Em);
% %Em=4*Em-3*Eme*ones(size(modes,1),1);
% %Em=Em-2*min(Em)*ones(size(modes,1),1);

% 
Ptr_max=0.045:0.01:0.33;Ptr_max=Ptr_max';
% 
for j=1:size(Ptr_max,1)

%Rmin=2*10^6:10^6:30*10^6; Rmin=Rmin.';
% 
%for j=1:size(Rmin,1)
    
Sm=log2(Nlm)+Nam.*log2(Nm);
Rm=B*Sm;

Ptot=Pc*Rm+ Nam*B*Pb + Pf*ones(size(modes,1),1) + Nam*Pc1 + Nam*B*Pc2 + (1/eta)*(Ptrm.*Km);
Em=Rm./Ptot;
    
[Smax]=max(Sm);
[Emax]=max(Em);

Sm_norm=Sm/Smax;
Em_norm=Em/Emax;
    
for i=1:size(modes,1)
    if Ptrm(i)>Ptr_max(j)
        Em(i)=0;
        Sm(i)=0;
        Sm_norm(i)=0;
        Em_norm(i)=0;
        
    elseif pe_m(i)>pe_max
        Em(i)=0;
        Sm(i)=0;
        Sm_norm(i)=0;
        Em_norm(i)=0;
        
    elseif Rm(i)<Rmin
        Em(i)=0;
        Sm(i)=0;
        Sm_norm(i)=0;
        Em_norm(i)=0;
    end
end


[Sm_max, m_sm]=max(Sm);
[Em_max, m_em]=max(Em);

Wm=k*Sm_norm+(1-k)*Em_norm;

[W_max(j), m_opt(j)]  = max(Wm);

Sm_opt(j)=Sm_norm(m_opt(j));
Em_opt(j)=Em_norm(m_opt(j));

Wm_in_Sm(j)=k*Sm_norm(m_sm)+(1-k)*Em_norm(m_sm);
       
Wm_in_Em(j)=k*Sm_norm(m_em)+(1-k)*Em_norm(m_em);

end


 figure(1)
plot(1000*Ptr_max,W_max,'b-','LineWidth',3);
hold on
plot(1000*Ptr_max,Wm_in_Sm,'g--','LineWidth',3);
plot(1000*Ptr_max,Wm_in_Em,'r-.','LineWidth',3);
axis([30 350 0.48 0.66])
% % % title('Convergence Plot','FontWeight','bold')
xlabel('Maximum transmission power (mW)','FontWeight','bold')
ylabel('WEF','FontWeight','bold')
legend('WEF-based','SE-based','EE-based');


