% Spatial Modulation in MIMO Systems

clc 
close all
clear all

%% Defining constants

er_max=10^-6;                                % max ser allowed
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
Ptr_max=0.251;                                  % Maximum allowed Txd Power = 24dBm


%% Generate Rayleigh fading channel coefficients

H = 1/sqrt(2)*(randn(Nr,Nt)+1i*randn(Nr,Nt));    % Channel matrix
Rr = eye(Nr,Nr);                                 % Receiver correlation matrix
Rt = eye(Nt,Nt);                                 % Transmit correlation matrix
P=rank(Rr);
k=0.5;                                         % fraction of SE
%% 

Nm=[2;4;8;16;32;64;128;256];                     
Ptm=0.01*[1;2;3;4;5;6;7;8];
chi_cm=10^(-8)*[3.5; 14; 58; 220; 930; 3700; 15000;  60000];    % Assuming Nam=2 initially
P_minm=[0.00316 0.00501187 0.00758577 0.0107152 0.01412537 0.02041738 0.0302 0.05012].';
% for i=1:8
% rho_minm(i,1)=(er_max.*Nm(i).*Nm(i).*16.*1.*factorial(P).*factorial(P-1))/(factorial(2*P-1).*chi_cm(i));
% end

Km = [1;1;1.6;1.6;1.57;2.14 ;2.4;2.9];     %PAPR
modes=[Nm Ptm chi_cm P_minm Km];

for i=1:8
    j=0;
    for Ptrm=P_minm(1):0.0001:Ptr_max
        j=j+1;
        if(Ptrm < P_minm(i))
            continue;
        end
        Na(i,j)=floor(Ptrm/Ptm(i));
        if(Na(i,j)<1) 
            continue;
        end
        if(Na(i,j)>Nr)
            Na(i,j)=Nr;
        end
        
        Power(i,j)=Ptrm;
        
        SE(i,j)=floor(log2(nchoosek(Nt,Na(i,j))))+Na(i,j)*log2(Nm(i));
        
        Ptotal(i,j) = Pc*B*SE(i,j)+ Na(i,j)*(B*Pb+Pc1+B*Pc2+((Ptm(i)*Km(i))/eta)) + Pf;

        EE(i,j) = SE(i,j)/Ptotal(i,j);                                      % Energy Efficiency
        
    end
    
    EE_avg(i) = mean(nonzeros(EE(i,:)));
end

for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
        if(EE(i,j)<EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        elseif(EE(i,j)>EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        end
    end
    SE_max(i)=max(SE(i,:));

    EE_max(i)=max(EE(i,:));

end
    SE_max_global=max(SE_max);
    EE_max_global=max(EE_max);


    figure(1)
    title('Variation of SE with Power')
    subplot(4,2,1)
    plot(nonzeros(Power(1,:)*1000),nonzeros(SE(1,:)),'LineWidth',2)
    title('Mode 1 : BPSK','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')
    
    subplot(4,2,3)
    plot(nonzeros(Power(2,:)*1000),nonzeros(SE(2,:)),'LineWidth',2)
    title('Mode 2 : QPSK','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')

    subplot(4,2,5)
    plot(nonzeros(Power(3,:)*1000),nonzeros(SE(3,:)),'LineWidth',2)
    title('Mode 3 : 8-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')
    
    subplot(4,2,7)
    plot(nonzeros(Power(4,:)*1000),nonzeros(SE(4,:)),'LineWidth',2)
    title('Mode 4 : 16-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')

    subplot(4,2,2)
    plot(nonzeros(Power(5,:)*1000),nonzeros(SE(5,:)),'LineWidth',2)
    title('Mode 5 : 32-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')
    
    subplot(4,2,4)
    plot(nonzeros(Power(6,:)*1000),nonzeros(SE(6,:)),'LineWidth',2)
    title('Mode 6 : 64-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')
    
    subplot(4,2,6)
    plot(nonzeros(Power(7,:)*1000),nonzeros(SE(7,:)),'LineWidth',2)
    title('Mode 7 : 128-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')
    
    subplot(4,2,8)
    plot(nonzeros(Power(8,:)*1000),nonzeros(SE(8,:)),'LineWidth',2)
    title('Mode 8 : 256-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('SE (bits/s/Hz)')
    hold on;

    
    figure(2)
    title('Variation of EE with Power')
    subplot(4,2,1)
    plot(nonzeros(Power(1,:)*1000),nonzeros(EE(1,:)),'LineWidth',2)
    title('Mode 1 : BPSK','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')
    
    subplot(4,2,3)
    plot(nonzeros(Power(2,:)*1000),nonzeros(EE(2,:)),'LineWidth',2)
    title('Mode 2 : QPSK','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')

    subplot(4,2,5)
    plot(nonzeros(Power(3,:)*1000),nonzeros(EE(3,:)),'LineWidth',2)
    title('Mode 3 : 8-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')
    
    subplot(4,2,7)
    plot(nonzeros(Power(4,:)*1000),nonzeros(EE(4,:)),'LineWidth',2)
    title('Mode 4 : 16-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')

    subplot(4,2,2)
    plot(nonzeros(Power(5,:)*1000),nonzeros(EE(5,:)),'LineWidth',2)
    title('Mode 5 : 32-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')
    
    subplot(4,2,4)
    plot(nonzeros(Power(6,:)*1000),nonzeros(EE(6,:)),'LineWidth',2)
    title('Mode 6 : 64-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')
    
    subplot(4,2,6)
    plot(nonzeros(Power(7,:)*1000),nonzeros(EE(7,:)),'LineWidth',2)
    title('Mode 7 : 128-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')
    
    subplot(4,2,8)
    plot(nonzeros(Power(8,:)*1000),nonzeros(EE(8,:)),'LineWidth',2)
    title('Mode 8 : 256-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('EE (bits/\muJ)')
    hold on;
    
for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
         SE(i,j)=SE(i,j)/SE_max(i);
         EE(i,j)=EE(i,j)/EE_max(i);
        
        f(i,j)= k*SE(i,j)+(1-k)*EE(i,j);
    end
end

figure(3)
title('Variation of WEF with Power')
    subplot(4,2,1)
    plot(nonzeros(Power(1,:)*1000),nonzeros(f(1,:)),'LineWidth',2)
    title('Mode 1 : BPSK','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')
    
    subplot(4,2,3)
    plot(nonzeros(Power(2,:)*1000),nonzeros(f(2,:)),'LineWidth',2)
    title('Mode 2 : QPSK','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')

    subplot(4,2,5)
    plot(nonzeros(Power(3,:)*1000),nonzeros(f(3,:)),'LineWidth',2)
    title('Mode 3 : 8-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')
    
    subplot(4,2,7)
    plot(nonzeros(Power(4,:)*1000),nonzeros(f(4,:)),'LineWidth',2)
    title('Mode 4 : 16-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')

    subplot(4,2,2)
    plot(nonzeros(Power(5,:)*1000),nonzeros(f(5,:)),'LineWidth',2)
    title('Mode 5 : 32-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')
    
    subplot(4,2,4)
    plot(nonzeros(Power(6,:)*1000),nonzeros(f(6,:)),'LineWidth',2)
    title('Mode 6 : 64-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')
    
    subplot(4,2,6)
    plot(nonzeros(Power(7,:)*1000),nonzeros(f(7,:)),'LineWidth',2)
    title('Mode 7 : 128-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')
    
    subplot(4,2,8)
    plot(nonzeros(Power(8,:)*1000),nonzeros(f(8,:)),'LineWidth',2)
    title('Mode 8 : 256-QAM','FontWeight','Bold')
    xlabel('Power (mWatts)')
    ylabel('F^{SE,EE}')
    hold on;

    k=0.4;
    for i=1:8
    j=0;
    for Ptrm=P_minm(1):0.0001:Ptr_max
        j=j+1;
        if(Ptrm < P_minm(i))
            continue;
        end
        Na(i,j)=floor(Ptrm/Ptm(i));
        if(Na(i,j)<1) 
            continue;
        end
        if(Na(i,j)>Nr)
            Na(i,j)=Nr;
        end
        
        Power(i,j)=Ptrm;
        
        SE(i,j)=floor(log2(nchoosek(Nt,Na(i,j))))+Na(i,j)*log2(Nm(i));
        
        Ptotal(i,j) = Pc*B*SE(i,j)+ Na(i,j)*(B*Pb+Pc1+B*Pc2+((Ptm(i)*Km(i))/eta)) + Pf;

        EE(i,j) = SE(i,j)/Ptotal(i,j);                                      % Energy Efficiency
        
    end
    
    EE_avg(i) = mean(nonzeros(EE(i,:)));
end

for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
        if(EE(i,j)<EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        elseif(EE(i,j)>EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        end
    end
    SE_max(i)=max(SE(i,:));

    EE_max(i)=max(EE(i,:));

end
    SE_max_global=max(SE_max);
    EE_max_global=max(EE_max);
    
    
    for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
         SE(i,j)=SE(i,j)/SE_max_global;
         EE(i,j)=EE(i,j)/EE_max_global;
        
        f(i,j)= k*SE(i,j)+(1-k)*EE(i,j);
    end
end
for j=1:2479
    [max_SE(j) I_SE(j)]=max(SE(:,j));
   [max_EE(j) I_EE(j)]=max(EE(:,j));
   [max_f(j) I_f(j)]=max(f(:,j));
end
   
figure(4)
subplot(4,2,1)
plot(nonzeros(Power(1,:))*1000,nonzeros(max_f(:)),'LineWidth',2)
axis([0 250 0.3 0.72])
title('WEF-based Optimal Mode Selection: k=0.35','FontWeight','Bold')
xlabel('Power (mWatts)')
ylabel('F^{SE,EE}_{max}');

subplot(4,2,3)
plot(nonzeros(Power(1,:))*1000,nonzeros(I_f(70:2479)),'LineWidth',2)
axis([0 250 0 8])
xlabel('Power (mWatts)')
ylabel('Mode index');


k=0.5;
    for i=1:8
    j=0;
    for Ptrm=P_minm(1):0.0001:Ptr_max
        j=j+1;
        if(Ptrm < P_minm(i))
            continue;
        end
        Na(i,j)=floor(Ptrm/Ptm(i));
        if(Na(i,j)<1) 
            continue;
        end
        if(Na(i,j)>Nr)
            Na(i,j)=Nr;
        end
        
        Power(i,j)=Ptrm;
        
        SE(i,j)=floor(log2(nchoosek(Nt,Na(i,j))))+Na(i,j)*log2(Nm(i));
        
        Ptotal(i,j) = Pc*B*SE(i,j)+ Na(i,j)*(B*Pb+Pc1+B*Pc2+((Ptm(i)*Km(i))/eta)) + Pf;

        EE(i,j) = SE(i,j)/Ptotal(i,j);                                      % Energy Efficiency
        
    end
    
    EE_avg(i) = mean(nonzeros(EE(i,:)));
end

for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
        if(EE(i,j)<EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        elseif(EE(i,j)>EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        end
    end
    SE_max(i)=max(SE(i,:));

    EE_max(i)=max(EE(i,:));

end
    SE_max_global=max(SE_max);
    EE_max_global=max(EE_max);
    
    
    for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
         SE(i,j)=SE(i,j)/SE_max_global;
         EE(i,j)=EE(i,j)/EE_max_global;
        
        f(i,j)= k*SE(i,j)+(1-k)*EE(i,j);
    end
end
for j=1:2479
    [max_SE(j) I_SE(j)]=max(SE(:,j));
   [max_EE(j) I_EE(j)]=max(EE(:,j));
   [max_f(j) I_f(j)]=max(f(:,j));
end
   

subplot(4,2,5)
plot(nonzeros(Power(1,:))*1000,nonzeros(max_f(:)),'LineWidth',2)
axis([0 250 0.3 0.72])
title('WEF-based Optimal Mode Selection: k=0.45','FontWeight','Bold')
xlabel('Power (mWatts)')
ylabel('F^{SE,EE}_{max}');

subplot(4,2,7)
plot(nonzeros(Power(1,:))*1000,nonzeros(I_f(70:2479)),'LineWidth',2)
axis([0 250 0 8])
xlabel('Power (mWatts)')
ylabel('Mode index');
    
    

k=0.6;
    for i=1:8
    j=0;
    for Ptrm=P_minm(1):0.0001:Ptr_max
        j=j+1;
        if(Ptrm < P_minm(i))
            continue;
        end
        Na(i,j)=floor(Ptrm/Ptm(i));
        if(Na(i,j)<1) 
            continue;
        end
        if(Na(i,j)>Nr)
            Na(i,j)=Nr;
        end
        
        Power(i,j)=Ptrm;
        
        SE(i,j)=floor(log2(nchoosek(Nt,Na(i,j))))+Na(i,j)*log2(Nm(i));
        
        Ptotal(i,j) = Pc*B*SE(i,j)+ Na(i,j)*(B*Pb+Pc1+B*Pc2+((Ptm(i)*Km(i))/eta)) + Pf;

        EE(i,j) = SE(i,j)/Ptotal(i,j);                                      % Energy Efficiency
        
    end
    
    EE_avg(i) = mean(nonzeros(EE(i,:)));
end

for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
        if(EE(i,j)<EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        elseif(EE(i,j)>EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        end
    end
    SE_max(i)=max(SE(i,:));

    EE_max(i)=max(EE(i,:));

end
    SE_max_global=max(SE_max);
    EE_max_global=max(EE_max);
    
    
    for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
         SE(i,j)=SE(i,j)/SE_max_global;
         EE(i,j)=EE(i,j)/EE_max_global;
        
        f(i,j)= k*SE(i,j)+(1-k)*EE(i,j);
    end
end
for j=1:2479
    [max_SE(j) I_SE(j)]=max(SE(:,j));
   [max_EE(j) I_EE(j)]=max(EE(:,j));
   [max_f(j) I_f(j)]=max(f(:,j));
end
   

subplot(4,2,2)
plot(nonzeros(Power(1,:))*1000,nonzeros(max_f(:)),'LineWidth',2)
axis([0 250 0.3 0.72])
title('WEF-based Optimal Mode Selection: k=0.55','FontWeight','Bold')
xlabel('Power (mWatts)')
ylabel('F^{SE,EE}_{max}');

subplot(4,2,4)
plot(nonzeros(Power(1,:))*1000,nonzeros(I_f(70:2479)),'LineWidth',2)
axis([0 250 0 8])
xlabel('Power (mWatts)')
ylabel('Mode index');


k=0.7;
    for i=1:8
    j=0;
    for Ptrm=P_minm(1):0.0001:Ptr_max
        j=j+1;
        if(Ptrm < P_minm(i))
            continue;
        end
        Na(i,j)=floor(Ptrm/Ptm(i));
        if(Na(i,j)<1) 
            continue;
        end
        if(Na(i,j)>Nr)
            Na(i,j)=Nr;
        end
        
        Power(i,j)=Ptrm;
        
        SE(i,j)=floor(log2(nchoosek(Nt,Na(i,j))))+Na(i,j)*log2(Nm(i));
        
        Ptotal(i,j) = Pc*B*SE(i,j)+ Na(i,j)*(B*Pb+Pc1+B*Pc2+((Ptm(i)*Km(i))/eta)) + Pf;

        EE(i,j) = SE(i,j)/Ptotal(i,j);                                      % Energy Efficiency
        
    end
    
    EE_avg(i) = mean(nonzeros(EE(i,:)));
end

for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
        if(EE(i,j)<EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        elseif(EE(i,j)>EE_avg(i))
            EE(i,j)=10*EE(i,j)-9*EE_avg(i);
        end
    end
    SE_max(i)=max(SE(i,:));

    EE_max(i)=max(EE(i,:));

end
    SE_max_global=max(SE_max);
    EE_max_global=max(EE_max);
    
    
    for i=1:8
    for j=1:length(Power(i,:))
        if(Power(i,j) < 0.0000001)
            continue;
        end
        
         SE(i,j)=SE(i,j)/SE_max_global;
         EE(i,j)=EE(i,j)/EE_max_global;
        
        f(i,j)= k*SE(i,j)+(1-k)*EE(i,j);
    end
end
for j=1:2479
    [max_SE(j) I_SE(j)]=max(SE(:,j));
   [max_EE(j) I_EE(j)]=max(EE(:,j));
   [max_f(j) I_f(j)]=max(f(:,j));
end
   

subplot(4,2,6)
plot(nonzeros(Power(1,:))*1000,nonzeros(max_f(:)),'LineWidth',2)
axis([0 250 0.3 0.72])
title('WEF-based Optimal Mode Selection: k=0.65','FontWeight','Bold')
xlabel('Power (mWatts)')
ylabel('F^{SE,EE}_{max}');

subplot(4,2,8)
plot(nonzeros(Power(1,:))*1000,nonzeros(I_f(70:2479)),'LineWidth',2)
axis([0 250 0 8])
xlabel('Power (mWatts)')
ylabel('Mode index');
hold on;

% [ppp qqq]=max(max_f(:));
% ppp
% I_f(qqq)
% Power(I_f(qqq),qqq)*1000
    
% for i=1:8
%     figure(3)
%     plot(nonzeros(Power(i,:)*1000),nonzeros(f(i,:)))
%     hold on;
% end

