%% test complex space
clear all;
close all;
%% due to the limit of file size on github, the data is just for 10 loops
KH = load('Channel_kernel_N_10.mat'); %N=10 case
hat_Phi = load('hat_phi_N_10.mat');
hat_Psi = load('hat_psi_N_10.mat');

KH = KH.Channel_kernel_N_10;
hat_Phi = hat_Phi.hat_phi_N_10;
hat_Psi = hat_Psi.hat_psi_N_10;
% 4x5x4x5 channel kernel
Nu = 4;
T = 5;
Bw = 20*10^6; % 20MHz bw
Te = T/Bw; % duration time for eigenfunctions.
Nloop = size(KH,5);
Nep = size(hat_Phi,3);
N = size(hat_Phi,2);
M = 4;
K = log2(M);
sq2 = sqrt(2);

SNRdB = [0:1:30];
taps_max = 5;
taps_min = 3;
packets = 1;

BER_total = zeros(length(SNRdB),Nep);
O = zeros(length(SNRdB),Nep);

%%
for m = 1:Nloop
   hat_phi = hat_Phi(:,:,:,m);
   hat_psi = hat_Psi(:,:,:,m);
   Hst = KH(:,:,:,:,m);
   for n = 1:length(SNRdB)
      sigma2 = 10^(-SNRdB(n)/10); sigma = sqrt(sigma2);
      ber_temp = zeros(Nep,1);
      O_temp = zeros(Nep,1);
  
      data = randi([0 1],N,K*packets);
      s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
      s = s';
      Tx = zeros(Nu,T,packets,Nep);  
        for k = 1:Nep
         O1 = soft_orth(hat_phi(:,:,k),N);
         O2 = soft_orth(hat_psi(:,:,k),N);   
         O_temp(k) = O_temp(k) + O2;

           for i = 1:N
              for j = 1:packets
                 temp_phi = hat_phi(:,i,k);
                 Phi = reshape(temp_phi, Nu,T);
                 Tx(:,:,j,k) = Tx(:,:,j,k) + Phi*s(i,j); %multiplexing
              end
           end
        end
   
        
        %%
           Rx = zeros(Nu, T, packets,Nep);
           for k = 1:Nep
               for i = 1:Nu
                  for j = 1:T
                     for l = 1:packets
                         Rx(i,j,l,k) = sum(sum(squeeze(Hst(i,j,:,:)) .* Tx(:,:,l,Nep))) + sigma*(randn()+randn()*1i)/sq2; %integral 
                     end
                  end
               end      
           end
        %%
        
        for k = 1:Nep
           r = zeros(N,packets);
           for i = 1:N
              for j = 1:packets
                 temp_psi = hat_psi(:,i,k);
                 Psi = reshape(temp_psi,Nu,T);
                 r(i,j) = sum(sum(conj(Psi).*Rx(:,:,j,k))); %demultiplexing
              end
           end
           p = mean(abs(r').^2);
           hat_s = qamdemod(r'./p,M,'OutputType','bit','UnitAveragePower',true);
           hat_s = hat_s';
           ber_temp(k) = ber_temp(k)+ length(find(round(data - hat_s)))/(N*K*packets);
        end
   BER_total(n,:) = BER_total(n,:) + ber_temp';
   O(n,:) = O(n,:) + O_temp';
   end
end
%%
BER_total = BER_total/Nloop;
O = O/Nloop;
%% BER
f1 = figure;
plot_ber = repmat(SNRdB',1,Nep);
plot_O = O;
plt = surf(plot_ber,plot_O,BER_total);
colorbar;
set(gca,'fontsize',26);
colormap('hot');
caxis([0,0.2]); 
ylim([0 0.05])
view(0,90);
plt.EdgeColor = 'none';
xlabel('SNR [dB]')
ylabel('Soft Orthogonality $$O(\hat{\Psi})$$','Interpreter','Latex')
zlabel('BER')
