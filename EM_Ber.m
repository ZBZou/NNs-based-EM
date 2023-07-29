%% test complex space
clear all;
close all;

% 4x5x4x5 channel kernel
Nu = 4;
T = 5;
Bw = 20*10^6; % 20MHz bw
Te = T/Bw; % duration time for eigenfunctions.

M = 4;
K = log2(M);
sq2 = sqrt(2);

Nloop = 100;
SNRdB = [0:1:30];
ep = [0:0.01:3];
taps_max = 5;
taps_min = 3;
packets = 100;
N = 15; 

BER_total = zeros(length(SNRdB),length(ep));
O = zeros(length(SNRdB),length(ep));

%%
for n = 1:length(SNRdB)
   sigma2 = 10^(-SNRdB(n)/10); sigma = sqrt(sigma2);
   sq2 = sqrt(2);
   Nep = length(ep);
   ber_temp = zeros(Nep,1);
   O_temp = zeros(Nep,1);
   n
   for m = 1:Nloop
     
      Taps = randi([taps_min,taps_max],Nu,T);
      data = randi([0 1],N,K*packets);
      s = qammod(data',M,'InputType','bit','UnitAveragePower',true);
      s = s';
      Hst = zeros(Nu,T, Nu, T); %4x10x4x10 kernel
        for u1 = 1:Nu
            for t1 = 1:T
                for u2 = 1:Nu
                   for t2 = 1:T
                      taps = Taps(u2,t1);
                      if t2 > t1 || t2 < t1 - taps
                         continue
                      end
                      vi = sqrt(t2/T)*rand();
                      ei = sqrt(t2/T)*(rand()+rand()*1i)/sq2;
                      Hst(u1,t1,u2,t2) = ei+vi*(randn()+randn()*1i)/sq2;
                   end
                end 
            end 
        end 
        
        H2d = reshape(Hst, Nu*T,Nu*T);
        [U,Sig,V] = svd(H2d);
        phi = V(:,1:N);
        psi = U(:,1:N);
        sig = diag(Sig(1:N));
        
        Tx = zeros(Nu,T,packets,Nep);
        hat_psi = zeros(Nu*T,N,Nep);
        for k = 1:Nep
            ep_k = ep(k);
            hat_phi = approx_eigenfunc(phi,ep_k);
            hat_psi(:,:,k) = dual(hat_phi,H2d);
            O1 = soft_orth(hat_phi,N);
            O2 = soft_orth(hat_psi(:,:,k),N);   
            O_temp(k) = O_temp(k) + O2;
           for i = 1:N
              for j = 1:packets
                 temp_phi = hat_phi(:,i);
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
   end
   BER_total(n,:) = ber_temp'/Nloop;
   O(n,:) = O_temp'/Nloop;
end
%%
data_size = ones(size(BER_total,1),size(BER_total,2));
Thp = (K*N*data_size-K*N*BER_total)/(Te*10^6); %throughput
%% BER
f1 = figure;
plot_ber = repmat(SNRdB',1,Nep);
plot_O = 0;
plt = surf(plot_ber,plot_O,BER_total);
colorbar;
set(gca,'fontsize',26);
colormap('hot');
ylim([0 0.05])
view(0,90);
plt.EdgeColor = 'none';
xlabel('SNR [dB]')
ylabel('Soft Orthogonality $$O(\hat{\Psi})$$','Interpreter','Latex')
zlabel('BER')
%%
savefile = 0;
if savefile ==1
   save Ber_4x5_N20_ep3_4QAM_2d.mat SNRdB plot_O BER_total -v7.3
end
%%
savefig = 0;
if savefig == 1
filename = 'Ber_4x5_N20_ep3_4QAM_2d';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f1, name1);
exportgraphics(f1, name2);
end

%% thp
f2 = figure;

plt = surf(plot_ber,plot_O,Thp);
plt.EdgeColor = 'none';
view(55,45);
ylim([0 0.05])
plt.EdgeColor = 'none';
set(gca,'fontsize',26);
xlabel('SNR [dB]','Rotation',-35,'VerticalAlignment','middle');
ylabel('Soft Orthogonality $$O(\hat{\Psi})$$','Interpreter','Latex','Rotation',15,'VerticalAlignment','middle')
zlabel('Throughput [Mbps]')

%%
savefile = 1;
if savefile ==1
   save Th_4x5_N20_ep3_4QAM_3d.mat SNRdB plot_O Thp -v7.3
end
%%
savefig = 0;
if savefig == 1
filename = 'Th_4x5_N20_ep3_4QAM_3d';
name1 = append(filename, '.fig');
name2 = append(filename, '.pdf');
saveas(f2, name1);
exportgraphics(f2, name2);
end


