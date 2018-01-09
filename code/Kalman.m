function [y_out,hs_out, v_out] = Kalman(x_TD, hs_TD, v_TD, fs, beta, process_noise )

% MVDR speech enhancement with a Kalman filter for RTF estimation
%Stijn Rutten

%INPUT:
%   x_TD            : microphone signal
%   hs_TD           : clean speech signal
%   v_TD            : noise signal
%   beta            : forgetting factor
%   process_noise   : added process noise in dB
%
% OUTPUT:
%   y_out       : output fileterd microphone signal
%   hs_out      : output fileterd speech signal
%   v_out       : output fileterd noise signal

%% STFT 

N_fft = 512;
% frame shift
R_fft = N_fft/2;
% analysis window
win = sqrt(hann(N_fft,'periodic'));
% number of bins in onsided FFT 
N_half = floor(N_fft/2)+1;

x_STFT = calc_STFT(hs_TD, fs, win, N_fft, N_fft/R_fft, 'onesided'); % speech component (not available in practice, but just to have a look at it)
v_STFT = calc_STFT(v_TD, fs, win, N_fft, N_fft/R_fft, 'onesided'); % noise component (not available in practice, but you can use it to compute the noise covariance matrix)
y_STFT = calc_STFT(x_TD, fs, win, N_fft, N_fft/R_fft, 'onesided'); % microphone signals

%figure; imagesc(mag2db(abs(x_STFT(end:-1:1,:,1))), [-65, 15]); colorbar; title('speech component, 1st mic');
%figure; imagesc(mag2db(abs(v_STFT(end:-1:1,:,1))), [-65, 15]); colorbar; title('noise component, 1st mic');
%figure; imagesc(mag2db(abs(y_STFT(end:-1:1,:,1))), [-65, 15]); colorbar; title('microphne signal, 1st mic');

bins = N_half; 
frames = size(y_STFT,2);
channels= size(y_STFT, 3);


%% algoritme
% var

Ryy = zeros(channels,channels,bins);
Rnn = zeros(channels,channels,bins);
epsilon = 1* 10^-6;

e1 = zeros(channels,1);
e1(1,1) = 1;
RTF = zeros( channels,bins);
RTF( 1,:) = 1;
y_out_STFT = zeros(bins,frames);
v_out_STFT = zeros(bins,frames);
x_out_STFT = zeros(bins,frames);
w = zeros(channels,bins);

P = zeros(channels-1,channels-1, bins);
for k= 1: bins
    P(:,:,k) = eye(channels-1,channels-1);
end
B = zeros ( channels, channels-1);

%init
Rww = db2pow(process_noise) .* eye(channels-1, channels-1);
RTF(2:channels ,:) = randn(channels-1,bins) + randn(channels-1,bins)*1i;

for i=2 :frames   
    for j =1 : bins
        
        if mag2db(abs(x_STFT(j,i,1))) > -17
            Ryy(:,:,j)= (beta.*squeeze( Ryy(:,:,j))) + ((1-beta).* squeeze(y_STFT(j,i,:))*squeeze(y_STFT(j,i,:))') ;
        else
            Rnn(:,:,j)= (beta.*squeeze( Rnn(:,:,j))) + ((1-beta).* squeeze(y_STFT(j,i,:))*squeeze(y_STFT(j,i,:))') ;
            
        end
        
        if mag2db(abs(x_STFT(j,i,1))) > -17 %speech
            %step1 Time update
            A =  Ryy(:,:,j)  / Rnn(:,:,j) ;
            kappa = (e1.' * A *  RTF(:,j));
            A= (1/kappa).* A;
            RTF(2:channels,j) =   A(2:channels,2:channels)*RTF(2:channels,j) + A(2:channels,1);
            
            P(:,:,j) = A(2:channels,2:channels) * P(:,:,j) * (A(2:channels,2:channels))'+ Rww(:,:);
                 
            
            %step2 Measurement update
            B(1,:) = (-(RTF(2:channels,j))');
            B( 2:channels,:) = eye(channels-1);
            u = B' * squeeze(y_STFT(j,i,:));
            Ru = ((y_STFT(j,i,1)* y_STFT(j,i,1)' )*  P(:,:,j)) + B' * Rnn(:,:,j)* B;
            K = (y_STFT(j,i,1)' * P(:,:,j))/ (Ru );
            RTF(2:channels,j) = RTF(2:channels,j) + K * u;
            P(:,:,j) =  P(:,:,j) - y_STFT(j,i,1)*K*P(:,:,j);
            
            P(:,:,j) = 0.5*(P(:,:,j) + P(:,:,j)');
        end
        

        w(:,j) = ((Rnn(:,:,j) +epsilon * eye(channels))\ RTF(:,j) ) /  ( (RTF(:,j)'/(Rnn(:,:,j) +epsilon * eye(channels))*RTF(:,j))); 
        y_out(j,i) = w(:,j)' * squeeze(y_STFT(j,i,:));
        hs_out(j,i) = w(:,j)' * squeeze(x_STFT(j,i,:));
        v_out(j,i) = w(:,j)' * squeeze(v_STFT(j,i,:));
     end
end
       
        
   

