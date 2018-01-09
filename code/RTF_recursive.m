function [y_out,hs_out, v_out] = RTF_recursive (x_TD, hs_TD, v_TD, fs, beta)

% MVDR beamformer speech enhancement with recursive algorithm for RTF
% estimation based on the GEVD 
% Stijn Rutten

%INPUT:
%   x_TD            : microphone signal
%   hs_TD           : clean speech signal
%   v_TD            : noise signal
%   beta            : forgetting factor
%
% OUTPUT:
%   y_out       : output fileterd microphone signal
%   hs_out       : output fileterd speech signal
%   v_out       : output fileterd noise signal

%% STFT 
% var of STFT
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
y_out = zeros(bins,frames);
v_out = zeros(bins,frames);
hs_out = zeros(bins,frames);
w = zeros(channels,bins);


%init

RTF(:,:) = randn(channels,bins) +randn(channels,bins)*1i;
RTF(1,1) =1;


for i=2 :frames   
    for j =1 : bins
        
       
        if mag2db(abs(x_STFT(j,i,1))) > -17
            Ryy(:,:,j)= (beta.*squeeze( Ryy(:,:,j))) + ((1-beta).* squeeze(y_STFT(j,i,:))*squeeze(y_STFT(j,i,:))') ;
     
        else
            Rnn(:,:,j)= (beta.*squeeze( Rnn(:,:,j))) + ((1-beta).* squeeze(y_STFT(j,i,:))*squeeze(y_STFT(j,i,:))') ;
            
        end
        
        if mag2db(abs(x_STFT(j,i,1))) > -17
            %step1 Time update
            A =  Ryy(:,:,j) / (Rnn(:,:,j)) ;
            kappa = (e1.' * A *  RTF(:,j));
            A= (1/kappa)*A;
            RTF(:,j) =   A*RTF(:,j);
       
        end
        
        w(:,j) = ((Rnn(:,:,j) +epsilon * eye(channels))\ RTF(:,j) ) /  ( (RTF(:,j)'/(Rnn(:,:,j) +epsilon * eye(channels))*RTF(:,j)));
        y_out(j,i) = w(:,j)' * squeeze(y_STFT(j,i,:));
        hs_out(j,i) = w(:,j)' * squeeze(x_STFT(j,i,:));
        v_out(j,i) = w(:,j)' * squeeze(v_STFT(j,i,:));
     end
end