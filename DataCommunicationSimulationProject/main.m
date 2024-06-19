clear all;
clc;

%Number of level
level=32;

%Loading voice
[x,fm]=audioread('bir.wav');

%Fundamental frequency
N=floor(0.02*fm);
C=xcorr(x,N,'coeff');
N1=floor(0.002*fm);
[x0,vmax]=max(C(N+N1:2*N+1));
t0=(vmax+N1)/fm;
f0=1/t0;
fundamental_frequency=strcat(num2str(f0),' Hz');

%PLOT
%Plotting input signals (voices)
figure(1)
plot(x)
axis([ 0 4500 min(x) max(x) ])
title('Input signal 1');

% Playing voices
disp('Playing input signals');
soundsc(x,fm);
pause(3);

%%
%Quantization
[y1, x2, errorcuantizacion] = quantize(x,level);

%Quantization error
quantization_error = strcat(num2str(errorcuantizacion),' %');

%Variables to plot
xg=x2; yq=y1;
xq=x; fmq=fm;

%PLOT
%Plotting input signal with level of quantization
figure(2)
subplot(2,1,1)
plot(x);
axis([ 0 4500 min(x) max(x) ])
hold on
for vv=1:level
    hold on
    plot(yq(vv,:));
end
grid on
xlabel('samples')
ylabel('x(t)')
title('Input signal with level')

%Ploting input signal quantized
subplot(2,1,2)
plot(xq)
axis([ 0 4500 min(xq) max(xq)])
grid on
xlabel('samples')
title('Input signal quantized')

%%
%--------------------- NO CODIFICATION ------------------------------------
%Transform decimal to bin
x3=x2-1;
bits=dec_bin(x3,log2(level));

%Matrix to vector
tem=[];
for i=1:size(bits,1)
    for j=1:size(bits,2)
        tem=[tem bits(i,j)];
    end
end
bitsc1=tem;

%--------------------- HAMMING 7,4 ----------------------------------------
% Dimensiones y matrices
n=7;
k=4;

%P matrix
P=[1 1 0; 0 1 1; 1 1 1; 1 0 1];

%Matrix generator
identity=eye(k);
G=[P identity];

%Variables to divide messages in packages
tamanio=size(bitsc1,2); %Size of message
div=1;
matrizc=[];

while(div<tamanio)
    %Divide message in packets
    m=bitsc1(div:div+3);
    div=div+4;

    %Matrix c=mG
    c=mod(m*G,2);
    matrizc=[matrizc c];
end



%%
% Array to plot constellation
        constellationArray=[-1,0;1,0];
        modulation_name = 'BPSK';

%PLOT
%Plot axis
figure(3)
plot([-2 2],[0 0],'b');
hold on
plot([0 0 ],[-2 2],'b');

%Plot constellation
for i=1:size(constellationArray,1)
    hold on
    p=plot(constellationArray(i,1),constellationArray(i,2),'*');
    set(p,'Color','red','LineWidth',2);
end

title(strcat('Constellation Diagram: ', modulation_name))
axis([-2 2 -2 2])

%-MODULATION WITH HAMMING CODE
bitsm2=2*matrizc-1;

%%
%----------------------- BER CURVES ---------------------------------------
%Probability of error
pet2=[]; %Hamming 7,4


%Data for Hamming codes
n=7;
k=4;

%Matriz x (Hamming)
P=[1 1 0; 0 1 1; 1 1 1; 1 0 1];

%Generator matrix (Hamming)
Identidad=eye(k);
G=[P Identidad];

%Matrix H (Hamming)
Identidad=eye(n-k);
H=[Identidad P'];

%Syndrome decoding table (Hamming)
Identidad=eye(n);
t1=[zeros(1,n); Identidad];
tabla=t1*H';

%Generator matrix for Convolutional codes
g=[1 0 1;1 1 1];

% Eb/N0 The energy per bit to noise power spectral density ratio
% 1 to 6 where 6 is the least noisy
ebn0db=0:1:6;

%--------------------- LOOP TO PLOT BER -----------------------------------
for eb=ebn0db
    %AWGN Channel
    ebn0=10^(eb/10);
    sigma=1/sqrt(2*ebn0);

    %Variables to caculate error (Hamming 7, 4)
    error2=0;
    pe2=[];
    

        for numiter=1:5
            %--------------- HAMMING 7,4 ----------------------------------
            N=size(bitsm2,2);

            % CALCULATE FOR BPSK MODULATION
                %AWGN Channel
                ruido=normrnd(0,sigma,1,N);
                y2=bitsm2+ruido;

                %BPSK demodulation
                bitsmr2=sign(y2);
                bitsr2=(bitsmr2+1)/2;

                %Variables to divide message in packets
                % message + parity bit
                tamanio=size(bitsr2,2);
                div=1; %message + parity bit

                bitsr=[]; %decoded matrix - only message

                while(div<tamanio)
                    %Divide message into packets
                    bits2=bitsr2(div:div+6);

                    %Syndrom
                    sindrome=mod(bits2*H',2);

                    %Correcction
                    for i=1:(n+1)
                        if tabla(i,:)==sindrome
                            if i~=1
                                bits2(i-1)=~bits2(i-1);
                            end
                        end
                    end

                    %Message decoded
                    b=bits2(4:7);
                    bitsr=[bitsr b];

                   %Increase div
                    div=div+7;
                end

                %BER - Hamming 7,4 
                error2=error2+sum(xor(bitsc1,bitsr));
                pe2=[pe2 sum(xor(bitsc1,bitsr))/N];
          
                
            end
        
            pet2=[pet2 mean(pe2)];
      
        
end
%Pe error (Probability of error)
    errorpe_bpsk_hamming=mean(pet2)
   

%PLOT
%BER according modulation
figure(4)
        semilogy(ebn0db,pet2,'rx-');
        hold on;
        xlabel('Eb/N0, dB')
        ylabel('Bit Error Rate')
        title('BER Curves (BPSK)')
        ylim([10^(-4) 10^(-1)]);
        hleg = legend('Hamming');
        set(hleg,'Location','EastOutside')
        grid on;
