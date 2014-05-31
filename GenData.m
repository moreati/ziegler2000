% Build encoder
%[n,k,dim,G,H] = Hamming(3)
%[n,k,dim,G,H] = Hamming(4)
%[n,k,M,G] = Golay23(1)          % g = [1 1 0 0 0 1 1 1 0 1 0 1], t=3
%[n,k,M,G] = Golay23(2)          % g = [1 0 1 0 1 1 1 0 0 0 1 1], t=3
%[n,k,M,q,G] = CRC(7,4,2,[1 0 1 1])    % note: also a Hamming & BCH code
%[n,k,M,q,G] = CRC(2047,2035,2,[1 1 0 0 0 0 0 0 0 1 1 1 1]) % CRC-12
%[n,k,M,q,G] = BCH(63,45,3,2,[1 1 1 1 0 0 0 0 0 1 0 1 1 0 0 1 1 1 1])
%[n,k,M,q,G] = BCH(21,15,1,4,[1 0 1 0 1 1 1]);
%[n,k,M,q,G] = BCH(21,12,2,4,[1 1 1 0 1 1 0 0 1 1]);
%[n,k,M,q,G] = BCH(85,73,3,16,[1 4 11 10 2 4 11 15 10 15 8 5 4]);
% Reed-Solomon
%[n,k,M,q,G] = BCH(63,57,3,64,[1 60 49 44 56 11 22]);
%[n,k,M,q,G] = BCH(255,245,5,256,[1 254 72 53 70 129 83 79 111 51 66]);
%[n,k,M,q,G] = ReedMuller(1,4)      % (16,5) binary Reed-Muller, t=3
%[n,k,M,q,G] = ReedMuller(1,5)      % (32,6) binary Reed-Muller, t=7
%G = [1 0 1; 1 1 1];                      % rate 1/2, K=3 conv code
%G = [1 1 0 1;1 1 1 1];                   % rate 1/2, K=4 conv code
%G = [1 0 0 1 1;1 1 1 0 1];               % rate 1/2, K=5 conv code
%G = [1 1 0 1 0 1;1 0 1 1 1 1];           % rate 1/2, K=6 conv code
%G = [1 0 1 1 0 1 1;1 1 1 1 0 0 1];       % rate 1/2, K=7 conv code
%G = [1 1 1 0 0 1 0 1;1 0 0 1 1 1 1 1];   % rate 1/2, K=8 conv code
%G = oct2gen([561;753]);                  % rate 1/2, K=9 conv code
%G = oct2gen([565;543]);                  % rate 1/2, K=9 (suboptimal)
%G = oct2gen([43572;56246])               % rate 1/2, K=14 conv code
%G = [1 0 1; 1 1 1; 1 1 1];               % rate 1/3, K=3 conv code
%G = [1 0 1 1; 1 0 0 1; 1 1 1 1];         % rate 1/3, K=4 (suboptimal)
%G = oct2gen([452;662;756])               % rate 1/3, K=8 conv code
%G = oct2gen([43512;73542;76266])         % rate 1/3, K=14 conv code
%G = [1 0 1; 1 1 1; 1 1 1; 1 1 1];        % rate 1/4, K=3 conv code
%G = oct2gen([564;564;634;714]);          % rate 1/4, K=7 conv code
%G = oct2gen([42226;46372;73256;73276])   % rate 1/4, K=14 conv code
%G = oct2gen([7;1;4; 2;5;7])           % rate 2/3, K=3, M=4 conv code
%G = oct2gen([63;15;46; 32;65;61])     % rate 2/3, K=6, M=10 conv code
% rate 3/4, K=4, M=9 conv code
G = oct2gen([40;14;34;60; 04;64;20;70; 34;00;60;64])

% Generate coded   bit stream
%[m,c,WordCount]   = BlockEnc('',100,1,n,k,q,G,0)    % no file output
%[m,c,WordCount]   = BlockEnc('Hamming3',1000,1000,n,k,q,G,1);
%[m,c,WordCount]   = BlockEnc('Hamming4',1000,100,n,k,q,G,1);
%[m,c,WordCount]   = BlockEnc('Golay23-1',1000,100,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('Golay23-2',1000,100,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('crc7-4-1',1000,1,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('bch63-45-3',1000,10,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('bch21-15-1-4',1000,5,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('bch21-12-2-4',1000,5,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('bch85-73-3-16',1000,5,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('rs63-57-3',100,14,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('rs255-245-5',100,14,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('rm16-5-3',1000,100,n,k,q,G,1);
%[m,c,WordCount] = BlockEnc('rm32-6-7',1000,10,n,k,q,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K3',100000,2,1,3,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K4',100000,2,1,4,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K5',100000,2,1,5,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K6',100000,2,1,6,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K7',100000,2,1,7,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K8',100000,2,1,8,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K9',10000,2,1,9,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K9b',10000,2,1,9,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-2_K14',10000,2,1,14,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-3_K3',10000,3,1,3,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-3_K4b',10000,3,1,4,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-3_K8',10000,3,1,8,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-3_K14',10000,3,1,14,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-4_K3',10000,4,1,3,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-4_K7',10000,4,1,7,G,1);
%[m,c,WordCount] = ConvEnc2('ConvR1-4_K14',10000,4,1,14,G,1);
%[m,c,WordCount] = ConvEnc3('ConvR2-3_K3_M4',10000,3,2,3,G,1);
%[m,c,WordCount] = ConvEnc3('ConvR2-3_K6_M10',10000,3,2,6,G,1);
[m,c,WordCount] = ConvEnc3('ConvR3-4_K4_M9',10000,4,3,4,G,1);




