function [G2] = FindGenPolys(FileName,n,K)

% FIND CONVOLUTIONAL CODE GENERATOR POLYNOMIALS VIA DATA CANCELLATION

WordCount = 20*K;
N = WordCount+(n-1)*(K-1);

EncodedBitCount = n*WordCount;
s = sprintf('%s.fec',FileName);
CodedFile = fopen(s,'r');

[r,WordsRead] = fread(CodedFile,EncodedBitCount,'ubit1');
if(WordsRead ~= EncodedBitCount)
   disp 'Error reading file'




   WordsRead
end

r2 = reshape(r,n,WordCount);

perm = prod([1:n])      % n! = # of possible orderings of gi(x)
PermLUT = perms(1:n);         % table of index permutations for gi(x)
%PermLUT = fliplr(perms(1:n)) % table of index permutations for gi(x)
comb = 2^(n*K)          % # of tap combinations across all polynomials
MaxComb = 1e6;          % maximum # of tap combinations to search


% Search based on best known polynomials first
if n==2

    switch K
       case 3,
          G = oct2gen([5;7]);
       case 4,
          G = oct2gen([64;74]);
%         G = oct2gen([74;64]);
       case 5,
          G = oct2gen([46;72]);
       case 6,
          G = oct2gen([65;57]);
       case 7,
          G = oct2gen([554;744]);
       case 8,
          G = oct2gen([712;476]);
       case 9,
          G = oct2gen([561;753]);
       case 10,
          G = oct2gen([4734;6624]);
       case 11,
          G = oct2gen([4672;7542]);
       case 12,
          G = oct2gen([4335;5723]);
       case 13,
          G = oct2gen([42554;77304]);
       case 14,
          G = oct2gen([43572;56246]);
    end

    % Try possible orderings until a winner is found
    fprintf('Searching for optimal generators\n')
    for i=1:perm
       a1 = G(i,:);
       a2 = G(perm+1-i,:);
       P1 = gfconv(a1,r2(1,:));
       P1 = [P1 zeros(1,N-length(P1))];
       P2 = gfconv(a2,r2(2,:));
       P2 = [P2 zeros(1,N-length(P2))];
       S = mod(P1+P2,2);
       S = S(1:WordCount);
       if sum(S)==0
          fprintf('Found generators at trial %d:\n',i)




          G2 = oct2gen([a2;a1],[1,n,K-1])
          break
       end
    end

    if (sum(S)~=0) & (comb<MaxComb)
       % Try all generator tap combinations, except when any gi(x)=0
       fprintf('Searching suboptimal generators\n')
       mask = 2^K-1;
       imin = 2^K;
       j = 0;
       for i=imin:comb-1
          g1o = bitand(bitshift(i,-K),mask);
          g2o = bitand(i,mask);
          if (g1o>0) & (g2o>0)
             j = j+1;
             a1 = de2bi(g1o,K);
             a2 = de2bi(g2o,K);
             P1 = gfconv(a1,r2(1,:));
             P1 = [P1 zeros(1,N-length(P1))];
             P2 = gfconv(a2,r2(2,:));
             P2 = [P2 zeros(1,N-length(P2))];
             S = mod(P1+P2,2);
             S = S(1:WordCount);
             if sum(S)==0
                fprintf('Found generators at i=%d (trial %d):\n',i,j)
                G2 = oct2gen([a2;a1],[1,n,K-1])
                break
             end
          end
       end
    end

    if sum(S)~=0
       fprintf('No generators found\n')
    end

elseif n==3

    switch K
       case 3,
          G = oct2gen([5;7;7]);
%          G = oct2gen([7;5;7]);
%          G = oct2gen([7;7;5]);
       case 4,
          G = oct2gen([54;64;74]);
       case 5,
          G = oct2gen([52;66;76]);
       case 6,
          G = oct2gen([47;53;75]);
       case 7,
          G = oct2gen([554;624;764]);
       case 8,
          G = oct2gen([452;662;756]);
       case 9,
          G = oct2gen([557;663;711]);




       case   10,
          G   = oct2gen([4474;5724;7154]);
       case   11,
          G   = oct2gen([4726;5562;6372]);
       case   12,
          G   = oct2gen([4767;5723;6265]);
       case   13,
          G   = oct2gen([42554;43364;77304]);
       case   14,
          G   = oct2gen([43512;73542;76266]);
%         G   = oct2gen([73542;43512;76266]);
    end

    % Try possible orderings until a winner is found
    fprintf('Searching for optimal generators\n')
    for i=1:perm
       g1 = G(PermLUT(i,1),:);
       g2 = G(PermLUT(i,2),:);
       g3 = G(PermLUT(i,3),:);
       a1 = gfconv(g2,g3);
       a2 = gfconv(g1,g3);
       a3 = gfconv(g1,g2);
       P1 = gfconv(a1,r2(1,:));
       P1 = [P1 zeros(1,N-length(P1))];
       P2 = gfconv(a2,r2(2,:));
       P2 = [P2 zeros(1,N-length(P2))];
       P3 = gfconv(a3,r2(3,:));
       P3 = [P3 zeros(1,N-length(P3))];
       S = mod(P1+P2+P3,3);
       S = S(1:WordCount);
       if sum(S)==0
          fprintf('Found generators at trial %d:\n',i)
          G2 = oct2gen([g1;g2;g3],[1,n,K-1])
          break
       end
    end

    if (sum(S)~=0) & (comb<MaxComb)
       % Try all generator tap combinations, except when any gi(x)=0
       fprintf('Searching suboptimal generators\n')
       mask = 2^K-1;
       imin = 2^((n-1)*K);
       j = 0;
       for i=imin:comb-1
          g1o = bitand(bitshift(i,-2*K),mask);
          g2o = bitand(bitshift(i,-K),mask);
          g3o = bitand(i,mask);
          if (g1o>0) & (g2o>0) & (g3o>0)
             j = j+1;
             g1 = de2bi(g1o,K);
             g2 = de2bi(g2o,K);
             g3 = de2bi(g3o,K);
             a1 = gfconv(g2,g3);
             a2 = gfconv(g1,g3);
             a3 = gfconv(g1,g2);
             P1 = gfconv(a1,r2(1,:));




             P1 = [P1 zeros(1,N-length(P1))];
             P2 = gfconv(a2,r2(2,:));
             P2 = [P2 zeros(1,N-length(P2))];
             P3 = gfconv(a3,r2(3,:));
             P3 = [P3 zeros(1,N-length(P3))];
             S = mod(P1+P2+P3,3);
             S = S(1:WordCount);
             if sum(S)==0
                fprintf('Found generators at i=%d (trial %d):\n',i,j)
                G2 = oct2gen([g1;g2;g3],[1,n,K-1])
                break
             end
          end
       end
    end

    if sum(S)~=0
       fprintf('No generators found\n')
    end

elseif n==4

    switch K
       case 3,
          G = oct2gen([5;7;7;7]);
%          G = oct2gen([7;7;5;7]);
       case 4,
          G = oct2gen([54;64;64;74]);
       case 5,
          G = oct2gen([52;56;66;76]);
       case 6,
          G = oct2gen([53;67;71;75]);
       case 7,
          G = oct2gen([564;564;634;714]);
%          G = oct2gen([564;634;714;564]);
       case 8,
          G = oct2gen([472;572;626;736]);
       case 9,
          G = oct2gen([463;535;733;745]);
       case 10,
          G = oct2gen([4474;5724;7154;7254]);
       case 11,
          G = oct2gen([4656;4726;5562;6372]);
       case 12,
          G = oct2gen([4767;5723;6265;7455]);
       case 13,
          G = oct2gen([44624;52374;66754;73534]);
       case 14,
          G = oct2gen([42226;46372;73256;73276]);
    end

    % Try possible orderings until a winner is found
    fprintf('Searching for optimal generators\n')
    for i=1:perm
       g1 = G(PermLUT(i,1),:);
       g2 = G(PermLUT(i,2),:);




   g3 = G(PermLUT(i,3),:);
   g4 = G(PermLUT(i,4),:);
   a12 = gfconv(g1,g2);
   a13 = gfconv(g1,g3);
   a23 = gfconv(g2,g3);
   a1 = gfconv(a23,g4);
   a2 = gfconv(a13,g4);
   a3 = gfconv(a12,g4);
   a4 = gfconv(a12,g3);
   P1 = gfconv(a1,r2(1,:));
   P1 = [P1 zeros(1,N-length(P1))];
   P2 = gfconv(a2,r2(2,:));
   P2 = [P2 zeros(1,N-length(P2))];
   P3 = gfconv(a3,r2(3,:));
   P3 = [P3 zeros(1,N-length(P3))];
   P4 = gfconv(a4,r2(4,:));
   P4 = [P4 zeros(1,N-length(P4))];
   S = mod(P1+P2+P3+P4,4);
   S = S(1:WordCount);
   if sum(S)==0
      fprintf('Found generators at trial %d:\n',i)
      G2 = oct2gen([g1;g2;g3;g4],[1,n,K-1])
      break
   end
end

if (sum(S)~=0) & (comb<MaxComb)
   % Try all generator tap combinations, except when any gi(x)=0
   fprintf('Searching suboptimal generators\n')
   mask = 2^K-1;
   imin = 2^((n-1)*K);
   j = 0;
   for i=imin:comb-1
      g1o = bitand(bitshift(i,-3*K),mask);
      g2o = bitand(bitshift(i,-2*K),mask);
      g3o = bitand(bitshift(i,-K),mask);
      g4o = bitand(i,mask);
      if (g1o>0) & (g2o>0) & (g3o>0) & (g4o>0)
         j = j+1;
         g1 = de2bi(g1o,K);
         g2 = de2bi(g2o,K);
         g3 = de2bi(g3o,K);
         g4 = de2bi(g4o,K);
         a12 = gfconv(g1,g2);
         a13 = gfconv(g1,g3);
         a23 = gfconv(g2,g3);
         a1 = gfconv(a23,g4);
         a2 = gfconv(a13,g4);
         a3 = gfconv(a12,g4);
         a4 = gfconv(a12,g3);
         P1 = gfconv(a1,r2(1,:));
         P1 = [P1 zeros(1,N-length(P1))];
         P2 = gfconv(a2,r2(2,:));
         P2 = [P2 zeros(1,N-length(P2))];
         P3 = gfconv(a3,r2(3,:));
         P3 = [P3 zeros(1,N-length(P3))];




            P4 = gfconv(a4,r2(4,:));
            P4 = [P4 zeros(1,N-length(P4))];
            S = mod(P1+P2+P3+P4,4);
            S = S(1:WordCount);
            if sum(S)==0
               fprintf('Found generators at i=%d (trial %d):\n',i,j)
               G2 = oct2gen([g1;g2;g3;g4],[1,n,K-1])
               break
            end
         end
      end
   end

   if sum(S)~=0
      fprintf('No generators found\n')
   end

end




