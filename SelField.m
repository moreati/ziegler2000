function [vectors,symbols] = SelField(qm,mx)

global   vectors8_m1;
global   symbols8_m1;
global   vectors8_m2;
global   symbols8_m2;
global   vectors16_m1;
global   symbols16_m1;
global   vectors32_m1;
global   symbols32_m1;
global   vectors64_m1;
global   symbols64_m1;
global   vectors128_m1;
global   symbols128_m1;
global   vectors256_m1;
global   symbols256_m1;
global   vectors512_m1;
global   symbols512_m1;
global   vectors1K_m1;
global   symbols1K_m1;
global   vectors2K_m1;
global   symbols2K_m1;
global   vectors2K_m2;




global   symbols2K_m2;
global   vectors4K_m1;
global   symbols4K_m1;
global   vectors16K_m1;
global   symbols16K_m1;
global   vectors32K_m1;
global   symbols32K_m1;
global   vectors64K_m1;
global   symbols64K_m1;
global   vectors256K_m1;
global   symbols256K_m1;
global   vectors1M_m1;
global   symbols1M_m1;


switch qm
   case 2,
      vectors = [0 1];
      symbols = [0 1];
   case 4,
      vectors = [0 1 2 3];
      symbols = [0 1 2 3];
   case 8,
      if mx==1
         vectors = vectors8_m1;
         symbols = symbols8_m1;
      else
         vectors = vectors8_m2;
         symbols = symbols8_m2;
      end
   case 16,
      vectors = vectors16_m1;
      symbols = symbols16_m1;
   case 32,
      vectors = vectors32_m1;
      symbols = symbols32_m1;
   case 64,
      vectors = vectors64_m1;
      symbols = symbols64_m1;
   case 128,
      vectors = vectors128_m1;
      symbols = symbols128_m1;
   case 256,
      vectors = vectors256_m1;
      symbols = symbols256_m1;
   case 512,
      vectors = vectors512_m1;
      symbols = symbols512_m1;
   case 1024,
      vectors = vectors1K_m1;
      symbols = symbols1K_m1;
   case 2048,
      if mx==1
         vectors = vectors2K_m1;
         symbols = symbols2K_m1;
      else




         vectors = vectors2K_m2;
         symbols = symbols2K_m2;
      end
   case 4096,
      vectors = vectors4K_m1;
      symbols = symbols4K_m1;
   case 16384,
      vectors = vectors16K_m1;
      symbols = symbols16K_m1;
   case 32768,
      vectors = vectors32K_m1;
      symbols = symbols32K_m1;
   case 65536,
      vectors = vectors64K_m1;
      symbols = symbols64K_m1;
   case 262144,
      vectors = vectors256K_m1;
      symbols = symbols256K_m1;
   case 1048576,
      vectors = vectors1M_m1;
      symbols = symbols1M_m1;
   otherwise,
      disp 'Invalid Galois Field'
      % IMPORTANT NOTE: FIELDS GF(q^m) FOR WHICH m = 13, 17, AND 19
      % HAVE NO VALID LENGTHS n, AND THEREFORE ARE NOT REQUIRED!
      % (e.g., 8K, 128K, 512K)
end




