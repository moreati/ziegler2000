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

% GF8 symbol to vector mapping and reverse mapping:
% {0 1 alpha alpha^2 ... alpha^6} <---> symbols = {0 1 2 3 ... 7}
%                                   <---> powers = {N/A 0 1 2 ... 6}
% and GF8 is defined by the minimal polynomial selected by mx:
%     mx=1: M(x) = 1 + x + x^3 --> alpha^3 = alpha + 1
%     mx=2: M(x) = 1 + x^2 + x^3 --> alpha^3 = alpha^2 + 1
% 3-bit vectors represent: {alpha^2 alpha 1}
% e.g., alpha^3 = 011 = 3 (mx=1), alpha^3 = 101 = 5 (mx=2)
[vectors8_m1,symbols8_m1] = GenGFx(8,3);
[vectors8_m2,symbols8_m2] = GenGFx(8,5);
[vectors16_m1,symbols16_m1] = GenGFx(16,3);
[vectors32_m1,symbols32_m1] = GenGFx(32,5);
[vectors64_m1,symbols64_m1] = GenGFx(64,3);
[vectors128_m1,symbols128_m1] = GenGFx(128,9);
[vectors256_m1,symbols256_m1] = GenGFx(256,29);
[vectors512_m1,symbols512_m1] = GenGFx(512,17);
[vectors1K_m1,symbols1K_m1] = GenGFx(1024,9);
[vectors2K_m1,symbols2K_m1] = GenGFx(2048,5);
[vectors4K_m1,symbols4K_m1] = GenGFx(4096,83);
[vectors16K_m1,symbols16K_m1] = GenGFx(16384,43);
[vectors32K_m1,symbols32K_m1] = GenGFx(32768,3);
[vectors64K_m1,symbols64K_m1] = GenGFx(65536,45);
[vectors256K_m1,symbols256K_m1] = GenGFx(262144,39);
[vectors1M_m1,symbols1M_m1] = GenGFx(1048576,83);

