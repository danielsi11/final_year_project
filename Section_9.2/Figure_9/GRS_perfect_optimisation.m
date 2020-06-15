function [rate, p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8, p_9, ...
     p_123, p_124, p_125, p_126, p_127, p_128, p_129, p_134, p_135, p_136, p_137, p_138, ...
     p_139, p_145, p_146, p_147, p_148, p_149, p_156, p_157, p_158, p_159, p_167, p_168, ...
     p_169, p_178, p_179, p_189, p_234, p_235, p_236, p_237, p_238, p_239, p_245, p_246, ...
     p_247, p_248, p_249, p_256, p_257, p_258, p_259, p_267, p_268, p_269, p_278, p_279, ...
     p_289, p_345, p_346, p_347, p_348, p_349, p_356, p_357, p_358, p_359, p_367, p_368, ...
     p_369, p_378, p_379, p_389, p_456, p_457, p_458, p_459, p_467, p_468, p_469, p_478, ...
     p_479, p_489, p_567, p_568, p_569, p_578, p_579, p_589, p_678, p_679, p_689, p_789, ...
     p_123456789] = GRS_perfect_optimisation(H,PAC,g,W)
    
[N,K] = size(H); 

D1 = zeros(N,N);
D2 = zeros(N,N);
D3 = zeros(N,N);
D4 = zeros(N,N);
D5 = zeros(N,N);
D6 = zeros(N,N);
D7 = zeros(N,N);
D8 = zeros(N,N);
D9 = zeros(N,N);

D1(1,1) = 1;
D2(2,2) = 1;
D3(3,3) = 1;
D4(4,4) = 1;
D5(5,5) = 1;
D6(6,6) = 1;
D7(7,7) = 1;
D8(8,8) = 1;
D9(9,9) = 1;

g_1_123456789 = g{1};
g_1_123 = g{2};
g_1_124 = g{3};
g_1_125 = g{4};
g_1_126 = g{5};
g_1_127 = g{6};
g_1_128 = g{7};
g_1_129 = g{8};
g_1_134 = g{9};
g_1_135 = g{10};
g_1_136 = g{11};
g_1_137 = g{12};
g_1_138 = g{13};
g_1_139 = g{14};
g_1_145 = g{15};
g_1_146 = g{16};
g_1_147 = g{17};
g_1_148 = g{18};
g_1_149 = g{19};
g_1_156 = g{20};
g_1_157 = g{21};
g_1_158 = g{22};
g_1_159 = g{23};
g_1_167 = g{24};
g_1_168 = g{25};
g_1_169 = g{26};
g_1_178 = g{27};
g_1_179 = g{28};
g_1_189 = g{29};
g_1 = g{30};

g_2_123456789 = g{31};
g_2_123 = g{32};
g_2_124 = g{33};
g_2_125 = g{34};
g_2_126 = g{35};
g_2_127 = g{36};
g_2_128 = g{37};
g_2_129 = g{38};
g_2_134 = g{39};
g_2_135 = g{40};
g_2_136 = g{41};
g_2_137 = g{42};
g_2_138 = g{43};
g_2_139 = g{44};
g_2_145 = g{45};
g_2_146 = g{46};
g_2_147 = g{47};
g_2_148 = g{48};
g_2_149 = g{49};
g_2_156 = g{50};
g_2_157 = g{51};
g_2_158 = g{52};
g_2_159 = g{53};
g_2_167 = g{54};
g_2_168 = g{55};
g_2_169 = g{56};
g_2_178 = g{57};
g_2_179 = g{58};
g_2_189 = g{59};
g_2 = g{60};

g_3_123456789 = g{61};
g_3_123 = g{62};
g_3_124 = g{63};
g_3_125 = g{64};
g_3_126 = g{65};
g_3_127 = g{66};
g_3_128 = g{67};
g_3_129 = g{68};
g_3_234 = g{69};
g_3_235 = g{70};
g_3_236 = g{71};
g_3_237 = g{72};
g_3_238 = g{73};
g_3_239 = g{74};
g_3_245 = g{75};
g_3_246 = g{76};
g_3_247 = g{77};
g_3_248 = g{78};
g_3_249 = g{79};
g_3_256 = g{80};
g_3_257 = g{81};
g_3_258 = g{82};
g_3_259 = g{83};
g_3_267 = g{84};
g_3_268 = g{85};
g_3_269 = g{86};
g_3_278 = g{87};
g_3_279 = g{88};
g_3_289 = g{89};
g_3 = g{90};

g_4_123456789 = g{91};
g_4_123 = g{92};
g_4_124 = g{93};
g_4_125 = g{94};
g_4_126 = g{95};
g_4_127 = g{96};
g_4_128 = g{97};
g_4_129 = g{98};
g_4_234 = g{99};
g_4_235 = g{100};
g_4_236 = g{101};
g_4_237 = g{102};
g_4_238 = g{103};
g_4_239 = g{104};
g_4_245 = g{105};
g_4_246 = g{106};
g_4_247 = g{107};
g_4_248 = g{108};
g_4_249 = g{109};
g_4_256 = g{110};
g_4_257 = g{111};
g_4_258 = g{112};
g_4_259 = g{113};
g_4_267 = g{114};
g_4_268 = g{115};
g_4_269 = g{116};
g_4_278 = g{117};
g_4_279 = g{118};
g_4_289 = g{119};
g_4 = g{120};

g_5_123456789 = g{121};
g_5_123 = g{122};
g_5_134 = g{123};
g_5_135 = g{124};
g_5_136 = g{125};
g_5_137 = g{126};
g_5_138 = g{127};
g_5_139 = g{128};
g_5_234 = g{129};
g_5_235 = g{130};
g_5_236 = g{131};
g_5_237 = g{132};
g_5_238 = g{133};
g_5_239 = g{134};
g_5_345 = g{135};
g_5_346 = g{136};
g_5_347 = g{137};
g_5_348 = g{138};
g_5_349 = g{139};
g_5_356 = g{140};
g_5_357 = g{141};
g_5_358 = g{142};
g_5_359 = g{143};
g_5_367 = g{144};
g_5_368 = g{145};
g_5_369 = g{146};
g_5_378 = g{147};
g_5_379 = g{148};
g_5_389 = g{149};
g_5 = g{150};

g_6_123456789 = g{151};
g_6_123 = g{152};
g_6_134 = g{153};
g_6_135 = g{154};
g_6_136 = g{155};
g_6_137 = g{156};
g_6_138 = g{157};
g_6_139 = g{158};
g_6_234 = g{159};
g_6_235 = g{160};
g_6_236 = g{161};
g_6_237 = g{162};
g_6_238 = g{163};
g_6_239 = g{164};
g_6_345 = g{165};
g_6_346 = g{166};
g_6_347 = g{167};
g_6_348 = g{168};
g_6_349 = g{169};
g_6_356 = g{170};
g_6_357 = g{171};
g_6_358 = g{172};
g_6_359 = g{173};
g_6_367 = g{174};
g_6_368 = g{175};
g_6_369 = g{176};
g_6_378 = g{177};
g_6_379 = g{178};
g_6_389 = g{179};
g_6 = g{180};

g_7_123456789 = g{181};
g_7_124 = g{182};
g_7_134 = g{183};
g_7_145 = g{184};
g_7_146 = g{185};
g_7_147 = g{186};
g_7_148 = g{187};
g_7_149 = g{188};
g_7_234 = g{189};
g_7_245 = g{190};
g_7_246 = g{191};
g_7_247 = g{192};
g_7_248 = g{193};
g_7_249 = g{194};
g_7_345 = g{195};
g_7_346 = g{196};
g_7_347 = g{197};
g_7_348 = g{198};
g_7_349 = g{199};
g_7_456 = g{200};
g_7_457 = g{201};
g_7_458 = g{202};
g_7_459 = g{203};
g_7_467 = g{204};
g_7_468 = g{205};
g_7_469 = g{206};
g_7_478 = g{207};
g_7_479 = g{208};
g_7_489 = g{209};
g_7 = g{210};

g_8_123456789 = g{211};
g_8_124 = g{212};
g_8_134 = g{213};
g_8_145 = g{214};
g_8_146 = g{215};
g_8_147 = g{216};
g_8_148 = g{217};
g_8_149 = g{218};
g_8_234 = g{219};
g_8_245 = g{220};
g_8_246 = g{221};
g_8_247 = g{222};
g_8_248 = g{223};
g_8_249 = g{224};
g_8_345 = g{225};
g_8_346 = g{226};
g_8_347 = g{227};
g_8_348 = g{228};
g_8_349 = g{229};
g_8_456 = g{230};
g_8_457 = g{231};
g_8_458 = g{232};
g_8_459 = g{233};
g_8_467 = g{234};
g_8_468 = g{235};
g_8_469 = g{236};
g_8_478 = g{237};
g_8_479 = g{238};
g_8_489 = g{239};
g_8 = g{240};

g_9_123456789 = g{241};
g_9_125 = g{242};
g_9_135 = g{243};
g_9_145 = g{244};
g_9_156 = g{245};
g_9_157 = g{246};
g_9_158 = g{247};
g_9_159 = g{248};
g_9_235 = g{249};
g_9_245 = g{250};
g_9_256 = g{251};
g_9_257 = g{252};
g_9_258 = g{253};
g_9_259 = g{254};
g_9_345 = g{255};
g_9_356 = g{256};
g_9_357 = g{257};
g_9_358 = g{258};
g_9_359 = g{259};
g_9_456 = g{260};
g_9_457 = g{261};
g_9_458 = g{262};
g_9_459 = g{263};
g_9_567 = g{264};
g_9_568 = g{265};
g_9_569 = g{266};
g_9_578 = g{267};
g_9_579 = g{268};
g_9_589 = g{269};
g_9 = g{270};

g_10_123456789 = g{271};
g_10_125 = g{272};
g_10_135 = g{273};
g_10_145 = g{274};
g_10_156 = g{275};
g_10_157 = g{276};
g_10_158 = g{277};
g_10_159 = g{278};
g_10_235 = g{279};
g_10_245 = g{280};
g_10_256 = g{281};
g_10_257 = g{282};
g_10_258 = g{283};
g_10_259 = g{284};
g_10_345 = g{285};
g_10_356 = g{286};
g_10_357 = g{287};
g_10_358 = g{288};
g_10_359 = g{289};
g_10_456 = g{290};
g_10_457 = g{291};
g_10_458 = g{292};
g_10_459 = g{293};
g_10_567 = g{294};
g_10_568 = g{295};
g_10_569 = g{296};
g_10_578 = g{297};
g_10_579 = g{298};
g_10_589 = g{299};
g_10 = g{300};

g_11_123456789 = g{301};
g_11_126 = g{302};
g_11_136 = g{303};
g_11_146 = g{304};
g_11_156 = g{305};
g_11_167 = g{306};
g_11_168 = g{307};
g_11_169 = g{308};
g_11_236 = g{309};
g_11_246 = g{310};
g_11_256 = g{311};
g_11_267 = g{312};
g_11_268 = g{313};
g_11_269 = g{314};
g_11_346 = g{315};
g_11_356 = g{316};
g_11_367 = g{317};
g_11_368 = g{318};
g_11_369 = g{319};
g_11_456 = g{320};
g_11_467 = g{321};
g_11_468 = g{322};
g_11_469 = g{323};
g_11_567 = g{324};
g_11_568 = g{325};
g_11_569 = g{326};
g_11_678 = g{327};
g_11_679 = g{328};
g_11_689 = g{329};
g_11 = g{330};

g_12_123456789 = g{331};
g_12_126 = g{332};
g_12_136 = g{333};
g_12_146 = g{334};
g_12_156 = g{335};
g_12_167 = g{336};
g_12_168 = g{337};
g_12_169 = g{338};
g_12_236 = g{339};
g_12_246 = g{340};
g_12_256 = g{341};
g_12_267 = g{342};
g_12_268 = g{343};
g_12_269 = g{344};
g_12_346 = g{345};
g_12_356 = g{346};
g_12_367 = g{347};
g_12_368 = g{348};
g_12_369 = g{349};
g_12_456 = g{350};
g_12_467 = g{351};
g_12_468 = g{352};
g_12_469 = g{353};
g_12_567 = g{354};
g_12_568 = g{355};
g_12_569 = g{356};
g_12_678 = g{357};
g_12_679 = g{358};
g_12_689 = g{359};
g_12 = g{360};

g_13_123456789 = g{361};
g_13_127 = g{362};
g_13_137 = g{363};
g_13_147 = g{364};
g_13_157 = g{365};
g_13_167 = g{366};
g_13_178 = g{367};
g_13_179 = g{368};
g_13_237 = g{369};
g_13_247 = g{370};
g_13_257 = g{371};
g_13_267 = g{372};
g_13_278 = g{373};
g_13_279 = g{374};
g_13_347 = g{375};
g_13_357 = g{376};
g_13_367 = g{377};
g_13_378 = g{378};
g_13_379 = g{379};
g_13_457 = g{380};
g_13_467 = g{381};
g_13_478 = g{382};
g_13_479 = g{383};
g_13_567 = g{384};
g_13_578 = g{385};
g_13_579 = g{386};
g_13_678 = g{387};
g_13_679 = g{388};
g_13_789 = g{389};
g_13 = g{390};

g_14_123456789 = g{391};
g_14_127 = g{392};
g_14_137 = g{393};
g_14_147 = g{394};
g_14_157 = g{395};
g_14_167 = g{396};
g_14_178 = g{397};
g_14_179 = g{398};
g_14_237 = g{399};
g_14_247 = g{400};
g_14_257 = g{401};
g_14_267 = g{402};
g_14_278 = g{403};
g_14_279 = g{404};
g_14_347 = g{405};
g_14_357 = g{406};
g_14_367 = g{407};
g_14_378 = g{408};
g_14_379 = g{409};
g_14_457 = g{410};
g_14_467 = g{411};
g_14_478 = g{412};
g_14_479 = g{413};
g_14_567 = g{414};
g_14_578 = g{415};
g_14_579 = g{416};
g_14_678 = g{417};
g_14_679 = g{418};
g_14_789 = g{419};
g_14 = g{420};

g_15_123456789 = g{421};
g_15_128 = g{422};
g_15_138 = g{423};
g_15_148 = g{424};
g_15_158 = g{425};
g_15_168 = g{426};
g_15_178 = g{427};
g_15_189 = g{428};
g_15_238 = g{429};
g_15_248 = g{430};
g_15_258 = g{431};
g_15_268 = g{432};
g_15_278 = g{433};
g_15_289 = g{434};
g_15_348 = g{435};
g_15_358 = g{436};
g_15_368 = g{437};
g_15_378 = g{438};
g_15_389 = g{439};
g_15_458 = g{440};
g_15_468 = g{441};
g_15_478 = g{442};
g_15_489 = g{443};
g_15_568 = g{444};
g_15_578 = g{445};
g_15_589 = g{446};
g_15_678 = g{447};
g_15_689 = g{448};
g_15_789 = g{449};
g_15 = g{450};

g_16_123456789 = g{451};
g_16_128 = g{452};
g_16_138 = g{453};
g_16_148 = g{454};
g_16_158 = g{455};
g_16_168 = g{456};
g_16_178 = g{457};
g_16_189 = g{458};
g_16_238 = g{459};
g_16_248 = g{460};
g_16_258 = g{461};
g_16_268 = g{462};
g_16_278 = g{463};
g_16_289 = g{464};
g_16_348 = g{465};
g_16_358 = g{466};
g_16_368 = g{467};
g_16_378 = g{468};
g_16_389 = g{469};
g_16_458 = g{470};
g_16_468 = g{471};
g_16_478 = g{472};
g_16_489 = g{473};
g_16_568 = g{474};
g_16_578 = g{475};
g_16_589 = g{476};
g_16_678 = g{477};
g_16_689 = g{478};
g_16_789 = g{479};
g_16 = g{480};

g_17_123456789 = g{481};
g_17_129 = g{482};
g_17_139 = g{483};
g_17_149 = g{484};
g_17_159 = g{485};
g_17_169 = g{486};
g_17_179 = g{487};
g_17_189 = g{488};
g_17_239 = g{489};
g_17_249 = g{490};
g_17_259 = g{491};
g_17_269 = g{492};
g_17_279 = g{493};
g_17_289 = g{494};
g_17_349 = g{495};
g_17_359 = g{496};
g_17_369 = g{497};
g_17_379 = g{498};
g_17_389 = g{499};
g_17_459 = g{500};
g_17_469 = g{501};
g_17_479 = g{502};
g_17_489 = g{503};
g_17_569 = g{504};
g_17_579 = g{505};
g_17_589 = g{506};
g_17_679 = g{507};
g_17_689 = g{508};
g_17_789 = g{509};
g_17 = g{510};

g_18_123456789 = g{511};
g_18_129 = g{512};
g_18_139 = g{513};
g_18_149 = g{514};
g_18_159 = g{515};
g_18_169 = g{516};
g_18_179 = g{517};
g_18_189 = g{518};
g_18_239 = g{519};
g_18_249 = g{520};
g_18_259 = g{521};
g_18_269 = g{522};
g_18_279 = g{523};
g_18_289 = g{524};
g_18_349 = g{525};
g_18_359 = g{526};
g_18_369 = g{527};
g_18_379 = g{528};
g_18_389 = g{529};
g_18_459 = g{530};
g_18_469 = g{531};
g_18_479 = g{532};
g_18_489 = g{533};
g_18_569 = g{534};
g_18_579 = g{535};
g_18_589 = g{536};
g_18_679 = g{537};
g_18_689 = g{538};
g_18_789 = g{539};
g_18 = g{540};

W_1_123456789 = W{1};
W_1_123 = W{2};
W_1_124 = W{3};
W_1_125 = W{4};
W_1_126 = W{5};
W_1_127 = W{6};
W_1_128 = W{7};
W_1_129 = W{8};
W_1_134 = W{9};
W_1_135 = W{10};
W_1_136 = W{11};
W_1_137 = W{12};
W_1_138 = W{13};
W_1_139 = W{14};
W_1_145 = W{15};
W_1_146 = W{16};
W_1_147 = W{17};
W_1_148 = W{18};
W_1_149 = W{19};
W_1_156 = W{20};
W_1_157 = W{21};
W_1_158 = W{22};
W_1_159 = W{23};
W_1_167 = W{24};
W_1_168 = W{25};
W_1_169 = W{26};
W_1_178 = W{27};
W_1_179 = W{28};
W_1_189 = W{29};
W_1 = W{30};

W_2_123456789 = W{31};
W_2_123 = W{32};
W_2_124 = W{33};
W_2_125 = W{34};
W_2_126 = W{35};
W_2_127 = W{36};
W_2_128 = W{37};
W_2_129 = W{38};
W_2_134 = W{39};
W_2_135 = W{40};
W_2_136 = W{41};
W_2_137 = W{42};
W_2_138 = W{43};
W_2_139 = W{44};
W_2_145 = W{45};
W_2_146 = W{46};
W_2_147 = W{47};
W_2_148 = W{48};
W_2_149 = W{49};
W_2_156 = W{50};
W_2_157 = W{51};
W_2_158 = W{52};
W_2_159 = W{53};
W_2_167 = W{54};
W_2_168 = W{55};
W_2_169 = W{56};
W_2_178 = W{57};
W_2_179 = W{58};
W_2_189 = W{59};
W_2 = W{60};

W_3_123456789 = W{61};
W_3_123 = W{62};
W_3_124 = W{63};
W_3_125 = W{64};
W_3_126 = W{65};
W_3_127 = W{66};
W_3_128 = W{67};
W_3_129 = W{68};
W_3_234 = W{69};
W_3_235 = W{70};
W_3_236 = W{71};
W_3_237 = W{72};
W_3_238 = W{73};
W_3_239 = W{74};
W_3_245 = W{75};
W_3_246 = W{76};
W_3_247 = W{77};
W_3_248 = W{78};
W_3_249 = W{79};
W_3_256 = W{80};
W_3_257 = W{81};
W_3_258 = W{82};
W_3_259 = W{83};
W_3_267 = W{84};
W_3_268 = W{85};
W_3_269 = W{86};
W_3_278 = W{87};
W_3_279 = W{88};
W_3_289 = W{89};
W_3 = W{90};

W_4_123456789 = W{91};
W_4_123 = W{92};
W_4_124 = W{93};
W_4_125 = W{94};
W_4_126 = W{95};
W_4_127 = W{96};
W_4_128 = W{97};
W_4_129 = W{98};
W_4_234 = W{99};
W_4_235 = W{100};
W_4_236 = W{101};
W_4_237 = W{102};
W_4_238 = W{103};
W_4_239 = W{104};
W_4_245 = W{105};
W_4_246 = W{106};
W_4_247 = W{107};
W_4_248 = W{108};
W_4_249 = W{109};
W_4_256 = W{110};
W_4_257 = W{111};
W_4_258 = W{112};
W_4_259 = W{113};
W_4_267 = W{114};
W_4_268 = W{115};
W_4_269 = W{116};
W_4_278 = W{117};
W_4_279 = W{118};
W_4_289 = W{119};
W_4 = W{120};

W_5_123456789 = W{121};
W_5_123 = W{122};
W_5_134 = W{123};
W_5_135 = W{124};
W_5_136 = W{125};
W_5_137 = W{126};
W_5_138 = W{127};
W_5_139 = W{128};
W_5_234 = W{129};
W_5_235 = W{130};
W_5_236 = W{131};
W_5_237 = W{132};
W_5_238 = W{133};
W_5_239 = W{134};
W_5_345 = W{135};
W_5_346 = W{136};
W_5_347 = W{137};
W_5_348 = W{138};
W_5_349 = W{139};
W_5_356 = W{140};
W_5_357 = W{141};
W_5_358 = W{142};
W_5_359 = W{143};
W_5_367 = W{144};
W_5_368 = W{145};
W_5_369 = W{146};
W_5_378 = W{147};
W_5_379 = W{148};
W_5_389 = W{149};
W_5 = W{150};

W_6_123456789 = W{151};
W_6_123 = W{152};
W_6_134 = W{153};
W_6_135 = W{154};
W_6_136 = W{155};
W_6_137 = W{156};
W_6_138 = W{157};
W_6_139 = W{158};
W_6_234 = W{159};
W_6_235 = W{160};
W_6_236 = W{161};
W_6_237 = W{162};
W_6_238 = W{163};
W_6_239 = W{164};
W_6_345 = W{165};
W_6_346 = W{166};
W_6_347 = W{167};
W_6_348 = W{168};
W_6_349 = W{169};
W_6_356 = W{170};
W_6_357 = W{171};
W_6_358 = W{172};
W_6_359 = W{173};
W_6_367 = W{174};
W_6_368 = W{175};
W_6_369 = W{176};
W_6_378 = W{177};
W_6_379 = W{178};
W_6_389 = W{179};
W_6 = W{180};

W_7_123456789 = W{181};
W_7_124 = W{182};
W_7_134 = W{183};
W_7_145 = W{184};
W_7_146 = W{185};
W_7_147 = W{186};
W_7_148 = W{187};
W_7_149 = W{188};
W_7_234 = W{189};
W_7_245 = W{190};
W_7_246 = W{191};
W_7_247 = W{192};
W_7_248 = W{193};
W_7_249 = W{194};
W_7_345 = W{195};
W_7_346 = W{196};
W_7_347 = W{197};
W_7_348 = W{198};
W_7_349 = W{199};
W_7_456 = W{200};
W_7_457 = W{201};
W_7_458 = W{202};
W_7_459 = W{203};
W_7_467 = W{204};
W_7_468 = W{205};
W_7_469 = W{206};
W_7_478 = W{207};
W_7_479 = W{208};
W_7_489 = W{209};
W_7 = W{210};

W_8_123456789 = W{211};
W_8_124 = W{212};
W_8_134 = W{213};
W_8_145 = W{214};
W_8_146 = W{215};
W_8_147 = W{216};
W_8_148 = W{217};
W_8_149 = W{218};
W_8_234 = W{219};
W_8_245 = W{220};
W_8_246 = W{221};
W_8_247 = W{222};
W_8_248 = W{223};
W_8_249 = W{224};
W_8_345 = W{225};
W_8_346 = W{226};
W_8_347 = W{227};
W_8_348 = W{228};
W_8_349 = W{229};
W_8_456 = W{230};
W_8_457 = W{231};
W_8_458 = W{232};
W_8_459 = W{233};
W_8_467 = W{234};
W_8_468 = W{235};
W_8_469 = W{236};
W_8_478 = W{237};
W_8_479 = W{238};
W_8_489 = W{239};
W_8 = W{240};

W_9_123456789 = W{241};
W_9_125 = W{242};
W_9_135 = W{243};
W_9_145 = W{244};
W_9_156 = W{245};
W_9_157 = W{246};
W_9_158 = W{247};
W_9_159 = W{248};
W_9_235 = W{249};
W_9_245 = W{250};
W_9_256 = W{251};
W_9_257 = W{252};
W_9_258 = W{253};
W_9_259 = W{254};
W_9_345 = W{255};
W_9_356 = W{256};
W_9_357 = W{257};
W_9_358 = W{258};
W_9_359 = W{259};
W_9_456 = W{260};
W_9_457 = W{261};
W_9_458 = W{262};
W_9_459 = W{263};
W_9_567 = W{264};
W_9_568 = W{265};
W_9_569 = W{266};
W_9_578 = W{267};
W_9_579 = W{268};
W_9_589 = W{269};
W_9 = W{270};

W_10_123456789 = W{271};
W_10_125 = W{272};
W_10_135 = W{273};
W_10_145 = W{274};
W_10_156 = W{275};
W_10_157 = W{276};
W_10_158 = W{277};
W_10_159 = W{278};
W_10_235 = W{279};
W_10_245 = W{280};
W_10_256 = W{281};
W_10_257 = W{282};
W_10_258 = W{283};
W_10_259 = W{284};
W_10_345 = W{285};
W_10_356 = W{286};
W_10_357 = W{287};
W_10_358 = W{288};
W_10_359 = W{289};
W_10_456 = W{290};
W_10_457 = W{291};
W_10_458 = W{292};
W_10_459 = W{293};
W_10_567 = W{294};
W_10_568 = W{295};
W_10_569 = W{296};
W_10_578 = W{297};
W_10_579 = W{298};
W_10_589 = W{299};
W_10 = W{300};

W_11_123456789 = W{301};
W_11_126 = W{302};
W_11_136 = W{303};
W_11_146 = W{304};
W_11_156 = W{305};
W_11_167 = W{306};
W_11_168 = W{307};
W_11_169 = W{308};
W_11_236 = W{309};
W_11_246 = W{310};
W_11_256 = W{311};
W_11_267 = W{312};
W_11_268 = W{313};
W_11_269 = W{314};
W_11_346 = W{315};
W_11_356 = W{316};
W_11_367 = W{317};
W_11_368 = W{318};
W_11_369 = W{319};
W_11_456 = W{320};
W_11_467 = W{321};
W_11_468 = W{322};
W_11_469 = W{323};
W_11_567 = W{324};
W_11_568 = W{325};
W_11_569 = W{326};
W_11_678 = W{327};
W_11_679 = W{328};
W_11_689 = W{329};
W_11 = W{330};

W_12_123456789 = W{331};
W_12_126 = W{332};
W_12_136 = W{333};
W_12_146 = W{334};
W_12_156 = W{335};
W_12_167 = W{336};
W_12_168 = W{337};
W_12_169 = W{338};
W_12_236 = W{339};
W_12_246 = W{340};
W_12_256 = W{341};
W_12_267 = W{342};
W_12_268 = W{343};
W_12_269 = W{344};
W_12_346 = W{345};
W_12_356 = W{346};
W_12_367 = W{347};
W_12_368 = W{348};
W_12_369 = W{349};
W_12_456 = W{350};
W_12_467 = W{351};
W_12_468 = W{352};
W_12_469 = W{353};
W_12_567 = W{354};
W_12_568 = W{355};
W_12_569 = W{356};
W_12_678 = W{357};
W_12_679 = W{358};
W_12_689 = W{359};
W_12 = W{360};

W_13_123456789 = W{361};
W_13_127 = W{362};
W_13_137 = W{363};
W_13_147 = W{364};
W_13_157 = W{365};
W_13_167 = W{366};
W_13_178 = W{367};
W_13_179 = W{368};
W_13_237 = W{369};
W_13_247 = W{370};
W_13_257 = W{371};
W_13_267 = W{372};
W_13_278 = W{373};
W_13_279 = W{374};
W_13_347 = W{375};
W_13_357 = W{376};
W_13_367 = W{377};
W_13_378 = W{378};
W_13_379 = W{379};
W_13_457 = W{380};
W_13_467 = W{381};
W_13_478 = W{382};
W_13_479 = W{383};
W_13_567 = W{384};
W_13_578 = W{385};
W_13_579 = W{386};
W_13_678 = W{387};
W_13_679 = W{388};
W_13_789 = W{389};
W_13 = W{390};

W_14_123456789 = W{391};
W_14_127 = W{392};
W_14_137 = W{393};
W_14_147 = W{394};
W_14_157 = W{395};
W_14_167 = W{396};
W_14_178 = W{397};
W_14_179 = W{398};
W_14_237 = W{399};
W_14_247 = W{400};
W_14_257 = W{401};
W_14_267 = W{402};
W_14_278 = W{403};
W_14_279 = W{404};
W_14_347 = W{405};
W_14_357 = W{406};
W_14_367 = W{407};
W_14_378 = W{408};
W_14_379 = W{409};
W_14_457 = W{410};
W_14_467 = W{411};
W_14_478 = W{412};
W_14_479 = W{413};
W_14_567 = W{414};
W_14_578 = W{415};
W_14_579 = W{416};
W_14_678 = W{417};
W_14_679 = W{418};
W_14_789 = W{419};
W_14 = W{420};

W_15_123456789 = W{421};
W_15_128 = W{422};
W_15_138 = W{423};
W_15_148 = W{424};
W_15_158 = W{425};
W_15_168 = W{426};
W_15_178 = W{427};
W_15_189 = W{428};
W_15_238 = W{429};
W_15_248 = W{430};
W_15_258 = W{431};
W_15_268 = W{432};
W_15_278 = W{433};
W_15_289 = W{434};
W_15_348 = W{435};
W_15_358 = W{436};
W_15_368 = W{437};
W_15_378 = W{438};
W_15_389 = W{439};
W_15_458 = W{440};
W_15_468 = W{441};
W_15_478 = W{442};
W_15_489 = W{443};
W_15_568 = W{444};
W_15_578 = W{445};
W_15_589 = W{446};
W_15_678 = W{447};
W_15_689 = W{448};
W_15_789 = W{449};
W_15 = W{450};

W_16_123456789 = W{451};
W_16_128 = W{452};
W_16_138 = W{453};
W_16_148 = W{454};
W_16_158 = W{455};
W_16_168 = W{456};
W_16_178 = W{457};
W_16_189 = W{458};
W_16_238 = W{459};
W_16_248 = W{460};
W_16_258 = W{461};
W_16_268 = W{462};
W_16_278 = W{463};
W_16_289 = W{464};
W_16_348 = W{465};
W_16_358 = W{466};
W_16_368 = W{467};
W_16_378 = W{468};
W_16_389 = W{469};
W_16_458 = W{470};
W_16_468 = W{471};
W_16_478 = W{472};
W_16_489 = W{473};
W_16_568 = W{474};
W_16_578 = W{475};
W_16_589 = W{476};
W_16_678 = W{477};
W_16_689 = W{478};
W_16_789 = W{479};
W_16 = W{480};

W_17_123456789 = W{481};
W_17_129 = W{482};
W_17_139 = W{483};
W_17_149 = W{484};
W_17_159 = W{485};
W_17_169 = W{486};
W_17_179 = W{487};
W_17_189 = W{488};
W_17_239 = W{489};
W_17_249 = W{490};
W_17_259 = W{491};
W_17_269 = W{492};
W_17_279 = W{493};
W_17_289 = W{494};
W_17_349 = W{495};
W_17_359 = W{496};
W_17_369 = W{497};
W_17_379 = W{498};
W_17_389 = W{499};
W_17_459 = W{500};
W_17_469 = W{501};
W_17_479 = W{502};
W_17_489 = W{503};
W_17_569 = W{504};
W_17_579 = W{505};
W_17_589 = W{506};
W_17_679 = W{507};
W_17_689 = W{508};
W_17_789 = W{509};
W_17 = W{510};

W_18_123456789 = W{511};
W_18_129 = W{512};
W_18_139 = W{513};
W_18_149 = W{514};
W_18_159 = W{515};
W_18_169 = W{516};
W_18_179 = W{517};
W_18_189 = W{518};
W_18_239 = W{519};
W_18_249 = W{520};
W_18_259 = W{521};
W_18_269 = W{522};
W_18_279 = W{523};
W_18_289 = W{524};
W_18_349 = W{525};
W_18_359 = W{526};
W_18_369 = W{527};
W_18_379 = W{528};
W_18_389 = W{529};
W_18_459 = W{530};
W_18_469 = W{531};
W_18_479 = W{532};
W_18_489 = W{533};
W_18_569 = W{534};
W_18_579 = W{535};
W_18_589 = W{536};
W_18_679 = W{537};
W_18_689 = W{538};
W_18_789 = W{539};
W_18 = W{540};


%CVX Optimisation Tool
cvx_begin quiet

variable p_1(N,1) complex
variable p_2(N,1) complex
variable p_3(N,1) complex
variable p_4(N,1) complex
variable p_5(N,1) complex
variable p_6(N,1) complex
variable p_7(N,1) complex
variable p_8(N,1) complex
variable p_9(N,1) complex

variable p_123(N,1) complex
variable p_124(N,1) complex
variable p_125(N,1) complex
variable p_126(N,1) complex
variable p_127(N,1) complex
variable p_128(N,1) complex
variable p_129(N,1) complex
variable p_134(N,1) complex
variable p_135(N,1) complex
variable p_136(N,1) complex
variable p_137(N,1) complex
variable p_138(N,1) complex
variable p_139(N,1) complex
variable p_145(N,1) complex
variable p_146(N,1) complex
variable p_147(N,1) complex
variable p_148(N,1) complex
variable p_149(N,1) complex
variable p_156(N,1) complex
variable p_157(N,1) complex
variable p_158(N,1) complex
variable p_159(N,1) complex
variable p_167(N,1) complex
variable p_168(N,1) complex
variable p_169(N,1) complex
variable p_178(N,1) complex
variable p_179(N,1) complex
variable p_189(N,1) complex
variable p_234(N,1) complex
variable p_235(N,1) complex
variable p_236(N,1) complex
variable p_237(N,1) complex
variable p_238(N,1) complex
variable p_239(N,1) complex
variable p_245(N,1) complex
variable p_246(N,1) complex
variable p_247(N,1) complex
variable p_248(N,1) complex
variable p_249(N,1) complex
variable p_256(N,1) complex
variable p_257(N,1) complex
variable p_258(N,1) complex
variable p_259(N,1) complex
variable p_267(N,1) complex
variable p_268(N,1) complex
variable p_269(N,1) complex
variable p_278(N,1) complex
variable p_279(N,1) complex
variable p_289(N,1) complex
variable p_345(N,1) complex
variable p_346(N,1) complex
variable p_347(N,1) complex
variable p_348(N,1) complex
variable p_349(N,1) complex
variable p_356(N,1) complex
variable p_357(N,1) complex
variable p_358(N,1) complex
variable p_359(N,1) complex
variable p_367(N,1) complex
variable p_368(N,1) complex
variable p_369(N,1) complex
variable p_378(N,1) complex
variable p_379(N,1) complex
variable p_389(N,1) complex
variable p_456(N,1) complex
variable p_457(N,1) complex
variable p_458(N,1) complex
variable p_459(N,1) complex
variable p_467(N,1) complex
variable p_468(N,1) complex
variable p_469(N,1) complex
variable p_478(N,1) complex
variable p_479(N,1) complex
variable p_489(N,1) complex
variable p_567(N,1) complex
variable p_568(N,1) complex
variable p_569(N,1) complex
variable p_578(N,1) complex
variable p_579(N,1) complex
variable p_589(N,1) complex
variable p_678(N,1) complex
variable p_679(N,1) complex
variable p_689(N,1) complex
variable p_789(N,1) complex

variable p_123456789(N,1) complex

variable r_1 
variable r_2 
variable r_3
variable r_4 
variable r_5 
variable r_6
variable r_7
variable r_8
variable r_9

variable r_1_123
variable r_2_123
variable r_3_123
variable r_1_124
variable r_2_124
variable r_4_124
variable r_1_125
variable r_2_125
variable r_5_125
variable r_1_126
variable r_2_126
variable r_6_126
variable r_1_127
variable r_2_127
variable r_7_127
variable r_1_128
variable r_2_128
variable r_8_128
variable r_1_129
variable r_2_129
variable r_9_129
variable r_1_134
variable r_3_134
variable r_4_134
variable r_1_135
variable r_3_135
variable r_5_135
variable r_1_136
variable r_3_136
variable r_6_136
variable r_1_137
variable r_3_137
variable r_7_137
variable r_1_138
variable r_3_138
variable r_8_138
variable r_1_139
variable r_3_139
variable r_9_139
variable r_1_145
variable r_4_145
variable r_5_145
variable r_1_146
variable r_4_146
variable r_6_146
variable r_1_147
variable r_4_147
variable r_7_147
variable r_1_148
variable r_4_148
variable r_8_148
variable r_1_149
variable r_4_149
variable r_9_149
variable r_1_156
variable r_5_156
variable r_6_156
variable r_1_157
variable r_5_157
variable r_7_157
variable r_1_158
variable r_5_158
variable r_8_158
variable r_1_159
variable r_5_159
variable r_9_159
variable r_1_167
variable r_6_167
variable r_7_167
variable r_1_168
variable r_6_168
variable r_8_168
variable r_1_169
variable r_6_169
variable r_9_169
variable r_1_178
variable r_7_178
variable r_8_178
variable r_1_179
variable r_7_179
variable r_9_179
variable r_1_189
variable r_8_189
variable r_9_189
variable r_2_234
variable r_3_234
variable r_4_234
variable r_2_235
variable r_3_235
variable r_5_235
variable r_2_236
variable r_3_236
variable r_6_236
variable r_2_237
variable r_3_237
variable r_7_237
variable r_2_238
variable r_3_238
variable r_8_238
variable r_2_239
variable r_3_239
variable r_9_239
variable r_2_245
variable r_4_245
variable r_5_245
variable r_2_246
variable r_4_246
variable r_6_246
variable r_2_247
variable r_4_247
variable r_7_247
variable r_2_248
variable r_4_248
variable r_8_248
variable r_2_249
variable r_4_249
variable r_9_249
variable r_2_256
variable r_5_256
variable r_6_256
variable r_2_257
variable r_5_257
variable r_7_257
variable r_2_258
variable r_5_258
variable r_8_258
variable r_2_259
variable r_5_259
variable r_9_259
variable r_2_267
variable r_6_267
variable r_7_267
variable r_2_268
variable r_6_268
variable r_8_268
variable r_2_269
variable r_6_269
variable r_9_269
variable r_2_278
variable r_7_278
variable r_8_278
variable r_2_279
variable r_7_279
variable r_9_279
variable r_2_289
variable r_8_289
variable r_9_289
variable r_3_345
variable r_4_345
variable r_5_345
variable r_3_346
variable r_4_346
variable r_6_346
variable r_3_347
variable r_4_347
variable r_7_347
variable r_3_348
variable r_4_348
variable r_8_348
variable r_3_349
variable r_4_349
variable r_9_349
variable r_3_356
variable r_5_356
variable r_6_356
variable r_3_357
variable r_5_357
variable r_7_357
variable r_3_358
variable r_5_358
variable r_8_358
variable r_3_359
variable r_5_359
variable r_9_359
variable r_3_367
variable r_6_367
variable r_7_367
variable r_3_368
variable r_6_368
variable r_8_368
variable r_3_369
variable r_6_369
variable r_9_369
variable r_3_378
variable r_7_378
variable r_8_378
variable r_3_379
variable r_7_379
variable r_9_379
variable r_3_389
variable r_8_389
variable r_9_389
variable r_4_456
variable r_5_456
variable r_6_456
variable r_4_457
variable r_5_457
variable r_7_457
variable r_4_458
variable r_5_458
variable r_8_458
variable r_4_459
variable r_5_459
variable r_9_459
variable r_4_467
variable r_6_467
variable r_7_467
variable r_4_468
variable r_6_468
variable r_8_468
variable r_4_469
variable r_6_469
variable r_9_469
variable r_4_478
variable r_7_478
variable r_8_478
variable r_4_479
variable r_7_479
variable r_9_479
variable r_4_489
variable r_8_489
variable r_9_489
variable r_5_567
variable r_6_567
variable r_7_567
variable r_5_568
variable r_6_568
variable r_8_568
variable r_5_569
variable r_6_569
variable r_9_569
variable r_5_578
variable r_7_578
variable r_8_578
variable r_5_579
variable r_7_579
variable r_9_579
variable r_5_589
variable r_8_589
variable r_9_589
variable r_6_678
variable r_7_678
variable r_8_678
variable r_6_679
variable r_7_679
variable r_9_679
variable r_6_689
variable r_8_689
variable r_9_689
variable r_7_789
variable r_8_789
variable r_9_789

variable r_1_123456789
variable r_2_123456789
variable r_3_123456789
variable r_4_123456789
variable r_5_123456789
variable r_6_123456789
variable r_7_123456789
variable r_8_123456789
variable r_9_123456789

variable r_g

expression constraints(1,819);


%Recieved Power
T_1 = square_abs(H(:,1)'*p_234) + square_abs(H(:,1)'*p_235) + square_abs(H(:,1)'*p_236) + square_abs(H(:,1)'*p_237) + ...
      square_abs(H(:,1)'*p_238) + square_abs(H(:,1)'*p_239) + square_abs(H(:,1)'*p_245) + square_abs(H(:,1)'*p_246) + ...
      square_abs(H(:,1)'*p_247) + square_abs(H(:,1)'*p_248) + square_abs(H(:,1)'*p_249) + square_abs(H(:,1)'*p_256) + ...
      square_abs(H(:,1)'*p_257) + square_abs(H(:,1)'*p_258) + square_abs(H(:,1)'*p_259) + square_abs(H(:,1)'*p_267) + ...
      square_abs(H(:,1)'*p_268) + square_abs(H(:,1)'*p_269) + square_abs(H(:,1)'*p_278) + square_abs(H(:,1)'*p_279) + ...
      square_abs(H(:,1)'*p_289) + square_abs(H(:,1)'*p_345) + square_abs(H(:,1)'*p_346) + square_abs(H(:,1)'*p_347) + ...
      square_abs(H(:,1)'*p_348) + square_abs(H(:,1)'*p_349) + square_abs(H(:,1)'*p_356) + square_abs(H(:,1)'*p_357) + ...
      square_abs(H(:,1)'*p_358) + square_abs(H(:,1)'*p_359) + square_abs(H(:,1)'*p_367) + square_abs(H(:,1)'*p_368) + ...
      square_abs(H(:,1)'*p_369) + square_abs(H(:,1)'*p_378) + square_abs(H(:,1)'*p_379) + square_abs(H(:,1)'*p_389) + ...
      square_abs(H(:,1)'*p_456) + square_abs(H(:,1)'*p_457) + square_abs(H(:,1)'*p_458) + square_abs(H(:,1)'*p_459) + ...
      square_abs(H(:,1)'*p_467) + square_abs(H(:,1)'*p_468) + square_abs(H(:,1)'*p_469) + square_abs(H(:,1)'*p_478) + ...
      square_abs(H(:,1)'*p_479) + square_abs(H(:,1)'*p_489) + square_abs(H(:,1)'*p_567) + square_abs(H(:,1)'*p_568) + ...
      square_abs(H(:,1)'*p_569) + square_abs(H(:,1)'*p_578) + square_abs(H(:,1)'*p_579) + square_abs(H(:,1)'*p_589) + ...
      square_abs(H(:,1)'*p_678) + square_abs(H(:,1)'*p_679) + square_abs(H(:,1)'*p_689) + square_abs(H(:,1)'*p_789) + ...
      square_abs(H(:,1)'*p_1) + square_abs(H(:,1)'*p_2) + square_abs(H(:,1)'*p_3) + ...
      square_abs(H(:,1)'*p_4) + square_abs(H(:,1)'*p_5) + square_abs(H(:,1)'*p_6) + ...
      square_abs(H(:,1)'*p_7) + square_abs(H(:,1)'*p_8) + square_abs(H(:,1)'*p_9) + 1;

T_1_189 = square_abs(H(:,1)'*p_189) + T_1;
T_1_179 = square_abs(H(:,1)'*p_179) + T_1_189;
T_1_178 = square_abs(H(:,1)'*p_178) + T_1_179;
T_1_169 = square_abs(H(:,1)'*p_169) + T_1_178;
T_1_168 = square_abs(H(:,1)'*p_168) + T_1_169;
T_1_167 = square_abs(H(:,1)'*p_167) + T_1_168;
T_1_159 = square_abs(H(:,1)'*p_159) + T_1_167;
T_1_158 = square_abs(H(:,1)'*p_158) + T_1_159;
T_1_157 = square_abs(H(:,1)'*p_157) + T_1_158;
T_1_156 = square_abs(H(:,1)'*p_156) + T_1_157;
T_1_149 = square_abs(H(:,1)'*p_149) + T_1_156;
T_1_148 = square_abs(H(:,1)'*p_148) + T_1_149;
T_1_147 = square_abs(H(:,1)'*p_147) + T_1_148;
T_1_146 = square_abs(H(:,1)'*p_146) + T_1_147;
T_1_145 = square_abs(H(:,1)'*p_145) + T_1_146;
T_1_139 = square_abs(H(:,1)'*p_139) + T_1_145;
T_1_138 = square_abs(H(:,1)'*p_138) + T_1_139;
T_1_137 = square_abs(H(:,1)'*p_137) + T_1_138;
T_1_136 = square_abs(H(:,1)'*p_136) + T_1_137;
T_1_135 = square_abs(H(:,1)'*p_135) + T_1_136;
T_1_134 = square_abs(H(:,1)'*p_134) + T_1_135;
T_1_129 = square_abs(H(:,1)'*p_129) + T_1_134;
T_1_128 = square_abs(H(:,1)'*p_128) + T_1_129;
T_1_127 = square_abs(H(:,1)'*p_127) + T_1_128;
T_1_126 = square_abs(H(:,1)'*p_126) + T_1_127;
T_1_125 = square_abs(H(:,1)'*p_125) + T_1_126;
T_1_124 = square_abs(H(:,1)'*p_124) + T_1_125;
T_1_123 = square_abs(H(:,1)'*p_123) + T_1_124;
T_1_123456789 = square_abs(H(:,1)'*p_123456789) + T_1_123;

T_2 = square_abs(H(:,2)'*p_234) + square_abs(H(:,2)'*p_235) + square_abs(H(:,2)'*p_236) + square_abs(H(:,2)'*p_237) + ...
      square_abs(H(:,2)'*p_238) + square_abs(H(:,2)'*p_239) + square_abs(H(:,2)'*p_245) + square_abs(H(:,2)'*p_246) + ...
      square_abs(H(:,2)'*p_247) + square_abs(H(:,2)'*p_248) + square_abs(H(:,2)'*p_249) + square_abs(H(:,2)'*p_256) + ...
      square_abs(H(:,2)'*p_257) + square_abs(H(:,2)'*p_258) + square_abs(H(:,2)'*p_259) + square_abs(H(:,2)'*p_267) + ...
      square_abs(H(:,2)'*p_268) + square_abs(H(:,2)'*p_269) + square_abs(H(:,2)'*p_278) + square_abs(H(:,2)'*p_279) + ...
      square_abs(H(:,2)'*p_289) + square_abs(H(:,2)'*p_345) + square_abs(H(:,2)'*p_346) + square_abs(H(:,2)'*p_347) + ...
      square_abs(H(:,2)'*p_348) + square_abs(H(:,2)'*p_349) + square_abs(H(:,2)'*p_356) + square_abs(H(:,2)'*p_357) + ...
      square_abs(H(:,2)'*p_358) + square_abs(H(:,2)'*p_359) + square_abs(H(:,2)'*p_367) + square_abs(H(:,2)'*p_368) + ...
      square_abs(H(:,2)'*p_369) + square_abs(H(:,2)'*p_378) + square_abs(H(:,2)'*p_379) + square_abs(H(:,2)'*p_389) + ...
      square_abs(H(:,2)'*p_456) + square_abs(H(:,2)'*p_457) + square_abs(H(:,2)'*p_458) + square_abs(H(:,2)'*p_459) + ...
      square_abs(H(:,2)'*p_467) + square_abs(H(:,2)'*p_468) + square_abs(H(:,2)'*p_469) + square_abs(H(:,2)'*p_478) + ...
      square_abs(H(:,2)'*p_479) + square_abs(H(:,2)'*p_489) + square_abs(H(:,2)'*p_567) + square_abs(H(:,2)'*p_568) + ...
      square_abs(H(:,2)'*p_569) + square_abs(H(:,2)'*p_578) + square_abs(H(:,2)'*p_579) + square_abs(H(:,2)'*p_589) + ...
      square_abs(H(:,2)'*p_678) + square_abs(H(:,2)'*p_679) + square_abs(H(:,2)'*p_689) + square_abs(H(:,2)'*p_789) + ...
      square_abs(H(:,2)'*p_1) + square_abs(H(:,2)'*p_2) + square_abs(H(:,2)'*p_3) + ...
      square_abs(H(:,2)'*p_4) + square_abs(H(:,2)'*p_5) + square_abs(H(:,2)'*p_6) + ...
      square_abs(H(:,2)'*p_7) + square_abs(H(:,2)'*p_8) + square_abs(H(:,2)'*p_9) + 1;

T_2_189 = square_abs(H(:,2)'*p_189) + T_2;
T_2_179 = square_abs(H(:,2)'*p_179) + T_2_189;
T_2_178 = square_abs(H(:,2)'*p_178) + T_2_179;
T_2_169 = square_abs(H(:,2)'*p_169) + T_2_178;
T_2_168 = square_abs(H(:,2)'*p_168) + T_2_169;
T_2_167 = square_abs(H(:,2)'*p_167) + T_2_168;
T_2_159 = square_abs(H(:,2)'*p_159) + T_2_167;
T_2_158 = square_abs(H(:,2)'*p_158) + T_2_159;
T_2_157 = square_abs(H(:,2)'*p_157) + T_2_158;
T_2_156 = square_abs(H(:,2)'*p_156) + T_2_157;
T_2_149 = square_abs(H(:,2)'*p_149) + T_2_156;
T_2_148 = square_abs(H(:,2)'*p_148) + T_2_149;
T_2_147 = square_abs(H(:,2)'*p_147) + T_2_148;
T_2_146 = square_abs(H(:,2)'*p_146) + T_2_147;
T_2_145 = square_abs(H(:,2)'*p_145) + T_2_146;
T_2_139 = square_abs(H(:,2)'*p_139) + T_2_145;
T_2_138 = square_abs(H(:,2)'*p_138) + T_2_139;
T_2_137 = square_abs(H(:,2)'*p_137) + T_2_138;
T_2_136 = square_abs(H(:,2)'*p_136) + T_2_137;
T_2_135 = square_abs(H(:,2)'*p_135) + T_2_136;
T_2_134 = square_abs(H(:,2)'*p_134) + T_2_135;
T_2_129 = square_abs(H(:,2)'*p_129) + T_2_134;
T_2_128 = square_abs(H(:,2)'*p_128) + T_2_129;
T_2_127 = square_abs(H(:,2)'*p_127) + T_2_128;
T_2_126 = square_abs(H(:,2)'*p_126) + T_2_127;
T_2_125 = square_abs(H(:,2)'*p_125) + T_2_126;
T_2_124 = square_abs(H(:,2)'*p_124) + T_2_125;
T_2_123 = square_abs(H(:,2)'*p_123) + T_2_124;
T_2_123456789 = square_abs(H(:,2)'*p_123456789) + T_2_123;

T_3 = square_abs(H(:,3)'*p_134) + square_abs(H(:,3)'*p_135) + square_abs(H(:,3)'*p_136) + square_abs(H(:,3)'*p_137) + ...
      square_abs(H(:,3)'*p_138) + square_abs(H(:,3)'*p_139) + square_abs(H(:,3)'*p_145) + square_abs(H(:,3)'*p_146) + ...
      square_abs(H(:,3)'*p_147) + square_abs(H(:,3)'*p_148) + square_abs(H(:,3)'*p_149) + square_abs(H(:,3)'*p_156) + ...
      square_abs(H(:,3)'*p_157) + square_abs(H(:,3)'*p_158) + square_abs(H(:,3)'*p_159) + square_abs(H(:,3)'*p_167) + ...
      square_abs(H(:,3)'*p_168) + square_abs(H(:,3)'*p_169) + square_abs(H(:,3)'*p_178) + square_abs(H(:,3)'*p_179) + ...
      square_abs(H(:,3)'*p_189) + square_abs(H(:,3)'*p_345) + square_abs(H(:,3)'*p_346) + square_abs(H(:,3)'*p_347) + ...
      square_abs(H(:,3)'*p_348) + square_abs(H(:,3)'*p_349) + square_abs(H(:,3)'*p_356) + square_abs(H(:,3)'*p_357) + ...
      square_abs(H(:,3)'*p_358) + square_abs(H(:,3)'*p_359) + square_abs(H(:,3)'*p_367) + square_abs(H(:,3)'*p_368) + ...
      square_abs(H(:,3)'*p_369) + square_abs(H(:,3)'*p_378) + square_abs(H(:,3)'*p_379) + square_abs(H(:,3)'*p_389) + ...
      square_abs(H(:,3)'*p_456) + square_abs(H(:,3)'*p_457) + square_abs(H(:,3)'*p_458) + square_abs(H(:,3)'*p_459) + ...
      square_abs(H(:,3)'*p_467) + square_abs(H(:,3)'*p_468) + square_abs(H(:,3)'*p_469) + square_abs(H(:,3)'*p_478) + ...
      square_abs(H(:,3)'*p_479) + square_abs(H(:,3)'*p_489) + square_abs(H(:,3)'*p_567) + square_abs(H(:,3)'*p_568) + ...
      square_abs(H(:,3)'*p_569) + square_abs(H(:,3)'*p_578) + square_abs(H(:,3)'*p_579) + square_abs(H(:,3)'*p_589) + ...
      square_abs(H(:,3)'*p_678) + square_abs(H(:,3)'*p_679) + square_abs(H(:,3)'*p_689) + square_abs(H(:,3)'*p_789) + ...
      square_abs(H(:,3)'*p_1) + square_abs(H(:,3)'*p_2) + square_abs(H(:,3)'*p_3) + ...
      square_abs(H(:,3)'*p_4) + square_abs(H(:,3)'*p_5) + square_abs(H(:,3)'*p_6) + ...
      square_abs(H(:,3)'*p_7) + square_abs(H(:,3)'*p_8) + square_abs(H(:,3)'*p_9) + 1;

T_3_289 = square_abs(H(:,3)'*p_289) + T_3;
T_3_279 = square_abs(H(:,3)'*p_279) + T_3_289;
T_3_278 = square_abs(H(:,3)'*p_278) + T_3_279;
T_3_269 = square_abs(H(:,3)'*p_269) + T_3_278;
T_3_268 = square_abs(H(:,3)'*p_268) + T_3_269;
T_3_267 = square_abs(H(:,3)'*p_267) + T_3_268;
T_3_259 = square_abs(H(:,3)'*p_259) + T_3_267;
T_3_258 = square_abs(H(:,3)'*p_258) + T_3_259;
T_3_257 = square_abs(H(:,3)'*p_257) + T_3_258;
T_3_256 = square_abs(H(:,3)'*p_256) + T_3_257;
T_3_249 = square_abs(H(:,3)'*p_249) + T_3_256;
T_3_248 = square_abs(H(:,3)'*p_248) + T_3_249;
T_3_247 = square_abs(H(:,3)'*p_247) + T_3_248;
T_3_246 = square_abs(H(:,3)'*p_246) + T_3_247;
T_3_245 = square_abs(H(:,3)'*p_245) + T_3_246;
T_3_239 = square_abs(H(:,3)'*p_239) + T_3_245;
T_3_238 = square_abs(H(:,3)'*p_238) + T_3_239;
T_3_237 = square_abs(H(:,3)'*p_237) + T_3_238;
T_3_236 = square_abs(H(:,3)'*p_236) + T_3_237;
T_3_235 = square_abs(H(:,3)'*p_235) + T_3_236;
T_3_234 = square_abs(H(:,3)'*p_234) + T_3_235;
T_3_129 = square_abs(H(:,3)'*p_129) + T_3_234;
T_3_128 = square_abs(H(:,3)'*p_128) + T_3_129;
T_3_127 = square_abs(H(:,3)'*p_127) + T_3_128;
T_3_126 = square_abs(H(:,3)'*p_126) + T_3_127;
T_3_125 = square_abs(H(:,3)'*p_125) + T_3_126;
T_3_124 = square_abs(H(:,3)'*p_124) + T_3_125;
T_3_123 = square_abs(H(:,3)'*p_123) + T_3_124;
T_3_123456789 = square_abs(H(:,3)'*p_123456789) + T_3_123;

T_4 = square_abs(H(:,4)'*p_134) + square_abs(H(:,4)'*p_135) + square_abs(H(:,4)'*p_136) + square_abs(H(:,4)'*p_137) + ...
      square_abs(H(:,4)'*p_138) + square_abs(H(:,4)'*p_139) + square_abs(H(:,4)'*p_145) + square_abs(H(:,4)'*p_146) + ...
      square_abs(H(:,4)'*p_147) + square_abs(H(:,4)'*p_148) + square_abs(H(:,4)'*p_149) + square_abs(H(:,4)'*p_156) + ...
      square_abs(H(:,4)'*p_157) + square_abs(H(:,4)'*p_158) + square_abs(H(:,4)'*p_159) + square_abs(H(:,4)'*p_167) + ...
      square_abs(H(:,4)'*p_168) + square_abs(H(:,4)'*p_169) + square_abs(H(:,4)'*p_178) + square_abs(H(:,4)'*p_179) + ...
      square_abs(H(:,4)'*p_189) + square_abs(H(:,4)'*p_345) + square_abs(H(:,4)'*p_346) + square_abs(H(:,4)'*p_347) + ...
      square_abs(H(:,4)'*p_348) + square_abs(H(:,4)'*p_349) + square_abs(H(:,4)'*p_356) + square_abs(H(:,4)'*p_357) + ...
      square_abs(H(:,4)'*p_358) + square_abs(H(:,4)'*p_359) + square_abs(H(:,4)'*p_367) + square_abs(H(:,4)'*p_368) + ...
      square_abs(H(:,4)'*p_369) + square_abs(H(:,4)'*p_378) + square_abs(H(:,4)'*p_379) + square_abs(H(:,4)'*p_389) + ...
      square_abs(H(:,4)'*p_456) + square_abs(H(:,4)'*p_457) + square_abs(H(:,4)'*p_458) + square_abs(H(:,4)'*p_459) + ...
      square_abs(H(:,4)'*p_467) + square_abs(H(:,4)'*p_468) + square_abs(H(:,4)'*p_469) + square_abs(H(:,4)'*p_478) + ...
      square_abs(H(:,4)'*p_479) + square_abs(H(:,4)'*p_489) + square_abs(H(:,4)'*p_567) + square_abs(H(:,4)'*p_568) + ...
      square_abs(H(:,4)'*p_569) + square_abs(H(:,4)'*p_578) + square_abs(H(:,4)'*p_579) + square_abs(H(:,4)'*p_589) + ...
      square_abs(H(:,4)'*p_678) + square_abs(H(:,4)'*p_679) + square_abs(H(:,4)'*p_689) + square_abs(H(:,4)'*p_789) + ...
      square_abs(H(:,4)'*p_1) + square_abs(H(:,4)'*p_2) + square_abs(H(:,4)'*p_3) + ...
      square_abs(H(:,4)'*p_4) + square_abs(H(:,4)'*p_5) + square_abs(H(:,4)'*p_6) + ...
      square_abs(H(:,4)'*p_7) + square_abs(H(:,4)'*p_8) + square_abs(H(:,4)'*p_9) + 1;

T_4_289 = square_abs(H(:,4)'*p_289) + T_4;
T_4_279 = square_abs(H(:,4)'*p_279) + T_4_289;
T_4_278 = square_abs(H(:,4)'*p_278) + T_4_279;
T_4_269 = square_abs(H(:,4)'*p_269) + T_4_278;
T_4_268 = square_abs(H(:,4)'*p_268) + T_4_269;
T_4_267 = square_abs(H(:,4)'*p_267) + T_4_268;
T_4_259 = square_abs(H(:,4)'*p_259) + T_4_267;
T_4_258 = square_abs(H(:,4)'*p_258) + T_4_259;
T_4_257 = square_abs(H(:,4)'*p_257) + T_4_258;
T_4_256 = square_abs(H(:,4)'*p_256) + T_4_257;
T_4_249 = square_abs(H(:,4)'*p_249) + T_4_256;
T_4_248 = square_abs(H(:,4)'*p_248) + T_4_249;
T_4_247 = square_abs(H(:,4)'*p_247) + T_4_248;
T_4_246 = square_abs(H(:,4)'*p_246) + T_4_247;
T_4_245 = square_abs(H(:,4)'*p_245) + T_4_246;
T_4_239 = square_abs(H(:,4)'*p_239) + T_4_245;
T_4_238 = square_abs(H(:,4)'*p_238) + T_4_239;
T_4_237 = square_abs(H(:,4)'*p_237) + T_4_238;
T_4_236 = square_abs(H(:,4)'*p_236) + T_4_237;
T_4_235 = square_abs(H(:,4)'*p_235) + T_4_236;
T_4_234 = square_abs(H(:,4)'*p_234) + T_4_235;
T_4_129 = square_abs(H(:,4)'*p_129) + T_4_234;
T_4_128 = square_abs(H(:,4)'*p_128) + T_4_129;
T_4_127 = square_abs(H(:,4)'*p_127) + T_4_128;
T_4_126 = square_abs(H(:,4)'*p_126) + T_4_127;
T_4_125 = square_abs(H(:,4)'*p_125) + T_4_126;
T_4_124 = square_abs(H(:,4)'*p_124) + T_4_125;
T_4_123 = square_abs(H(:,4)'*p_123) + T_4_124;
T_4_123456789 = square_abs(H(:,4)'*p_123456789) + T_4_123;

T_5 = square_abs(H(:,5)'*p_124) + square_abs(H(:,5)'*p_125) + square_abs(H(:,5)'*p_126) + square_abs(H(:,5)'*p_127) + ...
      square_abs(H(:,5)'*p_128) + square_abs(H(:,5)'*p_129) + square_abs(H(:,5)'*p_145) + square_abs(H(:,5)'*p_146) + ...
      square_abs(H(:,5)'*p_147) + square_abs(H(:,5)'*p_148) + square_abs(H(:,5)'*p_149) + square_abs(H(:,5)'*p_156) + ...
      square_abs(H(:,5)'*p_157) + square_abs(H(:,5)'*p_158) + square_abs(H(:,5)'*p_159) + square_abs(H(:,5)'*p_167) + ...
      square_abs(H(:,5)'*p_168) + square_abs(H(:,5)'*p_169) + square_abs(H(:,5)'*p_178) + square_abs(H(:,5)'*p_179) + ...
      square_abs(H(:,5)'*p_189) + square_abs(H(:,5)'*p_245) + square_abs(H(:,5)'*p_246) + square_abs(H(:,5)'*p_247) + ...
      square_abs(H(:,5)'*p_248) + square_abs(H(:,5)'*p_249) + square_abs(H(:,5)'*p_256) + square_abs(H(:,5)'*p_257) + ...
      square_abs(H(:,5)'*p_258) + square_abs(H(:,5)'*p_259) + square_abs(H(:,5)'*p_267) + square_abs(H(:,5)'*p_268) + ...
      square_abs(H(:,5)'*p_269) + square_abs(H(:,5)'*p_278) + square_abs(H(:,5)'*p_279) + square_abs(H(:,5)'*p_289) + ...
      square_abs(H(:,5)'*p_456) + square_abs(H(:,5)'*p_457) + square_abs(H(:,5)'*p_458) + square_abs(H(:,5)'*p_459) + ...
      square_abs(H(:,5)'*p_467) + square_abs(H(:,5)'*p_468) + square_abs(H(:,5)'*p_469) + square_abs(H(:,5)'*p_478) + ...
      square_abs(H(:,5)'*p_479) + square_abs(H(:,5)'*p_489) + square_abs(H(:,5)'*p_567) + square_abs(H(:,5)'*p_568) + ...
      square_abs(H(:,5)'*p_569) + square_abs(H(:,5)'*p_578) + square_abs(H(:,5)'*p_579) + square_abs(H(:,5)'*p_589) + ...
      square_abs(H(:,5)'*p_678) + square_abs(H(:,5)'*p_679) + square_abs(H(:,5)'*p_689) + square_abs(H(:,5)'*p_789) + ...
      square_abs(H(:,5)'*p_1) + square_abs(H(:,5)'*p_2) + square_abs(H(:,5)'*p_3) + ...
      square_abs(H(:,5)'*p_4) + square_abs(H(:,5)'*p_5) + square_abs(H(:,5)'*p_6) + ...
      square_abs(H(:,5)'*p_7) + square_abs(H(:,5)'*p_8) + square_abs(H(:,5)'*p_9) + 1;

T_5_389 = square_abs(H(:,5)'*p_389) + T_5;
T_5_379 = square_abs(H(:,5)'*p_379) + T_5_389;
T_5_378 = square_abs(H(:,5)'*p_378) + T_5_379;
T_5_369 = square_abs(H(:,5)'*p_369) + T_5_378;
T_5_368 = square_abs(H(:,5)'*p_368) + T_5_369;
T_5_367 = square_abs(H(:,5)'*p_367) + T_5_368;
T_5_359 = square_abs(H(:,5)'*p_359) + T_5_367;
T_5_358 = square_abs(H(:,5)'*p_358) + T_5_359;
T_5_357 = square_abs(H(:,5)'*p_357) + T_5_358;
T_5_356 = square_abs(H(:,5)'*p_356) + T_5_357;
T_5_349 = square_abs(H(:,5)'*p_349) + T_5_356;
T_5_348 = square_abs(H(:,5)'*p_348) + T_5_349;
T_5_347 = square_abs(H(:,5)'*p_347) + T_5_348;
T_5_346 = square_abs(H(:,5)'*p_346) + T_5_347;
T_5_345 = square_abs(H(:,5)'*p_345) + T_5_346;
T_5_239 = square_abs(H(:,5)'*p_239) + T_5_345;
T_5_238 = square_abs(H(:,5)'*p_238) + T_5_239;
T_5_237 = square_abs(H(:,5)'*p_237) + T_5_238;
T_5_236 = square_abs(H(:,5)'*p_236) + T_5_237;
T_5_235 = square_abs(H(:,5)'*p_235) + T_5_236;
T_5_234 = square_abs(H(:,5)'*p_234) + T_5_235;
T_5_139 = square_abs(H(:,5)'*p_139) + T_5_234;
T_5_138 = square_abs(H(:,5)'*p_138) + T_5_139;
T_5_137 = square_abs(H(:,5)'*p_137) + T_5_138;
T_5_136 = square_abs(H(:,5)'*p_136) + T_5_137;
T_5_135 = square_abs(H(:,5)'*p_135) + T_5_136;
T_5_134 = square_abs(H(:,5)'*p_134) + T_5_135;
T_5_123 = square_abs(H(:,5)'*p_123) + T_5_134;
T_5_123456789 = square_abs(H(:,5)'*p_123456789) + T_5_123;

T_6 = square_abs(H(:,6)'*p_124) + square_abs(H(:,6)'*p_125) + square_abs(H(:,6)'*p_126) + square_abs(H(:,6)'*p_127) + ...
      square_abs(H(:,6)'*p_128) + square_abs(H(:,6)'*p_129) + square_abs(H(:,6)'*p_145) + square_abs(H(:,6)'*p_146) + ...
      square_abs(H(:,6)'*p_147) + square_abs(H(:,6)'*p_148) + square_abs(H(:,6)'*p_149) + square_abs(H(:,6)'*p_156) + ...
      square_abs(H(:,6)'*p_157) + square_abs(H(:,6)'*p_158) + square_abs(H(:,6)'*p_159) + square_abs(H(:,6)'*p_167) + ...
      square_abs(H(:,6)'*p_168) + square_abs(H(:,6)'*p_169) + square_abs(H(:,6)'*p_178) + square_abs(H(:,6)'*p_179) + ...
      square_abs(H(:,6)'*p_189) + square_abs(H(:,6)'*p_245) + square_abs(H(:,6)'*p_246) + square_abs(H(:,6)'*p_247) + ...
      square_abs(H(:,6)'*p_248) + square_abs(H(:,6)'*p_249) + square_abs(H(:,6)'*p_256) + square_abs(H(:,6)'*p_257) + ...
      square_abs(H(:,6)'*p_258) + square_abs(H(:,6)'*p_259) + square_abs(H(:,6)'*p_267) + square_abs(H(:,6)'*p_268) + ...
      square_abs(H(:,6)'*p_269) + square_abs(H(:,6)'*p_278) + square_abs(H(:,6)'*p_279) + square_abs(H(:,6)'*p_289) + ...
      square_abs(H(:,6)'*p_456) + square_abs(H(:,6)'*p_457) + square_abs(H(:,6)'*p_458) + square_abs(H(:,6)'*p_459) + ...
      square_abs(H(:,6)'*p_467) + square_abs(H(:,6)'*p_468) + square_abs(H(:,6)'*p_469) + square_abs(H(:,6)'*p_478) + ...
      square_abs(H(:,6)'*p_479) + square_abs(H(:,6)'*p_489) + square_abs(H(:,6)'*p_567) + square_abs(H(:,6)'*p_568) + ...
      square_abs(H(:,6)'*p_569) + square_abs(H(:,6)'*p_578) + square_abs(H(:,6)'*p_579) + square_abs(H(:,6)'*p_589) + ...
      square_abs(H(:,6)'*p_678) + square_abs(H(:,6)'*p_679) + square_abs(H(:,6)'*p_689) + square_abs(H(:,6)'*p_789) + ...
      square_abs(H(:,6)'*p_1) + square_abs(H(:,6)'*p_2) + square_abs(H(:,6)'*p_3) + ...
      square_abs(H(:,6)'*p_4) + square_abs(H(:,6)'*p_5) + square_abs(H(:,6)'*p_6) + ...
      square_abs(H(:,6)'*p_7) + square_abs(H(:,6)'*p_8) + square_abs(H(:,6)'*p_9) + 1;

T_6_389 = square_abs(H(:,6)'*p_389) + T_6;
T_6_379 = square_abs(H(:,6)'*p_379) + T_6_389;
T_6_378 = square_abs(H(:,6)'*p_378) + T_6_379;
T_6_369 = square_abs(H(:,6)'*p_369) + T_6_378;
T_6_368 = square_abs(H(:,6)'*p_368) + T_6_369;
T_6_367 = square_abs(H(:,6)'*p_367) + T_6_368;
T_6_359 = square_abs(H(:,6)'*p_359) + T_6_367;
T_6_358 = square_abs(H(:,6)'*p_358) + T_6_359;
T_6_357 = square_abs(H(:,6)'*p_357) + T_6_358;
T_6_356 = square_abs(H(:,6)'*p_356) + T_6_357;
T_6_349 = square_abs(H(:,6)'*p_349) + T_6_356;
T_6_348 = square_abs(H(:,6)'*p_348) + T_6_349;
T_6_347 = square_abs(H(:,6)'*p_347) + T_6_348;
T_6_346 = square_abs(H(:,6)'*p_346) + T_6_347;
T_6_345 = square_abs(H(:,6)'*p_345) + T_6_346;
T_6_239 = square_abs(H(:,6)'*p_239) + T_6_345;
T_6_238 = square_abs(H(:,6)'*p_238) + T_6_239;
T_6_237 = square_abs(H(:,6)'*p_237) + T_6_238;
T_6_236 = square_abs(H(:,6)'*p_236) + T_6_237;
T_6_235 = square_abs(H(:,6)'*p_235) + T_6_236;
T_6_234 = square_abs(H(:,6)'*p_234) + T_6_235;
T_6_139 = square_abs(H(:,6)'*p_139) + T_6_234;
T_6_138 = square_abs(H(:,6)'*p_138) + T_6_139;
T_6_137 = square_abs(H(:,6)'*p_137) + T_6_138;
T_6_136 = square_abs(H(:,6)'*p_136) + T_6_137;
T_6_135 = square_abs(H(:,6)'*p_135) + T_6_136;
T_6_134 = square_abs(H(:,6)'*p_134) + T_6_135;
T_6_123 = square_abs(H(:,6)'*p_123) + T_6_134;
T_6_123456789 = square_abs(H(:,6)'*p_123456789) + T_6_123;

T_7 = square_abs(H(:,7)'*p_123) + square_abs(H(:,7)'*p_125) + square_abs(H(:,7)'*p_126) + square_abs(H(:,7)'*p_127) + ...
      square_abs(H(:,7)'*p_128) + square_abs(H(:,7)'*p_129) + square_abs(H(:,7)'*p_135) + square_abs(H(:,7)'*p_136) + ...
      square_abs(H(:,7)'*p_137) + square_abs(H(:,7)'*p_138) + square_abs(H(:,7)'*p_139) + square_abs(H(:,7)'*p_156) + ...
      square_abs(H(:,7)'*p_157) + square_abs(H(:,7)'*p_158) + square_abs(H(:,7)'*p_159) + square_abs(H(:,7)'*p_167) + ...
      square_abs(H(:,7)'*p_168) + square_abs(H(:,7)'*p_169) + square_abs(H(:,7)'*p_178) + square_abs(H(:,7)'*p_179) + ...
      square_abs(H(:,7)'*p_189) + square_abs(H(:,7)'*p_235) + square_abs(H(:,7)'*p_236) + square_abs(H(:,7)'*p_237) + ...
      square_abs(H(:,7)'*p_238) + square_abs(H(:,7)'*p_239) + square_abs(H(:,7)'*p_256) + square_abs(H(:,7)'*p_257) + ...
      square_abs(H(:,7)'*p_258) + square_abs(H(:,7)'*p_259) + square_abs(H(:,7)'*p_267) + square_abs(H(:,7)'*p_268) + ...
      square_abs(H(:,7)'*p_269) + square_abs(H(:,7)'*p_278) + square_abs(H(:,7)'*p_279) + square_abs(H(:,7)'*p_289) + ...
      square_abs(H(:,7)'*p_356) + square_abs(H(:,7)'*p_357) + square_abs(H(:,7)'*p_358) + square_abs(H(:,7)'*p_359) + ...
      square_abs(H(:,7)'*p_367) + square_abs(H(:,7)'*p_368) + square_abs(H(:,7)'*p_369) + square_abs(H(:,7)'*p_378) + ...
      square_abs(H(:,7)'*p_379) + square_abs(H(:,7)'*p_389) + square_abs(H(:,7)'*p_567) + square_abs(H(:,7)'*p_568) + ...
      square_abs(H(:,7)'*p_569) + square_abs(H(:,7)'*p_578) + square_abs(H(:,7)'*p_579) + square_abs(H(:,7)'*p_589) + ...
      square_abs(H(:,7)'*p_678) + square_abs(H(:,7)'*p_679) + square_abs(H(:,7)'*p_689) + square_abs(H(:,7)'*p_789) + ...
      square_abs(H(:,7)'*p_1) + square_abs(H(:,7)'*p_2) + square_abs(H(:,7)'*p_3) + ...
      square_abs(H(:,7)'*p_4) + square_abs(H(:,7)'*p_5) + square_abs(H(:,7)'*p_6) + ...
      square_abs(H(:,7)'*p_7) + square_abs(H(:,7)'*p_8) + square_abs(H(:,7)'*p_9) + 1;

T_7_489 = square_abs(H(:,7)'*p_489) + T_7;
T_7_479 = square_abs(H(:,7)'*p_479) + T_7_489;
T_7_478 = square_abs(H(:,7)'*p_478) + T_7_479;
T_7_469 = square_abs(H(:,7)'*p_469) + T_7_478;
T_7_468 = square_abs(H(:,7)'*p_468) + T_7_469;
T_7_467 = square_abs(H(:,7)'*p_467) + T_7_468;
T_7_459 = square_abs(H(:,7)'*p_459) + T_7_467;
T_7_458 = square_abs(H(:,7)'*p_458) + T_7_459;
T_7_457 = square_abs(H(:,7)'*p_457) + T_7_458;
T_7_456 = square_abs(H(:,7)'*p_456) + T_7_457;
T_7_349 = square_abs(H(:,7)'*p_349) + T_7_456;
T_7_348 = square_abs(H(:,7)'*p_348) + T_7_349;
T_7_347 = square_abs(H(:,7)'*p_347) + T_7_348;
T_7_346 = square_abs(H(:,7)'*p_346) + T_7_347;
T_7_345 = square_abs(H(:,7)'*p_345) + T_7_346;
T_7_249 = square_abs(H(:,7)'*p_249) + T_7_345;
T_7_248 = square_abs(H(:,7)'*p_248) + T_7_249;
T_7_247 = square_abs(H(:,7)'*p_247) + T_7_248;
T_7_246 = square_abs(H(:,7)'*p_246) + T_7_247;
T_7_245 = square_abs(H(:,7)'*p_245) + T_7_246;
T_7_234 = square_abs(H(:,7)'*p_234) + T_7_245;
T_7_149 = square_abs(H(:,7)'*p_149) + T_7_234;
T_7_148 = square_abs(H(:,7)'*p_148) + T_7_149;
T_7_147 = square_abs(H(:,7)'*p_147) + T_7_148;
T_7_146 = square_abs(H(:,7)'*p_146) + T_7_147;
T_7_145 = square_abs(H(:,7)'*p_145) + T_7_146;
T_7_134 = square_abs(H(:,7)'*p_134) + T_7_145;
T_7_124 = square_abs(H(:,7)'*p_124) + T_7_134;
T_7_123456789 = square_abs(H(:,7)'*p_123456789) + T_7_124;

T_8 = square_abs(H(:,8)'*p_123) + square_abs(H(:,8)'*p_125) + square_abs(H(:,8)'*p_126) + square_abs(H(:,8)'*p_127) + ...
      square_abs(H(:,8)'*p_128) + square_abs(H(:,8)'*p_129) + square_abs(H(:,8)'*p_135) + square_abs(H(:,8)'*p_136) + ...
      square_abs(H(:,8)'*p_137) + square_abs(H(:,8)'*p_138) + square_abs(H(:,8)'*p_139) + square_abs(H(:,8)'*p_156) + ...
      square_abs(H(:,8)'*p_157) + square_abs(H(:,8)'*p_158) + square_abs(H(:,8)'*p_159) + square_abs(H(:,8)'*p_167) + ...
      square_abs(H(:,8)'*p_168) + square_abs(H(:,8)'*p_169) + square_abs(H(:,8)'*p_178) + square_abs(H(:,8)'*p_179) + ...
      square_abs(H(:,8)'*p_189) + square_abs(H(:,8)'*p_235) + square_abs(H(:,8)'*p_236) + square_abs(H(:,8)'*p_237) + ...
      square_abs(H(:,8)'*p_238) + square_abs(H(:,8)'*p_239) + square_abs(H(:,8)'*p_256) + square_abs(H(:,8)'*p_257) + ...
      square_abs(H(:,8)'*p_258) + square_abs(H(:,8)'*p_259) + square_abs(H(:,8)'*p_267) + square_abs(H(:,8)'*p_268) + ...
      square_abs(H(:,8)'*p_269) + square_abs(H(:,8)'*p_278) + square_abs(H(:,8)'*p_279) + square_abs(H(:,8)'*p_289) + ...
      square_abs(H(:,8)'*p_356) + square_abs(H(:,8)'*p_357) + square_abs(H(:,8)'*p_358) + square_abs(H(:,8)'*p_359) + ...
      square_abs(H(:,8)'*p_367) + square_abs(H(:,8)'*p_368) + square_abs(H(:,8)'*p_369) + square_abs(H(:,8)'*p_378) + ...
      square_abs(H(:,8)'*p_379) + square_abs(H(:,8)'*p_389) + square_abs(H(:,8)'*p_567) + square_abs(H(:,8)'*p_568) + ...
      square_abs(H(:,8)'*p_569) + square_abs(H(:,8)'*p_578) + square_abs(H(:,8)'*p_579) + square_abs(H(:,8)'*p_589) + ...
      square_abs(H(:,8)'*p_678) + square_abs(H(:,8)'*p_679) + square_abs(H(:,8)'*p_689) + square_abs(H(:,8)'*p_789) + ...
      square_abs(H(:,8)'*p_1) + square_abs(H(:,8)'*p_2) + square_abs(H(:,8)'*p_3) + ...
      square_abs(H(:,8)'*p_4) + square_abs(H(:,8)'*p_5) + square_abs(H(:,8)'*p_6) + ...
      square_abs(H(:,8)'*p_7) + square_abs(H(:,8)'*p_8) + square_abs(H(:,8)'*p_9) + 1;

T_8_489 = square_abs(H(:,8)'*p_489) + T_8;
T_8_479 = square_abs(H(:,8)'*p_479) + T_8_489;
T_8_478 = square_abs(H(:,8)'*p_478) + T_8_479;
T_8_469 = square_abs(H(:,8)'*p_469) + T_8_478;
T_8_468 = square_abs(H(:,8)'*p_468) + T_8_469;
T_8_467 = square_abs(H(:,8)'*p_467) + T_8_468;
T_8_459 = square_abs(H(:,8)'*p_459) + T_8_467;
T_8_458 = square_abs(H(:,8)'*p_458) + T_8_459;
T_8_457 = square_abs(H(:,8)'*p_457) + T_8_458;
T_8_456 = square_abs(H(:,8)'*p_456) + T_8_457;
T_8_349 = square_abs(H(:,8)'*p_349) + T_8_456;
T_8_348 = square_abs(H(:,8)'*p_348) + T_8_349;
T_8_347 = square_abs(H(:,8)'*p_347) + T_8_348;
T_8_346 = square_abs(H(:,8)'*p_346) + T_8_347;
T_8_345 = square_abs(H(:,8)'*p_345) + T_8_346;
T_8_249 = square_abs(H(:,8)'*p_249) + T_8_345;
T_8_248 = square_abs(H(:,8)'*p_248) + T_8_249;
T_8_247 = square_abs(H(:,8)'*p_247) + T_8_248;
T_8_246 = square_abs(H(:,8)'*p_246) + T_8_247;
T_8_245 = square_abs(H(:,8)'*p_245) + T_8_246;
T_8_234 = square_abs(H(:,8)'*p_234) + T_8_245;
T_8_149 = square_abs(H(:,8)'*p_149) + T_8_234;
T_8_148 = square_abs(H(:,8)'*p_148) + T_8_149;
T_8_147 = square_abs(H(:,8)'*p_147) + T_8_148;
T_8_146 = square_abs(H(:,8)'*p_146) + T_8_147;
T_8_145 = square_abs(H(:,8)'*p_145) + T_8_146;
T_8_134 = square_abs(H(:,8)'*p_134) + T_8_145;
T_8_124 = square_abs(H(:,8)'*p_124) + T_8_134;
T_8_123456789 = square_abs(H(:,8)'*p_123456789) + T_8_124;    

T_9 = square_abs(H(:,9)'*p_123) + square_abs(H(:,9)'*p_124) + square_abs(H(:,9)'*p_126) + square_abs(H(:,9)'*p_127) + ...
      square_abs(H(:,9)'*p_128) + square_abs(H(:,9)'*p_129) + square_abs(H(:,9)'*p_134) + square_abs(H(:,9)'*p_136) + ...
      square_abs(H(:,9)'*p_137) + square_abs(H(:,9)'*p_138) + square_abs(H(:,9)'*p_139) + square_abs(H(:,9)'*p_146) + ...
      square_abs(H(:,9)'*p_147) + square_abs(H(:,9)'*p_148) + square_abs(H(:,9)'*p_149) + square_abs(H(:,9)'*p_167) + ...
      square_abs(H(:,9)'*p_168) + square_abs(H(:,9)'*p_169) + square_abs(H(:,9)'*p_178) + square_abs(H(:,9)'*p_179) + ...
      square_abs(H(:,9)'*p_189) + square_abs(H(:,9)'*p_234) + square_abs(H(:,9)'*p_236) + square_abs(H(:,9)'*p_237) + ...
      square_abs(H(:,9)'*p_238) + square_abs(H(:,9)'*p_239) + square_abs(H(:,9)'*p_246) + square_abs(H(:,9)'*p_247) + ...
      square_abs(H(:,9)'*p_248) + square_abs(H(:,9)'*p_249) + square_abs(H(:,9)'*p_267) + square_abs(H(:,9)'*p_268) + ...
      square_abs(H(:,9)'*p_269) + square_abs(H(:,9)'*p_278) + square_abs(H(:,9)'*p_279) + square_abs(H(:,9)'*p_289) + ...
      square_abs(H(:,9)'*p_346) + square_abs(H(:,9)'*p_347) + square_abs(H(:,9)'*p_348) + square_abs(H(:,9)'*p_349) + ...
      square_abs(H(:,9)'*p_367) + square_abs(H(:,9)'*p_368) + square_abs(H(:,9)'*p_369) + square_abs(H(:,9)'*p_378) + ...
      square_abs(H(:,9)'*p_379) + square_abs(H(:,9)'*p_389) + square_abs(H(:,9)'*p_467) + square_abs(H(:,9)'*p_468) + ...
      square_abs(H(:,9)'*p_469) + square_abs(H(:,9)'*p_478) + square_abs(H(:,9)'*p_479) + square_abs(H(:,9)'*p_489) + ...
      square_abs(H(:,9)'*p_678) + square_abs(H(:,9)'*p_679) + square_abs(H(:,9)'*p_689) + square_abs(H(:,9)'*p_789) + ...
      square_abs(H(:,9)'*p_1) + square_abs(H(:,9)'*p_2) + square_abs(H(:,9)'*p_3) + ...
      square_abs(H(:,9)'*p_4) + square_abs(H(:,9)'*p_5) + square_abs(H(:,9)'*p_6) + ...
      square_abs(H(:,9)'*p_7) + square_abs(H(:,9)'*p_8) + square_abs(H(:,9)'*p_9) + 1;

T_9_589 = square_abs(H(:,9)'*p_589) + T_9;  
T_9_579 = square_abs(H(:,9)'*p_579) + T_9_589;  
T_9_578 = square_abs(H(:,9)'*p_578) + T_9_579;  
T_9_569 = square_abs(H(:,9)'*p_569) + T_9_578;  
T_9_568 = square_abs(H(:,9)'*p_568) + T_9_569;  
T_9_567 = square_abs(H(:,9)'*p_567) + T_9_568;  
T_9_459 = square_abs(H(:,9)'*p_459) + T_9_567;  
T_9_458 = square_abs(H(:,9)'*p_458) + T_9_459;  
T_9_457 = square_abs(H(:,9)'*p_457) + T_9_458;  
T_9_456 = square_abs(H(:,9)'*p_456) + T_9_457;  
T_9_359 = square_abs(H(:,9)'*p_359) + T_9_456;  
T_9_358 = square_abs(H(:,9)'*p_358) + T_9_359;  
T_9_357 = square_abs(H(:,9)'*p_357) + T_9_358;  
T_9_356 = square_abs(H(:,9)'*p_356) + T_9_357;  
T_9_345 = square_abs(H(:,9)'*p_345) + T_9_356;  
T_9_259 = square_abs(H(:,9)'*p_259) + T_9_345;  
T_9_258 = square_abs(H(:,9)'*p_258) + T_9_259;  
T_9_257 = square_abs(H(:,9)'*p_257) + T_9_258;  
T_9_256 = square_abs(H(:,9)'*p_256) + T_9_257;  
T_9_245 = square_abs(H(:,9)'*p_245) + T_9_256;  
T_9_235 = square_abs(H(:,9)'*p_235) + T_9_245;
T_9_159 = square_abs(H(:,9)'*p_159) + T_9_235;  
T_9_158 = square_abs(H(:,9)'*p_158) + T_9_159;  
T_9_157 = square_abs(H(:,9)'*p_157) + T_9_158;  
T_9_156 = square_abs(H(:,9)'*p_156) + T_9_157;  
T_9_145 = square_abs(H(:,9)'*p_145) + T_9_156;  
T_9_135 = square_abs(H(:,9)'*p_135) + T_9_145;
T_9_125 = square_abs(H(:,9)'*p_125) + T_9_135;
T_9_123456789 = square_abs(H(:,9)'*p_123456789) + T_9_125;

T_10 = square_abs(H(:,10)'*p_123) + square_abs(H(:,10)'*p_124) + square_abs(H(:,10)'*p_126) + square_abs(H(:,10)'*p_127) + ...
       square_abs(H(:,10)'*p_128) + square_abs(H(:,10)'*p_129) + square_abs(H(:,10)'*p_134) + square_abs(H(:,10)'*p_136) + ...
       square_abs(H(:,10)'*p_137) + square_abs(H(:,10)'*p_138) + square_abs(H(:,10)'*p_139) + square_abs(H(:,10)'*p_146) + ...
       square_abs(H(:,10)'*p_147) + square_abs(H(:,10)'*p_148) + square_abs(H(:,10)'*p_149) + square_abs(H(:,10)'*p_167) + ...
       square_abs(H(:,10)'*p_168) + square_abs(H(:,10)'*p_169) + square_abs(H(:,10)'*p_178) + square_abs(H(:,10)'*p_179) + ...
       square_abs(H(:,10)'*p_189) + square_abs(H(:,10)'*p_234) + square_abs(H(:,10)'*p_236) + square_abs(H(:,10)'*p_237) + ...
       square_abs(H(:,10)'*p_238) + square_abs(H(:,10)'*p_239) + square_abs(H(:,10)'*p_246) + square_abs(H(:,10)'*p_247) + ...
       square_abs(H(:,10)'*p_248) + square_abs(H(:,10)'*p_249) + square_abs(H(:,10)'*p_267) + square_abs(H(:,10)'*p_268) + ...
       square_abs(H(:,10)'*p_269) + square_abs(H(:,10)'*p_278) + square_abs(H(:,10)'*p_279) + square_abs(H(:,10)'*p_289) + ...
       square_abs(H(:,10)'*p_346) + square_abs(H(:,10)'*p_347) + square_abs(H(:,10)'*p_348) + square_abs(H(:,10)'*p_349) + ...
       square_abs(H(:,10)'*p_367) + square_abs(H(:,10)'*p_368) + square_abs(H(:,10)'*p_369) + square_abs(H(:,10)'*p_378) + ...
       square_abs(H(:,10)'*p_379) + square_abs(H(:,10)'*p_389) + square_abs(H(:,10)'*p_467) + square_abs(H(:,10)'*p_468) + ...
       square_abs(H(:,10)'*p_469) + square_abs(H(:,10)'*p_478) + square_abs(H(:,10)'*p_479) + square_abs(H(:,10)'*p_489) + ...
       square_abs(H(:,10)'*p_678) + square_abs(H(:,10)'*p_679) + square_abs(H(:,10)'*p_689) + square_abs(H(:,10)'*p_789) + ...
       square_abs(H(:,10)'*p_1) + square_abs(H(:,10)'*p_2) + square_abs(H(:,10)'*p_3) + ...
       square_abs(H(:,10)'*p_4) + square_abs(H(:,10)'*p_5) + square_abs(H(:,10)'*p_6) + ...
       square_abs(H(:,10)'*p_7) + square_abs(H(:,10)'*p_8) + square_abs(H(:,10)'*p_9) + 1;

T_10_589 = square_abs(H(:,10)'*p_589) + T_10;  
T_10_579 = square_abs(H(:,10)'*p_579) + T_10_589;  
T_10_578 = square_abs(H(:,10)'*p_578) + T_10_579;  
T_10_569 = square_abs(H(:,10)'*p_569) + T_10_578;  
T_10_568 = square_abs(H(:,10)'*p_568) + T_10_569;  
T_10_567 = square_abs(H(:,10)'*p_567) + T_10_568;  
T_10_459 = square_abs(H(:,10)'*p_459) + T_10_567;  
T_10_458 = square_abs(H(:,10)'*p_458) + T_10_459;  
T_10_457 = square_abs(H(:,10)'*p_457) + T_10_458;  
T_10_456 = square_abs(H(:,10)'*p_456) + T_10_457;  
T_10_359 = square_abs(H(:,10)'*p_359) + T_10_456;  
T_10_358 = square_abs(H(:,10)'*p_358) + T_10_359;  
T_10_357 = square_abs(H(:,10)'*p_357) + T_10_358;  
T_10_356 = square_abs(H(:,10)'*p_356) + T_10_357;  
T_10_345 = square_abs(H(:,10)'*p_345) + T_10_356;  
T_10_259 = square_abs(H(:,10)'*p_259) + T_10_345;  
T_10_258 = square_abs(H(:,10)'*p_258) + T_10_259;  
T_10_257 = square_abs(H(:,10)'*p_257) + T_10_258;  
T_10_256 = square_abs(H(:,10)'*p_256) + T_10_257;  
T_10_245 = square_abs(H(:,10)'*p_245) + T_10_256;  
T_10_235 = square_abs(H(:,10)'*p_235) + T_10_245;
T_10_159 = square_abs(H(:,10)'*p_159) + T_10_235;  
T_10_158 = square_abs(H(:,10)'*p_158) + T_10_159;  
T_10_157 = square_abs(H(:,10)'*p_157) + T_10_158;  
T_10_156 = square_abs(H(:,10)'*p_156) + T_10_157;  
T_10_145 = square_abs(H(:,10)'*p_145) + T_10_156;  
T_10_135 = square_abs(H(:,10)'*p_135) + T_10_145;
T_10_125 = square_abs(H(:,10)'*p_125) + T_10_135;
T_10_123456789 = square_abs(H(:,10)'*p_123456789) + T_10_125;

T_11 = square_abs(H(:,11)'*p_123) + square_abs(H(:,11)'*p_124) + square_abs(H(:,11)'*p_125) + square_abs(H(:,11)'*p_127) + ...
       square_abs(H(:,11)'*p_128) + square_abs(H(:,11)'*p_129) + square_abs(H(:,11)'*p_134) + square_abs(H(:,11)'*p_135) + ...
       square_abs(H(:,11)'*p_137) + square_abs(H(:,11)'*p_138) + square_abs(H(:,11)'*p_139) + square_abs(H(:,11)'*p_145) + ...
       square_abs(H(:,11)'*p_147) + square_abs(H(:,11)'*p_148) + square_abs(H(:,11)'*p_149) + square_abs(H(:,11)'*p_157) + ...
       square_abs(H(:,11)'*p_158) + square_abs(H(:,11)'*p_159) + square_abs(H(:,11)'*p_178) + square_abs(H(:,11)'*p_179) + ...
       square_abs(H(:,11)'*p_189) + square_abs(H(:,11)'*p_234) + square_abs(H(:,11)'*p_235) + square_abs(H(:,11)'*p_237) + ...
       square_abs(H(:,11)'*p_238) + square_abs(H(:,11)'*p_239) + square_abs(H(:,11)'*p_245) + square_abs(H(:,11)'*p_247) + ...
       square_abs(H(:,11)'*p_248) + square_abs(H(:,11)'*p_249) + square_abs(H(:,11)'*p_257) + square_abs(H(:,11)'*p_258) + ...
       square_abs(H(:,11)'*p_259) + square_abs(H(:,11)'*p_278) + square_abs(H(:,11)'*p_279) + square_abs(H(:,11)'*p_289) + ...
       square_abs(H(:,11)'*p_345) + square_abs(H(:,11)'*p_347) + square_abs(H(:,11)'*p_348) + square_abs(H(:,11)'*p_349) + ...
       square_abs(H(:,11)'*p_357) + square_abs(H(:,11)'*p_358) + square_abs(H(:,11)'*p_359) + square_abs(H(:,11)'*p_378) + ...
       square_abs(H(:,11)'*p_379) + square_abs(H(:,11)'*p_389) + square_abs(H(:,11)'*p_457) + square_abs(H(:,11)'*p_458) + ...
       square_abs(H(:,11)'*p_459) + square_abs(H(:,11)'*p_478) + square_abs(H(:,11)'*p_479) + square_abs(H(:,11)'*p_489) + ...
       square_abs(H(:,11)'*p_578) + square_abs(H(:,11)'*p_579) + square_abs(H(:,11)'*p_589) + square_abs(H(:,11)'*p_789) + ...
       square_abs(H(:,11)'*p_1) + square_abs(H(:,11)'*p_2) + square_abs(H(:,11)'*p_3) + ...
       square_abs(H(:,11)'*p_4) + square_abs(H(:,11)'*p_5) + square_abs(H(:,11)'*p_6) + ...
       square_abs(H(:,11)'*p_7) + square_abs(H(:,11)'*p_8) + square_abs(H(:,11)'*p_9) + 1;

T_11_689 = square_abs(H(:,11)'*p_689) + T_11;
T_11_679 = square_abs(H(:,11)'*p_679) + T_11_689;
T_11_678 = square_abs(H(:,11)'*p_678) + T_11_679;
T_11_569 = square_abs(H(:,11)'*p_569) + T_11_678;
T_11_568 = square_abs(H(:,11)'*p_568) + T_11_569;
T_11_567 = square_abs(H(:,11)'*p_567) + T_11_568;
T_11_469 = square_abs(H(:,11)'*p_469) + T_11_567;
T_11_468 = square_abs(H(:,11)'*p_468) + T_11_469;
T_11_467 = square_abs(H(:,11)'*p_467) + T_11_468;
T_11_456 = square_abs(H(:,11)'*p_456) + T_11_467;
T_11_369 = square_abs(H(:,11)'*p_369) + T_11_456;
T_11_368 = square_abs(H(:,11)'*p_368) + T_11_369;
T_11_367 = square_abs(H(:,11)'*p_367) + T_11_368;
T_11_356 = square_abs(H(:,11)'*p_356) + T_11_367;
T_11_346 = square_abs(H(:,11)'*p_346) + T_11_356;
T_11_269 = square_abs(H(:,11)'*p_269) + T_11_346;
T_11_268 = square_abs(H(:,11)'*p_268) + T_11_269;
T_11_267 = square_abs(H(:,11)'*p_267) + T_11_268;
T_11_256 = square_abs(H(:,11)'*p_256) + T_11_267;
T_11_246 = square_abs(H(:,11)'*p_246) + T_11_256;
T_11_236 = square_abs(H(:,11)'*p_236) + T_11_246;
T_11_169 = square_abs(H(:,11)'*p_169) + T_11_236;
T_11_168 = square_abs(H(:,11)'*p_168) + T_11_169;
T_11_167 = square_abs(H(:,11)'*p_167) + T_11_168;
T_11_156 = square_abs(H(:,11)'*p_156) + T_11_167;
T_11_146 = square_abs(H(:,11)'*p_146) + T_11_156;
T_11_136 = square_abs(H(:,11)'*p_136) + T_11_146;
T_11_126 = square_abs(H(:,11)'*p_126) + T_11_136;
T_11_123456789 = square_abs(H(:,11)'*p_123456789) + T_11_126;

T_12 = square_abs(H(:,12)'*p_123) + square_abs(H(:,12)'*p_124) + square_abs(H(:,12)'*p_125) + square_abs(H(:,12)'*p_127) + ...
       square_abs(H(:,12)'*p_128) + square_abs(H(:,12)'*p_129) + square_abs(H(:,12)'*p_134) + square_abs(H(:,12)'*p_135) + ...
       square_abs(H(:,12)'*p_137) + square_abs(H(:,12)'*p_138) + square_abs(H(:,12)'*p_139) + square_abs(H(:,12)'*p_145) + ...
       square_abs(H(:,12)'*p_147) + square_abs(H(:,12)'*p_148) + square_abs(H(:,12)'*p_149) + square_abs(H(:,12)'*p_157) + ...
       square_abs(H(:,12)'*p_158) + square_abs(H(:,12)'*p_159) + square_abs(H(:,12)'*p_178) + square_abs(H(:,12)'*p_179) + ...
       square_abs(H(:,12)'*p_189) + square_abs(H(:,12)'*p_234) + square_abs(H(:,12)'*p_235) + square_abs(H(:,12)'*p_237) + ...
       square_abs(H(:,12)'*p_238) + square_abs(H(:,12)'*p_239) + square_abs(H(:,12)'*p_245) + square_abs(H(:,12)'*p_247) + ...
       square_abs(H(:,12)'*p_248) + square_abs(H(:,12)'*p_249) + square_abs(H(:,12)'*p_257) + square_abs(H(:,12)'*p_258) + ...
       square_abs(H(:,12)'*p_259) + square_abs(H(:,12)'*p_278) + square_abs(H(:,12)'*p_279) + square_abs(H(:,12)'*p_289) + ...
       square_abs(H(:,12)'*p_345) + square_abs(H(:,12)'*p_347) + square_abs(H(:,12)'*p_348) + square_abs(H(:,12)'*p_349) + ...
       square_abs(H(:,12)'*p_357) + square_abs(H(:,12)'*p_358) + square_abs(H(:,12)'*p_359) + square_abs(H(:,12)'*p_378) + ...
       square_abs(H(:,12)'*p_379) + square_abs(H(:,12)'*p_389) + square_abs(H(:,12)'*p_457) + square_abs(H(:,12)'*p_458) + ...
       square_abs(H(:,12)'*p_459) + square_abs(H(:,12)'*p_478) + square_abs(H(:,12)'*p_479) + square_abs(H(:,12)'*p_489) + ...
       square_abs(H(:,12)'*p_578) + square_abs(H(:,12)'*p_579) + square_abs(H(:,12)'*p_589) + square_abs(H(:,12)'*p_789) + ...
       square_abs(H(:,12)'*p_1) + square_abs(H(:,12)'*p_2) + square_abs(H(:,12)'*p_3) + ...
       square_abs(H(:,12)'*p_4) + square_abs(H(:,12)'*p_5) + square_abs(H(:,12)'*p_6) + ...
       square_abs(H(:,12)'*p_7) + square_abs(H(:,12)'*p_8) + square_abs(H(:,12)'*p_9) + 1;

T_12_689 = square_abs(H(:,12)'*p_689) + T_12;
T_12_679 = square_abs(H(:,12)'*p_679) + T_12_689;
T_12_678 = square_abs(H(:,12)'*p_678) + T_12_679;
T_12_569 = square_abs(H(:,12)'*p_569) + T_12_678;
T_12_568 = square_abs(H(:,12)'*p_568) + T_12_569;
T_12_567 = square_abs(H(:,12)'*p_567) + T_12_568;
T_12_469 = square_abs(H(:,12)'*p_469) + T_12_567;
T_12_468 = square_abs(H(:,12)'*p_468) + T_12_469;
T_12_467 = square_abs(H(:,12)'*p_467) + T_12_468;
T_12_456 = square_abs(H(:,12)'*p_456) + T_12_467;
T_12_369 = square_abs(H(:,12)'*p_369) + T_12_456;
T_12_368 = square_abs(H(:,12)'*p_368) + T_12_369;
T_12_367 = square_abs(H(:,12)'*p_367) + T_12_368;
T_12_356 = square_abs(H(:,12)'*p_356) + T_12_367;
T_12_346 = square_abs(H(:,12)'*p_346) + T_12_356;
T_12_269 = square_abs(H(:,12)'*p_269) + T_12_346;
T_12_268 = square_abs(H(:,12)'*p_268) + T_12_269;
T_12_267 = square_abs(H(:,12)'*p_267) + T_12_268;
T_12_256 = square_abs(H(:,12)'*p_256) + T_12_267;
T_12_246 = square_abs(H(:,12)'*p_246) + T_12_256;
T_12_236 = square_abs(H(:,12)'*p_236) + T_12_246;
T_12_169 = square_abs(H(:,12)'*p_169) + T_12_236;
T_12_168 = square_abs(H(:,12)'*p_168) + T_12_169;
T_12_167 = square_abs(H(:,12)'*p_167) + T_12_168;
T_12_156 = square_abs(H(:,12)'*p_156) + T_12_167;
T_12_146 = square_abs(H(:,12)'*p_146) + T_12_156;
T_12_136 = square_abs(H(:,12)'*p_136) + T_12_146;
T_12_126 = square_abs(H(:,12)'*p_126) + T_12_136;
T_12_123456789 = square_abs(H(:,12)'*p_123456789) + T_12_126;

T_13 = square_abs(H(:,13)'*p_123) + square_abs(H(:,13)'*p_124) + square_abs(H(:,13)'*p_125) + square_abs(H(:,13)'*p_126) + ...
       square_abs(H(:,13)'*p_128) + square_abs(H(:,13)'*p_129) + square_abs(H(:,13)'*p_134) + square_abs(H(:,13)'*p_135) + ...
       square_abs(H(:,13)'*p_136) + square_abs(H(:,13)'*p_138) + square_abs(H(:,13)'*p_139) + square_abs(H(:,13)'*p_145) + ...
       square_abs(H(:,13)'*p_146) + square_abs(H(:,13)'*p_148) + square_abs(H(:,13)'*p_149) + square_abs(H(:,13)'*p_156) + ...
       square_abs(H(:,13)'*p_158) + square_abs(H(:,13)'*p_159) + square_abs(H(:,13)'*p_168) + square_abs(H(:,13)'*p_169) + ...
       square_abs(H(:,13)'*p_189) + square_abs(H(:,13)'*p_234) + square_abs(H(:,13)'*p_235) + square_abs(H(:,13)'*p_236) + ...
       square_abs(H(:,13)'*p_238) + square_abs(H(:,13)'*p_239) + square_abs(H(:,13)'*p_245) + square_abs(H(:,13)'*p_246) + ...
       square_abs(H(:,13)'*p_248) + square_abs(H(:,13)'*p_249) + square_abs(H(:,13)'*p_256) + square_abs(H(:,13)'*p_258) + ...
       square_abs(H(:,13)'*p_259) + square_abs(H(:,13)'*p_268) + square_abs(H(:,13)'*p_269) + square_abs(H(:,13)'*p_289) + ...
       square_abs(H(:,13)'*p_345) + square_abs(H(:,13)'*p_346) + square_abs(H(:,13)'*p_348) + square_abs(H(:,13)'*p_349) + ...
       square_abs(H(:,13)'*p_356) + square_abs(H(:,13)'*p_358) + square_abs(H(:,13)'*p_359) + square_abs(H(:,13)'*p_368) + ...
       square_abs(H(:,13)'*p_369) + square_abs(H(:,13)'*p_389) + square_abs(H(:,13)'*p_456) + square_abs(H(:,13)'*p_458) + ...
       square_abs(H(:,13)'*p_459) + square_abs(H(:,13)'*p_468) + square_abs(H(:,13)'*p_469) + square_abs(H(:,13)'*p_489) + ...
       square_abs(H(:,13)'*p_568) + square_abs(H(:,13)'*p_569) + square_abs(H(:,13)'*p_589) + square_abs(H(:,13)'*p_689) + ...
       square_abs(H(:,13)'*p_1) + square_abs(H(:,13)'*p_2) + square_abs(H(:,13)'*p_3) + ...
       square_abs(H(:,13)'*p_4) + square_abs(H(:,13)'*p_5) + square_abs(H(:,13)'*p_6) + ...
       square_abs(H(:,13)'*p_7) + square_abs(H(:,13)'*p_8) + square_abs(H(:,13)'*p_9) + 1;

T_13_789 = square_abs(H(:,13)'*p_789) + T_13;
T_13_679 = square_abs(H(:,13)'*p_679) + T_13_789;
T_13_678 = square_abs(H(:,13)'*p_678) + T_13_679;
T_13_579 = square_abs(H(:,13)'*p_579) + T_13_678;
T_13_578 = square_abs(H(:,13)'*p_578) + T_13_579;
T_13_567 = square_abs(H(:,13)'*p_567) + T_13_578;
T_13_479 = square_abs(H(:,13)'*p_479) + T_13_567;
T_13_478 = square_abs(H(:,13)'*p_478) + T_13_479;
T_13_467 = square_abs(H(:,13)'*p_467) + T_13_478;
T_13_457 = square_abs(H(:,13)'*p_457) + T_13_467;
T_13_379 = square_abs(H(:,13)'*p_379) + T_13_457;
T_13_378 = square_abs(H(:,13)'*p_378) + T_13_379;
T_13_367 = square_abs(H(:,13)'*p_367) + T_13_378;
T_13_357 = square_abs(H(:,13)'*p_357) + T_13_367;
T_13_347 = square_abs(H(:,13)'*p_347) + T_13_357;
T_13_279 = square_abs(H(:,13)'*p_279) + T_13_347;
T_13_278 = square_abs(H(:,13)'*p_278) + T_13_279;
T_13_267 = square_abs(H(:,13)'*p_267) + T_13_278;
T_13_257 = square_abs(H(:,13)'*p_257) + T_13_267;
T_13_247 = square_abs(H(:,13)'*p_247) + T_13_257;
T_13_237 = square_abs(H(:,13)'*p_237) + T_13_247;
T_13_179 = square_abs(H(:,13)'*p_179) + T_13_237;
T_13_178 = square_abs(H(:,13)'*p_178) + T_13_179;
T_13_167 = square_abs(H(:,13)'*p_167) + T_13_178;
T_13_157 = square_abs(H(:,13)'*p_157) + T_13_167;
T_13_147 = square_abs(H(:,13)'*p_147) + T_13_157;
T_13_137 = square_abs(H(:,13)'*p_137) + T_13_147;
T_13_127 = square_abs(H(:,13)'*p_127) + T_13_137;
T_13_123456789 = square_abs(H(:,13)'*p_123456789) + T_13_127;

T_14 = square_abs(H(:,14)'*p_123) + square_abs(H(:,14)'*p_124) + square_abs(H(:,14)'*p_125) + square_abs(H(:,14)'*p_126) + ...
       square_abs(H(:,14)'*p_128) + square_abs(H(:,14)'*p_129) + square_abs(H(:,14)'*p_134) + square_abs(H(:,14)'*p_135) + ...
       square_abs(H(:,14)'*p_136) + square_abs(H(:,14)'*p_138) + square_abs(H(:,14)'*p_139) + square_abs(H(:,14)'*p_145) + ...
       square_abs(H(:,14)'*p_146) + square_abs(H(:,14)'*p_148) + square_abs(H(:,14)'*p_149) + square_abs(H(:,14)'*p_156) + ...
       square_abs(H(:,14)'*p_158) + square_abs(H(:,14)'*p_159) + square_abs(H(:,14)'*p_168) + square_abs(H(:,14)'*p_169) + ...
       square_abs(H(:,14)'*p_189) + square_abs(H(:,14)'*p_234) + square_abs(H(:,14)'*p_235) + square_abs(H(:,14)'*p_236) + ...
       square_abs(H(:,14)'*p_238) + square_abs(H(:,14)'*p_239) + square_abs(H(:,14)'*p_245) + square_abs(H(:,14)'*p_246) + ...
       square_abs(H(:,14)'*p_248) + square_abs(H(:,14)'*p_249) + square_abs(H(:,14)'*p_256) + square_abs(H(:,14)'*p_258) + ...
       square_abs(H(:,14)'*p_259) + square_abs(H(:,14)'*p_268) + square_abs(H(:,14)'*p_269) + square_abs(H(:,14)'*p_289) + ...
       square_abs(H(:,14)'*p_345) + square_abs(H(:,14)'*p_346) + square_abs(H(:,14)'*p_348) + square_abs(H(:,14)'*p_349) + ...
       square_abs(H(:,14)'*p_356) + square_abs(H(:,14)'*p_358) + square_abs(H(:,14)'*p_359) + square_abs(H(:,14)'*p_368) + ...
       square_abs(H(:,14)'*p_369) + square_abs(H(:,14)'*p_389) + square_abs(H(:,14)'*p_456) + square_abs(H(:,14)'*p_458) + ...
       square_abs(H(:,14)'*p_459) + square_abs(H(:,14)'*p_468) + square_abs(H(:,14)'*p_469) + square_abs(H(:,14)'*p_489) + ...
       square_abs(H(:,14)'*p_568) + square_abs(H(:,14)'*p_569) + square_abs(H(:,14)'*p_589) + square_abs(H(:,14)'*p_689) + ...
       square_abs(H(:,14)'*p_1) + square_abs(H(:,14)'*p_2) + square_abs(H(:,14)'*p_3) + ...
       square_abs(H(:,14)'*p_4) + square_abs(H(:,14)'*p_5) + square_abs(H(:,14)'*p_6) + ...
       square_abs(H(:,14)'*p_7) + square_abs(H(:,14)'*p_8) + square_abs(H(:,14)'*p_9) + 1;

T_14_789 = square_abs(H(:,14)'*p_789) + T_14;
T_14_679 = square_abs(H(:,14)'*p_679) + T_14_789;
T_14_678 = square_abs(H(:,14)'*p_678) + T_14_679;
T_14_579 = square_abs(H(:,14)'*p_579) + T_14_678;
T_14_578 = square_abs(H(:,14)'*p_578) + T_14_579;
T_14_567 = square_abs(H(:,14)'*p_567) + T_14_578;
T_14_479 = square_abs(H(:,14)'*p_479) + T_14_567;
T_14_478 = square_abs(H(:,14)'*p_478) + T_14_479;
T_14_467 = square_abs(H(:,14)'*p_467) + T_14_478;
T_14_457 = square_abs(H(:,14)'*p_457) + T_14_467;
T_14_379 = square_abs(H(:,14)'*p_379) + T_14_457;
T_14_378 = square_abs(H(:,14)'*p_378) + T_14_379;
T_14_367 = square_abs(H(:,14)'*p_367) + T_14_378;
T_14_357 = square_abs(H(:,14)'*p_357) + T_14_367;
T_14_347 = square_abs(H(:,14)'*p_347) + T_14_357;
T_14_279 = square_abs(H(:,14)'*p_279) + T_14_347;
T_14_278 = square_abs(H(:,14)'*p_278) + T_14_279;
T_14_267 = square_abs(H(:,14)'*p_267) + T_14_278;
T_14_257 = square_abs(H(:,14)'*p_257) + T_14_267;
T_14_247 = square_abs(H(:,14)'*p_247) + T_14_257;
T_14_237 = square_abs(H(:,14)'*p_237) + T_14_247;
T_14_179 = square_abs(H(:,14)'*p_179) + T_14_237;
T_14_178 = square_abs(H(:,14)'*p_178) + T_14_179;
T_14_167 = square_abs(H(:,14)'*p_167) + T_14_178;
T_14_157 = square_abs(H(:,14)'*p_157) + T_14_167;
T_14_147 = square_abs(H(:,14)'*p_147) + T_14_157;
T_14_137 = square_abs(H(:,14)'*p_137) + T_14_147;
T_14_127 = square_abs(H(:,14)'*p_127) + T_14_137;
T_14_123456789 = square_abs(H(:,14)'*p_123456789) + T_14_127;

T_15 = square_abs(H(:,15)'*p_123) + square_abs(H(:,15)'*p_124) + square_abs(H(:,15)'*p_125) + square_abs(H(:,15)'*p_126) + ...
       square_abs(H(:,15)'*p_127) + square_abs(H(:,15)'*p_129) + square_abs(H(:,15)'*p_134) + square_abs(H(:,15)'*p_135) + ...
       square_abs(H(:,15)'*p_136) + square_abs(H(:,15)'*p_137) + square_abs(H(:,15)'*p_139) + square_abs(H(:,15)'*p_145) + ...
       square_abs(H(:,15)'*p_146) + square_abs(H(:,15)'*p_147) + square_abs(H(:,15)'*p_149) + square_abs(H(:,15)'*p_156) + ...
       square_abs(H(:,15)'*p_157) + square_abs(H(:,15)'*p_159) + square_abs(H(:,15)'*p_167) + square_abs(H(:,15)'*p_169) + ...
       square_abs(H(:,15)'*p_179) + square_abs(H(:,15)'*p_234) + square_abs(H(:,15)'*p_235) + square_abs(H(:,15)'*p_236) + ...
       square_abs(H(:,15)'*p_237) + square_abs(H(:,15)'*p_239) + square_abs(H(:,15)'*p_245) + square_abs(H(:,15)'*p_246) + ...
       square_abs(H(:,15)'*p_247) + square_abs(H(:,15)'*p_249) + square_abs(H(:,15)'*p_256) + square_abs(H(:,15)'*p_257) + ...
       square_abs(H(:,15)'*p_259) + square_abs(H(:,15)'*p_267) + square_abs(H(:,15)'*p_269) + square_abs(H(:,15)'*p_279) + ...
       square_abs(H(:,15)'*p_345) + square_abs(H(:,15)'*p_346) + square_abs(H(:,15)'*p_347) + square_abs(H(:,15)'*p_349) + ...
       square_abs(H(:,15)'*p_356) + square_abs(H(:,15)'*p_357) + square_abs(H(:,15)'*p_359) + square_abs(H(:,15)'*p_367) + ...
       square_abs(H(:,15)'*p_369) + square_abs(H(:,15)'*p_379) + square_abs(H(:,15)'*p_456) + square_abs(H(:,15)'*p_457) + ...
       square_abs(H(:,15)'*p_459) + square_abs(H(:,15)'*p_467) + square_abs(H(:,15)'*p_469) + square_abs(H(:,15)'*p_479) + ...
       square_abs(H(:,15)'*p_567) + square_abs(H(:,15)'*p_569) + square_abs(H(:,15)'*p_579) + square_abs(H(:,15)'*p_679) + ...
       square_abs(H(:,15)'*p_1) + square_abs(H(:,15)'*p_2) + square_abs(H(:,15)'*p_3) + ...
       square_abs(H(:,15)'*p_4) + square_abs(H(:,15)'*p_5) + square_abs(H(:,15)'*p_6) + ...
       square_abs(H(:,15)'*p_7) + square_abs(H(:,15)'*p_8) + square_abs(H(:,15)'*p_9) + 1;

T_15_789 = square_abs(H(:,15)'*p_789) + T_15;
T_15_689 = square_abs(H(:,15)'*p_689) + T_15_789;
T_15_678 = square_abs(H(:,15)'*p_678) + T_15_689;
T_15_589 = square_abs(H(:,15)'*p_589) + T_15_678;
T_15_578 = square_abs(H(:,15)'*p_578) + T_15_589;
T_15_568 = square_abs(H(:,15)'*p_568) + T_15_578;
T_15_489 = square_abs(H(:,15)'*p_489) + T_15_568;
T_15_478 = square_abs(H(:,15)'*p_478) + T_15_489;
T_15_468 = square_abs(H(:,15)'*p_468) + T_15_478;
T_15_458 = square_abs(H(:,15)'*p_458) + T_15_468;
T_15_389 = square_abs(H(:,15)'*p_389) + T_15_458;
T_15_378 = square_abs(H(:,15)'*p_378) + T_15_389;
T_15_368 = square_abs(H(:,15)'*p_368) + T_15_378;
T_15_358 = square_abs(H(:,15)'*p_358) + T_15_368;
T_15_348 = square_abs(H(:,15)'*p_348) + T_15_358;
T_15_289 = square_abs(H(:,15)'*p_289) + T_15_348;
T_15_278 = square_abs(H(:,15)'*p_278) + T_15_289;
T_15_268 = square_abs(H(:,15)'*p_268) + T_15_278;
T_15_258 = square_abs(H(:,15)'*p_258) + T_15_268;
T_15_248 = square_abs(H(:,15)'*p_248) + T_15_258;
T_15_238 = square_abs(H(:,15)'*p_238) + T_15_248;
T_15_189 = square_abs(H(:,15)'*p_189) + T_15_238;
T_15_178 = square_abs(H(:,15)'*p_178) + T_15_189;
T_15_168 = square_abs(H(:,15)'*p_168) + T_15_178;
T_15_158 = square_abs(H(:,15)'*p_158) + T_15_168;
T_15_148 = square_abs(H(:,15)'*p_148) + T_15_158;
T_15_138 = square_abs(H(:,15)'*p_138) + T_15_148;
T_15_128 = square_abs(H(:,15)'*p_128) + T_15_138;
T_15_123456789 = square_abs(H(:,15)'*p_123456789) + T_15_128;

T_16 = square_abs(H(:,16)'*p_123) + square_abs(H(:,16)'*p_124) + square_abs(H(:,16)'*p_125) + square_abs(H(:,16)'*p_126) + ...
       square_abs(H(:,16)'*p_127) + square_abs(H(:,16)'*p_129) + square_abs(H(:,16)'*p_134) + square_abs(H(:,16)'*p_135) + ...
       square_abs(H(:,16)'*p_136) + square_abs(H(:,16)'*p_137) + square_abs(H(:,16)'*p_139) + square_abs(H(:,16)'*p_145) + ...
       square_abs(H(:,16)'*p_146) + square_abs(H(:,16)'*p_147) + square_abs(H(:,16)'*p_149) + square_abs(H(:,16)'*p_156) + ...
       square_abs(H(:,16)'*p_157) + square_abs(H(:,16)'*p_159) + square_abs(H(:,16)'*p_167) + square_abs(H(:,16)'*p_169) + ...
       square_abs(H(:,16)'*p_179) + square_abs(H(:,16)'*p_234) + square_abs(H(:,16)'*p_235) + square_abs(H(:,16)'*p_236) + ...
       square_abs(H(:,16)'*p_237) + square_abs(H(:,16)'*p_239) + square_abs(H(:,16)'*p_245) + square_abs(H(:,16)'*p_246) + ...
       square_abs(H(:,16)'*p_247) + square_abs(H(:,16)'*p_249) + square_abs(H(:,16)'*p_256) + square_abs(H(:,16)'*p_257) + ...
       square_abs(H(:,16)'*p_259) + square_abs(H(:,16)'*p_267) + square_abs(H(:,16)'*p_269) + square_abs(H(:,16)'*p_279) + ...
       square_abs(H(:,16)'*p_345) + square_abs(H(:,16)'*p_346) + square_abs(H(:,16)'*p_347) + square_abs(H(:,16)'*p_349) + ...
       square_abs(H(:,16)'*p_356) + square_abs(H(:,16)'*p_357) + square_abs(H(:,16)'*p_359) + square_abs(H(:,16)'*p_367) + ...
       square_abs(H(:,16)'*p_369) + square_abs(H(:,16)'*p_379) + square_abs(H(:,16)'*p_456) + square_abs(H(:,16)'*p_457) + ...
       square_abs(H(:,16)'*p_459) + square_abs(H(:,16)'*p_467) + square_abs(H(:,16)'*p_469) + square_abs(H(:,16)'*p_479) + ...
       square_abs(H(:,16)'*p_567) + square_abs(H(:,16)'*p_569) + square_abs(H(:,16)'*p_579) + square_abs(H(:,16)'*p_679) + ...
       square_abs(H(:,16)'*p_1) + square_abs(H(:,16)'*p_2) + square_abs(H(:,16)'*p_3) + ...
       square_abs(H(:,16)'*p_4) + square_abs(H(:,16)'*p_5) + square_abs(H(:,16)'*p_6) + ...
       square_abs(H(:,16)'*p_7) + square_abs(H(:,16)'*p_8) + square_abs(H(:,16)'*p_9) + 1;

T_16_789 = square_abs(H(:,16)'*p_789) + T_16;
T_16_689 = square_abs(H(:,16)'*p_689) + T_16_789;
T_16_678 = square_abs(H(:,16)'*p_678) + T_16_689;
T_16_589 = square_abs(H(:,16)'*p_589) + T_16_678;
T_16_578 = square_abs(H(:,16)'*p_578) + T_16_589;
T_16_568 = square_abs(H(:,16)'*p_568) + T_16_578;
T_16_489 = square_abs(H(:,16)'*p_489) + T_16_568;
T_16_478 = square_abs(H(:,16)'*p_478) + T_16_489;
T_16_468 = square_abs(H(:,16)'*p_468) + T_16_478;
T_16_458 = square_abs(H(:,16)'*p_458) + T_16_468;
T_16_389 = square_abs(H(:,16)'*p_389) + T_16_458;
T_16_378 = square_abs(H(:,16)'*p_378) + T_16_389;
T_16_368 = square_abs(H(:,16)'*p_368) + T_16_378;
T_16_358 = square_abs(H(:,16)'*p_358) + T_16_368;
T_16_348 = square_abs(H(:,16)'*p_348) + T_16_358;
T_16_289 = square_abs(H(:,16)'*p_289) + T_16_348;
T_16_278 = square_abs(H(:,16)'*p_278) + T_16_289;
T_16_268 = square_abs(H(:,16)'*p_268) + T_16_278;
T_16_258 = square_abs(H(:,16)'*p_258) + T_16_268;
T_16_248 = square_abs(H(:,16)'*p_248) + T_16_258;
T_16_238 = square_abs(H(:,16)'*p_238) + T_16_248;
T_16_189 = square_abs(H(:,16)'*p_189) + T_16_238;
T_16_178 = square_abs(H(:,16)'*p_178) + T_16_189;
T_16_168 = square_abs(H(:,16)'*p_168) + T_16_178;
T_16_158 = square_abs(H(:,16)'*p_158) + T_16_168;
T_16_148 = square_abs(H(:,16)'*p_148) + T_16_158;
T_16_138 = square_abs(H(:,16)'*p_138) + T_16_148;
T_16_128 = square_abs(H(:,16)'*p_128) + T_16_138;
T_16_123456789 = square_abs(H(:,16)'*p_123456789) + T_16_128;

T_17 = square_abs(H(:,17)'*p_123) + square_abs(H(:,17)'*p_124) + square_abs(H(:,17)'*p_125) + square_abs(H(:,17)'*p_126) + ...
       square_abs(H(:,17)'*p_127) + square_abs(H(:,17)'*p_128) + square_abs(H(:,17)'*p_134) + square_abs(H(:,17)'*p_135) + ...
       square_abs(H(:,17)'*p_136) + square_abs(H(:,17)'*p_137) + square_abs(H(:,17)'*p_138) + square_abs(H(:,17)'*p_145) + ...
       square_abs(H(:,17)'*p_146) + square_abs(H(:,17)'*p_147) + square_abs(H(:,17)'*p_148) + square_abs(H(:,17)'*p_156) + ...
       square_abs(H(:,17)'*p_157) + square_abs(H(:,17)'*p_158) + square_abs(H(:,17)'*p_167) + square_abs(H(:,17)'*p_168) + ...
       square_abs(H(:,17)'*p_178) + square_abs(H(:,17)'*p_234) + square_abs(H(:,17)'*p_235) + square_abs(H(:,17)'*p_236) + ...
       square_abs(H(:,17)'*p_237) + square_abs(H(:,17)'*p_238) + square_abs(H(:,17)'*p_245) + square_abs(H(:,17)'*p_246) + ...
       square_abs(H(:,17)'*p_247) + square_abs(H(:,17)'*p_248) + square_abs(H(:,17)'*p_256) + square_abs(H(:,17)'*p_257) + ...
       square_abs(H(:,17)'*p_258) + square_abs(H(:,17)'*p_267) + square_abs(H(:,17)'*p_268) + square_abs(H(:,17)'*p_278) + ...
       square_abs(H(:,17)'*p_345) + square_abs(H(:,17)'*p_346) + square_abs(H(:,17)'*p_347) + square_abs(H(:,17)'*p_348) + ...
       square_abs(H(:,17)'*p_356) + square_abs(H(:,17)'*p_357) + square_abs(H(:,17)'*p_358) + square_abs(H(:,17)'*p_367) + ...
       square_abs(H(:,17)'*p_368) + square_abs(H(:,17)'*p_378) + square_abs(H(:,17)'*p_456) + square_abs(H(:,17)'*p_457) + ...
       square_abs(H(:,17)'*p_458) + square_abs(H(:,17)'*p_467) + square_abs(H(:,17)'*p_468) + square_abs(H(:,17)'*p_478) + ...
       square_abs(H(:,17)'*p_567) + square_abs(H(:,17)'*p_568) + square_abs(H(:,17)'*p_578) + square_abs(H(:,17)'*p_678) + ...
       square_abs(H(:,17)'*p_1) + square_abs(H(:,17)'*p_2) + square_abs(H(:,17)'*p_3) + ...
       square_abs(H(:,17)'*p_4) + square_abs(H(:,17)'*p_5) + square_abs(H(:,17)'*p_6) + ...
       square_abs(H(:,17)'*p_7) + square_abs(H(:,17)'*p_8) + square_abs(H(:,17)'*p_9) + 1;

T_17_789 = square_abs(H(:,17)'*p_789) + T_17;
T_17_689 = square_abs(H(:,17)'*p_689) + T_17_789;
T_17_679 = square_abs(H(:,17)'*p_679) + T_17_689;
T_17_589 = square_abs(H(:,17)'*p_589) + T_17_679;
T_17_579 = square_abs(H(:,17)'*p_579) + T_17_589;
T_17_569 = square_abs(H(:,17)'*p_569) + T_17_579;
T_17_489 = square_abs(H(:,17)'*p_489) + T_17_569;
T_17_479 = square_abs(H(:,17)'*p_479) + T_17_489;
T_17_469 = square_abs(H(:,17)'*p_469) + T_17_479;
T_17_459 = square_abs(H(:,17)'*p_459) + T_17_469;
T_17_389 = square_abs(H(:,17)'*p_389) + T_17_459;
T_17_379 = square_abs(H(:,17)'*p_379) + T_17_389;
T_17_369 = square_abs(H(:,17)'*p_369) + T_17_379;
T_17_359 = square_abs(H(:,17)'*p_359) + T_17_369;
T_17_349 = square_abs(H(:,17)'*p_349) + T_17_359;
T_17_289 = square_abs(H(:,17)'*p_289) + T_17_349;
T_17_279 = square_abs(H(:,17)'*p_279) + T_17_289;
T_17_269 = square_abs(H(:,17)'*p_269) + T_17_279;
T_17_259 = square_abs(H(:,17)'*p_259) + T_17_269;
T_17_249 = square_abs(H(:,17)'*p_249) + T_17_259;
T_17_239 = square_abs(H(:,17)'*p_239) + T_17_249;
T_17_189 = square_abs(H(:,17)'*p_189) + T_17_239;
T_17_179 = square_abs(H(:,17)'*p_179) + T_17_189;
T_17_169 = square_abs(H(:,17)'*p_169) + T_17_179;
T_17_159 = square_abs(H(:,17)'*p_159) + T_17_169;
T_17_149 = square_abs(H(:,17)'*p_149) + T_17_159;
T_17_139 = square_abs(H(:,17)'*p_139) + T_17_149;
T_17_129 = square_abs(H(:,17)'*p_129) + T_17_139;
T_17_123456789 = square_abs(H(:,17)'*p_123456789) + T_17_129;

T_18 = square_abs(H(:,18)'*p_123) + square_abs(H(:,18)'*p_124) + square_abs(H(:,18)'*p_125) + square_abs(H(:,18)'*p_126) + ...
       square_abs(H(:,18)'*p_127) + square_abs(H(:,18)'*p_128) + square_abs(H(:,18)'*p_134) + square_abs(H(:,18)'*p_135) + ...
       square_abs(H(:,18)'*p_136) + square_abs(H(:,18)'*p_137) + square_abs(H(:,18)'*p_138) + square_abs(H(:,18)'*p_145) + ...
       square_abs(H(:,18)'*p_146) + square_abs(H(:,18)'*p_147) + square_abs(H(:,18)'*p_148) + square_abs(H(:,18)'*p_156) + ...
       square_abs(H(:,18)'*p_157) + square_abs(H(:,18)'*p_158) + square_abs(H(:,18)'*p_167) + square_abs(H(:,18)'*p_168) + ...
       square_abs(H(:,18)'*p_178) + square_abs(H(:,18)'*p_234) + square_abs(H(:,18)'*p_235) + square_abs(H(:,18)'*p_236) + ...
       square_abs(H(:,18)'*p_237) + square_abs(H(:,18)'*p_238) + square_abs(H(:,18)'*p_245) + square_abs(H(:,18)'*p_246) + ...
       square_abs(H(:,18)'*p_247) + square_abs(H(:,18)'*p_248) + square_abs(H(:,18)'*p_256) + square_abs(H(:,18)'*p_257) + ...
       square_abs(H(:,18)'*p_258) + square_abs(H(:,18)'*p_267) + square_abs(H(:,18)'*p_268) + square_abs(H(:,18)'*p_278) + ...
       square_abs(H(:,18)'*p_345) + square_abs(H(:,18)'*p_346) + square_abs(H(:,18)'*p_347) + square_abs(H(:,18)'*p_348) + ...
       square_abs(H(:,18)'*p_356) + square_abs(H(:,18)'*p_357) + square_abs(H(:,18)'*p_358) + square_abs(H(:,18)'*p_367) + ...
       square_abs(H(:,18)'*p_368) + square_abs(H(:,18)'*p_378) + square_abs(H(:,18)'*p_456) + square_abs(H(:,18)'*p_457) + ...
       square_abs(H(:,18)'*p_458) + square_abs(H(:,18)'*p_467) + square_abs(H(:,18)'*p_468) + square_abs(H(:,18)'*p_478) + ...
       square_abs(H(:,18)'*p_567) + square_abs(H(:,18)'*p_568) + square_abs(H(:,18)'*p_578) + square_abs(H(:,18)'*p_678) + ...
       square_abs(H(:,18)'*p_1) + square_abs(H(:,18)'*p_2) + square_abs(H(:,18)'*p_3) + ...
       square_abs(H(:,18)'*p_4) + square_abs(H(:,18)'*p_5) + square_abs(H(:,18)'*p_6) + ...
       square_abs(H(:,18)'*p_7) + square_abs(H(:,18)'*p_8) + square_abs(H(:,18)'*p_9) + 1;

T_18_789 = square_abs(H(:,18)'*p_789) + T_18;
T_18_689 = square_abs(H(:,18)'*p_689) + T_18_789;
T_18_679 = square_abs(H(:,18)'*p_679) + T_18_689;
T_18_589 = square_abs(H(:,18)'*p_589) + T_18_679;
T_18_579 = square_abs(H(:,18)'*p_579) + T_18_589;
T_18_569 = square_abs(H(:,18)'*p_569) + T_18_579;
T_18_489 = square_abs(H(:,18)'*p_489) + T_18_569;
T_18_479 = square_abs(H(:,18)'*p_479) + T_18_489;
T_18_469 = square_abs(H(:,18)'*p_469) + T_18_479;
T_18_459 = square_abs(H(:,18)'*p_459) + T_18_469;
T_18_389 = square_abs(H(:,18)'*p_389) + T_18_459;
T_18_379 = square_abs(H(:,18)'*p_379) + T_18_389;
T_18_369 = square_abs(H(:,18)'*p_369) + T_18_379;
T_18_359 = square_abs(H(:,18)'*p_359) + T_18_369;
T_18_349 = square_abs(H(:,18)'*p_349) + T_18_359;
T_18_289 = square_abs(H(:,18)'*p_289) + T_18_349;
T_18_279 = square_abs(H(:,18)'*p_279) + T_18_289;
T_18_269 = square_abs(H(:,18)'*p_269) + T_18_279;
T_18_259 = square_abs(H(:,18)'*p_259) + T_18_269;
T_18_249 = square_abs(H(:,18)'*p_249) + T_18_259;
T_18_239 = square_abs(H(:,18)'*p_239) + T_18_249;
T_18_189 = square_abs(H(:,18)'*p_189) + T_18_239;
T_18_179 = square_abs(H(:,18)'*p_179) + T_18_189;
T_18_169 = square_abs(H(:,18)'*p_169) + T_18_179;
T_18_159 = square_abs(H(:,18)'*p_159) + T_18_169;
T_18_149 = square_abs(H(:,18)'*p_149) + T_18_159;
T_18_139 = square_abs(H(:,18)'*p_139) + T_18_149;
T_18_129 = square_abs(H(:,18)'*p_129) + T_18_139;
T_18_123456789 = square_abs(H(:,18)'*p_123456789) + T_18_129;


%MSE
E_1_123456789 = abs(g_1_123456789)^2*T_1_123456789 - 2*real(g_1_123456789*H(:,1)'*p_123456789) + 1;
E_1_123 = abs(g_1_123)^2*T_1_123 - 2*real(g_1_123*H(:,1)'*p_123) + 1;
E_1_124 = abs(g_1_124)^2*T_1_124 - 2*real(g_1_124*H(:,1)'*p_124) + 1;
E_1_125 = abs(g_1_125)^2*T_1_125 - 2*real(g_1_125*H(:,1)'*p_125) + 1;
E_1_126 = abs(g_1_126)^2*T_1_126 - 2*real(g_1_126*H(:,1)'*p_126) + 1;
E_1_127 = abs(g_1_127)^2*T_1_127 - 2*real(g_1_127*H(:,1)'*p_127) + 1;
E_1_128 = abs(g_1_128)^2*T_1_128 - 2*real(g_1_128*H(:,1)'*p_128) + 1;
E_1_129 = abs(g_1_129)^2*T_1_129 - 2*real(g_1_129*H(:,1)'*p_129) + 1;
E_1_134 = abs(g_1_134)^2*T_1_134 - 2*real(g_1_134*H(:,1)'*p_134) + 1;
E_1_135 = abs(g_1_135)^2*T_1_135 - 2*real(g_1_135*H(:,1)'*p_135) + 1;
E_1_136 = abs(g_1_136)^2*T_1_136 - 2*real(g_1_136*H(:,1)'*p_136) + 1;
E_1_137 = abs(g_1_137)^2*T_1_137 - 2*real(g_1_137*H(:,1)'*p_137) + 1;
E_1_138 = abs(g_1_138)^2*T_1_138 - 2*real(g_1_138*H(:,1)'*p_138) + 1;
E_1_139 = abs(g_1_139)^2*T_1_139 - 2*real(g_1_139*H(:,1)'*p_139) + 1;
E_1_145 = abs(g_1_145)^2*T_1_145 - 2*real(g_1_145*H(:,1)'*p_145) + 1;
E_1_146 = abs(g_1_146)^2*T_1_146 - 2*real(g_1_146*H(:,1)'*p_146) + 1;
E_1_147 = abs(g_1_147)^2*T_1_147 - 2*real(g_1_147*H(:,1)'*p_147) + 1;
E_1_148 = abs(g_1_148)^2*T_1_148 - 2*real(g_1_148*H(:,1)'*p_148) + 1;
E_1_149 = abs(g_1_149)^2*T_1_149 - 2*real(g_1_149*H(:,1)'*p_149) + 1;
E_1_156 = abs(g_1_156)^2*T_1_156 - 2*real(g_1_156*H(:,1)'*p_156) + 1;
E_1_157 = abs(g_1_157)^2*T_1_157 - 2*real(g_1_157*H(:,1)'*p_157) + 1;
E_1_158 = abs(g_1_158)^2*T_1_158 - 2*real(g_1_158*H(:,1)'*p_158) + 1;
E_1_159 = abs(g_1_159)^2*T_1_159 - 2*real(g_1_159*H(:,1)'*p_159) + 1;
E_1_167 = abs(g_1_167)^2*T_1_167 - 2*real(g_1_167*H(:,1)'*p_167) + 1;
E_1_168 = abs(g_1_168)^2*T_1_168 - 2*real(g_1_168*H(:,1)'*p_168) + 1;
E_1_169 = abs(g_1_169)^2*T_1_169 - 2*real(g_1_169*H(:,1)'*p_169) + 1;
E_1_178 = abs(g_1_178)^2*T_1_178 - 2*real(g_1_178*H(:,1)'*p_178) + 1;
E_1_179 = abs(g_1_179)^2*T_1_179 - 2*real(g_1_179*H(:,1)'*p_179) + 1;
E_1_189 = abs(g_1_189)^2*T_1_189 - 2*real(g_1_189*H(:,1)'*p_189) + 1;
E_1 = abs(g_1)^2*T_1 - 2*real(g_1*H(:,1)'*p_1) + 1;

E_2_123456789 = abs(g_2_123456789)^2*T_2_123456789 - 2*real(g_2_123456789*H(:,2)'*p_123456789) + 1;
E_2_123 = abs(g_2_123)^2*T_2_123 - 2*real(g_2_123*H(:,2)'*p_123) + 1;
E_2_124 = abs(g_2_124)^2*T_2_124 - 2*real(g_2_124*H(:,2)'*p_124) + 1;
E_2_125 = abs(g_2_125)^2*T_2_125 - 2*real(g_2_125*H(:,2)'*p_125) + 1;
E_2_126 = abs(g_2_126)^2*T_2_126 - 2*real(g_2_126*H(:,2)'*p_126) + 1;
E_2_127 = abs(g_2_127)^2*T_2_127 - 2*real(g_2_127*H(:,2)'*p_127) + 1;
E_2_128 = abs(g_2_128)^2*T_2_128 - 2*real(g_2_128*H(:,2)'*p_128) + 1;
E_2_129 = abs(g_2_129)^2*T_2_129 - 2*real(g_2_129*H(:,2)'*p_129) + 1;
E_2_134 = abs(g_2_134)^2*T_2_134 - 2*real(g_2_134*H(:,2)'*p_134) + 1;
E_2_135 = abs(g_2_135)^2*T_2_135 - 2*real(g_2_135*H(:,2)'*p_135) + 1;
E_2_136 = abs(g_2_136)^2*T_2_136 - 2*real(g_2_136*H(:,2)'*p_136) + 1;
E_2_137 = abs(g_2_137)^2*T_2_137 - 2*real(g_2_137*H(:,2)'*p_137) + 1;
E_2_138 = abs(g_2_138)^2*T_2_138 - 2*real(g_2_138*H(:,2)'*p_138) + 1;
E_2_139 = abs(g_2_139)^2*T_2_139 - 2*real(g_2_139*H(:,2)'*p_139) + 1;
E_2_145 = abs(g_2_145)^2*T_2_145 - 2*real(g_2_145*H(:,2)'*p_145) + 1;
E_2_146 = abs(g_2_146)^2*T_2_146 - 2*real(g_2_146*H(:,2)'*p_146) + 1;
E_2_147 = abs(g_2_147)^2*T_2_147 - 2*real(g_2_147*H(:,2)'*p_147) + 1;
E_2_148 = abs(g_2_148)^2*T_2_148 - 2*real(g_2_148*H(:,2)'*p_148) + 1;
E_2_149 = abs(g_2_149)^2*T_2_149 - 2*real(g_2_149*H(:,2)'*p_149) + 1;
E_2_156 = abs(g_2_156)^2*T_2_156 - 2*real(g_2_156*H(:,2)'*p_156) + 1;
E_2_157 = abs(g_2_157)^2*T_2_157 - 2*real(g_2_157*H(:,2)'*p_157) + 1;
E_2_158 = abs(g_2_158)^2*T_2_158 - 2*real(g_2_158*H(:,2)'*p_158) + 1;
E_2_159 = abs(g_2_159)^2*T_2_159 - 2*real(g_2_159*H(:,2)'*p_159) + 1;
E_2_167 = abs(g_2_167)^2*T_2_167 - 2*real(g_2_167*H(:,2)'*p_167) + 1;
E_2_168 = abs(g_2_168)^2*T_2_168 - 2*real(g_2_168*H(:,2)'*p_168) + 1;
E_2_169 = abs(g_2_169)^2*T_2_169 - 2*real(g_2_169*H(:,2)'*p_169) + 1;
E_2_178 = abs(g_2_178)^2*T_2_178 - 2*real(g_2_178*H(:,2)'*p_178) + 1;
E_2_179 = abs(g_2_179)^2*T_2_179 - 2*real(g_2_179*H(:,2)'*p_179) + 1;
E_2_189 = abs(g_2_189)^2*T_2_189 - 2*real(g_2_189*H(:,2)'*p_189) + 1;
E_2 = abs(g_2)^2*T_2 - 2*real(g_2*H(:,2)'*p_1) + 1;

E_3_123456789 = abs(g_3_123456789)^2*T_3_123456789 - 2*real(g_3_123456789*H(:,3)'*p_123456789) + 1;
E_3_123 = abs(g_3_123)^2*T_3_123 - 2*real(g_3_123*H(:,3)'*p_123) + 1;
E_3_124 = abs(g_3_124)^2*T_3_124 - 2*real(g_3_124*H(:,3)'*p_124) + 1;
E_3_125 = abs(g_3_125)^2*T_3_125 - 2*real(g_3_125*H(:,3)'*p_125) + 1;
E_3_126 = abs(g_3_126)^2*T_3_126 - 2*real(g_3_126*H(:,3)'*p_126) + 1;
E_3_127 = abs(g_3_127)^2*T_3_127 - 2*real(g_3_127*H(:,3)'*p_127) + 1;
E_3_128 = abs(g_3_128)^2*T_3_128 - 2*real(g_3_128*H(:,3)'*p_128) + 1;
E_3_129 = abs(g_3_129)^2*T_3_129 - 2*real(g_3_129*H(:,3)'*p_129) + 1;
E_3_234 = abs(g_3_234)^2*T_3_234 - 2*real(g_3_234*H(:,3)'*p_234) + 1;
E_3_235 = abs(g_3_235)^2*T_3_235 - 2*real(g_3_235*H(:,3)'*p_235) + 1;
E_3_236 = abs(g_3_236)^2*T_3_236 - 2*real(g_3_236*H(:,3)'*p_236) + 1;
E_3_237 = abs(g_3_237)^2*T_3_237 - 2*real(g_3_237*H(:,3)'*p_237) + 1;
E_3_238 = abs(g_3_238)^2*T_3_238 - 2*real(g_3_238*H(:,3)'*p_238) + 1;
E_3_239 = abs(g_3_239)^2*T_3_239 - 2*real(g_3_239*H(:,3)'*p_239) + 1;
E_3_245 = abs(g_3_245)^2*T_3_245 - 2*real(g_3_245*H(:,3)'*p_245) + 1;
E_3_246 = abs(g_3_246)^2*T_3_246 - 2*real(g_3_246*H(:,3)'*p_246) + 1;
E_3_247 = abs(g_3_247)^2*T_3_247 - 2*real(g_3_247*H(:,3)'*p_247) + 1;
E_3_248 = abs(g_3_248)^2*T_3_248 - 2*real(g_3_248*H(:,3)'*p_248) + 1;
E_3_249 = abs(g_3_249)^2*T_3_249 - 2*real(g_3_249*H(:,3)'*p_249) + 1;
E_3_256 = abs(g_3_256)^2*T_3_256 - 2*real(g_3_256*H(:,3)'*p_256) + 1;
E_3_257 = abs(g_3_257)^2*T_3_257 - 2*real(g_3_257*H(:,3)'*p_257) + 1;
E_3_258 = abs(g_3_258)^2*T_3_258 - 2*real(g_3_258*H(:,3)'*p_258) + 1;
E_3_259 = abs(g_3_259)^2*T_3_259 - 2*real(g_3_259*H(:,3)'*p_259) + 1;
E_3_267 = abs(g_3_267)^2*T_3_267 - 2*real(g_3_267*H(:,3)'*p_267) + 1;
E_3_268 = abs(g_3_268)^2*T_3_268 - 2*real(g_3_268*H(:,3)'*p_268) + 1;
E_3_269 = abs(g_3_269)^2*T_3_269 - 2*real(g_3_269*H(:,3)'*p_269) + 1;
E_3_278 = abs(g_3_278)^2*T_3_278 - 2*real(g_3_278*H(:,3)'*p_278) + 1;
E_3_279 = abs(g_3_279)^2*T_3_279 - 2*real(g_3_279*H(:,3)'*p_279) + 1;
E_3_289 = abs(g_3_289)^2*T_3_289 - 2*real(g_3_289*H(:,3)'*p_289) + 1;
E_3 = abs(g_3)^2*T_3 - 2*real(g_3*H(:,3)'*p_2) + 1;

E_4_123456789 = abs(g_4_123456789)^2*T_4_123456789 - 2*real(g_4_123456789*H(:,4)'*p_123456789) + 1;
E_4_123 = abs(g_4_123)^2*T_4_123 - 2*real(g_4_123*H(:,4)'*p_123) + 1;
E_4_124 = abs(g_4_124)^2*T_4_124 - 2*real(g_4_124*H(:,4)'*p_124) + 1;
E_4_125 = abs(g_4_125)^2*T_4_125 - 2*real(g_4_125*H(:,4)'*p_125) + 1;
E_4_126 = abs(g_4_126)^2*T_4_126 - 2*real(g_4_126*H(:,4)'*p_126) + 1;
E_4_127 = abs(g_4_127)^2*T_4_127 - 2*real(g_4_127*H(:,4)'*p_127) + 1;
E_4_128 = abs(g_4_128)^2*T_4_128 - 2*real(g_4_128*H(:,4)'*p_128) + 1;
E_4_129 = abs(g_4_129)^2*T_4_129 - 2*real(g_4_129*H(:,4)'*p_129) + 1;
E_4_234 = abs(g_4_234)^2*T_4_234 - 2*real(g_4_234*H(:,4)'*p_234) + 1;
E_4_235 = abs(g_4_235)^2*T_4_235 - 2*real(g_4_235*H(:,4)'*p_235) + 1;
E_4_236 = abs(g_4_236)^2*T_4_236 - 2*real(g_4_236*H(:,4)'*p_236) + 1;
E_4_237 = abs(g_4_237)^2*T_4_237 - 2*real(g_4_237*H(:,4)'*p_237) + 1;
E_4_238 = abs(g_4_238)^2*T_4_238 - 2*real(g_4_238*H(:,4)'*p_238) + 1;
E_4_239 = abs(g_4_239)^2*T_4_239 - 2*real(g_4_239*H(:,4)'*p_239) + 1;
E_4_245 = abs(g_4_245)^2*T_4_245 - 2*real(g_4_245*H(:,4)'*p_245) + 1;
E_4_246 = abs(g_4_246)^2*T_4_246 - 2*real(g_4_246*H(:,4)'*p_246) + 1;
E_4_247 = abs(g_4_247)^2*T_4_247 - 2*real(g_4_247*H(:,4)'*p_247) + 1;
E_4_248 = abs(g_4_248)^2*T_4_248 - 2*real(g_4_248*H(:,4)'*p_248) + 1;
E_4_249 = abs(g_4_249)^2*T_4_249 - 2*real(g_4_249*H(:,4)'*p_249) + 1;
E_4_256 = abs(g_4_256)^2*T_4_256 - 2*real(g_4_256*H(:,4)'*p_256) + 1;
E_4_257 = abs(g_4_257)^2*T_4_257 - 2*real(g_4_257*H(:,4)'*p_257) + 1;
E_4_258 = abs(g_4_258)^2*T_4_258 - 2*real(g_4_258*H(:,4)'*p_258) + 1;
E_4_259 = abs(g_4_259)^2*T_4_259 - 2*real(g_4_259*H(:,4)'*p_259) + 1;
E_4_267 = abs(g_4_267)^2*T_4_267 - 2*real(g_4_267*H(:,4)'*p_267) + 1;
E_4_268 = abs(g_4_268)^2*T_4_268 - 2*real(g_4_268*H(:,4)'*p_268) + 1;
E_4_269 = abs(g_4_269)^2*T_4_269 - 2*real(g_4_269*H(:,4)'*p_269) + 1;
E_4_278 = abs(g_4_278)^2*T_4_278 - 2*real(g_4_278*H(:,4)'*p_278) + 1;
E_4_279 = abs(g_4_279)^2*T_4_279 - 2*real(g_4_279*H(:,4)'*p_279) + 1;
E_4_289 = abs(g_4_289)^2*T_4_289 - 2*real(g_4_289*H(:,4)'*p_289) + 1;
E_4 = abs(g_4)^2*T_4 - 2*real(g_4*H(:,4)'*p_2) + 1;

E_5_123456789 = abs(g_5_123456789)^2*T_5_123456789 - 2*real(g_5_123456789*H(:,5)'*p_123456789) + 1;
E_5_123 = abs(g_5_123)^2*T_5_123 - 2*real(g_5_123*H(:,5)'*p_123) + 1;
E_5_134 = abs(g_5_134)^2*T_5_134 - 2*real(g_5_134*H(:,5)'*p_134) + 1;
E_5_135 = abs(g_5_135)^2*T_5_135 - 2*real(g_5_135*H(:,5)'*p_135) + 1;
E_5_136 = abs(g_5_136)^2*T_5_136 - 2*real(g_5_136*H(:,5)'*p_136) + 1;
E_5_137 = abs(g_5_137)^2*T_5_137 - 2*real(g_5_137*H(:,5)'*p_137) + 1;
E_5_138 = abs(g_5_138)^2*T_5_138 - 2*real(g_5_138*H(:,5)'*p_138) + 1;
E_5_139 = abs(g_5_139)^2*T_5_139 - 2*real(g_5_139*H(:,5)'*p_139) + 1;
E_5_234 = abs(g_5_234)^2*T_5_234 - 2*real(g_5_234*H(:,5)'*p_234) + 1;
E_5_235 = abs(g_5_235)^2*T_5_235 - 2*real(g_5_235*H(:,5)'*p_235) + 1;
E_5_236 = abs(g_5_236)^2*T_5_236 - 2*real(g_5_236*H(:,5)'*p_236) + 1;
E_5_237 = abs(g_5_237)^2*T_5_237 - 2*real(g_5_237*H(:,5)'*p_237) + 1;
E_5_238 = abs(g_5_238)^2*T_5_238 - 2*real(g_5_238*H(:,5)'*p_238) + 1;
E_5_239 = abs(g_5_239)^2*T_5_239 - 2*real(g_5_239*H(:,5)'*p_239) + 1;
E_5_345 = abs(g_5_345)^2*T_5_345 - 2*real(g_5_345*H(:,5)'*p_345) + 1;
E_5_346 = abs(g_5_346)^2*T_5_346 - 2*real(g_5_346*H(:,5)'*p_346) + 1;
E_5_347 = abs(g_5_347)^2*T_5_347 - 2*real(g_5_347*H(:,5)'*p_347) + 1;
E_5_348 = abs(g_5_348)^2*T_5_348 - 2*real(g_5_348*H(:,5)'*p_348) + 1;
E_5_349 = abs(g_5_349)^2*T_5_349 - 2*real(g_5_349*H(:,5)'*p_349) + 1;
E_5_356 = abs(g_5_356)^2*T_5_356 - 2*real(g_5_356*H(:,5)'*p_356) + 1;
E_5_357 = abs(g_5_357)^2*T_5_357 - 2*real(g_5_357*H(:,5)'*p_357) + 1;
E_5_358 = abs(g_5_358)^2*T_5_358 - 2*real(g_5_358*H(:,5)'*p_358) + 1;
E_5_359 = abs(g_5_359)^2*T_5_359 - 2*real(g_5_359*H(:,5)'*p_359) + 1;
E_5_367 = abs(g_5_367)^2*T_5_367 - 2*real(g_5_367*H(:,5)'*p_367) + 1;
E_5_368 = abs(g_5_368)^2*T_5_368 - 2*real(g_5_368*H(:,5)'*p_368) + 1;
E_5_369 = abs(g_5_369)^2*T_5_369 - 2*real(g_5_369*H(:,5)'*p_369) + 1;
E_5_378 = abs(g_5_378)^2*T_5_378 - 2*real(g_5_378*H(:,5)'*p_378) + 1;
E_5_379 = abs(g_5_379)^2*T_5_379 - 2*real(g_5_379*H(:,5)'*p_379) + 1;
E_5_389 = abs(g_5_389)^2*T_5_389 - 2*real(g_5_389*H(:,5)'*p_389) + 1;
E_5 = abs(g_5)^2*T_5 - 2*real(g_5*H(:,5)'*p_3) + 1;

E_6_123456789 = abs(g_6_123456789)^2*T_6_123456789 - 2*real(g_6_123456789*H(:,6)'*p_123456789) + 1;
E_6_123 = abs(g_6_123)^2*T_6_123 - 2*real(g_6_123*H(:,6)'*p_123) + 1;
E_6_134 = abs(g_6_134)^2*T_6_134 - 2*real(g_6_134*H(:,6)'*p_134) + 1;
E_6_135 = abs(g_6_135)^2*T_6_135 - 2*real(g_6_135*H(:,6)'*p_135) + 1;
E_6_136 = abs(g_6_136)^2*T_6_136 - 2*real(g_6_136*H(:,6)'*p_136) + 1;
E_6_137 = abs(g_6_137)^2*T_6_137 - 2*real(g_6_137*H(:,6)'*p_137) + 1;
E_6_138 = abs(g_6_138)^2*T_6_138 - 2*real(g_6_138*H(:,6)'*p_138) + 1;
E_6_139 = abs(g_6_139)^2*T_6_139 - 2*real(g_6_139*H(:,6)'*p_139) + 1;
E_6_234 = abs(g_6_234)^2*T_6_234 - 2*real(g_6_234*H(:,6)'*p_234) + 1;
E_6_235 = abs(g_6_235)^2*T_6_235 - 2*real(g_6_235*H(:,6)'*p_235) + 1;
E_6_236 = abs(g_6_236)^2*T_6_236 - 2*real(g_6_236*H(:,6)'*p_236) + 1;
E_6_237 = abs(g_6_237)^2*T_6_237 - 2*real(g_6_237*H(:,6)'*p_237) + 1;
E_6_238 = abs(g_6_238)^2*T_6_238 - 2*real(g_6_238*H(:,6)'*p_238) + 1;
E_6_239 = abs(g_6_239)^2*T_6_239 - 2*real(g_6_239*H(:,6)'*p_239) + 1;
E_6_345 = abs(g_6_345)^2*T_6_345 - 2*real(g_6_345*H(:,6)'*p_345) + 1;
E_6_346 = abs(g_6_346)^2*T_6_346 - 2*real(g_6_346*H(:,6)'*p_346) + 1;
E_6_347 = abs(g_6_347)^2*T_6_347 - 2*real(g_6_347*H(:,6)'*p_347) + 1;
E_6_348 = abs(g_6_348)^2*T_6_348 - 2*real(g_6_348*H(:,6)'*p_348) + 1;
E_6_349 = abs(g_6_349)^2*T_6_349 - 2*real(g_6_349*H(:,6)'*p_349) + 1;
E_6_356 = abs(g_6_356)^2*T_6_356 - 2*real(g_6_356*H(:,6)'*p_356) + 1;
E_6_357 = abs(g_6_357)^2*T_6_357 - 2*real(g_6_357*H(:,6)'*p_357) + 1;
E_6_358 = abs(g_6_358)^2*T_6_358 - 2*real(g_6_358*H(:,6)'*p_358) + 1;
E_6_359 = abs(g_6_359)^2*T_6_359 - 2*real(g_6_359*H(:,6)'*p_359) + 1;
E_6_367 = abs(g_6_367)^2*T_6_367 - 2*real(g_6_367*H(:,6)'*p_367) + 1;
E_6_368 = abs(g_6_368)^2*T_6_368 - 2*real(g_6_368*H(:,6)'*p_368) + 1;
E_6_369 = abs(g_6_369)^2*T_6_369 - 2*real(g_6_369*H(:,6)'*p_369) + 1;
E_6_378 = abs(g_6_378)^2*T_6_378 - 2*real(g_6_378*H(:,6)'*p_378) + 1;
E_6_379 = abs(g_6_379)^2*T_6_379 - 2*real(g_6_379*H(:,6)'*p_379) + 1;
E_6_389 = abs(g_6_389)^2*T_6_389 - 2*real(g_6_389*H(:,6)'*p_389) + 1;
E_6 = abs(g_6)^2*T_6 - 2*real(g_6*H(:,6)'*p_3) + 1;

E_7_123456789 = abs(g_7_123456789)^2*T_7_123456789 - 2*real(g_7_123456789*H(:,7)'*p_123456789) + 1;
E_7_124 = abs(g_7_124)^2*T_7_124 - 2*real(g_7_124*H(:,7)'*p_124) + 1;
E_7_134 = abs(g_7_134)^2*T_7_134 - 2*real(g_7_134*H(:,7)'*p_134) + 1;
E_7_145 = abs(g_7_145)^2*T_7_145 - 2*real(g_7_145*H(:,7)'*p_145) + 1;
E_7_146 = abs(g_7_146)^2*T_7_146 - 2*real(g_7_146*H(:,7)'*p_146) + 1;
E_7_147 = abs(g_7_147)^2*T_7_147 - 2*real(g_7_147*H(:,7)'*p_147) + 1;
E_7_148 = abs(g_7_148)^2*T_7_148 - 2*real(g_7_148*H(:,7)'*p_148) + 1;
E_7_149 = abs(g_7_149)^2*T_7_149 - 2*real(g_7_149*H(:,7)'*p_149) + 1;
E_7_234 = abs(g_7_234)^2*T_7_234 - 2*real(g_7_234*H(:,7)'*p_234) + 1;
E_7_245 = abs(g_7_245)^2*T_7_245 - 2*real(g_7_245*H(:,7)'*p_245) + 1;
E_7_246 = abs(g_7_246)^2*T_7_246 - 2*real(g_7_246*H(:,7)'*p_246) + 1;
E_7_247 = abs(g_7_247)^2*T_7_247 - 2*real(g_7_247*H(:,7)'*p_247) + 1;
E_7_248 = abs(g_7_248)^2*T_7_248 - 2*real(g_7_248*H(:,7)'*p_248) + 1;
E_7_249 = abs(g_7_249)^2*T_7_249 - 2*real(g_7_249*H(:,7)'*p_249) + 1;
E_7_345 = abs(g_7_345)^2*T_7_345 - 2*real(g_7_345*H(:,7)'*p_345) + 1;
E_7_346 = abs(g_7_346)^2*T_7_346 - 2*real(g_7_346*H(:,7)'*p_346) + 1;
E_7_347 = abs(g_7_347)^2*T_7_347 - 2*real(g_7_347*H(:,7)'*p_347) + 1;
E_7_348 = abs(g_7_348)^2*T_7_348 - 2*real(g_7_348*H(:,7)'*p_348) + 1;
E_7_349 = abs(g_7_349)^2*T_7_349 - 2*real(g_7_349*H(:,7)'*p_349) + 1;
E_7_456 = abs(g_7_456)^2*T_7_456 - 2*real(g_7_456*H(:,7)'*p_456) + 1;
E_7_457 = abs(g_7_457)^2*T_7_457 - 2*real(g_7_457*H(:,7)'*p_457) + 1;
E_7_458 = abs(g_7_458)^2*T_7_458 - 2*real(g_7_458*H(:,7)'*p_458) + 1;
E_7_459 = abs(g_7_459)^2*T_7_459 - 2*real(g_7_459*H(:,7)'*p_459) + 1;
E_7_467 = abs(g_7_467)^2*T_7_467 - 2*real(g_7_467*H(:,7)'*p_467) + 1;
E_7_468 = abs(g_7_468)^2*T_7_468 - 2*real(g_7_468*H(:,7)'*p_468) + 1;
E_7_469 = abs(g_7_469)^2*T_7_469 - 2*real(g_7_469*H(:,7)'*p_469) + 1;
E_7_478 = abs(g_7_478)^2*T_7_478 - 2*real(g_7_478*H(:,7)'*p_478) + 1;
E_7_479 = abs(g_7_479)^2*T_7_479 - 2*real(g_7_479*H(:,7)'*p_479) + 1;
E_7_489 = abs(g_7_489)^2*T_7_489 - 2*real(g_7_489*H(:,7)'*p_489) + 1;
E_7 = abs(g_7)^2*T_7 - 2*real(g_7*H(:,7)'*p_4) + 1;

E_8_123456789 = abs(g_8_123456789)^2*T_8_123456789 - 2*real(g_8_123456789*H(:,8)'*p_123456789) + 1;
E_8_124 = abs(g_8_124)^2*T_8_124 - 2*real(g_8_124*H(:,8)'*p_124) + 1;
E_8_134 = abs(g_8_134)^2*T_8_134 - 2*real(g_8_134*H(:,8)'*p_134) + 1;
E_8_145 = abs(g_8_145)^2*T_8_145 - 2*real(g_8_145*H(:,8)'*p_145) + 1;
E_8_146 = abs(g_8_146)^2*T_8_146 - 2*real(g_8_146*H(:,8)'*p_146) + 1;
E_8_147 = abs(g_8_147)^2*T_8_147 - 2*real(g_8_147*H(:,8)'*p_147) + 1;
E_8_148 = abs(g_8_148)^2*T_8_148 - 2*real(g_8_148*H(:,8)'*p_148) + 1;
E_8_149 = abs(g_8_149)^2*T_8_149 - 2*real(g_8_149*H(:,8)'*p_149) + 1;
E_8_234 = abs(g_8_234)^2*T_8_234 - 2*real(g_8_234*H(:,8)'*p_234) + 1;
E_8_245 = abs(g_8_245)^2*T_8_245 - 2*real(g_8_245*H(:,8)'*p_245) + 1;
E_8_246 = abs(g_8_246)^2*T_8_246 - 2*real(g_8_246*H(:,8)'*p_246) + 1;
E_8_247 = abs(g_8_247)^2*T_8_247 - 2*real(g_8_247*H(:,8)'*p_247) + 1;
E_8_248 = abs(g_8_248)^2*T_8_248 - 2*real(g_8_248*H(:,8)'*p_248) + 1;
E_8_249 = abs(g_8_249)^2*T_8_249 - 2*real(g_8_249*H(:,8)'*p_249) + 1;
E_8_345 = abs(g_8_345)^2*T_8_345 - 2*real(g_8_345*H(:,8)'*p_345) + 1;
E_8_346 = abs(g_8_346)^2*T_8_346 - 2*real(g_8_346*H(:,8)'*p_346) + 1;
E_8_347 = abs(g_8_347)^2*T_8_347 - 2*real(g_8_347*H(:,8)'*p_347) + 1;
E_8_348 = abs(g_8_348)^2*T_8_348 - 2*real(g_8_348*H(:,8)'*p_348) + 1;
E_8_349 = abs(g_8_349)^2*T_8_349 - 2*real(g_8_349*H(:,8)'*p_349) + 1;
E_8_456 = abs(g_8_456)^2*T_8_456 - 2*real(g_8_456*H(:,8)'*p_456) + 1;
E_8_457 = abs(g_8_457)^2*T_8_457 - 2*real(g_8_457*H(:,8)'*p_457) + 1;
E_8_458 = abs(g_8_458)^2*T_8_458 - 2*real(g_8_458*H(:,8)'*p_458) + 1;
E_8_459 = abs(g_8_459)^2*T_8_459 - 2*real(g_8_459*H(:,8)'*p_459) + 1;
E_8_467 = abs(g_8_467)^2*T_8_467 - 2*real(g_8_467*H(:,8)'*p_467) + 1;
E_8_468 = abs(g_8_468)^2*T_8_468 - 2*real(g_8_468*H(:,8)'*p_468) + 1;
E_8_469 = abs(g_8_469)^2*T_8_469 - 2*real(g_8_469*H(:,8)'*p_469) + 1;
E_8_478 = abs(g_8_478)^2*T_8_478 - 2*real(g_8_478*H(:,8)'*p_478) + 1;
E_8_479 = abs(g_8_479)^2*T_8_479 - 2*real(g_8_479*H(:,8)'*p_479) + 1;
E_8_489 = abs(g_8_489)^2*T_8_489 - 2*real(g_8_489*H(:,8)'*p_489) + 1;
E_8 = abs(g_8)^2*T_8 - 2*real(g_8*H(:,8)'*p_4) + 1;

E_9_123456789 = abs(g_9_123456789)^2*T_9_123456789 - 2*real(g_9_123456789*H(:,9)'*p_123456789) + 1;
E_9_125 = abs(g_9_125)^2*T_9_125 - 2*real(g_9_125*H(:,9)'*p_125) + 1;
E_9_135 = abs(g_9_135)^2*T_9_135 - 2*real(g_9_135*H(:,9)'*p_135) + 1;
E_9_145 = abs(g_9_145)^2*T_9_145 - 2*real(g_9_145*H(:,9)'*p_145) + 1;
E_9_156 = abs(g_9_156)^2*T_9_156 - 2*real(g_9_156*H(:,9)'*p_156) + 1;
E_9_157 = abs(g_9_157)^2*T_9_157 - 2*real(g_9_157*H(:,9)'*p_157) + 1;
E_9_158 = abs(g_9_158)^2*T_9_158 - 2*real(g_9_158*H(:,9)'*p_158) + 1;
E_9_159 = abs(g_9_159)^2*T_9_159 - 2*real(g_9_159*H(:,9)'*p_159) + 1;
E_9_235 = abs(g_9_235)^2*T_9_235 - 2*real(g_9_235*H(:,9)'*p_235) + 1;
E_9_245 = abs(g_9_245)^2*T_9_245 - 2*real(g_9_245*H(:,9)'*p_245) + 1;
E_9_256 = abs(g_9_256)^2*T_9_256 - 2*real(g_9_256*H(:,9)'*p_256) + 1;
E_9_257 = abs(g_9_257)^2*T_9_257 - 2*real(g_9_257*H(:,9)'*p_257) + 1;
E_9_258 = abs(g_9_258)^2*T_9_258 - 2*real(g_9_258*H(:,9)'*p_258) + 1;
E_9_259 = abs(g_9_259)^2*T_9_259 - 2*real(g_9_259*H(:,9)'*p_259) + 1;
E_9_345 = abs(g_9_345)^2*T_9_345 - 2*real(g_9_345*H(:,9)'*p_345) + 1;
E_9_356 = abs(g_9_356)^2*T_9_356 - 2*real(g_9_356*H(:,9)'*p_356) + 1;
E_9_357 = abs(g_9_357)^2*T_9_357 - 2*real(g_9_357*H(:,9)'*p_357) + 1;
E_9_358 = abs(g_9_358)^2*T_9_358 - 2*real(g_9_358*H(:,9)'*p_358) + 1;
E_9_359 = abs(g_9_359)^2*T_9_359 - 2*real(g_9_359*H(:,9)'*p_359) + 1;
E_9_456 = abs(g_9_456)^2*T_9_456 - 2*real(g_9_456*H(:,9)'*p_456) + 1;
E_9_457 = abs(g_9_457)^2*T_9_457 - 2*real(g_9_457*H(:,9)'*p_457) + 1;
E_9_458 = abs(g_9_458)^2*T_9_458 - 2*real(g_9_458*H(:,9)'*p_458) + 1;
E_9_459 = abs(g_9_459)^2*T_9_459 - 2*real(g_9_459*H(:,9)'*p_459) + 1;
E_9_567 = abs(g_9_567)^2*T_9_567 - 2*real(g_9_567*H(:,9)'*p_567) + 1;
E_9_568 = abs(g_9_568)^2*T_9_568 - 2*real(g_9_568*H(:,9)'*p_568) + 1;
E_9_569 = abs(g_9_569)^2*T_9_569 - 2*real(g_9_569*H(:,9)'*p_569) + 1;
E_9_578 = abs(g_9_578)^2*T_9_578 - 2*real(g_9_578*H(:,9)'*p_578) + 1;
E_9_579 = abs(g_9_579)^2*T_9_579 - 2*real(g_9_579*H(:,9)'*p_579) + 1;
E_9_589 = abs(g_9_589)^2*T_9_589 - 2*real(g_9_589*H(:,9)'*p_589) + 1;
E_9 = abs(g_9)^2*T_9 - 2*real(g_9*H(:,9)'*p_5) + 1;

E_10_123456789 = abs(g_10_123456789)^2*T_10_123456789 - 2*real(g_10_123456789*H(:,10)'*p_123456789) + 1;
E_10_125 = abs(g_10_125)^2*T_10_125 - 2*real(g_10_125*H(:,10)'*p_125) + 1;
E_10_135 = abs(g_10_135)^2*T_10_135 - 2*real(g_10_135*H(:,10)'*p_135) + 1;
E_10_145 = abs(g_10_145)^2*T_10_145 - 2*real(g_10_145*H(:,10)'*p_145) + 1;
E_10_156 = abs(g_10_156)^2*T_10_156 - 2*real(g_10_156*H(:,10)'*p_156) + 1;
E_10_157 = abs(g_10_157)^2*T_10_157 - 2*real(g_10_157*H(:,10)'*p_157) + 1;
E_10_158 = abs(g_10_158)^2*T_10_158 - 2*real(g_10_158*H(:,10)'*p_158) + 1;
E_10_159 = abs(g_10_159)^2*T_10_159 - 2*real(g_10_159*H(:,10)'*p_159) + 1;
E_10_235 = abs(g_10_235)^2*T_10_235 - 2*real(g_10_235*H(:,10)'*p_235) + 1;
E_10_245 = abs(g_10_245)^2*T_10_245 - 2*real(g_10_245*H(:,10)'*p_245) + 1;
E_10_256 = abs(g_10_256)^2*T_10_256 - 2*real(g_10_256*H(:,10)'*p_256) + 1;
E_10_257 = abs(g_10_257)^2*T_10_257 - 2*real(g_10_257*H(:,10)'*p_257) + 1;
E_10_258 = abs(g_10_258)^2*T_10_258 - 2*real(g_10_258*H(:,10)'*p_258) + 1;
E_10_259 = abs(g_10_259)^2*T_10_259 - 2*real(g_10_259*H(:,10)'*p_259) + 1;
E_10_345 = abs(g_10_345)^2*T_10_345 - 2*real(g_10_345*H(:,10)'*p_345) + 1;
E_10_356 = abs(g_10_356)^2*T_10_356 - 2*real(g_10_356*H(:,10)'*p_356) + 1;
E_10_357 = abs(g_10_357)^2*T_10_357 - 2*real(g_10_357*H(:,10)'*p_357) + 1;
E_10_358 = abs(g_10_358)^2*T_10_358 - 2*real(g_10_358*H(:,10)'*p_358) + 1;
E_10_359 = abs(g_10_359)^2*T_10_359 - 2*real(g_10_359*H(:,10)'*p_359) + 1;
E_10_456 = abs(g_10_456)^2*T_10_456 - 2*real(g_10_456*H(:,10)'*p_456) + 1;
E_10_457 = abs(g_10_457)^2*T_10_457 - 2*real(g_10_457*H(:,10)'*p_457) + 1;
E_10_458 = abs(g_10_458)^2*T_10_458 - 2*real(g_10_458*H(:,10)'*p_458) + 1;
E_10_459 = abs(g_10_459)^2*T_10_459 - 2*real(g_10_459*H(:,10)'*p_459) + 1;
E_10_567 = abs(g_10_567)^2*T_10_567 - 2*real(g_10_567*H(:,10)'*p_567) + 1;
E_10_568 = abs(g_10_568)^2*T_10_568 - 2*real(g_10_568*H(:,10)'*p_568) + 1;
E_10_569 = abs(g_10_569)^2*T_10_569 - 2*real(g_10_569*H(:,10)'*p_569) + 1;
E_10_578 = abs(g_10_578)^2*T_10_578 - 2*real(g_10_578*H(:,10)'*p_578) + 1;
E_10_579 = abs(g_10_579)^2*T_10_579 - 2*real(g_10_579*H(:,10)'*p_579) + 1;
E_10_589 = abs(g_10_589)^2*T_10_589 - 2*real(g_10_589*H(:,10)'*p_589) + 1;
E_10 = abs(g_10)^2*T_10 - 2*real(g_10*H(:,10)'*p_5) + 1;

E_11_123456789 = abs(g_11_123456789)^2*T_11_123456789 - 2*real(g_11_123456789*H(:,11)'*p_123456789) + 1;
E_11_126 = abs(g_11_126)^2*T_11_126 - 2*real(g_11_126*H(:,11)'*p_126) + 1;
E_11_136 = abs(g_11_136)^2*T_11_136 - 2*real(g_11_136*H(:,11)'*p_136) + 1;
E_11_146 = abs(g_11_146)^2*T_11_146 - 2*real(g_11_146*H(:,11)'*p_146) + 1;
E_11_156 = abs(g_11_156)^2*T_11_156 - 2*real(g_11_156*H(:,11)'*p_156) + 1;
E_11_167 = abs(g_11_167)^2*T_11_167 - 2*real(g_11_167*H(:,11)'*p_167) + 1;
E_11_168 = abs(g_11_168)^2*T_11_168 - 2*real(g_11_168*H(:,11)'*p_168) + 1;
E_11_169 = abs(g_11_169)^2*T_11_169 - 2*real(g_11_169*H(:,11)'*p_169) + 1;
E_11_236 = abs(g_11_236)^2*T_11_236 - 2*real(g_11_236*H(:,11)'*p_236) + 1;
E_11_246 = abs(g_11_246)^2*T_11_246 - 2*real(g_11_246*H(:,11)'*p_246) + 1;
E_11_256 = abs(g_11_256)^2*T_11_256 - 2*real(g_11_256*H(:,11)'*p_256) + 1;
E_11_267 = abs(g_11_267)^2*T_11_267 - 2*real(g_11_267*H(:,11)'*p_267) + 1;
E_11_268 = abs(g_11_268)^2*T_11_268 - 2*real(g_11_268*H(:,11)'*p_268) + 1;
E_11_269 = abs(g_11_269)^2*T_11_269 - 2*real(g_11_269*H(:,11)'*p_269) + 1;
E_11_346 = abs(g_11_346)^2*T_11_346 - 2*real(g_11_346*H(:,11)'*p_346) + 1;
E_11_356 = abs(g_11_356)^2*T_11_356 - 2*real(g_11_356*H(:,11)'*p_356) + 1;
E_11_367 = abs(g_11_367)^2*T_11_367 - 2*real(g_11_367*H(:,11)'*p_367) + 1;
E_11_368 = abs(g_11_368)^2*T_11_368 - 2*real(g_11_368*H(:,11)'*p_368) + 1;
E_11_369 = abs(g_11_369)^2*T_11_369 - 2*real(g_11_369*H(:,11)'*p_369) + 1;
E_11_456 = abs(g_11_456)^2*T_11_456 - 2*real(g_11_456*H(:,11)'*p_456) + 1;
E_11_467 = abs(g_11_467)^2*T_11_467 - 2*real(g_11_467*H(:,11)'*p_467) + 1;
E_11_468 = abs(g_11_468)^2*T_11_468 - 2*real(g_11_468*H(:,11)'*p_468) + 1;
E_11_469 = abs(g_11_469)^2*T_11_469 - 2*real(g_11_469*H(:,11)'*p_469) + 1;
E_11_567 = abs(g_11_567)^2*T_11_567 - 2*real(g_11_567*H(:,11)'*p_567) + 1;
E_11_568 = abs(g_11_568)^2*T_11_568 - 2*real(g_11_568*H(:,11)'*p_568) + 1;
E_11_569 = abs(g_11_569)^2*T_11_569 - 2*real(g_11_569*H(:,11)'*p_569) + 1;
E_11_678 = abs(g_11_678)^2*T_11_678 - 2*real(g_11_678*H(:,11)'*p_678) + 1;
E_11_679 = abs(g_11_679)^2*T_11_679 - 2*real(g_11_679*H(:,11)'*p_679) + 1;
E_11_689 = abs(g_11_689)^2*T_11_689 - 2*real(g_11_689*H(:,11)'*p_689) + 1;
E_11 = abs(g_11)^2*T_11 - 2*real(g_11*H(:,11)'*p_6) + 1;

E_12_123456789 = abs(g_12_123456789)^2*T_12_123456789 - 2*real(g_12_123456789*H(:,12)'*p_123456789) + 1;
E_12_126 = abs(g_12_126)^2*T_12_126 - 2*real(g_12_126*H(:,12)'*p_126) + 1;
E_12_136 = abs(g_12_136)^2*T_12_136 - 2*real(g_12_136*H(:,12)'*p_136) + 1;
E_12_146 = abs(g_12_146)^2*T_12_146 - 2*real(g_12_146*H(:,12)'*p_146) + 1;
E_12_156 = abs(g_12_156)^2*T_12_156 - 2*real(g_12_156*H(:,12)'*p_156) + 1;
E_12_167 = abs(g_12_167)^2*T_12_167 - 2*real(g_12_167*H(:,12)'*p_167) + 1;
E_12_168 = abs(g_12_168)^2*T_12_168 - 2*real(g_12_168*H(:,12)'*p_168) + 1;
E_12_169 = abs(g_12_169)^2*T_12_169 - 2*real(g_12_169*H(:,12)'*p_169) + 1;
E_12_236 = abs(g_12_236)^2*T_12_236 - 2*real(g_12_236*H(:,12)'*p_236) + 1;
E_12_246 = abs(g_12_246)^2*T_12_246 - 2*real(g_12_246*H(:,12)'*p_246) + 1;
E_12_256 = abs(g_12_256)^2*T_12_256 - 2*real(g_12_256*H(:,12)'*p_256) + 1;
E_12_267 = abs(g_12_267)^2*T_12_267 - 2*real(g_12_267*H(:,12)'*p_267) + 1;
E_12_268 = abs(g_12_268)^2*T_12_268 - 2*real(g_12_268*H(:,12)'*p_268) + 1;
E_12_269 = abs(g_12_269)^2*T_12_269 - 2*real(g_12_269*H(:,12)'*p_269) + 1;
E_12_346 = abs(g_12_346)^2*T_12_346 - 2*real(g_12_346*H(:,12)'*p_346) + 1;
E_12_356 = abs(g_12_356)^2*T_12_356 - 2*real(g_12_356*H(:,12)'*p_356) + 1;
E_12_367 = abs(g_12_367)^2*T_12_367 - 2*real(g_12_367*H(:,12)'*p_367) + 1;
E_12_368 = abs(g_12_368)^2*T_12_368 - 2*real(g_12_368*H(:,12)'*p_368) + 1;
E_12_369 = abs(g_12_369)^2*T_12_369 - 2*real(g_12_369*H(:,12)'*p_369) + 1;
E_12_456 = abs(g_12_456)^2*T_12_456 - 2*real(g_12_456*H(:,12)'*p_456) + 1;
E_12_467 = abs(g_12_467)^2*T_12_467 - 2*real(g_12_467*H(:,12)'*p_467) + 1;
E_12_468 = abs(g_12_468)^2*T_12_468 - 2*real(g_12_468*H(:,12)'*p_468) + 1;
E_12_469 = abs(g_12_469)^2*T_12_469 - 2*real(g_12_469*H(:,12)'*p_469) + 1;
E_12_567 = abs(g_12_567)^2*T_12_567 - 2*real(g_12_567*H(:,12)'*p_567) + 1;
E_12_568 = abs(g_12_568)^2*T_12_568 - 2*real(g_12_568*H(:,12)'*p_568) + 1;
E_12_569 = abs(g_12_569)^2*T_12_569 - 2*real(g_12_569*H(:,12)'*p_569) + 1;
E_12_678 = abs(g_12_678)^2*T_12_678 - 2*real(g_12_678*H(:,12)'*p_678) + 1;
E_12_679 = abs(g_12_679)^2*T_12_679 - 2*real(g_12_679*H(:,12)'*p_679) + 1;
E_12_689 = abs(g_12_689)^2*T_12_689 - 2*real(g_12_689*H(:,12)'*p_689) + 1;
E_12 = abs(g_12)^2*T_12 - 2*real(g_12*H(:,12)'*p_6) + 1;

E_13_123456789 = abs(g_13_123456789)^2*T_13_123456789 - 2*real(g_13_123456789*H(:,13)'*p_123456789) + 1;
E_13_127 = abs(g_13_127)^2*T_13_127 - 2*real(g_13_127*H(:,13)'*p_127) + 1;
E_13_137 = abs(g_13_137)^2*T_13_137 - 2*real(g_13_137*H(:,13)'*p_137) + 1;
E_13_147 = abs(g_13_147)^2*T_13_147 - 2*real(g_13_147*H(:,13)'*p_147) + 1;
E_13_157 = abs(g_13_157)^2*T_13_157 - 2*real(g_13_157*H(:,13)'*p_157) + 1;
E_13_167 = abs(g_13_167)^2*T_13_167 - 2*real(g_13_167*H(:,13)'*p_167) + 1;
E_13_178 = abs(g_13_178)^2*T_13_178 - 2*real(g_13_178*H(:,13)'*p_178) + 1;
E_13_179 = abs(g_13_179)^2*T_13_179 - 2*real(g_13_179*H(:,13)'*p_179) + 1;
E_13_237 = abs(g_13_237)^2*T_13_237 - 2*real(g_13_237*H(:,13)'*p_237) + 1;
E_13_247 = abs(g_13_247)^2*T_13_247 - 2*real(g_13_247*H(:,13)'*p_247) + 1;
E_13_257 = abs(g_13_257)^2*T_13_257 - 2*real(g_13_257*H(:,13)'*p_257) + 1;
E_13_267 = abs(g_13_267)^2*T_13_267 - 2*real(g_13_267*H(:,13)'*p_267) + 1;
E_13_278 = abs(g_13_278)^2*T_13_278 - 2*real(g_13_278*H(:,13)'*p_278) + 1;
E_13_279 = abs(g_13_279)^2*T_13_279 - 2*real(g_13_279*H(:,13)'*p_279) + 1;
E_13_347 = abs(g_13_347)^2*T_13_347 - 2*real(g_13_347*H(:,13)'*p_347) + 1;
E_13_357 = abs(g_13_357)^2*T_13_357 - 2*real(g_13_357*H(:,13)'*p_357) + 1;
E_13_367 = abs(g_13_367)^2*T_13_367 - 2*real(g_13_367*H(:,13)'*p_367) + 1;
E_13_378 = abs(g_13_378)^2*T_13_378 - 2*real(g_13_378*H(:,13)'*p_378) + 1;
E_13_379 = abs(g_13_379)^2*T_13_379 - 2*real(g_13_379*H(:,13)'*p_379) + 1;
E_13_457 = abs(g_13_457)^2*T_13_457 - 2*real(g_13_457*H(:,13)'*p_457) + 1;
E_13_467 = abs(g_13_467)^2*T_13_467 - 2*real(g_13_467*H(:,13)'*p_467) + 1;
E_13_478 = abs(g_13_478)^2*T_13_478 - 2*real(g_13_478*H(:,13)'*p_478) + 1;
E_13_479 = abs(g_13_479)^2*T_13_479 - 2*real(g_13_479*H(:,13)'*p_479) + 1;
E_13_567 = abs(g_13_567)^2*T_13_567 - 2*real(g_13_567*H(:,13)'*p_567) + 1;
E_13_578 = abs(g_13_578)^2*T_13_578 - 2*real(g_13_578*H(:,13)'*p_578) + 1;
E_13_579 = abs(g_13_579)^2*T_13_579 - 2*real(g_13_579*H(:,13)'*p_579) + 1;
E_13_678 = abs(g_13_678)^2*T_13_678 - 2*real(g_13_678*H(:,13)'*p_678) + 1;
E_13_679 = abs(g_13_679)^2*T_13_679 - 2*real(g_13_679*H(:,13)'*p_679) + 1;
E_13_789 = abs(g_13_789)^2*T_13_789 - 2*real(g_13_789*H(:,13)'*p_789) + 1;
E_13 = abs(g_13)^2*T_13 - 2*real(g_13*H(:,13)'*p_7) + 1;

E_14_123456789 = abs(g_14_123456789)^2*T_14_123456789 - 2*real(g_14_123456789*H(:,14)'*p_123456789) + 1;
E_14_127 = abs(g_14_127)^2*T_14_127 - 2*real(g_14_127*H(:,14)'*p_127) + 1;
E_14_137 = abs(g_14_137)^2*T_14_137 - 2*real(g_14_137*H(:,14)'*p_137) + 1;
E_14_147 = abs(g_14_147)^2*T_14_147 - 2*real(g_14_147*H(:,14)'*p_147) + 1;
E_14_157 = abs(g_14_157)^2*T_14_157 - 2*real(g_14_157*H(:,14)'*p_157) + 1;
E_14_167 = abs(g_14_167)^2*T_14_167 - 2*real(g_14_167*H(:,14)'*p_167) + 1;
E_14_178 = abs(g_14_178)^2*T_14_178 - 2*real(g_14_178*H(:,14)'*p_178) + 1;
E_14_179 = abs(g_14_179)^2*T_14_179 - 2*real(g_14_179*H(:,14)'*p_179) + 1;
E_14_237 = abs(g_14_237)^2*T_14_237 - 2*real(g_14_237*H(:,14)'*p_237) + 1;
E_14_247 = abs(g_14_247)^2*T_14_247 - 2*real(g_14_247*H(:,14)'*p_247) + 1;
E_14_257 = abs(g_14_257)^2*T_14_257 - 2*real(g_14_257*H(:,14)'*p_257) + 1;
E_14_267 = abs(g_14_267)^2*T_14_267 - 2*real(g_14_267*H(:,14)'*p_267) + 1;
E_14_278 = abs(g_14_278)^2*T_14_278 - 2*real(g_14_278*H(:,14)'*p_278) + 1;
E_14_279 = abs(g_14_279)^2*T_14_279 - 2*real(g_14_279*H(:,14)'*p_279) + 1;
E_14_347 = abs(g_14_347)^2*T_14_347 - 2*real(g_14_347*H(:,14)'*p_347) + 1;
E_14_357 = abs(g_14_357)^2*T_14_357 - 2*real(g_14_357*H(:,14)'*p_357) + 1;
E_14_367 = abs(g_14_367)^2*T_14_367 - 2*real(g_14_367*H(:,14)'*p_367) + 1;
E_14_378 = abs(g_14_378)^2*T_14_378 - 2*real(g_14_378*H(:,14)'*p_378) + 1;
E_14_379 = abs(g_14_379)^2*T_14_379 - 2*real(g_14_379*H(:,14)'*p_379) + 1;
E_14_457 = abs(g_14_457)^2*T_14_457 - 2*real(g_14_457*H(:,14)'*p_457) + 1;
E_14_467 = abs(g_14_467)^2*T_14_467 - 2*real(g_14_467*H(:,14)'*p_467) + 1;
E_14_478 = abs(g_14_478)^2*T_14_478 - 2*real(g_14_478*H(:,14)'*p_478) + 1;
E_14_479 = abs(g_14_479)^2*T_14_479 - 2*real(g_14_479*H(:,14)'*p_479) + 1;
E_14_567 = abs(g_14_567)^2*T_14_567 - 2*real(g_14_567*H(:,14)'*p_567) + 1;
E_14_578 = abs(g_14_578)^2*T_14_578 - 2*real(g_14_578*H(:,14)'*p_578) + 1;
E_14_579 = abs(g_14_579)^2*T_14_579 - 2*real(g_14_579*H(:,14)'*p_579) + 1;
E_14_678 = abs(g_14_678)^2*T_14_678 - 2*real(g_14_678*H(:,14)'*p_678) + 1;
E_14_679 = abs(g_14_679)^2*T_14_679 - 2*real(g_14_679*H(:,14)'*p_679) + 1;
E_14_789 = abs(g_14_789)^2*T_14_789 - 2*real(g_14_789*H(:,14)'*p_789) + 1;
E_14 = abs(g_14)^2*T_14 - 2*real(g_14*H(:,14)'*p_7) + 1;

E_15_123456789 = abs(g_15_123456789)^2*T_15_123456789 - 2*real(g_15_123456789*H(:,15)'*p_123456789) + 1;
E_15_128 = abs(g_15_128)^2*T_15_128 - 2*real(g_15_128*H(:,15)'*p_128) + 1;
E_15_138 = abs(g_15_138)^2*T_15_138 - 2*real(g_15_138*H(:,15)'*p_138) + 1;
E_15_148 = abs(g_15_148)^2*T_15_148 - 2*real(g_15_148*H(:,15)'*p_148) + 1;
E_15_158 = abs(g_15_158)^2*T_15_158 - 2*real(g_15_158*H(:,15)'*p_158) + 1;
E_15_168 = abs(g_15_168)^2*T_15_168 - 2*real(g_15_168*H(:,15)'*p_168) + 1;
E_15_178 = abs(g_15_178)^2*T_15_178 - 2*real(g_15_178*H(:,15)'*p_178) + 1;
E_15_189 = abs(g_15_189)^2*T_15_189 - 2*real(g_15_189*H(:,15)'*p_189) + 1;
E_15_238 = abs(g_15_238)^2*T_15_238 - 2*real(g_15_238*H(:,15)'*p_238) + 1;
E_15_248 = abs(g_15_248)^2*T_15_248 - 2*real(g_15_248*H(:,15)'*p_248) + 1;
E_15_258 = abs(g_15_258)^2*T_15_258 - 2*real(g_15_258*H(:,15)'*p_258) + 1;
E_15_268 = abs(g_15_268)^2*T_15_268 - 2*real(g_15_268*H(:,15)'*p_268) + 1;
E_15_278 = abs(g_15_278)^2*T_15_278 - 2*real(g_15_278*H(:,15)'*p_278) + 1;
E_15_289 = abs(g_15_289)^2*T_15_289 - 2*real(g_15_289*H(:,15)'*p_289) + 1;
E_15_348 = abs(g_15_348)^2*T_15_348 - 2*real(g_15_348*H(:,15)'*p_348) + 1;
E_15_358 = abs(g_15_358)^2*T_15_358 - 2*real(g_15_358*H(:,15)'*p_358) + 1;
E_15_368 = abs(g_15_368)^2*T_15_368 - 2*real(g_15_368*H(:,15)'*p_368) + 1;
E_15_378 = abs(g_15_378)^2*T_15_378 - 2*real(g_15_378*H(:,15)'*p_378) + 1;
E_15_389 = abs(g_15_389)^2*T_15_389 - 2*real(g_15_389*H(:,15)'*p_389) + 1;
E_15_458 = abs(g_15_458)^2*T_15_458 - 2*real(g_15_458*H(:,15)'*p_458) + 1;
E_15_468 = abs(g_15_468)^2*T_15_468 - 2*real(g_15_468*H(:,15)'*p_468) + 1;
E_15_478 = abs(g_15_478)^2*T_15_478 - 2*real(g_15_478*H(:,15)'*p_478) + 1;
E_15_489 = abs(g_15_489)^2*T_15_489 - 2*real(g_15_489*H(:,15)'*p_489) + 1;
E_15_568 = abs(g_15_568)^2*T_15_568 - 2*real(g_15_568*H(:,15)'*p_568) + 1;
E_15_578 = abs(g_15_578)^2*T_15_578 - 2*real(g_15_578*H(:,15)'*p_578) + 1;
E_15_589 = abs(g_15_589)^2*T_15_589 - 2*real(g_15_589*H(:,15)'*p_589) + 1;
E_15_678 = abs(g_15_678)^2*T_15_678 - 2*real(g_15_678*H(:,15)'*p_678) + 1;
E_15_689 = abs(g_15_689)^2*T_15_689 - 2*real(g_15_689*H(:,15)'*p_689) + 1;
E_15_789 = abs(g_15_789)^2*T_15_789 - 2*real(g_15_789*H(:,15)'*p_789) + 1;
E_15 = abs(g_15)^2*T_15 - 2*real(g_15*H(:,15)'*p_8) + 1;

E_16_123456789 = abs(g_16_123456789)^2*T_16_123456789 - 2*real(g_16_123456789*H(:,16)'*p_123456789) + 1;
E_16_128 = abs(g_16_128)^2*T_16_128 - 2*real(g_16_128*H(:,16)'*p_128) + 1;
E_16_138 = abs(g_16_138)^2*T_16_138 - 2*real(g_16_138*H(:,16)'*p_138) + 1;
E_16_148 = abs(g_16_148)^2*T_16_148 - 2*real(g_16_148*H(:,16)'*p_148) + 1;
E_16_158 = abs(g_16_158)^2*T_16_158 - 2*real(g_16_158*H(:,16)'*p_158) + 1;
E_16_168 = abs(g_16_168)^2*T_16_168 - 2*real(g_16_168*H(:,16)'*p_168) + 1;
E_16_178 = abs(g_16_178)^2*T_16_178 - 2*real(g_16_178*H(:,16)'*p_178) + 1;
E_16_189 = abs(g_16_189)^2*T_16_189 - 2*real(g_16_189*H(:,16)'*p_189) + 1;
E_16_238 = abs(g_16_238)^2*T_16_238 - 2*real(g_16_238*H(:,16)'*p_238) + 1;
E_16_248 = abs(g_16_248)^2*T_16_248 - 2*real(g_16_248*H(:,16)'*p_248) + 1;
E_16_258 = abs(g_16_258)^2*T_16_258 - 2*real(g_16_258*H(:,16)'*p_258) + 1;
E_16_268 = abs(g_16_268)^2*T_16_268 - 2*real(g_16_268*H(:,16)'*p_268) + 1;
E_16_278 = abs(g_16_278)^2*T_16_278 - 2*real(g_16_278*H(:,16)'*p_278) + 1;
E_16_289 = abs(g_16_289)^2*T_16_289 - 2*real(g_16_289*H(:,16)'*p_289) + 1;
E_16_348 = abs(g_16_348)^2*T_16_348 - 2*real(g_16_348*H(:,16)'*p_348) + 1;
E_16_358 = abs(g_16_358)^2*T_16_358 - 2*real(g_16_358*H(:,16)'*p_358) + 1;
E_16_368 = abs(g_16_368)^2*T_16_368 - 2*real(g_16_368*H(:,16)'*p_368) + 1;
E_16_378 = abs(g_16_378)^2*T_16_378 - 2*real(g_16_378*H(:,16)'*p_378) + 1;
E_16_389 = abs(g_16_389)^2*T_16_389 - 2*real(g_16_389*H(:,16)'*p_389) + 1;
E_16_458 = abs(g_16_458)^2*T_16_458 - 2*real(g_16_458*H(:,16)'*p_458) + 1;
E_16_468 = abs(g_16_468)^2*T_16_468 - 2*real(g_16_468*H(:,16)'*p_468) + 1;
E_16_478 = abs(g_16_478)^2*T_16_478 - 2*real(g_16_478*H(:,16)'*p_478) + 1;
E_16_489 = abs(g_16_489)^2*T_16_489 - 2*real(g_16_489*H(:,16)'*p_489) + 1;
E_16_568 = abs(g_16_568)^2*T_16_568 - 2*real(g_16_568*H(:,16)'*p_568) + 1;
E_16_578 = abs(g_16_578)^2*T_16_578 - 2*real(g_16_578*H(:,16)'*p_578) + 1;
E_16_589 = abs(g_16_589)^2*T_16_589 - 2*real(g_16_589*H(:,16)'*p_589) + 1;
E_16_678 = abs(g_16_678)^2*T_16_678 - 2*real(g_16_678*H(:,16)'*p_678) + 1;
E_16_689 = abs(g_16_689)^2*T_16_689 - 2*real(g_16_689*H(:,16)'*p_689) + 1;
E_16_789 = abs(g_16_789)^2*T_16_789 - 2*real(g_16_789*H(:,16)'*p_789) + 1;
E_16 = abs(g_16)^2*T_16 - 2*real(g_16*H(:,16)'*p_8) + 1;

E_17_123456789 = abs(g_17_123456789)^2*T_17_123456789 - 2*real(g_17_123456789*H(:,17)'*p_123456789) + 1;
E_17_129 = abs(g_17_129)^2*T_17_129 - 2*real(g_17_129*H(:,17)'*p_129) + 1;
E_17_139 = abs(g_17_139)^2*T_17_139 - 2*real(g_17_139*H(:,17)'*p_139) + 1;
E_17_149 = abs(g_17_149)^2*T_17_149 - 2*real(g_17_149*H(:,17)'*p_149) + 1;
E_17_159 = abs(g_17_159)^2*T_17_159 - 2*real(g_17_159*H(:,17)'*p_159) + 1;
E_17_169 = abs(g_17_169)^2*T_17_169 - 2*real(g_17_169*H(:,17)'*p_169) + 1;
E_17_179 = abs(g_17_179)^2*T_17_179 - 2*real(g_17_179*H(:,17)'*p_179) + 1;
E_17_189 = abs(g_17_189)^2*T_17_189 - 2*real(g_17_189*H(:,17)'*p_189) + 1;
E_17_239 = abs(g_17_239)^2*T_17_239 - 2*real(g_17_239*H(:,17)'*p_239) + 1;
E_17_249 = abs(g_17_249)^2*T_17_249 - 2*real(g_17_249*H(:,17)'*p_249) + 1;
E_17_259 = abs(g_17_259)^2*T_17_259 - 2*real(g_17_259*H(:,17)'*p_259) + 1;
E_17_269 = abs(g_17_269)^2*T_17_269 - 2*real(g_17_269*H(:,17)'*p_269) + 1;
E_17_279 = abs(g_17_279)^2*T_17_279 - 2*real(g_17_279*H(:,17)'*p_279) + 1;
E_17_289 = abs(g_17_289)^2*T_17_289 - 2*real(g_17_289*H(:,17)'*p_289) + 1;
E_17_349 = abs(g_17_349)^2*T_17_349 - 2*real(g_17_349*H(:,17)'*p_349) + 1;
E_17_359 = abs(g_17_359)^2*T_17_359 - 2*real(g_17_359*H(:,17)'*p_359) + 1;
E_17_369 = abs(g_17_369)^2*T_17_369 - 2*real(g_17_369*H(:,17)'*p_369) + 1;
E_17_379 = abs(g_17_379)^2*T_17_379 - 2*real(g_17_379*H(:,17)'*p_379) + 1;
E_17_389 = abs(g_17_389)^2*T_17_389 - 2*real(g_17_389*H(:,17)'*p_389) + 1;
E_17_459 = abs(g_17_459)^2*T_17_459 - 2*real(g_17_459*H(:,17)'*p_459) + 1;
E_17_469 = abs(g_17_469)^2*T_17_469 - 2*real(g_17_469*H(:,17)'*p_469) + 1;
E_17_479 = abs(g_17_479)^2*T_17_479 - 2*real(g_17_479*H(:,17)'*p_479) + 1;
E_17_489 = abs(g_17_489)^2*T_17_489 - 2*real(g_17_489*H(:,17)'*p_489) + 1;
E_17_569 = abs(g_17_569)^2*T_17_569 - 2*real(g_17_569*H(:,17)'*p_569) + 1;
E_17_579 = abs(g_17_579)^2*T_17_579 - 2*real(g_17_579*H(:,17)'*p_579) + 1;
E_17_589 = abs(g_17_589)^2*T_17_589 - 2*real(g_17_589*H(:,17)'*p_589) + 1;
E_17_679 = abs(g_17_679)^2*T_17_679 - 2*real(g_17_679*H(:,17)'*p_679) + 1;
E_17_689 = abs(g_17_689)^2*T_17_689 - 2*real(g_17_689*H(:,17)'*p_689) + 1;
E_17_789 = abs(g_17_789)^2*T_17_789 - 2*real(g_17_789*H(:,17)'*p_789) + 1;
E_17 = abs(g_17)^2*T_17 - 2*real(g_17*H(:,17)'*p_9) + 1;

E_18_123456789 = abs(g_18_123456789)^2*T_18_123456789 - 2*real(g_18_123456789*H(:,18)'*p_123456789) + 1;
E_18_129 = abs(g_18_129)^2*T_18_129 - 2*real(g_18_129*H(:,18)'*p_129) + 1;
E_18_139 = abs(g_18_139)^2*T_18_139 - 2*real(g_18_139*H(:,18)'*p_139) + 1;
E_18_149 = abs(g_18_149)^2*T_18_149 - 2*real(g_18_149*H(:,18)'*p_149) + 1;
E_18_159 = abs(g_18_159)^2*T_18_159 - 2*real(g_18_159*H(:,18)'*p_159) + 1;
E_18_169 = abs(g_18_169)^2*T_18_169 - 2*real(g_18_169*H(:,18)'*p_169) + 1;
E_18_179 = abs(g_18_179)^2*T_18_179 - 2*real(g_18_179*H(:,18)'*p_179) + 1;
E_18_189 = abs(g_18_189)^2*T_18_189 - 2*real(g_18_189*H(:,18)'*p_189) + 1;
E_18_239 = abs(g_18_239)^2*T_18_239 - 2*real(g_18_239*H(:,18)'*p_239) + 1;
E_18_249 = abs(g_18_249)^2*T_18_249 - 2*real(g_18_249*H(:,18)'*p_249) + 1;
E_18_259 = abs(g_18_259)^2*T_18_259 - 2*real(g_18_259*H(:,18)'*p_259) + 1;
E_18_269 = abs(g_18_269)^2*T_18_269 - 2*real(g_18_269*H(:,18)'*p_269) + 1;
E_18_279 = abs(g_18_279)^2*T_18_279 - 2*real(g_18_279*H(:,18)'*p_279) + 1;
E_18_289 = abs(g_18_289)^2*T_18_289 - 2*real(g_18_289*H(:,18)'*p_289) + 1;
E_18_349 = abs(g_18_349)^2*T_18_349 - 2*real(g_18_349*H(:,18)'*p_349) + 1;
E_18_359 = abs(g_18_359)^2*T_18_359 - 2*real(g_18_359*H(:,18)'*p_359) + 1;
E_18_369 = abs(g_18_369)^2*T_18_369 - 2*real(g_18_369*H(:,18)'*p_369) + 1;
E_18_379 = abs(g_18_379)^2*T_18_379 - 2*real(g_18_379*H(:,18)'*p_379) + 1;
E_18_389 = abs(g_18_389)^2*T_18_389 - 2*real(g_18_389*H(:,18)'*p_389) + 1;
E_18_459 = abs(g_18_459)^2*T_18_459 - 2*real(g_18_459*H(:,18)'*p_459) + 1;
E_18_469 = abs(g_18_469)^2*T_18_469 - 2*real(g_18_469*H(:,18)'*p_469) + 1;
E_18_479 = abs(g_18_479)^2*T_18_479 - 2*real(g_18_479*H(:,18)'*p_479) + 1;
E_18_489 = abs(g_18_489)^2*T_18_489 - 2*real(g_18_489*H(:,18)'*p_489) + 1;
E_18_569 = abs(g_18_569)^2*T_18_569 - 2*real(g_18_569*H(:,18)'*p_569) + 1;
E_18_579 = abs(g_18_579)^2*T_18_579 - 2*real(g_18_579*H(:,18)'*p_579) + 1;
E_18_589 = abs(g_18_589)^2*T_18_589 - 2*real(g_18_589*H(:,18)'*p_589) + 1;
E_18_679 = abs(g_18_679)^2*T_18_679 - 2*real(g_18_679*H(:,18)'*p_679) + 1;
E_18_689 = abs(g_18_689)^2*T_18_689 - 2*real(g_18_689*H(:,18)'*p_689) + 1;
E_18_789 = abs(g_18_789)^2*T_18_789 - 2*real(g_18_789*H(:,18)'*p_789) + 1;
E_18 = abs(g_18)^2*T_18 - 2*real(g_18*H(:,18)'*p_9) + 1;


%Rate-WMMSE Relationship
X_1_123456789 = W_1_123456789*E_1_123456789 - log2(W_1_123456789);
X_1_123 = W_1_123*E_1_123 - log2(W_1_123);
X_1_124 = W_1_124*E_1_124 - log2(W_1_124);
X_1_125 = W_1_125*E_1_125 - log2(W_1_125);
X_1_126 = W_1_126*E_1_126 - log2(W_1_126);
X_1_127 = W_1_127*E_1_127 - log2(W_1_127);
X_1_128 = W_1_128*E_1_128 - log2(W_1_128);
X_1_129 = W_1_129*E_1_129 - log2(W_1_129);
X_1_134 = W_1_134*E_1_134 - log2(W_1_134);
X_1_135 = W_1_135*E_1_135 - log2(W_1_135);
X_1_136 = W_1_136*E_1_136 - log2(W_1_136);
X_1_137 = W_1_137*E_1_137 - log2(W_1_137);
X_1_138 = W_1_138*E_1_138 - log2(W_1_138);
X_1_139 = W_1_139*E_1_139 - log2(W_1_139);
X_1_145 = W_1_145*E_1_145 - log2(W_1_145);
X_1_146 = W_1_146*E_1_146 - log2(W_1_146);
X_1_147 = W_1_147*E_1_147 - log2(W_1_147);
X_1_148 = W_1_148*E_1_148 - log2(W_1_148);
X_1_149 = W_1_149*E_1_149 - log2(W_1_149);
X_1_156 = W_1_156*E_1_156 - log2(W_1_156);
X_1_157 = W_1_157*E_1_157 - log2(W_1_157);
X_1_158 = W_1_158*E_1_158 - log2(W_1_158);
X_1_159 = W_1_159*E_1_159 - log2(W_1_159);
X_1_167 = W_1_167*E_1_167 - log2(W_1_167);
X_1_168 = W_1_168*E_1_168 - log2(W_1_168);
X_1_169 = W_1_169*E_1_169 - log2(W_1_169);
X_1_178 = W_1_178*E_1_178 - log2(W_1_178);
X_1_179 = W_1_179*E_1_179 - log2(W_1_179);
X_1_189 = W_1_189*E_1_189 - log2(W_1_189);
X_1 = W_1*E_1 - log2(W_1);

X_2_123456789 = W_2_123456789*E_2_123456789 - log2(W_2_123456789);
X_2_123 = W_2_123*E_2_123 - log2(W_2_123);
X_2_124 = W_2_124*E_2_124 - log2(W_2_124);
X_2_125 = W_2_125*E_2_125 - log2(W_2_125);
X_2_126 = W_2_126*E_2_126 - log2(W_2_126);
X_2_127 = W_2_127*E_2_127 - log2(W_2_127);
X_2_128 = W_2_128*E_2_128 - log2(W_2_128);
X_2_129 = W_2_129*E_2_129 - log2(W_2_129);
X_2_134 = W_2_134*E_2_134 - log2(W_2_134);
X_2_135 = W_2_135*E_2_135 - log2(W_2_135);
X_2_136 = W_2_136*E_2_136 - log2(W_2_136);
X_2_137 = W_2_137*E_2_137 - log2(W_2_137);
X_2_138 = W_2_138*E_2_138 - log2(W_2_138);
X_2_139 = W_2_139*E_2_139 - log2(W_2_139);
X_2_145 = W_2_145*E_2_145 - log2(W_2_145);
X_2_146 = W_2_146*E_2_146 - log2(W_2_146);
X_2_147 = W_2_147*E_2_147 - log2(W_2_147);
X_2_148 = W_2_148*E_2_148 - log2(W_2_148);
X_2_149 = W_2_149*E_2_149 - log2(W_2_149);
X_2_156 = W_2_156*E_2_156 - log2(W_2_156);
X_2_157 = W_2_157*E_2_157 - log2(W_2_157);
X_2_158 = W_2_158*E_2_158 - log2(W_2_158);
X_2_159 = W_2_159*E_2_159 - log2(W_2_159);
X_2_167 = W_2_167*E_2_167 - log2(W_2_167);
X_2_168 = W_2_168*E_2_168 - log2(W_2_168);
X_2_169 = W_2_169*E_2_169 - log2(W_2_169);
X_2_178 = W_2_178*E_2_178 - log2(W_2_178);
X_2_179 = W_2_179*E_2_179 - log2(W_2_179);
X_2_189 = W_2_189*E_2_189 - log2(W_2_189);
X_2 = W_2*E_2 - log2(W_2);

X_3_123456789 = W_3_123456789*E_3_123456789 - log2(W_3_123456789);
X_3_123 = W_3_123*E_3_123 - log2(W_3_123);
X_3_124 = W_3_124*E_3_124 - log2(W_3_124);
X_3_125 = W_3_125*E_3_125 - log2(W_3_125);
X_3_126 = W_3_126*E_3_126 - log2(W_3_126);
X_3_127 = W_3_127*E_3_127 - log2(W_3_127);
X_3_128 = W_3_128*E_3_128 - log2(W_3_128);
X_3_129 = W_3_129*E_3_129 - log2(W_3_129);
X_3_234 = W_3_234*E_3_234 - log2(W_3_234);
X_3_235 = W_3_235*E_3_235 - log2(W_3_235);
X_3_236 = W_3_236*E_3_236 - log2(W_3_236);
X_3_237 = W_3_237*E_3_237 - log2(W_3_237);
X_3_238 = W_3_238*E_3_238 - log2(W_3_238);
X_3_239 = W_3_239*E_3_239 - log2(W_3_239);
X_3_245 = W_3_245*E_3_245 - log2(W_3_245);
X_3_246 = W_3_246*E_3_246 - log2(W_3_246);
X_3_247 = W_3_247*E_3_247 - log2(W_3_247);
X_3_248 = W_3_248*E_3_248 - log2(W_3_248);
X_3_249 = W_3_249*E_3_249 - log2(W_3_249);
X_3_256 = W_3_256*E_3_256 - log2(W_3_256);
X_3_257 = W_3_257*E_3_257 - log2(W_3_257);
X_3_258 = W_3_258*E_3_258 - log2(W_3_258);
X_3_259 = W_3_259*E_3_259 - log2(W_3_259);
X_3_267 = W_3_267*E_3_267 - log2(W_3_267);
X_3_268 = W_3_268*E_3_268 - log2(W_3_268);
X_3_269 = W_3_269*E_3_269 - log2(W_3_269);
X_3_278 = W_3_278*E_3_278 - log2(W_3_278);
X_3_279 = W_3_279*E_3_279 - log2(W_3_279);
X_3_289 = W_3_289*E_3_289 - log2(W_3_289);
X_3 = W_3*E_3 - log2(W_3);

X_4_123456789 = W_4_123456789*E_4_123456789 - log2(W_4_123456789);
X_4_123 = W_4_123*E_4_123 - log2(W_4_123);
X_4_124 = W_4_124*E_4_124 - log2(W_4_124);
X_4_125 = W_4_125*E_4_125 - log2(W_4_125);
X_4_126 = W_4_126*E_4_126 - log2(W_4_126);
X_4_127 = W_4_127*E_4_127 - log2(W_4_127);
X_4_128 = W_4_128*E_4_128 - log2(W_4_128);
X_4_129 = W_4_129*E_4_129 - log2(W_4_129);
X_4_234 = W_4_234*E_4_234 - log2(W_4_234);
X_4_235 = W_4_235*E_4_235 - log2(W_4_235);
X_4_236 = W_4_236*E_4_236 - log2(W_4_236);
X_4_237 = W_4_237*E_4_237 - log2(W_4_237);
X_4_238 = W_4_238*E_4_238 - log2(W_4_238);
X_4_239 = W_4_239*E_4_239 - log2(W_4_239);
X_4_245 = W_4_245*E_4_245 - log2(W_4_245);
X_4_246 = W_4_246*E_4_246 - log2(W_4_246);
X_4_247 = W_4_247*E_4_247 - log2(W_4_247);
X_4_248 = W_4_248*E_4_248 - log2(W_4_248);
X_4_249 = W_4_249*E_4_249 - log2(W_4_249);
X_4_256 = W_4_256*E_4_256 - log2(W_4_256);
X_4_257 = W_4_257*E_4_257 - log2(W_4_257);
X_4_258 = W_4_258*E_4_258 - log2(W_4_258);
X_4_259 = W_4_259*E_4_259 - log2(W_4_259);
X_4_267 = W_4_267*E_4_267 - log2(W_4_267);
X_4_268 = W_4_268*E_4_268 - log2(W_4_268);
X_4_269 = W_4_269*E_4_269 - log2(W_4_269);
X_4_278 = W_4_278*E_4_278 - log2(W_4_278);
X_4_279 = W_4_279*E_4_279 - log2(W_4_279);
X_4_289 = W_4_289*E_4_289 - log2(W_4_289);
X_4 = W_4*E_4 - log2(W_4);

X_5_123456789 = W_5_123456789*E_5_123456789 - log2(W_5_123456789);
X_5_123 = W_5_123*E_5_123 - log2(W_5_123);
X_5_134 = W_5_134*E_5_134 - log2(W_5_134);
X_5_135 = W_5_135*E_5_135 - log2(W_5_135);
X_5_136 = W_5_136*E_5_136 - log2(W_5_136);
X_5_137 = W_5_137*E_5_137 - log2(W_5_137);
X_5_138 = W_5_138*E_5_138 - log2(W_5_138);
X_5_139 = W_5_139*E_5_139 - log2(W_5_139);
X_5_234 = W_5_234*E_5_234 - log2(W_5_234);
X_5_235 = W_5_235*E_5_235 - log2(W_5_235);
X_5_236 = W_5_236*E_5_236 - log2(W_5_236);
X_5_237 = W_5_237*E_5_237 - log2(W_5_237);
X_5_238 = W_5_238*E_5_238 - log2(W_5_238);
X_5_239 = W_5_239*E_5_239 - log2(W_5_239);
X_5_345 = W_5_345*E_5_345 - log2(W_5_345);
X_5_346 = W_5_346*E_5_346 - log2(W_5_346);
X_5_347 = W_5_347*E_5_347 - log2(W_5_347);
X_5_348 = W_5_348*E_5_348 - log2(W_5_348);
X_5_349 = W_5_349*E_5_349 - log2(W_5_349);
X_5_356 = W_5_356*E_5_356 - log2(W_5_356);
X_5_357 = W_5_357*E_5_357 - log2(W_5_357);
X_5_358 = W_5_358*E_5_358 - log2(W_5_358);
X_5_359 = W_5_359*E_5_359 - log2(W_5_359);
X_5_367 = W_5_367*E_5_367 - log2(W_5_367);
X_5_368 = W_5_368*E_5_368 - log2(W_5_368);
X_5_369 = W_5_369*E_5_369 - log2(W_5_369);
X_5_378 = W_5_378*E_5_378 - log2(W_5_378);
X_5_379 = W_5_379*E_5_379 - log2(W_5_379);
X_5_389 = W_5_389*E_5_389 - log2(W_5_389);
X_5 = W_5*E_5 - log2(W_5);

X_6_123456789 = W_6_123456789*E_6_123456789 - log2(W_6_123456789);
X_6_123 = W_6_123*E_6_123 - log2(W_6_123);
X_6_134 = W_6_134*E_6_134 - log2(W_6_134);
X_6_135 = W_6_135*E_6_135 - log2(W_6_135);
X_6_136 = W_6_136*E_6_136 - log2(W_6_136);
X_6_137 = W_6_137*E_6_137 - log2(W_6_137);
X_6_138 = W_6_138*E_6_138 - log2(W_6_138);
X_6_139 = W_6_139*E_6_139 - log2(W_6_139);
X_6_234 = W_6_234*E_6_234 - log2(W_6_234);
X_6_235 = W_6_235*E_6_235 - log2(W_6_235);
X_6_236 = W_6_236*E_6_236 - log2(W_6_236);
X_6_237 = W_6_237*E_6_237 - log2(W_6_237);
X_6_238 = W_6_238*E_6_238 - log2(W_6_238);
X_6_239 = W_6_239*E_6_239 - log2(W_6_239);
X_6_345 = W_6_345*E_6_345 - log2(W_6_345);
X_6_346 = W_6_346*E_6_346 - log2(W_6_346);
X_6_347 = W_6_347*E_6_347 - log2(W_6_347);
X_6_348 = W_6_348*E_6_348 - log2(W_6_348);
X_6_349 = W_6_349*E_6_349 - log2(W_6_349);
X_6_356 = W_6_356*E_6_356 - log2(W_6_356);
X_6_357 = W_6_357*E_6_357 - log2(W_6_357);
X_6_358 = W_6_358*E_6_358 - log2(W_6_358);
X_6_359 = W_6_359*E_6_359 - log2(W_6_359);
X_6_367 = W_6_367*E_6_367 - log2(W_6_367);
X_6_368 = W_6_368*E_6_368 - log2(W_6_368);
X_6_369 = W_6_369*E_6_369 - log2(W_6_369);
X_6_378 = W_6_378*E_6_378 - log2(W_6_378);
X_6_379 = W_6_379*E_6_379 - log2(W_6_379);
X_6_389 = W_6_389*E_6_389 - log2(W_6_389);
X_6 = W_6*E_6 - log2(W_6);

X_7_123456789 = W_7_123456789*E_7_123456789 - log2(W_7_123456789);
X_7_124 = W_7_124*E_7_124 - log2(W_7_124);
X_7_134 = W_7_134*E_7_134 - log2(W_7_134);
X_7_145 = W_7_145*E_7_145 - log2(W_7_145);
X_7_146 = W_7_146*E_7_146 - log2(W_7_146);
X_7_147 = W_7_147*E_7_147 - log2(W_7_147);
X_7_148 = W_7_148*E_7_148 - log2(W_7_148);
X_7_149 = W_7_149*E_7_149 - log2(W_7_149);
X_7_234 = W_7_234*E_7_234 - log2(W_7_234);
X_7_245 = W_7_245*E_7_245 - log2(W_7_245);
X_7_246 = W_7_246*E_7_246 - log2(W_7_246);
X_7_247 = W_7_247*E_7_247 - log2(W_7_247);
X_7_248 = W_7_248*E_7_248 - log2(W_7_248);
X_7_249 = W_7_249*E_7_249 - log2(W_7_249);
X_7_345 = W_7_345*E_7_345 - log2(W_7_345);
X_7_346 = W_7_346*E_7_346 - log2(W_7_346);
X_7_347 = W_7_347*E_7_347 - log2(W_7_347);
X_7_348 = W_7_348*E_7_348 - log2(W_7_348);
X_7_349 = W_7_349*E_7_349 - log2(W_7_349);
X_7_456 = W_7_456*E_7_456 - log2(W_7_456);
X_7_457 = W_7_457*E_7_457 - log2(W_7_457);
X_7_458 = W_7_458*E_7_458 - log2(W_7_458);
X_7_459 = W_7_459*E_7_459 - log2(W_7_459);
X_7_467 = W_7_467*E_7_467 - log2(W_7_467);
X_7_468 = W_7_468*E_7_468 - log2(W_7_468);
X_7_469 = W_7_469*E_7_469 - log2(W_7_469);
X_7_478 = W_7_478*E_7_478 - log2(W_7_478);
X_7_479 = W_7_479*E_7_479 - log2(W_7_479);
X_7_489 = W_7_489*E_7_489 - log2(W_7_489);
X_7 = W_7*E_7 - log2(W_7);

X_8_123456789 = W_8_123456789*E_8_123456789 - log2(W_8_123456789);
X_8_124 = W_8_124*E_8_124 - log2(W_8_124);
X_8_134 = W_8_134*E_8_134 - log2(W_8_134);
X_8_145 = W_8_145*E_8_145 - log2(W_8_145);
X_8_146 = W_8_146*E_8_146 - log2(W_8_146);
X_8_147 = W_8_147*E_8_147 - log2(W_8_147);
X_8_148 = W_8_148*E_8_148 - log2(W_8_148);
X_8_149 = W_8_149*E_8_149 - log2(W_8_149);
X_8_234 = W_8_234*E_8_234 - log2(W_8_234);
X_8_245 = W_8_245*E_8_245 - log2(W_8_245);
X_8_246 = W_8_246*E_8_246 - log2(W_8_246);
X_8_247 = W_8_247*E_8_247 - log2(W_8_247);
X_8_248 = W_8_248*E_8_248 - log2(W_8_248);
X_8_249 = W_8_249*E_8_249 - log2(W_8_249);
X_8_345 = W_8_345*E_8_345 - log2(W_8_345);
X_8_346 = W_8_346*E_8_346 - log2(W_8_346);
X_8_347 = W_8_347*E_8_347 - log2(W_8_347);
X_8_348 = W_8_348*E_8_348 - log2(W_8_348);
X_8_349 = W_8_349*E_8_349 - log2(W_8_349);
X_8_456 = W_8_456*E_8_456 - log2(W_8_456);
X_8_457 = W_8_457*E_8_457 - log2(W_8_457);
X_8_458 = W_8_458*E_8_458 - log2(W_8_458);
X_8_459 = W_8_459*E_8_459 - log2(W_8_459);
X_8_467 = W_8_467*E_8_467 - log2(W_8_467);
X_8_468 = W_8_468*E_8_468 - log2(W_8_468);
X_8_469 = W_8_469*E_8_469 - log2(W_8_469);
X_8_478 = W_8_478*E_8_478 - log2(W_8_478);
X_8_479 = W_8_479*E_8_479 - log2(W_8_479);
X_8_489 = W_8_489*E_8_489 - log2(W_8_489);
X_8 = W_8*E_8 - log2(W_8);

X_9_123456789 = W_9_123456789*E_9_123456789 - log2(W_9_123456789);
X_9_125 = W_9_125*E_9_125 - log2(W_9_125);
X_9_135 = W_9_135*E_9_135 - log2(W_9_135);
X_9_145 = W_9_145*E_9_145 - log2(W_9_145);
X_9_156 = W_9_156*E_9_156 - log2(W_9_156);
X_9_157 = W_9_157*E_9_157 - log2(W_9_157);
X_9_158 = W_9_158*E_9_158 - log2(W_9_158);
X_9_159 = W_9_159*E_9_159 - log2(W_9_159);
X_9_235 = W_9_235*E_9_235 - log2(W_9_235);
X_9_245 = W_9_245*E_9_245 - log2(W_9_245);
X_9_256 = W_9_256*E_9_256 - log2(W_9_256);
X_9_257 = W_9_257*E_9_257 - log2(W_9_257);
X_9_258 = W_9_258*E_9_258 - log2(W_9_258);
X_9_259 = W_9_259*E_9_259 - log2(W_9_259);
X_9_345 = W_9_345*E_9_345 - log2(W_9_345);
X_9_356 = W_9_356*E_9_356 - log2(W_9_356);
X_9_357 = W_9_357*E_9_357 - log2(W_9_357);
X_9_358 = W_9_358*E_9_358 - log2(W_9_358);
X_9_359 = W_9_359*E_9_359 - log2(W_9_359);
X_9_456 = W_9_456*E_9_456 - log2(W_9_456);
X_9_457 = W_9_457*E_9_457 - log2(W_9_457);
X_9_458 = W_9_458*E_9_458 - log2(W_9_458);
X_9_459 = W_9_459*E_9_459 - log2(W_9_459);
X_9_567 = W_9_567*E_9_567 - log2(W_9_567);
X_9_568 = W_9_568*E_9_568 - log2(W_9_568);
X_9_569 = W_9_569*E_9_569 - log2(W_9_569);
X_9_578 = W_9_578*E_9_578 - log2(W_9_578);
X_9_579 = W_9_579*E_9_579 - log2(W_9_579);
X_9_589 = W_9_589*E_9_589 - log2(W_9_589);
X_9 = W_9*E_9 - log2(W_9);

X_10_123456789 = W_10_123456789*E_10_123456789 - log2(W_10_123456789);
X_10_125 = W_10_125*E_10_125 - log2(W_10_125);
X_10_135 = W_10_135*E_10_135 - log2(W_10_135);
X_10_145 = W_10_145*E_10_145 - log2(W_10_145);
X_10_156 = W_10_156*E_10_156 - log2(W_10_156);
X_10_157 = W_10_157*E_10_157 - log2(W_10_157);
X_10_158 = W_10_158*E_10_158 - log2(W_10_158);
X_10_159 = W_10_159*E_10_159 - log2(W_10_159);
X_10_235 = W_10_235*E_10_235 - log2(W_10_235);
X_10_245 = W_10_245*E_10_245 - log2(W_10_245);
X_10_256 = W_10_256*E_10_256 - log2(W_10_256);
X_10_257 = W_10_257*E_10_257 - log2(W_10_257);
X_10_258 = W_10_258*E_10_258 - log2(W_10_258);
X_10_259 = W_10_259*E_10_259 - log2(W_10_259);
X_10_345 = W_10_345*E_10_345 - log2(W_10_345);
X_10_356 = W_10_356*E_10_356 - log2(W_10_356);
X_10_357 = W_10_357*E_10_357 - log2(W_10_357);
X_10_358 = W_10_358*E_10_358 - log2(W_10_358);
X_10_359 = W_10_359*E_10_359 - log2(W_10_359);
X_10_456 = W_10_456*E_10_456 - log2(W_10_456);
X_10_457 = W_10_457*E_10_457 - log2(W_10_457);
X_10_458 = W_10_458*E_10_458 - log2(W_10_458);
X_10_459 = W_10_459*E_10_459 - log2(W_10_459);
X_10_567 = W_10_567*E_10_567 - log2(W_10_567);
X_10_568 = W_10_568*E_10_568 - log2(W_10_568);
X_10_569 = W_10_569*E_10_569 - log2(W_10_569);
X_10_578 = W_10_578*E_10_578 - log2(W_10_578);
X_10_579 = W_10_579*E_10_579 - log2(W_10_579);
X_10_589 = W_10_589*E_10_589 - log2(W_10_589);
X_10 = W_10*E_10 - log2(W_10);

X_11_123456789 = W_11_123456789*E_11_123456789 - log2(W_11_123456789);
X_11_126 = W_11_126*E_11_126 - log2(W_11_126);
X_11_136 = W_11_136*E_11_136 - log2(W_11_136);
X_11_146 = W_11_146*E_11_146 - log2(W_11_146);
X_11_156 = W_11_156*E_11_156 - log2(W_11_156);
X_11_167 = W_11_167*E_11_167 - log2(W_11_167);
X_11_168 = W_11_168*E_11_168 - log2(W_11_168);
X_11_169 = W_11_169*E_11_169 - log2(W_11_169);
X_11_236 = W_11_236*E_11_236 - log2(W_11_236);
X_11_246 = W_11_246*E_11_246 - log2(W_11_246);
X_11_256 = W_11_256*E_11_256 - log2(W_11_256);
X_11_267 = W_11_267*E_11_267 - log2(W_11_267);
X_11_268 = W_11_268*E_11_268 - log2(W_11_268);
X_11_269 = W_11_269*E_11_269 - log2(W_11_269);
X_11_346 = W_11_346*E_11_346 - log2(W_11_346);
X_11_356 = W_11_356*E_11_356 - log2(W_11_356);
X_11_367 = W_11_367*E_11_367 - log2(W_11_367);
X_11_368 = W_11_368*E_11_368 - log2(W_11_368);
X_11_369 = W_11_369*E_11_369 - log2(W_11_369);
X_11_456 = W_11_456*E_11_456 - log2(W_11_456);
X_11_467 = W_11_467*E_11_467 - log2(W_11_467);
X_11_468 = W_11_468*E_11_468 - log2(W_11_468);
X_11_469 = W_11_469*E_11_469 - log2(W_11_469);
X_11_567 = W_11_567*E_11_567 - log2(W_11_567);
X_11_568 = W_11_568*E_11_568 - log2(W_11_568);
X_11_569 = W_11_569*E_11_569 - log2(W_11_569);
X_11_678 = W_11_678*E_11_678 - log2(W_11_678);
X_11_679 = W_11_679*E_11_679 - log2(W_11_679);
X_11_689 = W_11_689*E_11_689 - log2(W_11_689);
X_11 = W_11*E_11 - log2(W_11);

X_12_123456789 = W_12_123456789*E_12_123456789 - log2(W_12_123456789);
X_12_126 = W_12_126*E_12_126 - log2(W_12_126);
X_12_136 = W_12_136*E_12_136 - log2(W_12_136);
X_12_146 = W_12_146*E_12_146 - log2(W_12_146);
X_12_156 = W_12_156*E_12_156 - log2(W_12_156);
X_12_167 = W_12_167*E_12_167 - log2(W_12_167);
X_12_168 = W_12_168*E_12_168 - log2(W_12_168);
X_12_169 = W_12_169*E_12_169 - log2(W_12_169);
X_12_236 = W_12_236*E_12_236 - log2(W_12_236);
X_12_246 = W_12_246*E_12_246 - log2(W_12_246);
X_12_256 = W_12_256*E_12_256 - log2(W_12_256);
X_12_267 = W_12_267*E_12_267 - log2(W_12_267);
X_12_268 = W_12_268*E_12_268 - log2(W_12_268);
X_12_269 = W_12_269*E_12_269 - log2(W_12_269);
X_12_346 = W_12_346*E_12_346 - log2(W_12_346);
X_12_356 = W_12_356*E_12_356 - log2(W_12_356);
X_12_367 = W_12_367*E_12_367 - log2(W_12_367);
X_12_368 = W_12_368*E_12_368 - log2(W_12_368);
X_12_369 = W_12_369*E_12_369 - log2(W_12_369);
X_12_456 = W_12_456*E_12_456 - log2(W_12_456);
X_12_467 = W_12_467*E_12_467 - log2(W_12_467);
X_12_468 = W_12_468*E_12_468 - log2(W_12_468);
X_12_469 = W_12_469*E_12_469 - log2(W_12_469);
X_12_567 = W_12_567*E_12_567 - log2(W_12_567);
X_12_568 = W_12_568*E_12_568 - log2(W_12_568);
X_12_569 = W_12_569*E_12_569 - log2(W_12_569);
X_12_678 = W_12_678*E_12_678 - log2(W_12_678);
X_12_679 = W_12_679*E_12_679 - log2(W_12_679);
X_12_689 = W_12_689*E_12_689 - log2(W_12_689);
X_12 = W_12*E_12 - log2(W_12);

X_13_123456789 = W_13_123456789*E_13_123456789 - log2(W_13_123456789);
X_13_127 = W_13_127*E_13_127 - log2(W_13_127);
X_13_137 = W_13_137*E_13_137 - log2(W_13_137);
X_13_147 = W_13_147*E_13_147 - log2(W_13_147);
X_13_157 = W_13_157*E_13_157 - log2(W_13_157);
X_13_167 = W_13_167*E_13_167 - log2(W_13_167);
X_13_178 = W_13_178*E_13_178 - log2(W_13_178);
X_13_179 = W_13_179*E_13_179 - log2(W_13_179);
X_13_237 = W_13_237*E_13_237 - log2(W_13_237);
X_13_247 = W_13_247*E_13_247 - log2(W_13_247);
X_13_257 = W_13_257*E_13_257 - log2(W_13_257);
X_13_267 = W_13_267*E_13_267 - log2(W_13_267);
X_13_278 = W_13_278*E_13_278 - log2(W_13_278);
X_13_279 = W_13_279*E_13_279 - log2(W_13_279);
X_13_347 = W_13_347*E_13_347 - log2(W_13_347);
X_13_357 = W_13_357*E_13_357 - log2(W_13_357);
X_13_367 = W_13_367*E_13_367 - log2(W_13_367);
X_13_378 = W_13_378*E_13_378 - log2(W_13_378);
X_13_379 = W_13_379*E_13_379 - log2(W_13_379);
X_13_457 = W_13_457*E_13_457 - log2(W_13_457);
X_13_467 = W_13_467*E_13_467 - log2(W_13_467);
X_13_478 = W_13_478*E_13_478 - log2(W_13_478);
X_13_479 = W_13_479*E_13_479 - log2(W_13_479);
X_13_567 = W_13_567*E_13_567 - log2(W_13_567);
X_13_578 = W_13_578*E_13_578 - log2(W_13_578);
X_13_579 = W_13_579*E_13_579 - log2(W_13_579);
X_13_678 = W_13_678*E_13_678 - log2(W_13_678);
X_13_679 = W_13_679*E_13_679 - log2(W_13_679);
X_13_789 = W_13_789*E_13_789 - log2(W_13_789);
X_13 = W_13*E_13 - log2(W_13);

X_14_123456789 = W_14_123456789*E_14_123456789 - log2(W_14_123456789);
X_14_127 = W_14_127*E_14_127 - log2(W_14_127);
X_14_137 = W_14_137*E_14_137 - log2(W_14_137);
X_14_147 = W_14_147*E_14_147 - log2(W_14_147);
X_14_157 = W_14_157*E_14_157 - log2(W_14_157);
X_14_167 = W_14_167*E_14_167 - log2(W_14_167);
X_14_178 = W_14_178*E_14_178 - log2(W_14_178);
X_14_179 = W_14_179*E_14_179 - log2(W_14_179);
X_14_237 = W_14_237*E_14_237 - log2(W_14_237);
X_14_247 = W_14_247*E_14_247 - log2(W_14_247);
X_14_257 = W_14_257*E_14_257 - log2(W_14_257);
X_14_267 = W_14_267*E_14_267 - log2(W_14_267);
X_14_278 = W_14_278*E_14_278 - log2(W_14_278);
X_14_279 = W_14_279*E_14_279 - log2(W_14_279);
X_14_347 = W_14_347*E_14_347 - log2(W_14_347);
X_14_357 = W_14_357*E_14_357 - log2(W_14_357);
X_14_367 = W_14_367*E_14_367 - log2(W_14_367);
X_14_378 = W_14_378*E_14_378 - log2(W_14_378);
X_14_379 = W_14_379*E_14_379 - log2(W_14_379);
X_14_457 = W_14_457*E_14_457 - log2(W_14_457);
X_14_467 = W_14_467*E_14_467 - log2(W_14_467);
X_14_478 = W_14_478*E_14_478 - log2(W_14_478);
X_14_479 = W_14_479*E_14_479 - log2(W_14_479);
X_14_567 = W_14_567*E_14_567 - log2(W_14_567);
X_14_578 = W_14_578*E_14_578 - log2(W_14_578);
X_14_579 = W_14_579*E_14_579 - log2(W_14_579);
X_14_678 = W_14_678*E_14_678 - log2(W_14_678);
X_14_679 = W_14_679*E_14_679 - log2(W_14_679);
X_14_789 = W_14_789*E_14_789 - log2(W_14_789);
X_14 = W_14*E_14 - log2(W_14);

X_15_123456789 = W_15_123456789*E_15_123456789 - log2(W_15_123456789);
X_15_128 = W_15_128*E_15_128 - log2(W_15_128);
X_15_138 = W_15_138*E_15_138 - log2(W_15_138);
X_15_148 = W_15_148*E_15_148 - log2(W_15_148);
X_15_158 = W_15_158*E_15_158 - log2(W_15_158);
X_15_168 = W_15_168*E_15_168 - log2(W_15_168);
X_15_178 = W_15_178*E_15_178 - log2(W_15_178);
X_15_189 = W_15_189*E_15_189 - log2(W_15_189);
X_15_238 = W_15_238*E_15_238 - log2(W_15_238);
X_15_248 = W_15_248*E_15_248 - log2(W_15_248);
X_15_258 = W_15_258*E_15_258 - log2(W_15_258);
X_15_268 = W_15_268*E_15_268 - log2(W_15_268);
X_15_278 = W_15_278*E_15_278 - log2(W_15_278);
X_15_289 = W_15_289*E_15_289 - log2(W_15_289);
X_15_348 = W_15_348*E_15_348 - log2(W_15_348);
X_15_358 = W_15_358*E_15_358 - log2(W_15_358);
X_15_368 = W_15_368*E_15_368 - log2(W_15_368);
X_15_378 = W_15_378*E_15_378 - log2(W_15_378);
X_15_389 = W_15_389*E_15_389 - log2(W_15_389);
X_15_458 = W_15_458*E_15_458 - log2(W_15_458);
X_15_468 = W_15_468*E_15_468 - log2(W_15_468);
X_15_478 = W_15_478*E_15_478 - log2(W_15_478);
X_15_489 = W_15_489*E_15_489 - log2(W_15_489);
X_15_568 = W_15_568*E_15_568 - log2(W_15_568);
X_15_578 = W_15_578*E_15_578 - log2(W_15_578);
X_15_589 = W_15_589*E_15_589 - log2(W_15_589);
X_15_678 = W_15_678*E_15_678 - log2(W_15_678);
X_15_689 = W_15_689*E_15_689 - log2(W_15_689);
X_15_789 = W_15_789*E_15_789 - log2(W_15_789);
X_15 = W_15*E_15 - log2(W_15);

X_16_123456789 = W_16_123456789*E_16_123456789 - log2(W_16_123456789);
X_16_128 = W_16_128*E_16_128 - log2(W_16_128);
X_16_138 = W_16_138*E_16_138 - log2(W_16_138);
X_16_148 = W_16_148*E_16_148 - log2(W_16_148);
X_16_158 = W_16_158*E_16_158 - log2(W_16_158);
X_16_168 = W_16_168*E_16_168 - log2(W_16_168);
X_16_178 = W_16_178*E_16_178 - log2(W_16_178);
X_16_189 = W_16_189*E_16_189 - log2(W_16_189);
X_16_238 = W_16_238*E_16_238 - log2(W_16_238);
X_16_248 = W_16_248*E_16_248 - log2(W_16_248);
X_16_258 = W_16_258*E_16_258 - log2(W_16_258);
X_16_268 = W_16_268*E_16_268 - log2(W_16_268);
X_16_278 = W_16_278*E_16_278 - log2(W_16_278);
X_16_289 = W_16_289*E_16_289 - log2(W_16_289);
X_16_348 = W_16_348*E_16_348 - log2(W_16_348);
X_16_358 = W_16_358*E_16_358 - log2(W_16_358);
X_16_368 = W_16_368*E_16_368 - log2(W_16_368);
X_16_378 = W_16_378*E_16_378 - log2(W_16_378);
X_16_389 = W_16_389*E_16_389 - log2(W_16_389);
X_16_458 = W_16_458*E_16_458 - log2(W_16_458);
X_16_468 = W_16_468*E_16_468 - log2(W_16_468);
X_16_478 = W_16_478*E_16_478 - log2(W_16_478);
X_16_489 = W_16_489*E_16_489 - log2(W_16_489);
X_16_568 = W_16_568*E_16_568 - log2(W_16_568);
X_16_578 = W_16_578*E_16_578 - log2(W_16_578);
X_16_589 = W_16_589*E_16_589 - log2(W_16_589);
X_16_678 = W_16_678*E_16_678 - log2(W_16_678);
X_16_689 = W_16_689*E_16_689 - log2(W_16_689);
X_16_789 = W_16_789*E_16_789 - log2(W_16_789);
X_16 = W_16*E_16 - log2(W_16);

X_17_123456789 = W_17_123456789*E_17_123456789 - log2(W_17_123456789);
X_17_129 = W_17_129*E_17_129 - log2(W_17_129);
X_17_139 = W_17_139*E_17_139 - log2(W_17_139);
X_17_149 = W_17_149*E_17_149 - log2(W_17_149);
X_17_159 = W_17_159*E_17_159 - log2(W_17_159);
X_17_169 = W_17_169*E_17_169 - log2(W_17_169);
X_17_179 = W_17_179*E_17_179 - log2(W_17_179);
X_17_189 = W_17_189*E_17_189 - log2(W_17_189);
X_17_239 = W_17_239*E_17_239 - log2(W_17_239);
X_17_249 = W_17_249*E_17_249 - log2(W_17_249);
X_17_259 = W_17_259*E_17_259 - log2(W_17_259);
X_17_269 = W_17_269*E_17_269 - log2(W_17_269);
X_17_279 = W_17_279*E_17_279 - log2(W_17_279);
X_17_289 = W_17_289*E_17_289 - log2(W_17_289);
X_17_349 = W_17_349*E_17_349 - log2(W_17_349);
X_17_359 = W_17_359*E_17_359 - log2(W_17_359);
X_17_369 = W_17_369*E_17_369 - log2(W_17_369);
X_17_379 = W_17_379*E_17_379 - log2(W_17_379);
X_17_389 = W_17_389*E_17_389 - log2(W_17_389);
X_17_459 = W_17_459*E_17_459 - log2(W_17_459);
X_17_469 = W_17_469*E_17_469 - log2(W_17_469);
X_17_479 = W_17_479*E_17_479 - log2(W_17_479);
X_17_489 = W_17_489*E_17_489 - log2(W_17_489);
X_17_569 = W_17_569*E_17_569 - log2(W_17_569);
X_17_579 = W_17_579*E_17_579 - log2(W_17_579);
X_17_589 = W_17_589*E_17_589 - log2(W_17_589);
X_17_679 = W_17_679*E_17_679 - log2(W_17_679);
X_17_689 = W_17_689*E_17_689 - log2(W_17_689);
X_17_789 = W_17_789*E_17_789 - log2(W_17_789);
X_17 = W_17*E_17 - log2(W_17);

X_18_123456789 = W_18_123456789*E_18_123456789 - log2(W_18_123456789);
X_18_129 = W_18_129*E_18_129 - log2(W_18_129);
X_18_139 = W_18_139*E_18_139 - log2(W_18_139);
X_18_149 = W_18_149*E_18_149 - log2(W_18_149);
X_18_159 = W_18_159*E_18_159 - log2(W_18_159);
X_18_169 = W_18_169*E_18_169 - log2(W_18_169);
X_18_179 = W_18_179*E_18_179 - log2(W_18_179);
X_18_189 = W_18_189*E_18_189 - log2(W_18_189);
X_18_239 = W_18_239*E_18_239 - log2(W_18_239);
X_18_249 = W_18_249*E_18_249 - log2(W_18_249);
X_18_259 = W_18_259*E_18_259 - log2(W_18_259);
X_18_269 = W_18_269*E_18_269 - log2(W_18_269);
X_18_279 = W_18_279*E_18_279 - log2(W_18_279);
X_18_289 = W_18_289*E_18_289 - log2(W_18_289);
X_18_349 = W_18_349*E_18_349 - log2(W_18_349);
X_18_359 = W_18_359*E_18_359 - log2(W_18_359);
X_18_369 = W_18_369*E_18_369 - log2(W_18_369);
X_18_379 = W_18_379*E_18_379 - log2(W_18_379);
X_18_389 = W_18_389*E_18_389 - log2(W_18_389);
X_18_459 = W_18_459*E_18_459 - log2(W_18_459);
X_18_469 = W_18_469*E_18_469 - log2(W_18_469);
X_18_479 = W_18_479*E_18_479 - log2(W_18_479);
X_18_489 = W_18_489*E_18_489 - log2(W_18_489);
X_18_569 = W_18_569*E_18_569 - log2(W_18_569);
X_18_579 = W_18_579*E_18_579 - log2(W_18_579);
X_18_589 = W_18_589*E_18_589 - log2(W_18_589);
X_18_679 = W_18_679*E_18_679 - log2(W_18_679);
X_18_689 = W_18_689*E_18_689 - log2(W_18_689);
X_18_789 = W_18_789*E_18_789 - log2(W_18_789);
X_18 = W_18*E_18 - log2(W_18);


%Objective Function
object_func = r_g;

%Optimisation
maximize(object_func)

%Constraints    
constraints(1) = r_1_123456789 + r_1_123 + r_1_124 + r_1_125 + r_1_126 + r_1_127 + r_1_128 + r_1_129 + ...
                 r_1_134 + r_1_135 + r_1_136 + r_1_137 + r_1_138 + r_1_139 + r_1_145 + r_1_146 + ...
                 r_1_147 + r_1_148 + r_1_149 + r_1_156 + r_1_157 + r_1_158 + r_1_159 + r_1_167 + ...
                 r_1_168 + r_1_169 + r_1_178 + r_1_179 + r_1_189 + r_1 - r_g;
             
constraints(2) = r_2_123456789 + r_2_123 + r_2_124 + r_2_125 + r_2_126 + r_2_127 + r_2_128 + r_2_129 + ...
                 r_2_234 + r_2_235 + r_2_236 + r_2_237 + r_2_238 + r_2_239 + r_2_245 + r_2_246 + ...
                 r_2_247 + r_2_248 + r_2_249 + r_2_256 + r_2_257 + r_2_258 + r_2_259 + r_2_267 + ...
                 r_2_268 + r_2_269 + r_2_278 + r_2_279 + r_2_289 + r_2 - r_g;
             
constraints(3) = r_3_123456789 + r_3_123 + r_3_134 + r_3_135 + r_3_136 + r_3_137 + r_3_138 + r_3_139 + ...
                 r_3_234 + r_3_235 + r_3_236 + r_3_237 + r_3_238 + r_3_239 + r_3_345 + r_3_346 + ...
                 r_3_347 + r_3_348 + r_3_349 + r_3_356 + r_3_357 + r_3_358 + r_3_359 + r_3_367 + ...
                 r_3_368 + r_3_369 + r_3_378 + r_3_379 + r_3_389 + r_3 - r_g;

constraints(4) = r_4_123456789 + r_4_124 + r_4_134 + r_4_145 + r_4_146 + r_4_147 + r_4_148 + r_4_149 + ...
                 r_4_234 + r_4_245 + r_4_246 + r_4_247 + r_4_248 + r_4_249 + r_4_345 + r_4_346 + ...
                 r_4_347 + r_4_348 + r_4_349 + r_4_456 + r_4_457 + r_4_458 + r_4_459 + r_4_467 + ...
                 r_4_468 + r_4_469 + r_4_478 + r_4_479 + r_4_489 + r_4 - r_g;

constraints(5) = r_5_123456789 + r_5_125 + r_5_135 + r_5_145 + r_5_156 + r_5_157 + r_5_158 + r_5_159 + ...
                 r_5_235 + r_5_245 + r_5_256 + r_5_257 + r_5_258 + r_5_259 + r_5_345 + r_5_356 + ...
                 r_5_357 + r_5_358 + r_5_359 + r_5_456 + r_5_457 + r_5_458 + r_5_459 + r_5_567 + ...
                 r_5_568 + r_5_569 + r_5_578 + r_5_579 + r_5_589 + r_5 - r_g;

constraints(6) = r_6_123456789 + r_6_126 + r_6_136 + r_6_146 + r_6_156 + r_6_167 + r_6_168 + r_6_169 + ...
                 r_6_236 + r_6_246 + r_6_256 + r_6_267 + r_6_268 + r_6_269 + r_6_346 + r_6_356 + ...
                 r_6_367 + r_6_368 + r_6_369 + r_6_456 + r_6_467 + r_6_468 + r_6_469 + r_6_567 + ...
                 r_6_568 + r_6_569 + r_6_678 + r_6_679 + r_6_689 + r_6 - r_g;

constraints(7) = r_7_123456789 + r_7_127 + r_7_137 + r_7_147 + r_7_157 + r_7_167 + r_7_178 + r_7_179 + ...
                 r_7_237 + r_7_247 + r_7_257 + r_7_267 + r_7_278 + r_7_279 + r_7_347 + r_7_357 + ...
                 r_7_367 + r_7_378 + r_7_379 + r_7_457 + r_7_467 + r_7_478 + r_7_479 + r_7_567 + ...
                 r_7_578 + r_7_579 + r_7_678 + r_7_679 + r_7_789 + r_7 - r_g;

constraints(8) = r_8_123456789 + r_8_128 + r_8_138 + r_8_148 + r_8_158 + r_8_168 + r_8_178 + r_8_189 + ...
                 r_8_238 + r_8_248 + r_8_258 + r_8_268 + r_8_278 + r_8_289 + r_8_348 + r_8_358 + ...
                 r_8_368 + r_8_378 + r_8_389 + r_8_458 + r_8_468 + r_8_478 + r_8_489 + r_8_568 + ...
                 r_8_578 + r_8_589 + r_8_678 + r_8_689 + r_8_789 + r_8 - r_g;

constraints(9) = r_9_123456789 + r_9_129 + r_9_139 + r_9_149 + r_9_159 + r_9_169 + r_9_179 + r_9_189 + ...
                 r_9_239 + r_9_249 + r_9_259 + r_9_269 + r_9_279 + r_9_289 + r_9_349 + r_9_359 + ...
                 r_9_369 + r_9_379 + r_9_389 + r_9_459 + r_9_469 + r_9_479 + r_9_489 + r_9_569 + ...
                 r_9_579 + r_9_589 + r_9_679 + r_9_689 + r_9_789 + r_9 - r_g;

constraints(10) = 1 - X_1 - r_1;
constraints(11) = 1 - X_2 - r_1;
constraints(12) = 1 - X_3 - r_2;
constraints(13) = 1 - X_4 - r_2;
constraints(14) = 1 - X_5 - r_3;
constraints(15) = 1 - X_6 - r_3;
constraints(16) = 1 - X_7 - r_4;
constraints(17) = 1 - X_8 - r_4;
constraints(18) = 1 - X_9 - r_5;
constraints(19) = 1 - X_10 - r_5;
constraints(20) = 1 - X_11 - r_6;
constraints(21) = 1 - X_12 - r_6;
constraints(22) = 1 - X_13 - r_7;
constraints(23) = 1 - X_14 - r_7;
constraints(24) = 1 - X_15 - r_8;
constraints(25) = 1 - X_16 - r_8;
constraints(26) = 1 - X_17 - r_9;
constraints(27) = 1 - X_18 - r_9;             
             
constraints(28) = 1 - X_1_123 - r_1_123 - r_2_123 - r_3_123;
constraints(29) = 1 - X_2_123 - r_1_123 - r_2_123 - r_3_123;
constraints(30) = 1 - X_3_123 - r_1_123 - r_2_123 - r_3_123;
constraints(31) = 1 - X_4_123 - r_1_123 - r_2_123 - r_3_123;
constraints(32) = 1 - X_5_123 - r_1_123 - r_2_123 - r_3_123;
constraints(33) = 1 - X_6_123 - r_1_123 - r_2_123 - r_3_123;

constraints(34) = 1 - X_1_124 - r_1_124 - r_2_124 - r_4_124;
constraints(35) = 1 - X_2_124 - r_1_124 - r_2_124 - r_4_124;
constraints(36) = 1 - X_3_124 - r_1_124 - r_2_124 - r_4_124;
constraints(37) = 1 - X_4_124 - r_1_124 - r_2_124 - r_4_124;
constraints(38) = 1 - X_7_124 - r_1_124 - r_2_124 - r_4_124;
constraints(39) = 1 - X_8_124 - r_1_124 - r_2_124 - r_4_124;

constraints(40) = 1 - X_1_125 - r_1_125 - r_2_125 - r_5_125;
constraints(41) = 1 - X_2_125 - r_1_125 - r_2_125 - r_5_125;
constraints(42) = 1 - X_3_125 - r_1_125 - r_2_125 - r_5_125;
constraints(43) = 1 - X_4_125 - r_1_125 - r_2_125 - r_5_125;
constraints(44) = 1 - X_9_125 - r_1_125 - r_2_125 - r_5_125;
constraints(45) = 1 - X_10_125 - r_1_125 - r_2_125 - r_5_125;

constraints(46) = 1 - X_1_126 - r_1_126 - r_2_126 - r_6_126;
constraints(47) = 1 - X_2_126 - r_1_126 - r_2_126 - r_6_126;
constraints(48) = 1 - X_3_126 - r_1_126 - r_2_126 - r_6_126;
constraints(49) = 1 - X_4_126 - r_1_126 - r_2_126 - r_6_126;
constraints(50) = 1 - X_11_126 - r_1_126 - r_2_126 - r_6_126;
constraints(51) = 1 - X_12_126 - r_1_126 - r_2_126 - r_6_126;

constraints(52) = 1 - X_1_127 - r_1_127 - r_2_127 - r_7_127;
constraints(53) = 1 - X_2_127 - r_1_127 - r_2_127 - r_7_127;
constraints(54) = 1 - X_3_127 - r_1_127 - r_2_127 - r_7_127;
constraints(55) = 1 - X_4_127 - r_1_127 - r_2_127 - r_7_127;
constraints(56) = 1 - X_13_127 - r_1_127 - r_2_127 - r_7_127;
constraints(57) = 1 - X_14_127 - r_1_127 - r_2_127 - r_7_127;

constraints(58) = 1 - X_1_128 - r_1_128 - r_2_128 - r_8_128;
constraints(59) = 1 - X_2_128 - r_1_128 - r_2_128 - r_8_128;
constraints(60) = 1 - X_3_128 - r_1_128 - r_2_128 - r_8_128;
constraints(61) = 1 - X_4_128 - r_1_128 - r_2_128 - r_8_128;
constraints(62) = 1 - X_15_128 - r_1_128 - r_2_128 - r_8_128;
constraints(63) = 1 - X_16_128 - r_1_128 - r_2_128 - r_8_128;

constraints(64) = 1 - X_1_129 - r_1_129 - r_2_129 - r_9_129;
constraints(65) = 1 - X_2_129 - r_1_129 - r_2_129 - r_9_129;
constraints(66) = 1 - X_3_129 - r_1_129 - r_2_129 - r_9_129;
constraints(67) = 1 - X_4_129 - r_1_129 - r_2_129 - r_9_129;
constraints(68) = 1 - X_17_129 - r_1_129 - r_2_129 - r_9_129;
constraints(69) = 1 - X_18_129 - r_1_129 - r_2_129 - r_9_129;

constraints(70) = 1 - X_1_134 - r_1_134 - r_3_134 - r_4_134;
constraints(71) = 1 - X_2_134 - r_1_134 - r_3_134 - r_4_134;
constraints(72) = 1 - X_5_134 - r_1_134 - r_3_134 - r_4_134;
constraints(73) = 1 - X_6_134 - r_1_134 - r_3_134 - r_4_134;
constraints(74) = 1 - X_7_134 - r_1_134 - r_3_134 - r_4_134;
constraints(75) = 1 - X_8_134 - r_1_134 - r_3_134 - r_4_134;

constraints(76) = 1 - X_1_135 - r_1_135 - r_3_135 - r_5_135;			
constraints(77) = 1 - X_2_135 - r_1_135 - r_3_135 - r_5_135;			
constraints(78) = 1 - X_5_135 - r_1_135 - r_3_135 - r_5_135;			
constraints(79) = 1 - X_6_135 - r_1_135 - r_3_135 - r_5_135;			
constraints(80) = 1 - X_9_135 - r_1_135 - r_3_135 - r_5_135;			
constraints(81) = 1 - X_10_135 - r_1_135 - r_3_135 - r_5_135;	

constraints(82) = 1 - X_1_136 - r_1_136 - r_3_136 - r_6_136;			
constraints(83) = 1 - X_2_136 - r_1_136 - r_3_136 - r_6_136;			
constraints(84) = 1 - X_5_136 - r_1_136 - r_3_136 - r_6_136;			
constraints(85) = 1 - X_6_136 - r_1_136 - r_3_136 - r_6_136;			
constraints(86) = 1 - X_11_136 - r_1_136 - r_3_136 - r_6_136;			
constraints(87) = 1 - X_12_136 - r_1_136 - r_3_136 - r_6_136;			

constraints(88) = 1 - X_1_137 - r_1_137 - r_3_137 - r_7_137;			
constraints(89) = 1 - X_2_137 - r_1_137 - r_3_137 - r_7_137;			
constraints(90) = 1 - X_5_137 - r_1_137 - r_3_137 - r_7_137;			
constraints(91) = 1 - X_6_137 - r_1_137 - r_3_137 - r_7_137;			
constraints(92) = 1 - X_13_137 - r_1_137 - r_3_137 - r_7_137;			
constraints(93) = 1 - X_14_137 - r_1_137 - r_3_137 - r_7_137;			

constraints(94) = 1 - X_1_138 - r_1_138 - r_3_138 - r_8_138;			
constraints(95) = 1 - X_2_138 - r_1_138 - r_3_138 - r_8_138;			
constraints(96) = 1 - X_5_138 - r_1_138 - r_3_138 - r_8_138;			
constraints(97) = 1 - X_6_138 - r_1_138 - r_3_138 - r_8_138;			
constraints(98) = 1 - X_15_138 - r_1_138 - r_3_138 - r_8_138;			
constraints(99) = 1 - X_16_138 - r_1_138 - r_3_138 - r_8_138;	

constraints(100) = 1 - X_1_139 - r_1_139 - r_3_139 - r_9_139;			
constraints(101) = 1 - X_2_139 - r_1_139 - r_3_139 - r_9_139;			
constraints(102) = 1 - X_5_139 - r_1_139 - r_3_139 - r_9_139;			
constraints(103) = 1 - X_6_139 - r_1_139 - r_3_139 - r_9_139;			
constraints(104) = 1 - X_17_139 - r_1_139 - r_3_139 - r_9_139;			
constraints(105) = 1 - X_18_139 - r_1_139 - r_3_139 - r_9_139;			

constraints(106) = 1 - X_1_145 - r_1_145 - r_4_145 - r_5_145;			
constraints(107) = 1 - X_2_145 - r_1_145 - r_4_145 - r_5_145;			
constraints(108) = 1 - X_7_145 - r_1_145 - r_4_145 - r_5_145;			
constraints(109) = 1 - X_8_145 - r_1_145 - r_4_145 - r_5_145;			
constraints(110) = 1 - X_9_145 - r_1_145 - r_4_145 - r_5_145;			
constraints(111) = 1 - X_10_145 - r_1_145 - r_4_145 - r_5_145;		

constraints(112) = 1 - X_1_146 - r_1_146 - r_4_146 - r_6_146;			
constraints(113) = 1 - X_2_146 - r_1_146 - r_4_146 - r_6_146;			
constraints(114) = 1 - X_7_146 - r_1_146 - r_4_146 - r_6_146;			
constraints(115) = 1 - X_8_146 - r_1_146 - r_4_146 - r_6_146;			
constraints(116) = 1 - X_11_146 - r_1_146 - r_4_146 - r_6_146;			
constraints(117) = 1 - X_12_146 - r_1_146 - r_4_146 - r_6_146;		

constraints(118) = 1 - X_1_147 - r_1_147 - r_4_147 - r_7_147;
constraints(119) = 1 - X_2_147 - r_1_147 - r_4_147 - r_7_147;
constraints(120) = 1 - X_7_147 - r_1_147 - r_4_147 - r_7_147;
constraints(121) = 1 - X_8_147 - r_1_147 - r_4_147 - r_7_147;
constraints(122) = 1 - X_13_147 - r_1_147 - r_4_147 - r_7_147;
constraints(123) = 1 - X_14_147 - r_1_147 - r_4_147 - r_7_147;

constraints(124) = 1 - X_1_148 - r_1_148 - r_4_148 - r_8_148;
constraints(125) = 1 - X_2_148 - r_1_148 - r_4_148 - r_8_148;
constraints(126) = 1 - X_7_148 - r_1_148 - r_4_148 - r_8_148;
constraints(127) = 1 - X_8_148 - r_1_148 - r_4_148 - r_8_148;
constraints(128) = 1 - X_15_148 - r_1_148 - r_4_148 - r_8_148;
constraints(129) = 1 - X_16_148 - r_1_148 - r_4_148 - r_8_148;

constraints(130) = 1 - X_1_149 - r_1_149 - r_4_149 - r_9_149;
constraints(131) = 1 - X_2_149 - r_1_149 - r_4_149 - r_9_149;
constraints(132) = 1 - X_7_149 - r_1_149 - r_4_149 - r_9_149;
constraints(133) = 1 - X_8_149 - r_1_149 - r_4_149 - r_9_149;
constraints(134) = 1 - X_17_149 - r_1_149 - r_4_149 - r_9_149;
constraints(135) = 1 - X_18_149 - r_1_149 - r_4_149 - r_9_149;

constraints(136) = 1 - X_1_156 - r_1_156 - r_5_156 - r_6_156;
constraints(137) = 1 - X_2_156 - r_1_156 - r_5_156 - r_6_156;
constraints(138) = 1 - X_9_156 - r_1_156 - r_5_156 - r_6_156;
constraints(139) = 1 - X_10_156 - r_1_156 - r_5_156 - r_6_156;
constraints(140) = 1 - X_11_156 - r_1_156 - r_5_156 - r_6_156;
constraints(141) = 1 - X_12_156 - r_1_156 - r_5_156 - r_6_156;

constraints(142) = 1 - X_1_157 - r_1_157 - r_5_157 - r_7_157;
constraints(143) = 1 - X_2_157 - r_1_157 - r_5_157 - r_7_157;
constraints(144) = 1 - X_9_157 - r_1_157 - r_5_157 - r_7_157;
constraints(145) = 1 - X_10_157 - r_1_157 - r_5_157 - r_7_157;
constraints(146) = 1 - X_13_157 - r_1_157 - r_5_157 - r_7_157;
constraints(147) = 1 - X_14_157 - r_1_157 - r_5_157 - r_7_157;

constraints(148) = 1 - X_1_158 - r_1_158 - r_5_158 - r_8_158;
constraints(149) = 1 - X_2_158 - r_1_158 - r_5_158 - r_8_158;
constraints(150) = 1 - X_9_158 - r_1_158 - r_5_158 - r_8_158;
constraints(151) = 1 - X_10_158 - r_1_158 - r_5_158 - r_8_158;
constraints(152) = 1 - X_15_158 - r_1_158 - r_5_158 - r_8_158;
constraints(153) = 1 - X_16_158 - r_1_158 - r_5_158 - r_8_158;

constraints(154) = 1 - X_1_159 - r_1_159 - r_5_159 - r_9_159;
constraints(155) = 1 - X_2_159 - r_1_159 - r_5_159 - r_9_159;
constraints(156) = 1 - X_9_159 - r_1_159 - r_5_159 - r_9_159;
constraints(157) = 1 - X_10_159 - r_1_159 - r_5_159 - r_9_159;
constraints(158) = 1 - X_17_159 - r_1_159 - r_5_159 - r_9_159;
constraints(159) = 1 - X_18_159 - r_1_159 - r_5_159 - r_9_159;

constraints(160) = 1 - X_1_167 - r_1_167 - r_6_167 - r_7_167;
constraints(161) = 1 - X_2_167 - r_1_167 - r_6_167 - r_7_167;
constraints(162) = 1 - X_11_167 - r_1_167 - r_6_167 - r_7_167;
constraints(163) = 1 - X_12_167 - r_1_167 - r_6_167 - r_7_167;
constraints(164) = 1 - X_13_167 - r_1_167 - r_6_167 - r_7_167;
constraints(165) = 1 - X_14_167 - r_1_167 - r_6_167 - r_7_167;

constraints(166) = 1 - X_1_168 - r_1_168 - r_6_168 - r_8_168;
constraints(167) = 1 - X_2_168 - r_1_168 - r_6_168 - r_8_168;
constraints(168) = 1 - X_11_168 - r_1_168 - r_6_168 - r_8_168;
constraints(169) = 1 - X_12_168 - r_1_168 - r_6_168 - r_8_168;
constraints(170) = 1 - X_15_168 - r_1_168 - r_6_168 - r_8_168;
constraints(171) = 1 - X_16_168 - r_1_168 - r_6_168 - r_8_168;

constraints(172) = 1 - X_1_169 - r_1_169 - r_6_169 - r_9_169;
constraints(173) = 1 - X_2_169 - r_1_169 - r_6_169 - r_9_169;
constraints(174) = 1 - X_11_169 - r_1_169 - r_6_169 - r_9_169;
constraints(175) = 1 - X_12_169 - r_1_169 - r_6_169 - r_9_169;
constraints(176) = 1 - X_17_169 - r_1_169 - r_6_169 - r_9_169;
constraints(177) = 1 - X_18_169 - r_1_169 - r_6_169 - r_9_169;

constraints(178) = 1 - X_1_178 - r_1_178 - r_7_178 - r_8_178;
constraints(179) = 1 - X_2_178 - r_1_178 - r_7_178 - r_8_178;
constraints(180) = 1 - X_13_178 - r_1_178 - r_7_178 - r_8_178;
constraints(181) = 1 - X_14_178 - r_1_178 - r_7_178 - r_8_178;
constraints(182) = 1 - X_15_178 - r_1_178 - r_7_178 - r_8_178;
constraints(183) = 1 - X_16_178 - r_1_178 - r_7_178 - r_8_178;

constraints(184) = 1 - X_1_179 - r_1_179 - r_7_179 - r_9_179;
constraints(185) = 1 - X_2_179 - r_1_179 - r_7_179 - r_9_179;
constraints(186) = 1 - X_13_179 - r_1_179 - r_7_179 - r_9_179;
constraints(187) = 1 - X_14_179 - r_1_179 - r_7_179 - r_9_179;
constraints(188) = 1 - X_17_179 - r_1_179 - r_7_179 - r_9_179;
constraints(189) = 1 - X_18_179 - r_1_179 - r_7_179 - r_9_179;

constraints(190) = 1 - X_1_189 - r_1_189 - r_8_189 - r_9_189;
constraints(191) = 1 - X_2_189 - r_1_189 - r_8_189 - r_9_189;
constraints(192) = 1 - X_15_189 - r_1_189 - r_8_189 - r_9_189;
constraints(193) = 1 - X_16_189 - r_1_189 - r_8_189 - r_9_189;
constraints(194) = 1 - X_17_189 - r_1_189 - r_8_189 - r_9_189;
constraints(195) = 1 - X_18_189 - r_1_189 - r_8_189 - r_9_189;

constraints(196) = 1 - X_3_234 - r_2_234 - r_3_234 - r_4_234;
constraints(197) = 1 - X_4_234 - r_2_234 - r_3_234 - r_4_234;
constraints(198) = 1 - X_5_234 - r_2_234 - r_3_234 - r_4_234;
constraints(199) = 1 - X_6_234 - r_2_234 - r_3_234 - r_4_234;
constraints(200) = 1 - X_7_234 - r_2_234 - r_3_234 - r_4_234;
constraints(201) = 1 - X_8_234 - r_2_234 - r_3_234 - r_4_234;

constraints(202) = 1 - X_3_235 - r_2_235 - r_3_235 - r_5_235;
constraints(203) = 1 - X_4_235 - r_2_235 - r_3_235 - r_5_235;
constraints(204) = 1 - X_5_235 - r_2_235 - r_3_235 - r_5_235;
constraints(205) = 1 - X_6_235 - r_2_235 - r_3_235 - r_5_235;
constraints(206) = 1 - X_9_235 - r_2_235 - r_3_235 - r_5_235;
constraints(207) = 1 - X_10_235 - r_2_235 - r_3_235 - r_5_235;

constraints(208) = 1 - X_3_236 - r_2_236 - r_3_236 - r_6_236;
constraints(209) = 1 - X_4_236 - r_2_236 - r_3_236 - r_6_236;
constraints(210) = 1 - X_5_236 - r_2_236 - r_3_236 - r_6_236;
constraints(211) = 1 - X_6_236 - r_2_236 - r_3_236 - r_6_236;
constraints(212) = 1 - X_11_236 - r_2_236 - r_3_236 - r_6_236;
constraints(213) = 1 - X_12_236 - r_2_236 - r_3_236 - r_6_236;

constraints(214) = 1 - X_3_237 - r_2_237 - r_3_237 - r_7_237;
constraints(215) = 1 - X_4_237 - r_2_237 - r_3_237 - r_7_237;
constraints(216) = 1 - X_5_237 - r_2_237 - r_3_237 - r_7_237;
constraints(217) = 1 - X_6_237 - r_2_237 - r_3_237 - r_7_237;
constraints(218) = 1 - X_13_237 - r_2_237 - r_3_237 - r_7_237;
constraints(219) = 1 - X_14_237 - r_2_237 - r_3_237 - r_7_237;

constraints(220) = 1 - X_3_238 - r_2_238 - r_3_238 - r_8_238;
constraints(221) = 1 - X_4_238 - r_2_238 - r_3_238 - r_8_238;
constraints(222) = 1 - X_5_238 - r_2_238 - r_3_238 - r_8_238;
constraints(223) = 1 - X_6_238 - r_2_238 - r_3_238 - r_8_238;
constraints(224) = 1 - X_15_238 - r_2_238 - r_3_238 - r_8_238;
constraints(225) = 1 - X_16_238 - r_2_238 - r_3_238 - r_8_238;

constraints(226) = 1 - X_3_239 - r_2_239 - r_3_239 - r_9_239;
constraints(227) = 1 - X_4_239 - r_2_239 - r_3_239 - r_9_239;
constraints(228) = 1 - X_5_239 - r_2_239 - r_3_239 - r_9_239;
constraints(229) = 1 - X_6_239 - r_2_239 - r_3_239 - r_9_239;
constraints(230) = 1 - X_17_239 - r_2_239 - r_3_239 - r_9_239;
constraints(231) = 1 - X_18_239 - r_2_239 - r_3_239 - r_9_239;

constraints(232) = 1 - X_3_245 - r_2_245 - r_4_245 - r_5_245;
constraints(233) = 1 - X_4_245 - r_2_245 - r_4_245 - r_5_245;
constraints(234) = 1 - X_7_245 - r_2_245 - r_4_245 - r_5_245;
constraints(235) = 1 - X_8_245 - r_2_245 - r_4_245 - r_5_245;
constraints(236) = 1 - X_9_245 - r_2_245 - r_4_245 - r_5_245;
constraints(237) = 1 - X_10_245 - r_2_245 - r_4_245 - r_5_245;

constraints(238) = 1 - X_3_246 - r_2_246 - r_4_246 - r_6_246;
constraints(239) = 1 - X_4_246 - r_2_246 - r_4_246 - r_6_246;
constraints(240) = 1 - X_7_246 - r_2_246 - r_4_246 - r_6_246;
constraints(241) = 1 - X_8_246 - r_2_246 - r_4_246 - r_6_246;
constraints(242) = 1 - X_11_246 - r_2_246 - r_4_246 - r_6_246;
constraints(243) = 1 - X_12_246 - r_2_246 - r_4_246 - r_6_246;

constraints(244) = 1 - X_3_247 - r_2_247 - r_4_247 - r_7_247;
constraints(245) = 1 - X_4_247 - r_2_247 - r_4_247 - r_7_247;
constraints(246) = 1 - X_7_247 - r_2_247 - r_4_247 - r_7_247;
constraints(247) = 1 - X_8_247 - r_2_247 - r_4_247 - r_7_247;
constraints(248) = 1 - X_13_247 - r_2_247 - r_4_247 - r_7_247;
constraints(249) = 1 - X_14_247 - r_2_247 - r_4_247 - r_7_247;

constraints(250) = 1 - X_3_248 - r_2_248 - r_4_248 - r_8_248;
constraints(251) = 1 - X_4_248 - r_2_248 - r_4_248 - r_8_248;
constraints(252) = 1 - X_7_248 - r_2_248 - r_4_248 - r_8_248;
constraints(253) = 1 - X_8_248 - r_2_248 - r_4_248 - r_8_248;
constraints(254) = 1 - X_15_248 - r_2_248 - r_4_248 - r_8_248;
constraints(255) = 1 - X_16_248 - r_2_248 - r_4_248 - r_8_248;

constraints(256) = 1 - X_3_249 - r_2_249 - r_4_249 - r_9_249;
constraints(257) = 1 - X_4_249 - r_2_249 - r_4_249 - r_9_249;
constraints(258) = 1 - X_7_249 - r_2_249 - r_4_249 - r_9_249;
constraints(259) = 1 - X_8_249 - r_2_249 - r_4_249 - r_9_249;
constraints(260) = 1 - X_17_249 - r_2_249 - r_4_249 - r_9_249;
constraints(261) = 1 - X_18_249 - r_2_249 - r_4_249 - r_9_249;

constraints(262) = 1 - X_3_256 - r_2_256 - r_5_256 - r_6_256;
constraints(263) = 1 - X_4_256 - r_2_256 - r_5_256 - r_6_256;
constraints(264) = 1 - X_9_256 - r_2_256 - r_5_256 - r_6_256;
constraints(265) = 1 - X_10_256 - r_2_256 - r_5_256 - r_6_256;
constraints(266) = 1 - X_11_256 - r_2_256 - r_5_256 - r_6_256;
constraints(267) = 1 - X_12_256 - r_2_256 - r_5_256 - r_6_256;

constraints(268) = 1 - X_3_257 - r_2_257 - r_5_257 - r_7_257;
constraints(269) = 1 - X_4_257 - r_2_257 - r_5_257 - r_7_257;
constraints(270) = 1 - X_9_257 - r_2_257 - r_5_257 - r_7_257;
constraints(271) = 1 - X_10_257 - r_2_257 - r_5_257 - r_7_257;
constraints(272) = 1 - X_13_257 - r_2_257 - r_5_257 - r_7_257;
constraints(273) = 1 - X_14_257 - r_2_257 - r_5_257 - r_7_257;

constraints(274) = 1 - X_3_258 - r_2_258 - r_5_258 - r_8_258;
constraints(275) = 1 - X_4_258 - r_2_258 - r_5_258 - r_8_258;
constraints(276) = 1 - X_9_258 - r_2_258 - r_5_258 - r_8_258;
constraints(277) = 1 - X_10_258 - r_2_258 - r_5_258 - r_8_258;
constraints(278) = 1 - X_15_258 - r_2_258 - r_5_258 - r_8_258;
constraints(279) = 1 - X_16_258 - r_2_258 - r_5_258 - r_8_258;

constraints(280) = 1 - X_3_259 - r_2_259 - r_5_259 - r_9_259;
constraints(281) = 1 - X_4_259 - r_2_259 - r_5_259 - r_9_259;
constraints(282) = 1 - X_9_259 - r_2_259 - r_5_259 - r_9_259;
constraints(283) = 1 - X_10_259 - r_2_259 - r_5_259 - r_9_259;
constraints(284) = 1 - X_17_259 - r_2_259 - r_5_259 - r_9_259;
constraints(285) = 1 - X_18_259 - r_2_259 - r_5_259 - r_9_259;

constraints(286) = 1 - X_3_267 - r_2_267 - r_6_267 - r_7_267;
constraints(287) = 1 - X_4_267 - r_2_267 - r_6_267 - r_7_267;
constraints(288) = 1 - X_11_267 - r_2_267 - r_6_267 - r_7_267;
constraints(289) = 1 - X_12_267 - r_2_267 - r_6_267 - r_7_267;
constraints(290) = 1 - X_13_267 - r_2_267 - r_6_267 - r_7_267;
constraints(291) = 1 - X_14_267 - r_2_267 - r_6_267 - r_7_267;

constraints(292) = 1 - X_3_268 - r_2_268 - r_6_268 - r_8_268;
constraints(293) = 1 - X_4_268 - r_2_268 - r_6_268 - r_8_268;
constraints(294) = 1 - X_11_268 - r_2_268 - r_6_268 - r_8_268;
constraints(295) = 1 - X_12_268 - r_2_268 - r_6_268 - r_8_268;
constraints(296) = 1 - X_15_268 - r_2_268 - r_6_268 - r_8_268;
constraints(297) = 1 - X_16_268 - r_2_268 - r_6_268 - r_8_268;

constraints(298) = 1 - X_3_269 - r_2_269 - r_6_269 - r_9_269;
constraints(299) = 1 - X_4_269 - r_2_269 - r_6_269 - r_9_269;
constraints(300) = 1 - X_11_269 - r_2_269 - r_6_269 - r_9_269;
constraints(301) = 1 - X_12_269 - r_2_269 - r_6_269 - r_9_269;
constraints(302) = 1 - X_17_269 - r_2_269 - r_6_269 - r_9_269;
constraints(303) = 1 - X_18_269 - r_2_269 - r_6_269 - r_9_269;

constraints(304) = 1 - X_3_278 - r_2_278 - r_7_278 - r_8_278;
constraints(305) = 1 - X_4_278 - r_2_278 - r_7_278 - r_8_278;
constraints(306) = 1 - X_13_278 - r_2_278 - r_7_278 - r_8_278;
constraints(307) = 1 - X_14_278 - r_2_278 - r_7_278 - r_8_278;
constraints(308) = 1 - X_15_278 - r_2_278 - r_7_278 - r_8_278;
constraints(309) = 1 - X_16_278 - r_2_278 - r_7_278 - r_8_278;

constraints(310) = 1 - X_3_279 - r_2_279 - r_7_279 - r_9_279;
constraints(311) = 1 - X_4_279 - r_2_279 - r_7_279 - r_9_279;
constraints(312) = 1 - X_13_279 - r_2_279 - r_7_279 - r_9_279;
constraints(313) = 1 - X_14_279 - r_2_279 - r_7_279 - r_9_279;
constraints(314) = 1 - X_17_279 - r_2_279 - r_7_279 - r_9_279;
constraints(315) = 1 - X_18_279 - r_2_279 - r_7_279 - r_9_279;

constraints(316) = 1 - X_3_289 - r_2_289 - r_8_289 - r_9_289;
constraints(317) = 1 - X_4_289 - r_2_289 - r_8_289 - r_9_289;
constraints(318) = 1 - X_15_289 - r_2_289 - r_8_289 - r_9_289;
constraints(319) = 1 - X_16_289 - r_2_289 - r_8_289 - r_9_289;
constraints(320) = 1 - X_17_289 - r_2_289 - r_8_289 - r_9_289;
constraints(321) = 1 - X_18_289 - r_2_289 - r_8_289 - r_9_289;

constraints(322) = 1 - X_5_345 - r_3_345 - r_4_345 - r_5_345;
constraints(323) = 1 - X_6_345 - r_3_345 - r_4_345 - r_5_345;
constraints(324) = 1 - X_7_345 - r_3_345 - r_4_345 - r_5_345;
constraints(325) = 1 - X_8_345 - r_3_345 - r_4_345 - r_5_345;
constraints(326) = 1 - X_9_345 - r_3_345 - r_4_345 - r_5_345;
constraints(327) = 1 - X_10_345 - r_3_345 - r_4_345 - r_5_345;

constraints(328) = 1 - X_5_346 - r_3_346 - r_4_346 - r_6_346;
constraints(329) = 1 - X_6_346 - r_3_346 - r_4_346 - r_6_346;
constraints(330) = 1 - X_7_346 - r_3_346 - r_4_346 - r_6_346;
constraints(331) = 1 - X_8_346 - r_3_346 - r_4_346 - r_6_346;
constraints(332) = 1 - X_11_346 - r_3_346 - r_4_346 - r_6_346;
constraints(333) = 1 - X_12_346 - r_3_346 - r_4_346 - r_6_346;

constraints(334) = 1 - X_5_347 - r_3_347 - r_4_347 - r_7_347;
constraints(335) = 1 - X_6_347 - r_3_347 - r_4_347 - r_7_347;
constraints(336) = 1 - X_7_347 - r_3_347 - r_4_347 - r_7_347;
constraints(337) = 1 - X_8_347 - r_3_347 - r_4_347 - r_7_347;
constraints(338) = 1 - X_13_347 - r_3_347 - r_4_347 - r_7_347;
constraints(339) = 1 - X_14_347 - r_3_347 - r_4_347 - r_7_347;

constraints(340) = 1 - X_5_348 - r_3_348 - r_4_348 - r_8_348;
constraints(341) = 1 - X_6_348 - r_3_348 - r_4_348 - r_8_348;
constraints(342) = 1 - X_7_348 - r_3_348 - r_4_348 - r_8_348;
constraints(343) = 1 - X_8_348 - r_3_348 - r_4_348 - r_8_348;
constraints(344) = 1 - X_15_348 - r_3_348 - r_4_348 - r_8_348;
constraints(345) = 1 - X_16_348 - r_3_348 - r_4_348 - r_8_348;

constraints(346) = 1 - X_5_349 - r_3_349 - r_4_349 - r_9_349;
constraints(347) = 1 - X_6_349 - r_3_349 - r_4_349 - r_9_349;
constraints(348) = 1 - X_7_349 - r_3_349 - r_4_349 - r_9_349;
constraints(349) = 1 - X_8_349 - r_3_349 - r_4_349 - r_9_349;
constraints(350) = 1 - X_17_349 - r_3_349 - r_4_349 - r_9_349;
constraints(351) = 1 - X_18_349 - r_3_349 - r_4_349 - r_9_349;

constraints(352) = 1 - X_5_356 - r_3_356 - r_5_356 - r_6_356;
constraints(353) = 1 - X_6_356 - r_3_356 - r_5_356 - r_6_356;
constraints(354) = 1 - X_9_356 - r_3_356 - r_5_356 - r_6_356;
constraints(355) = 1 - X_10_356 - r_3_356 - r_5_356 - r_6_356;
constraints(356) = 1 - X_11_356 - r_3_356 - r_5_356 - r_6_356;
constraints(357) = 1 - X_12_356 - r_3_356 - r_5_356 - r_6_356;

constraints(358) = 1 - X_5_357 - r_3_357 - r_5_357 - r_7_357;
constraints(359) = 1 - X_6_357 - r_3_357 - r_5_357 - r_7_357;
constraints(360) = 1 - X_9_357 - r_3_357 - r_5_357 - r_7_357;
constraints(361) = 1 - X_10_357 - r_3_357 - r_5_357 - r_7_357;
constraints(362) = 1 - X_13_357 - r_3_357 - r_5_357 - r_7_357;
constraints(363) = 1 - X_14_357 - r_3_357 - r_5_357 - r_7_357;

constraints(364) = 1 - X_5_358 - r_3_358 - r_5_358 - r_8_358;
constraints(365) = 1 - X_6_358 - r_3_358 - r_5_358 - r_8_358;
constraints(366) = 1 - X_9_358 - r_3_358 - r_5_358 - r_8_358;
constraints(367) = 1 - X_10_358 - r_3_358 - r_5_358 - r_8_358;
constraints(368) = 1 - X_15_358 - r_3_358 - r_5_358 - r_8_358;
constraints(369) = 1 - X_16_358 - r_3_358 - r_5_358 - r_8_358;

constraints(370) = 1 - X_5_359 - r_3_359 - r_5_359 - r_9_359;
constraints(371) = 1 - X_6_359 - r_3_359 - r_5_359 - r_9_359;
constraints(372) = 1 - X_9_359 - r_3_359 - r_5_359 - r_9_359;
constraints(373) = 1 - X_10_359 - r_3_359 - r_5_359 - r_9_359;
constraints(374) = 1 - X_17_359 - r_3_359 - r_5_359 - r_9_359;
constraints(375) = 1 - X_18_359 - r_3_359 - r_5_359 - r_9_359;

constraints(376) = 1 - X_5_367 - r_3_367 - r_6_367 - r_7_367;
constraints(377) = 1 - X_6_367 - r_3_367 - r_6_367 - r_7_367;
constraints(378) = 1 - X_11_367 - r_3_367 - r_6_367 - r_7_367;
constraints(379) = 1 - X_12_367 - r_3_367 - r_6_367 - r_7_367;
constraints(380) = 1 - X_13_367 - r_3_367 - r_6_367 - r_7_367;
constraints(381) = 1 - X_14_367 - r_3_367 - r_6_367 - r_7_367;

constraints(382) = 1 - X_5_368 - r_3_368 - r_6_368 - r_8_368;
constraints(383) = 1 - X_6_368 - r_3_368 - r_6_368 - r_8_368;
constraints(384) = 1 - X_11_368 - r_3_368 - r_6_368 - r_8_368;
constraints(385) = 1 - X_12_368 - r_3_368 - r_6_368 - r_8_368;
constraints(386) = 1 - X_15_368 - r_3_368 - r_6_368 - r_8_368;
constraints(387) = 1 - X_16_368 - r_3_368 - r_6_368 - r_8_368;

constraints(388) = 1 - X_5_369 - r_3_369 - r_6_369 - r_9_369;
constraints(389) = 1 - X_6_369 - r_3_369 - r_6_369 - r_9_369;
constraints(390) = 1 - X_11_369 - r_3_369 - r_6_369 - r_9_369;
constraints(391) = 1 - X_12_369 - r_3_369 - r_6_369 - r_9_369;
constraints(392) = 1 - X_17_369 - r_3_369 - r_6_369 - r_9_369;
constraints(393) = 1 - X_18_369 - r_3_369 - r_6_369 - r_9_369;

constraints(394) = 1 - X_5_378 - r_3_378 - r_7_378 - r_8_378;
constraints(395) = 1 - X_6_378 - r_3_378 - r_7_378 - r_8_378;
constraints(396) = 1 - X_13_378 - r_3_378 - r_7_378 - r_8_378;
constraints(397) = 1 - X_14_378 - r_3_378 - r_7_378 - r_8_378;
constraints(398) = 1 - X_15_378 - r_3_378 - r_7_378 - r_8_378;
constraints(399) = 1 - X_16_378 - r_3_378 - r_7_378 - r_8_378;

constraints(400) = 1 - X_5_379 - r_3_379 - r_7_379 - r_9_379;
constraints(401) = 1 - X_6_379 - r_3_379 - r_7_379 - r_9_379;
constraints(402) = 1 - X_13_379 - r_3_379 - r_7_379 - r_9_379;
constraints(403) = 1 - X_14_379 - r_3_379 - r_7_379 - r_9_379;
constraints(404) = 1 - X_17_379 - r_3_379 - r_7_379 - r_9_379;
constraints(405) = 1 - X_18_379 - r_3_379 - r_7_379 - r_9_379;

constraints(406) = 1 - X_5_389 - r_3_389 - r_8_389 - r_9_389;
constraints(407) = 1 - X_6_389 - r_3_389 - r_8_389 - r_9_389;
constraints(408) = 1 - X_15_389 - r_3_389 - r_8_389 - r_9_389;
constraints(409) = 1 - X_16_389 - r_3_389 - r_8_389 - r_9_389;
constraints(410) = 1 - X_17_389 - r_3_389 - r_8_389 - r_9_389;
constraints(411) = 1 - X_18_389 - r_3_389 - r_8_389 - r_9_389;

constraints(412) = 1 - X_7_456 - r_4_456 - r_5_456 - r_6_456;
constraints(413) = 1 - X_8_456 - r_4_456 - r_5_456 - r_6_456;
constraints(414) = 1 - X_9_456 - r_4_456 - r_5_456 - r_6_456;
constraints(415) = 1 - X_10_456 - r_4_456 - r_5_456 - r_6_456;
constraints(416) = 1 - X_11_456 - r_4_456 - r_5_456 - r_6_456;
constraints(417) = 1 - X_12_456 - r_4_456 - r_5_456 - r_6_456;

constraints(418) = 1 - X_7_457 - r_4_457 - r_5_457 - r_7_457;
constraints(419) = 1 - X_8_457 - r_4_457 - r_5_457 - r_7_457;
constraints(420) = 1 - X_9_457 - r_4_457 - r_5_457 - r_7_457;
constraints(421) = 1 - X_10_457 - r_4_457 - r_5_457 - r_7_457;
constraints(422) = 1 - X_13_457 - r_4_457 - r_5_457 - r_7_457;
constraints(423) = 1 - X_14_457 - r_4_457 - r_5_457 - r_7_457;

constraints(424) = 1 - X_7_458 - r_4_458 - r_5_458 - r_8_458;
constraints(425) = 1 - X_8_458 - r_4_458 - r_5_458 - r_8_458;
constraints(426) = 1 - X_9_458 - r_4_458 - r_5_458 - r_8_458;
constraints(427) = 1 - X_10_458 - r_4_458 - r_5_458 - r_8_458;
constraints(428) = 1 - X_15_458 - r_4_458 - r_5_458 - r_8_458;
constraints(429) = 1 - X_16_458 - r_4_458 - r_5_458 - r_8_458;

constraints(430) = 1 - X_7_459 - r_4_459 - r_5_459 - r_9_459;
constraints(431) = 1 - X_8_459 - r_4_459 - r_5_459 - r_9_459;
constraints(432) = 1 - X_9_459 - r_4_459 - r_5_459 - r_9_459;
constraints(433) = 1 - X_10_459 - r_4_459 - r_5_459 - r_9_459;
constraints(434) = 1 - X_17_459 - r_4_459 - r_5_459 - r_9_459;
constraints(435) = 1 - X_18_459 - r_4_459 - r_5_459 - r_9_459;

constraints(436) = 1 - X_7_467 - r_4_467 - r_6_467 - r_7_467;
constraints(437) = 1 - X_8_467 - r_4_467 - r_6_467 - r_7_467;
constraints(438) = 1 - X_11_467 - r_4_467 - r_6_467 - r_7_467;
constraints(439) = 1 - X_12_467 - r_4_467 - r_6_467 - r_7_467;
constraints(440) = 1 - X_13_467 - r_4_467 - r_6_467 - r_7_467;
constraints(441) = 1 - X_14_467 - r_4_467 - r_6_467 - r_7_467;

constraints(442) = 1 - X_7_468 - r_4_468 - r_6_468 - r_8_468;
constraints(443) = 1 - X_8_468 - r_4_468 - r_6_468 - r_8_468;
constraints(444) = 1 - X_11_468 - r_4_468 - r_6_468 - r_8_468;
constraints(445) = 1 - X_12_468 - r_4_468 - r_6_468 - r_8_468;
constraints(446) = 1 - X_15_468 - r_4_468 - r_6_468 - r_8_468;
constraints(447) = 1 - X_16_468 - r_4_468 - r_6_468 - r_8_468;

constraints(448) = 1 - X_7_469 - r_4_469 - r_6_469 - r_9_469;
constraints(449) = 1 - X_8_469 - r_4_469 - r_6_469 - r_9_469;
constraints(450) = 1 - X_11_469 - r_4_469 - r_6_469 - r_9_469;
constraints(451) = 1 - X_12_469 - r_4_469 - r_6_469 - r_9_469;
constraints(452) = 1 - X_17_469 - r_4_469 - r_6_469 - r_9_469;
constraints(453) = 1 - X_18_469 - r_4_469 - r_6_469 - r_9_469;

constraints(454) = 1 - X_7_478 - r_4_478 - r_7_478 - r_8_478;
constraints(455) = 1 - X_8_478 - r_4_478 - r_7_478 - r_8_478;
constraints(456) = 1 - X_13_478 - r_4_478 - r_7_478 - r_8_478;
constraints(457) = 1 - X_14_478 - r_4_478 - r_7_478 - r_8_478;
constraints(458) = 1 - X_15_478 - r_4_478 - r_7_478 - r_8_478;
constraints(459) = 1 - X_16_478 - r_4_478 - r_7_478 - r_8_478;

constraints(460) = 1 - X_7_479 - r_4_479 - r_7_479 - r_9_479;
constraints(461) = 1 - X_8_479 - r_4_479 - r_7_479 - r_9_479;
constraints(462) = 1 - X_13_479 - r_4_479 - r_7_479 - r_9_479;
constraints(463) = 1 - X_14_479 - r_4_479 - r_7_479 - r_9_479;
constraints(464) = 1 - X_17_479 - r_4_479 - r_7_479 - r_9_479;
constraints(465) = 1 - X_18_479 - r_4_479 - r_7_479 - r_9_479;

constraints(466) = 1 - X_7_489 - r_4_489 - r_8_489 - r_9_489;
constraints(467) = 1 - X_8_489 - r_4_489 - r_8_489 - r_9_489;
constraints(468) = 1 - X_15_489 - r_4_489 - r_8_489 - r_9_489;
constraints(469) = 1 - X_16_489 - r_4_489 - r_8_489 - r_9_489;
constraints(470) = 1 - X_17_489 - r_4_489 - r_8_489 - r_9_489;
constraints(471) = 1 - X_18_489 - r_4_489 - r_8_489 - r_9_489;

constraints(472) = 1 - X_9_567 - r_5_567 - r_6_567 - r_7_567;
constraints(473) = 1 - X_10_567 - r_5_567 - r_6_567 - r_7_567;
constraints(474) = 1 - X_11_567 - r_5_567 - r_6_567 - r_7_567;
constraints(475) = 1 - X_12_567 - r_5_567 - r_6_567 - r_7_567;
constraints(476) = 1 - X_13_567 - r_5_567 - r_6_567 - r_7_567;
constraints(477) = 1 - X_14_567 - r_5_567 - r_6_567 - r_7_567;

constraints(478) = 1 - X_9_568 - r_5_568 - r_6_568 - r_8_568;
constraints(479) = 1 - X_10_568 - r_5_568 - r_6_568 - r_8_568;
constraints(480) = 1 - X_11_568 - r_5_568 - r_6_568 - r_8_568;
constraints(481) = 1 - X_12_568 - r_5_568 - r_6_568 - r_8_568;
constraints(482) = 1 - X_15_568 - r_5_568 - r_6_568 - r_8_568;
constraints(483) = 1 - X_16_568 - r_5_568 - r_6_568 - r_8_568;

constraints(484) = 1 - X_9_569 - r_5_569 - r_6_569 - r_9_569;
constraints(485) = 1 - X_10_569 - r_5_569 - r_6_569 - r_9_569;
constraints(486) = 1 - X_11_569 - r_5_569 - r_6_569 - r_9_569;
constraints(487) = 1 - X_12_569 - r_5_569 - r_6_569 - r_9_569;
constraints(488) = 1 - X_17_569 - r_5_569 - r_6_569 - r_9_569;
constraints(489) = 1 - X_18_569 - r_5_569 - r_6_569 - r_9_569;

constraints(490) = 1 - X_9_578 - r_5_578 - r_7_578 - r_8_578;
constraints(491) = 1 - X_10_578 - r_5_578 - r_7_578 - r_8_578;
constraints(492) = 1 - X_13_578 - r_5_578 - r_7_578 - r_8_578;
constraints(493) = 1 - X_14_578 - r_5_578 - r_7_578 - r_8_578;
constraints(494) = 1 - X_15_578 - r_5_578 - r_7_578 - r_8_578;
constraints(495) = 1 - X_16_578 - r_5_578 - r_7_578 - r_8_578;

constraints(496) = 1 - X_9_579 - r_5_579 - r_7_579 - r_9_579;
constraints(497) = 1 - X_10_579 - r_5_579 - r_7_579 - r_9_579;
constraints(498) = 1 - X_13_579 - r_5_579 - r_7_579 - r_9_579;
constraints(499) = 1 - X_14_579 - r_5_579 - r_7_579 - r_9_579;
constraints(500) = 1 - X_17_579 - r_5_579 - r_7_579 - r_9_579;
constraints(501) = 1 - X_18_579 - r_5_579 - r_7_579 - r_9_579;

constraints(502) = 1 - X_9_589 - r_5_589 - r_8_589 - r_9_589;
constraints(503) = 1 - X_10_589 - r_5_589 - r_8_589 - r_9_589;
constraints(504) = 1 - X_15_589 - r_5_589 - r_8_589 - r_9_589;
constraints(505) = 1 - X_16_589 - r_5_589 - r_8_589 - r_9_589;
constraints(506) = 1 - X_17_589 - r_5_589 - r_8_589 - r_9_589;
constraints(507) = 1 - X_18_589 - r_5_589 - r_8_589 - r_9_589;

constraints(508) = 1 - X_11_678 - r_6_678 - r_7_678 - r_8_678;
constraints(509) = 1 - X_12_678 - r_6_678 - r_7_678 - r_8_678;
constraints(510) = 1 - X_13_678 - r_6_678 - r_7_678 - r_8_678;
constraints(511) = 1 - X_14_678 - r_6_678 - r_7_678 - r_8_678;
constraints(512) = 1 - X_15_678 - r_6_678 - r_7_678 - r_8_678;
constraints(513) = 1 - X_16_678 - r_6_678 - r_7_678 - r_8_678;

constraints(514) = 1 - X_11_679 - r_6_679 - r_7_679 - r_9_679;
constraints(515) = 1 - X_12_679 - r_6_679 - r_7_679 - r_9_679;
constraints(516) = 1 - X_13_679 - r_6_679 - r_7_679 - r_9_679;
constraints(517) = 1 - X_14_679 - r_6_679 - r_7_679 - r_9_679;
constraints(518) = 1 - X_17_679 - r_6_679 - r_7_679 - r_9_679;
constraints(519) = 1 - X_18_679 - r_6_679 - r_7_679 - r_9_679;

constraints(520) = 1 - X_11_689 - r_6_689 - r_8_689 - r_9_689;
constraints(521) = 1 - X_12_689 - r_6_689 - r_8_689 - r_9_689;
constraints(522) = 1 - X_15_689 - r_6_689 - r_8_689 - r_9_689;
constraints(523) = 1 - X_16_689 - r_6_689 - r_8_689 - r_9_689;
constraints(524) = 1 - X_17_689 - r_6_689 - r_8_689 - r_9_689;
constraints(525) = 1 - X_18_689 - r_6_689 - r_8_689 - r_9_689;

constraints(526) = 1 - X_13_789 - r_7_789 - r_8_789 - r_9_789;
constraints(527) = 1 - X_14_789 - r_7_789 - r_8_789 - r_9_789;
constraints(528) = 1 - X_15_789 - r_7_789 - r_8_789 - r_9_789;
constraints(529) = 1 - X_16_789 - r_7_789 - r_8_789 - r_9_789;
constraints(530) = 1 - X_17_789 - r_7_789 - r_8_789 - r_9_789;
constraints(531) = 1 - X_18_789 - r_7_789 - r_8_789 - r_9_789;

constraints(532) = 1 - X_1_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(533) = 1 - X_2_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(534) = 1 - X_3_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(535) = 1 - X_4_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(536) = 1 - X_5_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(537) = 1 - X_6_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(538) = 1 - X_7_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(539) = 1 - X_8_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(540) = 1 - X_9_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(541) = 1 - X_10_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(542) = 1 - X_11_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(543) = 1 - X_12_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(544) = 1 - X_13_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(545) = 1 - X_14_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(546) = 1 - X_15_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(547) = 1 - X_16_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(548) = 1 - X_17_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(549) = 1 - X_18_123456789 - r_1_123456789 - r_2_123456789 - r_3_123456789 - ...
                   r_4_123456789 - r_5_123456789 - r_6_123456789 - r_7_123456789 - r_8_123456789 - r_9_123456789;

constraints(550) = r_1_123456789;
constraints(551) = r_2_123456789;
constraints(552) = r_3_123456789;
constraints(553) = r_4_123456789;
constraints(554) = r_5_123456789;
constraints(555) = r_6_123456789;
constraints(556) = r_7_123456789;
constraints(557) = r_8_123456789;
constraints(558) = r_9_123456789;
constraints(559) = r_1_123;
constraints(560) = r_2_123;
constraints(561) = r_3_123;
constraints(562) = r_1_124;
constraints(563) = r_2_124;
constraints(564) = r_4_124;
constraints(565) = r_1_125;
constraints(566) = r_2_125;
constraints(567) = r_5_125;
constraints(568) = r_1_126;
constraints(569) = r_2_126;
constraints(570) = r_6_126;
constraints(571) = r_1_127;
constraints(572) = r_2_127;
constraints(573) = r_7_127;
constraints(574) = r_1_128;
constraints(575) = r_2_128;
constraints(576) = r_8_128;
constraints(577) = r_1_129;
constraints(578) = r_2_129;
constraints(579) = r_9_129;
constraints(580) = r_1_134;
constraints(581) = r_3_134;
constraints(582) = r_4_134;
constraints(583) = r_1_135;
constraints(584) = r_3_135;
constraints(585) = r_5_135;
constraints(586) = r_1_136;
constraints(587) = r_3_136;
constraints(588) = r_6_136;
constraints(589) = r_1_137;
constraints(590) = r_3_137;
constraints(591) = r_7_137;
constraints(592) = r_1_138;
constraints(593) = r_3_138;
constraints(594) = r_8_138;
constraints(595) = r_1_139;
constraints(596) = r_3_139;
constraints(597) = r_9_139;
constraints(598) = r_1_145;
constraints(599) = r_4_145;
constraints(600) = r_5_145;
constraints(601) = r_1_146;
constraints(602) = r_4_146;
constraints(603) = r_6_146;
constraints(604) = r_1_147;
constraints(605) = r_4_147;
constraints(606) = r_7_147;
constraints(607) = r_1_148;
constraints(608) = r_4_148;
constraints(609) = r_8_148;
constraints(610) = r_1_149;
constraints(611) = r_4_149;
constraints(612) = r_9_149;
constraints(613) = r_1_156;
constraints(614) = r_5_156;
constraints(615) = r_6_156;
constraints(616) = r_1_157;
constraints(617) = r_5_157;
constraints(618) = r_7_157;
constraints(619) = r_1_158;
constraints(620) = r_5_158;
constraints(621) = r_8_158;
constraints(622) = r_1_159;
constraints(623) = r_5_159;
constraints(624) = r_9_159;
constraints(625) = r_1_167;
constraints(626) = r_6_167;
constraints(627) = r_7_167;
constraints(628) = r_1_168;
constraints(629) = r_6_168;
constraints(630) = r_8_168;
constraints(631) = r_1_169;
constraints(632) = r_6_169;
constraints(633) = r_9_169;
constraints(634) = r_1_178;
constraints(635) = r_7_178;
constraints(636) = r_8_178;
constraints(637) = r_1_179;
constraints(638) = r_7_179;
constraints(639) = r_9_179;
constraints(640) = r_1_189;
constraints(641) = r_8_189;
constraints(642) = r_9_189;
constraints(643) = r_2_234;
constraints(644) = r_3_234;
constraints(645) = r_4_234;
constraints(646) = r_2_235;
constraints(647) = r_3_235;
constraints(648) = r_5_235;
constraints(649) = r_2_236;
constraints(650) = r_3_236;
constraints(651) = r_6_236;
constraints(652) = r_2_237;
constraints(653) = r_3_237;
constraints(654) = r_7_237;
constraints(655) = r_2_238;
constraints(656) = r_3_238;
constraints(657) = r_8_238;
constraints(658) = r_2_239;
constraints(659) = r_3_239;
constraints(660) = r_9_239;
constraints(661) = r_2_245;
constraints(662) = r_4_245;
constraints(663) = r_5_245;
constraints(664) = r_2_246;
constraints(665) = r_4_246;
constraints(666) = r_6_246;
constraints(667) = r_2_247;
constraints(668) = r_4_247;
constraints(669) = r_7_247;
constraints(670) = r_2_248;
constraints(671) = r_4_248;
constraints(672) = r_8_248;
constraints(673) = r_2_249;
constraints(674) = r_4_249;
constraints(675) = r_9_249;
constraints(676) = r_2_256;
constraints(677) = r_5_256;
constraints(678) = r_6_256;
constraints(679) = r_2_257;
constraints(680) = r_5_257;
constraints(681) = r_7_257;
constraints(682) = r_2_258;
constraints(683) = r_5_258;
constraints(684) = r_8_258;
constraints(685) = r_2_259;
constraints(686) = r_5_259;
constraints(687) = r_9_259;
constraints(688) = r_2_267;
constraints(689) = r_6_267;
constraints(690) = r_7_267;
constraints(691) = r_2_268;
constraints(692) = r_6_268;
constraints(693) = r_8_268;
constraints(694) = r_2_269;
constraints(695) = r_6_269;
constraints(696) = r_9_269;
constraints(697) = r_2_278;
constraints(698) = r_7_278;
constraints(699) = r_8_278;
constraints(700) = r_2_279;
constraints(701) = r_7_279;
constraints(702) = r_9_279;
constraints(703) = r_2_289;
constraints(704) = r_8_289;
constraints(705) = r_9_289;
constraints(706) = r_3_345;
constraints(707) = r_4_345;
constraints(708) = r_5_345;
constraints(709) = r_3_346;
constraints(710) = r_4_346;
constraints(711) = r_6_346;
constraints(712) = r_3_347;
constraints(713) = r_4_347;
constraints(714) = r_7_347;
constraints(715) = r_3_348;
constraints(716) = r_4_348;
constraints(717) = r_8_348;
constraints(718) = r_3_349;
constraints(719) = r_4_349;
constraints(720) = r_9_349;
constraints(721) = r_3_356;
constraints(722) = r_5_356;
constraints(723) = r_6_356;
constraints(724) = r_3_357;
constraints(725) = r_5_357;
constraints(726) = r_7_357;
constraints(727) = r_3_358;
constraints(728) = r_5_358;
constraints(729) = r_8_358;
constraints(730) = r_3_359;
constraints(731) = r_5_359;
constraints(732) = r_9_359;
constraints(733) = r_3_367;
constraints(734) = r_6_367;
constraints(735) = r_7_367;
constraints(736) = r_3_368;
constraints(737) = r_6_368;
constraints(738) = r_8_368;
constraints(739) = r_3_369;
constraints(740) = r_6_369;
constraints(741) = r_9_369;
constraints(742) = r_3_378;
constraints(743) = r_7_378;
constraints(744) = r_8_378;
constraints(745) = r_3_379;
constraints(746) = r_7_379;
constraints(747) = r_9_379;
constraints(748) = r_3_389;
constraints(749) = r_8_389;
constraints(750) = r_9_389;
constraints(751) = r_4_456;
constraints(752) = r_5_456;
constraints(753) = r_6_456;
constraints(754) = r_4_457;
constraints(755) = r_5_457;
constraints(756) = r_7_457;
constraints(757) = r_4_458;
constraints(758) = r_5_458;
constraints(759) = r_8_458;
constraints(760) = r_4_459;
constraints(761) = r_5_459;
constraints(762) = r_9_459;
constraints(763) = r_4_467;
constraints(764) = r_6_467;
constraints(765) = r_7_467;
constraints(766) = r_4_468;
constraints(767) = r_6_468;
constraints(768) = r_8_468;
constraints(769) = r_4_469;
constraints(770) = r_6_469;
constraints(771) = r_9_469;
constraints(772) = r_4_478;
constraints(773) = r_7_478;
constraints(774) = r_8_478;
constraints(775) = r_4_479;
constraints(776) = r_7_479;
constraints(777) = r_9_479;
constraints(778) = r_4_489;
constraints(779) = r_8_489;
constraints(780) = r_9_489;
constraints(781) = r_5_567;
constraints(782) = r_6_567;
constraints(783) = r_7_567;
constraints(784) = r_5_568;
constraints(785) = r_6_568;
constraints(786) = r_8_568;
constraints(787) = r_5_569;
constraints(788) = r_6_569;
constraints(789) = r_9_569;
constraints(790) = r_5_578;
constraints(791) = r_7_578;
constraints(792) = r_8_578;
constraints(793) = r_5_579;
constraints(794) = r_7_579;
constraints(795) = r_9_579;
constraints(796) = r_5_589;
constraints(797) = r_8_589;
constraints(798) = r_9_589;
constraints(799) = r_6_678;
constraints(800) = r_7_678;
constraints(801) = r_8_678;
constraints(802) = r_6_679;
constraints(803) = r_7_679;
constraints(804) = r_9_679;
constraints(805) = r_6_689;
constraints(806) = r_8_689;
constraints(807) = r_9_689;
constraints(808) = r_7_789;
constraints(809) = r_8_789;
constraints(810) = r_9_789;

constraints(811) = PAC - (p_1'*D1*p_1 + p_2'*D1*p_2 + p_3'*D1*p_3 + p_4'*D1*p_4 + ...
                   p_5'*D1*p_5 + p_6'*D1*p_6 + p_7'*D1*p_7 + p_8'*D1*p_8 + p_9'*D1*p_9 + ...
                   p_123'*D1*p_123 + p_124'*D1*p_124 + p_125'*D1*p_125 + p_126'*D1*p_126 + ...
                   p_127'*D1*p_127 + p_128'*D1*p_128 + p_129'*D1*p_129 + p_134'*D1*p_134 + ...
                   p_135'*D1*p_135 + p_136'*D1*p_136 + p_137'*D1*p_137 + p_138'*D1*p_138 + ...
                   p_139'*D1*p_139 + p_145'*D1*p_145 + p_146'*D1*p_146 + p_147'*D1*p_147 + ...
                   p_148'*D1*p_148 + p_149'*D1*p_149 + p_156'*D1*p_156 + p_157'*D1*p_157 + ...
                   p_158'*D1*p_158 + p_159'*D1*p_159 + p_167'*D1*p_167 + p_168'*D1*p_168 + ...
                   p_169'*D1*p_169 + p_178'*D1*p_178 + p_179'*D1*p_179 + p_189'*D1*p_189 + ...
                   p_234'*D1*p_234 + p_235'*D1*p_235 + p_236'*D1*p_236 + p_237'*D1*p_237 + ...
                   p_238'*D1*p_238 + p_239'*D1*p_239 + p_245'*D1*p_245 + p_246'*D1*p_246 + ...
                   p_247'*D1*p_247 + p_248'*D1*p_248 + p_249'*D1*p_249 + p_256'*D1*p_256 + ...
                   p_257'*D1*p_257 + p_258'*D1*p_258 + p_259'*D1*p_259 + p_267'*D1*p_267 + ...
                   p_268'*D1*p_268 + p_269'*D1*p_269 + p_278'*D1*p_278 + p_279'*D1*p_279 + ...
                   p_289'*D1*p_289 + p_345'*D1*p_345 + p_346'*D1*p_346 + p_347'*D1*p_347 + ...
                   p_348'*D1*p_348 + p_349'*D1*p_349 + p_356'*D1*p_356 + p_357'*D1*p_357 + ...
                   p_358'*D1*p_358 + p_359'*D1*p_359 + p_367'*D1*p_367 + p_368'*D1*p_368 + ...
                   p_369'*D1*p_369 + p_378'*D1*p_378 + p_379'*D1*p_379 + p_389'*D1*p_389 + ...
                   p_456'*D1*p_456 + p_457'*D1*p_457 + p_458'*D1*p_458 + p_459'*D1*p_459 + ...
                   p_467'*D1*p_467 + p_468'*D1*p_468 + p_469'*D1*p_469 + p_478'*D1*p_478 + ...
                   p_479'*D1*p_479 + p_489'*D1*p_489 + p_567'*D1*p_567 + p_568'*D1*p_568 + ...
                   p_569'*D1*p_569 + p_578'*D1*p_578 + p_579'*D1*p_579 + p_589'*D1*p_589 + ...
                   p_678'*D1*p_678 + p_679'*D1*p_679 + p_689'*D1*p_689 + p_789'*D1*p_789 + ...
                   p_123456789'*D1*p_123456789);

constraints(812) = PAC - (p_1'*D2*p_1 + p_2'*D2*p_2 + p_3'*D2*p_3 + p_4'*D2*p_4 + ...
                   p_5'*D2*p_5 + p_6'*D2*p_6 + p_7'*D2*p_7 + p_8'*D2*p_8 + p_9'*D2*p_9 + ...
                   p_123'*D2*p_123 + p_124'*D2*p_124 + p_125'*D2*p_125 + p_126'*D2*p_126 + ...
                   p_127'*D2*p_127 + p_128'*D2*p_128 + p_129'*D2*p_129 + p_134'*D2*p_134 + ...
                   p_135'*D2*p_135 + p_136'*D2*p_136 + p_137'*D2*p_137 + p_138'*D2*p_138 + ...
                   p_139'*D2*p_139 + p_145'*D2*p_145 + p_146'*D2*p_146 + p_147'*D2*p_147 + ...
                   p_148'*D2*p_148 + p_149'*D2*p_149 + p_156'*D2*p_156 + p_157'*D2*p_157 + ...
                   p_158'*D2*p_158 + p_159'*D2*p_159 + p_167'*D2*p_167 + p_168'*D2*p_168 + ...
                   p_169'*D2*p_169 + p_178'*D2*p_178 + p_179'*D2*p_179 + p_189'*D2*p_189 + ...
                   p_234'*D2*p_234 + p_235'*D2*p_235 + p_236'*D2*p_236 + p_237'*D2*p_237 + ...
                   p_238'*D2*p_238 + p_239'*D2*p_239 + p_245'*D2*p_245 + p_246'*D2*p_246 + ...
                   p_247'*D2*p_247 + p_248'*D2*p_248 + p_249'*D2*p_249 + p_256'*D2*p_256 + ...
                   p_257'*D2*p_257 + p_258'*D2*p_258 + p_259'*D2*p_259 + p_267'*D2*p_267 + ...
                   p_268'*D2*p_268 + p_269'*D2*p_269 + p_278'*D2*p_278 + p_279'*D2*p_279 + ...
                   p_289'*D2*p_289 + p_345'*D2*p_345 + p_346'*D2*p_346 + p_347'*D2*p_347 + ...
                   p_348'*D2*p_348 + p_349'*D2*p_349 + p_356'*D2*p_356 + p_357'*D2*p_357 + ...
                   p_358'*D2*p_358 + p_359'*D2*p_359 + p_367'*D2*p_367 + p_368'*D2*p_368 + ...
                   p_369'*D2*p_369 + p_378'*D2*p_378 + p_379'*D2*p_379 + p_389'*D2*p_389 + ...
                   p_456'*D2*p_456 + p_457'*D2*p_457 + p_458'*D2*p_458 + p_459'*D2*p_459 + ...
                   p_467'*D2*p_467 + p_468'*D2*p_468 + p_469'*D2*p_469 + p_478'*D2*p_478 + ...
                   p_479'*D2*p_479 + p_489'*D2*p_489 + p_567'*D2*p_567 + p_568'*D2*p_568 + ...
                   p_569'*D2*p_569 + p_578'*D2*p_578 + p_579'*D2*p_579 + p_589'*D2*p_589 + ...
                   p_678'*D2*p_678 + p_679'*D2*p_679 + p_689'*D2*p_689 + p_789'*D2*p_789 + ...
                   p_123456789'*D2*p_123456789);

constraints(813) = PAC - (p_1'*D3*p_1 + p_2'*D3*p_2 + p_3'*D3*p_3 + p_4'*D3*p_4 + ...
                   p_5'*D3*p_5 + p_6'*D3*p_6 + p_7'*D3*p_7 + p_8'*D3*p_8 + p_9'*D3*p_9 + ...
                   p_123'*D3*p_123 + p_124'*D3*p_124 + p_125'*D3*p_125 + p_126'*D3*p_126 + ...
                   p_127'*D3*p_127 + p_128'*D3*p_128 + p_129'*D3*p_129 + p_134'*D3*p_134 + ...
                   p_135'*D3*p_135 + p_136'*D3*p_136 + p_137'*D3*p_137 + p_138'*D3*p_138 + ...
                   p_139'*D3*p_139 + p_145'*D3*p_145 + p_146'*D3*p_146 + p_147'*D3*p_147 + ...
                   p_148'*D3*p_148 + p_149'*D3*p_149 + p_156'*D3*p_156 + p_157'*D3*p_157 + ...
                   p_158'*D3*p_158 + p_159'*D3*p_159 + p_167'*D3*p_167 + p_168'*D3*p_168 + ...
                   p_169'*D3*p_169 + p_178'*D3*p_178 + p_179'*D3*p_179 + p_189'*D3*p_189 + ...
                   p_234'*D3*p_234 + p_235'*D3*p_235 + p_236'*D3*p_236 + p_237'*D3*p_237 + ...
                   p_238'*D3*p_238 + p_239'*D3*p_239 + p_245'*D3*p_245 + p_246'*D3*p_246 + ...
                   p_247'*D3*p_247 + p_248'*D3*p_248 + p_249'*D3*p_249 + p_256'*D3*p_256 + ...
                   p_257'*D3*p_257 + p_258'*D3*p_258 + p_259'*D3*p_259 + p_267'*D3*p_267 + ...
                   p_268'*D3*p_268 + p_269'*D3*p_269 + p_278'*D3*p_278 + p_279'*D3*p_279 + ...
                   p_289'*D3*p_289 + p_345'*D3*p_345 + p_346'*D3*p_346 + p_347'*D3*p_347 + ...
                   p_348'*D3*p_348 + p_349'*D3*p_349 + p_356'*D3*p_356 + p_357'*D3*p_357 + ...
                   p_358'*D3*p_358 + p_359'*D3*p_359 + p_367'*D3*p_367 + p_368'*D3*p_368 + ...
                   p_369'*D3*p_369 + p_378'*D3*p_378 + p_379'*D3*p_379 + p_389'*D3*p_389 + ...
                   p_456'*D3*p_456 + p_457'*D3*p_457 + p_458'*D3*p_458 + p_459'*D3*p_459 + ...
                   p_467'*D3*p_467 + p_468'*D3*p_468 + p_469'*D3*p_469 + p_478'*D3*p_478 + ...
                   p_479'*D3*p_479 + p_489'*D3*p_489 + p_567'*D3*p_567 + p_568'*D3*p_568 + ...
                   p_569'*D3*p_569 + p_578'*D3*p_578 + p_579'*D3*p_579 + p_589'*D3*p_589 + ...
                   p_678'*D3*p_678 + p_679'*D3*p_679 + p_689'*D3*p_689 + p_789'*D3*p_789 + ...
                   p_123456789'*D3*p_123456789);

constraints(814) = PAC - (p_1'*D4*p_1 + p_2'*D4*p_2 + p_3'*D4*p_3 + p_4'*D4*p_4 + ...
                   p_5'*D4*p_5 + p_6'*D4*p_6 + p_7'*D4*p_7 + p_8'*D4*p_8 + p_9'*D4*p_9 + ...
                   p_123'*D4*p_123 + p_124'*D4*p_124 + p_125'*D4*p_125 + p_126'*D4*p_126 + ...
                   p_127'*D4*p_127 + p_128'*D4*p_128 + p_129'*D4*p_129 + p_134'*D4*p_134 + ...
                   p_135'*D4*p_135 + p_136'*D4*p_136 + p_137'*D4*p_137 + p_138'*D4*p_138 + ...
                   p_139'*D4*p_139 + p_145'*D4*p_145 + p_146'*D4*p_146 + p_147'*D4*p_147 + ...
                   p_148'*D4*p_148 + p_149'*D4*p_149 + p_156'*D4*p_156 + p_157'*D4*p_157 + ...
                   p_158'*D4*p_158 + p_159'*D4*p_159 + p_167'*D4*p_167 + p_168'*D4*p_168 + ...
                   p_169'*D4*p_169 + p_178'*D4*p_178 + p_179'*D4*p_179 + p_189'*D4*p_189 + ...
                   p_234'*D4*p_234 + p_235'*D4*p_235 + p_236'*D4*p_236 + p_237'*D4*p_237 + ...
                   p_238'*D4*p_238 + p_239'*D4*p_239 + p_245'*D4*p_245 + p_246'*D4*p_246 + ...
                   p_247'*D4*p_247 + p_248'*D4*p_248 + p_249'*D4*p_249 + p_256'*D4*p_256 + ...
                   p_257'*D4*p_257 + p_258'*D4*p_258 + p_259'*D4*p_259 + p_267'*D4*p_267 + ...
                   p_268'*D4*p_268 + p_269'*D4*p_269 + p_278'*D4*p_278 + p_279'*D4*p_279 + ...
                   p_289'*D4*p_289 + p_345'*D4*p_345 + p_346'*D4*p_346 + p_347'*D4*p_347 + ...
                   p_348'*D4*p_348 + p_349'*D4*p_349 + p_356'*D4*p_356 + p_357'*D4*p_357 + ...
                   p_358'*D4*p_358 + p_359'*D4*p_359 + p_367'*D4*p_367 + p_368'*D4*p_368 + ...
                   p_369'*D4*p_369 + p_378'*D4*p_378 + p_379'*D4*p_379 + p_389'*D4*p_389 + ...
                   p_456'*D4*p_456 + p_457'*D4*p_457 + p_458'*D4*p_458 + p_459'*D4*p_459 + ...
                   p_467'*D4*p_467 + p_468'*D4*p_468 + p_469'*D4*p_469 + p_478'*D4*p_478 + ...
                   p_479'*D4*p_479 + p_489'*D4*p_489 + p_567'*D4*p_567 + p_568'*D4*p_568 + ...
                   p_569'*D4*p_569 + p_578'*D4*p_578 + p_579'*D4*p_579 + p_589'*D4*p_589 + ...
                   p_678'*D4*p_678 + p_679'*D4*p_679 + p_689'*D4*p_689 + p_789'*D4*p_789 + ...
                   p_123456789'*D4*p_123456789);

constraints(815) = PAC - (p_1'*D5*p_1 + p_2'*D5*p_2 + p_3'*D5*p_3 + p_4'*D5*p_4 + ...
                   p_5'*D5*p_5 + p_6'*D5*p_6 + p_7'*D5*p_7 + p_8'*D5*p_8 + p_9'*D5*p_9 + ...
                   p_123'*D5*p_123 + p_124'*D5*p_124 + p_125'*D5*p_125 + p_126'*D5*p_126 + ...
                   p_127'*D5*p_127 + p_128'*D5*p_128 + p_129'*D5*p_129 + p_134'*D5*p_134 + ...
                   p_135'*D5*p_135 + p_136'*D5*p_136 + p_137'*D5*p_137 + p_138'*D5*p_138 + ...
                   p_139'*D5*p_139 + p_145'*D5*p_145 + p_146'*D5*p_146 + p_147'*D5*p_147 + ...
                   p_148'*D5*p_148 + p_149'*D5*p_149 + p_156'*D5*p_156 + p_157'*D5*p_157 + ...
                   p_158'*D5*p_158 + p_159'*D5*p_159 + p_167'*D5*p_167 + p_168'*D5*p_168 + ...
                   p_169'*D5*p_169 + p_178'*D5*p_178 + p_179'*D5*p_179 + p_189'*D5*p_189 + ...
                   p_234'*D5*p_234 + p_235'*D5*p_235 + p_236'*D5*p_236 + p_237'*D5*p_237 + ...
                   p_238'*D5*p_238 + p_239'*D5*p_239 + p_245'*D5*p_245 + p_246'*D5*p_246 + ...
                   p_247'*D5*p_247 + p_248'*D5*p_248 + p_249'*D5*p_249 + p_256'*D5*p_256 + ...
                   p_257'*D5*p_257 + p_258'*D5*p_258 + p_259'*D5*p_259 + p_267'*D5*p_267 + ...
                   p_268'*D5*p_268 + p_269'*D5*p_269 + p_278'*D5*p_278 + p_279'*D5*p_279 + ...
                   p_289'*D5*p_289 + p_345'*D5*p_345 + p_346'*D5*p_346 + p_347'*D5*p_347 + ...
                   p_348'*D5*p_348 + p_349'*D5*p_349 + p_356'*D5*p_356 + p_357'*D5*p_357 + ...
                   p_358'*D5*p_358 + p_359'*D5*p_359 + p_367'*D5*p_367 + p_368'*D5*p_368 + ...
                   p_369'*D5*p_369 + p_378'*D5*p_378 + p_379'*D5*p_379 + p_389'*D5*p_389 + ...
                   p_456'*D5*p_456 + p_457'*D5*p_457 + p_458'*D5*p_458 + p_459'*D5*p_459 + ...
                   p_467'*D5*p_467 + p_468'*D5*p_468 + p_469'*D5*p_469 + p_478'*D5*p_478 + ...
                   p_479'*D5*p_479 + p_489'*D5*p_489 + p_567'*D5*p_567 + p_568'*D5*p_568 + ...
                   p_569'*D5*p_569 + p_578'*D5*p_578 + p_579'*D5*p_579 + p_589'*D5*p_589 + ...
                   p_678'*D5*p_678 + p_679'*D5*p_679 + p_689'*D5*p_689 + p_789'*D5*p_789 + ...
                   p_123456789'*D5*p_123456789);

constraints(816) = PAC - (p_1'*D6*p_1 + p_2'*D6*p_2 + p_3'*D6*p_3 + p_4'*D6*p_4 + ...
                   p_5'*D6*p_5 + p_6'*D6*p_6 + p_7'*D6*p_7 + p_8'*D6*p_8 + p_9'*D6*p_9 + ...
                   p_123'*D6*p_123 + p_124'*D6*p_124 + p_125'*D6*p_125 + p_126'*D6*p_126 + ...
                   p_127'*D6*p_127 + p_128'*D6*p_128 + p_129'*D6*p_129 + p_134'*D6*p_134 + ...
                   p_135'*D6*p_135 + p_136'*D6*p_136 + p_137'*D6*p_137 + p_138'*D6*p_138 + ...
                   p_139'*D6*p_139 + p_145'*D6*p_145 + p_146'*D6*p_146 + p_147'*D6*p_147 + ...
                   p_148'*D6*p_148 + p_149'*D6*p_149 + p_156'*D6*p_156 + p_157'*D6*p_157 + ...
                   p_158'*D6*p_158 + p_159'*D6*p_159 + p_167'*D6*p_167 + p_168'*D6*p_168 + ...
                   p_169'*D6*p_169 + p_178'*D6*p_178 + p_179'*D6*p_179 + p_189'*D6*p_189 + ...
                   p_234'*D6*p_234 + p_235'*D6*p_235 + p_236'*D6*p_236 + p_237'*D6*p_237 + ...
                   p_238'*D6*p_238 + p_239'*D6*p_239 + p_245'*D6*p_245 + p_246'*D6*p_246 + ...
                   p_247'*D6*p_247 + p_248'*D6*p_248 + p_249'*D6*p_249 + p_256'*D6*p_256 + ...
                   p_257'*D6*p_257 + p_258'*D6*p_258 + p_259'*D6*p_259 + p_267'*D6*p_267 + ...
                   p_268'*D6*p_268 + p_269'*D6*p_269 + p_278'*D6*p_278 + p_279'*D6*p_279 + ...
                   p_289'*D6*p_289 + p_345'*D6*p_345 + p_346'*D6*p_346 + p_347'*D6*p_347 + ...
                   p_348'*D6*p_348 + p_349'*D6*p_349 + p_356'*D6*p_356 + p_357'*D6*p_357 + ...
                   p_358'*D6*p_358 + p_359'*D6*p_359 + p_367'*D6*p_367 + p_368'*D6*p_368 + ...
                   p_369'*D6*p_369 + p_378'*D6*p_378 + p_379'*D6*p_379 + p_389'*D6*p_389 + ...
                   p_456'*D6*p_456 + p_457'*D6*p_457 + p_458'*D6*p_458 + p_459'*D6*p_459 + ...
                   p_467'*D6*p_467 + p_468'*D6*p_468 + p_469'*D6*p_469 + p_478'*D6*p_478 + ...
                   p_479'*D6*p_479 + p_489'*D6*p_489 + p_567'*D6*p_567 + p_568'*D6*p_568 + ...
                   p_569'*D6*p_569 + p_578'*D6*p_578 + p_579'*D6*p_579 + p_589'*D6*p_589 + ...
                   p_678'*D6*p_678 + p_679'*D6*p_679 + p_689'*D6*p_689 + p_789'*D6*p_789 + ...
                   p_123456789'*D6*p_123456789);

constraints(817) = PAC - (p_1'*D7*p_1 + p_2'*D7*p_2 + p_3'*D7*p_3 + p_4'*D7*p_4 + ...
                   p_5'*D7*p_5 + p_6'*D7*p_6 + p_7'*D7*p_7 + p_8'*D7*p_8 + p_9'*D7*p_9 + ...
                   p_123'*D7*p_123 + p_124'*D7*p_124 + p_125'*D7*p_125 + p_126'*D7*p_126 + ...
                   p_127'*D7*p_127 + p_128'*D7*p_128 + p_129'*D7*p_129 + p_134'*D7*p_134 + ...
                   p_135'*D7*p_135 + p_136'*D7*p_136 + p_137'*D7*p_137 + p_138'*D7*p_138 + ...
                   p_139'*D7*p_139 + p_145'*D7*p_145 + p_146'*D7*p_146 + p_147'*D7*p_147 + ...
                   p_148'*D7*p_148 + p_149'*D7*p_149 + p_156'*D7*p_156 + p_157'*D7*p_157 + ...
                   p_158'*D7*p_158 + p_159'*D7*p_159 + p_167'*D7*p_167 + p_168'*D7*p_168 + ...
                   p_169'*D7*p_169 + p_178'*D7*p_178 + p_179'*D7*p_179 + p_189'*D7*p_189 + ...
                   p_234'*D7*p_234 + p_235'*D7*p_235 + p_236'*D7*p_236 + p_237'*D7*p_237 + ...
                   p_238'*D7*p_238 + p_239'*D7*p_239 + p_245'*D7*p_245 + p_246'*D7*p_246 + ...
                   p_247'*D7*p_247 + p_248'*D7*p_248 + p_249'*D7*p_249 + p_256'*D7*p_256 + ...
                   p_257'*D7*p_257 + p_258'*D7*p_258 + p_259'*D7*p_259 + p_267'*D7*p_267 + ...
                   p_268'*D7*p_268 + p_269'*D7*p_269 + p_278'*D7*p_278 + p_279'*D7*p_279 + ...
                   p_289'*D7*p_289 + p_345'*D7*p_345 + p_346'*D7*p_346 + p_347'*D7*p_347 + ...
                   p_348'*D7*p_348 + p_349'*D7*p_349 + p_356'*D7*p_356 + p_357'*D7*p_357 + ...
                   p_358'*D7*p_358 + p_359'*D7*p_359 + p_367'*D7*p_367 + p_368'*D7*p_368 + ...
                   p_369'*D7*p_369 + p_378'*D7*p_378 + p_379'*D7*p_379 + p_389'*D7*p_389 + ...
                   p_456'*D7*p_456 + p_457'*D7*p_457 + p_458'*D7*p_458 + p_459'*D7*p_459 + ...
                   p_467'*D7*p_467 + p_468'*D7*p_468 + p_469'*D7*p_469 + p_478'*D7*p_478 + ...
                   p_479'*D7*p_479 + p_489'*D7*p_489 + p_567'*D7*p_567 + p_568'*D7*p_568 + ...
                   p_569'*D7*p_569 + p_578'*D7*p_578 + p_579'*D7*p_579 + p_589'*D7*p_589 + ...
                   p_678'*D7*p_678 + p_679'*D7*p_679 + p_689'*D7*p_689 + p_789'*D7*p_789 + ...
                   p_123456789'*D7*p_123456789);

constraints(818) = PAC - (p_1'*D8*p_1 + p_2'*D8*p_2 + p_3'*D8*p_3 + p_4'*D8*p_4 + ...
                   p_5'*D8*p_5 + p_6'*D8*p_6 + p_7'*D8*p_7 + p_8'*D8*p_8 + p_9'*D8*p_9 + ...
                   p_123'*D8*p_123 + p_124'*D8*p_124 + p_125'*D8*p_125 + p_126'*D8*p_126 + ...
                   p_127'*D8*p_127 + p_128'*D8*p_128 + p_129'*D8*p_129 + p_134'*D8*p_134 + ...
                   p_135'*D8*p_135 + p_136'*D8*p_136 + p_137'*D8*p_137 + p_138'*D8*p_138 + ...
                   p_139'*D8*p_139 + p_145'*D8*p_145 + p_146'*D8*p_146 + p_147'*D8*p_147 + ...
                   p_148'*D8*p_148 + p_149'*D8*p_149 + p_156'*D8*p_156 + p_157'*D8*p_157 + ...
                   p_158'*D8*p_158 + p_159'*D8*p_159 + p_167'*D8*p_167 + p_168'*D8*p_168 + ...
                   p_169'*D8*p_169 + p_178'*D8*p_178 + p_179'*D8*p_179 + p_189'*D8*p_189 + ...
                   p_234'*D8*p_234 + p_235'*D8*p_235 + p_236'*D8*p_236 + p_237'*D8*p_237 + ...
                   p_238'*D8*p_238 + p_239'*D8*p_239 + p_245'*D8*p_245 + p_246'*D8*p_246 + ...
                   p_247'*D8*p_247 + p_248'*D8*p_248 + p_249'*D8*p_249 + p_256'*D8*p_256 + ...
                   p_257'*D8*p_257 + p_258'*D8*p_258 + p_259'*D8*p_259 + p_267'*D8*p_267 + ...
                   p_268'*D8*p_268 + p_269'*D8*p_269 + p_278'*D8*p_278 + p_279'*D8*p_279 + ...
                   p_289'*D8*p_289 + p_345'*D8*p_345 + p_346'*D8*p_346 + p_347'*D8*p_347 + ...
                   p_348'*D8*p_348 + p_349'*D8*p_349 + p_356'*D8*p_356 + p_357'*D8*p_357 + ...
                   p_358'*D8*p_358 + p_359'*D8*p_359 + p_367'*D8*p_367 + p_368'*D8*p_368 + ...
                   p_369'*D8*p_369 + p_378'*D8*p_378 + p_379'*D8*p_379 + p_389'*D8*p_389 + ...
                   p_456'*D8*p_456 + p_457'*D8*p_457 + p_458'*D8*p_458 + p_459'*D8*p_459 + ...
                   p_467'*D8*p_467 + p_468'*D8*p_468 + p_469'*D8*p_469 + p_478'*D8*p_478 + ...
                   p_479'*D8*p_479 + p_489'*D8*p_489 + p_567'*D8*p_567 + p_568'*D8*p_568 + ...
                   p_569'*D8*p_569 + p_578'*D8*p_578 + p_579'*D8*p_579 + p_589'*D8*p_589 + ...
                   p_678'*D8*p_678 + p_679'*D8*p_679 + p_689'*D8*p_689 + p_789'*D8*p_789 + ...
                   p_123456789'*D8*p_123456789);

constraints(819) = PAC - (p_1'*D9*p_1 + p_2'*D9*p_2 + p_3'*D9*p_3 + p_4'*D9*p_4 + ...
                   p_5'*D9*p_5 + p_6'*D9*p_6 + p_7'*D9*p_7 + p_8'*D9*p_8 + p_9'*D9*p_9 + ...
                   p_123'*D9*p_123 + p_124'*D9*p_124 + p_125'*D9*p_125 + p_126'*D9*p_126 + ...
                   p_127'*D9*p_127 + p_128'*D9*p_128 + p_129'*D9*p_129 + p_134'*D9*p_134 + ...
                   p_135'*D9*p_135 + p_136'*D9*p_136 + p_137'*D9*p_137 + p_138'*D9*p_138 + ...
                   p_139'*D9*p_139 + p_145'*D9*p_145 + p_146'*D9*p_146 + p_147'*D9*p_147 + ...
                   p_148'*D9*p_148 + p_149'*D9*p_149 + p_156'*D9*p_156 + p_157'*D9*p_157 + ...
                   p_158'*D9*p_158 + p_159'*D9*p_159 + p_167'*D9*p_167 + p_168'*D9*p_168 + ...
                   p_169'*D9*p_169 + p_178'*D9*p_178 + p_179'*D9*p_179 + p_189'*D9*p_189 + ...
                   p_234'*D9*p_234 + p_235'*D9*p_235 + p_236'*D9*p_236 + p_237'*D9*p_237 + ...
                   p_238'*D9*p_238 + p_239'*D9*p_239 + p_245'*D9*p_245 + p_246'*D9*p_246 + ...
                   p_247'*D9*p_247 + p_248'*D9*p_248 + p_249'*D9*p_249 + p_256'*D9*p_256 + ...
                   p_257'*D9*p_257 + p_258'*D9*p_258 + p_259'*D9*p_259 + p_267'*D9*p_267 + ...
                   p_268'*D9*p_268 + p_269'*D9*p_269 + p_278'*D9*p_278 + p_279'*D9*p_279 + ...
                   p_289'*D9*p_289 + p_345'*D9*p_345 + p_346'*D9*p_346 + p_347'*D9*p_347 + ...
                   p_348'*D9*p_348 + p_349'*D9*p_349 + p_356'*D9*p_356 + p_357'*D9*p_357 + ...
                   p_358'*D9*p_358 + p_359'*D9*p_359 + p_367'*D9*p_367 + p_368'*D9*p_368 + ...
                   p_369'*D9*p_369 + p_378'*D9*p_378 + p_379'*D9*p_379 + p_389'*D9*p_389 + ...
                   p_456'*D9*p_456 + p_457'*D9*p_457 + p_458'*D9*p_458 + p_459'*D9*p_459 + ...
                   p_467'*D9*p_467 + p_468'*D9*p_468 + p_469'*D9*p_469 + p_478'*D9*p_478 + ...
                   p_479'*D9*p_479 + p_489'*D9*p_489 + p_567'*D9*p_567 + p_568'*D9*p_568 + ...
                   p_569'*D9*p_569 + p_578'*D9*p_578 + p_579'*D9*p_579 + p_589'*D9*p_589 + ...
                   p_678'*D9*p_678 + p_679'*D9*p_679 + p_689'*D9*p_689 + p_789'*D9*p_789 + ...
                   p_123456789'*D9*p_123456789);        
               
subject to
    constraints >= zeros(size(constraints))

cvx_end

%Result
rate = object_func;