function [xp,xm,HIp,HIm,Dp,Dm,Dpp,Dpm,Dmp,Dmm] = sbp_staggered_upwind_84th(n,h)

assert(n >= 12,'Not enough grid points');  

np=n+1;
nm=n+2;

bp=8; % Number of boundary points


%H- norm
Hm_U=[0.27058177e8 / 0.194594400e9 0 0 0 0 0 0 0; 0 0.378196601e9 / 0.464486400e9 0 0 0 0 0 0; 0 0 0.250683799e9 / 0.232243200e9 0 0 0 0 0; 0 0 0 0.48820787e8 / 0.51609600e8 0 0 0 0; 0 0 0 0 0.23953873e8 / 0.23224320e8 0 0 0; 0 0 0 0 0 0.55018021e8 / 0.55738368e8 0 0; 0 0 0 0 0 0 0.284772253e9 / 0.283852800e9 0; 0 0 0 0 0 0 0 0.6036097241e10 / 0.6038323200e10;];

%H+ norm
Hp_U=[0.273925927e9 / 0.928972800e9 0 0 0 0 0 0 0; 0 0.1417489391e10 / 0.928972800e9 0 0 0 0 0 0; 0 0 0.26528603e8 / 0.103219200e9 0 0 0 0 0; 0 0 0 0.66843235e8 / 0.37158912e8 0 0 0 0; 0 0 0 0 0.76542481e8 / 0.185794560e9 0 0 0; 0 0 0 0 0 0.132009637e9 / 0.103219200e9 0 0; 0 0 0 0 0 0 0.857580049e9 / 0.928972800e9 0; 0 0 0 0 0 0 0 0.937663193e9 / 0.928972800e9;];

% Upper part of Qpp. Notice that the B matrix is not included here.
% Qpp+Qpp^T=S/2,Qpp+Qpm^T=0 
Qpp_U=[-0.5053e4 / 0.43253760e8 0.299126491231e12 / 0.442810368000e12 -0.562745123e9 / 0.8856207360e10 -0.59616938177e11 / 0.398529331200e12 -0.5489464381e10 / 0.177124147200e12 0.6994759151e10 / 0.88562073600e11 0.712609553e9 / 0.181149696000e12 -0.14167e5 / 0.998400e6; -0.74652297589e11 / 0.110702592000e12 -0.995441e6 / 0.302776320e9 0.662285893e9 / 0.2236416000e10 0.2364992837e10 / 0.5535129600e10 0.8893339297e10 / 0.49816166400e11 -0.2166314077e10 / 0.8945664000e10 -0.78813191e8 / 0.2952069120e10 0.88810243397e11 / 0.1992646656000e13; 0.16939159e8 / 0.276756480e9 -0.2302477691e10 / 0.8200192000e10 -0.197067e6 / 0.9175040e7 0.46740527413e11 / 0.88562073600e11 -0.4827794723e10 / 0.11070259200e11 0.155404199e9 / 0.1230028800e10 0.15110006581e11 / 0.295206912000e12 -0.2491856531e10 / 0.88562073600e11; 0.7568509969e10 / 0.49816166400e11 -0.19762404121e11 / 0.44281036800e11 -0.2554553593e10 / 0.5535129600e10 -0.19550057e8 / 0.302776320e9 0.145135201e9 / 0.210862080e9 0.1596117533e10 / 0.11070259200e11 0.104983477e9 / 0.3558297600e10 -0.1087937837e10 / 0.25303449600e11; 0.5282544031e10 / 0.177124147200e12 -0.131784830977e12 / 0.797058662400e12 0.8310607171e10 / 0.22140518400e11 -0.56309573e8 / 0.105431040e9 -0.36477607e8 / 0.302776320e9 0.57443289107e11 / 0.88562073600e11 -0.220576549e9 / 0.691891200e9 0.39942520697e11 / 0.398529331200e12; -0.1743516779e10 / 0.22140518400e11 0.7784403199e10 / 0.32800768000e11 -0.459036521e9 / 0.4920115200e10 -0.5749179391e10 / 0.22140518400e11 -0.562460513e9 / 0.1383782400e10 -0.1500741e7 / 0.9175040e7 0.70986504791e11 / 0.73801728000e11 -0.930160589e9 / 0.3406233600e10; -0.712609553e9 / 0.181149696000e12 0.180761e6 / 0.6589440e7 -0.18016744831e11 / 0.295206912000e12 0.18651112477e11 / 0.797058662400e12 0.7155507361e10 / 0.44281036800e11 -0.49358401541e11 / 0.73801728000e11 -0.55032223e8 / 0.302776320e9 0.38275518139e11 / 0.40255488000e11; 0.14167e5 / 0.998400e6 -0.88810243397e11 / 0.1992646656000e13 0.650307179e9 / 0.22140518400e11 0.474654619e9 / 0.16102195200e11 -0.7267828861e10 / 0.199264665600e12 0.42227833e8 / 0.425779200e9 -0.71245142351e11 / 0.110702592000e12 -0.7998899e7 / 0.43253760e8;];

% Upper part of Qmp. Notice that the B matrix is not included here.
% Qmp+Qmp^T=S/2,Qmp+Qmm^T=0 
Qmp_U=[-0.92025012754706822244637e23 / 0.73350939131274317328275670e26 0.930302337374620084855601123690977e33 / 0.1359398284732080668181467336976000e34 -0.727851704787797291000113057457371e33 / 0.2718796569464161336362934673952000e34 0.7360201789777281105403766444579e31 / 0.90626552315472044545431155798400e32 0.13716157770579047943179985700303e32 / 0.543759313892832267272586934790400e33 -0.553550686983724599329589830113e30 / 0.21750372555713290690903477391616e32 0.39889728754596585193681401053e29 / 0.20596943708061828305779808136000e32 0.4597826214543803453588826224861e31 / 0.2718796569464161336362934673952000e34; -0.460817194517403953445488602736801e33 / 0.679699142366040334090733668488000e33 -0.2251237342833501172927156567e28 / 0.267062619272621870023659683840e30 0.190175037040902421814031988319053e33 / 0.202800676510147232549216572416000e33 -0.280264041304585424221475951435491e33 / 0.1081603608054118573595821719552000e34 -0.2948540748840639854083293613843e31 / 0.81120270604058893019686628966400e32 0.10373234521893519136055741251159e32 / 0.194688649449741343247247909519360e33 -0.391946441371315598560753809131e30 / 0.59488198442976521547770194575360e32 -0.15471439859462263133402977429037e32 / 0.6026077244872946338605292437504000e34; 0.702728186606917922841148808013121e33 / 0.2718796569464161336362934673952000e34 -0.738964213741742941634289388463587e33 / 0.811202706040588930196866289664000e33 -0.6938440904692701484464706797e28 / 0.267062619272621870023659683840e30 0.13587935221144065448123508546711e32 / 0.16900056375845602712434714368000e32 -0.398581023027193812136121833051e30 / 0.3244810824162355720787465158656e31 -0.4390487184318690815376655433561e31 / 0.486721623624353358118119773798400e33 0.39375020562052974775725982794643e32 / 0.5948819844297652154777019457536000e34 -0.337431206994896862016299739723e30 / 0.1054563517852765609255926176563200e34; -0.1626484714285090813572964936931e31 / 0.22656638078868011136357788949600e32 0.247081446636583523913555235517491e33 / 0.1081603608054118573595821719552000e34 -0.98606689532786656855690277893813e32 / 0.135200451006764821699477714944000e33 -0.54412910204947725598668603049e29 / 0.801187857817865610070979051520e30 0.44527159063104584858762191285733e32 / 0.54080180402705928679791085977600e32 -0.4170461618635408764831166558231e31 / 0.18541776138070604118785515192320e32 0.2733053052508781643445069010989e31 / 0.55081665224978260692379809792000e32 -0.62697761224771731704532896739409e32 / 0.7030423452351770728372841177088000e34; -0.16732447706243527499842078942603e32 / 0.543759313892832267272586934790400e33 0.166114324967929244010272382781e30 / 0.2897152521573531893560236748800e31 0.45101984053799299547057123965e29 / 0.811202706040588930196866289664e30 -0.17968772701532140224971787578029e32 / 0.27040090201352964339895542988800e32 -0.48635228792591105226576208151e29 / 0.400593928908932805035489525760e30 0.89480251714879473157102672858403e32 / 0.97344324724870671623623954759680e32 -0.42143863277048975845235502012773e32 / 0.148720496107441303869425486438400e33 0.4389611982769118812600028148913e31 / 0.52728175892638280462796308828160e32; 0.590475930631333482244268674627e30 / 0.21750372555713290690903477391616e32 -0.1702402353852012839500776925027e31 / 0.27812664207105906178178272788480e32 0.5538920825569170589993365115709e31 / 0.121680405906088339529529943449600e33 0.13866428236712167588433125984777e32 / 0.129792432966494228831498606346240e33 -0.16470253483258635888811378530287e32 / 0.24336081181217667905905988689920e32 -0.391812951622286986994545412081e30 / 0.2403563573453596830212937154560e31 0.345131812155051195787190432571781e33 / 0.356929190657859129286621167452160e33 -0.8126635967231661246920267342728147e34 / 0.25309524428466374622142228237516800e35; -0.88170259036818455465427409231e29 / 0.41193887416123656611559616272000e32 0.943459254925453305858100463659e30 / 0.118976396885953043095540389150720e33 -0.104249970715879362499981962934393e33 / 0.5948819844297652154777019457536000e34 0.78160374875076232344526852587e29 / 0.18360555074992753564126603264000e32 0.37536035579237342566706229296521e32 / 0.297440992214882607738850972876800e33 -0.240812609144054591374890622842541e33 / 0.356929190657859129286621167452160e33 -0.145390363516271219220036634153e30 / 0.801187857817865610070979051520e30 0.38093040928245760105862644889779897e35 / 0.38667328987934739006050626473984000e35; -0.4596580943680048141920105029111e31 / 0.2718796569464161336362934673952000e34 0.15262190991038620305994595494037e32 / 0.6026077244872946338605292437504000e34 0.1782508142749448850000523441843e31 / 0.1054563517852765609255926176563200e34 -0.33510338065544202178903991034341e32 / 0.7030423452351770728372841177088000e34 -0.8226959851063633384855147175189e31 / 0.421825407141106243702370470625280e33 0.3729744974003821244717412110773747e34 / 0.25309524428466374622142228237516800e35 -0.6554525408845613610278033864711443e34 / 0.9666832246983684751512656618496000e34 -0.49383212966905764275823531781e29 / 0.267062619272621870023659683840e30;];

% The staggered + operator, upper part. Notice that the B matrix is not included here.Qp+Qm^T=0
Qp_U=[-0.6661444046602902130086192779621746771e37 / 0.23342839720855518078613888837079040000e38 0.12367666586033683088530161144944922123649e41 / 0.14797137255429936031548004199961722880000e41 0.4395151478839780503265910147432452357e37 / 0.186518536833150454179176523528929280000e39 -0.87963490365651280757669626389462038003e38 / 0.1345194295948176002868000381814702080000e40 -0.224187853266917333808369723332131873619e39 / 0.5178998039400477611041801469986603008000e40 0.55314264508663245255306547445101748003e38 / 0.2959427451085987206309600839992344576000e40 0.232303691019191510116997121080032335473e39 / 0.7398568627714968015774002099980861440000e40 -0.2816079756494683118054172509551114287e37 / 0.182680706857159704093185237036564480000e39; -0.10353150194859080046224306922028068797e38 / 0.43350988053017390717425793554575360000e38 -0.41810120772068966379630018024826564868789e41 / 0.44391411766289808094644012599885168640000e41 0.14548989645937048285445132077600344947e38 / 0.16118885899161150361163403267932160000e38 0.1353060140543056426043612771323929599939e40 / 0.6341630252327115442092001799983595520000e40 0.52584560489710238297754746193308329991e38 / 0.317081512616355772104600089999179776000e39 -0.322880558646947920788817445027309369183e39 / 0.8878282353257961618928802519977033728000e40 -0.2523050363123545986203815545725359199053e40 / 0.22195705883144904047322006299942584320000e41 0.2171073153815036601208579021528549178587e40 / 0.44391411766289808094644012599885168640000e41; -0.891108828407068615176254357157902263e36 / 0.14450329351005796905808597851525120000e38 0.9207480733237809022786174965361787199767e40 / 0.44391411766289808094644012599885168640000e41 -0.42541158717297365539935073107511004261e38 / 0.62172845611050151393058841176309760000e38 0.3752589973871193536179209598104238527889e40 / 0.4932379085143312010516001399987240960000e40 -0.503978909830206767329311700636686423011e39 / 0.2219570588314490404732200629994258432000e40 -0.25729240076595732139825342259542295873e38 / 0.328825272342887467367733426665816064000e39 0.283553624900223723890865927423934504751e39 / 0.2466189542571656005258000699993620480000e40 -0.129077863120107719239717930927898530811e39 / 0.4035582887844528008604001145444106240000e40; 0.150348887575956035991597139943595613e36 / 0.2000814833216187263881190471749632000e37 -0.20271017027971837857434511014725441049e38 / 0.422775350155141029472800119998906368000e39 -0.28593401740189393456434271848012721991e38 / 0.87041983855470211950282377646833664000e38 -0.3107901953269564253446821657614993689969e40 / 0.2959427451085987206309600839992344576000e40 0.163990240717967340033770840655582198979e39 / 0.147971372554299360315480041999617228800e39 0.1282390462809880153203810439898877805771e40 / 0.5326969411954776971357281511986220236800e40 0.55264142141446425784998699697993554089e38 / 0.1479713725542993603154800419996172288000e40 -0.11464469241412665723941542481946601159e38 / 0.328825272342887467367733426665816064000e39; 0.523088046329055388999216632374482217e36 / 0.8670197610603478143485158710915072000e37 -0.1060078188350823496391913491020744512571e40 / 0.8878282353257961618928802519977033728000e40 0.2489758378222482566046953325470732791e37 / 0.87041983855470211950282377646833664000e38 0.2223322627542429729644405110485179851147e40 / 0.8878282353257961618928802519977033728000e40 -0.405229856535104846314030690968355556777e39 / 0.443914117662898080946440125998851686400e39 0.1490365162783652350291746570529147907183e40 / 0.1775656470651592323785760503995406745600e40 -0.768972334347761129534143069144613163467e39 / 0.4439141176628980809464401259988516864000e40 0.34866718651627637034229483614924933299e38 / 0.1268326050465423088418400359996719104000e40; -0.522920230014765612739540250860177277e36 / 0.14450329351005796905808597851525120000e38 0.1453599225900263767517831305438030179313e40 / 0.44391411766289808094644012599885168640000e41 0.39466073050115777575909309244731209107e38 / 0.435209919277351059751411888234168320000e39 -0.501532616891797921922210650058733239369e39 / 0.4932379085143312010516001399987240960000e40 -0.276688880738967796756224778348832938429e39 / 0.2219570588314490404732200629994258432000e40 -0.346414321340583244146983258114429082007e39 / 0.328825272342887467367733426665816064000e39 0.12835585968364300662881260070481636919e38 / 0.10676145205937904784666669696942080000e38 -0.823464330534329517579044232163467356639e39 / 0.44391411766289808094644012599885168640000e41; -0.52802811745049309885720368328150617e35 / 0.1688999534533145092886719229399040000e37 0.306081878756402097710283710502974342821e39 / 0.4932379085143312010516001399987240960000e40 -0.43281676756908448328630499995327205907e38 / 0.1305629757832053179254235664702504960000e40 -0.136777517026276610492513752852976312597e39 / 0.4932379085143312010516001399987240960000e40 0.7290267185112688983569547959866321207e37 / 0.246618954257165600525800069999362048000e39 0.351639553804386998900032061955116226307e39 / 0.3804978151396269265255201079990157312000e40 -0.2840607921448362125113163437018746783083e40 / 0.2466189542571656005258000699993620480000e40 0.1859197240206073478058984196606793676119e40 / 0.1644126361714437336838667133329080320000e40; 0.5412884950805861225701675068383948563e37 / 0.303456916371121735021980554882027520000e39 -0.1279848124286614098666833963684894093627e40 / 0.44391411766289808094644012599885168640000e41 0.17218770500911534891873213313315117e35 / 0.39564538116122823613764717112197120000e38 0.904771800018341402297740826650985365019e39 / 0.44391411766289808094644012599885168640000e41 0.65377500469171784914136353007658339737e38 / 0.15536994118201432833125404409959809024000e41 -0.211014619053462658139013021424372225649e39 / 0.8878282353257961618928802519977033728000e40 0.195325945159072852492812627815477159003e39 / 0.3170815126163557721046000899991797760000e40 -0.52260858238454846311625894178508086247499e41 / 0.44391411766289808094644012599885168640000e41;];

Hp=spdiags(ones(np,1),0,np,np);
Hp(1:bp,1:bp)=Hp_U;
Hp(np-bp+1:np,np-bp+1:np)=fliplr(flipud(Hp_U));
Hp=Hp*h;
HIp=inv(Hp);

Hm=spdiags(ones(nm,1),0,nm,nm);
Hm(1:bp,1:bp)=Hm_U;
Hm(nm-bp+1:nm,nm-bp+1:nm)=fliplr(flipud(Hm_U));
Hm=Hm*h;
HIm=inv(Hm);

tt=[0.1447e4 / 0.688128e6 -0.17119e5 / 0.688128e6 0.8437e4 / 0.57344e5 -0.5543e4 / 0.8192e4 -0.15159e5 / 0.81920e5 0.16139e5 / 0.16384e5 -0.2649e4 / 0.8192e4 0.15649e5 / 0.172032e6 -0.3851e4 / 0.229376e6 0.5053e4 / 0.3440640e7;];

Qpp=spdiags(repmat(tt,[np,1]),-4:5,np,np);
Qpp(1:bp,1:bp)=Qpp_U;
Qpp(np-bp+1:np,np-bp+1:np)=flipud( fliplr(Qpp_U ) )'; 

Qpm=-Qpp';

Qmp=spdiags(repmat(tt,[nm,1]),-4:5,nm,nm);
Qmp(1:bp,1:bp)=Qmp_U;
Qmp(nm-bp+1:nm,nm-bp+1:nm)=flipud( fliplr(Qmp_U ) )'; 

Qmm=-Qmp';


Bpp=spalloc(np,np,2);Bpp(1,1)=-1;Bpp(np,np)=1;
Bmp=spalloc(nm,nm,2);Bmp(1,1)=-1;Bmp(nm,nm)=1;

Dpp=HIp*(Qpp+1/2*Bpp) ;
Dpm=HIp*(Qpm+1/2*Bpp) ;


Dmp=HIm*(Qmp+1/2*Bmp) ;
Dmm=HIm*(Qmm+1/2*Bmp) ;


%%% Start with the staggered
Qp=spdiags(repmat([0.5e1 / 0.7168e4 -0.49e2 / 0.5120e4 0.245e3 / 0.3072e4 -0.1225e4 / 0.1024e4 0.1225e4 / 0.1024e4 -0.245e3 / 0.3072e4 0.49e2 / 0.5120e4 -0.5e1 / 0.7168e4;],[np,1]),-3:4,np,nm);
Qp(1:bp,1:bp)=Qp_U;
Qp(np-bp+1:np,nm-bp+1:nm)=flipud( fliplr(-Qp_U ) ); 
Qm=-Qp';

Bp=spalloc(np,nm,2);Bp(1,1)=-1;Bp(np,nm)=1;
Bm=Bp';

Dp=HIp*(Qp+1/2*Bp) ;

Dm=HIm*(Qm+1/2*Bm) ;

% grids
xp = h*[0:n]';
xm = h*[0 1/2+0:n n]';  


