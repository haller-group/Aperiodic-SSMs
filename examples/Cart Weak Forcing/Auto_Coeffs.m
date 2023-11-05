function [h030,h120,h210,h300,f030,f120,f210,f300,hc030,hc120,hc210 ...
    ,hc300,fc030,fc120,fc210,fc300,h500,h050,h320,h230,h410,h140,f500,f050,f320,f230,f410,f140,hc500,hc050,hc320,hc230,hc410,hc140,fc500,fc050,fc320,fc230,fc410,fc140...
    ] = Auto_Coeffs(gamma,l2,l3,l1c,l1)


h300 = (-0.002000098274315038 - 0.0012278056114135897.*1i).*gamma./(l2-3.*l1);

h030  = (-0.0023436771712469356 - 0.00012278861230853068.*1i).*gamma./(l2-3.*l1c);


h210 = (-0.00652664370974273 - 0.0026408368793066138.*1i).*gamma ./(l2-l1c-2.*l1);

h120 = (-0.006873403125525224 - 0.0015255907088896161.*1i).*gamma./(l2-2.*l1c-l1);

f300 = (-0.002518381494372588 - 0.00045485569140714453.*1i).*gamma./(l3-3.*l1);

f030 = (-0.0024295966989206184 + 0.0008038649960023223.*1i).*gamma./(l3-3.*l1c);

f210 = (-0.0076767746924559154 - 0.00009685959994975245.*1i).*gamma./(l3-2.*l1-l1c);

f120 = (-0.007587168013826304 + 0.0011735130978412925.*1i).*gamma./(l3-l1-2.*l1c);


% h030 = -(0.013813008195636451 + 0.000996033944032752*1i)*gamma/(l2 - 3*l1c);
% h120 = -(0.04050205667290954 + 0.009257694990541146*1i)*gamma/(l2 - l1-2*l1c);
% h210 = -(0.038622163372940145 + 0.015311760436585583*1i)*gamma/(l2 - 2*l1 - l1c);
% h300 = -(0.011947703447961826 + 0.007003117973441996*1i)*gamma/(l2 - 3*l1);
% 
% f030 = -(0.014143339101206014 - 0.004483904619644111*1i)*gamma/(l3 - 3*l1c);
% f120 = -(0.04398260498843236 - 0.0068399877601822145*1i)*gamma/(l3 - l1-2*l1c);
% f210 = -(0.04451123689935621 - 0.00006902041225595448*1i)*gamma/(l3 - 2*l1 - l1c);
% f300 = -(0.01466786866532256 + 0.0022345179315944852*1i)*gamma/(l3 - 3*l1);

hc030 = conj(h030);
hc120 = conj(h120);
hc210 = conj(h210);
hc300 = conj(h300);

fc030 = conj(f030);
fc120 = conj(f120);
fc210 = conj(f210);
fc300 = conj(f300);

 h500 = ((-0.00023622803615410554 - 0.009277726390043492.*1i).*f300.*gamma - (0.007605359729204068 - 0.005318882781466016.*1i).*fc300.*gamma + ...
  (0.000027701429951862793 - 0.00016440759757741346.*1i).*h210.*gamma + (0.009326917463392308 - 0.02394550101253702.*1i).*h300.*gamma - ...
  (0.014012952364953033 - 0.022134856413765693.*1i).*hc300.*gamma ...
)./(l2-4.*l1) ;

h050 =( (-0.003248499258299676 - 0.008693633487926733.*1i).*f030.*gamma - (0.005455370746887958 - 0.007508058382003553.*1i).*fc030.*gamma + (0.0008461542230817063 - 0.026651154132919654.*1i).*h030.*gamma + ...
  (0.000027701429951862753 + 0.0001644075975774134.*1i).*h120.*gamma - (0.006029588580057364 - 0.025494288850558457.*1i).*hc030.*gamma ...
)./(l2-5.*l1c);

h320 = ((-0.00023622803615410532 - 0.009277726390043496.*1i).*f120.*gamma - (0.0035333396150121726 + 0.018222062281182238.*1i).*f210.*gamma - (0.003248499258299676 + 0.008693633487926736.*1i).*f300.*gamma - ...
  (0.007605359729204069 - 0.005318882781466017.*1i).*fc120.*gamma - (0.013242929070983806 - 0.013005878372304993.*1i).*fc210.*gamma - ...
  (0.005455370746887958 - 0.0075080583820035535.*1i).*fc300.*gamma + (0.00008310428985558839 - 0.0004932227927322404.*1i).*h030.*gamma + (0.009436229796470096 - 0.0252611689177229.*1i).*h120.*gamma + ...
  (0.010398115069578722 - 0.050809120042499656.*1i).*h210.*gamma + (0.0010101227226983847 - 0.024677652275140828.*1i).*h300.*gamma - (0.014012952364953037 - 0.0221348564137657.*1i).*hc120.*gamma - ...
  (0.020322136546874044 - 0.04829357696352668.*1i).*hc210.*gamma - (0.006029588580057366 - 0.02549428885055845.*1i).*hc300.*gamma ...
)./(l2 - 3.*l1 - 2.*l1c) ;

h230 =( (-0.00023622803615410554 - 0.00927772639004349.*1i).*f030.*gamma - (0.003533339615012173 + 0.01822206228118223.*1i).*f120.*gamma - (0.0032484992582996774 + 0.008693633487926731.*1i).*f210.*gamma - ...
  (0.007605359729204069 - 0.005318882781466015.*1i).*fc030.*gamma - (0.013242929070983806 - 0.013005878372304987.*1i).*fc120.*gamma - (0.005455370746887958 - 0.007508058382003553.*1i).*fc210.*gamma + ...
  (0.009490885963008986 - 0.025919002870315844.*1i).*h030.*gamma + (0.010398115069578724 - 0.05179584781084494.*1i).*h120.*gamma + (0.0009554665561594903 - 0.02533548622773376.*1i).*h210.*gamma + ...
  (0.00008310428985558828 + 0.0004932227927322403.*1i).*h300.*gamma - (0.014012952364953033 - 0.022134856413765693.*1i).*hc030.*gamma - ...
  (0.020322136546874044 - 0.04829357696352668.*1i).*hc120.*gamma - (0.006029588580057366 - 0.02549428885055845.*1i).*hc210.*gamma ...
) ./(l2 - 2.*l1 - 3.*l1c) ;


h410 =( (-0.00023622803615410554 - 0.009277726390043492.*1i).*f210.*gamma - (0.0035333396150121696 + 0.018222062281182234.*1i).*f300.*gamma - ...
  (0.007605359729204068 - 0.005318882781466016.*1i).*fc210.*gamma - (0.013242929070983806 - 0.013005878372304994.*1i).*fc300.*gamma + ...
  (0.000055402859903725586 - 0.0003288151951548269.*1i).*h120.*gamma + (0.009381573629931198 - 0.024603334965129956.*1i).*h210.*gamma + ...
  (0.010398115069578722 - 0.049822392274154366.*1i).*h300.*gamma - (0.014012952364953033 - 0.022134856413765693.*1i).*hc210.*gamma - (0.020322136546874047 - 0.04829357696352668.*1i).*hc300.*gamma ...
)./(l2-4.*l1-l1c) ;

h140 =( (-0.003533339615012171 - 0.01822206228118224.*1i).*f030.*gamma - (0.003248499258299678 + 0.008693633487926733.*1i).*f120.*gamma - (0.013242929070983806 - 0.013005878372304993.*1i).*fc030.*gamma - ...
  (0.005455370746887958 - 0.0075080583820035535.*1i).*fc120.*gamma + (0.010398115069578724 - 0.05278257557919024.*1i).*h030.*gamma + (0.000900810389620599 - 0.025993320180326712.*1i).*h120.*gamma + ...
  (0.00005540285990372551 + 0.00032881519515482686.*1i).*h210.*gamma - (0.020322136546874044 - 0.048293576963526666.*1i).*hc030.*gamma - ...
  (0.006029588580057364 - 0.02549428885055845.*1i).*hc120.*gamma ...
) ./(l2-l1-4.*l1c) ;


f500 = ((0.000027701429951862793 - 0.00016440759757741346.*1i).*f210.*gamma - (0.004078671263120068 + 0.00885868321014472.*1i).*f300.*gamma - ...
  (0.005618863087675488 - 0.008416839639608065.*1i).*fc300.*gamma - (0.00003452421379762649 + 0.02856672324055655.*1i).*h300.*gamma - (0.005465562695104693 - 0.028039017333882654.*1i).*hc300.*gamma ...
)./(l3 - 5.*l1);

f050 =( (-0.006907414912072124 - 0.008017742041373812.*1i).*f030.*gamma + (0.00002770142995186278 + 0.00016440759757741343.*1i).*f120.*gamma - ...
  (0.0025672956478306514 - 0.009788963426918578.*1i).*fc030.*gamma - (0.009347369737784578 + 0.026994176179484367.*1i).*h030.*gamma + (0.003975816651248213 + 0.0282887212609469.*1i).*hc030.*gamma ...
)./(l3-5.*l1c);


f320=( (0.00008310428985558839 - 0.0004932227927322404.*1i).*f030.*gamma - (0.003969358930042283 + 0.0101743511153306.*1i).*f120.*gamma - (0.011056215320036504 + 0.016618489314752922.*1i).*f210.*gamma - ...
  (0.0067434464124554476 + 0.006044240183594991.*1i).*f300.*gamma - (0.00561886308767549 - 0.008416839639608069.*1i).*fc120.*gamma - (0.008300356530330903 - 0.018459775977434624.*1i).*fc210.*gamma - ...
  (0.0025672956478306535 - 0.009788963426918582.*1i).*fc300.*gamma - (0.00003452421379762649 + 0.028566723240556556.*1i).*h120.*gamma - (0.009512772381278846 + 0.05633597994281762.*1i).*h210.*gamma - ...
  (0.009347369737784583 + 0.02699417617948437.*1i).*h300.*gamma - (0.005465562695104693 - 0.028039017333882654.*1i).*hc120.*gamma - (0.001510528161398305 - 0.05711351660657216.*1i).*hc210.*gamma + ...
  (0.0039758166512482195 + 0.028288721260946913.*1i).*hc300.*gamma ...
)./(l3-3.*l1-2.*l1c);


f230=( (-0.00391470276350339 - 0.010832185067923543.*1i).*f030.*gamma - (0.011056215320036505 + 0.017605217083098205.*1i).*f120.*gamma - (0.006798102578994342 + 0.0067020741361879296.*1i).*f210.*gamma + ...
  (0.0000831042898555883 + 0.0004932227927322403.*1i).*f300.*gamma - (0.00561886308767549 - 0.008416839639608065.*1i).*fc030.*gamma - (0.0083003565303309 - 0.018459775977434617.*1i).*fc120.*gamma - ...
  (0.0025672956478306527 - 0.009788963426918582.*1i).*fc210.*gamma - (0.00003452421379762302 + 0.028566723240556553.*1i).*h030.*gamma - (0.009512772381278855 + 0.05633597994281763.*1i).*h120.*gamma - ...
  (0.009347369737784585 + 0.026994176179484373.*1i).*h210.*gamma - (0.005465562695104693 - 0.02803901733388266.*1i).*hc030.*gamma - (0.0015105281613982946 - 0.05711351660657216.*1i).*hc120.*gamma + ...
  (0.003975816651248217 + 0.028288721260946913.*1i).*hc210.*gamma ...
)./(l3-2.*l1-3.*l1c);


f410 =( (0.000055402859903725586 - 0.0003288151951548269.*1i).*f120.*gamma - (0.004024015096581174 + 0.00951651716273766.*1i).*f210.*gamma - (0.011056215320036502 + 0.015631761546407625.*1i).*f300.*gamma - ...
  (0.005618863087675488 - 0.008416839639608065.*1i).*fc210.*gamma - (0.0083003565303309 - 0.018459775977434617.*1i).*fc300.*gamma - (0.00003452421379762649 + 0.02856672324055655.*1i).*h210.*gamma - ...
  (0.009512772381278846 + 0.056335979942817596.*1i).*h300.*gamma - (0.005465562695104693 - 0.028039017333882654.*1i).*hc210.*gamma - (0.001510528161398298 - 0.05711351660657212.*1i).*hc300.*gamma ...
)./(l3 - 4.*l1-l1c);


f140=( (-0.011056215320036502 - 0.018591944851443498.*1i).*f030.*gamma - (0.006852758745533234 + 0.0073599080887808725.*1i).*f120.*gamma + ...
  (0.00005540285990372548 + 0.00032881519515482686.*1i).*f210.*gamma - (0.008300356530330904 - 0.018459775977434617.*1i).*fc030.*gamma - ...
  (0.002567295647830653 - 0.009788963426918578.*1i).*fc120.*gamma - (0.009512772381278848 + 0.05633597994281762.*1i).*h030.*gamma - (0.009347369737784583 + 0.026994176179484373.*1i).*h120.*gamma - ...
  (0.0015105281613983015 - 0.05711351660657217.*1i).*hc030.*gamma + (0.003975816651248215 + 0.028288721260946913.*1i).*hc120.*gamma ...
)./(l3-4.*l1c-l1) ;




%  h500 = ((-0.0012617843553779954 - 0.03030346023774042.*sqrt(-1)).*f300.*gamma - (0.024657040863375276 - 0.017661317561315525.*sqrt(-1)).*fc300.*gamma + (0.0001701727879495008 - 0.0010946704818697492.*sqrt(-1)).*h210.*gamma + ...
%   (0.028608293871942997 - 0.07658099680681464.*sqrt(-1)).*h300.*gamma - (0.04462316572619118 - 0.07242588638683141.*sqrt(-1)).*hc300.*gamma )./(l2-4.*l1) ;
% 
% h050 =( (-0.010316551347688505 - 0.02852123017035375.*sqrt(-1)).*f030.*gamma - (0.018204347823182657 - 0.02425888543343193.*sqrt(-1)).*fc030.*gamma + (0.0032540536354470134 - 0.08813393552782134.*sqrt(-1)).*h030.*gamma + ...
%   (0.0001701727879495008 + 0.0010946704818697492.*sqrt(-1)).*h120.*gamma - (0.020776845226240983 - 0.08249277932516373.*sqrt(-1)).*hc030.*gamma)./(l2-5.*l1c);
% 
% h320 = ((-0.0012617843553779972 - 0.030303460237740422.*sqrt(-1)).*f120.*gamma - (0.011714700237991642 + 0.05951750168558122.*sqrt(-1)).*f210.*gamma - (0.010316551347688503 + 0.02852123017035376.*sqrt(-1)).*f300.*gamma - ... 
%   (0.024657040863375286 - 0.01766131756131553.*sqrt(-1)).*fc120.*gamma - (0.043366191231965265 - 0.042413920669890734.*sqrt(-1)).*fc210.*gamma - (0.018204347823182657 - 0.024258885433431934.*sqrt(-1)).*fc300.*gamma + ...
%   (0.0005105183638485024 - 0.003284011445609248.*sqrt(-1)).*h030.*gamma + (0.02928136751686007 - 0.08534196695071036.*sqrt(-1)).*h120.*gamma + (0.0327483588007034 - 0.16336935198692218.*sqrt(-1)).*h210.*gamma + ...
%   (0.004263664102822608 - 0.07499248031197782.*sqrt(-1)).*h300.*gamma - (0.04462316572619121 - 0.07242588638683144.*sqrt(-1)).*hc120.*gamma - (0.06617026345731211 - 0.1567432294785688.*sqrt(-1)).*hc210.*gamma - ...
%   (0.020776845226240986 - 0.08249277932516373.*sqrt(-1)).*hc300.*gamma)./(l2 - 3.*l1 - 2.*l1c) ;
% 
% h230 = ((-0.0012617843553779998 - 0.030303460237740415.*sqrt(-1)).*f030.*gamma - (0.011714700237991638 + 0.05951750168558122.*sqrt(-1)).*f120.*gamma - (0.010316551347688503 + 0.02852123017035376.*sqrt(-1)).*f210.*gamma - ...
%   (0.024657040863375286 - 0.01766131756131553.*sqrt(-1)).*fc030.*gamma - (0.04336619123196525 - 0.04241392066989072.*sqrt(-1)).*fc120.*gamma - (0.018204347823182657 - 0.024258885433431937.*sqrt(-1)).*fc210.*gamma + ...
%   (0.029617904339318598 - 0.08972245202265822.*sqrt(-1)).*h030.*gamma + (0.0327483588007034 - 0.16994039243693151.*sqrt(-1)).*h120.*gamma + (0.003927127280364075 - 0.07937296538392563.*sqrt(-1)).*h210.*gamma + ...
%   (0.0005105183638485026 + 0.0032840114456092473.*sqrt(-1)).*h300.*gamma - (0.04462316572619119 - 0.07242588638683142.*sqrt(-1)).*hc030.*gamma - (0.06617026345731211 - 0.15674322947856875.*sqrt(-1)).*hc120.*gamma - ...
%   (0.02077684522624098 - 0.08249277932516376.*sqrt(-1)).*hc210.*gamma) ./(l2 - 2.*l1 - 3.*l1c) ;
% 
% 
% h410 = ((-0.0012617843553779954 - 0.03030346023774042.*sqrt(-1)).*f210.*gamma - (0.01171470023799163 + 0.05951750168558122.*sqrt(-1)).*f300.*gamma - (0.024657040863375276 - 0.017661317561315525.*sqrt(-1)).*fc210.*gamma - ...
%   (0.04336619123196524 - 0.042413920669890734.*sqrt(-1)).*fc300.*gamma + (0.0003403455758990016 - 0.0021893409637394985.*sqrt(-1)).*h120.*gamma + (0.028944830694401527 - 0.0809614818787625.*sqrt(-1)).*h210.*gamma + ...
%   (0.032748358800703384 - 0.15679831153691273.*sqrt(-1)).*h300.*gamma - (0.04462316572619118 - 0.07242588638683141.*sqrt(-1)).*hc210.*gamma - (0.06617026345731208 - 0.1567432294785687.*sqrt(-1)).*hc300.*gamma)./(l2-4.*l1-l1c) ;
% 
% h140 = ((-0.01171470023799164 - 0.05951750168558123.*sqrt(-1)).*f030.*gamma - (0.010316551347688505 + 0.028521230170353764.*sqrt(-1)).*f120.*gamma - (0.04336619123196525 - 0.04241392066989075.*sqrt(-1)).*fc030.*gamma - ...
%   (0.018204347823182657 - 0.024258885433431944.*sqrt(-1)).*fc120.*gamma + (0.0327483588007034 - 0.17651143288694096.*sqrt(-1)).*h030.*gamma + (0.0035905904579055393 - 0.08375345045587348.*sqrt(-1)).*h120.*gamma + ...
%   (0.00034034557589900186 + 0.0021893409637394985.*sqrt(-1)).*h210.*gamma - (0.06617026345731211 - 0.15674322947856875.*sqrt(-1)).*hc030.*gamma - (0.02077684522624099 - 0.08249277932516375.*sqrt(-1)).*hc120.*gamma) ./(l2-l1-4.*l1c) ;
% 
% 
% f500 =((0.0001701727879495008 - 0.0010946704818697492.*sqrt(-1)).*f210.*gamma - (0.014262508377861841 + 0.02649074195329167.*sqrt(-1)).*f300.*gamma - (0.01754122781913827 - 0.02735258501656721.*sqrt(-1)).*fc300.*gamma - ...
%   (0.002147814562160637 + 0.09111402560420669.*sqrt(-1)).*h300.*gamma - (0.015706213622955933 - 0.08977579641990115.*sqrt(-1)).*hc300.*gamma )./(l3 - 5.*l1);
% 
% f050 = ((-0.022568553119841402 - 0.027469841557874904.*sqrt(-1)).*f030.*gamma + (0.00017017278794950096 + 0.0010946704818697492.*sqrt(-1)).*f120.*gamma - ...
%   (0.008503471694552253 - 0.031361593544196646.*sqrt(-1)).*fc030.*gamma - (0.02944916082924804 + 0.08625036635081934.*sqrt(-1)).*h030.*gamma + (0.012019164418627895 + 0.0903433365327817.*sqrt(-1)).*hc030.*gamma)./(l3-5.*l1c);
% 
% 
% f320=((0.0005105183638485024 - 0.003284011445609248.*sqrt(-1)).*f030.*gamma - (0.013589434732944789 + 0.03525171209718738.*sqrt(-1)).*f120.*gamma - (0.03675409090001605 + 0.051310587256402335.*sqrt(-1)).*f210.*gamma - ...
%   (0.02155894265246581 + 0.01432838634203135.*sqrt(-1)).*f300.*gamma - (0.01754122781913827 - 0.02735258501656721.*sqrt(-1)).*fc120.*gamma - (0.026351442505734617 - 0.059405688278419556.*sqrt(-1)).*fc210.*gamma - ...
%   (0.008503471694552253 - 0.031361593544196646.*sqrt(-1)).*fc300.*gamma - (0.0021478145621606493 + 0.09111402560420669.*sqrt(-1)).*h120.*gamma - (0.0319691106416545 + 0.1794533115926578.*sqrt(-1)).*h210.*gamma - ...
%   (0.029449160829248047 + 0.08625036635081935.*sqrt(-1)).*h300.*gamma - (0.01570621362295594 - 0.08977579641990116.*sqrt(-1)).*hc120.*gamma - (0.003730473644842497 - 0.18224049671567236.*sqrt(-1)).*hc210.*gamma + ...
%   (0.0120191644186279 + 0.09034333653278172.*sqrt(-1)).*hc300.*gamma)./(l3-3.*l1-2.*l1c);
% 
% 
% f230 =((-0.013252897910486257 - 0.03963219716913523.*sqrt(-1)).*f030.*gamma - (0.03675409090001604 + 0.05788162770641174.*sqrt(-1)).*f120.*gamma - (0.021895479474924342 + 0.0187088714139792.*sqrt(-1)).*f210.*gamma + ...
%   (0.000510518363848502 + 0.003284011445609247.*sqrt(-1)).*f300.*gamma - (0.017541227819138278 - 0.02735258501656721.*sqrt(-1)).*fc030.*gamma - (0.02635144250573461 - 0.05940568827841955.*sqrt(-1)).*fc120.*gamma - ...
%   (0.008503471694552253 - 0.031361593544196646.*sqrt(-1)).*fc210.*gamma - (0.002147814562160646 + 0.09111402560420669.*sqrt(-1)).*h030.*gamma - (0.03196911064165448 + 0.17945331159265773.*sqrt(-1)).*h120.*gamma - ...
%   (0.029449160829248037 + 0.08625036635081934.*sqrt(-1)).*h210.*gamma - (0.015706213622955927 - 0.08977579641990115.*sqrt(-1)).*hc030.*gamma - (0.0037304736448425177 - 0.18224049671567236.*sqrt(-1)).*hc120.*gamma + ...
%   (0.012019164418627892 + 0.09034333653278172.*sqrt(-1)).*hc210.*gamma)./(l3-2.*l1-3.*l1c);
% 
% 
% f410 = ((0.0003403455758990016 - 0.0021893409637394985.*sqrt(-1)).*f120.*gamma - (0.01392597155540332 + 0.030871227025239516.*sqrt(-1)).*f210.*gamma - (0.03675409090001603 + 0.04473954680639291.*sqrt(-1)).*f300.*gamma - ...
%   (0.01754122781913827 - 0.02735258501656721.*sqrt(-1)).*fc210.*gamma - (0.026351442505734617 - 0.05940568827841955.*sqrt(-1)).*fc300.*gamma - (0.002147814562160637 + 0.09111402560420669.*sqrt(-1)).*h210.*gamma - ...
%   (0.031969110641654495 + 0.17945331159265768.*sqrt(-1)).*h300.*gamma - (0.015706213622955933 - 0.08977579641990115.*sqrt(-1)).*hc210.*gamma - (0.003730473644842504 - 0.18224049671567236.*sqrt(-1)).*hc300.*gamma)./(l3 - 4.*l1-l1c);
% 
% 
% f140 =((-0.036754090900016044 - 0.06445266815642114.*sqrt(-1)).*f030.*gamma - (0.02223201629738288 + 0.023089356485927054.*sqrt(-1)).*f120.*gamma + (0.0003403455758990014 + 0.0021893409637394985.*sqrt(-1)).*f210.*gamma - ...
%   (0.026351442505734617 - 0.05940568827841955.*sqrt(-1)).*fc030.*gamma - (0.008503471694552263 - 0.03136159354419665.*sqrt(-1)).*fc120.*gamma - (0.031969110641654495 + 0.17945331159265773.*sqrt(-1)).*h030.*gamma - ...
%   (0.029449160829248054 + 0.08625036635081934.*sqrt(-1)).*h120.*gamma - (0.003730473644842469 - 0.18224049671567238.*sqrt(-1)).*hc030.*gamma + (0.012019164418627905 + 0.09034333653278173.*sqrt(-1)).*hc120.*gamma)./(l3-4.*l1c-l1) ;


hc500 = conj(h500);

hc050 = conj(h050);

hc320 = conj(h320);

hc230 = conj(h230);

hc410 = conj(h410);

hc140 = conj(h140);


fc500 = conj(f500);

fc050 = conj(f050);

fc320 = conj(f320);

fc230 = conj(f230);

fc410 = conj(f410);

fc140 = conj(f140);

end