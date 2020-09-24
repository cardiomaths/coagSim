function[]=coagSim(donors, inhib)

format long;
%tic

fid=fopen('coagData.csv');
df = textscan(fid,'%s %f %f %f %f %f %f %f %f %f','HeaderLines',1,'Delimiter',',');
fclose(fid);

%initial conditions (set to zero)
u1=0; u2=0; u7=0;
u10=0; u11=0; u12=0; u13=0; u14=0; u15=0; u16=0; u17=0; u18=0; u19=0; u20=0; u21=0;
u22=0; u23=0; u24=0; u25=0; u26=0; u27=0; u28=0; u29=0; u30=0; u31=0;
u33=0; u37=0; u38=0;
u41=0; u42=0; u43=0;

%baseline procoagulants in a 'normal' population
u3=1.6e-7; %x 160nM
u4=2e-8; %v 20nM
u5=1.4e-6; %Thrombin 1400nM
u6=3.4e-6;%ATIII 3400nM
u8=9e-8; %ix 90nM
u9=7e-10; %viii 0.7nM
u32=1.39e-9; %TFPI 2.4nM
u36=1e-8;%VII 10nM
u35=1e-10;%VIIa 0.1nM
u34=1e-12;

if donors=="Case"
    p5=(df{4}(1)/100)*u5; %FII
    p4=(df{5}(1)/100)*u4; %FV
    p36=(df{6}(1)/100)*u36; %FVII
    p35=(df{6}(1)/100)*u35; %FVII
    p3=(df{7}(1)/100)*u3; %FX
    p9=(df{8}(1)/100)*u9; %FVIII
    p8=(df{9}(1)/100)*u8; %FIX
    p6=(df{10}(1)/100)*u6; %AT
    p32=df{3}(1)*u32; %TFPI
    p34=df{2}(1)*u34;%TF
end

if donors=="Control"
    p5=(df{4}(2)/100)*u5; %FII
    p4=(df{5}(2)/100)*u4; %FV
    p36=(df{6}(2)/100)*u36; %FVII
    p35=(df{6}(2)/100)*u35; %FVII
    p3=(df{7}(2)/100)*u3; %FX
    p9=(df{8}(2)/100)*u9; %FVIII
    p8=(df{9}(2)/100)*u8; %FIX
    p6=(df{10}(2)/100)*u6; %AT
    p32=df{3}(2)*u32; %TFPI
    p34=df{2}(2)*u34;%TF
end

if donors=="MaleCase"
    p5=(df{4}(3)/100)*u5; %FII
    p4=(df{5}(3)/100)*u4; %FV
    p36=(df{6}(3)/100)*u36; %FVII
    p35=(df{6}(3)/100)*u35; %FVII
    p3=(df{7}(3)/100)*u3; %FX
    p9=(df{8}(3)/100)*u9; %FVIII
    p8=(df{9}(3)/100)*u8; %FIX
    p6=(df{10}(3)/100)*u6; %AT
    p32=df{3}(3)*u32; %TFPI
    p34=df{2}(3)*u34;%TF
end

if donors=="MaleControl"
    p5=(df{4}(4)/100)*u5; %FII
    p4=(df{5}(4)/100)*u4; %FV
    p36=(df{6}(4)/100)*u36; %FVII
    p35=(df{6}(4)/100)*u35; %FVII
    p3=(df{7}(4)/100)*u3; %FX
    p9=(df{8}(4)/100)*u9; %FVIII
    p8=(df{9}(4)/100)*u8; %FIX
    p6=(df{10}(4)/100)*u6; %AT
    p32=df{3}(4)*u32; %TFPI
    p34=df{2}(4)*u34;%TF
end

if donors=="FemaleCase"
    p5=(df{4}(5)/100)*u5; %FII
    p4=(df{5}(5)/100)*u4; %FV
    p36=(df{6}(5)/100)*u36; %FVII
    p35=(df{6}(5)/100)*u35; %FVII
    p3=(df{7}(5)/100)*u3; %FX
    p9=(df{8}(5)/100)*u9; %FVIII
    p8=(df{9}(5)/100)*u8; %FIX
    p6=(df{10}(5)/100)*u6; %AT
    p32=df{3}(5)*u32; %TFPI
    p34=df{2}(6)*u34;%TF
end

if donors=="FemaleControl"
    p5=(df{4}(6)/100)*u5; %FII
    p4=(df{5}(6)/100)*u4; %FV
    p36=(df{6}(6)/100)*u36; %FVII
    p35=(df{6}(6)/100)*u35; %FVII
    p3=(df{7}(6)/100)*u3; %FX
    p9=(df{8}(6)/100)*u9; %FVIII
    p8=(df{9}(6)/100)*u8; %FIX
    p6=(df{10}(6)/100)*u6; %AT
    p32=df{3}(6)*u32; %TFPI
    p34=df{2}(6)*u34;%TF
end

 if inhib=="none"
    XaI=0;
    ThI=0;
 end
 
if strcmp(inhib, 'Warfarin25')
    p3=u3*0.25;
    p5=u5*0.25;
    p8=u8*0.25;
    p36=u36*0.25;
    p35=u35*0.25;
    XaI=0;
    ThI=0;
end;

if inhib=="Xa"
    XaI=1;
    u11=4e-9;
    ThI=0;
end;

if inhib=="Th"
    XaI=0;
    u41=3e-7;
    ThI=1;
end;

% Now put all the initial conditions together
u0=[u1 u2 p3 p4 p5 p6 u7 p8 p9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 u31 p32 u33 p34 p35 p36 u37 u38 0 0 u41 u42 u43];

%Parameters - taken from (Hockin, 2002)
a1=3.1e-3; a2=3.2e6;
a3=3.1e-3; a4=2.3e7;
a5=4.4e5;
a6=1.3e7;
a7=2.3e4;
a8=1.05; a9=2.5e7; a10=6;
a11=19; a12=2.2e7;
a13=2.4; a14=1e7; a15=1.8;
a16=7.5e3;
a17=2e7;
a18=5e-3; a19=1e7;
a20=1e-3; a21=1e8; a22=8.2;
a23=2.2e4; a24=6e-3;
a25=1e-3;
a26=2e7;
a27=0.2; a28=4e8;
a29=103; a30=1e8;
a31=63.5;
a32=1.5e7;
a33=3.6e-4; a34=9e5;
a35=1.1e-4; a36=3.2e8;
a37=5e7;
a38=1.5e3;
a39=7.1e3;
a40=4.9e2;
a41=7.1e3;
a42=2.3e2;

%introduced in later papers
a43=5.7e3;
a44=3e6;

a32=2.3e8;
a38=4.2e3;

%Time, the plots in the manscript go out to 20mins (20*60=1200)
etim=10000;
tspan=[0:1:etim];

options=odeset('RelTol',1e-13,'AbsTol',1e-60);
%tspan = span;

[t,u] = ode15s(@eq,tspan,u0,options);
d(:,1)=t;
d(:,2:44)=u(:,1:43);

fig1 = figure;
subplot(3,3,1);
plot(d(:,1)/60,d(:,35)*1e9,'-','LineWidth',2');
xlim([0,20]);
title('TF');
subplot(3,3,2);
plot(d(:,1)/60,d(:,6)*1e9,'-','LineWidth',2');
xlim([0,20]);
title('II');
subplot(3,3,3);
plot(d(:,1)/60,d(:,5)*1e9,'-','LineWidth',2');
xlim([0,20]);
title('V (nM)');
subplot(3,3,4);
plot(d(:,1)/60,d(:,37)*1e9,'-','LineWidth',2');
xlim([0,60]);
title('VII (nM)');
subplot(3,3,5);
plot(d(:,1)/60,d(:,10)*1e9,'-','LineWidth',2');
xlim([0,60]);
title('VIII (nM)');
subplot(3,3,6);
plot(d(:,1)/60,d(:,9)*1e9,'-','LineWidth',2');
xlim([0,100]);
title('IX (nM)');
subplot(3,3,7);
plot(d(:,1)/60,d(:,4)*1e9,'-','LineWidth',2');
xlim([0,100]);
title('X (nM)');
subplot(3,3,8);
plot(d(:,1)/60,d(:,33)*1e9,'-','LineWidth',2');
xlim([0,100]);
title('TFPI (nM)');
subplot(3,3,9);
plot(d(:,1)/60,d(:,7)*1e9,'-','LineWidth',2');
xlim([0,100]);
title('AT (nM)');

saveas(fig1,sprintf('Fig_ProcoagM%s%s.png',donors,inhib));

fig2 = figure;
subplot(4,4,1);
plot(d(:,1)/60,(d(:,17)+d(:,18))*1e9);
title('Thrombin (nM)');
xlim([0,20]);
ylim([0,750]);
subplot(4,4,2);
plot(d(:,1)/60,d(:,17)*1e9);
title('Thrombin (nM)');
xlim([0,20]);
ylim([0,750]);
subplot(4,4,3);
plot(d(:,1)/60,d(:,15)*1e9);
title('Va (nM)');
xlim([0,100]);
ylim([0,25]);
subplot(4,4,4);
plot(d(:,1)/60,d(:,36)*1e9);
title('VIIa (nM)');
xlim([0,100]);
ylim([0,25]);
subplot(4,4,5);
plot(d(:,1)/60,d(:,20)*1e9);
title('VIIIa (nM)');
xlim([0,100]);
ylim([0,2]);
subplot(4,4,6);
plot(d(:,1)/60,(d(:,22)+d(:,24))*1e9);
title('VIIIa other (nM)');
xlim([0,100]);
ylim([0,3]);
subplot(4,4,7);
plot(d(:,1)/60,d(:,19)*1e9);
title('IXa (nM)');
xlim([0,100]);
ylim([0,0.1]);
subplot(4,4,8);
plot(d(:,1)/60,d(:,14)*1e9);
title('Xa (nM)');
xlim([0,100]);
ylim([0,20]);
subplot(4,4,9);
plot(d(:,1)/60,d(:,13)*1e9);
title('TF:VIIa (nM)');
xlim([0,100]);
ylim([0,0.0003]);
subplot(4,4,10);
plot(d(:,1)/60,d(:,29)*1e9);
title('TF:VIIa:X (nM)');
xlim([0,100]);
ylim([0,0.0003]);
subplot(4,4,11);
plot(d(:,1)/60,d(:,39)*1e9);
title('TF:VIIa:IX (nM)');
xlim([0,100]);
ylim([0,0.0003]);
subplot(4,4,12);
plot(d(:,1)/60,(d(:,16)+d(:,11))*1e9);
title('Xa:Va (nM)');
xlim([0,100]);
ylim([0,25]);
subplot(4,4,13);
plot(d(:,1)/60,(d(:,21)+d(:,28))*1e9);
title('IX:VIIIa (nM)');
subplot(4,4,14);
plot(d(:,1)/60,(d(:,23)+d(:,26)+d(:,27)+d(:,25)+d(:,28))*1e9);
title('AT usage (nM)');
subplot(4,4,15);
plot(d(:,1)/60,d(:,32)*1e9);
title('TFPI usage by Xa (nM)');
xlim([0,100]);
ylim([0,1.5]);
subplot(4,4,16);
plot(d(:,1)/60,d(:,34)*1e9);
title('TFPI (nM)');
xlim([0,100]);
ylim([0,0.1]);

saveas(fig2,sprintf('Fig_ActiveM%s%s.png',donors,inhib));

    function du=eq(t,u)
        
        du=zeros(43,1);
        
        if XaI==0
            du(1)=0;
            du(2)=0;
            du(11)=0;
            du(30)=0;
        end
        
        if ThI==0
            du(41)=0;
            du(42)=0;
            du(43)=0;
        end
        
        %Extrinsic tenase
        %%%%%%%%%%%%%%%%%
        %TF(34), VII(36), VIIa(35), TF:VIIa(12), TF:VII(37), TF:VII:X(28), TF:VII:Xa(29),  TF:VII:IX(38)
        
        %TF
        du(34)= - a2*u(34)*u(36) + a1*u(37) - a4*u(34)*u(35) + a3*u(12);
        
        %VII
        du(36)= - a2*u(34)*u(36) + a1*u(37) - a5*u(12)*u(36) - a6*u(13)*u(36) - a7*u(16)*u(36);
        
        %TF:VII
        du(37)= -a1*u(37) + a2*u(34)*u(36);
        
        %VIIa
        du(35)= - a4*u(34)*u(35) + a3*u(12) + a5*u(12)*u(36) + a6*u(13)*u(36) + a7*u(16)*u(36);
        
        %TF:VIIa Extrinsic tenase (n)
        du(12)= - a3*u(12) +a4*u(34)*u(35) -a9*u(12)*u(3) +a8*u(28) - a12*u(12)*u(13) +a11*u(29) - a14*u(12)*u(8) + a13*u(38) + a15*u(38) - a37*u(12)*u(31) - a42*u(12)*u(6);
        
        %TF:VIIa:X
        du(28)=a9*u(12)*u(3) -a10*u(28) -a8*u(28);
        
        %TF:VIIa:Xa
        du(29)= a10*u(28) +a12*u(12)*u(13) -a11*u(29) -a36*u(29)*u(32) +a35*u(33);
        
        %TF:VIIa:IX
        du(38)= a14*u(12)*u(8) -a13*u(38) -a15*u(38);
        
        
        %TFPI Intercations
        %%%%%%%%%%%%%%%%%%
        %TFPI(32), Xa:TFPI(31), VIIa:TF:Xa:TFPI(33)
        
        %TFPI
        du(32)= - a34*u(13)*u(32) + a33*u(31) - a36*u(29)*u(32) + a35*u(33);
        
        %b_xa p
        du(31)=a34*u(13)*u(32) - a33*u(31) -a37*u(12)*u(31);
        
        %VIIa:TF:Xa:TFPI
        du(33)=a36*u(29)*u(32) -a35*u(33) +a37*u(12)*u(31);
        
        
        %ATIII
        %%%%%%
        %AT(6), Xa:AT(22), mIIa:AT(25), IXa:AT(26), IIa:AT(24), TF:VIIa:AT(7)
        
        %AT
        du(6)=-a38*u(13)*u(6) -a39*u(17)*u(6) -a40*u(18)*u(6) -a41*u(16)*u(6) -a42*u(12)*u(6);
        
        %Xa:AT
        du(22)=a38*u(13)*u(6);
        
        %mIIa:AT
        du(25)=a39*u(17)*u(6);
        
        %IXa:AT
        du(26)=a40*u(18)*u(6);
        
        %IIa:AT
        du(24)=a41*u(16)*u(6);
        
        %TF:VIIa:AT (Brummel eq 34)
        du(7)=a42*u(12)*u(6);
        
        
        %Common pathway
        %%%%%%%%%%%%%%%
        %x(3), V(4), II(5), Xa(13), Va(14), Xa:Va(15), Xa:Va:II(10), IIa(16), mIIa(17)
        
        %X
        du(3)= -a9*u(12)*u(3) +a8*u(28) -a21*u(20)*u(3) +a20*u(27) +a25*u(27) -a43*u(18)*u(3);
        
        %V
        du(4)= - a26*u(16)*u(4) - a44*u(17)*u(4);
        
        %prothrombin (Brummel eq 14)
        du(5)= - a16*u(13)*u(5) - a30*u(15)*u(5) + a29*u(10);
        
        %Xa
        du(13)= -a12*u(12)*u(13) +a11*u(29) -a28*u(13)*u(14) +a27*u(15) +a22*u(27) -a34*u(13)*u(32) +a33*u(31) -a38*u(13)*u(6) +a43*u(18)*u(3) - XaI*1e8*u(13)*u(11) + XaI*0.04*u(30);
        
        %Va
        du(14)= a26*u(16)*u(4) - a28*u(13)*u(14) + a27*u(15) + a44*u(17)*u(4);
        
        %Prothrombinase
        du(15)= a28*u(13)*u(14) - a27*u(15) - a30*u(15)*u(5) + a29*u(10) + a31*u(10) - XaI*1e8*u(15)*u(11) + XaI*0.21*u(1);
        
        %Xa:Va:II
        du(10)=a30*u(15)*u(5) -a29*u(10) -a31*u(10);% - XaI*1e8*u(10)*u(11) + XaI*0.21*u(2);
        
        %Thrombin
        du(16)= a16*u(13)*u(5) +a32*u(17)*u(15) -a41*u(16)*u(6) - ThI*2.5*1e6*u(41)*u(16) + 0.05*u(42);
        
        %mIIa
        du(17)=a31*u(10) -a32*u(17)*u(15) -a39*u(17)*u(6) - ThI*8.5*1e4*u(41)*u(17) - 0.05*u(43);
        
        
        %Intrinsic pathway
        %%%%%%%%%%%%%%%%%%
        %IX(8), VIII(9), IXa(18), VIIIa(19), VIIIa1(21), VIIIa2(23), VIIIa:IXa(20) VIIIa:IXa:X(27)
        
        %IX
        du(8)= - a14*u(12)*u(8) + a13*u(38);

        du(9)= - a17*u(16)*u(9);
        
        du(18)= a15*u(38) -a19*u(19)*u(18) +a18*u(20) +a25*u(27) +a25*u(20) -a40*u(18)*u(6);
        
        %VIIIa
        du(19)= a17*u(16)*u(9) - a19*u(19)*u(18) + a18*u(20) - a24*u(19) + a23*u(21)*u(23);
        
        %VIIIa1
        du(21)=a24*u(19)+a25*u(27)+a25*u(20)-a23*u(21)*u(23);
        
        %VIIIa2
        du(23)=a24*u(19)+a25*u(27)+a25*u(20)-a23*u(21)*u(23);
        
        %Intrinsic tenase (q)
        du(20)= a19*u(19)*u(18) - a18*u(20) - a21*u(20)*u(3) + a20*u(27) + a22*u(27) - a25*u(20);
        
        %intrinsic tenase bound to X
        du(27)=a21*u(20)*u(3) - a20*u(27) - a22*u(27) - a25*u(27);
        
        %Calculate ETP
        %%%%%%%%%%%%%%
        
        du(39)=u(16)+u(17);
        
        du(40)=u(13);
        
        %Factor Xa inhibitor (Rivaroxaban)
        du(11)= - XaI*1e8*u(13)*u(11) + XaI*0.04*u(30) - XaI*1e8*u(15)*u(11) + XaI*0.21*u(1);% - XaI*1e8*u(10)*u(11) + XaI*0.21*u(2);
        du(30)= XaI*1e8*u(13)*u(11) - XaI*0.04*u(30);
        du(1)=XaI*1e8*u(15)*u(11) - XaI*0.21*u(1);
        du(2)=0;%XaI*1e8*u(10)*u(11) - XaI*0.21*u(2);
        
        %Factor IIa inhibitor (DAPA)
        du(41)= - ThI*2.5*1e6*u(41)*u(16) + 0.05*u(42) - ThI*8.5*1e4*u(41)*u(17) + 0.05*u(43);
        du(42)= ThI*2.5*1e6*u(41)*u(16) - 0.05*u(42);
        du(43)= ThI*8.5*1e4*u(41)*u(17) - 0.05*u(43);
        
    end

end

