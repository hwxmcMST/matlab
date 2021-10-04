function [model,distType,distPara,nLSF]=exp_in(exID)

if exID==1
   
    nLSF=2;               % Number of limit state functions
    model=['exp1_1';'exp1_2';];
    distType=char('norm', 'Lognormal');
    distPara=[4, 0.7;4, 1];
    
    
elseif exID==2
    model=['exp2_1';'exp2_2';'exp2_3'];
    nLSF=3;
    distType=char('norm','norm','norm','norm','Lognormal','Lognormal','norm','norm');
    distPara=[14000,   14000*0.1;
              12,       12*0.01;
              9.00e-4,  9.00e-4*0.06;
              0.05,     0.05*0.08;
              2e11,     2e11*0.06;
              3e10,     3e10*0.06;
              3.35e8    3.35e8*0.18;
              1.34e7    1.34e7*0.18;];
   
elseif exID==3
    nLSF=3;               % Number of limit state functions
    model=['exp3_1';'exp3_2';'exp3_3'];
    distType=char('weib','norm','Lognormal','norm','norm','norm','norm','norm','norm','norm');
    distPara=[2E3,2E2;
        3E4,3E3;
        18E3,2E3;
        3E4,1E3;
        3E4,1E3;
        2E4,1E3;
        2E4,1E3;
        1E3,10;
        4.5E7,4.5E6;
        3.5E6,5E5];
 
elseif exID==4
    model = ['exp4_01';'exp4_02';'exp4_03';'exp4_04';'exp4_05';'exp4_06';'exp4_07';...
        'exp4_08';'exp4_09';'exp4_10'];
    nLSF = 10;
    distType = char('norm','norm','norm','norm','norm','norm','norm','norm'...
        ,'norm','norm','norm');
    distPara =  [0.5,   0.03;
                 1.31,  0.03;
                 0.5,   0.03;
                 1.395, 0.03;
                 0.875, 0.03;
                 1.2,   0.03;
                 0.4,   0.03;
                 0.345,  0.006;
                 0.192,  0.006;
                 0,    10;
                 0,    10;];
elseif exID==5             
    model = ['exp5_1'; 'exp5_2'; 'exp5_3';'exp5_4';'exp5_5';'exp5_6';'exp5_7'];
    nLSF = 7;
    distType = char('norm','Lognormal','Lognormal','Lognormal','Lognormal','Lognormal','Lognormal',...
       'norm','norm','norm','norm','norm','norm',...
        'norm','norm','norm','norm','norm','norm','norm','norm','norm');
    distPara =  [
                 16E3,    16E3*0.1;            % F
                 300E6,   300E6*0.19;          % plate strenght 
                 300E6,   300E6*0.19;          % plate strenght 
                 310E6,   310E6*0.19;           % bolt strenght
                 310E6,   310E6*0.19;           % bolt strenght 
                 310E6,   310E6*0.19;           % bolt strenght 
                 310E6,   310E6*0.19;           % bolt strenght 
                 10E-3,   10E-3*0.02;          % t1
                 10E-3,   10E-3*0.02;          % t2
                 16E-3,   16E-3*0.02;          % bolt length
                 16E-3,   16E-3*0.02;          % bolt length
                 16E-3,   16E-3*0.02;          % bolt length
                 16E-3,   16E-3*0.02;          % bolt length
                 320E-3,  320E-3*0.02;         % L1
                 50E-3,   50E-3*0.02;          % L2
                 200E-3,  200E-3*0.02;         % L3
                 75E-3,   75E-3*0.02;          % L4
                 60E-3,   60E-3*0.02;          % L5
                 144E-6,  144E-6*0.02;         % As
                 144E-6,  144E-6*0.02;         % As
                 144E-6,  144E-6*0.02;         % As
                 144E-6,  144E-6*0.02;];       % As
             
end