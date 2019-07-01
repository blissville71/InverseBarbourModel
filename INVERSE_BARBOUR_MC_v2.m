function INVERSE_BARBOUR_MC_v2(INPUT,OUTPUT,SIMNUM,MACHINE_ERROR)
%
%Code name: INVERSE_BARBOUR_MC_v2
%Authors: Michael Singer; Cristina Evans; Chris Sargeant, last updated, 2019
%Date created: 13/11/2014
%
%Inverse Barbour model to compute source water based on measured tree ring
%cellulose with error bars determined by Monte Carlo sampling of climate
%and physiological variables.
%
%Based on Barbour, M. M., J. S. Roden, G. D. Farquhar, and J. R. Ehleringer. 2004. 
%Expressing leaf water and cellulose oxygen isotope ratios as enrichment above source
%water reveals evidence of a Péclet effect. Oecologia 138:426-435, 10.2307/40005850.
%
%USAGE: INVERSE_BARBOUR_MC_v1(INPUT,OUTPUT,SIMNUM,MACHINE_ERROR)
%
%INPUT is a string containing the name of the input textfile
%The input file should contain 12 columns, most of which are paired columns of means and standard deviations of variables used for Monte Carlo simulation:
%sample code, measured cellulose d18O (‰), atmospheric water vapour d18O (‰), atmospheric water vapour d18O SD (‰), RH (%), RH SD (%), airT(°C), airT SD (°C), stomatal conductance (mol m-2 s-1), stomatal conductance SD (mol m-2 s-1), transpiration rate (mmol m-2 s-1), transpiration rate SD (mmol m-2 s-1)
 
%
%OUTPUT is a string containing the name of the output textfile
%
%SIMNUM is number of Monte Carlo simulations to run by varying input parameters
%
%MACHINE_ERROR is decimal containing the average analytical error report from the IRMS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%model parameters

%number of Monte Carlo simulations
sim_num = SIMNUM;
machine_error = MACHINE_ERROR; %enter the appropriate value here for analytical error based on machine output
%confidence_lim = 1.96; %based on 95% confidence. use 2.575 for 99%; 1.96 for 95%; 1.645 for 90%

%Standard Mean Ocean Water d18O
SMOW=0.0020052;

%concentration of d18O in water (mol m-3)
concwater=55500;

%diffusivity of d18O in water (m2 s-1)
diffhwater=2.66e-9;

%diffusive fractionation of d18O through stomata (per mil)
diffstomata=28;%changed 32 to 28 following Luz 2009 study that merlivat values more accurate than cappa 2003

%diffusive fractionation through boundary layer (per mil)
diffboundary=19;%changed from 21 to 19 as above

%equilibrium fractionation between C=O and water (per mil)
equilfrac=27;
%----------------------------------------------------------------------------
%open input table

A=importdata(INPUT);
data=A.data;
header=A.textdata;
%----------------------------------------------------------------------------
%determine the accuracy of the approximation

choose = 2;
chosenaccuracy=choose+2;
for p=1:chosenaccuracy
    accuracy(p)=10^(p-1);
end
[j,~] = size(data);

%assign site number

for n=1:j
   
    %choose a starting value for the source water
    %lets make it the actual cellulose value
    
    source=data(n,1);
    
    %measured d18O in tree ring cellulose
    actualcellulose=data(n,1);
    
    %atmospheric water vapour d18O (per mil)
    atmwat_mu=data(n,2);
    atmwat_sd=data(n,3);
    atmwat_PD = makedist('normal','mu',atmwat_mu,'sigma',atmwat_sd); %create normal distribution from mean and SD of measured values that is sampled below
      
    %relative humidity (%)
    relhum_mu=data(n,4);
    relhum_sd=data(n,5);
    relhum_PD = makedist('normal','mu',relhum_mu,'sigma',relhum_sd);

    %air temp (oC)
    airtemp_mu=data(n,6);
    airtemp_sd=data(n,7);
    airtemp_PD = makedist('normal','mu',airtemp_mu,'sigma',airtemp_sd);

    %leaf temperature difference (oC)
    leafdiff=1;
    
    %stomatal conductance (mol m-2 s-1)
    stomcond_mu=data(n,8);
    stomcond_sd=data(n,9);
    stomcond_PD = makedist('normal','mu',stomcond_mu,'sigma',stomcond_sd);

   
    %transpiration rate (mmol m-2 s-1)
    transp_mu=data(n,10);
    transp_sd=data(n,11);
    transp_PD = makedist('normal','mu',transp_mu,'sigma',transp_sd);
    
    %proportion of exchangeable oxygen in cellulose via Roden et al., 2000
    exchoxy=0.42;
    
    %proportion of xylem water in meristem
    propmeristem=1;
    %----------------------------------------------------------------------------
    %calculate new values under Monte Carlo scheme sampling climatic and
    %tree variables randomly for # of simulations defined in sim_num in order to create error
    %estimates around our computed values of back-calculated source waters
    
    for MC = 1:1:sim_num
        
        transp = random(transp_PD,1);
                
        stomcond = random(stomcond_PD,1);
                        
        airtemp = random(airtemp_PD,1);
                
        relhum = random(relhum_PD,1);
                
        atmwat = random(atmwat_PD,1);
        
        %effective path length (m) via Song et al., 2014
        transp2 = transp/1000;
        %pathlength = 0.0000236*transp2^(-1.2);
        path_mu=0.0000236*transp2^(-1.2);
        path_sd = 0.19*path_mu; %this error corresponds to the R2 value reported in the paper (0.81).
        path_PD = makedist('normal','mu',path_mu,'sigma',path_sd);
        pathlength = random(path_PD,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UNFINISHED
        %boundary layer conductance (mol m-2 s-1)
        %boundcond=1;
        boundcond_mu = 0.75;
        boundcond_sd = 0.15;
        boundcond_PD = makedist('normal','mu',boundcond_mu,'sigma',boundcond_sd);
        boundcond = random(boundcond_PD,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %leaf temperature
        leaftemp=airtemp+leafdiff;
        
        %equilibrium vapour pressure fractionation (per mil) Majoube,1971
        equilvapfrac=((exp((1137/((airtemp+273.16)^2))-(0.4156/(airtemp+273.16))-0.0020667))-1)*1000;
        
        %total diffusive fractionation (per mil)
        totdifffrac=(((1/stomcond)*diffstomata)+((1/boundcond)*diffboundary))/((1/boundcond)+(1/stomcond));
        
        %18O/16O water vapour (per mil)
        ratiovapour=((atmwat/1000)+1)*SMOW;
        
        %18O/16O source water (per mil)
        ratiosource=((source/1000)+1)*SMOW;
        
        %Delta18O water vapour (per mil)
        Deltavapour=((ratiovapour/ratiosource)-1)*1000;
        
        %saturated vapour pressure of the air (mbar)
        satvappress=(6.13753*exp(airtemp*((18.564-(airtemp/254.4)))/(airtemp+255.57)));
        
        %actual vapour pressure of the air (mbar)
        actvappress=(relhum/100)*satvappress;
        
        %leaf internal vapour pressure (mbar)
        leafvappress=(6.13753*exp(leaftemp*((18.564-(leaftemp/254.4)))/(leaftemp+255.57)));
        
        %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
        CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
        
        %Peclet number
        peclet=((transp/1000)*pathlength)/(concwater*diffhwater);
        
        %bulk leaf water (Delta18O) (per mil)
        bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
        
        %sucrose Delta18O (per mil)
        %sucroseD=bulkleafD+equilfrac;
        
        %cellulose Delta18O (per mil)
        celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
        
        %Craig-Gordon evaporation site water (d18Oe) (per mil)
        %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
        
        %bulk leaf water (d18O) (per mil)
        %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
        
        %sucrose d18O (per mil)
        %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
        
        %cellulose d18O (per mil)
        cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
        %----------------------------------------------------------------------------
        
        for y=1:chosenaccuracy
            
            %create a loop to calculate a more accurate value
            if mod(y,2)~=0
                while cellulose*accuracy(y)>actualcellulose*accuracy(y)
                    source=source-(1/accuracy(y));
                    %18O/16O source water (per mil)
                    ratiosource=((source/1000)+1)*SMOW;
                    
                    %Delta18O water vapour (per mil)
                    Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                    
                    %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                    CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                    
                    %bulk leaf water (Delta18O) (per mil)
                    bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                    
                    %sucrose Delta18O (per mil)
                    %sucroseD=bulkleafD+equilfrac;
                    
                    %cellulose Delta18O (per mil)
                    celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                    
                    %Craig-Gordon evaporation site water (d18Oe) (per mil)
                    %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %bulk leaf water (d18O) (per mil)
                    %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %sucrose d18O (per mil)
                    %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %cellulose d18O (per mil)
                    cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                end
            else
                while cellulose*accuracy(y)<actualcellulose*accuracy(y)
                    source=source+(1/accuracy(y));
                    
                    %18O/16O source water (per mil)
                    ratiosource=((source/1000)+1)*SMOW;
                    
                    %Delta18O water vapour (per mil)
                    Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                    
                    %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                    CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                    
                    %bulk leaf water (Delta18O) (per mil)
                    bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                    
                    %sucrose Delta18O (per mil)
                    %sucroseD=bulkleafD+equilfrac;
                    
                    %cellulose Delta18O (per mil)
                    celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                    
                    %Craig-Gordon evaporation site water (d18Oe) (per mil)
                    %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %bulk leaf water (d18O) (per mil)
                    %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %sucrose d18O (per mil)
                    %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                    
                    %cellulose d18O (per mil)
                    cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                end
            end
        end
        
        %----------------------------------------------------------------------------
        %store values to plot in a graph
        
        sourcestor(MC,n)=source;
        
    end
    std_source(n) = std(sourcestor(:,n));
    [~,~] = size(sourcestor);
    %----------------------------------------------------------------------------
    %redo the calculations for the measured (mean) values for climatic and
    %physiological variables
    
    atmwat=atmwat_mu;
    transp=transp_mu;
    airtemp=airtemp_mu;
    stomcond=stomcond_mu;
    relhum=relhum_mu;
    %barpress=barpress_real;
    
    %effective path length (m) via Song et al., 2014
    pathlength=0.0000236*transp^(-1.2);
    
    %leaf temperature
    leaftemp=airtemp+leafdiff;
    
    %equilibrium vapour pressure fractionation (per mil)
    equilvapfrac=((exp((1137/((airtemp+273.16)^2))-(0.4156/(airtemp+273.16))-0.0020667))-1)*1000;
    
    %total diffusive fractionation (per mil)
    totdifffrac=(((1/stomcond)*diffstomata)+((1/boundcond)*diffboundary))/((1/boundcond)+(1/stomcond));
    
    %18O/16O water vapour (per mil)
    ratiovapour=((atmwat/1000)+1)*SMOW;
    
    %18O/16O source water (per mil)
    ratiosource=((source/1000)+1)*SMOW;
    
    %Delta18O water vapour (per mil)
    Deltavapour=((ratiovapour/ratiosource)-1)*1000;
    
    %saturated vapour pressure of the air (mbar)
    satvappress=(6.13753*exp(airtemp*((18.564-(airtemp/254.4)))/(airtemp+255.57)));
    
    %actual vapour pressure of the air (mbar)
    actvappress=(relhum/100)*satvappress;
    
    %leaf internal vapour pressure (mbar)
    leafvappress=(6.13753*exp(leaftemp*((18.564-(leaftemp/254.4)))/(leaftemp+255.57)));
    
    %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
    CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
    
    %Peclet number
    peclet=((transp/1000)*pathlength)/(concwater*diffhwater);
    
    %bulk leaf water (Delta18O) (per mil)
    bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
    
    %sucrose Delta18O (per mil)
    %sucroseD=bulkleafD+equilfrac;
    
    %cellulose Delta18O (per mil)
    celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
    
    %Craig-Gordon evaporation site water (d18Oe) (per mil)
    %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
    
    %bulk leaf water (d18O) (per mil)
    %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
    
    %sucrose d18O (per mil)
    %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
    
    %cellulose d18O (per mil)
    cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
    %----------------------------------------------------------------------------
    
    for y=1:chosenaccuracy
        
        %create a loop to calculate a more accurate value
        if mod(y,2)~=0
            while cellulose*accuracy(y)>actualcellulose*accuracy(y)
                source=source-(1/accuracy(y));
                %18O/16O source water (per mil)
                ratiosource=((source/1000)+1)*SMOW;
                
                %Delta18O water vapour (per mil)
                Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                
                %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                
                %bulk leaf water (Delta18O) (per mil)
                bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                
                %sucrose Delta18O (per mil)
                %sucroseD=bulkleafD+equilfrac;
                
                %cellulose Delta18O (per mil)
                celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                
                %Craig-Gordon evaporation site water (d18Oe) (per mil)
                %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %bulk leaf water (d18O) (per mil)
                %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %sucrose d18O (per mil)
                %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %cellulose d18O (per mil)
                cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
            end
        else
            while cellulose*accuracy(y)<actualcellulose*accuracy(y)
                source=source+(1/accuracy(y));
                
                %18O/16O source water (per mil)
                ratiosource=((source/1000)+1)*SMOW;
                
                %Delta18O water vapour (per mil)
                Deltavapour=((ratiovapour/ratiosource)-1)*1000;
                
                %Craig-Gordon evaporation site water (Delta18Oe) (per mil)
                CGevapwaterD=(((1+(equilvapfrac/1000))*(1+(totdifffrac/1000)+(((Deltavapour/1000)-(totdifffrac/1000))*(actvappress/leafvappress))))-1)*1000;
                
                %bulk leaf water (Delta18O) (per mil)
                bulkleafD=(CGevapwaterD*(1-(exp(-peclet))))/peclet;
                
                %sucrose Delta18O (per mil)
                %sucroseD=bulkleafD+equilfrac;
                
                %cellulose Delta18O (per mil)
                celluloseD=bulkleafD*(1-(propmeristem*exchoxy))+equilfrac;
                
                %Craig-Gordon evaporation site water (d18Oe) (per mil)
                %CGevapwater=((((CGevapwaterD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %bulk leaf water (d18O) (per mil)
                %bulkleaf=((((bulkleafD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %sucrose d18O (per mil)
                %sucrose=((((sucroseD/1000)+1)*ratiosource)/SMOW-1)*1000;
                
                %cellulose d18O (per mil)
                cellulose=((((celluloseD/1000)+1)*ratiosource)/SMOW-1)*1000;
            end
        end
    end
    sourcestore(n)=source;
    cellulosestore(n)=cellulose;
    actualcellulosestore(n)=actualcellulose;
    
end
%----------------------------------------------------------------------------
%save the output

source = sourcestore;
upper_source_std = source+abs(std_source)+machine_error; %added in the average analytical error here defined above.
lower_source_std = source-abs(std_source)-machine_error;
Out=[sourcestore;upper_source_std;lower_source_std];


%plot the results (with SD error bounds)
figure
x = 1:j;
X=[x,fliplr(x)];
Y=[upper_source_std,fliplr(lower_source_std)];
fill(X,Y,'g')
hold on
h = plot(source,'r');
set(h, 'LineWidth', 0.5)
ylabel('Modeled Source Water d18O (per mil)');
xlabel('Samples plotted as a time series');
hold off
M=Out';

%open the output file
inverseoutput=fopen(OUTPUT,'w');
fprintf(inverseoutput,'%11s\t %25s\t %21s\t %21s\t %22s\t %28s\n','Sample_Code','Modeled_Source_Water_d18O','Source_Water_SD+_d18O', 'Source_Water_SD-_d18O', 'Modeled Tree Ring d18O', 'Measured Tree d18O');

for row = 1:j
    fprintf(inverseoutput,'%s\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\n',header{row+1,1}, M(row,:), cellulosestore(row), actualcellulosestore(row));
end

fclose(inverseoutput);

