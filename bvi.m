%% The Biological Value Index (Sanders, 1960; Loya-Salinas & Escofet, 1990) revisited.
%
%
%% Syntax
%
%   BVIs = BVI(x,p)
%   BVIs = BVI(x,p,sp,loc)
%
%   [BVIs, BVIl] = BVI(x,p,...)
%   [BVIs,BVIl,X] = BVI(x,p,...)
%
%   BVIs = BVI(__,plot)
% 
%   [~] = BVI(__,plot);
%
%% Description
%
% BVIs = BVI(x,p) takes an abundance or density table x and a cutoff
% percentage p. and outputs the Biological Value Index (Sanders, 1960;
% Loya-Salinas & Escofet, 1990) for each species that cummulative
% contribute up to p. Output BVIs is a table. x has species as rows and
% samples as columns. P may be expressed as a percentage (0-100) or a ratio
% (0-1).
%
% BVIs = BVI(x,p,sp,loc) takes an abundance or density matrix x, a cutoff
% percentage p, a cell array with species names sp (rows of x), and a cell
% array of sample names loc (columns of x). It first creates a table x and
% then calculates the Index of Bioloical Value (Sanders, 1960; Loya-Salinas
% & Escofet, 1990) for each species that cummulative contribute up to p.
% Output BVIs is a table.
%
% [BVIs, BVIl] = BVI(x,p,...) also outputs a table BVIs witht he scores
% of each species on each sample. BVIl has the same dimensions as x.
%
% [BVIs,BVIl,X] = BVI(x,p,...) also outputs X, the table used to run the
% function. Usefull if inputs were a matrix and two vectors containing
% species and locations. If x was a table, then X is the same table.
%
% BVIs = BVI(__,plot) takes the information either in table or matrix
% (see above) and calculates the Biological Value Index (Sanders,
% 1960; Loya-Salinas & Escofet, 1990) for each species that cummulative
% contribute up to p. plot is a logical of size 1x1. If plot = 1, the
% function also produces a stacked pecentages bar graph.
% 
% [~] = BVI(__,plot) Does not provide any output, and produces a stacked
% percentages bar gaph.
%
%% Examples
%
%   % Data is from Loya-Salinas & Escofet, 1990
%   
%   % sp contains the species names, loc contains the sampling locations
%
%   sp = {'Synchelidium spp.'
%        'Tridentella spp.'
%        'Nerine cirratulus'
%        'Nephtys californiensis'
%        'Glycera tenuis'
%        'Donax gouldii'
%        'Orchestoidea benedicti'
%        'Archaeomysis spp.'
%        'Armadillium spp.'
%        'Megalopus spp.'
%        'Emerita analoga'
%        'Pontharpinia spp.'
%        'Euzonus mucronata'
%        'Lepidopa californica'
%        'Magelona californica'
%        'Hanstorina spp.'
%        'Glycera dibranchiata'
%        'Archaeomysis maculata'};
%
%   loc = {'S1' 'S2' 'S3' 'S4' 'S5' 'S6'};
%
%   % x has the rows of sp and columns of loc
%
%  x = [  2398        1626         811        1275        1343        7079
%         2048        1125         528        1990        1098        1274
%           37         165        1141        1540         118          53
%          544         875         404         170          58          90
%          265         566         106         646         133         118
%          914          75          42           5           5          15
%           11         251         133          79         162         245
%          958          90         160          37         522         111
%           59         155          91          16         208         283
%          149         266          48          30           0           0
%          101          16          96          27           0          10
%           69          96          11          16           0           0
%            0           0          37         341           0           0
%            5          69          16           0          16           0
%            0           5           5           0          42          10
%            0          16           5           0           0           0
%            0           0           0           0           5           0
%            0           0           0           5           0           0];
% 
%   % our cutoff percentage is 95%
%   p=95;
%
% We will calculate the Biological Value Index for the matrix, for each
% sample, and also obtain a table X with data in sp, loc and x in it. We
% will also obtain a plot.
%
%   [BVIs,BVIl,X]= BVI(x,p,sp,loc,1);
%
%   disp(BVIs)
%
%                               BVI     RBVI 
%                               ___    ______
% 
%     Synchelidium spp.         57     17.221
%     Tridentella spp.          54     16.314
%     Glycera tenuis            34     10.272
%     Nephtys californiensis    33     9.9698
%     Archaeomysis spp.         31     9.3656
%     Nerine cirratulus         30     9.0634
%     Orchestoidea benedicti    27     8.1571
%     Armadillium spp.          21     6.3444
%     Megalopus spp.            13     3.9275
%     Donax gouldii              9      2.719
%     Others (8)                22     6.6465
%   
% These are the species that consistently contribute to 95% of the
% abundances. The rest of the species are labeled as Others (8) indicating
% that the 8 spp left only contributed to 5% of the abundances and were
% less consistent.
%
%   disp(BVIl)
%
%                               S1    S2    S3    S4    S5    S6    BVI     RBVI  
%                               __    __    __    __    __    __    ___    _______
% 
%     Synchelidium spp.         10    10     9     8    10    10    57      17.221
%     Tridentella spp.           9     9     8    10     9     9    54      16.314
%     Glycera tenuis             5     7     4     7     5     6    34      10.272
%     Nephtys californiensis     6     8     7     5     3     4    33      9.9698
%     Archaeomysis spp.          8     1     6     3     8     5    31      9.3656
%     Nerine cirratulus          0     4    10     9     4     3    30      9.0634
%     Orchestoidea benedicti     0     5     5     4     6     7    27      8.1571
%     Armadillium spp.           1     3     2     0     7     8    21      6.3444
%     Megalopus spp.             4     6     1     2     0     0    13      3.9275
%     Donax gouldii              7     0     0     0     0     2     9       2.719
%     Emerita analoga            3     0     3     1     0     1     8      2.4169
%     Euzonus mucronata          0     0     0     6     0     0     6      1.8127
%     Pontharpinia spp.          2     2     0     0     0     0     4      1.2085
%     Magelona californica       0     0     0     0     2     1     3     0.90634
%     Lepidopa californica       0     0     0     0     1     0     1     0.30211
%     Glycera dibranchiata       0     0     0     0     0     0     0           0
%     Archaeomysis maculata      0     0     0     0     0     0     0           0
%     Hanstorina spp.            0     0     0     0     0     0     0           0
%
% Here we see the scores for each species and each sample (location). This
% values are the same as the ones published by Loya-Salinas & Escofet,
% (1990; see Table III) and include an extra column with the new %BVI
% (relative BVI) that allows for comparisons.
% 
%   disp(X)
%
%                                S1      S2      S3      S4      S5      S6 
%                               ____    ____    ____    ____    ____    ____
% 
%     Synchelidium spp.         2398    1626     811    1275    1343    7079
%     Tridentella spp.          2048    1125     528    1990    1098    1274
%     Nerine cirratulus           37     165    1141    1540     118      53
%     Nephtys californiensis     544     875     404     170      58      90
%     Glycera tenuis             265     566     106     646     133     118
%     Donax gouldii              914      75      42       5       5      15
%     Orchestoidea benedicti      11     251     133      79     162     245
%     Archaeomysis spp.          958      90     160      37     522     111
%     Armadillium spp.            59     155      91      16     208     283
%     Megalopus spp.             149     266      48      30       0       0
%     Emerita analoga            101      16      96      27       0      10
%     Pontharpinia spp.           69      96      11      16       0       0
%     Euzonus mucronata            0       0      37     341       0       0
%     Lepidopa californica         5      69      16       0      16       0
%     Magelona californica         0       5       5       0      42      10
%     Hanstorina spp.              0      16       5       0       0       0
%     Glycera dibranchiata         0       0       0       0       5       0
%     Archaeomysis maculata        0       0       0       5       0       0
% 
% Created by Villaseñor-Derbez, J.C. Bren School of Environmental Science &
% Management, University of California, Santa Barbara, Santa Barbara, CA
% 93106, U.S.A. <jvillasenor@bren.ucsb.edu>
%
% Copyright (c) November 29, 2015. To cite this file, this would be an
% appropriate format:
% 
% Villaseñor-Derbez, J.C. 2015. The Biological Value Index (Sanders, 1960;
% Loya-Salinas & Escofet, 1990) revisited. A  MATLAB file.  [WWW document].
%
% DOI: 10.13140/RG.2.1.2082.8241
% 
% URL: 
% 
%% References
% [1] Sanders, H.L. 1960. Benthic Studies in Buzzards Bay III. The
% Structure of the Soft-Bottom Comminity. Limnology and Oceanography,
% 5(2): 138-153.
% 
% [2] Loya-Salinas, D.H. & Escofet, A. 1990. Contributions to the
% calculation of the biological value index (Sanders, 1960). Ciencias
% Marinas, 16(2):97-115.
% 
%
function [BVIs,BVIl,X] = bvi(x,p,varargin)

[l,w]=size(x);

%%
% Handling p
% Evaluate if p is in the range 0-1 or 0-100
if p>1
    p=p/100;
end

if p==100;
    error('Cut-off percentage can not be 100%')
end


%%
% Evaluate if x is a table
% If x is nt a table, then it makes sure that cell arrays are provided to
% create a table with row and variable names (species and samples).

if ~istable(x)
    if length(varargin)<2
        error('x is a matrix, not a table. If x is a matrix, you must provide cell array of species and a cell array of samples')
    else
        sp=varargin{1};
        loc=varargin{2};
        
        if ~iscell(sp)
            error('sp must be a cell array')
        end
        
        if ~iscell(loc)
            error('loc must be a cell array')
        end
        
        if length(sp)~=l
            error('The number of rows in x (species) must be the same as the number of elements in sp')
        end
        if length(loc)~=w
            error('The number of columns in x (samples) must be the same as the number of elements in loc')
        end
        X=table(x(:,1));
        for i=2:w
            X{:,i}=x(:,i);
        end
        X.Properties.RowNames=sp;
        X.Properties.VariableNames=loc;
        x=X;
    end
end
% Table has been created and stored in x, given that the user wishes to
% widraw it as an output


%%
% Evaluate for all-zero samples

check=sum((x{:,:}>0));
check2=sum((x{:,:}>0),2);

if any(check==0)
    warning('Your dataset contains at least one all-zero sample')
end

if any(check2==0,2)
    warning('Your dataset contains at least one all-zero species')
end


%%
% Calculating relative abundances
% The cicle creates a mtrix with repeted elements of sum(x)

ntot=zeros(l,w); %Dimensioning variable to save space
for i=1:l
    ntot(i,:)=sum(x{:,:}); %Calculates total abundances for each column
end

prop=x; %prop is now a table with the same size as x
prop{:,:}=x{:,:}./ntot; %We put the relative abundances into the table


%%
% Calculating the sum of relative abundances
% This calculates the sum of relative abundances across samples
t=sum(prop{:,:},2);

count=prop;
count{:,w+1}=t;
str=horzcat('Var', num2str(w+1));
count.Properties.VariableNames{str}='SumNi';

count=sortrows(count,'SumNi','descend');
% Relative abundances are now sorted according to the column with the
% summed abundances


%%
% Calculate the number of species (S) needed to reach p in each sample
% The cicle sorts one column at a time and calculates the number of species
% needed to reach a p value of p, or instantly above it.
y=count{:,1:w};
S=zeros(w,1);
for i=1:w
    j=1;
    sample=sort(y(:,i),'descend');
    pi=sample(j);
    while pi<p
        pi=pi+sample(j+1);
        j=j+1;
    end
    S(i)=j;
end

% S has the number of species in each sample that contribute to p

%%
% Allocation of points
% Allocation of points is done as described by Loya-Salinas & Escofet,
% (1990), taking into account possible ties in relative abundances.

location=count.Properties.VariableNames;
scores=count(:,1:w);

for i=1:w
    scores=sortrows(scores,location(i),'descend');
    column=zeros(l,1);
    un=sort(unique(scores{:,location(i)}),'descend');
    Ss=max(S); % Here Ss contains the largest number of species needed to achieve p
    for j=1:min(length(un)-1,max(S));
        column(scores{:,location(i)}==un(j))=Ss;
        Ss=max(S)-sum(column>0); %This operation modifies Ss to take ties into account
        if Ss<0
            Ss=0;
        end
    end
    scores{:,location(i)}=column;
end


%%
% Handling scores to output BVIs and BVIl

finalscores=sum(scores{:,:},2); % his is the BVI, sum of scores across samples
scoresl=scores; %This contains all the scores
scoresl{:,w+1}=finalscores; %BVI is added to the table
scoresl{:,w+2}=(finalscores*100/sum(finalscores)); %And the relative BVI too
scoresl=sortrows(scoresl,str,'descend');
scoresl.Properties.VariableNames{str}='BVI';
str=horzcat('Var', num2str(w+2));
scoresl.Properties.VariableNames{str}='RBVI';

finalscores=table(finalscores,'RowNames',scores.Properties.RowNames,'VariableNames',{'A'});
finalscores=sortrows(finalscores,'A','descend');

finalscoresS=finalscores{1:max(S),'A'};
finalscoresS(max(S)+1)=sum(finalscores{max(S)+1:end,'A'});
R=finalscoresS/sum(finalscoresS)*100;

species=finalscores.Properties.RowNames(1:max(S));
species(max(S)+1)={horzcat('Others (',num2str(l-max(S)),')')};


%%
% Outputs
% BVIs contains the BVI and Relative BVI scores
BVIs=table(finalscoresS, R,'RowNames',species,'VariableNames',{'BVI','RBVI'});
% BVIl contains scores for all spp and all samples, BVI and relative BVI
BVIl=scoresl;
% X contains the table generated or fed to BVI()
X=x;


%%
% Plotting
% If plotting was indicated it evaluates it in two ways. First, if inputs
% were a table, p and 1 and then if inputs were a matriz, p, sp, locations
% and 1

if length(varargin)==1
    if varargin{1}==1
        Ss=max(S);
        BVIll=zeros(Ss+1,w);
        BVIsl=BVIl{1:Ss,1:w};
        BVIsl(Ss+1,:)=sum(BVIl{Ss+1:end,1:w});
        for i=1:Ss+1
            BVIll(i,:)=sum(BVIsl);
        end
        BVIprop=BVIsl./BVIll;
        figure1=figure;
        Xticks=BVIl.Properties.VariableNames(1:end-1);
        axes1 = axes('Parent',figure1,'XTickLabel',Xticks);
        box(axes1,'on');
        hold(axes1,'all');
        bar(BVIprop',1,'stacked','EdgeColor','k','LineWidth',1.5)
        legend('\it%s',BVIs.Properties.RowNames(:),'Location','EastOutside')
        legend('boxoff')
        axis tight
        ylabel ('% BVI', 'FontSize',12)
        xlabel ('Samples', 'FontSize',12)
    end
elseif length(varargin)==3
    if varargin{3}==3
                Ss=max(S);
        BVIll=zeros(Ss+1,w);
        BVIsl=BVIl{1:Ss,1:w};
        BVIsl(Ss+1,:)=sum(BVIl{Ss+1:end,1:w});
        for i=1:Ss+1
            BVIll(i,:)=sum(BVIsl);
        end
        BVIprop=BVIsl./BVIll;
        figure1=figure;
        Xticks=BVIl.Properties.VariableNames(1:end-1);
        axes1 = axes('Parent',figure1,'XTickLabel',Xticks);
        box(axes1,'on');
        hold(axes1,'all');
        bar(BVIprop',1,'stacked','EdgeColor','k','LineWidth',1.5)
        legend('\it%s',BVIs.Properties.RowNames(:),'Location','EastOutside')
        legend('boxoff')
        axis tight
        ylabel ('% BVI', 'FontSize',12)
        xlabel ('Samples', 'FontSize',12)
    end
end
end