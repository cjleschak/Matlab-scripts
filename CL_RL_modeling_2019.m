%%
%Based on parameter sampling distributions in Davidow 2016 (Neuron) paper
close

%Number of random parameters to draw for each
numStarts=20;
itersPerSub=2000;

%Distribution parameters
betaShape= 1.1;
betaScale= 1.1;

gamShape=1.2;
gamScale=5;

%Plot distributions parameters are sampled from?
plotDists=1;


%%%%%%%%%%%%%%%%%%%%%%%%%

%used for plotting only
X_bet=0:0.001:1;    %LR; note P(selecting estimate of 0 or 1)=0.
X_gam=0:0.001:30;   %IT; note P(selecting an estimate >30) < .001 under this distribution.

%Learning rate (alphaParam) distribution (orig 1.1, 1.1)
betaDist=makedist('beta',betaShape,betaScale);
betaPlot=betapdf(X_bet',betaShape,betaScale);


%Temperature (betaParam) distribution (orig. 1.2, 5)
gamDist=makedist('gamma',gamShape,gamScale);
gamPlot=gampdf(X_gam',gamShape,gamScale);

%Plot the two sampling distributions
if plotDists==1
    plot(X_bet,betaPlot,'LineWidth',2);
    hold on
    plot(X_gam,gamPlot,'LineWidth',2);
    hold on
    legend('Learning rate','Temperature');
    axis([0 5 0 1.2]);
    ylabel('Density');
    xlabel('Possible Values');
    title('Sampling Distributions for Parameters');

    %Plot histogram of actually sampled parameters.
    %histogram(randTemp_params,500,'Normalization','pdf')
    %histogram(randLR_params,12,'Normalization','pdf')
end



%Read data.
fname='/Users/cjleschak/Documents/UCLA/Lab/Lab Resources/RL Workshop Duke/SSSNN/sampData_CJL.csv';
%fname='/Users/cjleschak/Documents/UCLA/Lab/Lab Resources/RL Workshop Duke/SSSNN/SampleData1_forCJLScript.csv';
data = csvread([fname],1); %Import data, ignoring column headers 
%trial num, cue, flower chosen, correct? = csv column headings


%How many different cues?
nCues=length(unique(data(:,2)));
%nTrialsPerCue=length(data)/nCues; %(assumes equal num of trials per cue).
trialfit = 1;
partialfit = 0;
fit=0;
bestfit=10000;

for iter=1:itersPerSub
    %Set priors.
    expected1=zeros(15,nCues);   %starting expected value (EV) for option 1 for each cue (15 trials per cue)
    expected2=zeros(15,nCues);   %starting expected value (EV) for option 2 for each cue (15 trials per cue)
    pr1 = zeros(15,nCues);       %starting probability of choosing option 1 for each cue (15 trials per cue)
    pr2 = zeros(15,nCues);       %starting probability of choosing option 2 for each cue (15 trials per cue)

    randLR_params=random(betaDist,numStarts,1);  %draw random numbers from distribution.
    randTemp_params=random(gamDist,numStarts,1); %draw random numbers from distribution.

    
    %For each set of parameters...
    for n_param=1:length(numStarts)
        cur_LR=randLR_params(numStarts,1);
        cur_TEMP=randTemp_params(numStarts,1);

        cueTrialNum=zeros(1,nCues);

        %For each trial...
        for trialnum = 1:length(data)
            thisCue=data(trialnum,2);
            flowerChoice = data(trialnum,3);
            reward=data(trialnum,4); %outcome on this trial (1=reward; 0=no reward)

            cueTrialNum(thisCue)=cueTrialNum(thisCue)+1;

            if thisCue<3
                optimalChoice=1;
            else
                optimalChoice=2;
            end

            %softmax function
            pr1(cueTrialNum(thisCue),thisCue)= exp(expected1(cueTrialNum(thisCue),thisCue)/cur_TEMP)/(exp(expected1(cueTrialNum(thisCue),thisCue)/cur_TEMP) + exp(expected2(cueTrialNum(thisCue),thisCue)/cur_TEMP));
            pr2(cueTrialNum(thisCue),thisCue)= 1 - pr1(cueTrialNum(thisCue),thisCue);

            %Update expected values
            if flowerChoice==1
              trialfit(trialnum) = pr1(cueTrialNum(thisCue),thisCue);     %probability of choosing option 1 for this cue
              expected1(cueTrialNum(thisCue)+1,thisCue) = expected1(cueTrialNum(thisCue),thisCue) + cur_LR * (reward-expected1(cueTrialNum(thisCue),thisCue)); %EV for option 1, in context of this cue
              expected2(cueTrialNum(thisCue)+1,thisCue) = expected2(cueTrialNum(thisCue),thisCue); %keep EV of unchosen option same.
           elseif flowerChoice==2
              trialfit(trialnum) = pr2(cueTrialNum(thisCue),thisCue);     %probability of choosing option 2 for this cue
              expected2(cueTrialNum(thisCue)+1,thisCue) = expected2(cueTrialNum(thisCue),thisCue) + cur_LR * (reward-expected2(cueTrialNum(thisCue),thisCue)); %EV for option 2, in context of this cue
              expected1(cueTrialNum(thisCue)+1,thisCue) = expected1(cueTrialNum(thisCue),thisCue); %keep EV of unchosen option same.
            end

           if flowerChoice==optimalChoice
               optimalSelected(cueTrialNum(thisCue),thisCue)=1;
           else
               optimalSelected(cueTrialNum(thisCue),thisCue)=0;
           end


           partialfit=partialfit + log(trialfit(trialnum)); %By end of loop, "fit" is fit of all trials summed together.
           fit = -1*partialfit; %we want to find the maximum but using fmin functions, so must multiply by -1

           
        end
    
    if fit<=bestfit       %takes the best fitting parameters if the fit is the best fit so far
        bestfit=fit;
        bestLR=cur_LR;
        bestTEMP=cur_TEMP;
    end
    end
end
fprintf('Best learning rate: %f\n',bestLR)
fprintf('Best temperature: %f\n',bestTEMP)
fprintf('Fit (lower better): %f\n',bestfit)


for thisCue=1:nCues
    block_optAvg(1,thisCue)=mean(optimalSelected(1:3,thisCue));
    block_optAvg(2,thisCue)=mean(optimalSelected(4:6,thisCue));
    block_optAvg(3,thisCue)=mean(optimalSelected(7:9,thisCue));
    block_optAvg(4,thisCue)=mean(optimalSelected(10:12,thisCue));
    block_optAvg(5,thisCue)=mean(optimalSelected(13:15,thisCue));    
end

figure
hold on
plot(1:length(block_optAvg),block_optAvg(:,1),'LineWidth',1)
plot(1:length(block_optAvg),block_optAvg(:,2),'LineWidth',1)
plot(1:length(block_optAvg),block_optAvg(:,3),'LineWidth',1)
plot(1:length(block_optAvg),block_optAvg(:,4),'LineWidth',1)

block_optAvgAll(1)=mean(block_optAvg(1,:));
block_optAvgAll(2)=mean(block_optAvg(2,:));
block_optAvgAll(3)=mean(block_optAvg(3,:));
block_optAvgAll(4)=mean(block_optAvg(4,:));
block_optAvgAll(5)=mean(block_optAvg(5,:));

plot(1:length(block_optAvg),block_optAvgAll,'LineWidth',4,'Color','k')

legend('Butterfly 1','Butterfly 2', 'Butterfly 3', 'Butterfly 4', 'Average','Location','southeast');

axis([1 5 0 1.05]);
ylabel('% Optimal Chosen');
xlabel('Learning Block');
title('Learning Accuracy');

