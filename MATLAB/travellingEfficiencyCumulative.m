clear;


noOfGames=textread('F:\\project8\\GT Charge\\Put_NO_OF_GAMES.txt','%d');
noOfNodes=textread('F:\\project8\\node location\\Put_NO_OF_NODES.txt','%d');


axis([0 noOfGames 0 2]);
hold on;


distanceGTCharge=textread('F:\\project8\\GT Charge\\Put_TRAVELLING_EFFICIENCY_CUMULATIVE.txt','%f');
distanceGTChargeModified=textread('F:\\project8\\GT Charge Modified\\Put_TRAVELLING_EFFICIENCY_CUMULATIVE.txt','%f');
distanceCCGK=textread('F:\\project8\\CCGK\\Put_TRAVELLING_EFFICIENCY_CUMULATIVE.txt','%f');

plot(distanceGTCharge,'r-','LineWidth',0.5);
plot(distanceGTChargeModified,'g-','LineWidth',1);
plot(distanceCCGK,'b-','LineWidth',1);

