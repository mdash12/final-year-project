clear;


noOfGames=textread('F:\\project8\\GT Charge\\Put_NO_OF_GAMES.txt','%d');
noOfNodes=textread('F:\\project8\\node location\\Put_NO_OF_NODES.txt','%d');


axis([0 noOfGames 0 100]);
hold on;


distanceGTCharge=textread('F:\\project8\\GT Charge\\Put_DISTANCE_TRAVELLED.txt','%f');
distanceGTChargeModified=textread('F:\\project8\\GT Charge Modified\\Put_DISTANCE_TRAVELLED.txt','%f');
distanceCCGK=textread('F:\\project8\\CCGK\\Put_DISTANCE_TRAVELLED.txt','%f');

plot(distanceGTCharge,'r-','LineWidth',1.5);
plot(distanceGTChargeModified,'g-','LineWidth',0.5);
plot(distanceCCGK,'b-','LineWidth',0.5);


   