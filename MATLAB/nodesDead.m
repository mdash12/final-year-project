clear;


noOfGames=textread('F:\\project8\\GT Charge\\Put_NO_OF_GAMES.txt','%d');
noOfNodes=textread('F:\\project8\\node location\\Put_NO_OF_NODES.txt','%d');


axis([10 noOfGames 0 noOfNodes]);
hold on;


distanceGTCharge=textread('F:\\project8\\GT Charge\\Put_NO_OF_DEAD_NODES.txt','%d');
distanceGTChargeModified=textread('F:\\project8\\GT Charge Modified\\Put_NO_OF_DEAD_NODES.txt','%d');
distanceCCGK=textread('F:\\project8\\CCGK\\Put_NODES_DEAD.txt','%d');

plot(distanceGTCharge,'r-','LineWidth',1.5);
plot(distanceGTChargeModified,'m-','LineWidth',0.5);
plot(distanceCCGK,'b-','LineWidth',0.5);


   