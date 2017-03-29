clear;


noOfGames=textread('F:\\project8\\GT Charge\\Put_NO_OF_GAMES.txt','%d');
noOfNodes=textread('F:\\project8\\node location\\Put_NO_OF_NODES.txt','%d');


axis([0 noOfGames 0 10]);
hold on;


avgEnergyTransferredGTCharge=textread('F:\\project8\\GT Charge\\Put_AVG_ENERGY_TRANSFERRED.txt','%f');
avgEnergyTransferredGTChargeModified=textread('F:\\project8\\GT Charge Modified\\Put_AVG_ENERGY_TRANSFERRED.txt','%f');
avgEnergyTransferredCCGK=textread('F:\\project8\\CCGK\\Put_AVG_ENERGY_TRANSFERRED.txt','%f');

plot(avgEnergyTransferredGTCharge,'r-','LineWidth',1.5);
plot(avgEnergyTransferredGTChargeModified,'g-','LineWidth',0.5);
plot(avgEnergyTransferredCCGK,'b-','LineWidth',0.5);


   