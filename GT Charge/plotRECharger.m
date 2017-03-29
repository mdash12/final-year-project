clear;

noOfChargers=textread('Put_NO_OF_CHARGERS.txt','%d');
disp(noOfChargers);


no_of_games=textread('Put_NO_OF_GAMES.txt','%d');
disp(no_of_games);

axis([0 noOfChargers 0 1000]);
hold on;

x=textread('Put_REMENERGY_CHARGER.txt');

no_of_games=no_of_games+1;

for i=1:no_of_games
    temp=x((i-1)+1:i,:);
    k=(i-1)*1/no_of_games;
    plot(temp,'-','Color',[k 1-k k]);
    hold on;
       
end
   