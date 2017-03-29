clear;

no_of_games=textread('Put_NO_OF_GAMES.txt','%d');
disp(no_of_games);



no_of_nodes=textread('Put_NO_OF_NODES.txt','%d');
disp(no_of_nodes);

axis([0 no_of_games 0 no_of_nodes]);
hold on;


requests=textread('Put_REQUESTS.txt','%d');
reqServed=textread('Put_REQ_SERVED.txt','%d');

disp(requests);

plot(requests,'g-','LineWidth',0.5);
hold on;
plot(reqServed,'b-','LineWidth',0.5);
