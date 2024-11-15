% equivalent to view and commit
system("git add *");
msg = input("commit message: ",'s');
system("git commit -m "" "+ msg + " ""  ");

% equivalent to !git push in terminal
system("git push");

%%% Note: type '!git pull' in terminal to pull