system("git add *");
msg = input("commit message: ",'s');
system("git commit -m "" "+ msg + " ""  ");
system("git push");