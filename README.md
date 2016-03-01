A LAMMPS fix for parameterizing Tersoff potentials with long-range charge and LJ forces. 

To use, add these files to the LAMMPS src directory and compile. A makefile is provided in MAKE, specialized for the computer systems of the Clancy research group at Cornell. 

Getting started:

1.	Open a terminal in the src directory in your lammps folder

2.	Generate an SSH key and add it to the ssh-agent

		ssh-keygen -t rsa -b 4096 -C "your_email@example.com"  #Creates an ssh key, using your GitHub e-mail as a label
		
	When prompted to "Enter a file in which to save the key," press Enter
	If asked to Overwrite, enter 'y'
	At the prompt, type a secure passphrase
	Retype your secure passphrase
	
3.	Adding a new SSH key to your GitHub account

		gedit ~/.ssh/id_rsa.pub
	
	In the top right corner of any GitHub page in your browser, click on your profile photo, the click 'Settings'
	In the user settings sidebar, click 'SSH keys'
	Click 'New SSH key'
	In the "Title" field, add a descriptive label for the new key
	Copy and Paste the contents from the 'id_rsa.pub' file into the "Key" field
	Click 'Add SSH key'
	
4.	Load your keys into your SSH agent
	
		eval "$(ssh-agent -s)"
		
		ssh-add
		
	Enter passphrase
	
5.	Test your SSH connection

		ssh -T git@github.com
		
	You should see the message "Hi 'username'! You've successfully authenicated, but GitHub does not provide shell access."

6.	Clone repository into src folder
	
		git init
		
		git remote add origin PATH/TO/REPO	 #get the PATH/TO/REPO from the SSH clone URL on the repository page
		
		git fetch
		
		git checkout -t origin/master
		
		
When you start work:
git pull

To save changes:
git add -u
git commit -m "Comment about the changes"
git push
