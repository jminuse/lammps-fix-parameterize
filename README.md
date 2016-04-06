A LAMMPS fix for parameterizing Tersoff potentials. Optional Coulombic and Lennard-Jones contributions are available for the range outside the Tersoff cutoffs. 

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
		
		git remote add origin git@github.com:jminuse/lammps-min-params.git
		
		git fetch
		
		git checkout -t origin/master

	"git checkout -t origin/master" may give the error "Untracked working tree file 'FILE_NAME_HERE' would be overwritten by merge."

	In this case, the simplest option is to delete the offending file, then repeat

		rm FILE_NAME_HERE

		git checkout -t origin/master 

	until "Untracked working tree file" is replaced by a message like "Branch master set up to track remote branch master from origin."

7. You're done! Now each time you sit down to work, just update your local copy via:

		git pull

	And when you save your changes:

		git add -u

		git commit -m "Comment about the changes"

		git push

	This will make your changes appear to others when they use "git pull" later.

