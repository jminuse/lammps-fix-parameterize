A LAMMPS fix for parameterizing Tersoff potentials. Optional Coulombic and Lennard-Jones contributions are available for the range outside the Tersoff cutoffs. 
This fix works for the December 7th, 2015 version of LAMMPs which can be downloaded [here](http://lammps.sandia.gov/tars/lammps-7Dec15.tar.gz).  

# Getting started
	
1.	If you don't have an SSH key, generate one like this:

		ssh-keygen -t rsa -b 4096 -C "your_email@example.com"  #Creates an ssh key, using your GitHub e-mail as a label
		
	When prompted to "Enter a file in which to save the key," press Enter
	If asked to Overwrite, enter 'y'
	At the prompt, type a secure passphrase
	Retype your secure passphrase

	- Add the SSH key to your GitHub account  

		`gedit ~/.ssh/id_rsa.pub`
	
	In the top right corner of any GitHub page in your browser, click on your profile photo, the click 'Settings'
	In the user settings sidebar, click 'SSH keys'
	Click 'New SSH key'
	In the "Title" field, add a descriptive label for the new key
	Copy and Paste the contents from the 'id_rsa.pub' file into the "Key" field
	Click 'Add SSH key'

	- Load your keys into your SSH agent  
	
		```
		eval "$(ssh-agent -s)"
		
		ssh-add
		```
		
	Enter passphrase
	
	- Test your SSH connection  

		`ssh -T git@github.com`
		
	You should see the message "Hi 'username'! You've successfully authenicated, but GitHub does not provide shell access."

2.	From the lammps source folder, download the repository from github:
	
		git init
		
		git remote add origin git@github.com:jminuse/lammps-min-params.git
		
		git fetch
		
		git checkout -t origin/master

	"git checkout -t origin/master" may give the error "Untracked working tree file 'FILE_NAME_HERE' would be overwritten by merge."

	In this case, the simplest option is to delete the offending file, then repeat

		rm FILE_NAME_HERE

		git checkout -t origin/master 

	until "Untracked working tree file" is replaced by a message like "Branch master set up to track remote branch master from origin."

3.	Make the new version of lammps by running
	
		make serial

4. 	You're done! Now each time you sit down to work, just update your local copy via:

		git pull

	And when you save your changes:

		git add -u

		git commit -m "Comment about the changes"

		git push

	This will make your changes appear to others when they use "git pull" later.

# To parameterize a force field

1.	Build example structures for yoru system in Avogadro (including any OPLS bonds)

2.	Optimize with orca (e.g. Opt B98-D3 def2-TZVP)

3.	Analyze at higher level (e.g. SP RI-B2PLYP D3BJ def2-TZVP)

4.	Put molecule file system.cml in orca output directory

5.	Run `python run_min_params.py RUN_NAME`. The output parameters will start to appear in lammps/RUN_NAME_best.tersoff.

7.	When error is reasonable (e.g. below 10%) `cp lammps/RUN_NAME_best.tersoff lammps/md_input.tersoff`. This allows you to run an annealing job with these parameters using `python run_normal.py RUN_NAME_2`

9.	See what went wrong and adjust parameter guess (lammps/input.tersoff), bounds, or data set. Frames from `view lammps/RUN_NAME_2` can be used as example structures. 
