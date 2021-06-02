# BIRDIE
This repository contains code and documentation related to the South Africa Biodiversity Data Pipeline for Wetlands and Waterbirds (BIRDIE) project: https://jrsbiodiversity.org/grants/sanbi-2020/

## Basic structure
There are two basic elements in the BIRDIE repo: i) a package that contains functions and data that is installed as an extension to base R and ii) an analysis folder that contains analysis scripts and data that runs in R but that is not integrated in R. In other words, the package adds functionality and the analysis uses this functionality.

The folders *R*, *data* and *man* contain functions, data and documentation for the BIRDIE package respectively. See the file *main_dev.R* and https://r-pkgs.org/, to get and idea of how to include new functions, data and modify existing ones. The folder *data_prep* stores scripts that are used to create data that is contained in the package, so that we can keep track of how data were produced.

The folder *analysis* contains data, scripts and outputs from analyses related to the BIRDIE project. This folder is not part of the package, so nothing in there is meant to be installed.

Data that is meant to be included in and exported with the package should go to the *data* folder, and it should be properly documented (see *main_dev.R* and https://r-pkgs.org/). These data is usually meant for running examples or it may also be data that is used very frequently. Data that is only meant for running anaysis and not to exported with the package, should go to the *analysis/data* folder. It is not necessary to document these data.

## INSTRUCTION TO INSTALL PACKAGE

1. Clone the repository to your local machine:
   - In RStudio, create a new project
   - In the 'Create project' menu, select 'Version Control'/'Git'
   - Copy the repository URL (click on the 'Code' green button and copy the link)
   - Choose the appropiate directory and 'Create project'
2. Install the package 'devtools' in case you don´t have it and run devtools::install() from the project directory
3. Remember to pull the latest version regularly

## INSTRUCTIONS TO CONTRIBUTE CODE TO THE PACKAGE

For site owners:

There is the danger of multiple people working simultaneously on the project code. If you make changes locally on your computer and, before you push your changes, others push theirs, there might be conflicts. This is because the HEAD pointer in the main branch has moved since you started working. 

To deal with these lurking issues, I would suggest opening and working on a topic branch. This is a just a regular branch that has a short lifespan. In steps:

- Open a branch at your local machine
- Push to the remote repo
- Make your changes in your local machine
- Commit and push to remote
- Merge your changes:
  - In the GitHub repo you will now see an option that notifies of changes in a branch: click compare and pull request.
  - If there are no conflicts 'merge pull request'
- Delete the branch. You will have to delete it in the remote repo (GitHub) and also in your local machine. In your local machine you have to use Git directly, apparently RStudio doesn´t do it:
  - In your local machine, change to master branch.
  - Either use the Git GUI (go to branches/delete/select branch/push).
  - Or use the console typing 'git branch -d your_branch_name'.

Opening branches is quick and easy, so there is no harm in opening multiple branches a day. However, it is important to merge and delete them often to keep things tidy. Git provides functionality to deal with conflicting branches. More about branches here:

https://git-scm.com/book/en/v2/Git-Branching-Branches-in-a-Nutshell

Another idea is to use the 'issues' tab that you find in the project header. There, we can identify issues with the package, assign tasks and warn other contributors that we will be working on the code.
