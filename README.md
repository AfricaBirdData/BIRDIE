# BIRDIE
This repository contains code and documentation related to the South Africa 
Biodiversity Data Pipeline for Wetlands and Waterbirds (BIRDIE) 
project: https://jrsbiodiversity.org/grants/sanbi-2020/.

The aim of the BIRDIE project is to develop a data-to-decision pipeline for wetlands
and waterbirds that uses use statistical analyses to compute policy-relevant indicators
using citizen science data. The data pipeline initiates pulling data from different databases
and repositories, data then moves through cleaning and validation steps, and into statistical
analyses and summaries that result in indicators that are useful for decision makers.
These indicators are then post-processed and presented on website. These indicators
are useful to support South Africa’s national and international reporting requirements,
management of wetland sites, and broader interest in wetlands and waterbirds.

This repository contains the code necessary for data acquisition and statistical
modelling, which are conducted in R. Modelling outputs are then exported to a
MySQL database that is hosted on a different system together with the website and
the API services that allow communication between the two. However, here we will not
be concerned with anything beyond the computation of basic indicators in R.

There are two basic elements in the BIRDIE repo: i) a package that contains functions
and data that is installed as an extension to base R and ii) an analysis folder that
contains analysis scripts that run in R. In other words, the package adds functionality
to R and the analysis scripts use this functionality. This means that to use the
pipeline we don't only have to install the package, but we also need to clone the repository
to have access to the scripts.

Inside of the package there are some vignettes that explain how to run the 
different modules of the pipeline. These vignettes are also available as articles
in the website africabirddata.github.io/BIRDIE/.

## INSTRUCTIONS TO CLONE AND INSTALL PACKAGE

1. Clone the repository to your local machine:
   - In RStudio, create a new project
   - In the 'Create project' menu, select 'Version Control'/'Git'
   - Copy the repository URL (click on the 'Code' green button and copy the link)
   - Choose the appropiate directory and 'Create project'
2. Install the package 'devtools' in case you don´t have it and run `devtools::install(build_vignettes = TRUE)` from the project directory
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
