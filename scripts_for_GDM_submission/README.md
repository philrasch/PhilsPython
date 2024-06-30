
This directory houses the data and scripts for submitted version of the paper titled
"A protocol for model intercomparison of impacts of Marine Cloud
Brightening Climate Intervention" by
Philip J. Rasch, Haruki Hirasawa, Mingxuan Wu, Sarah J. Doherty, Robert Wood, Hailong Wang,
Andy Jones, James Haywood, and Hansi Singh

For questions contact [Phil](mailto:philip.j.rasch@gmail.com) or [Haruki](mailto:harukih@uw.edu)

Examples from simulations with three models are described in the paper (CESM2, E3SMv2, and UKESM1). The data used for figures in the paper is stored in subdirectories with the same names. The python scripts to manipulate and plot that data are stored in a folder called analysis_scripts. Python scripts created by Haruki and Phil are stored in the form of jupyter notebooks with "\*.ipynb" file extensions. Haruki named his scripts according the figure numbers appearing in the paper (for example "Rasch_et_al_Figure_2_3_4.ipynb"). Phil named his scripts with names describing the figure content (for example "gavg_timeseries.ipynb" for figures plotting timeseries of global averages of simulations). There is also a file called "pjrlib.py" containing functions that Phil frequently uses in those scripts. It is included when needed.

The organizational structure of the data won't correspond exactly with the organization of the data being distributed to the outside world, since we have other data (cases, variables, ensemble members, time periods) that might be used in future analyses so the organization is quite elaborate, but we are trying to keep things simple for people wanting just the data used for the paper. So we have simplified the structure to minimize work for others.

But the filenames of the simple organizational structure, and the complex structure we used are unique, and descriptive. So if you just focus on the filenames, and rewrite the directories used in the scripts, they will work fine and the ipynb files will show you the expected results.

