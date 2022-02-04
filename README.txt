README

This dataset is from the paper titled
    "The evolution of marine larval dispersal kernels in spatially structured habitats: analytical models, individual-based simulations, and comparisons with empirical estimates"
    By: Allison K Shaw, Cassidy C D'Aloia, Peter M Buston
    Published in: The American Naturalist

Contact ashaw@umn.edu for assistance.

The following files are included:

(1) Folder (.tar.gz) "files_runsims"
     * Contains Matlab code (.m files) for running model simulation data
       used for Figures 2-4 in the paper.
     * Contents:
	> create_env_bounded.m: creates bounded, homogeneous environment

	> create_env_heterogeneous.m: creates bounded, heterogeneous environment

	> create_env_reef.m: creates seascape environment (homogeneous or
           heterogeneous)

	> create_env_unbounded.m: creates unbounded, homogenous environment

        > IBM_dispersal_2D.m: runs a simulation of the individual-based model

        > IBM_dispersal_2D_wK.m: runs alternative version of simulation of the
           individual-based model with more than 1 individual allowed per site

        > main_runcode_step1.m: creates the environments

        > main_runcode_step2.m: runs IBM simulations in each environment

        > main_runcode_multpersite.m: runs alternative IBM simulations in
           seascape environment allowing 2 and 4 individuals per site

        > trimmed_ascii.txt: data file of 0s and 1s of the seascape environment,
           where unsuitable habitat is '0' and forereef is '1'.


(2) Folder (.tar.gz) "output_environments"
     * Contains Matlab data files (.mat files) for each type of simulation
         environment.
     * Contents:
	> env_bounded_sx*.mat (3 types of bounded, homogenous environment)

	> env_bounded_het_sx*.mat (3 types of bounded, heterogeneous environment)

	> env_reef_*.mat (4 types of seascape environment)

	> env_unbounded_sx*.mat (1 type of unbounded, homogenous environment)


(3) Folder (.tar.gz) "output_simulations"
     * Contains Matlab data files (.mat files) for each individual-based
         simulation.
     * Contents:
 	> IBM_unbounded_sx*.mat (30 IBM output for different survival values in
           unbounded, homogenous env)

	> IBM_bounded_het_sx*.mat (30 IBM output for different survival values in
           bounded, heterogeneous env)

	> IBM_bounded_sx*.mat (30 IBM output for different survival values in
           bounded, homogenous env)

	> IBM_reef_*.mat (6 IBM output in seascape environment)


(4) Folder (.tar.gz) "files_plotsims"
     * Contains Matlab code (.m files) for plotting Figures 2-4 in the paper,
         using the above .mat files.
     * Contents:
	> plot_fig2.m: generates figure 2

	> plot_fig3.m: generates figure 3

	> plot_fig4.m: generates figure 4


(5) Folder (.tar.gz) "files_plotcomparisons"
     * Contains Matlab code (.m files), R code (.R files) and data files (.csv)
          to compare simulated and empirical data and to plot Figure 5 in the
          paper.

	> get_counts.m: opens the IBM output for the reef simulations and counts
           the number of individuals allocated to each dispersal bin, for
           comparison with empirical data.

        > plot_fig5.r: R code to compare the empirical and simulated dispersal
           data, and to produce figures 5 and E.3. 

        > dispersal_frequencies.csv: The relative frequencies of dispersal
           distances for E. lori empirical data, and various simulations. The
           simulations include four different scenarios of resource
           heterogeneity ("simulated_high; "simulated_med"; "simulated_low";
           "simulated_none"); and three different carrying capacities per site
           ("capacity_1"; "capacity_2"; "capacity_3"). Data are formatted for
           plotting with ggplot2.

        > distances.csv: Observed dispersal distances for E. lori empirical data
           ("dist_empirical") and simulated data under four scenarios of
           resource heterogeneity (none: "dist_no"; low: "dist_low"; medium:
           "dist_med"; and high: "dist_high"). 


