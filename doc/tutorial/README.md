---
author: Remco Bouckaert
level: Intermediate
title: Transmission Tree Tutorial
subtitle: Inferring who-infected-who
beastversion: 2.7.5
---


# Background

----

# Programs used in this Exercise

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees {% cite Bouckaert2014,bouckaert2019beast --file master-refs.bib %}. This tutorial uses the BEAST2 version 2.7.5.

### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### Tracer

Tracer is used to summarise the posterior estimates of the various parameters sampled by the Markov Chain. This program can be used for visual inspection and to assess convergence. It helps to quickly view median estimates and 95% highest posterior density intervals of the parameters, and calculates the effective sample sizes (ESS) of parameters. It can also be used to investigate potential parameter correlations. We will be using Tracer v1.7.0.

----

# Practical: Transmission Tree Analysis

We will set up an analysis in BEAUti using a transmission tree prior of a tuberculosis outbreak in Hamburg, Germany, earlier analysed in {% cite roetzer2013whole %}.
To reduce run-time, we only analyse a subset of 40 samples, then run BEAST and analyse the results. 
We will be using the `transmission` package, so make sure it is installed, like so:

> * Start BEAUti
> * Click to the `File => Manage packages` menu item.
> * Select `transmission` in the list of packages and the click `Install` button.
> * Close BEAUti -- it needs to restart to pick up the new packages.

## Set up in BEAUti

> * Start BEAUti -- it should open with the Standard template.
> * Next, select the `File => Import Alignment` menu.

A dialog is shown where you can select a file containing a tree in NEXUS format.

> Select the file `roetzer40.tree` that comes with this tutorial in the data section.

<figure>
	<a id="fig:BEAUti1"></a>
	<img style="width:30%;" src="figures/BEAUti-import.png" alt="">
	<img style="width:65%;" src="figures/BEAUti-partitions.png" alt="">
	<figcaption>Figure: Add partition through the `File => Import Imprt` menu.</figcaption>
</figure>

In the partition panel, a new partition will be added with the name roetzer40. 
Next, we will set up tip dates.

> * Click on the `Tip Dates` tab.
> * Select the `Use tip dates` checkbox.
> * Click the `as dates with format` and select `yyyy-M-dd` from the drop down box.
> * Click the `Auto-configure` button
> * A dialog pops up. We want everything after the first colon, so change the underscore in the first entry to `:` and click the `OK` button.

<figure>
	<a id="fig:BEAUti2"></a>
	<img style="width:40%;" src="figures/BEAUti-dates.png" alt="">
	<img style="width:55%;" src="figures/BEAUti-dates2.png" alt="">
	<figcaption>Figure: Set up tip dates in the tip-dates panel.</figcaption>
</figure>

Now, we can set up the site model. We will use HKY+4G+I in this analysis.

> * Click the `Site Model` panel.
> * Select `HKY` from the drop down box with label `Subst Model`.
> * Change `Gamma category count` to 4.
> * Change `Proportion invariant` to 0.1 and select the `estimate` check box next to it.

The site model panel should look similar to this:

<figure>
	<a id="fig:BEAUti3"></a>
	<img style="width:70%;" src="figures/BEAUti-site-model.png" alt="">
	<figcaption>Figure: Set up site model to HKY+4G+I.</figcaption>
</figure>

We will leave the clock model to a strict clock. Because we use tip dates, the clock rate is estimated by default. 
Next, we set up the transmission tree prior.

> * Click the `Priors` panel.
> * Change the default Yule Model for tree prior to `Transmission`.
> * New priors appear for the block-count, block-start and end, and for transmission tree population size. Click on the triangle next to `Tree.t:roetzer40` to show the parameter so the `Transmission` tree likelihood.
> * Go to the population size prior (at the bottom), open the distribution by clicking the triangle next to the prior, and set the lower bound to 0.1 (this prevents the tree collapsing).

<figure>
	<a id="fig:BEAUti4"></a>
	<img style="width:70%;" src="figures/BEAUti-priors.png" alt="">
	<figcaption>Figure: Select the `Transmission` tree prior, and this is how the priors panel looks like.</figcaption>
</figure>


<figure>
	<a id="fig:BEAUti5"></a>
	<img style="width:45%;" src="figures/BEAUti-priors1.png" alt="">
	<img style="width:45%;" src="figures/BEAUti-priors2.png" alt="">
	<figcaption>Figure: Shows all options of the transmission tree prior.</figcaption>
</figure>


The transmission tree likelihood has the following components:
* samplingHazard: determines the hazard of being sampled. It has a sampling probability (`C` in the priors tab) and a `shape` and `rate` parameter for a Gamma distribution that determine the time of sampling after a host got infected.
* transmissionHazard: determines the hazard of transmitting an infection. It has an average number of transmissions `C` and a `shape` and `rate` parameter for a Gamma distribution that determine the time from infection to time of infecting another host. In general, the average transmission time should be larger than the average sampling time (so shape/rate of transmission should be larger than shape/rate of the sampling hazard).
* endTime: time at which the study finished relative to the latest sample. So, if the units of time is years, and the study stopped collecting samples 3 months after the latest sample, it means the endTime is 1/4 year after the latest sample, and endTime=-0.25.
* deltaStartTime: time at which the study start till root of tree (optional, default: 0). 
* allowTransmissionsAfterSampling: flag to indicate sampling does not affect the probability of onwards transmissions. If false, no onwards transmissions are allowed (not clear how this affects the unknown unknowns though). (optional, default: true)
* includeCoalescent: flag for debugging that includes contribution from coalescent to posterior if true (default: true).
* population model: A constant size population is assumed for the coalescent process inside each host. The prior on the population size can be set.

For more details, we refer to the paper.


Since we don't want the analysis to take too long, we only run for 5 million samples.

> * In the MCMC panel, set the chainLength to 5 million samples.
> * Optionally, you might want to reduce the log frequency of the screen logger to 100000.
> * Safe the file to `roetzer.xml`

<figure>
	<a id="fig:BEAUti6"></a>
	<img style="width:45%;" src="figures/BEAUti-mcmc.png" alt="">
	<figcaption>Figure: MCMC settings.</figcaption>
</figure>

## Run with BEAST

> Run BEAST on `roetzer.xml`

This should not take more than 5 minutes, but if you don't want to wait that long you can use the data in the `precooked_runs` directory that comes with this tutorial.
The longer runs (marked with `long` in their name) are 

## Check convergence

> Run `Tracer`, and make sure all parameters have sufficiently large ESSs

<figure>
	<a id="fig:Tracer"></a>
	<img style="width:65%;" src="figures/tracer.png" alt="">
	<figcaption>Figure: Convergence of MCMC in Tracer.</figcaption>
</figure>

The short 5 million sample run has not quite converged yet, but the longer run does.


Inspect the tree file, for example in DensiTree. It shows that there is a lot of uncertainty in the tree distribution

<figure>
	<a id="fig:DensiTree"></a>
	<img style="width:65%;" src="figures/densitree.png" alt="">
	<figcaption>Figure: Inspect the tree distribution in DensiTree.</figcaption>
</figure>

## Visualising Who-Infected-Who

The `WIWVisualiser` app creates an SVG files that visualises who-infected-who. 
To start the `WIWVisualiser` app, 

> * Select the `File =>> Launch apps` menu in BEAUti.
> * Select `WIWVisualiser` from the list of apps, and click the `Launch` button.
> * In the menu that pops up, select the tree file (or select `roetzer-roetzer40.trees` from the pre-cooked runs)
> * Select an appropriate output file name
> * Set the `Partition` to `roetzer40`
> * Click the `OK` button, and after a little while a message in the terminal appears that the svg file has been written.
> * Close the terminal window.

<figure>
	<a id="fig:WIWVisualiser"></a>
	<img style="width:65%;" src="figures/WIWVisualiser.png" alt="">
	<figcaption>Figure: Who-infected-who visualiser options</figcaption>
</figure>


Alternatively, you can run it from the command line like so:

```
/path/to/applauncher WIWVisualiser -tree roetzer-roetzer40.trees -partition roetzer40 -out /tmp/roetzer.svg
```

You can open the svg file in a web browser for closer inspection, or in a vector drawing program like Inkscape or Illustrator.


<figure>
	<a id="fig:WIWVisualiser2"></a>
	<img style="width:65%;" src="figures/WIWVisualiser2.png" alt="">
	<figcaption>Figure: Part of the who-infected-who network shown in a browser.
The numbers on edges are the probability that the host at the tail infected the host at the head.
Numbers in brackets represent the total probability a host was infected by another sampled host (so one minus that number is the probability it was infected by a host that was not sampled).
	</figcaption>
</figure>


WIWVisualiser has the following options:

* trees (TreeFile): tree file file with transmission trees. (optional, default: [[none]])
* log (LogFile): trace file containing infectorOf log. Ignored if tree file is specified (optional, default: [[none]])
* burnin (Integer): percentage of trees to used as burn-in (and will be ignored) (optional, default: 10)
* out (OutFile): output file, or stdout if not specified (optional, default: /tmp/wiw.svg)
* prefix (String): prefix of infectorOf entry, e.g., infectorOf (optional, default: infectorOf)
* threshold (Double): probability threshold below which edges will be ignored. (optional, default: 0.1)
* partition (String): name of the partition appended to `blockcount, blockend and blockstart` (optional)
* suppressSingleton (Boolean): do not show taxa that are not connected to any otehr taxa (optional, default: true)

## Transmission Tree Statistics

The `TransmissionTreeStats` app provide statistics of set of transmission trees. It has the following options:

* trees (TreeFile): tree file file with transmission trees. (optional, default: [[none]])
* burnin (Integer): percentage of trees to used as burn-in (and will be ignored). NB default 0 (optional, default: 0)
* partition (String): name of the partition appended to `blockcount, blockend and blockstart` (optional)
* out (String): directory where to put files with tranmsision & sampling time stats (optional, default: /tmp)

This is a sample output of `TransmissionTreeStats`:
```
average internal Branch Length = 0.62
average leaf transmission count = 1
average internal transmission count = 1
average transmission count per tree = 170
average unsampled hosts per tree = 85
average maximum block count = 1
```

It also produces two log files in the specified output directory called `timeTillSampling.dat` and `timeTillTransmission.dat`.
The `timeTillSampling.dat` can be opened in `Tracer` and reveals for each sampled host the time from first infection to sampling.
It can be revealing in that for some hosts the distribution is multi-modal.

<figure>
        <a id="fig:TransmissionTreeStats"></a>
        <img style="width:45%;" src="figures/TransmissionTreeStats.png" alt="">
        <img style="width:45%;" src="figures/TransmissionTreeStats2.png" alt="">
        <figcaption>Figure: TransmissionTreeStats GUI version (left) and distributions of time-till-sampling (right) for three selected hosts.</figcaption>
</figure>

The `timeTillTransmission.dat` file contains information about the time it takes between the first infection of a host and the time the host infects another host. 
It also records the time till infects a second host.
If no host(s) are infected, -1 is output -- this shows up when you open the file in Tracer as a peak around -1.
Checking the time-till-first-tranmission and time-till-second-transmission can be useful in verifying the transmission hazard is properly parameterised.


## Transmission Tree Simulation

The `tranmission` package comes with a tranmission tree simulator that can be run as an app via the `applauncher`. It simulates transmission trees with colouring and block counts.

`TransmissionTreeSimulator` has the following options:
* endTime (Function): end time of the study (optional, default: constant)
* popSize (Function): population size governing the coalescent process (optional, default: constant)
* sampleShape (Function): shape parameter of the sampling intensity function (optional, default: constant)
* sampleRate (Function): rate parameter of the sampling intensity function (optional, default: constant)
* sampleConstant (Function): constant multiplier of the sampling intensity function (optional, default: constant)
* transmissionShape (Function): shape parameter of the transmission intensity function (optional, default: constant)
* transmissionRate (Function): rate parameter of the transmission intensity function (optional, default: constant)
* transmissionConstant (Function): constant multiplier of the transmission intensity function (optional, default: constant)
* out (OutFile): output file. Print to stdout if not specified (optional)
* trace (OutFile): trace output file, or stdout if not specified (optional, default: [[none]])
* seed (Long): random number seed used to initialise the random number generator (optional)
* maxAttempts (Integer): maximum number of attempts to generate coalescent sub-trees (optional, default: 1000)
* taxonCount (Integer): generate tree with taxonCount number of taxa. Ignored if negative (optional, default: -1)
* maxTaxonCount (Integer): reject any tree with more than this number of taxa. Ignored if negative (optional, default: -1)
* treeCount (Integer): generate treeCount number of trees (optional, default: 1)
* directOnly (Boolean): consider direct infections only, if false block counts are ignored (optional, default: true)
* quiet (Boolean): suppress some screen output (optional, default: false)


----

# Useful Links

- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file master-refs.bib %}
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users)

----

# Relevant References

{% bibliography --cited --file master-refs.bib %}
