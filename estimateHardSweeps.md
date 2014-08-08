# A rough step-by-step for evaluating the prevalence of hard sweeps during maize domestication 

##### 8-7-2014

## Analysis I: Maize and teosinte vs. tripsicum

### Step 1
Compute the diversity per site for every bp in the maize (and teosinte)  genome with data. It may be necessary to "correct for systematic variation in the mutation rate" by dividing by the maize-tripsicum/teosinte-tripsicum divergence, similar to what was done by Hernandez et al (2011, pg 920).


### Step 2
Identify potentially functional substitutions (amino acid substitutions) and potentially neutral substitutions (e.g. synonymous, four-fold degenerate?). Do this for both maize and teosinte by comparing to tripsicum.


### Step 3
For every bp in maize and teosinte with data, compute the distance from that bp to the nearest substitution of each class above (ie potentially functional, potentially neutral).


### Step 4
Calculate, for each class above, the average diversity as a function of distance to the nearest substitution of each class. Do this for both maize and teosinte (so each taxa is comparable to a population in the Hernandez (2011) paper).


### Step 5
Interpret the results--If diversity is lower around potentially functional substitutions for both maize and teosinte, this can be interpreted as evidence for wide-spread hard sweeps, at least during the evolution of teosinte. If diversity is not different surrounding potentially functional and potentially neutral substitutions for either maize or teosinte, this can be interpreted as evidence for a lack of wide-spread hard sweeps. If maize and teosinte differ, things get interesting. If maize shows evidence of reduced diversity around potentially functional sites but teosinte doesn't, this suggests that hard sweeps may have played a substantial role in the evolution of maize but not for teosinte. If teosinte shows evidence of reduced diversity around potentially functional sites but maize doesn't, this is hard to interpret... For any of the above outcomes, the next analysis will be useful.


## Analysis II: Maize vs teosinte

### Step 1
Compute the diversity per site for every bp in the maize genome with data. It may be necessary to "correct for systematic variation in the mutation rate" by dividing by the maize-tripsicum  divergence, similar to what was done by Hernandez et al (2011, pg 920). (you have already done this in I.1)


### Step 2
Identify potentially functional substitutions (amino acid substitutions) and potentially neutral substitutions (e.g. synonymous, four-fold degenerate?). Do this for **maize** by comparing to **teosinte**.


### Step 3
For every bp in maize with data, compute the distance from that bp to the nearest substitution of each class above (ie potentially functional, potentially neutral).


### Step 4
Calculate, for each class above, the average diversity as a function of distance to the nearest substitution of each class. Notice that substitutions were computed by comparing **maize** and **teosinte**, so all substitutions are recent.

### Step 5
Interpret results. If there is a reduction in diversity around potentially functional substitutions in maize, this provides evidence that hard sweeps were wide-spread _during_ maize domestication. If there is not a reduction, this provides evidence that hard sweeps were not wide-spread _during_ maize domestication. Coupled with **Analysis I**, conclusions drawn should be strong.
