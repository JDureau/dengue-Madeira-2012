dengue-Madeira-2012
===================

Until last year, dengue had disappeared from the European continent. The last epidemic goes back to 1927-1928, 
in Greece.
However, concerns of a return of dengue in Europe had started to rise in the recent years, due to the dissemination 
of Aedes albopictus across European countries. This mosquito plays a central role in dengue transmission, as it serves as a vector for the 
virus. In reaction to this threat, severeal European Research programs have been launched to improve our understanding of dengue, on several aspects
that are key to predicting and preventing coming epidemics. 


In September 2012, a first epidemic occured in Europe. [2159 cases were recorded over 3 months][1] in the Portughese 
island of Madeira. Among these case, a few individuals were hospitalised for mild symptoms of fever but no severe case
has been recorded. 

![data](https://raw.github.com/JDureau/dengue-Madeira-2012/master/images/data.png?login=JDureau&token=c5b1e3d648591265b128978f10a0bcee)


Although multiple and crucial aspects of dengue transmission are still to be explored, some epidemiologists argue that
severe cases [are more likely to correspond to secondary infections][3]. After having previously been infected with one of
the 4 dengue strains, an individual that is re-infected with another strain would have a much higher probability of developing 
severe symptoms as hemorragic dengue fever. If we follow this assumption, and consider that all infections that occured
in 2012 were primary infections, there is a risk for severe cases in 2013 if any of the primary infected gets re-infected.


We illustrate here how mechanistic models can be used to forecast future epidemics while reflecting the different sources of uncertainty. We have extended a
multi-strain model that had been introduced by [Aguiar et al][2] to study dengue dynamics in South-East Asia. Our model 
accounts for 

* *strain competition*: we know that four strains of dengue coexist
* *cross immunity*: after recovering from a dengue infection, individuals are resistant to all strains for a short period
* *under-reporting*: a high but unknown proportion of individuals infected with dengue do not develop symptoms. They are 
generally not recorded, but they are nevertheless infectious.
* *seasonality*: dengue transmission is highly infuenced by temperature and humidity through the concentration of 
mosquitoes.
* *immigration*: travelers being infected in other parts of the world are regularly identified in Maidera.
* *demographic stochasticity*: disease transmissions are a random process, and this randomness needs to be reflected when 
construcing past and future scenarios of the epidemic. Specially in a small island like Madeira.
* *environmental stochasticty*: all factors playing a role on dengue transmissions, and their variability, may not be
perfectly incorporated in the model. An additional source of stochasticity called environmental is meant to reflect 
this uncertainty.

Under these assumptions, the data from the 2012 epidemic can be used to reconstruct the current state of immunity of 
the population of Madeira, and to project its evolution. A Bayesian approach is followed, to reflect the 
available information on the respective lengths of the infectivity and cross-immunity periods as well as the uncertainty 
on the proportion of asymptomatics and initial state of the population immunity. We nonetheless consider that only less
than 5% of the population had already been infected with dengue before September 2012. Under these assumptions, 
the predicted number of sever dengue cases occurring each week is the following:

![data](https://raw.github.com/JDureau/dengue-Madeira-2012/master/images/forecast.png?login=JDureau&token=e66b78f7f11574ef08f2b064073d0c67)


For the sake of transparency, and to foster further improvements of this preliminary exploration, this repository provides
the means to easily reproduce the presented results. It relies on the [plom-pipe library of inference methods][4]
developed in collaboration with Sébastien Ballesteros as part the the [PLOM.IO project][6].

Reproducing the results:
------------------------

Data is contained in the data folder, in the csv format. Additionally, the folder 2-strains contains json files that 
define a model and link it to the data, following the [PLOM.IO grammar][6]. To generate the code and play with the model
yourself, simply [install the package][7], and compile the model with:

    plom build -t map.json --local

The joint posterior density of paths and parameters can be explored with:

    plom pipe map.json | ./pmcmc --full -M 10000 -J 2000 -a 0.98 -N 8 -n 200
    
From these sampled trajectories, forecasts can be simulated with:

    plom predict map.json -n 303 -X X_0.csv -T trace_0.csv | ./simul sde -o 303 -D 470  --traj 
    


[1]: http://www.ecdc.europa.eu/en/press/news/Lists/News/ECDC_DispForm.aspx?List=32e43ee8-e230-4424-a783-85742124029a&ID=845        "Dengue epidemic in Madeira"
[2]: http://www.epiwork.eu/wp-content/uploads/2010/03/role.pdf "Aguiar et al."
[3]: http://www.ncbi.nlm.nih.gov/pubmed/20639791 "Dengue hemorrhagic fever and shock syndromes."
[4]: https://github.com/plom-io/plom-pipe "plom-pipe"
[6]: http://plom.io/cli/grammar "PLOM.IO grammar"
[7]: http://plom.io/cli "workflow"
