# Genetic Early Warning Signals

Good wildlife management is based on timely action. However, changing regulations and putting in place conservation effort takes time. Detecting population with a high risk of colapse in time to apply these effert before the population is doomed would allow  for much more efficent management. However, good predictive modeling of ecological dynmics is neraly imposible given the complexity of the underlying processes and the high level of sotchasticity. This is futher complicated by the fact that the transitions from one state to another (viable population to extinct) in these systems are often non-linear. As conditions deteriorate, a thershold is passed leading to the rapid collapse of thepopulation.


These rapid shift in state triggered by a small change in conditions are termed critical transition ( or catastrophic bifurcations). Catastrophic bifurcations arise in systems with alternative stable states (or, in general, alternative attractors). The intersting aspect of critical transitions, is that they have a distinct mathematical signature that is detectable even before the transition occcurs. This mathematical signature (a symptom of the approaching collapse) can be identifed even without complete understanding of the complexe underlying mechanism leading to the alternative stable states. These mathematical signatures could therefore be used as early warning signals (EWS) of collapse even in system with incomplete  knowledge of the machanistic causing the dynamics.


EWS have been detected in a wide range of complex system going from physics, astonomy, climate to ecology. **expand**


The use of EWS has attracted significant attention in ecology from both theoretical fields (due to it's interesting methematical properties) and from the wildlife conservation fields (due to its possible impact of management efficiency). Its has recently been suggested that basing EWS not only in trends in abundacne but also on phenotypic trends could significantly improve the effectiveness of EWS (Clemments & Ozgul 2016). 

Increase phenotypic monitoring may not be sufficient to properly quantify intraspecific diversity and its impact on population process. Genetic monitoring may also be  important (Mimura et al. 2016). 

While EWS seem to bee extremely promissing tools, their usefullness in nature has not been demonstrated. The aim is to see if we can detect early warning signals of a population collapse based on the genetic diversity data of a population of wild ungulate which underwent an important collapse. This population of bighorn sheep has undergone intense monitoring both in term of demography, phenotype and genetics so if any sort of signal is to be found, it should be possible to do so in this population.

## Getting Started

This R code allows one to use the demographic, phenotypic and genetic information on the population of Ram Mountain, Alberta, to calculate and detect EWS. It was mostly stolen from the supplementary material of (Clemets & Ozgul 2016). The cleaned data file sould be place in the "cache" folder and R working directory should be set to the root folder. Raw data go in the "data" folder. Raw and cleaned data are available on demand only.


## References and further reading

* [Clements & Ozgul 2016](http://www.nature.com/doifinder/10.1038/ncomms10984)
* [Mimura et al. 2016](http://doi.wiley.com/10.1111/eva.12436)
* [Early Warning Signals Toolbox](http://www.early-warning-signals.org/)


## To do !

* Figure out if (and best approach to use) correctionfor multiple test is needed.
    + At the moment, it's all based on multiplying the sd of critical limits by a constant determined by a normal probability distribution  
* More phenotypic signals
* determine theoretical justifications for alternative stable states

## Author
* [Gabriel Pigeon](https://github.com/GabrielPigeon)
