# TSO-DSO Coordination for re-dispatch with flexible DERs

Pearson, Simon; Wellnitz, Sonja; Crespo del Granado, Pedro; Hashemipour, Naser

The value of TSO-DSO coordination in re-dispatch with flexible Decentralized Energy Sources: Insights for Germany in 2030

For the open access publication, we provide our model code based on JuliaLang (v1.5.4) under the MIT licence.

# Abstract

As Germany plans to raise the share of energy consumption satisfied by \ac{res} up to 65\%, congestion in the transmission grids will drastically increase in the power system unless the grids are substantially upgraded or new flexibility options are considered. In this paper, we explore possible integration mechanisms of a potentially powerful source of flexibility: \acfp{der}. Currently, there is no advanced regulatory system for including them into the established electricity markets. We investigate the application of load flexibility DERs can provide for assisting the re-dispatch necessary in electricity markets that employ a zonal pricing mechanism. We implement two different cases with varying levels of involvement of the DSOs and compare their performance with a business-as-usual case in one scenario from 2015 and one prediction for 2030. Findings include that while both cases facilitate the system-wide re-dispatch concerning volume and cost, the average value of optimal load-shifting is not high enough in 2015 to incentivize investment in this area. However, at the higher percentages of generation from \ac{res} in the future scenario, this value becomes promising and using \acp{der} for this purpose may provide long-term benefits to the system operators and owners of assets alike.

# Links

- This model is based on the source code published unter the MIT licence by Xiong et al. (2020) for their paper "Spatial flexibility in redispatch: Supporting low carbon energy systems with Power-to-Gas" (https://doi.org/10.1016/j.apenergy.2020.116201).
- The technology class definitions are based on [Joulia.jl](https://github.com/JuliaEnergy/Joulia.jl/) by J. Weibezahn and M. Kendziorski.
- Like Xiong et al., we use the open source electricity system data set [ELMOD](https://ideas.repec.org/p/diw/diwddc/dd83.html) for Germany (2015).
- The projected 2030 electricity system data set is [The German Electricity System in 2030: Data on Consumption, Generation, and the Grid](https://bwdatadiss.kit.edu/dataset/254) by vom Scheidt et al. (2020).