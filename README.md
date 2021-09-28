# Lahontan-MPVA
This repository includes source code for a multi-population viability analysis conducted for 155 Lahontan cutthroat trout populations across the sub-species' range in the Great Basin, USA.  

The methods were published in *Ecology* ([Leasure et al. 2018](https://doi.org/10.1002/ecy.2538)).

## Repository Structure
`./scripts/` includes all of the scripts used to process data and fit statistical models.  

`./scripts/02 MPVA.JAGS.R` is the JAGS code for the Bayesian statistical model.  

`./wd/` is a working directory that contains input and output data. It is included in .gitignore so the data do not appear on GitHub (file sizes are prohibitive).  

## Contributing
This work was funded by the NASA Ecological Forecasting program (grant no NNX14AC91G), U.S. Bureau of Land Management (cooperative agreement no L14AC0009), U.S. Fish and Wildlife Service (grant no F16AC01280), National Fish and Wildlife Foundation (grant no 45345, 50018, and 53844), and partial salary support from Trout Unlimited. Survey data for Lahontan cutthroat trout and nonnative trout were contributed by Nevada Department of Wildlife, Oregon Department of Fish and Wildlife, California Department of Fish and Wildlife, U.S. Fish and Wildlife Service, University of Nevada–Reno, and Trout Unlimited. Tom Hobbs and Mevin Hooten provided important statistical training during their Bayesian modeling workshop at Colorado State University funded by the National
Science Foundation (grant no DEB1145200). The research team who developed this statistical model included Doug Leasure, Seth Wenger, Nate Chelgren, Helen Neville, Dan Dauwalter, Robin Bjork, Kurt Fesenmyer, Jason Dunham, Mary Peacock, Charlie Luce, Abby Lute, and Dan Isaak. Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.

## License
You are free to reuse and redistribute this code under the terms of a [GNU GPLv3.0](https://www.gnu.org/licenses/gpl-3.0.en.html) license.
