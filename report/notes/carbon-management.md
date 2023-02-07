# Carbon Management

## Scenarios

primary:

- co2 network y/n
- H2 network y/n

secondary:

- vary allowed co2 sequestration between 200 Mt and 1000 Mt
- with/without onshore sequestration potentials (default only offshore)
- vary recycling / circular economy potentials (plastics, steel, aluminum)
- myopic pathway to observe CO2/H2 pipeline development
- maybe: vary allowed biomass potentials (Low, Med, High)
- maybe: net-zero and negative emissions scenarios? (leave to Chalmers)

setup:

- 181 nodes, 25-hourly/4-hourly, without gas network, with biomass transport

## Graphics

- map with CO2 network, sequestration as choropleth, sources as pie charts (Fabian H)
- carbon sankey diagrams (Fabian N)
- carbon balance, origin of carbon emissions
- area chart with costs as function of sequestration potentials
- mismatch of biomass potentials and biomass usage

## Abstract

- necessary to establish an appropriate carbon management strategy
- the required carbon (especially from the chemical industry) must be managed in a cycle in order to avoid emissions
- implement the required production and demand processes of the chemical industry, the CO2 sources in question and the necessary logistics (transport and distribution)
- industry as the "carbon manager"
- closing the loop for carbon cycles
- more recycling, circular economy
- keep carbon stock in chemicals sector
- "carbon backbone" and detailed carbon management required
- CO2 infrastructure is "final component" for climate neutrality
- "Carbon management strategy"
- under what circumstances there is co2 transport?
- cement process emissions widely distributed
- truck or ship rather than pipeline? see below!
- or local methanation instead? (deliver H2/elec to sites with process emissions)
- recycling of plastics to strengthen resource efficiency
- many studies do not account for non-energy feedstocks
- and they assume stable domestic basic materials production
- projections for rising steel (2x), platics (3x), cement (2x) demand
- threat of industry relocation -> "green leakage" in analogy to "carbon leakage"
- BECCS, DACCS and green polymers offset unavoidable residual emissions
- climate change mitigation requires a lot of new infrastructure (e.g. networks)
- chicken-egg problem between supply, demand, infrastructure in industry (Open Grid Europe, OGE)
- "carbon backbone" and detailed carbon management required
- CO2 infrastructure is "final component" for climate neutrality
- OGE starts with CO2 pipeline infrastructure in North West Germany 2025/2030 with TES
- OGE calls for CO2 strategy in analogy to H2 strategy
- It is criticised that CO2 infrastructure is not the highest priority in policymaking
- industry transition also requires a system perspective (Philipp Hauser, Agora EW)
- DG Climate says all infrastructure will need to be cross-border (also CO2)
- CO2 as raw material (use in produce synthetic fuels and feedstocks)
- CO2 sources from industry (e.g. calcination of cement)
- scenario: carbon sequestration only in North Sea

## Model Development (not for IEW!)

roughly in order of priority

- revised CO2 pipeline connection options (currently where gas/electricity network)
- integrate CO2 geological sequestration potential properly (e.g. configuration options, upstream processing scripts, converting absolute potentials to annual values, https://github.com/PyPSA/pypsa-eur-sec/issues/156)
- oxyfuel option in cement (maybe: https://publications.jrc.ec.europa.eu/repository/handle/JRC131246)
- add Allam cycle plants
- hydrogen turbines
- allow local intermediate above-ground storage storage of CO2 (needs separate co2 sequestration buses and a unidirectional link such that sequestered co2 cannot be extracted)
- endogenise choices for industrial process heat, then don't need to allow biomass to be transported for feasibility!
- Add waste-to-energy (WtE) plants w/wo CC to consume non-recycled plastics (https://github.com/PyPSA/pypsa-eur-sec/issues/173)
- Allow process emissions from mechanical and chemical recycling of plastics to be captured (https://github.com/PyPSA/pypsa-eur-sec/issues/170)
- non-energy fossil supply outside of the chemical industry (https://github.com/PyPSA/pypsa-eur-sec/issues/169)
- carbon dioxide consumption for urea, food and drink (https://github.com/PyPSA/pypsa-eur-sec/issues/168)
- Fix the electricity and heat consumption of carbon capture in industry (https://github.com/PyPSA/pypsa-eur-sec/issues/82)
- add import options, in particular TES (carbon cycling)
- capture CO2 from biogas upgrading (https://github.com/PyPSA/pypsa-eur-sec/pull/291, https://github.com/PyPSA/pypsa-eur-sec/issues/49)
- https://github.com/PyPSA/pypsa-eur-sec/issues/272

## Notes on CO2 transport options

Why transport?
- separated geological storage potential and CO2 utilisation sites
- modes: pipeline, ship, road
- for large volumes and long distances only pipeline or ship are viable options

CO2 properties
- compress CO2 to high-pressure gas/fluid (80-160 bar) or liquefaction
- density 800-1000 kg/m³

CO2 leakage (emissions per transported volume):
- pipeline 0.05%
- ship 0.4%
- truck 1.6%

By ship:
- more flexible as routes can be changed/upgraded
- but requires costly CO2 terminals
- for medium to long transport and large amounts
- not in operation to date

By pipeline:
- favoured for transport distances up to 500-700km
- CO2 pielines exist in USA, Canada, Norway, Turkey for enhanced oil recovery (EOR) and in the Netherlands from gas processing plants to greenhouses
- diameters 4" to 20"

Gas pipeline retrofitting:
- not suitable for dense phase CO2 transport because max operating pressure of gas transmission lines is 80 bars, 40 bars for main distribution lines

## Links

https://www.catf.us/2023/01/germany-makes-significant-progress-in-advancing-carbon-capture/

https://www.catf.us/2022/04/why-germany-needs-carbon-management-strategy/

https://langfristszenarien.de/enertile-explorer-wAssets/docs/LFSIII_Webinar16.11.2022_Industrie_final.pdf (page 43)

https://ens.dk/sites/ens.dk/files/Analyser/technology_data_for_energy_transport.pdf

https://www.wirtschaft.nrw/carbon-management-strategie-nrw

https://background.tagesspiegel.de/energie-klima/deutschland-braucht-eine-carbon-management-strategie

https://cdn.catf.us/wp-content/uploads/2022/07/07113504/Fu%CC%88nf-Kernkomponente-einer-nationalen-Carbon-Management-Strategie-in-Deutschland.pdf

https://co2-netz.de/en

https://oge.net/en/press-releases/2022/oge-and-tes-join-forces-to-develop-a-1-000-km-co-2-transmission-system

https://www.bveg.de/wp-content/uploads/2022/12/20221220_Positionspapier_Carbon-Management-Strategie.pdf

https://www.agora-energiewende.de/veroeffentlichungen/klimaneutrale-industrie-hauptstudie/

https://network.bellona.org/content/uploads/sites/3/2022/12/Letter-on-CCS-Strategy-1-1.pdf

https://wintershalldea.com/en/newsroom/wintershall-dea-and-equinor-partner-large-scale-ccs-value-chain-north-sea
