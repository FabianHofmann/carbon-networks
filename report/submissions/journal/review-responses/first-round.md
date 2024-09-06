## Reviewer #1 (Remarks to the Author):


@FabianHofmann
> The article presents an interesting analysis of the foreseen H2 and CO2 transport infrastructures.
> Literature is not accurate and there are many references cited without going in their details.
> The study lacks of realistic elements such as leakage issues, location of compression stations in relation with the capture location.


Thank you for your comments. We appreciate your feedback and will address your concerns in the revised manuscript.

We have expanded the literature review section to provide more context for the study. Since the study touches a wide range of topics in the political and technological landscape, we have tried to find a good balance between providing a comprehensive overview and keeping the text concise. (#TODO)

To simulate a more realistic transport system, the following assumptions have been integrated into the energy system model for transport efficiencies and compression requirements. For hydrogen pipelines an average loss of 1.2% per 1000 km, aligning with industry estimates that account for hydrogen's high permeability and small molecular size, which can lead to higher losses compared to other gases​. The compression energy required is assumed to be 1.8% of the transported energy, which is consistent with the significant energy demands for compressing hydrogen due to its low energy density​ [Environmental Defense Fund​](https://www.edf.org/sites/default/files/documents/Pipeline%20Methane%20Leaks%20Report.pdf). For natural gas pipelines, the model assumes a 1% loss rate per 1000 km, accounting for typical leakage and operational losses in well-maintained natural gas infrastructure​ as reported in [Global Energy Monitor](https://www.gem.wiki/Natural_gas_transmission_leakage_rates). The compression requirement is set at 1%, which falls within the common range for natural gas pipeline operations​ [Environmental Defense Fund​](https://www.edf.org/sites/default/files/documents/Pipeline%20Methane%20Leaks%20Report.pdf)​. For CO2 pipelines an average loss rate of 0.8% per 1000 km is reflecting a generally low leakage rate observed in supercritical or dense phase CO2 transport​ [SpringerLink](https://link.springer.com/article/10.1007/s11356-023-27289-3) [MDPI](https://www.mdpi.com/1996-1073/14/15/4601)​. The compression energy needed is set to 250 KWh per 1000km and tonne CO2 reflecting the energy required to compress and maintain CO2 in a supercritical state for efficient transport​.

---
@FabianHofmann
> There are several graphs that make the reader confused rather than explain the amount of work done by the authors.

We have revised the figures to make them more accessible to the reader. However, we would appreciate it if you could provide more specific feedback on which figures you found confusing so that we can address these issues more effectively. Most of the figures are contain a lot of information which we find to be necessary to underline the argumentation. We have tried to make them as clear as possible.

---

@nworbmot
> The results are not understandable in terms of transferring into policies.

We have written a short policy outlook in the conclusion as well as linking in the introduction to the recently-announced industrial carbon management strategy of the European Commission. Our addition to the conclusion reads:

"For policy makers these results show the need for coordinated planning across sectors like carbon, hydrogen and synthetic fuels since they are strongly interacting, as well as across borders because of the localized nature of good renewable resources and sequestration potentials. Deploying multiple networks offers some robustness should unforeseen problems arise with one of them, and our results show that the system can even cope with neither a CO2 nor H2 network for a cost penalty.  Many of the technologies have not yet been deployed at scale before and will need support for financing, regulation and gaining the acceptance of local populations."


---


## Reviewer #2 (Remarks to the Author):
@FabianHofmann

> Overall, a very impressive paper, with a very clear and concise structure, easy to follow with many interesting insights. The research gap and necessity of the paper clearly are demonstrated, the topic is higly relevant and the chosen methodology are well suited for the task.

> There are a couple of aspects I would like to understand better, i.e. aspects in which I think additional information would be beneficial.

> 1: CCS is a challenge for scenario analyses like yours, since it low cost makes it a backstop technology. The model will try to solve many things with CCS if one does not prevent somehow. Your approach is to limit CS. “This is in order to avoid offsetting fossil emissions with carbon removal if they can technologically be avoided in the first place.” I was wondering about this throughout the paper. So, first of all I could be helpful to point out the importance and weight of this setting early in the text. I think the implications of this should be discussed more, because I believe this to be crucial: At your point of operation, your model is able to provide and sequester CO2 for roughly160 €/t if I’m correct. That basically means that all mitigation measures costlier than that are only part of the solution because of the CS limit. It would be interesting to discuss the dual value/shadow price of tis cap, to understand how much economic pressure could be on relying more on CCS. The fossil lock-in your CS cap prevents could be an economically very attractive. So, if compensating fossil fuel usage is cheaper than the chosen path of FT and other PtX options, I think this would deserve a discussion in the paper. Especially because the compensation would not need to take place in Europe but could be done in regions with even lower DACS/BECCS costs.


Thank you for your insightful comment. In response, we have added an explanation of for our choice to the introduction and have added a new subchapter in the Appendix "Higher Sequestration Limits" which we refer to in the main text. The added chapter discusses new model results for higher seqestration limits and their implications.

Our analysis shows that increasing the limit from 200 Mt/year to 800 Mt/year reduces system costs quite significantly, by up to 9.1%. This reduction is primarily due to a shift towards greater reliance on fossil fuels with CCS, which in turn reduces the need for CU and renewable investments by about a third. At the same time, the hydrogen network loses importance whereas the carbon network quadruples its capacity. Point-sources of carbon in the inland where fossil fuels are used, like gas for heating or naphta for industry, are nearly all captured and tranpsorted to be sequesterd.


Another interesting finding is the drop of the shadow price from €225/t to €51/t when increasing the sequestration limit. Additional sequestration capacity provides economic benefits, these benefits, however, diminish as the limit increases​.


However, this also introduces the risk of reinforcing fossil fuel dependence, potentially undermining long-term decarbonization efforts. The discussion in the appendix now highlights these trade-offs, emphasizing the importance of carefully balancing sequestration limits to avoid economic over-reliance on CCS​.



---

@chrstphtrs
> 2. The H2 provision in the Baseline case confuses me a bit, namely that it works surprisingly well. How is a constant supply possible in regions without geological storage potential? Does this mean that in hours of low/no RES power plants in the regions with storages produce electricity, which is then brought to H2 sinks and turned into H2 again? In our scenarios this led to much bigger costs and almost impossible necessities for electricity network expansions. Is this because your H2 demand is relatively low, and you transform it into other carriers? Or do you have some other form of storage option at the consumption sites. Overall, I would just like to understand better how you supply H2 steadily without a network.

Thank you for the comment, this is indeed a very interesting question.

The effect you describe would indeed be possible in our model, although the model avoids it as much as possible. Finally, the geological potential for hydrogen storage is lower in some regions, but never zero. So the model always has the ability to store some hydrogen in the ground if needed. The bulk of the H2 demand (FT, methanolisation, methanation) is placed by the model and the spatially fixed hydrogen demand from industry is comparatively low.

Therefore, if there is a situation where a fixed hydrogen demand is critical to be met in a particular region, the model would first consider increasing the local hydrogen storage capacity to meet that hour. If this is not possible, the model would consider increasing the electricity grid capacity to meet the local electrolysis demand. However, the imported electricity would mainly come from non-renewable sources nearby or from renewable sources far away (note that there is no hour when no renewable energy is available in Europe). In the baseline scenario we don't see any hydrogen re-electrification in regions with large hydrogen storage, indicating that the effect you describe is not pronounced in our model.

---
@fneum
> 3: In the introduction: “It is assumed that Europe is self-sufficient in energy and does not import any fuels.” I would like to see the consequences of this setting discussed, especially regarding the large liquid quantities it entails. Just a suggestion, but a variation in which half of the quantities are imported discussed briefliy in the Annex would be very valuable. Much of the infrastructure you model and discuss discuss hinges on this, and importing a major share is clearly an realistic and attractive option in reality.



---
@chrstphtrs
> Your H2 prices seem comparatively low to similar exercises. What interest rate/ WACC did you use? I think that should always be made explicit, as it has huge a impact on many results.

We use technology cost assumptions from the Technology Data database referenced in the Methods section. The database applies a discount rate of 7% to derive annualized investment costs for all major technologies, including hydrogen electrolyzers. We have added this information to the Methods section.

---
@FabianHofmann
> Maybe this would be to much, but while seeing your H2 prices I was very interested in seeing your electricity prices. If a map/maps would be too much, maybe state the ranges somewhere for reference?

We have added a chapter to the Appendix that showes maps for other energy carrier (electricity, heat, gas, biomass) in the hybrid scenario. These show the average production, consumption, transport and prices per carrier. We hope this will provide the information you are looking for.

---
@chrstphtrs
> I think an overview table with differences of the scenario settings in the beginning would be helpful, and maybe with some key results in the end.



---
@FabianHofmann
> The colour coding in Figure 6 could be better, the greens are almost indistinguishable and too many similar browns.


---
@chrstphtrs
> In the end of section 2: “In the network layouts, the slight shift in CU deployment results in a partial dismantling of hydrogen pipelines, especially from the Iberian Peninsula to Central Europe (0.5 bn€/a).” Is “dismantling” really the right word here? The model simply does not build those pipelines in this scenario, and you cannot dismantle what you never built. Or maybe you can, and this is a tricky koan.

You are absolutely correct. We rephrased this passage, and similar ones, to accurately reflect the net-zero and net-negative scenarios as mutually exclusive futures, not as chronological events.

---
@chrstphtrs
> “High-cost Direct Air Capture”: I may be nitpicking, but I was wondering whether this expression is really fair in combination with the CS cap. I think you force the model to pursue many even more expensive option like FT, this is just the most expensive carbon option for sequestration. Again, really just a thought.

You are correct, the relative statement of low and high cost only applies to elementary carbon dioxide needed for sequestration. We amended the statement to make it more clear.
