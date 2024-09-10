
# Response to Reviewer 2


We thank both reviewers for their helpful comments. We have rerun all scenarios in the paper with improved assumptions as suggested by the authors (improved accounting for leakage and compression losses) as well as new scenarios requested by the reviewers for a sensitivity to sequestration rates.


Please let us reply to the reviewers' comments below.

Comment:

> Overall, a very impressive paper, with a very clear and concise structure, easy to follow with many interesting insights. The research gap and necessity of the paper clearly are demonstrated, the topic is higly relevant and the chosen methodology are well suited for the task.

> There are a couple of aspects I would like to understand better, i.e. aspects in which I think additional information would be beneficial.

> 1: CCS is a challenge for scenario analyses like yours, since it low cost makes it a backstop technology. The model will try to solve many things with CCS if one does not prevent somehow. Your approach is to limit CS. “This is in order to avoid offsetting fossil emissions with carbon removal if they can technologically be avoided in the first place.” I was wondering about this throughout the paper. So, first of all I could be helpful to point out the importance and weight of this setting early in the text. I think the implications of this should be discussed more, because I believe this to be crucial: At your point of operation, your model is able to provide and sequester CO2 for roughly160 €/t if I’m correct. That basically means that all mitigation measures costlier than that are only part of the solution because of the CS limit. It would be interesting to discuss the dual value/shadow price of tis cap, to understand how much economic pressure could be on relying more on CCS. The fossil lock-in your CS cap prevents could be an economically very attractive. So, if compensating fossil fuel usage is cheaper than the chosen path of FT and other PtX options, I think this would deserve a discussion in the paper. Especially because the compensation would not need to take place in Europe but could be done in regions with even lower DACS/BECCS costs.


Thank you for your insightful comment. In response, we have added an explanation for our choice to the introduction and have added new modelling runs in a new subchapter in the Appendix "Higher Sequestration Limits" which we refer to in the main text. The added chapter discusses new model results for higher seqestration limits and their implications.

Our analysis shows that increasing the limit from 200 Mt/year to 800 Mt/year reduces system costs quite significantly, by up to 9.1%. This reduction is primarily due to a shift towards greater reliance on fossil fuels either with CCS or unabated and compensated with CDR, which in turn reduces the need for CU and renewable investments by about a third. At the same time, the hydrogen network loses importance whereas the carbon network quadruples its capacity. Point-sources of carbon in the inland where fossil fuels are used, like gas for heating or naphtha for industry, are nearly all captured and tranpsorted to be sequesterd.


Another interesting finding is the drop of the shadow price from €225/t to €51/t when increasing the sequestration limit. Additional sequestration capacity provides economic benefits, these benefits, however, diminish as the limit increases​.


However, this also introduces the risk of reinforcing fossil fuel dependence, potentially undermining long-term decarbonization efforts. The discussion in the appendix now highlights these trade-offs, emphasizing the importance of carefully balancing sequestration limits to avoid economic over-reliance on CCS​.



Comment:

> 2. The H2 provision in the Baseline case confuses me a bit, namely that it works surprisingly well. How is a constant supply possible in regions without geological storage potential? Does this mean that in hours of low/no RES power plants in the regions with storages produce electricity, which is then brought to H2 sinks and turned into H2 again? In our scenarios this led to much bigger costs and almost impossible necessities for electricity network expansions. Is this because your H2 demand is relatively low, and you transform it into other carriers? Or do you have some other form of storage option at the consumption sites. Overall, I would just like to understand better how you supply H2 steadily without a network.

Thank you for the comment, this is indeed a very interesting question.

The spatially-fixed hydrogen demand from industry (primarily for primary steel but also for ammonia) is comparatively low, while the bulk of the H2 demand (FT, methanolisation, methanation) can be moved by the model to locations with low-cost renewables and geological H2 storage. This avoids over-burdening the electricity grid.

In addition, to provide hydrogen in some regions with poor renewable resources the model resorts to SMR with either green or fossil methane transported in the gas grid.

Comment:

> 3: In the introduction: “It is assumed that Europe is self-sufficient in energy and does not import any fuels.” I would like to see the consequences of this setting discussed, especially regarding the large liquid quantities it entails. Just a suggestion, but a variation in which half of the quantities are imported discussed briefliy in the Annex would be very valuable. Much of the infrastructure you model and discuss discuss hinges on this, and importing a major share is clearly an realistic and attractive option in reality.


This is indeed an important point that we have not discussed in detail in the paper. The assumption of Europe's self-sufficiency in energy was made to simplify the system design and focus on internal infrastructure. It also has the advantage of making good use of Europe's sustainable biomass sources for sourcing CO2, and avoiding the difficulty of regulating sustainable CO2 provision for carbonaceous fuels synthesised outside of Europe. However, we agree that importing a significant share of energy carriers, particularly liquid fuels, is a realistic and potentially cost-beneficial option.

We have expanded the discussion of the limitations in the appendix to explore the potential consequences of allowing imports. In the higher sequestration scenario, for example, our model shows that increased imports of fossil oil significantly reduce the need for domestic Fischer-Tropsch (FT) fuel production. While this reduces hydrogen transport and carbon utilization requirements in Central Europe, it also necessitates substantial carbon sequestration within Europe to offset emissions from fossil-based imports. If these imports were carbon-neutral, the need for dedicated sequestration infrastructure would be significantly reduced, as the emissions would already be compensated.
Thus, the carbon neutrality of imported fuels has a profound effect on the optimal infrastructure layout, especially for hydrogen transport and carbon sequestration. A scenario that allows for partial energy imports could lead to lower investments in the carbon and hydrogen networks as well as direct air capture.

However, a critical political and strategic consideration is the control over the CO2 supply chain. There may not be sufficient sustainable CO2 available abroad, and Europe may prefer to maintain control of the CO2 supply chain to ensure security and reliability. Relying on external sources for carbon-neutral fuels could introduce geopolitical risks, similar to those faced with fossil fuel imports.


Comment:

> Your H2 prices seem comparatively low to similar exercises. What interest rate/ WACC did you use? I think that should always be made explicit, as it has huge a impact on many results.

We use technology cost assumptions from the Danish Energy Agency referenced in the Methods section. The database applies a discount rate of 7% to derive annualized investment costs for all major technologies, including hydrogen electrolyzers. We have added this information to the Methods section.

Comment:

> Maybe this would be to much, but while seeing your H2 prices I was very interested in seeing your electricity prices. If a map/maps would be too much, maybe state the ranges somewhere for reference?

We have added a chapter to the Appendix that showes maps for other energy carrier (electricity, heat, gas, biomass) in the hybrid scenario. These show the average production, consumption, transport and prices per carrier. We hope this will provide the information you are looking for.

Comment:

> I think an overview table with differences of the scenario settings in the beginning would be helpful, and maybe with some key results in the end.

We have added an overview table in the introduction that shows the main differences between the scenarios. Since the difference of the scenarios is important on the qualitative level, we have not added a table with key results in the end. We hope that the reader can follow the argumentation and the results are clear enough.

Comment:

> The colour coding in Figure 6 could be better, the greens are almost indistinguishable and too many similar browns.


Thank you for your feedback. We have revised the color scheme in Figure 6 to make it more accessible to the reader. We hope that the new color scheme will make the figure easier to read.


Comment:

> In the end of section 2: “In the network layouts, the slight shift in CU deployment results in a partial dismantling of hydrogen pipelines, especially from the Iberian Peninsula to Central Europe (0.5 bn€/a).” Is “dismantling” really the right word here? The model simply does not build those pipelines in this scenario, and you cannot dismantle what you never built. Or maybe you can, and this is a tricky koan.

You are absolutely correct. We rephrased this passage, and similar ones, that the the net-zero and net-negative scenarios are modelled as different futures, not as chronological events.



Comment:

> “High-cost Direct Air Capture”: I may be nitpicking, but I was wondering whether this expression is really fair in combination with the CS cap. I think you force the model to pursue many even more expensive option like FT, this is just the most expensive carbon option for sequestration. Again, really just a thought.

You are correct, the relative statement of low and high cost only applies to elementary carbon dioxide needed for sequestration. We amended the statement to make it more clear.
