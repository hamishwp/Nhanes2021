---
title: "NHANES Blood Pressure-Based Mortality Risk - Appendix"
author: "Rscripts by Hamish Patten, DW Bester and David Steinsaltz"
date: "03/08/2024"
output: 
  bookdown::pdf_document2:
    keep_tex: true
    toc: true
    toc_depth: 3
    number_sections: true
    table_caption: true
    latex_engine: xelatex  # Change the LaTeX engine to xelatex
header-includes:
  - "\\usepackage{lscape}"  # Use lscape package to change page orientation
bibliography: nhanesBP.bib
---


# Appendix A -- The data







\begin{figure}

{\centering \includegraphics{NHANES-BP-Appendix_files/figure-latex/Venn1-1} 

}

\caption{Venn diagram of subjects excluded from the analysis.}(\#fig:Venn1)
\end{figure}

## Exclusions

There were 19592 subjects in the initial data set.
Of these 4573 were excluded because they had missing data or were not followed up, or belonged to the "Other" ethnic group.
This left 15019 subjects for further consideration.
A small number of subjects were excluded because their blood pressure measurements were outside the normal range, as described below in section \ref{sec:BPrange}.
As our method depends on estimating the mortality rates for each demographic group (ethnicity and sex),
we removed the small number of subjects whose ethnic group was given as "Other" (n=751).
(The three included ethnic groups were Mexican American (n=5150), Black (n=5336), and White (n=8355).
In the end there were 14654 subjects in the analysis data set.
A Venn diagram of the different causes of exclusion is given in figure \ref{fig:Venn1}.
We will refer to this as the "full population".
Of these, 9008 had a computable FRS score.
We call this the ``FRS population''.




## Exploratory data analysis
The empirical means of the home and clinic measures in population B are tabulated in Table \ref{tab:summaries}. We note that the home measures are systematically higher than the clinic measures, within every demographic group, with greater differences for subjects who are white or Mexican, and female. 
The average difference is about 2.2 for diastolic and 2.7 for systolic, which is small compared with the general range of the differences, which have SD of 
10.5 (diastolic) and 14.8 (systolic).

\begin{table}[!h]
\centering
\caption{(\#tab:summaries)Summary data for blood pressure}
\centering
\begin{tabular}[t]{llllrr}
\toprule
Place & Sys/Dias & Sex & Ethnicity & Mean & SD\\
\midrule
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{128.8} & \cellcolor{gray!10}{2.7}\\
Home & Systolic & Male & White & 131.2 & 2.8\\
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{127.1} & \cellcolor{gray!10}{2.5}\\
Home & Systolic & Female & Black & 123.2 & 2.7\\
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{127.2} & \cellcolor{gray!10}{2.9}\\
Home & Systolic & Female & Mexican & 121.2 & 2.7\\
\addlinespace
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Diastolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{78.7} & \cellcolor{gray!10}{2.4}\\
Home & Diastolic & Male & White & 77.1 & 2.3\\
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Diastolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{76.8} & \cellcolor{gray!10}{2.4}\\
Home & Diastolic & Female & Black & 75.0 & 2.4\\
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Diastolic} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{73.7} & \cellcolor{gray!10}{2.3}\\
Home & Diastolic & Female & Mexican & 72.0 & 2.4\\
\addlinespace
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{127.7} & \cellcolor{gray!10}{3.5}\\
Clinic & Systolic & Male & White & 128.0 & 4.2\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{123.5} & \cellcolor{gray!10}{3.5}\\
Clinic & Systolic & Female & Black & 122.4 & 3.6\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{123.5} & \cellcolor{gray!10}{4.1}\\
Clinic & Systolic & Female & Mexican & 117.8 & 3.4\\
\addlinespace
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Diastolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{77.8} & \cellcolor{gray!10}{3.1}\\
Clinic & Diastolic & Male & White & 75.4 & 3.1\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Diastolic} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{74.9} & \cellcolor{gray!10}{3.2}\\
Clinic & Diastolic & Female & Black & 72.4 & 3.0\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Diastolic} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{70.5} & \cellcolor{gray!10}{3.0}\\
Clinic & Diastolic & Female & Mexican & 69.2 & 3.0\\
\bottomrule
\end{tabular}
\end{table}

### Correlations between measurements {#sec:correlations}
In Figure \ref{fig:SD-mean} see that there is relatively little correlation between empirical SD and empirical mean SD for the different BP types and places. This is reassuring, as it avoids the possibility of a collinearity effect confounding the sampling of mean and SD, which are being treated as independent covariates in the model.


![(\#fig:SD-mean)Scatterplot of individual mean BP against individual SD of BP](NHANES-BP-Appendix_files/figure-latex/SD-mean-1.pdf) 

<!--#```{r Delta-mean,fig.pos='H', fig.cap='Scatterplot of individual mean BP against individual absolute difference between Clinic and Home mean', echo=FALSE,results='asis',message=FALSE,warning=FALSE}
#bp.sub <- data.frame(Mean=bp.data$TotalMean[bp.data$BPplace=='Clinic'],Delta=bp.data$Delta[bp.data$BPplace=='Clinic'], #BPtype=bp.data$BPtype[bp.data$BPplace=='Clinic']) %>%   # Home and clinic are identical
#  group_by(BPtype)
  
# bp.data.cor <- bp.sub %>%
#   summarize(correlationDelta = cor(Mean, Delta, use = "complete.obs"))
# # Scatterplot of mean BP against SD of BP from bp.data
#     ggplot(bp.sub %>% group_by(BPtype), aes(x = Mean, y = Delta)) + 
#     geom_point(alpha=.03, color = 'navyblue' ) +
#     labs(title = "Scatterplot of Mean against |Delta|",
#          x = "Mean",y = "|Delta|") + ylim(0,20)+
#       geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add regression line
#     facet_grid(~ BPtype, scales = 'free_x') +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     geom_text(data = bp.data.cor2, aes(label = sprintf("Cor: %.2f", correlationDelta), x = Inf, y = Inf), 
#                  vjust = "top", hjust = "right", inherit.aes = FALSE)
-->
\begin{table}

\caption{(\#tab:Delta-mean)Correlation between mean and Delta. Rows correspond to type of Delta, columns to type of mean.}
\centering
\begin{tabular}[t]{lrr}
\toprule
  & SysMean & DiasMean\\
\midrule
SysDelta & 0.137 & 0.014\\
DiasDelta & 0.347 & 0.128\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\centering
\caption{(\#tab:SD-mean2)Correlation between mean and SD. Rows correspond to type and location of SD, columns to type and location of mean.}
\centering
\begin{tabular}[t]{lrrrr}
\toprule
  & Clinic Sys Mean & Home Sys Mean & Clinic Dias Mean & Home Dias Mean\\
\midrule
\cellcolor{gray!10}{Clinic Sys SD} & \cellcolor{gray!10}{0.269} & \cellcolor{gray!10}{0.254} & \cellcolor{gray!10}{0.095} & \cellcolor{gray!10}{0.106}\\
Home Sys SD & 0.157 & 0.177 & 0.047 & 0.073\\
\cellcolor{gray!10}{Clinic Dias SD} & \cellcolor{gray!10}{0.041} & \cellcolor{gray!10}{0.030} & \cellcolor{gray!10}{-0.076} & \cellcolor{gray!10}{-0.037}\\
Home Dias SD & 0.045 & 0.060 & 0.010 & 0.003\\
\bottomrule
\end{tabular}
\end{table}

In Table \ref{tab:Delta-mean} we show the correlations between overall mean and absolute difference ($|\Delta|$) between clinic and home measurements.
The results are given as a $2\times2$ table, showing correlations within systolic and diastolic BP, and between the two.
The only moderately high correlation is between Systolic mean and Diastolic absolute Delta, which would correspond to a Variance Inflation Factor of 1.14.
While this is not directly relevant to the present Bayesian methodology, it suggests that this correlation should not substantially affect the estimation of the model coefficients.

In Table \ref{tab:SD-mean2} we show the correlations between mean and standard deviation for the three BP measures, considering all pairs of (Clinic,Home) and (Systolic,Diastolic). 
Finally, Table \ref{tab:SD-mean3} shows the correlations between systolic and diastolic, ranging over (Clinic,Home) and (Mean,SD).
(Some of the numbers here of course duplicate those in Table \ref{tab:SD-mean2}.)
Again, the correlations are too low to require any special treatment.

\begin{table}[!h]
\centering
\caption{(\#tab:SD-mean3)Correlation between diastolic and systolic summary statistics. Rows correspond to variables and locations for diastolic, columns to variables and locations for systolic.}
\centering
\begin{tabular}[t]{lrrrr}
\toprule
  & Clinic Sys SD & Clinic Sys Mean & Home Sys SD & Home Sys Mean\\
\midrule
\cellcolor{gray!10}{Clinic Dias SD} & \cellcolor{gray!10}{0.150} & \cellcolor{gray!10}{0.041} & \cellcolor{gray!10}{0.007} & \cellcolor{gray!10}{0.030}\\
Clinic Dias Mean & 0.095 & 0.497 & 0.047 & 0.378\\
\cellcolor{gray!10}{Home Dias SD} & \cellcolor{gray!10}{0.022} & \cellcolor{gray!10}{0.045} & \cellcolor{gray!10}{0.224} & \cellcolor{gray!10}{0.060}\\
Home Dias Mean & 0.106 & 0.405 & 0.073 & 0.547\\
\bottomrule
\end{tabular}
\end{table}

## Errors in blood pressure measurement or recording {#sec:errors}

The blood pressure measurement or recording errors were found particularly in the home measurements. 
While these did not destroy the usefulness of the home measurements, they did require some attention and decisions for how to work with these defects. 
We also consider them inherently interesting, and worth registering for future researchers working on these or similar data. 
In particular, the problem we have called “dependent replication” was entirely unexpected, although not unprecedented, and is of particular concern to researchers trying to estimate individual variation in clinically relevant measures.

### Last-digit preference {#sec:lastdigit}

Mild tendency for observers to prefer certain last digits in reporting BP measurements has been reported in other studies, though an analysis of the 1999 wave of NHANES reported no last-digit preference [@ostchega2003national].  
The last-digit preference in NHANES III, on the other hand, is substantial, with about 26.7% of all the clinic-measured systolic BP measurements ending in 0, but only about 31.9% ending in 4 or 6. Because the shifts due to last-digit preference are presumably small, we expect them to have little effect on the main effects that we are examining in this paper, but they do increase the probability of two measurements being rounded to the same value, something that needs to be taken into account in examining the problem of dependent replication.

\begin{table}[!h]
\centering
\caption{(\#tab:digit-summary)Summary data for BP end digits}
\centering
\begin{tabular}[t]{llrrrrr}
\toprule
Place & Sys/Dias & 0 & 2 & 4 & 6 & 8\\
\midrule
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{0.240} & \cellcolor{gray!10}{0.199} & \cellcolor{gray!10}{0.159} & \cellcolor{gray!10}{0.169} & \cellcolor{gray!10}{0.233}\\
Home & Diastolic & 0.186 & 0.179 & 0.198 & 0.217 & 0.219\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{0.267} & \cellcolor{gray!10}{0.188} & \cellcolor{gray!10}{0.160} & \cellcolor{gray!10}{0.159} & \cellcolor{gray!10}{0.226}\\
Clinic & Diastolic & 0.192 & 0.189 & 0.209 & 0.212 & 0.198\\
\bottomrule
\end{tabular}
\end{table}

### Dependent replicates {#sec:pseudorep}

While the protocol calls for each subject to have three independent BP measures taken, it is not impossible that the observers may have been influenced by one measure in recording the next.
This could happen in either direction: later measurements could be pulled closer to the first, or there could be an inclination to avoid repeated measures.
This is relevant, because erroneously repeated measures would artificially decrease the variance of the three measurements, and avoiding repeated measures would have the opposite effect.

The end-digit bias may be expected to have an effect here, since it influences the probability of two measurements being rounded to the same value.
We begin by noting the standard deviations for measurements of individual subjects as given in the column 'Mean of SD' in Table \ref{tab:sd-summary}.
The column 'Prob all rep' gives the theoretical probability that two of the three measurements for a subject would have the same value, if the measurements were independent and normally distributed with the given standard deviation (adjusted for the rounding), and assuming that rounding to particular digits is done in proportion to the fractions listed in Table \ref{tab:digit-summary}.
The column 'Prob 2 rep' gives the probability that two of the three measurements would have the same value, under the same conditions.
The column 'Frac all rep' gives the observed fraction of subjects for whom all three measurements were equal, and 'Frac 2 rep' gives the fraction for whom two of the three measurements were equal.
The observed fractions for three equal measurements are all very close to the theoretical probabilities, but the observed fractions for two equal measurements are substantially lower than the theoretical probabilities.
(For comparison, a 95\% probability range for the fraction of subjects with two equal measurements is about $\pm 0.008$.)

In Figure \ref{fig:examinerPlot}, we show the fraction of subjects with two equal measurements, by examiner, blocked by place and type.
We see that the fraction of subjects with two equal measurements varies substantially by examiner, and that the variation is greater for the systolic than for the diastolic measurements.



\begin{table}[!h]
\centering
\caption{(\#tab:sd-summary)Summary data for repeated measures}
\centering
\begin{tabular}[t]{llrrrrr}
\toprule
Place & Sys/Dias & Mean of SD & Frac all rep & Prob all rep & Frac 2 rep & Prob 2 rep\\
\midrule
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{2.739} & \cellcolor{gray!10}{0.050} & \cellcolor{gray!10}{0.048} & \cellcolor{gray!10}{0.432} & \cellcolor{gray!10}{0.510}\\
Home & Diastolic & 2.343 & 0.063 & 0.061 & 0.488 & 0.564\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{3.775} & \cellcolor{gray!10}{0.024} & \cellcolor{gray!10}{0.028} & \cellcolor{gray!10}{0.355} & \cellcolor{gray!10}{0.400}\\
Clinic & Diastolic & 3.082 & 0.028 & 0.036 & 0.414 & 0.459\\
\bottomrule
\end{tabular}
\end{table}


We show the fraction of subjects with two equal measurements in Figure {fig:examinerPlot}, split by examiner, blocked by place and type.
We see that the fraction of subjects with two equal measurements varies substantially by examiner, and that the variation is greater for the systolic than for the diastolic measurements.






![(\#fig:examinerPlot)Number of subjects with 2 equal measurements by examiner, blocked by place and type. Red band shows 95% probability range. Vertical green dashed line shows expected fraction; blue dotted line shows observed fraction over all examiners.](NHANES-BP-Appendix_files/figure-latex/examinerPlot-1.pdf) 

![(\#fig:examinerPlot3)Number of subjects with 3 equal measurements by examiner, blocked by place and type. Red band shows 95% probability range. Vertical green dashed line shows expected fraction; blue dotted line shows observed fraction over all examiners.](NHANES-BP-Appendix_files/figure-latex/examinerPlot3-1.pdf) 

In Figure \ref{fig:examinerPlot3}, we show the fraction of subjects with three equal measurements, by examiner, blocked by place and type.
Relative to the expected random fluctuations, we see that there is even more variation among the examiners.
One examiner (3001) produced consistently excessive numbers of triple repeats in Home measurements, and a deficit of triple repeats in Clinic measurements.

One further point to explore is the position of the two equal measures in a group of three.
If there are three independent measures, with two equal, each of the three has equal probability of being the odd one out.
On the other hand, if there is a trend in the measurements, then the second is least likely to be the odd one out.

In fact, what we observe is that it is the third measurement that is least likely to differ from the other two, while the first is most likely.
This is what we would expect if examiners sometimes either intentionally copied the second measurement into the space for the third, or unintentionally allowed themselves to be influenced into observing the same number.
The proportions are listed in Table \ref{tab:proportionChisq}, together with chi-squared tests for difference from the expected equal proportions for each site and type.
On the other hand, if there is a trend in the measurements, then the second is least likely to be the odd one out, which is also not what we see.

We see that there is a huge deviation from the expected proportions in the Home measurements, but less in the Home measurements, and more deviation in Systolic than in Diastolic measurements.

\begin{table}[!h]
\centering
\caption{(\#tab:proportionChisq)Chi-square test for difference between observed proportions (all examiners), stratified by place and type}
\centering
\begin{tabular}[t]{llrrrrl}
\toprule
Place & Sys/Dias & Freq1 & Freq2 & Freq3 & ChiSq & p-value\\
\midrule
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{2657} & \cellcolor{gray!10}{2149} & \cellcolor{gray!10}{1522} & \cellcolor{gray!10}{306.0} & \cellcolor{gray!10}{3.57e-67}\\
Home & Diastolic & 2864 & 2541 & 1746 & 278.0 & 4.3e-61\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{1905} & \cellcolor{gray!10}{1702} & \cellcolor{gray!10}{1594} & \cellcolor{gray!10}{28.8} & \cellcolor{gray!10}{5.57e-07}\\
Clinic & Diastolic & 2172 & 1992 & 1906 & 18.2 & 1.12e-04\\
\bottomrule
\end{tabular}
\end{table}
To explore this further, we can look at the proportions of first, second and third measurements from each examiner that are different from the other two.
The results of a chi-squared test for each examiner (stratified by site and type of BP) for difference from the expected equal proportions are shown in Figure \ref{fig:proportionChisq3}.
The dashed line represents a p-value of 0.001.
Here we see that the Home measurements are extremely variable, while the Clinic measurements are quite consistent with the expected proportions, with the single exception of examiner 3004, who is far from the expected equal proportions in all categories of measurement.

![(\#fig:proportionChisq3)Proportions of first, second and third measurements from each examiner that are different from the other two, by place and type. Chi-squared value for difference from expected proportions. Dashed line represents p-value 0.001.](NHANES-BP-Appendix_files/figure-latex/proportionChisq3-1.pdf) 

Given that the position of the differing measure clearly differs from the expected equal proportions, we might ask whether the examiners agree on a common proportion, suggesting that there might be some underlying systematic (observer-independent) reason for the differing measurements.
In Table \ref{tab:ternaryChi2} we show the results of a chi-squared test for equality of observed proportions among the examiners, stratified by place and type.
Interestingly, we see here that the examiners are fairly consistent in their proportions for the Home measures, but not for the Clinic measures.


\begin{table}[!h]
\centering
\caption{(\#tab:ternaryChi2)Chi-square test for difference between observed proportions among the examiners, stratified by place and type}
\centering
\begin{tabular}[t]{llrl}
\toprule
Place & Sys/Dias & ChiSq & p-value\\
\midrule
\cellcolor{gray!10}{Home} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{32.8} & \cellcolor{gray!10}{1.69e-01}\\
Home & Diastolic & 38.9 & 5.01e-02\\
\cellcolor{gray!10}{Clinic} & \cellcolor{gray!10}{Systolic} & \cellcolor{gray!10}{60.0} & \cellcolor{gray!10}{1.67e-04}\\
Clinic & Diastolic & 110.0 & 3.05e-12\\
\bottomrule
\end{tabular}
\end{table}

![(\#fig:ternaryPlot)Ternary plot of the position of the measurement that is unique, among subjects with 2 equal measurements. V1 is the fraction with the first distinct, V2 is the fraction with the second distinct, V3 is the fraction with the third distinct.](NHANES-BP-Appendix_files/figure-latex/ternaryPlot-1.pdf) 

Looking at a ternary plot Figure \ref{fig:ternaryPlot} for the proportions from the 13 different examiners, we see very clearly the bias toward having the last two measures agree, for almost all examiners, and examiner 3004 (marked larger) standing out as a clear outlier.


Overall, we can only conclude that there are clearly some irregularities in the BP measurement process, but we cannot identify a specific structure to them, or propose a remedy.
As the irregularities are not very large, we will proceed with the analysis without attempting to correct for them.


### Missing or implausible measurements {#sec:BPrange}

Some of the reported measures were extremely implausible, particularly for diastolic BP. NAsubjects had at least one diastolic BP measure recorded as 0, in addition to the 3916 subjects who were missing at least one measurement. We excluded all of these subjects, and indeed any subject who had at least one measurement recorded outside the ranges (40,140) for diastolic and (60,250) for systolic BP, as recommended by the CDC [@CDCBP]. 
There was just one subject with systolic BP measures that were too low, but NA subjects with low diastolic BP (in addition to those with measures recorded as 0). 
One subject was excluded for diastolic BP 156, and three were excluded for systolic BP that was too high, with the maximum being 264. 

\newpage

<!-- ## Exploratory data analysis -->
<!-- The empirical means of the home and clinic measures in population B are tabulated in Table 1. We note that the home measures are systematically higher than the clinic measures, within every demographic group, with greater differences for subjects who are white or Mexican, and female.  -->
<!-- The average difference is about 2.2 for diastolic and 2.7 for systolic, which is small compared with the general range of the differences, which have SD of  -->
<!-- 10.5 (diastolic) and 14.8 (systolic). -->

<!-- ```{r summaries, include=TRUE, echo=FALSE,results='asis',message=FALSE,warning=FALSE} -->
<!-- # Printing the table using kable -->
<!--     cat(kable(mean_sd_summary,format="latex", escape = F,booktabs = T, digits=1, -->
<!--               linesep = rep(c(rep("",5),"\\addlinespace"),4), -->
<!--                   col.names=c("Sys/Dias", "Place","Sex", "Ethnicity", 'Mean', 'SD'),caption = 'Summary data for blood pressure') %>%  kable_styling(latex_options = "hold_position") ) -->
<!--     cat('\n') -->
<!-- ``` -->

<!-- ### Correlations between measurements {#sec:correlations} -->
<!-- We see that there is relatively little correlation between empirical SD and empirical mean SD for the different BP types and places. This is reassuring, as it avoids the possibility of a collinearity effect confounding the sampling of mean and SD, which are being treated as independent covariates in the model. -->

<!-- ```{r SD-mean, dpi = 1,fig.pos='H', fig.width=8,fig.height=6, fig.cap='2D density plot of individual mean BP against individual SD of BP', echo=FALSE,results='asis',message=FALSE,warning=FALSE} -->
<!-- # 2D density plot of mean BP against SD of BP from bp.data -->
<!-- ggplot(bp.data, aes(x = Mean, y = SD)) + -->
<!--   geom_bin2d() + -->
<!--   labs(title = "2D Density Plot of Mean against SD", -->
<!--        x = "Mean",y = "SD") + ylim(0,15)+ -->
<!--   scale_fill_continuous(type = "viridis") + -->
<!--   theme_bw()+ -->
<!-- geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add regression line -->
<!--   facet_grid(BPtype ~ BPplace, scales = 'free_x') + -->
<!--   theme(plot.title = element_text(hjust = 0.5)) + -->
<!--   geom_text(data = bp.data.cor, aes(label = sprintf("Cor: %.2f", correlation), x = Inf, y = Inf),  -->
<!--             vjust = "top", hjust = "right", inherit.aes = FALSE) -->

<!-- ``` -->





# Appendix B -- Model details
This appendix aims to add more detail about the numerical modelling than was provided in the article. This is to ensure that the research methods are transparent and entirely reproducible. The numerical modelling presented in this paper was performed using R combined with Rstan. More detail will be provided here about the model, about the specific methodology used to parameterize the model, and more results are provided that were not included in the main text.

## The Statistical Model
The model used in this research is built from the theory of joint modelling of longitudinal and time-to-event data. This will be described in detail later on in this section, however, in brief, this allows the simultaneous modelling of both longitudinal observation data (in this article, this is blood pressure measurements) and also the time-to-event outcome. 
In this research the event of interest is either death from any cause, or death from specifically cardiovascular or cerebrovascular causes. We henceforth will refer to this latter mortality as CVD.
In the latter case, death from a different cause is treated as a noninformative censoring event.

### Survival Analysis (Time-to-Event)

The basic survival model is a Gompertz hazard rate with proportional hazards influences of the blood pressure covariates. 
The Gompertz equation 
\begin{equation}\label{gompertz}
h_0(t)=B\exp{\left(\theta(x+T)\right)},
\end{equation}
describes the baseline hazard of the population to a particular risk, which, for this article, investigates CVD mortality specifically, as well as studying mortality risk in general. $x\in\mathbb{N^N}$ is the age of the individual at the initial interview time, for $N$ the number of individuals, and $T\in\mathbb{R}^{+,N}$ the time since the individual entered the survey. 
Note that both $B$ and $\theta$ have 6 different values, depending on the sex reported at the initial interview -- female or male --- or the race --- black, white or 'other'. 
Note that 'other' in the race category is a combination of all non-black or non-white racial identities, such as Hispanic populations. 
The log-linear proportional hazards model links the covariates of the model (mean systolic blood pressure, variance in the diastolic blood pressure, etc) to the survival outcome of the individual via the equation
\begin{equation}\label{prophaz}
h(t)=h_0(t)\exp{\left(\boldsymbol{\beta}\cdot(\boldsymbol{X}-\hat{\boldsymbol{X}})\right)},
\end{equation}
where $\boldsymbol{X}\in\mathbb{R}^{+,N\times d}$ is a vector of summary statistics of the blood pressure measurements of individual covariates in our model, $\hat{\boldsymbol{X}}\in\mathbb{R}^{+,d}$ is the centering of the covariates such that the equation $\sum_i^N \exp{(\boldsymbol{\beta}\cdot(\boldsymbol{X}-\hat{\boldsymbol{X}}))}=0$ is approximately satisfied (more on this later), and $\boldsymbol{\beta}\in\mathbb{R}^d$ implies the strength of the influence of the covariate on the mortality risk. 
The majority of mortality events are censored --- not yet known at the time of data collection --- the censoring indicator being notated as $\delta\in \{0,1\}$.
When CVD mortality is the event being analysed, deaths due to other causes are treated as noninformative censoring events.
In this study, we explored the following covariates:



\begin{table}[!h]
\centering
\caption{(\#tab:runnumers)Explanations of the different models simulated in this work, according to run number.}
\centering
\begin{tabular}[t]{lll}
\toprule
Variable Name & Support & Description\\
\midrule
\cellcolor{gray!10}{$FRS-1998$} & \cellcolor{gray!10}{$R^N$} & \cellcolor{gray!10}{1998 version of the FRS score}\\
$FRS-ATP$ & $R^N$ & ATP version of the FRS score\\
\cellcolor{gray!10}{$M_S$} & \cellcolor{gray!10}{$R^{+,N}$} & \cellcolor{gray!10}{Mean systolic blood pressure}\\
$M_D$ & $R^{+,N}$ & Mean diastolic blood pressure\\
\cellcolor{gray!10}{$\Delta_S$} & \cellcolor{gray!10}{$R^{+,N}$} & \cellcolor{gray!10}{Semi-difference between Home and Clinic mean systolic blood pressure}\\
$\Delta_D$ & $R^{+,N}$ & Semi-difference between Home and Clinic mean diastolic blood pressure\\
\cellcolor{gray!10}{$\sigma_{\{S,H\}}$} & \cellcolor{gray!10}{$R^{+,N}$} & \cellcolor{gray!10}{Standard deviation of the systolic blood pressure taken at home}\\
$\sigma_{\{D,H\}}$ & $R^{+,N}$ & Standard deviation of the diastolic blood pressure taken at home\\
\cellcolor{gray!10}{$\sigma_{\{S,C\}}$} & \cellcolor{gray!10}{$R^{+,N}$} & \cellcolor{gray!10}{Standard deviation of the systolic blood pressure taken at the clinic}\\
$\sigma_{\{D,C\}}$ & $R^{+,N}$ & Standard deviation of the diastolic blood pressure taken at the clinic\\
\cellcolor{gray!10}{$\tau_{\{S,H\}}$} & \cellcolor{gray!10}{$R^{+,N}$} & \cellcolor{gray!10}{Precision of the systolic blood pressure taken at home}\\
$\tau_{\{D,H\}}$ & $R^{+,N}$ & Precision of the diastolic blood pressure taken at home\\
\cellcolor{gray!10}{$\tau_{\{S,C\}}$} & \cellcolor{gray!10}{$R^{+,N}$} & \cellcolor{gray!10}{Precision of the systolic blood pressure taken at the clinic}\\
$\tau_{\{D,C\}}$ & $R^{+,N}$ & Precision of the diastolic blood pressure taken at the clinic\\
\bottomrule
\end{tabular}
\end{table}

Please note that the last four elements of this list, the precision values, were only carried out to ensure model consistency with the use of standard deviation instead. 
Note as well that the $\Delta$ covariates, representing the medium-term variability, enter into the log relative risk sum as an **absolute value**.

For the parametrization of this model, we assume that the Gompertz parameters and the parameters in the linear predictor term have prior distributions as follows:
\begin{equation}\label{priorsS}
\begin{aligned}
  \boldsymbol{B}\sim\mathcal{N}(\mu_B,\sigma_B),\\
  \boldsymbol{\theta}\sim\mathcal{N}(\mu_\theta,\sigma_\theta),\\
  \boldsymbol{\beta}\sim \mathcal{N}(\mu_\beta,\sigma_\beta),
\end{aligned}
\end{equation}

The likelihood for this Gompertz proportional hazards model, over all individuals in the census, is as follows:
\begin{equation}\label{likesurv}
L_S(\boldsymbol{v},\boldsymbol{\delta})=\prod_i^N f(v_i,\delta_i|B_i,\theta_i,\beta_i,\boldsymbol{X},\hat{\boldsymbol{X}})=\prod_i^N h(v_i|B_i,\theta_i,\beta_i,\boldsymbol{X},\hat{\boldsymbol{X}})^{\delta_i} \exp{\left( -\sum_i^N H(v_i|B_i,\theta_i,\beta_i,\boldsymbol{X},\hat{\boldsymbol{X}}) \right)},
\end{equation}
with $H(v)=\int_0^v h(w) \mathrm{d}w$ the cumulative hazard.

### Longitudinal Modelling

The mortality hazard rates are assumed to be influenced by individual-level blood pressure means and variability characteristics.
These characteristics are not directly observed, but are inferred from their influence on the individual blood pressure measurements, which have been observed.
Let $Y_i(t_j)$ be the observed blood pressure for patient $i$ at time $t_j$, for the individual $i\in 1,2,...,N$ and the number of blood pressure measurements per individual $j\in 1,2,...,k$. Due to the fact that the blood pressure measurement data was taken at both the home and clinic (written using subscripts H and C, respectively), with approximately 6 months between these two measurements, we model the blood pressure using the following model, assuming the diastolic $Y_{i}^D$ and systolic $Y_{i}^S$ blood pressure to be Gaussian-distributed:
\begin{equation}\label{bp}
\begin{aligned}
  (Y_{i}^D)_{H} \sim \mathcal{N}(M_i^D+\Delta_i^D,(\sigma_i^D)_H),\\
  (Y_{i}^D)_{C} \sim \mathcal{N}(M_i^D-\Delta_i^D,(\sigma_i^D)_C),\\
  (Y_{i}^S)_{H} \sim \mathcal{N}(M_i^S+\Delta_i^S,(\sigma_i^S)_H),\\
  (Y_{i}^S)_{C} \sim \mathcal{N}(M_i^S-\Delta_i^S,(\sigma_i^S)_C),
\end{aligned}
\end{equation}
where superscripts $D$ and $S$ refer to diastolic and systolic blood pressure, respectively. 

The blood pressure characteristics --- the individual-level parameters --- are themselves distributed according to a hierarchical model, determined by population-level parameters (also called ``hyperparameters''):
\begin{equation}\label{priorsL}
\begin{aligned}
  M_i^{\{D,S\}}\sim \mathcal{N}(\mu_M^{\{D,S\}},\sigma_M^{\{D,S\}}),\\
  \Delta_i^{\{D,S\}}\sim \mathcal{N}(\mu_D^{\{D,S\}},\sigma_D^{\{D,S\}}),\\
  \sigma_{i,C}^{\{D,S\}}\sim \Gamma(r_C^{\{D,S\}},\lambda_C^{\{D,S\}}),\\
  \sigma_{i,H}^{\{D,S\}}\sim \Gamma(r_H^{\{D,S\}},\lambda_H^{\{D,S\}}).
\end{aligned}
\end{equation}
The longitudinal outcome modelling therefore aims to infer these hyperparameters
\begin{equation}
  \Theta=\left\{\mu_M^{\{D,S\}},\mu_D^{\{D,S\}},\sigma_M^{\{D,S\}},\sigma_D^{\{D,S\}},r_C^{\{D,S\}},\lambda_C^{\{D,S\}},r_H^{\{D,S\}},\lambda_H^{\{D,S\}}\right\},
\end{equation}
and to use the implied uncertainty about the individual-level parameters to inform the inference about the survival parameters.
The likelihood for the longitudinal measurements is therefore (combining the systolic and diastolic into a single parameter for simplicity):
\begin{equation}\label{likelong}
  L_L(\Theta|Y)=\prod_{i=1}^N\left(\prod_{j=1}^{k}f(y_{ij}|M_i,\Delta_i,\sigma_i)\right)f(M_i|\mu_M,\sigma_M)f(\Delta_i|\mu_D,\sigma_D)f(\tau_{i,C}|r_C,\lambda_C)f(\tau_{i,H}|r_H,\lambda_H)
\end{equation}

### Combined Hierarchical Model

Combining the longitudinal outcome and time-to-event partial likelihoods, and for a given parameter space value of $\Omega=\{\beta,B,\theta\}\cup \Theta$, the joint likelihood is
\begin{equation}
\begin{split}
  L(\Omega|Y)=\prod_{i=1}^N\left(\prod_{j=1}^{k}f(y_{ij}|M_i,\Delta_i,\sigma_i)\right)f&(v_i,\delta_i|B_i,\theta_i,\beta_i,\boldsymbol{X},\hat{\boldsymbol{X}})f(M_i|\mu_M,\sigma_M)\\
  &f(\Delta_i|\mu_D,\sigma_D)f(\tau_{i,C}|r_C,\lambda_C)f(\tau_{i,H}|r_H,\lambda_H).
  \end{split}
\end{equation}
One approach to estimating the complete set of hyperparameters
\begin{equation}
  \Omega_H=\{\mu_B,\sigma_B,\mu_\theta,\sigma_\theta,\mu_\beta,\sigma_\beta,\mu_M^{\{D,S\}},\sigma_M^{\{D,S\}},\mu_D^{\{D,S\}},\sigma_D^{\{D,S\}},r_C^{\{D,S\}},\lambda_C^{\{D,S\}},r_H^{\{D,S\}},\lambda_H^{\{D,S\}}\}
\end{equation}
is to impose a higher-level prior distribution, and use the machinery of Bayesian inference to produce posteriors for everything.
This approach runs into computational difficulties, which have led us to a two-stage `empirical Bayes' approach, where the hyperparameters for the longitudinal model are first fixed by a maximum-likelihood calculation, after which the remaining hyperparameters and individual-level parameters can be estimated with Bayesian machinery. 
For the time-to-event parameters we choose flat hyperpriors, selecting the hyperparameters $\mu_B=\mu_\theta=\mu_\beta=0$,  $\sigma_B=\sigma_\theta=2$, and $\sigma_\beta=100$.

### The modelling variants

In this article, we researched into 16 variants of the model-fitting problem, but focussed mainly on 8 of them. 
The 8 main models use the standard deviation, $\sigma$, as the measure of the influence of blood-pressure variability on mortality. 
We also produced the same 8 models but using precision, $\tau=1/\sigma^2$, as the measure of the influence of blood-pressure variability on mortality. 
However, this was only to ensure that there were no differences between the use of one over the other. 
Throughout the remainder of this appendix, we refer to the 8 main models using the following run numbers:

\begin{enumerate}
\item All participants (14,654), using mean systolic and diastolic blood pressure (not FRS) in the linear predictor term, with the outcome data as death specifically from CVD.
\item All participants (14,654), using mean systolic and diastolic blood pressure (not FRS) in the linear predictor term, with the outcome data as all-causes of death.
\item Only participants that had data from which FRS values could be computed (N=9,008) --- the ``FRS population'' but using mean systolic and diastolic blood pressure (not FRS) in the linear predictor term, with the outcome data as death specifically from CVD.
\item FRS population, but using mean systolic and diastolic blood pressure (not FRS) in the linear predictor term, with the outcome data as all-causes of death.
\item FRS population, and using the FRS ATP-III value in the linear predictor term, with the outcome data as death specifically from CVD.
\item FRS population, and using the FRS ATP-III value in the linear predictor term, with the outcome data as all-causes of death.
\item FRS population, and using the FRS 1998-version value in the linear predictor term, with the outcome data as death specifically from CVD.
\item FRS population, and using the FRS 1998-version value in the linear predictor term, with the outcome data as all-causes of death.
\label{tab:runnums}
\end{enumerate}

We also include Directed Acyclical Graph (DAG) sketches to help visualize the different models, as shown in figures \ref{fig:DAGmean} and \ref{fig:DAGFRS}. 
In order to read the DAGs, note that each square background layer that appears as a stack of layers represents different measured outcomes that were made in the first wave of the survey. 
The outcome variables measured are represented by a square-shaped text box, and a parameter of the model is represented by a circular-shaped text box. If either a square or circular text box is placed on top of a stacked rectangular layer, it means that multiple values of that variable (as many as there are layers to the stack) are either measured (for outcome variables) or simulated (for parameters of the model). Please note that the number of layers in the stack is written in the text box that does not contain a frame which is intentionally displayed on top of the stacked layer that it represents. For example, $i=1,...,N$. 
Finally, the direction of the arrows implies causality assumed in the model.

The distribution of the blood pressure parameters in the population are derived from the model, and are summarised with other outputs of the model in table \ref{tab:Mean-SD} of Appendix C.

![An illustration of the DAG of the mean blood pressure-based model presented in this article. ](./DAG_Mean2.png){#fig:DAGmean}

![An illustration of the DAG of the FRS-based model presented in this article. ](./DAG_FRS2.png){#fig:DAGFRS}

<!--
For convenience, we provide an overview of the different blood pressure values in the full and FRS populations in tables \ref{tab:BloodFull} and \ref{tab:BloodFRS}.


\begin{table}

\caption{(\#tab:BloodFull)Parameters for distribution of blood pressure, for the full population}
\centering
\begin{tabular}[t]{llrr}
\toprule
BP & Variable & Mean & SD\\
\midrule
Systolic & Overall Mean & 125.36 & 19.48\\
Systolic & Delta & 5.24 & 4.86\\
Systolic & Home SD & 2.74 & 2.05\\
Systolic & Clinic SD & 3.78 & 2.61\\
\addlinespace\\
Diastolic & Overall Mean & 74.34 & 10.29\\
\addlinespace
Diastolic & Delta & 3.90 & 3.14\\
Diastolic & Home SD & 2.34 & 1.75\\
Diastolic & Clinic SD & 3.08 & 2.08\\
\bottomrule
\end{tabular}
\end{table}


\begin{table}

\caption{(\#tab:BloodFRS)Parameters for distribution of blood pressure, for the FRS population}
\centering
\begin{tabular}[t]{llrr}
\toprule
BP & Variable & Mean & SD\\
\midrule
Systolic & Overall Mean & 125.89 & 18.33\\
Systolic & Delta & 5.23 & 4.71\\
Systolic & Home SD & 2.78 & 2.06\\
Systolic & Clinic SD & 3.78 & 2.52\\
\addlinespace\\
Diastolic & Overall Mean & 76.41 & 10.00\\
\addlinespace
Diastolic & Delta & 3.84 & 3.10\\
Diastolic & Home SD & 2.28 & 1.71\\
Diastolic & Clinic SD & 2.92 & 1.97\\
\bottomrule
\end{tabular}
\end{table}
-->


## Computational methodology

The methodology for this research can be split into three main sections: 1) calculating the empirical Bayes' parameters, 2) parameterizing the model using Hamiltonian Monte Carlo (HMC) and 3) re-centering the variables in the linear predictor equation. By applying empirical Bayes', Maximum Likelihood Estimates (MLEs) of some of the parameter distributions are provided. Note that the parameters estimated here are only the prior distribution of the global (not individual) blood pressure means and the variances, for both systolic and diastolic and home and clinic measurements. These estimates are then provided as prior distributions for the Stan MCMC simulations using HMC, where estimates can be made for all the parameter distributions of the model, given the specific centering applied. Finally, section (3) recalculates the centering values based on the previous MCMC iteration, and sets of the next iteration, while simultaneously checking for convergence in both the MCMC simulations and the centering values.

## Empirical Bayes Parameters

First, we extract the intervals for the digits in the blood pressure measurement recordings. Suppose the fractions of digits 0,2,4,6,8 are $b_0,b_2,b_4,b_6,b_8$.
Letting $B_0=0$ and $B_k=10\sum_{j=0}^{k-1}b_{2j}$ for $k=1,\dots,5$,
we want to choose a positive $a$ and place breaks at $-a+B_k$, so that measurements between $-a+B_k$ and $-a+B_{k+1}$ modulo
10 are assigned the final digit $2k$, for $k=0,\dots,4$.
We choose $a$ to minimise the total distance of the intervals from the rounded value:
$$
  \sum_{k=0}^4 \int_{-a+B_k}^{-a+B_{k+1}} \bigl| x-2k\bigr|\mathrm{d} x=\frac12\sum_{k=0}^4 \left(-a+B_k-2k\right)^2 + \left(-a+B_{k+1}-2k\right)^2,
$$
as long as $2k$ is in the appropriate interval. This is minimized at
$$
  a= \frac{1}{5}\left(B_1+B_2+B_3+B_4 - 15\right)=\sum_{j=0}^3 (8-2j) b_{2j} \, - 3.
$$



Next step, we fit the BP distribution parameters. We suppose that each individual has BP measures
$\tilde{y}_{ij}^l$ for $i=1,\dots,n$, $j=1,\dots,k$ (default $k=3$),
and $l=1,2$, which are rounded versions of 
$$
  y_{ij}^l \sim \mathcal{N}\bigl( \mu_i^l,(\tau_i^l)^{-1}\bigr),
$$
where 
\begin{align*}
\mu_i^1&=(M_i+\Delta_i)/2,\\
\mu_i^2&=(M_i-\Delta_i)/2,\\
M_i&\sim \mathcal{N}\bigl(m_M,\sigma^2_M) \text{ and }
  \Delta_i\sim \mathcal{N}\bigl(m_\Delta,\sigma^2_\Delta) \text{ independent,}\\
  \tau_i^l &\sim \mathrm{Gamma}(\alpha^l,\alpha^l/\theta^l ).
\end{align*}
(Note that $\alpha^l$ is the usual shape parameter,
while $\theta^l$ is the expectation.)

We wish to estimate the eight parameters 
$$
(m_M,m_\Delta,\sigma^2_M,\sigma^2_\Delta,\alpha^1,\theta^1,\alpha^2,\theta^2)
$$
We begin by assuming $y_{ij}^l$ observed directly. We estimate
by maximising the partial likelihood on the observations
\begin{align*}
  \bar{y}_{i+}&:= \frac{1}{2k} \sum_{j=1}^k \bigl( y_{ij}^1 + y_{ij}^2\bigr),\\
  \bar{y}_{i-}&:= \frac{1}{2k} \sum_{j=1}^k \bigl( y_{ij}^1 - y_{ij}^2 \bigr),\\
  s_i^l&:=  \frac{1}{k-1}\sum_{j=1}^k \Bigl( y_{ij}^l - \frac{1}{k} \sum_{j=1}^k y_{ij}^l \Bigr)^2.
\end{align*}
Note that 
$$
(k-1)s_i^l \tau_i^l =\sum_{j=1}^k \Bigl( z_{ij}^l - \frac{1}{k} \sum_{j=1}^k z_{ij}^l \Bigr)^2.
$$ 
where $z_{ij}^l$ are i.i.d.\ standard normal
is independent of $\tau_i^l$, thus has a chi-squared distribution
with $k-1$ degrees of freedom --- hence $\frac{k-1}{2}\cdot s_i^l\tau_i^l$ is
gamma distributed with parameters $(\frac{k-1}{2},1)$. Since $\frac{\alpha}{\theta}\tau_i^l$ is independent of $s_i^l\tau_i^l$, with $\mathrm{Gamma}(\alpha,1)$ distribution, we see that $\frac{(k-1)\theta}{2\alpha}s_i^l$ is the ratio of two independent gamma random variables, hence has beta-prime distribution with parameters $\left(\frac{k-1}{2}, \alpha \right)$, so log partial likelihood
$$
  \ell_{\operatorname{Beta}}(\alpha,\theta;s^l_\cdot)=n\alpha \log\frac{\alpha}{\theta}+n\log\Gamma\left(\alpha+\frac{k-1}{2}\right)-n\log\Gamma(\alpha)
  + \frac{k-1}{2} \sum_{i=1}^n \log s_i^l -\left(\alpha+\frac{k-1}{2}\right) \sum_{i=1}^n \log \left(s_i^l+\frac\alpha\theta\right).
$$
Note as well that these quantities $(k-1)s_i^l$ should correspond to empirically observed individual variances; hence we will compare these empirical variances (with imputed fractional parts) divided by the normalization factor $2\alpha/(k-1)\theta$ to the beta-prime distribution below as a goodness-of-fit test.

The partial Fisher Information has entries
\begin{align*}
 -\frac{\partial^2 \ell}{\partial \alpha^2} &=
   n\psi_1\left(\alpha\right) - n\psi_1\left(\alpha+\frac{k-1}{2}\right)
  - \frac{n}{\alpha} +\sum_{i=1}^n \frac{2\theta s_i^l + \alpha-(k-1)/2}{(\theta s_i^l + \alpha)^2}\\
-\frac{\partial^2 \ell}{\partial \theta^2} &=
   -\frac{n \alpha}{\theta^2} +\frac{\alpha}{\theta^2}\left(\alpha+\frac{k-1}{2}\right)\sum_{i=1}^n \frac{2\theta s_i^l + \alpha}{(\theta s_i^l + \alpha)^2}\\
-\frac{\partial^2 \ell}{\partial \theta\partial\alpha} &= \frac{n}{\theta}-
   \frac1\theta \sum_{i=1}^n \frac{\alpha^2+2\alpha\theta s_i^l+\frac{k-1}{2}\theta s_i^l}{(\theta s_i^l + \alpha)^2}.
\end{align*}
where $\psi_1$ is the trigamma function.

Let $(\hat\alpha^l,\hat\theta^l)$ be the maximum partial likelihood estimators. Conditioned on $(\tau_i^l)$ we have
\begin{align*}
  \bar{y}_{i+}&\sim \mathcal{N}\left(m_M, \sigma^2_M + \frac{1}{4k}\left( \frac{1}{\tau_i^1}+\frac{1}{\tau_i^2}\right)\right),\\
  \bar{y}_{i-}&\sim \mathcal{N}\left(m_\Delta,\sigma^2_\Delta + \frac{1}{4k}\left( \frac{1}{\tau_i^1}+\frac{1}{\tau_i^2}\right)\right).
\end{align*}
We would then have MLEs
\begin{align*}
  \hat{m}_M&= \frac{1}{n} \sum_{i=1}^n \bar{y}_{i+},\\
  \hat{m}_\Delta&= \frac{1}{n} \sum_{i=1}^n \bar{y}_{i-},
\end{align*}
which are approximately normally distributed, with means $m_M$ and $m_\Delta$ respectively, and conditional on $\tau_i^l$ standard errors
$$
  \frac{\sigma_M^2}{n}+\frac{1}{4kn^2} \sum_{i=1}^n (\tau_i^1)^{-1} + (\tau_i^2)^{-1} \quad \text{ and } \quad
  \frac{\sigma_\Delta^2}{n}+\frac{1}{4kn^2} \sum_{i=1}^n (\tau_i^1)^{-1} + (\tau_i^2)^{-1},
$$
which we may approximate --- with error on the order of $n^{-3/2}$ --- replacing the mean of $(\tau_i^l)^{-1}$ by its expected value $\theta^l/(\alpha^l-1)$ to obtain
\begin{align*}
  \mathrm{Var}(\hat{m}_M) &\approx \frac{\sigma_M^2}{n}+\frac{1}{4kn}\left( \frac{\theta^1}{\alpha^1-1}+ \frac{\theta^2}{\alpha^2-1}\right) \\
  \mathrm{Var}(\hat{m}_\Delta) &\approx \frac{\sigma_\Delta^2}{n}+\frac{1}{4kn}\left( \frac{\theta^1}{\alpha^1-1}+ \frac{\theta^2}{\alpha^2-1}\right) 
\end{align*}
Finally, conditioned on the $\tau_i^l$ we have that the random variables $\bar{y}_{i+}$ are normal with variance
$$
  \sigma_M^2+\frac{1}{4k}\left((\tau_i^1)^{-1} + (\tau_i^2)^{-1} \right),
$$
so the unconditional variance is the expected value, or
$$
  \sigma_M^2+\frac{1}{4k}\left(\frac{\theta_1}{\alpha_1-1}+ \frac{\theta_2}{\alpha_2-1} \right).
$$
This yields the estimators
\begin{align*}
  \hat\sigma_M^2 &=\frac{1}{n-1}\sum_{i=1}^n\left(\bar{y}_{i+}-n^{-1}\sum_{i=1}^n \bar{y}_{i+}\right)^2 - \frac{1}{4k}\left(\frac{\hat\theta_1}{\hat\alpha_1-1}+ \frac{\hat\theta_2}{\hat\alpha_2-1} \right),\\
  \hat\sigma_\Delta^2 &=\frac{1}{n-1}\sum_{i=1}^n\left(\bar{y}_{i-}-n^{-1}\sum_{i=1}^n \bar{y}_{i-}\right)^2 - \frac{1}{4k}\left(\frac{\hat\theta_1}{\hat\alpha_1-1}+ \frac{\hat\theta_2}{\hat\alpha_2-1} \right).
\end{align*}
Using the delta method, and the fact that the correlation between $\hat\alpha$ and $\hat\theta$ is small, we see that the variance of $\hat\theta/(\hat\alpha-1)$ is approximately
$$
  \frac{\sigma_\theta^2}{(\hat\alpha-1)^2} + \frac{\hat\theta^2\sigma_\alpha^2}{(\hat\alpha-1)^4},
$$
where $\sigma_\alpha$ and $\sigma_\theta$ are the standard errors for $\hat\alpha$ and $\hat\theta$ respectively. Define
$$
  \hat\sigma_{\alpha\theta}^2 := \frac{1}{16k^2} \left(\frac{\hat\sigma_{\theta_1}^2}{(\hat\alpha_1-1)^2} + \frac{(\hat\theta_1)^2\hat\sigma_{\alpha_1}^2}{(\hat\alpha_1-1)^4} + \frac{\hat\sigma_{\theta_2}^2}{(\hat\alpha_2-1)^2} + \frac{(\hat\theta_2)^2\hat\sigma_{\alpha_2}}{(\hat\alpha_2-1)^4}\right)
$$
so the standard errors for $\hat\sigma_M^2$ and $\hat\sigma_\Delta^2$ are approximately
\begin{align*}
  \operatorname{SE}\left(\hat\sigma_M^2\right)&\approx \Bigl(\frac{2\hat\sigma_M^4}{n} + \hat\sigma_{\alpha\theta}^2 \Bigr)^{1/2},\\
  \operatorname{SE}\left(\hat\sigma_\Delta^2\right)&\approx \Bigl(\frac{2\hat\sigma_\Delta^4}{n} + \hat\sigma_{\alpha\theta}^2 \Bigr)^{1/2}.
\end{align*}


<!-- ## Test whether parameters are being fit correctly -->
<!-- We simulate bootstrap data sets, find the average parameter estimate, and compare to the "true" parameters. -->
<!-- We also compare the average estimated SE to the "true" SE (which is the SD of the estimates). -->
\begin{table}[!h]
\centering
\caption{(\#tab:testparameters)Results of estimating parameters from simulated data from the whole population. 
    First column on top is the average parameter estimate from the simulations, second is the true parameter from which the simulations were made, third is the relative error.
    On bottom are the standard errors for the parameters: True is the theoretically computed standard error, SimAverage is the SD of the simulated parameter estimates, and RelError is the relative error.}
\centering
\begin{tabular}[t]{lllllll}
\toprule
\multicolumn{1}{c}{ } & \multicolumn{3}{c}{Systolic} & \multicolumn{3}{c}{Diastolic} \\
\cmidrule(l{3pt}r{3pt}){2-4} \cmidrule(l{3pt}r{3pt}){5-7}
  & SimAverage & True & RelError & SimAverage & True & RelError\\
\midrule
\addlinespace[0.3em]
\multicolumn{7}{l}{\textbf{Parameter estimates}}\\
\cellcolor{gray!10}{\hspace{1em}$m_M$} & \cellcolor{gray!10}{123} & \cellcolor{gray!10}{123} & \cellcolor{gray!10}{1.26e-04} & \cellcolor{gray!10}{72.3} & \cellcolor{gray!10}{72.3} & \cellcolor{gray!10}{-4.34e-04}\\
\hspace{1em}$m_\Delta$ & 1.33 & 1.35 & -0.0152 & 1.11 & 1.1 & 0.0131\\
\cellcolor{gray!10}{\hspace{1em}$\sigma^2_M$} & \cellcolor{gray!10}{377} & \cellcolor{gray!10}{377} & \cellcolor{gray!10}{7.77e-04} & \cellcolor{gray!10}{104} & \cellcolor{gray!10}{104} & \cellcolor{gray!10}{0.00298}\\
\hspace{1em}$\sigma^2_\Delta$ & 46.4 & 46.4 & -6.50e-04 & 21.9 & 21.9 & 0.00149\\
\cellcolor{gray!10}{\hspace{1em}$\alpha_H$} & \cellcolor{gray!10}{2.17} & \cellcolor{gray!10}{2.16} & \cellcolor{gray!10}{0.00492} & \cellcolor{gray!10}{2.31} & \cellcolor{gray!10}{2.32} & \cellcolor{gray!10}{-0.0036}\\
\hspace{1em}$\theta_H$ & 0.149 & 0.15 & -0.00411 & 0.195 & 0.194 & 0.00181\\
\cellcolor{gray!10}{\hspace{1em}$\alpha_C$} & \cellcolor{gray!10}{2.59} & \cellcolor{gray!10}{2.57} & \cellcolor{gray!10}{0.0085} & \cellcolor{gray!10}{2.79} & \cellcolor{gray!10}{2.76} & \cellcolor{gray!10}{0.0111}\\
\hspace{1em}$\theta_C$ & 0.0753 & 0.0751 & 0.00227 & 0.109 & 0.109 & -0.00403\\
\addlinespace[0.3em]
\multicolumn{7}{l}{\textbf{Parameter SE estimates}}\\
\cellcolor{gray!10}{\hspace{1em}$m_M$} & \cellcolor{gray!10}{0.161} & \cellcolor{gray!10}{0.167} & \cellcolor{gray!10}{-0.0352} & \cellcolor{gray!10}{0.0851} & \cellcolor{gray!10}{0.0768} & \cellcolor{gray!10}{0.109}\\
\hspace{1em}$m_\Delta$ & 0.058 & 0.0244 & 1.37 & 0.0404 & 0.0311 & 0.296\\
\cellcolor{gray!10}{\hspace{1em}$\sigma^2_M$} & \cellcolor{gray!10}{4.41} & \cellcolor{gray!10}{3.32} & \cellcolor{gray!10}{0.325} & \cellcolor{gray!10}{1.22} & \cellcolor{gray!10}{1.42} & \cellcolor{gray!10}{-0.139}\\
\hspace{1em}$\sigma^2_\Delta$ & 0.551 & 0.686 & -0.197 & 0.265 & 0.284 & -0.0655\\
\cellcolor{gray!10}{\hspace{1em}$\alpha_H$} & \cellcolor{gray!10}{0.0568} & \cellcolor{gray!10}{0.0488} & \cellcolor{gray!10}{0.166} & \cellcolor{gray!10}{0.0633} & \cellcolor{gray!10}{0.0634} & \cellcolor{gray!10}{-0.00127}\\
\hspace{1em}$\theta_H$ & 0.00211 & 0.00173 & 0.214 & 0.00272 & 0.00286 & -0.0485\\
\cellcolor{gray!10}{\hspace{1em}$\alpha_C$} & \cellcolor{gray!10}{0.0769} & \cellcolor{gray!10}{0.104} & \cellcolor{gray!10}{-0.258} & \cellcolor{gray!10}{0.087} & \cellcolor{gray!10}{0.0534} & \cellcolor{gray!10}{0.628}\\
\hspace{1em}$\theta_C$ & 0.00104 & 0.00137 & -0.245 & 0.00148 & 0.00133 & 0.113\\
\bottomrule
\end{tabular}
\end{table}


<!-- ## Multiple imputation for the real data -->


Now we compute the combined variance. For a parameter like $\alpha$ we estimate the variance of $\hat\alpha$ by
\newcommand{\E}{\mathbb{E}}
\renewcommand{\P}{\mathbb{P}}
$$
  \mathrm{Var}(\hat\alpha) = \mathbb{E}\bigl[ \mathrm{Var}\left(\hat\alpha\, |\, I\right)\bigr] + \mathrm{Var}\left(\mathbb{E} \left[ \hat\alpha\, |\, I \right]\right).
$$
Here $I$ represents the randomly imputed fractional part. 
We can estimate the first term by averaging the estimated variance (from Fisher Information) over all random imputations.
We estimate the second term by the variance of the $\alpha$ estimates over imputations. Note that this is not quite right, since what we really
want the variance of is $\alpha_0(I)$ --- effectively, the "true" parameter consistent with the imputation. This is a plug-in estimate,
as is the Fisher Information estimate of the variance.  

### Estimates for whole population
The estimates of the empirical Bayes parameters together with their standard errors are given in the column labelled "True" in Table \ref{tab:testparameters}.

<!--
We then compute the residuals. We define the deviance for an individual $i$ with observations $(Y_i)$
given the hyperparameters $h=(m_M,m_\Delta,\sigma^2_M,\sigma^2_\Delta,\alpha^H,\theta^H,\alpha^C,\theta^C)$
$$
  D= \sum_{i=1}^n \log \mathbb{P}\left\{ \mathbf{Y}_{i}\,|\, \text{hyperparameters}=h\right\}.
$$
\newcommand{\wtb}{\widetilde\mathbf}
Since the $\mathbf{Y}_i$ are independent conditioned on $h$,
\begin{align*}
D&= \sum_{i=1}^n \log \E_h\left[ \P\left\{ \mathbf{Y}_i \, |\, M_i,\Delta_i,\tau_i^{C},\tau_i^H \right\} \right]\\
    &\approx \sum_{i=1}^n \log \frac1R\sum_{r=1}^R \left[ \P\left\{ \mathbf{Y}_i \, |\, M_{i,r},\Delta_{i,r},\tau_{i,r}^{C},\tau_{i,r}^{H} \right\}\right] \frac{\pi_h(M_{i,r},\Delta_{i,r},\tau_{i,r}^{C},\tau_{i,r}^{H} )}{q(M_{i,r},\Delta_{i,r},\tau_{i,r}^{C},\tau_{i,r}^{H} \, | \, h,\, \mathbf{Y}_i)},
\end{align*}
where $(M_{i,r},\Delta_{i,r},\tau_{i,r}^{C},\tau_{i,r}^{H})$ are independent samples from a distribution $q$ that may depend
on $\mathbf{Y}_i$ and $h$, and $\pi_h$ is the true density of those individual parameters given hyperparameters $h$.
-->
<!-- 
We can try estimating this simply by Monte Carlo sampling of the individual parameters.
-->
<!--
We estimate this by importance sampling on the four individual parameters $(M_i,\Delta_i,\tau_i^H, \tau_i^C)$ from an approximate posterior.
We have, conditioned on the observations of sample variances for clinical and home measures $S_{Ci}^2$ and $S_{Hi}^2$,
\begin{align*}
  \tau_i^H &\sim \mathrm{Gamma}\left( \alpha^H+\frac{k-1}{2},\, \beta^H+\frac{k-1}{2} S_{Hi}^2 \right),\\
  \tau_i^C &\sim \mathrm{Gamma}\left( \alpha^H+\frac{k-1}{2},\, \beta^H+\frac{k-1}{2} S_{Hi}^2 \right).
\end{align*}
Then, recalling the definitions of $y_{i+}$ and $y_{i-}$,
conditioned on $\tau_i^H$ and $\tau_i^C$ we have
\begin{align*}
  M_i &= \mathcal{N} \left( (\tau_i^M)^{-1} \left( \frac{m_M}{\sigma_M^2} +\y_{i+}\cdot \frac{4k \tau_i^C \tau_i^H}{\tau_i^C + \tau_i^H}  \right) \, ,  (\tau_i^M)^{-1} \right),\\
  \Delta_i &= \mathcal{N} \left( (\tau_i^Delta)^{-1} \left( \frac{m_\Delta}{\sigma_\Delta^2} +\y_{i+}\cdot \frac{4k \tau_i^C \tau_i^H}{\tau_i^C + \tau_i^H}  \right) \, ,  (\tau_i^\Delta)^{-1} \right),\\
\end{align*}
where
\begin{align*}
  \tau_i^M & = \frac{1}{\sigma_M^2} + \frac{4k \tau^C_i \tau^H_i}}{\tau^C_i + \tau^H_i},\\
  \tau_i^\Delta & = \frac{1}{\sigma_\Delta^2} + \frac{4k \tau^C_i \tau^H_i}}{\tau^C_i + \tau^H_i},\\
\end{align*}
Because of rounding, the observed $S^2_i$ will be too small. We approximate the true $S^2_i$ by adding $\frac13$, the variance
of a uniform random variable on $[-1,1]$.
-->



Finally, we check the variance distribution empirically, to check whether the continuous distribution we have fit for individual variances describes the true distribution of variances in the population reasonably well.
The first thing we do is to compare the empirical
variances (with fractional parts imputed according to the observed proportions for the unequal digit preference, as discussed in section \ref{sec:lastdigit}) to the theoretical beta-prime distribution.
To match the standard distribution, the variances are normalized by being divided by the factor $\alpha/\theta$.
We show histograms of these "unrounded" empirical variances and the theoretical beta-prime distribution in Figure \ref{fig:histograms}.
Note that the distribution has a very long tail, and we have truncated about 2\% of the data to make the figures more readable.

In Figure \ref{fig:QQplots} we show essentially the same data in the form of Q--Q plots.
Here we have extended the plot far out into the tails of the distribution, including values in the range $[0,10]$, covering around 99.7\% of the data.
We generate data from the inferred model that mimic the true data, with three systolic and three diastolic BP measures per person.
As before, we impute the fractional parts to the real data.
This gives us a set of true variances and a set of simulated variances, which we hope will have approximately the same distribution.
We see some deviation here, but it is slight, and quite deep into the tails.
Furthermore, the deviation is in the direction of the simulated data having slightly fatter tails than the true data, which is the direction we would wish to err in for the sake of making conservative inferences.

The estimates of the empirical Bayes parameters together with their standard errors are given in the column labelled "True" in Table \ref{tab:testparameters}.
These parameters (and SEs) are accompanied by the results of 10 estimates of data simulated from the model with the parameters inferred from the data, and then fitted by the same procedure.
Note that the errors for the estimates are consistent with the stated standard errors ($\pm \sqrt{\operatorname{SE}}$), and the relative errors for the SE are small, confirming that the estimation procedure is reliable.

\begin{landscape}
\begin{figure}

{\centering \includegraphics{NHANES-BP-Appendix_files/figure-latex/histograms-1} 

}

\caption{Comparison of the distribution of empirical variances, normalized by dividing by $\beta=\alpha/\theta$, to the fitted beta-prime distribution.}(\#fig:histograms)
\end{figure}

\end{landscape}


\begin{figure}

{\centering \includegraphics{NHANES-BP-Appendix_files/figure-latex/QQplots-1} 

}

\caption{Q--Q plots of the variances of the observed data with imputed fractional parts (x-axis) against the variances of the simulated data (y-axis).}(\#fig:QQplots)
\end{figure}
<!--
Note: The empirical Bayes priors are now included in Table VarianceTest
To finish this section, we include a table of the empirical Bayes priors, shown in tables \ref{tab:empestsS} and \ref{tab:empestsD}.

\begin{table}

\caption{(\#tab:empestsS)Empirical Bayes prior hyperparameter estimates for the systolic blood pressure for the NHANES III, full population.}
\centering
\begin{tabular}[t]{ll}
\toprule
Parameter & Estimate\\
\midrule
$\alpha$ for Clinic SD & 2.57$ \pm $0.0775\\
$\alpha$ for Home SD & 2.16$ \pm $0.057\\
$\mu$ for $M$ & 123$ \pm $0.401\\
$\mu$ for $\Delta$ & 1.35$ \pm $0.241\\
$\sigma^2$ for $M$ & 19.4$ \pm $1.45\\
\addlinespace
$\sigma^2$ for $\Delta$ & 6.81$ \pm $0.862\\
$\theta$ for Clinic SD & 0.0751$ \pm $0.00105\\
$\theta$ for Home SD & 0.15$ \pm $0.00213\\
\bottomrule
\end{tabular}
\end{table}


\begin{table}

\caption{(\#tab:empestsD)Empirical Bayes prior hyperparameter estimates for the diastolic blood pressure for the NHANES III, full population.}
\centering
\begin{tabular}[t]{ll}
\toprule
Parameter & Estimate\\
\midrule
$\alpha$ for Clinic SD & 2.76$ \pm $0.0867\\
$\alpha$ for Home SD & 2.32$ \pm $0.0637\\
$\mu$ for $M$ & 72.3$ \pm $0.292\\
$\mu$ for $\Delta$ & 1.1$ \pm $0.201\\
$\sigma^2$ for $M$ & 10.2$ \pm $1.05\\
\addlinespace
$\sigma^2$ for $\Delta$ & 4.68$ \pm $0.718\\
$\theta$ for Clinic SD & 0.109$ \pm $0.00151\\
$\theta$ for Home SD & 0.194$ \pm $0.00272\\
\bottomrule
\end{tabular}
\end{table}

### Estimates for the FRS population
The estimates of the empirical Bayes parameters together with their standard errors are given in the column labelled "True" in Table \ref{tab:testparameters}.




Finally, we check the variance distribution empirically, to check whether the continuous distribution we have fit for individual variances describes the true distribution of variances in the population reasonably well.
The first thing we do is to compare the empirical
variances (with fractional parts imputed according to the observed proportions for the unequal digit preference, as discussed in section \ref{sec:lastdigit}) to the theoretical beta-prime distribution.
To match the standard distribution, the variances are normalized by being divided by the factor $\alpha/\theta$.
We show histograms of these "unrounded" empirical variances and the theoretical beta-prime distribution in Figure \ref{fig:histograms}.
Note that the distribution has a very long tail, and we have truncated about 2\% of the data to make the figures more readable.

In Figure \ref{fig:QQplots} we show essentially the same data in the form of Q--Q plots.
Here we have extended the plot far out into the tails of the distribution, including values in the range $[0,10]$, covering around 99.7\% of the data.
We generate data from the inferred model that mimic the true data, with three systolic and three diastolic BP measures per person.
As before, we impute the fractional parts to the real data.
This gives us a set of true variances and a set of simulated variances, which we hope will have approximately the same distribution.
We see some deviation here, but it is slight, and quite deep into the tails.
Furthermore, the deviation is in the direction of the simulated data having slightly fatter tails than the true data, which is the direction we would wish to err in for the sake of making conservative inferences.

The estimates of the empirical Bayes parameters together with their standard errors are given in the column labelled "True" in Table \ref{tab:testparameters}.
These parameters (and SEs) are accompanied by the results of 10 estimates of data simulated from the model with the parameters inferred from the data, and then fitted by the same procedure.
Note that the errors for the estimates are consistent with the stated standard errors ($\pm \sqrt{\operatorname{SE}}$), and the relative errors for the SE are small, confirming that the estimation procedure is reliable.

\begin{landscape}
\begin{figure}

{\centering \includegraphics{NHANES-BP-Appendix_files/figure-latex/histogramsFRS-1} 

}

\caption{Comparison of the distribution of empirical variances, normalized by dividing by $\beta=\alpha/\theta$, to the fitted beta-prime distribution.}(\#fig:histogramsFRS)
\end{figure}

\end{landscape}


\begin{figure}

{\centering \includegraphics{NHANES-BP-Appendix_files/figure-latex/QQplotsFRS-1} 

}

\caption{Q--Q plots of the variances of the observed data with imputed fractional parts (x-axis) against the variances of the simulated data (y-axis).}(\#fig:QQplotsFRS)
\end{figure}
-->

## Hamiltonian Monte Carlo (HMC)

The model, as described in the article, is a Bayesian hierarchical model. In order to parameterize such an intricate model, traditional Maximum Likelihood Estimation methods can no longer be applied. Therefore, we apply the Hamiltonian Monte Carlo (HMC) method. HMC is a form of Markov Chain Monte Carlo methods, which samples potential parameter space values of the model, then calculates directly the likelihood function based on that choice of parameters. The derivative of the likelihood function, $\phi$, guides parameter space exploration in $\theta$ towards the modal value of the joint posterior distribution. This method is ideal for complicated, non-Gaussian distribution forms. The three steps of HMC are:
\begin{enumerate}
\item Draw a sample of the derivative $\phi$ using the posterior distribution of $\phi$, which is the same as its prior.
\item Update the values of $\theta^*$ and $\phi^*$ using
  \begin{equation}
    \theta^*\leftarrow \theta+\epsilon M^{-1}\phi,
  \end{equation}
  and
  \begin{equation}
    \phi\leftarrow \phi+\epsilon\frac{1}{2}\frac{\mathrm{d}\log\{p(\theta|y)\}}{\mathrm{d}\theta},
  \end{equation}
where $M$ is the Jacobian of the parameters. This can be set to a diagonal matrix for no correlation between parameters, and is updated pointwise throughout the calculation. This is the leapfrog method, whereby $\epsilon$ dictates the scale size of the step to ensure convergence on the correct point is made, and L is the number of steps to be `leaped'.
\item Compute the rejection parameter:
  \begin{equation}
    r=\frac{p(\theta^*|y)p(\phi^*)}{p(\theta^{t-1}|y)p(\phi^{t-1})}
  \end{equation}
\item Set $\theta^t$ to $\theta^*$ with probability $\min\{1,r\}$, or otherwise keep $\theta^{t-1}$.
\end{enumerate}  
The tuning parameters $\epsilon$ and L should be chosen according to a desired acceptance rate. The No-U-Turn Sampler of Stan automates the calculation of these tuning parameters. A more detailed overview of HMC and the NUTS algorithm integrated into the Stan package, see [@NUTS].

### Centering the Linear Predictor

During the MCMC simulations, the centering values play a non-negligible role in shaping the model parameterization. If the centering parameters are held constant throughout all of the MCMC simulations, then the equation $\sum_i^N \exp{(\boldsymbol{\beta}\cdot(\boldsymbol{X}-\hat{X}))}=0$ is no longer guaranteed. However, automatically defining the centering values based on the model parameters sampled at the current MCMC iteration is not advisable as it can lead to poor parameter convergence. This is because it modifies the likelihood function at every MCMC iteration. Therefore, we iterate the MCMC algorithm multiple times. At every iteration, we recalculate the centering parameters to satisfy the requirement that the average of the linear predictor term going to zero, based on the posterior distributions of the previous MCMC simulation. This iteration is carried out until the centering parameters converge. Convergence is defined by optimising on two factors. The first is that the sum of the linear predictor term across all MCMC samples needs to tend to negligible values (we define this as the average difference being less than $10^{-7}$), see figure \ref{fig:linpred_conv}. The second convergence criteria is that the average Root Mean-Squared Error (RMSE) of the model predictions on the survival outcomes in the MCMC simulations needs to also decrease towards zero, see figure \ref{fig:linpred_conv} (top). For the second criteria, we stopped the simulations when either the difference in the RMSE stopped decreasing (below a threshold of $1\%$), or the RMSE value was less than 20, see figure \ref{fig:linpred_conv} (bottom). Illustration of the convergence is shown in figure \ref{fig:linpred_conv}.

\newpage
## Code Description
The code will be made available, but detailed references have been removed to preserve anonymity for the review process.

<!--
The code can be found at [https://github.com/hamishwp/NHANES_HPOX](https://github.com/hamishwp/NHANES_HPOX).
The numerical code has been built in multiple stages. Below, we explain the principal files required to replicate the entire analysis presented in the article. There are 5 main groups for the code:

1. Data cleaning scripts
2. Main file
3. Stan files for HMC
4. Centering recalculation scripts
5. Post-processing analysis

We provide a brief description of each of these below.

### Data Cleaning
This is found in the file `Dataclean2021.R`. Provided the raw NHANES dataset (in CSV format), it extracts all the data required for the simulations, and stores it in a structure that can be directly read in to the main file (`MCMC_DiasSyst_v3.R`) of this research. 

### Main
The main file is `MCMC_DiasSyst_v3.R`. It reads in the cleaned NHANES data, the specific choice of simulation parameters (for example, whether to use the FRS number or mean systolic & diastolic blood pressure), and runs the correct RStan scripts for that specific selection of simulation parameters. This script is intended for use on computing clusters.

### Stan

There are eight Stan files: 

1. `mystanmodel_DS_sigma_v2_autopred.stan` 
2. `mystanmodel_DS_tau_v2_autopred.stan` 
3. `mystanmodelFRS_DS_sigma_v2_autopred.stan` 
4. `mystanmodelFRS_DS_tau_v2_autopred.stan`
5. `mystanmodel_DS_sigma_v2.stan` 
6. `mystanmodel_DS_tau_v2.stan`
7. `mystanmodelFRS_DS_sigma_v2.stan`
8. `mystanmodelFRS_DS_tau_v2.stan`

These correspond to the following alternative simulation parameters:

* For the blood-pressure variability, choosing to use the standard-deviation $\sigma$ or the precision $\tau=1/\sigma$]}
* Using the FRS score or the mean diastolic and systolic blood pressure as a covariate in the analysis
* Whether the centering parameters, $\hat{X}$, in the linear predictor term are automatically calculated to satisfy $\sum_i^N \exp{(\boldsymbol{\beta}\cdot(\boldsymbol{X}-\hat{X}))}=0$ for every MCMC iteration, or whether the centering is held constant across all iterations

### Centering
The centering of the linear predictors, which is required as input to every MCMC simulation iteration, is recalculated in the files `AutoPred_Recalc.R` and `ManPred_Recalc.R`. This is then provided to the Main script, `MCMC_DiasSyst_v3.R`, which provides these centering values to the Stan code for the MCMC simulations.


### Post-processing
The post-processing script is called `PostProcessing.R`, which heavily relies on the `Functions.R` script which contains all the necessary functions to analyse the data. The post-processing script generates many useful plots of the MCMC posterior distribution for the user, including Bayes' factors, violin plots of the normalised beta and gompertz posteriors, and more.
-->


# Appendix C -- Further Results

In this section, we add some additional detail to the results section covered in the article. Extra information is given to explain how convergence of the simulations was ensured, and to also include more visualisations of the converged model parameterizations. The authors feel that this is particularly useful to provide confidence in the model parameterization and the predictions.

## Convergence of Simulations

Convergence of the simulations required to parameterize the model presented in this work is required for the MCMC simulations performed by Stan, as well as convergence in the centering values that requires repeating the Stan calculations several times. Convergence of the latter is shown in figure \ref{fig:linpred_conv}. The upper plot in figure \ref{fig:linpred_conv} illustrates convergence in the average Root Mean-Squared Error (RMSE) of the model predictions on the survival outcomes in the MCMC simulations. The lower plot in figure \ref{fig:linpred_conv} illustrates convergence in the average sum of the linear predictor terms over all MCMC chain iterations.

![Illustration of the convergence of the centering parameters of the model. ](./Rmarkdown_Plots/RMSE-Linpred_Convergence.png){#fig:linpred_conv}

With respect to convergence of the MCMC simulations, defining convergence first involves discarding the burn-in period of the simulations. When the time-evolution marker chain has a large number of samples, sequence thinning is used to reduce the amount of data storage - after convergence, take only the kth value of the simulations (after having discarded the burn-in phase values) and discard the rest. One measure of convergence is to bin similar markers and check that for each bin, the variation of the individual marker movement over a few time steps is larger than the variation of the ensemble markers in-between one-another. Other methods of convergence are stationarity and mixing. The former occurs by ensuring that the gradients of movements in the chains in time are in the same direction, the latter ensures that the amplitude of the movements in the chains are similar. To calculate the mixing and stationarity, one can do the following:
\begin{itemize}
\item Take the proposedly converged marker population, where there are N markers in total each of index length $\tau$ (thus of total physical time quantity $t\tau$). Split it k times, where k is a common denominator of $\tau$.
\item Now you have $kN$ MCMC chains each of length $|\tau/k|$
\item For the marker $\psi_{ij}$ with i and j the chain length (time) and marker number indices respectively, then the mean marker value over the chain length (time) is
  \begin{equation}
    \bar{\psi}_{|,j}=\frac{k}{\tau}\sum_{i=1}^{\tau/k}\psi_{ij}
  \end{equation}
  and the total average quantity of $\psi$ over all markers, over all chain lengths is therefore
  \begin{equation}
    \bar{\psi}_{||}=\frac{1}{kN}\sum_{j=1}^{kN}\bar{\psi}_{|j}
  \end{equation}  
\item Stationarity: compare the inter-marker variance (between sequence B):
  \begin{equation}
    B = \frac{\tau}{k(kN-1)}\sum_{j=1}^{kN}(\bar{\psi}_{|,j}-\bar{\psi}_{||})^2
  \end{equation}
\item Mixing: compare the variance along each markers chain length (within-sequence W):
  \begin{equation}
    W = \frac{1}{n(\tau-k)}\sum_{j=1}^{kN}\sum_{i=1}^{\tau/k}(\psi_{i,j}-\bar{\psi}_{|j})^2
  \end{equation}
\item Therefore, to estimate the marginal posterior variance of $p(\psi|y)$, then we use a weighted average
  \begin{equation}
    \hat{\text{Var}}^+(\psi|y)=\frac{\tau-k}{N}W+\frac{1}{Nk}B
  \end{equation}
  Note that this quantity overestimates the marginal posterior variance, but it is unbiased under stationarity: this can be used to infer convergence. When the varation in
  \begin{equation}
    \hat{R}=\sqrt{\frac{\hat{\text{Var}}^+(\psi|y)}{W}}
  \end{equation}
  should approach close to 1 for converged simulations.
\end{itemize}
Another convergence parameter is the number of effective independent marker draws. Upon convergence, the time evolution of each marker should be uncorrelated and independent to previous time steps. To find the average time-correlation over all particles, we use the variogram $V_t$:
\begin{equation}
  V_t=\frac{1}{Nk(\tau/k-\tilde{t})}\sum_{j=1}^{kN}\sum_{i=1}^{\tau/k}(\psi_{i,j}-\psi_{i-\tilde{t},j})^2,
\end{equation}
where $\tilde{t}\in 1,2,...,\tau/k$ is a time index. Then we get the time-correlations:
\begin{equation}
  \hat{\rho}_t=1-\frac{V_t}{2\hat{\text{Var}}^+}
\end{equation}
This comes from the expectation of the variance $E[(\psi_i-\psi_{i-t})^2]=2(1-\rho_t)\text{Var}(\psi)$. This can be used to infer the effective number of independent marker draws:
\begin{equation}
  \hat{n}_{eff}=\frac{mn}{1+2\sum_{\tilde{t}=1}^T\hat{\rho}_t}
\end{equation}
Where T is the index at which the sum of the autocorrelation estimates $\hat{\rho}_{t'}+\hat{\rho}_{t'+1}$ is negative. As a general guide, we should have $\hat{n}_{eff}\sim 10N/k$ effective independent marker draws and that $\hat{R}\to 1\sim 1.1$. In this research, we continued running the MCMC simulations until these two criteria were met (and went beyond: $\hat{R}<1.05$ for all parameters in all models and that $\hat{n}_{eff}>750$ for all parameters in all models).

## Results - Model Parameterization

We remind the reader of the list of numbers of the different models explored in this research, provided in the list found in section 'Proposed Models'. The authors will use the numbers in the list, referred to as the run number, in the following plots. One of the most important set of parameters of the model is the vector $\beta$ of covariates in the Cox' proportional hazards model. When the $\beta$ vector is normalised, the larger (in absolute terms) the value of $\beta$, the larger the correlation between that specific covariate and the risk of mortality. Positive values of $\beta$ imply a higher risk of mortality, and the inverse for negative values of $\beta$. As we can see from the violin plots of the MCMC posterior samples of the $\beta$ parameters in figure \ref{fig:betas}, the parameter that correlated the highest with both the mortality risk of HA-CVD-CeVD and for all mortalities, in absolute terms, was the 1998 version of the FRS score, shown in the top-right plot under run numbers 7 and 8. The FRS-1998 score correlated, on average over all the MCMC iterations, approximately $25\%$ more with mortality risk of HA-CVD-CeVD than the (more recently developed) FRS ATP III score. A similar, but slightly weaker, correlation was found between the two FRS scores for all mortality-based risk. The middle-left plot in figure \ref{fig:betas} shows that the mean diastolic blood pressure acts to decrease mortality risk. Finally, the influence of the longer-term difference in the mean blood pressure, displayed in the top-left and top-middle plots of figure \ref{fig:betas}, is also shown to increase mortality risk across all run numbers. The influence of the blood-pressure variability on mortality is illustrated to not be consistent across simulations, whereby the statistical significance of the effect is lower than for the other parameters in the linear predictor term. 



\begin{table}[!h]
\centering
\caption{(\#tab:iqrnormbetaRL1)Beta parameters for all-cause mortality, full NHANES III population, normalised by the interquartile range instead of the standard deviation.}
\centering
\begin{tabular}[t]{lcc}
\toprule
Covariate & Beta Normalised (IQR) & Range\\
\midrule
\cellcolor{gray!10}{Diastolic Mean} & \cellcolor{gray!10}{-0.077} & \cellcolor{gray!10}{(-0.097,-0.057)}\\
Systolic Mean & 0.240 & (0.22,0.26)\\
\cellcolor{gray!10}{Diastolic $|\Delta|$} & \cellcolor{gray!10}{0.052} & \cellcolor{gray!10}{(0.026,0.078)}\\
Systolic $|\Delta|$ & 0.074 & (0.055,0.092)\\
\cellcolor{gray!10}{Diastolic Clinic Stand Dev} & \cellcolor{gray!10}{-0.050} & \cellcolor{gray!10}{(-0.079,-0.022)}\\
Systolic Clinic Stand Dev & -0.011 & (-0.028,0.0052)\\
\cellcolor{gray!10}{Diastolic Home Stand Dev} & \cellcolor{gray!10}{0.019} & \cellcolor{gray!10}{(-0.0059,0.044)}\\
Systolic Home Stand Dev & -0.014 & (-0.032,0.004)\\
\bottomrule
\end{tabular}
\end{table}


\begin{table}[!h]
\centering
\caption{(\#tab:iqrnormbetaRL2)Beta parameters for cardiovascular mortality, full NHANES III population, normalised by the interquartile range instead of the standard deviation.}
\centering
\begin{tabular}[t]{lcc}
\toprule
Covariate & Beta Normalised (IQR) & Range\\
\midrule
\cellcolor{gray!10}{Diastolic Mean} & \cellcolor{gray!10}{-0.049} & \cellcolor{gray!10}{(-0.061,-0.038)}\\
Systolic Mean & 0.140 & (0.13,0.15)\\
\cellcolor{gray!10}{Diastolic $|\Delta|$} & \cellcolor{gray!10}{0.053} & \cellcolor{gray!10}{(0.038,0.068)}\\
Systolic $|\Delta|$ & 0.045 & (0.034,0.056)\\
\cellcolor{gray!10}{Diastolic Clinic Stand Dev} & \cellcolor{gray!10}{-0.019} & \cellcolor{gray!10}{(-0.034,-0.0044)}\\
Systolic Clinic Stand Dev & -0.009 & (-0.017,1.7e-05)\\
\cellcolor{gray!10}{Diastolic Home Stand Dev} & \cellcolor{gray!10}{0.012} & \cellcolor{gray!10}{(-6.0e-04,0.025)}\\
Systolic Home Stand Dev & 0.000 & (-0.0082,0.0088)\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\centering
\caption{(\#tab:iqrnormbetaRL3)Beta parameters for all-cause mortality, FRS NHANES III population, normalised by the interquartile range instead of the standard deviation.}
\centering
\begin{tabular}[t]{lcc}
\toprule
Covariate & Beta Normalised (IQR) & Range\\
\midrule
\cellcolor{gray!10}{Diastolic Mean} & \cellcolor{gray!10}{-0.022} & \cellcolor{gray!10}{(-0.05,0.0065)}\\
Systolic Mean & 0.260 & (0.23,0.28)\\
\cellcolor{gray!10}{Diastolic $|\Delta|$} & \cellcolor{gray!10}{0.060} & \cellcolor{gray!10}{(0.023,0.096)}\\
Systolic $|\Delta|$ & 0.069 & (0.039,0.098)\\
\cellcolor{gray!10}{Diastolic Clinic Stand Dev} & \cellcolor{gray!10}{0.003} & \cellcolor{gray!10}{(-0.037,0.043)}\\
Systolic Clinic Stand Dev & -0.027 & (-0.055,0.0013)\\
\cellcolor{gray!10}{Diastolic Home Stand Dev} & \cellcolor{gray!10}{0.021} & \cellcolor{gray!10}{(-0.01,0.052)}\\
Systolic Home Stand Dev & 0.032 & (0.0063,0.059)\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\centering
\caption{(\#tab:iqrnormbetaRL4)Beta parameters for cardiovascular mortality, FRS NHANES III population, normalised by the interquartile range instead of the standard deviation.}
\centering
\begin{tabular}[t]{lcc}
\toprule
Covariate & Beta Normalised (IQR) & Range\\
\midrule
\cellcolor{gray!10}{Diastolic Mean} & \cellcolor{gray!10}{-0.019} & \cellcolor{gray!10}{(-0.034,-0.0036)}\\
Systolic Mean & 0.140 & (0.13,0.16)\\
\cellcolor{gray!10}{Diastolic $|\Delta|$} & \cellcolor{gray!10}{0.054} & \cellcolor{gray!10}{(0.035,0.072)}\\
Systolic $|\Delta|$ & 0.043 & (0.028,0.059)\\
\cellcolor{gray!10}{Diastolic Clinic Stand Dev} & \cellcolor{gray!10}{0.007} & \cellcolor{gray!10}{(-0.012,0.026)}\\
Systolic Clinic Stand Dev & -0.013 & (-0.027,4.9e-04)\\
\cellcolor{gray!10}{Diastolic Home Stand Dev} & \cellcolor{gray!10}{0.005} & \cellcolor{gray!10}{(-0.011,0.021)}\\
Systolic Home Stand Dev & 0.017 & (0.0028,0.031)\\
\bottomrule
\end{tabular}
\end{table}

![Violin plots of the normalised $\beta$ parameters of the different models. ](./Plots/beta/Beta_parameter_normalised.png){#fig:betas}

With respect to the time-independent Gompertz parameter, denoted $B$ in this article, the results between all models that simulate CVD mortality risk, and all the models that simulation all-cause mortality risk are consistent with one-another. This is illustrated by the similarity between plots on the left hand side and the right hand side of figure \ref{fig:gompB}. The consistency appears across sex assigned at birth and race.

![Violin plots of the normalised B parameter (from the Gompertz equation) of the different models.](./Plots/gompertz/B_parameter.png){#fig:gompB}

Figure \ref{fig:gompt} reflects the same level of consistency for the Gompertz parameter that influences the temporal evolution of the mortality risk. It is worth noting that both figures \ref{fig:gompB} and \ref{fig:gompt} have inverse trends between the values of B and theta for each demographic group. This makes it difficult to imagine, based on these two plots, what the mortality risk is at different ages across demographics, yet it is evident that the form of the change in the mortality risk curve in time is different for each demographic group. Women are observed to have lower initial values of risk, but mortality risk later in life begins to increase much faster than for men. Additionally, hispanic populations are shown to have a larger initial mortality risk than black populations who are shown to have a larger initial mortality risk than white populations in the USA. However, mortality risk increases at a faster rate for white populations than for black populations, for which it increases faster than hispanic populations in the USA. For ease of comparison, we also present here tables of the mean and standard deviation values of the time dependent and independent Gompertz parameters in tables \ref{tab:RL12} to \ref{tab:RL78}.

![Violin plots of the normalised $\theta$ parameter (from the Gompertz equation) of the different models.](./Plots/gompertz/theta_parameter.png){#fig:gompt}


\begin{table}[!h]
\centering
\caption{(\#tab:Mean-SD)Summary of the posterior estimate of the blood pressure distribution.}
\centering
\begin{tabular}[t]{lllll}
\toprule
\multicolumn{1}{c}{ } & \multicolumn{2}{c}{Systolic} & \multicolumn{2}{c}{Diastolic} \\
\cmidrule(l{3pt}r{3pt}){2-3} \cmidrule(l{3pt}r{3pt}){4-5}
\addlinespace[0.3em]
\multicolumn{5}{l}{\textbf{Full population}}\\
\cellcolor{gray!10}{\hspace{1em}Overall Mean} & \cellcolor{gray!10}{125.4} & \cellcolor{gray!10}{19.5} & \cellcolor{gray!10}{74.3} & \cellcolor{gray!10}{10.3}\\
\hspace{1em}$|\Delta|=|\mathrm{Home}-\mathrm{Clinic}|/2$ & 5.24 & 4.83 & 3.90 & 3.14\\
\cellcolor{gray!10}{\hspace{1em}Home Stand Dev} & \cellcolor{gray!10}{2.74} & \cellcolor{gray!10}{2.05} & \cellcolor{gray!10}{2.34} & \cellcolor{gray!10}{1.75}\\
\hspace{1em}Clinic Stand Dev & 3.78 & 2.61 & 3.08 & 2.08\\
\addlinespace[0.3em]
\multicolumn{5}{l}{\textbf{FRS population}}\\
\cellcolor{gray!10}{\hspace{1em}Overall Mean} & \cellcolor{gray!10}{125.9} & \cellcolor{gray!10}{18.3} & \cellcolor{gray!10}{76.4} & \cellcolor{gray!10}{10.0}\\
\hspace{1em}$|\Delta|=|\mathrm{Home}-\mathrm{Clinic}|/2$ & 5.23 & 4.71 & 3.84 & 3.10\\
\cellcolor{gray!10}{\hspace{1em}Home Stand Dev} & \cellcolor{gray!10}{2.78} & \cellcolor{gray!10}{2.06} & \cellcolor{gray!10}{2.28} & \cellcolor{gray!10}{1.71}\\
\hspace{1em}Clinic Stand Dev & 3.78 & 2.52 & 2.92 & 1.97\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\centering
\caption{(\#tab:RL12)Parameters for survival model, NHANES III, Full population, using the systolic and diastolic mean model.}
\centering
\begin{tabular}[t]{llrrrr}
\toprule
Sex & Race/Ethnicity & B–Mean & B–SD & $\theta$-Mean & $\theta$–SD\\
\midrule
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{All-Cause Mortality}}\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{1.62e-04} & \cellcolor{gray!10}{3.14e-05} & \cellcolor{gray!10}{0.0720} & \cellcolor{gray!10}{0.00264}\\
\hspace{1em}Female & White & 3.12e-05 & 6.20e-06 & 0.0906 & 0.00238\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{2.66e-04} & \cellcolor{gray!10}{5.45e-05} & \cellcolor{gray!10}{0.0625} & \cellcolor{gray!10}{0.00279}\\
\hspace{1em}Male & Black & 4.29e-04 & 7.54e-05 & 0.0641 & 0.00243\\
\cellcolor{gray!10}{\hspace{1em}Male} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{5.08e-05} & \cellcolor{gray!10}{9.70e-06} & \cellcolor{gray!10}{0.0900} & \cellcolor{gray!10}{0.00236}\\
\hspace{1em}Male & Mexican & 5.27e-04 & 8.88e-05 & 0.0590 & 0.00245\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{CVD Mortality}}\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{1.36e-05} & \cellcolor{gray!10}{5.70e-06} & \cellcolor{gray!10}{0.0883} & \cellcolor{gray!10}{0.00534}\\
\hspace{1em}Female & White & 2.20e-06 & 9.00e-07 & 0.1080 & 0.00488\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{2.44e-05} & \cellcolor{gray!10}{1.07e-05} & \cellcolor{gray!10}{0.0775} & \cellcolor{gray!10}{0.00565}\\
\hspace{1em}Male & Black & 7.20e-05 & 2.49e-05 & 0.0712 & 0.00459\\
\cellcolor{gray!10}{\hspace{1em}Male} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{8.70e-06} & \cellcolor{gray!10}{3.40e-06} & \cellcolor{gray!10}{0.0973} & \cellcolor{gray!10}{0.00466}\\
\hspace{1em}Male & Mexican & 8.17e-05 & 3.00e-05 & 0.0678 & 0.00498\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\centering
\caption{(\#tab:RL34)Parameters for survival model, NHANES III, FRS-population only, using the systolic and diastolic mean model.}
\centering
\begin{tabular}[t]{llrrrr}
\toprule
Sex & Race/Ethnicity & B–Mean & B–SD & $\theta$-Mean & $\theta$–SD\\
\midrule
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{All-Cause Mortality}}\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{1.77e-04} & \cellcolor{gray!10}{4.66e-05} & \cellcolor{gray!10}{0.0699} & \cellcolor{gray!10}{0.00359}\\
\hspace{1em}Female & White & 3.12e-05 & 8.80e-06 & 0.0902 & 0.00355\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{1.56e-04} & \cellcolor{gray!10}{4.83e-05} & \cellcolor{gray!10}{0.0696} & \cellcolor{gray!10}{0.00427}\\
\hspace{1em}Male & Black & 3.41e-04 & 8.54e-05 & 0.0662 & 0.00352\\
\cellcolor{gray!10}{\hspace{1em}Male} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{4.27e-05} & \cellcolor{gray!10}{1.10e-05} & \cellcolor{gray!10}{0.0902} & \cellcolor{gray!10}{0.00335}\\
\hspace{1em}Male & Mexican & 4.51e-04 & 1.17e-04 & 0.0594 & 0.00369\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{CVD Mortality}}\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{8.60e-06} & \cellcolor{gray!10}{5.40e-06} & \cellcolor{gray!10}{0.0937} & \cellcolor{gray!10}{0.00795}\\
\hspace{1em}Female & White & 1.70e-06 & 1.30e-06 & 0.1100 & 0.00844\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{4.06e-05} & \cellcolor{gray!10}{2.65e-05} & \cellcolor{gray!10}{0.0708} & \cellcolor{gray!10}{0.00879}\\
\hspace{1em}Male & Black & 6.36e-05 & 3.53e-05 & 0.0697 & 0.00765\\
\cellcolor{gray!10}{\hspace{1em}Male} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{1.08e-05} & \cellcolor{gray!10}{6.50e-06} & \cellcolor{gray!10}{0.0896} & \cellcolor{gray!10}{0.00737}\\
\hspace{1em}Male & Mexican & 1.94e-04 & 9.69e-05 & 0.0539 & 0.00709\\
\bottomrule
\end{tabular}
\end{table}


\begin{table}[!h]
\centering
\caption{(\#tab:RL78)Parameters for survival model, NHANES III, FRS-population only, using the 1998 FRS-based model.}
\centering
\begin{tabular}[t]{llrrrr}
\toprule
Sex & Race/Ethnicity & B–Mean & B–SD & $\theta$-Mean & $\theta$–SD\\
\midrule
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{All-Cause Mortality}}\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{2.34e-04} & \cellcolor{gray!10}{7.08e-05} & \cellcolor{gray!10}{0.0684} & \cellcolor{gray!10}{0.00404}\\
\hspace{1em}Female & White & 3.66e-05 & 1.16e-05 & 0.0897 & 0.00392\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{1.88e-04} & \cellcolor{gray!10}{6.25e-05} & \cellcolor{gray!10}{0.0694} & \cellcolor{gray!10}{0.00442}\\
\hspace{1em}Male & Black & 4.94e-04 & 1.42e-04 & 0.0640 & 0.00384\\
\cellcolor{gray!10}{\hspace{1em}Male} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{5.01e-05} & \cellcolor{gray!10}{1.46e-05} & \cellcolor{gray!10}{0.0903} & \cellcolor{gray!10}{0.00365}\\
\hspace{1em}Male & Mexican & 6.00e-04 & 1.74e-04 & 0.0580 & 0.00391\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{CVD Mortality}}\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{2.62e-05} & \cellcolor{gray!10}{1.82e-05} & \cellcolor{gray!10}{0.0843} & \cellcolor{gray!10}{0.00859}\\
\hspace{1em}Female & White & 3.70e-06 & 2.70e-06 & 0.1030 & 0.00861\\
\cellcolor{gray!10}{\hspace{1em}Female} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{1.02e-04} & \cellcolor{gray!10}{7.46e-05} & \cellcolor{gray!10}{0.0641} & \cellcolor{gray!10}{0.00896}\\
\hspace{1em}Male & Black & 2.14e-04 & 1.34e-04 & 0.0608 & 0.00793\\
\cellcolor{gray!10}{\hspace{1em}Male} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{2.22e-05} & \cellcolor{gray!10}{1.46e-05} & \cellcolor{gray!10}{0.0859} & \cellcolor{gray!10}{0.00778}\\
\hspace{1em}Male & Mexican & 5.58e-04 & 3.25e-04 & 0.0461 & 0.00720\\
\bottomrule
\end{tabular}
\end{table}

## Results - Model Performance

To measure the performance of the model to predict the survival outcome of individuals in the population, figure \ref{fig:cumpred} shows, ordered by individual age, the cumulative hazard $H(t)$ predicted against the cumulative number of deaths in the populations, for each model explored in this research. Each model is shown to predict survival outcomes reliably, across the entire age range of the population.

![Predicted cumulative hazard against cumulative number of deaths in the population, ordered by the age of the individual. ](./Plots/Survival/redlinpred_Cumulative_haz-death_age.png){#fig:cumpred}

A common metric that is used to evaluate the performance of models such as presented in this article is called the Receiver Operating Characteristic (ROC) curve. With continuous predictor values such as cumulative hazard $H(T_i)$, a threshold can be defined whereby any individual who has a cumulative risk larger than the threshold $H(T_i)>\epsilon$ is predicted to die. The ratio of the number individuals that were predicted to die compared to the total number who die corresponds is referred to as the True Positive Ratio (TPR)
\begin{equation}
  TPR(\epsilon)=\frac{\sum_i\big(\mathbb{I}(H(T_i)>\epsilon \ \ \& \ \ \delta_i=1)\big)}{\sum_i\big(\mathbb{I}(\delta_i=1)\big)}.
\end{equation}
Note that TPR is also referred to as the recall or sensitivity. Conversely, the ratio of the number of individuals predicted to die but survive compared to the total number of individuals that survived is referred to as the False Positive Ratio (FPR)
\begin{equation}
  FPR(\epsilon)=\frac{\sum_i\big(\mathbb{I}(H(T_i)>\epsilon \ \ \& \ \ \delta_i=0)\big)}{\sum_i\big(\mathbb{I}(\delta_i=0)\big)}.
\end{equation}
Note that the FPR is also referred to as $1-\mathrm{specificity}$. An ROC curve is produced by varying the threshold value that is then used to calculate both the TPR and FPR, and plotting them against one another.
The area under this curve is a metric that indicates performance of the model to predict survival outcomes. AUROC=1 implies perfect predictions and AUROC=0.5 implies the contrary. However, our model is formulated such that the variables age and time since starting the survey both form part of Cox's proportional hazard. Furthermore, the Gompertz model is stratified by demographic group. Therefore, in this work, we present a modified ROC curve, which calculates the individuals cumulative hazard at a given time since the start of the survey, $T_{surv} \in {5, 10, 15}$ years, and calculates whether the model correctly predicted an event to occur before or after this time. 
Note that to do this, we split the ROC population by ages: 45-64 and 65-84. The modified TPR is then calculated via:
\begin{equation}
  \operatorname{TPR}(\epsilon)=\frac{\sum_i\big(\mathbb{I}(\delta_i=1 \ \ \& \ \ \ H(T_i)\geq \epsilon \ \ \& \ \ T_i<T_{surv})\big)}{\sum_i\big(\mathbb{I}(\delta_i=1 \ \ \& \ \ \ T_i<T_{surv})\big)},
\end{equation}
and the modified FPR:
\begin{equation}
  \operatorname{FPR}(\epsilon)=\frac{\sum_i\big(\mathbb{I}(H(T_i)\geq \epsilon \ \ \& \ \ T_i\geq T_{surv})\big)}{\sum_i\big(\mathbb{I}(T_i\geq T_{surv})\big)}.
\end{equation}

One problem with standard ROC curves is that they are stretched out along the diagonal of the plot, making it impossible to scale the actual differences to make the significant deviations from the diagonal more legible.
For that reason it is often better to plot a variant ROC curve, where the ordinate (the vertical coordinate) is the difference TPR--FPR, rather than TPR itself.
This modified plot --- also known as the Youden index plot --- rotates what was formerly the diagonal onto the abscissa, allowing the deviations to be plotted on the appropriate scale.
The information content is, of course, the same, and the standard AUC is now simply $0.5+$the area under this modified plot.

Before presenting any plots or AUC values, we first present the density distributions of the median posterior systolic $\Delta$ values for all individuals, split by demographic. 

![Density of the (median posterior) systolic $\Delta$ values, per demographic. ](./Rmarkdown_Plots/SysDelta_Densities_Demography.png){#fig:DeltaDens}

Using Welch's ANOVA test, we calculated that $p<1\times 10^{-9}$ for all demographics, including when split between the male and female populations. Figures \ref{fig:ROC_MeanBP} and \ref{fig:ROC_FRS} show the Youden index plots (modified ROC curves) --- including the AUC values --- of the model, the former using the mean systolic and diastolic blood pressure as covariates in the linear predictor term (in the Cox's proportional hazards component) and the latter using the FRS value instead. By making predictions of the 5, 10 and 15 year survival between the middle aged and old aged sub-groups, for the three different mortality causes, we start to build a picture of the performance of the model. For figure \ref{fig:ROC_MeanBP}, we notice that the AUC value of the predictions for the middle aged compared to the older aged population is higher, independent of the survival year prediction or the different mortality causes. The highest AUC is for the 45-64 year old population with a focus on CVD and heart attack-related mortality, for all three survival year periods. We also note that the TPR seems to start increasing at a faster rate for the population aged 45-64 than the 65-84 group, implying that it is possible to choose a threshold level, $\epsilon$, for the survival predictions that could correctly identify people at risk without incorrectly predicting as many people to be at risk of mortality as for the group aged 65-84. The results also reflect that the influence of choosing a 5, 10 or 15 year prediction period does not seem to significantly influence the results. Figure \ref{fig:ROC_FRS} displays similar results to the mean systolic and diastolic model when using the FRS value instead, with the main difference that the predictions of the middle aged group for CVD and heart attack-related mortality for 5 year survival seems to be lower than the equivalent in the older group or as compared to the mean blood pressure model. This is caused by a reduced mortality before 5 years for the middle aged population that had their FRS value calculated, where 36, 83 and 213 CVD and heart attack-related deaths occurred before 5, 10 and 15 years in this sub-group, respectively. This can be compared to 141, 356 and 952 all-cause deaths in this same sub-group (red-curve in figure \ref{fig:ROC_FRS} and \ref{fig:ROC_MeanBP}). Alternatively, when compared to the full-population (not just those who had the FRS value), the CVD and heart attack-related deaths for the middle aged population are 55, 115 and 282 over the 5, 10 and 15 year range, respectively.

![Youden index plots (modified ROC curves) for the model that used mean systolic and diastolic blood pressure as covariates in the linear predictor, stratified by the event type (cause of mortality). The columns split two groups in the population: those who start the survey aged between 45 to 64 and 65-84 years old. The rows split the model predictions between 5, 10 and 15 year survival. ](./Rmarkdown_Plots/ROC_MeanBPModel_CAx-EventType.png){#fig:ROC_MeanBP}

![Youden index plots (modified ROC curves) for the model that used the FRS value as covariates in the linear predictor, stratified by age group and the number of years the survival outcome was predicted since participant starting the survey. ](./Rmarkdown_Plots/ROC_FRSModel_CAx-EventType.png){#fig:ROC_FRS}

Comparison of the Youden index plots (modified ROC curves) and AUC values is also presented for the different demographic groups, see figure \ref{fig:ROC_Demog}. This figure shows the differences in the prediction performance (w.r.t. the ROC and AUC values)  of the full-population model for the 45-64 age population for their 10 year survival outcome, using the mean systolic and diastolic blood pressure model. This plot illustrates that potentially only the all-cause mortality has enough outcomes in each demographic group to separate the curves. The model seems to most accurately predict the 10 year survival outcome of the black and other ethnic groups, as well as the white female demographic as compared to the black female, other female and white male population. To provide insight into this, we also provide the frequency table of deaths for each demographic group, mortality cause and survival year, see table \ref{tab:DeathFreq1} and \ref{tab:DeathFreq2}.


![Youden index plots (modified ROC curves) stratified by the different demographic groups used in this research. The point and line colours represent the different event types that were used to predict on. ](./Rmarkdown_Plots/ROC_CAx-EventType_Demog_10Yr_45-64.png){#fig:ROC_Demog}

\begin{table}[!h]
\centering
\caption{(\#tab:DeathFreq1)Frequency table of the population aged 45-64 for the N-year survival outcomes as separated by demographic group and mortality cause.}
\centering
\begin{tabular}[t]{llllr}
\toprule
Year & EventType & Ethnicity & Gender & Deaths\\
\midrule
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{15}\\
5-Year Mortality & All Deaths & Black & Female & 44\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{9}\\
5-Year Mortality & All Deaths & White & Female & 33\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{9}\\
5-Year Mortality & All Deaths & Mexican & Female & 24\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{6}\\
5-Year Mortality & All Deaths & Black & Male & 25\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{7}\\
5-Year Mortality & All Deaths & White & Male & 24\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{6}\\
5-Year Mortality & All Deaths & Mexican & Male & 21\\
\addlinespace\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{23}\\
10-Year Mortality & All Deaths & Black & Female & 87\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{21}\\
10-Year Mortality & All Deaths & White & Female & 67\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{19}\\
10-Year Mortality & All Deaths & Mexican & Female & 63\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{17}\\
10-Year Mortality & All Deaths & Black & Male & 67\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{13}\\
10-Year Mortality & All Deaths & White & Male & 62\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{10}\\
10-Year Mortality & All Deaths & Mexican & Male & 44\\
\addlinespace\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{37}\\
15-Year Mortality & All Deaths & Black & Female & 132\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{44}\\
15-Year Mortality & All Deaths & White & Female & 139\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{35}\\
15-Year Mortality & All Deaths & Mexican & Female & 105\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{26}\\
15-Year Mortality & All Deaths & Black & Male & 113\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{29}\\
15-Year Mortality & All Deaths & White & Male & 117\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{18}\\
15-Year Mortality & All Deaths & Mexican & Male & 67\\
\bottomrule
\end{tabular}
\end{table}


\begin{table}[!h]
\centering
\caption{(\#tab:DeathFreq2)Frequency table of the population aged 65-84 for the N-year survival outcomes as separated by demographic group and mortality cause.}
\centering
\begin{tabular}[t]{llllr}
\toprule
Year & EventType & Ethnicity & Gender & Deaths\\
\midrule
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{35}\\
5-Year Mortality & All Deaths & Black & Female & 87\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{92}\\
5-Year Mortality & All Deaths & White & Female & 253\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{15}\\
5-Year Mortality & All Deaths & Mexican & Female & 45\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{19}\\
5-Year Mortality & All Deaths & Black & Male & 49\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{55}\\
5-Year Mortality & All Deaths & White & Male & 138\\
\cellcolor{gray!10}{5-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{10}\\
5-Year Mortality & All Deaths & Mexican & Male & 29\\
\addlinespace\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{56}\\
10-Year Mortality & All Deaths & Black & Female & 153\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{187}\\
10-Year Mortality & All Deaths & White & Female & 501\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{37}\\
10-Year Mortality & All Deaths & Mexican & Female & 105\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{50}\\
10-Year Mortality & All Deaths & Black & Male & 121\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{144}\\
10-Year Mortality & All Deaths & White & Male & 357\\
\cellcolor{gray!10}{10-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{44}\\
10-Year Mortality & All Deaths & Mexican & Male & 91\\
\addlinespace\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{79}\\
15-Year Mortality & All Deaths & Black & Female & 231\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{251}\\
15-Year Mortality & All Deaths & White & Female & 694\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Female} & \cellcolor{gray!10}{56}\\
15-Year Mortality & All Deaths & Mexican & Female & 170\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Black} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{76}\\
15-Year Mortality & All Deaths & Black & Male & 185\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{White} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{225}\\
15-Year Mortality & All Deaths & White & Male & 589\\
\cellcolor{gray!10}{15-Year Mortality} & \cellcolor{gray!10}{CVD} & \cellcolor{gray!10}{Mexican} & \cellcolor{gray!10}{Male} & \cellcolor{gray!10}{57}\\
15-Year Mortality & All Deaths & Mexican & Male & 139\\
\bottomrule
\end{tabular}
\end{table}


To finalise the section on the performance of the models using ROC and AUC values, we present a series of figures that provide Youden index plots and AUC values for different linear predictor (Cox's proportional hazards model) covariate configurations. By setting the covariate-specific $\beta$ parameter values to zero, we can measure the additional prediction performance that adding different covariates provides to the model. In figures \ref{fig:ROC_RL1} to \ref{fig:ROC_RL8oth}, we present three main formulations: using only the systolic and diastolic $\Delta$ terms, using the mean systolic and diastolic (figures \ref{fig:ROC_RL1}-\ref{fig:ROC_RL2oth}) or the FRS value (figures \ref{fig:ROC_RL7}-\ref{fig:ROC_RL8oth}) terms as well as the systolic and diastolic $\Delta$ terms (figures \ref{fig:ROC_RL1}-\ref{fig:ROC_RL8oth}) and, finally, using the systolic mean (figures \ref{fig:ROC_RL1}-\ref{fig:ROC_RL2oth}) or FRS value only (figures \ref{fig:ROC_RL7}-\ref{fig:ROC_RL8oth}). Figures \ref{fig:ROC_RL1}, \ref{fig:ROC_RL2} and \ref{fig:ROC_RL2oth} apply to the full-population with the models trained on CVD and heart-attack related mortality, all-cause, and other mortality, respectively. Figures \ref{fig:ROC_RL7}, \ref{fig:ROC_RL8} and \ref{fig:ROC_RL8oth} apply to the population with an FRS value, with the models trained on CVD and heart-attack related mortality, all-cause, and other mortality, respectively. The first thing to note as a commonality between all these different figures is that the use of the long-term variability, $\Delta$, in the model has comparable performance with that of using the systolic mean or FRS values only. Additionally, where the number of deaths permits for prediction, the use of both the $\Delta$ and mean/FRS values results in higher AUC values than using models that use one or the other. Finally, the use of the FRS value consistently under-performs the mean diastolic and systolic blood pressure-based model.

![Youden index plots (modified ROC curves) for the mean systolic and diastolic model, looking specifically at CVD and heart attack-related deaths, stratified by age group and the number of years the survival outcome was predicted since participant starting the survey. The colour of the points and lines represents the different linear predictor covariate models possible. ](./Rmarkdown_Plots/ROC_CAx-Covariates_EventType_RL1.png){#fig:ROC_RL1}

![Youden index plots (modified ROC curves) for the mean systolic and diastolic model, looking at all-cause deaths, stratified by age group and the number of years the survival outcome was predicted since participant starting the survey. The colour of the points and lines represents the different linear predictor covariate models possible. ](./Rmarkdown_Plots/ROC_CAx-Covariates_EventType_RL2.png){#fig:ROC_RL2}

![Youden index plots (modified ROC curves) for the mean systolic and diastolic model, looking non-CVD and heart attack-related deaths, stratified by age group and the number of years the survival outcome was predicted since participant starting the survey. The colour of the points and lines represents the different linear predictor covariate models possible. ](./Rmarkdown_Plots/ROC_CAx-Covariates_EventType_RL2oth.png){#fig:ROC_RL2oth}

![Youden index plots (modified ROC curves) for for the FRS-based model, looking specifically at CVD and heart attack-related deaths, stratified by age group and the number of years the survival outcome was predicted since participant starting the survey. The colour of the points and lines represents the different linear predictor covariate models possible. ](./Rmarkdown_Plots/ROC_CAx-Covariates_EventType_RL7.png){#fig:ROC_RL7}

![Youden index plots (modified ROC curves) for for the FRS-based model, looking at all-cause deaths, stratified by age group and the number of years the survival outcome was predicted since participant starting the survey. The colour of the points and lines represents the different linear predictor covariate models possible. ](./Rmarkdown_Plots/ROC_CAx-Covariates_EventType_RL8.png){#fig:ROC_RL8}

![Youden index plots (modified ROC curves) for for the FRS-based model, looking non-CVD and heart attack-related deaths, stratified by age group and the number of years the survival outcome was predicted since participant starting the survey. The colour of the points and lines represents the different linear predictor covariate models possible. ](./Rmarkdown_Plots/ROC_CAx-Covariates_EventType_RL8oth.png){#fig:ROC_RL8oth}

## Results - Excluded Population

As there were 3,916 individuals that were excluded from the research, it is important to verify that there is no bias in the survival outcomes linked to each of the demographic groups for the included and excluded populations. This section demonstrates this by considering differences in the Kaplan-Meier curves, shown in figure \ref{fig:excpop}. No significant difference in survival outcomes is observed between the different demographic groups for the excluded and included populations. However, it is important to note that bias has been introduced by the fact that the 751 individuals that were included in the 'other' ethnicity will therefore not be represented in this research.


![Kaplan-Meier plots of the full-population for all-cause mortality, split by demographic group and whether the individual was included in this research. ](./Rmarkdown_Plots/Excluded_SurvProbKM.png){#fig:excpop}

## Results - Frequentist Comparison

A comparison with the frequentist Cox proportional hazards model is presented in this section, to allow the reader to percieve the differences between the Bayesian hierarchical model and the standard frequentist approach. Tables \ref{tab:freqCVDF} to \ref{tab:freqALLNF} show the model results from using only the blood pressure covariates to predict survival outcomes.
The baseline mortality is stratified by sex and ethnicity.






``` r
cox_NHANES <- coxph( Surv(age, Time+age, eventCVDHrt) ~ FRS.1998 +
                     D_i_S + D_i_D + sigma_C_S + sigma_H_S +
                     sigma_C_D + sigma_H_D + strata(race) + strata(female), data = DF_nhanesFRS)
```

\begin{table}[!h]
\centering
\caption{(\#tab:freqCVDNF)Cox-PH parameter estimates for cardiovascular mortality, NHANES III, full population.}
\centering
\begin{tabular}[t]{lrll}
\toprule
Covariate & Beta Normalised & P-Value & Range\\
\midrule
\cellcolor{gray!10}{Systolic Mean} & \cellcolor{gray!10}{0.312} & \cellcolor{gray!10}{<0.001} & \cellcolor{gray!10}{(0.25,0.37)}\\
Diastolic Mean & -0.097 & 0.0028 & (-0.16,-0.033)\\
\cellcolor{gray!10}{Systolic $|\Delta|$} & \cellcolor{gray!10}{0.065} & \cellcolor{gray!10}{0.0025} & \cellcolor{gray!10}{(0.023,0.11)}\\
Diastolic $|\Delta|$ & 0.053 & 0.038 & (0.0029,0.11)\\
\cellcolor{gray!10}{Systolic Clinic Stand Dev} & \cellcolor{gray!10}{-0.008} & \cellcolor{gray!10}{0.72} & \cellcolor{gray!10}{(-0.054,0.036)}\\
Systolic Home Stand Dev & -0.017 & 0.46 & (-0.062,0.028)\\
\cellcolor{gray!10}{Diastolic Clinic Stand Dev} & \cellcolor{gray!10}{-0.027} & \cellcolor{gray!10}{0.3} & \cellcolor{gray!10}{(-0.08,0.025)}\\
Diastolic Home Stand Dev & 0.025 & 0.34 & (-0.027,0.08)\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\centering
\caption{(\#tab:freqALLNF)Cox-PH parameter estimates for all-cause mortality, NHANES III, FRS population.}
\centering
\begin{tabular}[t]{lrll}
\toprule
Covariate & Beta Normalised & P-Value & Range\\
\midrule
\cellcolor{gray!10}{Systolic Mean} & \cellcolor{gray!10}{0.166} & \cellcolor{gray!10}{<0.001} & \cellcolor{gray!10}{(0.13,0.19)}\\
Diastolic Mean & -0.051 & 0.0048 & (-0.088,-0.015)\\
\cellcolor{gray!10}{Systolic $|\Delta|$} & \cellcolor{gray!10}{0.047} & \cellcolor{gray!10}{<0.001} & \cellcolor{gray!10}{(0.022,0.07)}\\
Diastolic $|\Delta|$ & 0.053 & <0.001 & (0.026,0.083)\\
\cellcolor{gray!10}{Systolic Clinic Stand Dev} & \cellcolor{gray!10}{-0.009} & \cellcolor{gray!10}{0.49} & \cellcolor{gray!10}{(-0.036,0.017)}\\
Systolic Home Stand Dev & -0.004 & 0.74 & (-0.03,0.022)\\
\cellcolor{gray!10}{Diastolic Clinic Stand Dev} & \cellcolor{gray!10}{-0.005} & \cellcolor{gray!10}{0.72} & \cellcolor{gray!10}{(-0.033,0.023)}\\
Diastolic Home Stand Dev & 0.013 & 0.38 & (-0.017,0.044)\\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\centering
\caption{(\#tab:freqCVDF)Cox-PH parameter estimates for cardiovascular mortality, NHANES III, FRS population.}
\centering
\begin{tabular}[t]{lrll}
\toprule
Covariate & Beta Normalised & P-Value & Range\\
\midrule
\cellcolor{gray!10}{FRS (1998)} & \cellcolor{gray!10}{0.512} & \cellcolor{gray!10}{<0.001} & \cellcolor{gray!10}{(0.34,0.67)}\\
Systolic $|\Delta|$ & 0.108 & <0.001 & (0.045,0.18)\\
\cellcolor{gray!10}{Diastolic $|\Delta|$} & \cellcolor{gray!10}{0.102} & \cellcolor{gray!10}{0.0057} & \cellcolor{gray!10}{(0.029,0.18)}\\
Systolic Clinic Stand Dev & 0.019 & 0.59 & (-0.05,0.089)\\
\cellcolor{gray!10}{Systolic Home Stand Dev} & \cellcolor{gray!10}{0.060} & \cellcolor{gray!10}{0.085} & \cellcolor{gray!10}{(-0.0082,0.13)}\\
Diastolic Clinic Stand Dev & 0.025 & 0.49 & (-0.046,0.096)\\
\cellcolor{gray!10}{Diastolic Home Stand Dev} & \cellcolor{gray!10}{0.028} & \cellcolor{gray!10}{0.46} & \cellcolor{gray!10}{(-0.045,0.1)}\\
\bottomrule
\end{tabular}
\end{table}


\begin{table}[!h]
\centering
\caption{(\#tab:freqALLF)Cox-PH parameter estimates for all-cause mortality, NHANES III, FRS population.}
\centering
\begin{tabular}[t]{lrll}
\toprule
Covariate & Beta Normalised & P-Value & Range\\
\midrule
\cellcolor{gray!10}{FRS (1998)} & \cellcolor{gray!10}{0.223} & \cellcolor{gray!10}{<0.001} & \cellcolor{gray!10}{(0.14,0.3)}\\
Systolic $|\Delta|$ & 0.081 & <0.001 & (0.05,0.12)\\
\cellcolor{gray!10}{Diastolic $|\Delta|$} & \cellcolor{gray!10}{0.067} & \cellcolor{gray!10}{<0.001} & \cellcolor{gray!10}{(0.032,0.11)}\\
Systolic Clinic Stand Dev & 0.007 & 0.68 & (-0.027,0.042)\\
\cellcolor{gray!10}{Systolic Home Stand Dev} & \cellcolor{gray!10}{0.017} & \cellcolor{gray!10}{0.35} & \cellcolor{gray!10}{(-0.019,0.054)}\\
Diastolic Clinic Stand Dev & 0.023 & 0.19 & (-0.012,0.06)\\
\cellcolor{gray!10}{Diastolic Home Stand Dev} & \cellcolor{gray!10}{0.012} & \cellcolor{gray!10}{0.54} & \cellcolor{gray!10}{(-0.027,0.05)}\\
\bottomrule
\end{tabular}
\end{table}





## Results - Exploring $\Delta$ Directionality

This section of the appendix is to explore whether the directionality of the difference in clinic-home blood pressure (represented through the non-absolute value of the $\Delta$ covariate) may have an influence on the survival outcome in the population. In the work presented in this article, $\Delta$ is the absolute value of the differences in the means of the blood pressure measurements at the clinic and at home, respectively, for both diastolic and systolic blood pressure. By 'directionality', we refer to whether the difference between the clinic and home mean measurements are positive or negative. Figure \ref{fig:DeltaDensities} shows the clinic-home directionalities, split by demographic group, indicating no significant difference between the different demographic groups. There is a general trend that the directionality for systolic and diastolic blood pressure is more likely to be the same than opposite.

![The range of the non-absolute $\Delta$ values in the systolic and diastolic blood pressure measurements, split by demographic group. This reflects the differences between the average measurements at the clinic and at home. ](./Rmarkdown_Plots/Delta_plusminus_Demography.png){#fig:DeltaDensities}

In order to explore whether the directionality of the clinic-home measurements influences survival outcome, we will use a combination of Kaplan-Meier curves and Cox's proportional hazards regression. The latter will implement a simple Maximum Likelihood Estimation (MLE) method based on summary statistics of the Bayesian posterior blood pressure values, not the Bayesian-HMC method applied elsewhere in this article. The Kaplan-Meier curve is a plot of the change in survival probability of a population in time since the start of a survey/census. The survival distribution is calculated using
\begin{equation}\label{survKM}
\hat{S}(t)=\prod_{t_j \le t}\left(1-\frac{d_j}{r_j} \right),
\end{equation}
for $d_j$ the number of individuals who die within the time interval $t_j$ and $r_j$ the population that are alive (at risk of death) and not censored. Greenwood's formula is used to calculate the variance of the Kaplan-Meier estimation
\begin{equation}\label{sigKM}
\hat{\sigma}(t)^2=\hat{S}(t)^2\sum_{t_j \le t}\left(\frac{d_j}{r_j(r_j-d_j)} \right).
\end{equation}
The 100(1-$\alpha$)\% confidence intervals of the Kaplan-Meier estimate are assumed to be normally distributed
\begin{equation}\label{CIKM}
\hat{S}(t) \pm z_{1-\alpha/2}\hat{\sigma}(t).
\end{equation}

Figures \ref{fig:KM45tot} and \ref{fig:KM65tot} show the Kaplan-Meier estimates for the full NHANES population for CVD mortality, split by demographic group, for the age range 45-64 and 65-84, respectively. The survival probability of the older population decreases faster in time than the middle-aged (45-64) population.

![Kaplan-Meier plots of the full-population for CVD mortality, for ages between 45-64, split by demographic group. ](./Rmarkdown_Plots/SurvProbKM_45-64.png){#fig:KM45tot}

![Kaplan-Meier plots of the full-population for CVD mortality, for ages between 65-84, split by demographic group. ](./Rmarkdown_Plots/SurvProbKM_65-84.png){#fig:KM65tot}

By splitting the populations into the respective regions of the $\Delta$ directionality, we can use the different Kaplan-Meier plots to try to identify differences in the survival outcomes. Figures \ref{fig:KM45_deltaregion} and \ref{fig:KM65_deltaregion} show the Kaplan-Meier estimates for the full NHANES population for CVD mortality, split by $\Delta$ directionality and demographic group, for the age range 45-64 and 65-84, respectively. With the broad confidence intervals, all of the different $\Delta$ directionality regions overlap, for all demographic groups and both age range groups (where $\hat{S}(t)\neq 1$).

![Kaplan-Meier plots of the full-population for CVD mortality, for ages between 45-64, split by demographic group and region in systolic-diastolic $\Delta$ space. ](./Rmarkdown_Plots/SurvProbKM_Delta_45-65.png){#fig:KM45_deltaregion}

![Kaplan-Meier plots of the full-population for CVD mortality, for ages between 65-84, split by demographic group and region in systolic-diastolic $\Delta$ space. ](./Rmarkdown_Plots/SurvProbKM_Delta_65-85.png){#fig:KM65_deltaregion}

We further quantify this insignificant relationship between $\Delta$ directionality and survival outcome via the use of a Cox's Proportional Hazards (CPH) model. Via the use of the 'coxph' function in the 'survival' R package, we fit (using MLE) a CPH model. The covariates used in the model are $\Delta$ directionality region, gender, ethnicity and age. Table \ref{tab:DeltaDir} shows the summary of the model fit, which reflects that being in different $\Delta$ directionality regions has a non-statistically significant influence on survival outcomes. As shown in the remainder of this article, ethnicity, gender and age are shown to have statistically significant effects on survival outcome.

\begin{table}[!h]
\centering
\caption{(\#tab:DeltaDir)Parameters for distribution of blood pressure, for the full population}
\centering
\begin{tabular}[t]{lrrrrr}
\toprule
covariate & coef & exp(coef) & se(coef) & z & Pr(>|z|)\\
\midrule
\cellcolor{gray!10}{DeltaRegionSys. -ve, Dys. +ve} & \cellcolor{gray!10}{0.047} & \cellcolor{gray!10}{1.048} & \cellcolor{gray!10}{0.074} & \cellcolor{gray!10}{0.634} & \cellcolor{gray!10}{0.526}\\
DeltaRegionSys. +ve, Dys. -ve & -0.002 & 0.998 & 0.069 & -0.031 & 0.975\\
\cellcolor{gray!10}{DeltaRegionSys. +ve, Dys. +ve} & \cellcolor{gray!10}{-0.056} & \cellcolor{gray!10}{0.946} & \cellcolor{gray!10}{0.063} & \cellcolor{gray!10}{-0.881} & \cellcolor{gray!10}{0.378}\\
GenderMale & -0.413 & 0.662 & 0.050 & -8.282 & 0.000\\
\cellcolor{gray!10}{EthnicityWhite} & \cellcolor{gray!10}{-0.248} & \cellcolor{gray!10}{0.780} & \cellcolor{gray!10}{0.061} & \cellcolor{gray!10}{-4.043} & \cellcolor{gray!10}{0.000}\\
EthnicityMexican & -0.170 & 0.844 & 0.075 & -2.272 & 0.023\\
\cellcolor{gray!10}{age} & \cellcolor{gray!10}{0.100} & \cellcolor{gray!10}{1.105} & \cellcolor{gray!10}{0.002} & \cellcolor{gray!10}{50.573} & \cellcolor{gray!10}{0.000}\\
\bottomrule
\end{tabular}
\end{table}

Finally, we wish to confirm that the perfomance of the model does not depend on the directionality of $\Delta$. Figure \ref{fig:DeltaAUCs} plots the AUC values of 10-year CVD mortality for the all-covariate mean blood pressure-based model trained on the full NHANES population, split by the two age-ranges (45-64 and 65-84) and demographic groups. There is no clear trend between mode performance for the different regions of $\Delta$ directionality.

![AUC values of the full-population, all covariate model with the mean blood pressure model, based on CVD 10-year mortality. ](./Rmarkdown_Plots/DeltaDirection_AUCs.png){#fig:DeltaAUCs}

# References