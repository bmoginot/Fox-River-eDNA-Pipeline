# Stats
## Kruskal-Wallis Test
nonparametic (does not assume normality) one-way ANOVA

### Assumptions
1. ordinal or continuous response variable
2. independent
3. distributions in each group have a similar shape (???)

### One-Way ANOVA
compares the means of three or more independent groups to determine if there is a statistically significant difference between population means

#### Assumptions
1. normality (Shapiro-Wilk test i think??)
2. equal variance (Bartlett's test; variance test for multiple groups, i.e. ANOVA)
3. independence (we have random samples)

## Wilcoxon Rank Sum Test
nonparametric equivalent to two-sample independent t-test. compare between two groups whose sampling distribution is not normal and the sample sizes are small. it has the same assumptions as the Kruskal-Wallis test.

### Sum the Ranks
order values for both groups from least to greatest. assign ranks to values. lowest value is 1, highest is n1 + n2 where n is sample size. for equal values, give them the half-way point between the ranks, e.g., rank 7 and 8 are both 2.5, so they are both ranked 7.5. sum the ranks for values from group 1 and do the same for group 2. plug each R into test stat formula and take the lower of the two as the test stat.

## Bonferroni Test for Multiple Testing
alter alpha value based on the number of tests being performed

### Formula
new_alpha = original_alpha / n
original_alpha is usually 0.5. n is the number of tests being performed.

## Workflow (???)
check normality of data using Shapiro-Wilk (dr. picq's was not, so i'm assuming my data will not be normal either)
one-way ANOVA (Kruskal-Wallis) on three categories (tributary, park river, outside park river)
Wilcoxon rank sum test (pairwise, each category against against each other category)
Bonferroni correction to determine new alpha to consider for each pairwise Wilcoxon test