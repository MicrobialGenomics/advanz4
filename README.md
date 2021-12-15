## ADVANZ-4 PROJECT

### End Points to Include:
****

1. CD4+ count changes at week 48.  

2. Changes in CD4+ counts at 96 weeks, 

3. Proportion of subjects with >200 CD4+ counts/mm3 at weeks 48 and 96

4. Changes in immune activation (CD8+CD38+, HLA-DR), Bacterial translocation (sCD14), inflammation (IL-6, TNF-alpha and ultrasensitive PCR) and coagulation (D-Dimer) markers at weeks 0, 24, 48 and 96. 

6. Percentage of subjects with detectable fecal ADV, CMV, and enterovirus overall, and by study arm, at weeks 0, 24, 48 and 96 

7. The association of the gut microbiome composition and function and intestinal shedding of ADV, CMV, and enterovirus, with the percent of subjects at week 48 and 96 with: 

    - CD4+ count increases >50cells/mm3
    
    - CD4+ counts > 200 cells/mm3
    
    - CD4+ counts > 500 cells/mm3
    
    - Ratio CD4+/CD8+ > 1 
    
8. Association of the gut microbiome composition and function and intestinal shedding of ADV, CMV and Enterovirus with immune activation, inflammation, coagulation and bacterial translocation markers. 

### Analysis Plan:
****

1. **Cohort Description:**

Summarize follow up, dropouts and possible Adverse Events. Make sure there are no errors in the metadata. The immune state analysis have been mostly done and already published. We'll focus more on the microbiome part. 

2. **Microbiome analysis:**

    - Elaborate a concrete definition for the many different comparison the project asks for:
    
        a. CD4+ count increases: >50 at week 48 &	>50 at week 96 
        
        b. Absolute CD2: two new vairables. Week 48/96: <200 | >200 & < =500 | > 500 
        
    - **IDEA:** Longitudinal part of the analysis has been done already, and immune/inflammation markers have already been proven to strongly correlate with TimePoint. This means we could transform the TimePoint into a variable, which consists on a factor composed by the cumulativa number of patients who meet each criteria (CD4 increase, CD4count, etc), for each timepoint. This variable could be inputted into an LMM against viral counts/species counts etc. The idea is assuming immune markers are strongly correlated with treatment time hence, if there is any microbiomic/viromic change being affected by such markers, they should be affected as well, together with immune markers. However, it should be wise to check the correlation between TimePoint and immune markers to be sure our assumption is correct. 
    
        a. Species Composition: **Metaphlan3** 
        
            - Stacked Barplots 
            
            - NMDS and Bray-curtis first Differences 
            
            - Gene Richness 
            
            - DA (Wilcoxon), w48/46 

        b. Pathways: **HumaNn3** 
        
            - Stacked Barplots 
            
            - Beta diversity: NMDS and Bray-Curtis First Differences  
            
            - Gene Richness 
            
            - DA (Wilcoxon), w48/46 
