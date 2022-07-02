# Knowledge-Based Adjustment for Multiple Colinearity

Some of the clinical variables should be transformed and combined rather than modeling on its own. Here are the scientific evidence to justify these transformations.

## BMI
`X_imp$BMI <- X_imp$B_WEIGHT / (X_imp$B_HEIGHT / 100)^2`

## MBP
`X_imp$MBP <- X_imp$Diastolic_blood_pressure + 1/3 * (X_imp$Systolic_blood_pressure - X_imp$Diastolic_blood_pressure)`

[Body surface area, height, and body fat percentage as more sensitive risk factors of cancer and cardiovascular disease](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7300397/)

## NLR
`X_imp$NLR <- X_imp$Total_Neutrophils / X_imp$Lymphocytes`

or alternatively
`X_imp$NLR <- X_imp$ANC / X_imp$WBC`

[Neutrophil to lymphocyte ratio and breast cancer risk: analysis by subtype and potential interactions](https://www.nature.com/articles/s41598-020-70077-z)

[High Neutrophil to Lymphocyte Ratio (NLR) is Associated with Treatment Failure and Death in Melanoma Patients Treated with PD-1 Inhibitor Monotherapy](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6906249/)

[Neutrophil to lymphocyte ratio and cancer prognosis: an umbrella review of systematic reviews and meta-analyses of observational studies](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01817-1)

[Prognostic value of neutrophil-to-lymphocyte ratio and platelet-to-lymphocyte ratio in gastric cancer](https://journals.lww.com/md-journal/Fulltext/2018/03230/Prognostic_value_of_neutrophil_to_lymphocyte_ratio.14.aspx)

[A High Preoperative Neutrophil–to–lymphocyte Ratio Is Associated with Poor Survival in Patients with Colorectal Cancer](https://ar.iiarjournals.org/content/33/8/3291)

[The neutrophil-to-lymphocyte ratio: a narrative review](https://ecancer.org/en/journal/article/702-the-neutrophil-to-lymphocyte-ratio-a-narrative-review)

## PWR

`X_imp$PWR <- X_imp$`

[Platelet-to-lymphocyte ratio and lymphocyte-to-white blood cell ratio predict the efficacy of neoadjuvant chemotherapy and the prognosis of locally advanced gastric cancer patients treated with the oxaliplatin and capecitabine regimen](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6200072/)

[White blood cell and platelet indices as prognostic markers in patients with invasive ductal breast carcinoma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950562/)

## AGR
`X_imp$AGR <- X_imp$Albumin / X_imp$Hemoglobin`

[Prognostic value of the albumin–globulin ratio and albumin–globulin score in patients with multiple myeloma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7944530/ )

## HALP composite
Composite index, HALP, calculated as hemoglobin (g/L) × albumin (g/L) × lymphocytes (/L) / platelets (/L) [oncotarget](https://www.oncotarget.com/article/12271/text/)

## BUN/Cre
[The blood urea nitrogen/creatinine (BUN/cre) ratio was U-shaped associated with all-cause mortality in general population](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8856064/)


## UC(/Cre)
First report, not reliable marker

[Α higher ratio of serum uric acid to serum creatinine could predict the risk of total and cause specific mortality- insight from a US national survey](https://pubmed.ncbi.nlm.nih.gov/32535029/)

It's a good idea to use uric acid on its own. Hyperuricemia is an indication of systematic inflammation.

[Friend or Foe? An Unrecognized Role of Uric Acid in Cancer Development and the Potential Anticancer Effects of Uric Acid-lowering Drugs](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7378935/)

[Contribution of uric acid to cancer risk, recurrence, and mortality](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3560981/)

## AST/LDH

[468P Prognostic role of aspartate aminotransferase-to-alanine aminotransferase ratio and lactate dehydrogenase levels in resectable colorectal cancer](https://www.annalsofoncology.org/article/S0923-7534(21)03218-X/fulltext)

[Extremely Elevated Lactate Dehydrogenase and Aspartate Aminotransferase Level as an Indicator of Tumor Lysis during Bortezomib Treatment in a Multiple Myeloma Patient: Case Report](https://ashpublications.org/blood/article/108/11/5123/129444/Extremely-Elevated-Lactate-Dehydrogenase-and)


## AST/ALT
A classic marker, should be deployed before AST/LDH
[Preoperative serum liver enzyme markers for predicting early recurrence after curative resection of hepatocellular carcinoma](https://pubmed.ncbi.nlm.nih.gov/25865691/)

## Cl:Na ratio
Indicator of metabolic acidosis
`X_imp$`

[The value of the chloride: Sodium ratio in differentiating the aetiology of metabolic acidosis](https://link.springer.com/article/10.1007/s001340100915)

## SERMFE/HGB
[Prognostic Value of Ferritin-to-Hemoglobin Ratio in Patients with Advanced Non-Small-Cell Lung Cancer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6548010/)