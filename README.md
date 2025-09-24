# mantodea_supmat_coue
Supplementary material for "Analyzing potential niche shifts in alien pairs of mantis species (Insecta, Mantodea), with comments on the current taxonomic and ecological knowledge"

Equivalent of:
https://zenodo.org/records/10101384

angui_intra: intra-species analyses for Tenodera angustipennis, with both complete and trimmed dataset.

contrhierodula: intra-species analyses for the Hierodula species with controversial status (H. transcaucasica/H. tenuidentata). It includes analyses between H. transcaucasica and H. tenuidentata sensu stricto and the populations treated as separated species.

hierodula_europe: inter-species analysis with the two Hierodula present in Europe.

kendall: scripts and table used for calculating Kendall's rank correlation.

pate_intra: intra-species analyses for Hierodula patellifera, with both complete and trimmed (no India, Nepal and Philippines points) datasets. There is also the analyses between Asian distribution and India, Nepal and Philippines points

sinensis_intra: intra-species analyses for Tenodera sinesis, with both complete and trimmed datasets.

tenodera_europe: inter-species analysis with the two Tenodera present in North America, with both complete and trimmed datasets.

####NOTE

I did not include the raw presence points (the elaborated ones are together with the values in a csv file per each folder) and the wordlclim variables. For the latter, you can get them from the WorldClim website. If you want the raw presence points, I will give them on request

WARNING: AN ITALIAN LOCALE HAS BEEN USED FOR THIS ANALYSES; IF YOU HAVE ANY FORMAT ISSUE WITH THE CSV, CONSIDER USING THIS FOR IMPORTING INTO R:

read_csv("namefile.csv", locale = locale(decimal_mark = ",", grouping_mark = "."))

 

IF THAT DOES NOT WORK, PLEASE CONTACT ME
