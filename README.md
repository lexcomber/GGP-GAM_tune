# Spatially varying coefficient modelling with a Geographical Gaussian Process GAM (GGP-GAM)
Alexis Comber<sup>1*</sup>, Paul Harris<sup>2</sup> , Daisuke Murakami<sup>3</sup> , Tomoki Nakaya<sup>4</sup> , Naru Tsutsumida<sup>5</sup>, Takahiro Yoshida<sup>6</sup>  and Chris Brunsdon<sup>7</sup> 

<sup>1</sup> School of Geography, University of Leeds, Leeds, UK.\
<sup>2</sup> Department of Information and Computer Sciences, Saitama University, Japan\

<sup>2</sup> Sustainable Agriculture Sciences, Rothamsted Research, North Wyke, UK\
<sup>3</sup> Institute of Statistical Mathematics, Japan\
<sup>4</sup> Graduate School of Environmental Studies, Tohoku University, Japan\
<sup>5</sup> Department of Information and Computer Sciences, Saitama University, Japan\
<sup>6</sup> Center for Spatial Information Science, University of Tokyo, Japan\
<sup>7</sup> National Centre for Geocomputation, Maynooth University, Ireland\
<sup>*</sup> contact author: a.comber@leeds.ac.uk

## Abstract
This paper describes a novel spatially varying coefficient (SVC) regression model constructed from a Generalized Additive Model (GAM) with Gaussian Process (GP) splines parameterised at observation locations: a GGP-GAM. The ability of the GGP-GAM to estimate true spatially varying coefficients was compared with that of Multiscale Geographically Weighted Regression (MGWR) using simulated data with complex spatial heterogeneities. The GGP-GAM out-performs MGWR. The MGWR poorer performance was investigated and found to be due to relatively higher residuals. This may be due to the inability of MGWR and other kernel based approaches, to handle highly localised patterns of spatial heterogeneity. One of the GGP-GAM models was investigated in detail to illustrate GAM diagnostics, model checks, spline smooth convergence and basis evaluations, and the number of spline knots. A larger case study ($n = 250000$) was investigated to determine how acceptable trade-offs between GGP-GAM model complexity, performance and computational can be identified. A number of areas of further work are discussed.

This paper will be submitted to Geographic Analysis. 

## Code
To run the analysis in this paper you should download the the R script `GGP-GAM_tune.R`, the 5 data files unless you want to run the analysis from fresh (about 24 hrs). Package and other info is below. The data files and supporting scripts will need will need to be locally available . The code recreates the results in the same sequence as in the paper. 

If you have any problems with data / code / versions etc please contact Lex Comber at the email above.
```{r}
sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.6.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] spdep_1.2-5       doParallel_1.0.17 iterators_1.0.14  foreach_1.5.2     GWmodel_2.2-9     spatialreg_1.2-5 
 [7] Matrix_1.5-1      spData_2.2.0      Rcpp_1.0.10       robustbase_0.95-0 maptools_1.1-4    sp_1.6-0         
[13] mgcv_1.8-40       nlme_3.1-159      cols4all_0.6      cowplot_1.1.1     sf_1.0-9          spmoran_0.2.2.6  
[19] forcats_0.5.2     stringr_1.5.0     dplyr_1.1.0       purrr_1.0.1       readr_2.1.2       tidyr_1.3.0      
[25] tibble_3.1.8      ggplot2_3.4.1     tidyverse_1.3.2  

loaded via a namespace (and not attached):
 [1] googledrive_2.0.0   colorspace_2.1-0    deldir_1.0-6        ellipsis_0.3.2      class_7.3-20        fs_1.6.1           
 [7] proxy_0.4-27        farver_2.1.1        RSpectra_0.16-1     fansi_1.0.4         lubridate_1.8.0     xml2_1.3.3         
[13] codetools_0.2-18    splines_4.2.0       knitr_1.42          spam_2.9-1          jsonlite_1.8.4      broom_1.0.1        
[19] cluster_2.1.4       dbplyr_2.2.1        png_0.1-8           compiler_4.2.0      httr_1.4.4          backports_1.4.1    
[25] assertthat_0.2.1    gargle_1.2.1        cli_3.6.0           s2_1.1.2            tools_4.2.0         dotCall64_1.0-1    
[31] coda_0.19-4         gtable_0.3.1        glue_1.6.2          wk_0.7.1            maps_3.4.0          gmodels_2.18.1.1   
[37] cellranger_1.1.0    raster_3.6-14       vctrs_0.5.2         gdata_2.18.0.1      xfun_0.37           rvest_1.0.3        
[43] lifecycle_1.0.3     gtools_3.9.3        googlesheets4_1.0.1 terra_1.7-3         DEoptimR_1.0-11     zoo_1.8-11         
[49] LearnBayes_2.15.1   MASS_7.3-58.1       scales_1.2.1        hms_1.1.2           expm_0.999-6        RColorBrewer_1.1-3 
[55] fields_14.1         gridExtra_2.3       stringi_1.7.12      e1071_1.7-13        permute_0.9-7       boot_1.3-28        
[61] intervals_0.15.2    rlang_1.0.6         pkgconfig_2.0.3     lattice_0.20-45     labeling_0.4.2      tidyselect_1.2.0   
[67] magrittr_2.0.3      R6_2.5.1            generics_0.1.3      DBI_1.1.3           pillar_1.8.1        haven_2.5.1        
[73] foreign_0.8-82      withr_2.5.0         xts_0.12.1          units_0.8-1         abind_1.4-5         spacetime_1.2-8    
[79] modelr_0.1.9        crayon_1.5.2        rARPACK_0.11-0      KernSmooth_2.23-20  utf8_1.2.3          tzdb_0.3.0         
[85] viridis_0.6.2       grid_4.2.0          readxl_1.4.1        FNN_1.1.3.1         vegan_2.6-4         reprex_2.0.2       
[91] classInt_0.4-8      munsell_0.5.0       viridisLite_0.4.1  
```
