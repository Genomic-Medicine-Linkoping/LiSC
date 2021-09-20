#!/usr/bin/env Rscript

R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS:   /opt/sw/bioinfo-tools/sources/R-4.0.3/lib/libRblas.so
LAPACK: /opt/sw/bioinfo-tools/sources/R-4.0.3/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=sv_SE.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=sv_SE.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=sv_SE.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=sv_SE.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] stringr_1.4.0      data.table_1.14.0  cowplot_1.1.1      RColorBrewer_1.1-2
 [5] ggplot2_3.3.5      patchwork_1.1.1    scCATCH_2.1        SeuratObject_4.0.2
 [9] Seurat_4.0.4       dplyr_1.0.7        docopt_0.7.1      

loaded via a namespace (and not attached):
  [1] nlme_3.1-152          matrixStats_0.60.0    spatstat.sparse_2.0-0
  [4] progress_1.2.2        RcppAnnoy_0.0.19      httr_1.4.2           
  [7] sctransform_0.3.2     tools_4.0.3           utf8_1.2.2           
 [10] R6_2.5.1              irlba_2.3.3           rpart_4.1-15         
 [13] KernSmooth_2.23-20    uwot_0.1.10           mgcv_1.8-36          
 [16] DBI_1.1.1             lazyeval_0.2.2        colorspace_2.0-2     
 [19] withr_2.4.2           prettyunits_1.1.1     tidyselect_1.1.1     
 [22] gridExtra_2.3         compiler_4.0.3        cli_3.0.1            
 [25] plotly_4.9.4.1        labeling_0.4.2        scales_1.1.1         
 [28] lmtest_0.9-38         spatstat.data_2.1-0   ggridges_0.5.3       
 [31] pbapply_1.4-3         goftest_1.2-2         digest_0.6.27        
 [34] spatstat.utils_2.2-0  pkgconfig_2.0.3       htmltools_0.5.1.1    
 [37] parallelly_1.27.0     limma_3.46.0          fastmap_1.1.0        
 [40] htmlwidgets_1.5.3     rlang_0.4.11          rstudioapi_0.13      
 [43] shiny_1.6.0           farver_2.1.0          generics_0.1.0       
 [46] zoo_1.8-9             jsonlite_1.7.2        ica_1.0-2            
 [49] magrittr_2.0.1        Matrix_1.3-4          Rcpp_1.0.7           
 [52] munsell_0.5.0         fansi_0.5.0           abind_1.4-5          
 [55] reticulate_1.20       lifecycle_1.0.0       stringi_1.7.3        
 [58] MASS_7.3-54           Rtsne_0.15            plyr_1.8.6           
 [61] grid_4.0.3            parallel_4.0.3        listenv_0.8.0        
 [64] promises_1.2.0.1      ggrepel_0.9.1         crayon_1.4.1         
 [67] miniUI_0.1.1.1        deldir_0.2-10         lattice_0.20-44      
 [70] splines_4.0.3         tensor_1.5            hms_1.1.0            
 [73] pillar_1.6.2          igraph_1.2.6          spatstat.geom_2.2-2  
 [76] future.apply_1.8.1    reshape2_1.4.4        codetools_0.2-18     
 [79] leiden_0.3.9          glue_1.4.2            png_0.1-7            
 [82] vctrs_0.3.8           httpuv_1.6.2          gtable_0.3.0         
 [85] RANN_2.6.1            purrr_0.3.4           spatstat.core_2.3-0  
 [88] polyclip_1.10-0       tidyr_1.1.3           scattermore_0.7      
 [91] future_1.21.0         assertthat_0.2.1      mime_0.11            
 [94] xtable_1.8-4          RSpectra_0.16-0       later_1.3.0          
 [97] survival_3.2-12       viridisLite_0.4.0     tibble_3.1.3         
[100] cluster_2.1.2         globals_0.14.0        fitdistrplus_1.1-5   
[103] ellipsis_0.3.2        ROCR_1.0-11          

